"""
Anchor Graph for monomer decomposition.

Key idea: Build a graph where:
- Nodes = conserved parts (anchors) from FS-tree hints
- Edges = transitions between anchors with support (copy number)

This solves two problems:
1. Multiple conserved parts per monomer (overcutting)
2. Mutations in anchor (missing monomers) - handled by using shorter anchors
   where longer ones are absent

The graph structure reveals the monomer architecture.
"""

from collections import defaultdict, Counter
from typing import List, Tuple, Dict, Set, Optional
from dataclasses import dataclass, field


@dataclass
class AnchorHit:
    """Single occurrence of an anchor in sequence."""
    anchor: str
    position: int
    length: int

    def __lt__(self, other):
        return self.position < other.position


def get_top_candidates(candidates: List[Dict], top_k: int = 10) -> List[Dict]:
    """
    Get top-k candidates from compute_cuts output.

    Candidates are already scored by compute_cuts logic.
    """
    # Sort by adjusted_score (or base_score if adjusted not available)
    sorted_candidates = sorted(
        candidates,
        key=lambda x: x.get('adjusted_score', x.get('base_score', 0)),
        reverse=True
    )
    return sorted_candidates[:top_k]


def find_anchor_hits_hierarchical(
    sequence: str,
    anchors: List[str],
    verbose: bool = False
) -> List[AnchorHit]:
    """
    Find all anchor occurrences, using longer anchors first.

    Logic:
    1. Sort anchors by length (longest first)
    2. Find all positions of longest anchor
    3. For shorter anchors - only find in regions NOT covered by longer ones
    4. This handles mutations: if long anchor has mutation, short variant will be found

    Args:
        sequence: The array sequence
        anchors: List of anchor sequences (will be sorted by length)
        verbose: Print debug info

    Returns:
        List of AnchorHit sorted by position
    """
    if not anchors:
        return []

    # Sort by length descending
    sorted_anchors = sorted(anchors, key=len, reverse=True)

    # Track covered positions
    covered = set()
    all_hits = []

    for anchor in sorted_anchors:
        anchor_len = len(anchor)

        # Find all occurrences
        start = 0
        hits_for_anchor = 0

        while True:
            pos = sequence.find(anchor, start)
            if pos == -1:
                break

            # Check if this position is already covered by a longer anchor
            # A position is "covered" if ANY part of this anchor overlaps with covered region
            is_covered = any(p in covered for p in range(pos, pos + anchor_len))

            if not is_covered:
                all_hits.append(AnchorHit(anchor=anchor, position=pos, length=anchor_len))
                # Mark positions as covered
                for p in range(pos, pos + anchor_len):
                    covered.add(p)
                hits_for_anchor += 1

            start = pos + 1

        if verbose and hits_for_anchor > 0:
            print(f"  '{anchor}' (len={anchor_len}): {hits_for_anchor} hits")

    # Sort by position
    all_hits.sort()

    return all_hits


def build_transition_graph(hits: List[AnchorHit]) -> Dict:
    """
    Build graph of transitions between consecutive anchor hits.

    Returns dict with:
    - edges: {(from_anchor, to_anchor): count}
    - distances: {(from_anchor, to_anchor): [list of distances]}
    - sequences: {(from_anchor, to_anchor): [list of sequences between]}
    """
    edges = Counter()
    distances = defaultdict(list)

    for i in range(len(hits) - 1):
        curr = hits[i]
        next_hit = hits[i + 1]

        edge = (curr.anchor, next_hit.anchor)
        distance = next_hit.position - curr.position

        edges[edge] += 1
        distances[edge].append(distance)

    return {
        "edges": dict(edges),
        "distances": dict(distances),
    }


def find_monomer_cycle(graph: Dict, verbose: bool = False) -> List[str]:
    """
    Find the dominant cycle in the transition graph.

    This represents the sequence of anchors in one monomer.

    Strategy:
    1. First check for high-quality self-loops (single anchor per monomer, consistent distance)
    2. If no good self-loop, try multi-anchor cycles (multiple conserved parts per monomer)
    3. For multi-anchor cycles, ensure anchors are not overlapping variants
    """
    edges = graph["edges"]
    distances = graph.get("distances", {})

    if not edges:
        return []

    # Separate self-loops from cross-anchor edges
    self_loops = {e: c for e, c in edges.items() if e[0] == e[1]}
    cross_edges = {e: c for e, c in edges.items() if e[0] != e[1]}

    # First, evaluate best self-loop quality
    best_self_loop = None
    best_self_score = 0
    best_self_cv = float('inf')

    for (anchor, _), count in self_loops.items():
        edge = (anchor, anchor)
        if edge in distances:
            dists = distances[edge]
            if len(dists) >= 3:  # Need enough data points
                mean_dist = sum(dists) / len(dists)
                variance = sum((d - mean_dist) ** 2 for d in dists) / len(dists)
                cv = (variance ** 0.5) / mean_dist if mean_dist > 0 else float('inf')
                # Score: high count + low CV (consistency)
                score = count * (1 / (1 + cv * 2))

                if score > best_self_score:
                    best_self_score = score
                    best_self_cv = cv
                    best_self_loop = anchor

    # If we have a good self-loop (consistent distances), use it
    # CV < 0.3 means distances are very consistent (good monomer periodicity)
    if best_self_loop and best_self_cv < 0.3:
        if verbose:
            print(f"  High-quality self-loop: {best_self_loop[:20]}... (CV={best_self_cv:.3f})")
        return [best_self_loop]

    # Otherwise, try to find multi-anchor cycle (different conserved parts)
    if cross_edges:
        # Filter out edges between overlapping anchors (>50% overlap)
        filtered_cross_edges = {}
        for (from_a, to_a), count in cross_edges.items():
            # Check if anchors share a common substring
            min_len = min(len(from_a), len(to_a))
            is_overlapping = False
            for offset in range(max(1, min_len - 5)):
                if from_a[offset:] == to_a[:len(from_a)-offset]:
                    is_overlapping = True
                    break
                if to_a[offset:] == from_a[:len(to_a)-offset]:
                    is_overlapping = True
                    break
            if not is_overlapping:
                filtered_cross_edges[(from_a, to_a)] = count

        if filtered_cross_edges:
            # Find multi-anchor cycle using filtered edges
            best_cross_edge = max(filtered_cross_edges.items(), key=lambda x: x[1])[0]
            start_anchor = best_cross_edge[0]

            cycle = [start_anchor]
            current = start_anchor
            visited_anchors = {start_anchor}

            for _ in range(100):
                best_next = None
                best_count = 0

                for (from_a, to_a), count in filtered_cross_edges.items():
                    if from_a == current and count > best_count:
                        if to_a == start_anchor and len(cycle) > 1:
                            if verbose:
                                print(f"  Multi-anchor cycle: {' -> '.join([a[:15]+'...' for a in cycle])}")
                            return cycle
                        if to_a not in visited_anchors:
                            best_next = to_a
                            best_count = count

                if best_next is None:
                    for (from_a, to_a), count in filtered_cross_edges.items():
                        if from_a == current and to_a == start_anchor and len(cycle) > 1:
                            if verbose:
                                print(f"  Multi-anchor cycle: {' -> '.join([a[:15]+'...' for a in cycle])}")
                            return cycle
                    break

                cycle.append(best_next)
                visited_anchors.add(best_next)
                current = best_next

            if len(cycle) > 1:
                if verbose:
                    print(f"  Partial multi-anchor: {' -> '.join([a[:15]+'...' for a in cycle])}")
                return cycle

    # Fall back to best self-loop even if not super consistent
    if best_self_loop:
        if verbose:
            print(f"  Self-loop fallback: {best_self_loop[:20]}... (CV={best_self_cv:.3f})")
        return [best_self_loop]

    # Last resort: most frequent edge
    if edges:
        start_edge = max(edges.items(), key=lambda x: x[1])[0]
        if verbose:
            print(f"  Fallback to most frequent: {start_edge[0][:20]}...")
        return [start_edge[0]]

    return []


def estimate_monomer_length(graph: Dict, cycle: List[str]) -> Tuple[float, float]:
    """
    Estimate monomer length from cycle and transition distances.

    Returns (mean_length, variance)
    """
    if not cycle:
        return 0, 0

    distances = graph["distances"]

    # Sum distances along the cycle
    cycle_lengths = []

    # We need to find instances where the full cycle occurs
    # For now, estimate from individual edge distances
    total_distances = []

    for i in range(len(cycle)):
        from_anchor = cycle[i]
        to_anchor = cycle[(i + 1) % len(cycle)]
        edge = (from_anchor, to_anchor)

        if edge in distances:
            total_distances.extend(distances[edge])

    if not total_distances:
        return 0, 0

    # Rough estimate: sum of mean distances along cycle
    edge_means = []
    for i in range(len(cycle)):
        from_anchor = cycle[i]
        to_anchor = cycle[(i + 1) % len(cycle)]
        edge = (from_anchor, to_anchor)

        if edge in distances:
            edge_means.append(sum(distances[edge]) / len(distances[edge]))

    if edge_means:
        estimated_length = sum(edge_means)
        variance = sum((d - estimated_length/len(cycle))**2 for d in total_distances) / len(total_distances)
        return estimated_length, variance

    return 0, 0


def decompose_by_cycle(
    sequence: str,
    hits: List[AnchorHit],
    cycle: List[str],
    verbose: bool = False
) -> List[str]:
    """
    Decompose sequence into monomers based on the anchor cycle.

    Strategy:
    1. Find positions where the cycle starts (first anchor of cycle)
    2. Cut at those positions

    Args:
        sequence: The array sequence
        hits: Sorted list of AnchorHit
        cycle: List of anchors representing one monomer

    Returns:
        List of monomer sequences
    """
    if not cycle or not hits:
        return [sequence] if sequence else []

    start_anchor = cycle[0]

    # Find all positions where cycle starts
    cut_positions = []

    i = 0
    while i < len(hits):
        if hits[i].anchor == start_anchor:
            # Check if this is start of a valid cycle
            # (following anchors match the cycle pattern)
            is_cycle_start = True

            # We don't strictly require full cycle match -
            # just that it starts with the right anchor
            # This handles cases with mutations

            cut_positions.append(hits[i].position)
        i += 1

    if not cut_positions:
        return [sequence]

    if verbose:
        print(f"  Cut positions: {cut_positions[:10]}{'...' if len(cut_positions) > 10 else ''}")

    # Cut sequence
    monomers = []

    # Left flank
    if cut_positions[0] > 0:
        monomers.append(sequence[:cut_positions[0]])

    # Monomers
    for i in range(len(cut_positions) - 1):
        monomers.append(sequence[cut_positions[i]:cut_positions[i + 1]])

    # Last part (may be partial)
    if cut_positions[-1] < len(sequence):
        monomers.append(sequence[cut_positions[-1]:])

    return monomers


class AnchorGraphDecomposer:
    """
    Main class for anchor graph-based decomposition.

    Usage:
        decomposer = AnchorGraphDecomposer()
        decomposer.build_from_candidates(sequence, candidates)
        monomers = decomposer.decompose()
    """

    def __init__(self):
        self.sequence: Optional[str] = None
        self.anchors: List[str] = []
        self.hits: List[AnchorHit] = []
        self.graph: Dict = {}
        self.cycle: List[str] = []

    def build_from_candidates(
        self,
        sequence: str,
        candidates: List[Dict],
        top_k: int = 10,
        verbose: bool = False
    ) -> None:
        """
        Build anchor graph from compute_cuts candidates.

        Args:
            sequence: The array sequence
            candidates: List of candidate dicts from compute_cuts
            top_k: Number of top candidates to use
            verbose: Print debug info
        """
        self.sequence = sequence

        # Get top candidates
        top = get_top_candidates(candidates, top_k)
        self.anchors = [c['cut'] for c in top]

        if verbose:
            print(f"Top {len(self.anchors)} anchors:")
            for c in top:
                print(f"  '{c['cut']}' score={c.get('adjusted_score', c.get('base_score', 0)):.3f}")

        # Find hits hierarchically (longer first)
        if verbose:
            print(f"\nFinding anchor hits (longer anchors first):")
        self.hits = find_anchor_hits_hierarchical(sequence, self.anchors, verbose)

        if verbose:
            print(f"\nTotal hits: {len(self.hits)}")

        # Build transition graph
        self.graph = build_transition_graph(self.hits)

        if verbose:
            print(f"\nTransition graph:")
            for edge, count in sorted(self.graph["edges"].items(), key=lambda x: x[1], reverse=True)[:10]:
                dists = self.graph["distances"][edge]
                mean_dist = sum(dists) / len(dists)
                print(f"  {edge[0][:15]}... -> {edge[1][:15]}...: count={count}, mean_dist={mean_dist:.0f}")

        # Find monomer cycle
        if verbose:
            print(f"\nFinding monomer cycle:")
        self.cycle = find_monomer_cycle(self.graph, verbose)

        if verbose and self.cycle:
            length, var = estimate_monomer_length(self.graph, self.cycle)
            print(f"\nEstimated monomer length: {length:.0f} bp (var={var:.0f})")

    def decompose(self, verbose: bool = False) -> List[str]:
        """
        Decompose sequence into monomers.

        Returns:
            List of monomer sequences
        """
        if not self.cycle:
            if verbose:
                print("No cycle found, returning whole sequence")
            return [self.sequence] if self.sequence else []

        if verbose:
            print(f"\nDecomposing by cycle: {' -> '.join(a[:10]+'...' for a in self.cycle)}")

        monomers = decompose_by_cycle(self.sequence, self.hits, self.cycle, verbose)

        # Verify reconstruction
        reconstructed = "".join(monomers)
        if reconstructed != self.sequence:
            if verbose:
                print(f"WARNING: Reconstruction mismatch! {len(self.sequence)} vs {len(reconstructed)}")
        elif verbose:
            print(f"Reconstruction: PERFECT")

        return monomers

    def get_stats(self) -> Dict:
        """Get decomposer statistics."""
        length, var = estimate_monomer_length(self.graph, self.cycle)
        return {
            "num_anchors": len(self.anchors),
            "anchors": self.anchors,
            "num_hits": len(self.hits),
            "cycle": self.cycle,
            "estimated_monomer_length": length,
            "length_variance": var,
            "edge_counts": self.graph.get("edges", {}),
        }
