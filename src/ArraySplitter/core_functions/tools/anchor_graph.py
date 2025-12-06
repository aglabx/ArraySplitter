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


def find_monomer_cycle(graph: Dict, hits: List[AnchorHit], sequence: str, verbose: bool = False) -> Tuple[List[str], float, float, int]:
    """
    Find the best cycle in the transition graph based on CV of ACTUAL decomposition.

    This represents the sequence of anchors in one monomer.

    Strategy:
    1. Collect candidate anchors (from self-loops and multi-anchor cycles)
    2. For each anchor, try different step sizes (1, 2, 3) to handle A-B patterns
    3. Calculate CV for each (anchor, step) combination
    4. Select the combination with minimum CV (best uniformity)

    The key insight: An anchor may appear multiple times per monomer (A-B pattern).
    Using step=2 means we cut at every 2nd occurrence, effectively merging A+B.

    Returns:
        Tuple of (cycle, mean_period, cv, step)
    """
    edges = graph["edges"]
    distances = graph.get("distances", {})

    if not edges:
        return [], 0, float('inf'), 1

    # Collect all candidate anchors with their actual decomposition CV
    candidates = []  # List of (cycle, cv, mean_length, count)

    # Get unique anchors from graph
    anchors_to_try = set()
    for (from_a, to_a), count in edges.items():
        if count >= 3:  # Need enough occurrences
            anchors_to_try.add(from_a)
            if from_a != to_a:
                anchors_to_try.add(to_a)

    # Evaluate each anchor by ACTUAL decomposition
    for anchor in anchors_to_try:
        # Find all positions of this anchor
        positions = []
        for hit in hits:
            if hit.anchor == anchor:
                positions.append(hit.position)

        if len(positions) < 3:
            continue

        positions.sort()

        # Try different step sizes (1=direct, 2=A-B pattern, etc.)
        # This handles cases where anchor appears multiple times per monomer
        for step in [1, 2, 3]:
            if len(positions) <= step:
                continue

            # Calculate lengths with this step size
            lengths = []
            for i in range(len(positions) - step):
                length = positions[i + step] - positions[i]
                lengths.append(length)

            if len(lengths) < 2:
                continue

            mean_len = sum(lengths) / len(lengths)
            variance = sum((x - mean_len) ** 2 for x in lengths) / len(lengths)
            cv = (variance ** 0.5) / mean_len if mean_len > 0 else float('inf')

            # Store anchor, step, and CV
            candidates.append(([anchor], cv, mean_len, len(positions), step))
            if verbose:
                step_info = f"step={step}" if step > 1 else ""
                print(f"  Candidate: {anchor[:20]}... period={mean_len:.0f}, CV={cv:.3f}, cuts={len(positions)} {step_info}")

    # Select best candidate by CV (lower is better)
    if not candidates:
        # Last resort: most frequent edge
        if edges:
            start_edge = max(edges.items(), key=lambda x: x[1])[0]
            if verbose:
                print(f"  Fallback to most frequent: {start_edge[0][:20]}...")
            return [start_edge[0]], 0, float('inf'), 1
        return [], 0, float('inf'), 1

    # Sort by CV (ascending) - best uniformity first
    candidates.sort(key=lambda x: x[1])

    best_cycle, best_cv, best_mean, best_count, best_step = candidates[0]

    if verbose:
        step_info = f" step={best_step}" if best_step > 1 else ""
        print(f"  Selected: {best_cycle[0][:20]}... period={best_mean:.0f}, CV={best_cv:.3f}{step_info}")
        if len(candidates) > 1:
            print(f"  (rejected {len(candidates)-1} other candidates with higher CV)")

    return best_cycle, best_mean, best_cv, best_step


def _calculate_cycle_distances(cycle: List[str], distances: Dict) -> List[float]:
    """Calculate full cycle distances from edge distances."""
    if len(cycle) == 1:
        edge = (cycle[0], cycle[0])
        return distances.get(edge, [])

    # For multi-anchor cycle, sum distances along the cycle
    # This is approximate - we'd need to track actual cycle instances
    cycle_dists = []
    edge = (cycle[0], cycle[1])
    if edge in distances:
        # Use first edge distances as proxy
        # TODO: properly track full cycle instances
        cycle_dists = distances[edge]

    return cycle_dists


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
    step: int = 1,
    verbose: bool = False
) -> List[str]:
    """
    Decompose sequence into monomers based on the anchor cycle.

    Strategy:
    1. Find positions where the cycle starts (first anchor of cycle)
    2. If step > 1, only cut at every step-th position (A-B pattern handling)
    3. Cut at those positions

    Args:
        sequence: The array sequence
        hits: Sorted list of AnchorHit
        cycle: List of anchors representing one monomer
        step: Cut at every step-th occurrence (1=every, 2=every 2nd for A-B pattern)

    Returns:
        List of monomer sequences
    """
    if not cycle or not hits:
        return [sequence] if sequence else []

    start_anchor = cycle[0]

    # Find all positions where anchor appears
    all_positions = []

    i = 0
    while i < len(hits):
        if hits[i].anchor == start_anchor:
            all_positions.append(hits[i].position)
        i += 1

    if not all_positions:
        return [sequence]

    # Apply step: only use every step-th position
    cut_positions = all_positions[::step]

    if verbose:
        if step > 1:
            print(f"  All anchor positions: {len(all_positions)}, using every {step}th: {len(cut_positions)}")
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
        self.estimated_period: float = 0
        self.estimated_cv: float = float('inf')
        self.step: int = 1  # Step size for A-B pattern handling

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

        # Get top candidates, filtering by max anchor length
        # Anchors should be short (<=11bp) - they are graph nodes, not sequences
        MAX_ANCHOR_LENGTH = 11
        top = get_top_candidates(candidates, top_k)
        self.anchors = [c['cut'] for c in top if len(c['cut']) <= MAX_ANCHOR_LENGTH]

        if verbose:
            print(f"Top {len(self.anchors)} anchors (max length {MAX_ANCHOR_LENGTH}bp):")
            for c in top:
                if len(c['cut']) <= MAX_ANCHOR_LENGTH:
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
        self.cycle, self.estimated_period, self.estimated_cv, self.step = find_monomer_cycle(
            self.graph, self.hits, sequence, verbose
        )

        if verbose and self.cycle:
            step_info = f", step={self.step}" if self.step > 1 else ""
            print(f"\nEstimated monomer length: {self.estimated_period:.0f} bp (CV={self.estimated_cv:.3f}{step_info})")

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
            step_info = f" (step={self.step})" if self.step > 1 else ""
            print(f"\nDecomposing by cycle: {' -> '.join(a[:10]+'...' for a in self.cycle)}{step_info}")

        monomers = decompose_by_cycle(self.sequence, self.hits, self.cycle, self.step, verbose)

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
        return {
            "num_anchors": len(self.anchors),
            "anchors": self.anchors,
            "num_hits": len(self.hits),
            "cycle": self.cycle,
            "estimated_monomer_length": self.estimated_period,
            "estimated_cv": self.estimated_cv,
            "edge_counts": self.graph.get("edges", {}),
        }
