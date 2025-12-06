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
from statistics import median


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


# =============================================================================
# Distance Clustering Functions
# =============================================================================

def simple_kmeans(data: List[float], k: int, max_iter: int = 100) -> Tuple[List[List[float]], List[float]]:
    """
    Simple k-means clustering without external dependencies.

    Args:
        data: List of values to cluster
        k: Number of clusters
        max_iter: Maximum iterations

    Returns:
        (clusters, centers) - list of clusters and their centers
    """
    if not data or k <= 0:
        return [[]], [0]

    if k == 1:
        return [data], [sum(data) / len(data)]

    # Initialize centers evenly across the data range
    data_sorted = sorted(data)
    n = len(data_sorted)
    centers = [data_sorted[min(i * n // k, n - 1)] for i in range(k)]

    for _ in range(max_iter):
        # Assign points to nearest center
        clusters = [[] for _ in range(k)]
        for x in data:
            nearest = min(range(k), key=lambda i: abs(x - centers[i]))
            clusters[nearest].append(x)

        # Update centers
        new_centers = []
        for i, c in enumerate(clusters):
            if c:
                new_centers.append(sum(c) / len(c))
            else:
                # Empty cluster - keep old center
                new_centers.append(centers[i])

        if new_centers == centers:
            break
        centers = new_centers

    return clusters, centers


def find_optimal_clusters(data: List[float], max_k: int = 4) -> Tuple[int, List[List[float]], List[float]]:
    """
    Find optimal number of clusters (1 to max_k).

    Criteria:
    - All clusters must be significant (>15% of data)
    - Cluster centers must be well separated (gap > 20% of max center)

    Returns:
        (n_clusters, clusters, centers)
    """
    if len(data) < 4:
        return 1, [data], [median(data)]

    best_k = 1
    best_score = 0
    best_clusters = [data]
    best_centers = [median(data)]

    for k in range(2, min(max_k + 1, len(data) // 2 + 1)):
        clusters, centers = simple_kmeans(data, k)

        # Filter out empty clusters
        non_empty = [(c, ctr) for c, ctr in zip(clusters, centers) if c]
        if len(non_empty) < k:
            continue

        clusters = [x[0] for x in non_empty]
        centers = [x[1] for x in non_empty]

        # Check that all clusters are significant (>15% of data)
        min_size = min(len(c) for c in clusters)
        if min_size < len(data) * 0.15:
            continue

        # Check that centers are well separated
        centers_sorted = sorted(centers)
        if len(centers_sorted) < 2:
            continue

        min_gap = min(centers_sorted[i + 1] - centers_sorted[i]
                      for i in range(len(centers_sorted) - 1))
        max_center = max(centers_sorted)

        if min_gap < 0.2 * max_center:
            continue

        # Score: separation quality
        score = min_gap / max_center
        if score > best_score:
            best_score = score
            best_k = k
            best_clusters = clusters
            best_centers = centers

    return best_k, best_clusters, best_centers


def analyze_distance_distribution(distances: List[float], verbose: bool = False) -> Tuple[int, List[float], float]:
    """
    Analyze distance distribution to detect multi-region patterns.

    If anchor appears N times per monomer, distances will have N clusters.
    True period = sum of cluster centers.

    Args:
        distances: List of distances between consecutive anchor occurrences
        verbose: Print debug info

    Returns:
        (n_clusters, cluster_centers, true_period)
    """
    if len(distances) < 4:
        med = median(distances) if distances else 0
        return 1, [med], med

    n_clusters, clusters, centers = find_optimal_clusters(distances, max_k=4)

    # Sort centers by value for consistent ordering
    centers_sorted = sorted(centers)
    true_period = sum(centers_sorted)

    if verbose and n_clusters > 1:
        print(f"  Distance clustering: {n_clusters} clusters")
        for i, (c, ctr) in enumerate(sorted(zip(clusters, centers), key=lambda x: x[1])):
            print(f"    Cluster {i}: center={ctr:.0f}, count={len(c)}")
        print(f"  True period: {true_period:.0f}")

    return n_clusters, centers_sorted, true_period


def classify_positions_by_distance(
    positions: List[int],
    centers: List[float],
    verbose: bool = False
) -> List[int]:
    """
    Classify each position by the cluster of its distance to the next position.

    Args:
        positions: Sorted list of anchor positions
        centers: Cluster centers from analyze_distance_distribution

    Returns:
        List of cluster labels (0, 1, 2, ...) for each position.
        Last position gets label -1 (no next position).
    """
    if len(positions) < 2:
        return [-1] * len(positions)

    labels = []
    for i in range(len(positions) - 1):
        dist = positions[i + 1] - positions[i]
        # Find nearest cluster center
        nearest = min(range(len(centers)), key=lambda j: abs(dist - centers[j]))
        labels.append(nearest)

    # Last position has no "next", mark with -1
    labels.append(-1)

    if verbose:
        label_counts = Counter(l for l in labels if l >= 0)
        print(f"  Position labels: {dict(label_counts)}")

    return labels


def select_cut_positions_by_cluster(
    positions: List[int],
    labels: List[int],
    verbose: bool = False
) -> List[int]:
    """
    Select cut positions based on cluster analysis.

    Strategy: Cut at positions with the same label as the first position.
    This preserves the monomer structure starting from the left flank.

    Args:
        positions: Sorted list of anchor positions
        labels: Cluster labels for each position

    Returns:
        List of positions to cut at
    """
    if not positions or not labels:
        return positions

    # Use the label of the first position as the "start" label
    # This defines which region (A, B, C...) is the monomer start
    start_label = labels[0]

    if start_label < 0:
        # Only one position, return it
        return positions

    # Select all positions with the start label
    cut_positions = [pos for pos, label in zip(positions, labels) if label == start_label]

    # Also include the last position if it wasn't included
    # (to capture the final partial monomer)
    if positions[-1] not in cut_positions:
        cut_positions.append(positions[-1])

    if verbose:
        print(f"  Start label: {start_label}")
        print(f"  Cut positions: {len(cut_positions)} of {len(positions)} total")

    return sorted(cut_positions)


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


def find_monomer_cycle(graph: Dict, hits: List[AnchorHit], sequence: str, verbose: bool = False) -> Tuple[List[str], float, float, int, List[float]]:
    """
    Find the best anchor for decomposition using distance clustering.

    Strategy:
    1. For each anchor with self-loop, analyze distance distribution
    2. Cluster distances to detect multi-region patterns (A-B, A-B-C, etc.)
    3. True period = sum of cluster centers
    4. Select anchor with lowest CV after cluster-based correction

    Returns:
        Tuple of (cycle, mean_period, cv, n_clusters, cluster_centers)
    """
    edges = graph["edges"]
    distances = graph.get("distances", {})

    if not edges:
        return [], 0, float('inf'), 1, []

    # Collect candidate anchors with cluster analysis
    # Format: (cycle, cv, mean_period, count, n_clusters, centers)
    candidates = []

    # Get unique anchors from graph with self-loops
    anchors_to_try = set()
    for (from_a, to_a), count in edges.items():
        if count >= 3:  # Need enough occurrences
            anchors_to_try.add(from_a)
            if from_a != to_a:
                anchors_to_try.add(to_a)

    for anchor in anchors_to_try:
        # Find ALL positions for this anchor (not just self-loop)
        positions = [hit.position for hit in hits if hit.anchor == anchor]
        if len(positions) < 3:
            continue
        positions.sort()

        # Calculate distances between ALL consecutive positions of this anchor
        # This captures the full pattern even if other anchors appear in between
        dist_list = []
        for i in range(len(positions) - 1):
            dist_list.append(positions[i + 1] - positions[i])

        if len(dist_list) < 3:
            continue

        # Analyze distance distribution - detect multi-region patterns
        n_clusters, centers, true_period = analyze_distance_distribution(dist_list, verbose=False)

        # Classify positions and get cut positions
        labels = classify_positions_by_distance(positions, centers)
        cut_positions = select_cut_positions_by_cluster(positions, labels)

        if len(cut_positions) < 2:
            continue

        # Calculate CV of resulting monomer lengths
        lengths = []
        for i in range(len(cut_positions) - 1):
            lengths.append(cut_positions[i + 1] - cut_positions[i])

        if len(lengths) < 2:
            continue

        mean_len = sum(lengths) / len(lengths)
        variance = sum((x - mean_len) ** 2 for x in lengths) / len(lengths)
        cv = (variance ** 0.5) / mean_len if mean_len > 0 else float('inf')

        candidates.append(([anchor], cv, true_period, len(cut_positions), n_clusters, centers))

        if verbose:
            cluster_info = f" ({n_clusters} clusters)" if n_clusters > 1 else ""
            print(f"  Candidate: {anchor[:20]}... period={true_period:.0f}, CV={cv:.3f}, cuts={len(cut_positions)}{cluster_info}")

    # Select best candidate by CV
    if not candidates:
        # Fallback: most frequent edge
        if edges:
            start_edge = max(edges.items(), key=lambda x: x[1])[0]
            if verbose:
                print(f"  Fallback to most frequent: {start_edge[0][:20]}...")
            return [start_edge[0]], 0, float('inf'), 1, []
        return [], 0, float('inf'), 1, []

    # Sort by CV (ascending), then by n_clusters (ascending) - prefer simpler patterns
    candidates.sort(key=lambda x: (x[1], x[4]))

    best_cycle, best_cv, best_period, best_count, best_n_clusters, best_centers = candidates[0]

    if verbose:
        cluster_info = f" ({best_n_clusters} clusters)" if best_n_clusters > 1 else ""
        print(f"  Selected: {best_cycle[0][:20]}... period={best_period:.0f}, CV={best_cv:.3f}{cluster_info}")
        if len(candidates) > 1:
            print(f"  (rejected {len(candidates)-1} other candidates)")

    return best_cycle, best_period, best_cv, best_n_clusters, best_centers


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
    cluster_centers: List[float] = None,
    verbose: bool = False
) -> List[str]:
    """
    Decompose sequence into monomers based on the anchor cycle.

    Strategy:
    1. Find all positions of the anchor
    2. If cluster_centers provided (multi-region pattern), use cluster analysis
       to select only the "start" positions
    3. Cut at selected positions

    Args:
        sequence: The array sequence
        hits: Sorted list of AnchorHit
        cycle: List of anchors representing one monomer
        cluster_centers: If multi-region pattern detected, the cluster centers
                        for classifying positions. None = cut at every position.

    Returns:
        List of monomer sequences
    """
    if not cycle or not hits:
        return [sequence] if sequence else []

    start_anchor = cycle[0]

    # Find all positions where anchor appears
    all_positions = []
    for hit in hits:
        if hit.anchor == start_anchor:
            all_positions.append(hit.position)

    if not all_positions:
        return [sequence]

    all_positions.sort()

    # Determine cut positions based on cluster analysis
    if cluster_centers and len(cluster_centers) > 1:
        # Multi-region pattern: use cluster-based selection
        labels = classify_positions_by_distance(all_positions, cluster_centers)
        cut_positions = select_cut_positions_by_cluster(all_positions, labels, verbose)

        if verbose:
            print(f"  Multi-region pattern: {len(cluster_centers)} clusters")
            print(f"  All positions: {len(all_positions)}, cut positions: {len(cut_positions)}")
    else:
        # Single region: cut at every position
        cut_positions = all_positions

    if verbose:
        print(f"  Cut positions: {cut_positions[:10]}{'...' if len(cut_positions) > 10 else ''}")

    if not cut_positions:
        return [sequence]

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
        self.n_clusters: int = 1  # Number of distance clusters (1 = normal, >1 = multi-region)
        self.cluster_centers: List[float] = []  # Cluster centers for multi-region patterns

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

        # Find monomer cycle with cluster analysis
        if verbose:
            print(f"\nFinding monomer cycle:")
        self.cycle, self.estimated_period, self.estimated_cv, self.n_clusters, self.cluster_centers = find_monomer_cycle(
            self.graph, self.hits, sequence, verbose
        )

        if verbose and self.cycle:
            cluster_info = f", {self.n_clusters} clusters" if self.n_clusters > 1 else ""
            print(f"\nEstimated monomer length: {self.estimated_period:.0f} bp (CV={self.estimated_cv:.3f}{cluster_info})")

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
            cluster_info = f" ({self.n_clusters} clusters)" if self.n_clusters > 1 else ""
            print(f"\nDecomposing by cycle: {' -> '.join(a[:10]+'...' for a in self.cycle)}{cluster_info}")

        monomers = decompose_by_cycle(self.sequence, self.hits, self.cycle, self.cluster_centers, verbose)

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
            "n_clusters": self.n_clusters,
            "cluster_centers": self.cluster_centers,
            "edge_counts": self.graph.get("edges", {}),
        }
