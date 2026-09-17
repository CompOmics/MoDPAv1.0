"""Compare two thresholded PTM association networks.

The two networks share most of their nodes but differ in their edges. The
functions here test whether the differing edges are explained by indirect
(triangle- or square-closing) structure rather than by genuine rewiring.

Networks are undirected `networkx.Graph` objects. Edge weights, the signed
distance correlation (SDCor), are read from the attribute named by
`weight_attr`, which defaults to `Score`.

Every comparison is expressed in two directions:

* `G1_minus_G2`: edges present in `G1` but absent from `G2`, evaluated inside `G2`;
* `G2_minus_G1`: edges present in `G2` but absent from `G1`, evaluated inside `G1`.

Pass `directions=("G2_minus_G1",)` to compute only one of them, which halves the
runtime of the two expensive functions.

Typical use
-----------
    import networkx as nx
    from topology_compare import (
        compare_sets, indirectness_report, common_neighbor_enrichment,
        threshold_proximity, partition_stability,
    )

    G1 = nx.from_pandas_edgelist(edges_a, "nodeA", "nodeB", edge_attr="Score")
    G2 = nx.from_pandas_edgelist(edges_b, "nodeA", "nodeB", edge_attr="Score")

    # Reuse one comparison rather than recomputing the edge sets each time.
    comparison = compare_sets(G1, G2)
    print(comparison["node_jaccard"], comparison["edge_jaccard"])

    # 1. Are the differing edges close (triangle or square) in the network that
    #    lacks them, relative to random non-edges?
    report = indirectness_report(G1, G2, comparison=comparison,
                                 directions=("G2_minus_G1",))

    # 2. Do the differing edges have many common neighbours there?
    enrichment = common_neighbor_enrichment(G1, G2, comparison=comparison,
                                            directions=("G2_minus_G1",))

    # 3. Do the differing edges sit close to the score cutoff?
    print(threshold_proximity(G1, G2, comparison=comparison))

    # 4. Does the community structure survive?
    print(partition_stability(G1, G2, resolution=0.5))
"""

from __future__ import annotations

from collections import Counter
import random
from typing import Iterable, Sequence

import igraph as ig
import leidenalg as la
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.figure import Figure
from scipy.stats import mannwhitneyu
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
)

__all__ = [
    "DEFAULT_WEIGHT_ATTR",
    "DIRECTIONS",
    "jaccard",
    "compare_sets",
    "sample_non_edges",
    "endpoint_distances",
    "summarise_distances",
    "indirectness_report",
    "common_neighbor_counts",
    "common_neighbor_enrichment",
    "cn_effect_size",
    "plot_cn_comparison",
    "threshold_proximity",
    "partition_stability",
]

DEFAULT_WEIGHT_ATTR = "Score"
DIRECTIONS = ("G1_minus_G2", "G2_minus_G1")

Pair = frozenset


# ---------------------------------------------------------------------------
# Basic set comparisons
# ---------------------------------------------------------------------------

def _edge_set(G: nx.Graph) -> set[Pair]:
    return {frozenset((u, v)) for u, v in G.edges()}


def jaccard(a: Iterable, b: Iterable) -> float:
    """Jaccard index of two sets. Returns NaN when both are empty."""
    a, b = set(a), set(b)
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def compare_sets(G1: nx.Graph, G2: nx.Graph) -> dict:
    """Node and edge Jaccard indices plus the symmetric-difference edge sets.

    This materialises the full edge sets, so it costs O(|E1| + |E2|) time and
    memory. Compute it once and pass the result to the other functions through
    their `comparison` argument.
    """
    n1, n2 = set(G1.nodes()), set(G2.nodes())
    e1, e2 = _edge_set(G1), _edge_set(G2)
    return {
        "node_jaccard": jaccard(n1, n2),
        "edge_jaccard": jaccard(e1, e2),
        "edges_only_in_G1": e1 - e2,
        "edges_only_in_G2": e2 - e1,
        "shared_edges": e1 & e2,
    }


def _resolve(G1, G2, comparison, directions):
    """Return the comparison dict and the (label, G_absent, missing) triples."""
    if comparison is None:
        comparison = compare_sets(G1, G2)
    unknown = set(directions) - set(DIRECTIONS)
    if unknown:
        raise ValueError(f"unknown direction(s) {sorted(unknown)}; expected {DIRECTIONS}")
    available = {
        "G1_minus_G2": (G2, comparison["edges_only_in_G1"]),
        "G2_minus_G1": (G1, comparison["edges_only_in_G2"]),
    }
    return comparison, [(label, *available[label]) for label in directions]


# ---------------------------------------------------------------------------
# Random non-edge baseline
# ---------------------------------------------------------------------------

def sample_non_edges(G: nx.Graph, n: int, seed: int = 0) -> list[tuple]:
    """Sample up to `n` distinct random node pairs that are not adjacent in `G`.

    Sampling stops after `100 * n` attempts, so fewer than `n` pairs may be
    returned for a very dense graph. Pairs are not deduplicated, matching the
    behaviour of an independent draw with replacement.
    """
    rng = random.Random(seed)
    nodes = list(G.nodes())
    if len(nodes) < 2:
        return []
    pairs: list[tuple] = []
    tries, limit = 0, max(n * 100, 1000)
    while len(pairs) < n and tries < limit:
        tries += 1
        u, v = rng.sample(nodes, 2)
        if G.has_edge(u, v):
            continue
        pairs.append((u, v))
    return pairs


# ---------------------------------------------------------------------------
# 1. Are the differing edges indirect in the network that lacks them?
#    distance 2 -> the edge would close a triangle
#    distance 3 -> the edge would close a square (4-cycle)
#    In a dense graph most pairs are already at distance 2, so only the
#    enrichment over the random non-edge baseline is meaningful.
# ---------------------------------------------------------------------------

def endpoint_distances(
    G_absent: nx.Graph,
    pairs: Iterable,
    cutoff: int = 3,
) -> tuple[list[float], int]:
    """Shortest-path distance between the endpoints of each pair, inside `G_absent`.

    `np.inf` means the endpoints are farther apart than `cutoff` or lie in
    different components. Pairs with an endpoint missing from `G_absent` are
    skipped, since the two networks do not have identical node sets.

    Breadth-first searches are grouped by source node, so the cost is one BFS
    per distinct source rather than one per pair.

    Returns the list of distances and the number of skipped pairs.
    """
    by_source: dict = {}
    skipped = 0
    for pair in pairs:
        u, v = tuple(pair)
        if u in G_absent and v in G_absent:
            by_source.setdefault(u, []).append(v)
        else:
            skipped += 1

    distances: list[float] = []
    for u, targets in by_source.items():
        reachable = nx.single_source_shortest_path_length(G_absent, u, cutoff=cutoff)
        distances.extend(reachable.get(v, np.inf) for v in targets)
    return distances, skipped


def summarise_distances(distances: Sequence[float], n_skipped: int = 0) -> dict:
    """Summarise a list of endpoint distances as fractions per distance class."""
    counts = Counter(d for d in distances if np.isfinite(d))
    n = len(distances)
    nan = float("nan")
    return {
        "n": n,
        "n_skipped": n_skipped,
        "frac_triangle_d2": counts.get(2, 0) / n if n else nan,
        "frac_square_d3": counts.get(3, 0) / n if n else nan,
        "frac_unreachable_or_far": (
            sum(1 for d in distances if not np.isfinite(d)) / n if n else nan
        ),
        "distance_counts": dict(sorted(counts.items())),
    }


def indirectness_report(
    G1: nx.Graph,
    G2: nx.Graph,
    n_baseline: int = 100_000,
    cutoff: int = 3,
    seed: int = 0,
    comparison: dict | None = None,
    directions: Sequence[str] = DIRECTIONS,
) -> dict:
    """Endpoint-distance profile of the differing edges versus random non-edges.

    Read it as: are the differing edges enriched at distance 2 (triangle) or 3
    (square) relative to random non-edges in the network that lacks them?
    """
    _, targets = _resolve(G1, G2, comparison, directions)
    report = {}
    for label, G_absent, missing in targets:
        differing, skipped = endpoint_distances(G_absent, missing, cutoff=cutoff)
        baseline, _ = endpoint_distances(
            G_absent, sample_non_edges(G_absent, n_baseline, seed=seed), cutoff=cutoff
        )
        report[label] = {
            "differing_edges": summarise_distances(differing, skipped),
            "random_nonedges": summarise_distances(baseline),
        }
    return report


# ---------------------------------------------------------------------------
# 2. Do the differing edges have many common neighbours in the network that
#    lacks them? This is the number of triangles each edge would close.
# ---------------------------------------------------------------------------

def common_neighbor_counts(G: nx.Graph, pairs: Iterable) -> tuple[list[int], int]:
    """Number of common neighbours of each pair in `G`, plus the number skipped."""
    counts, skipped = [], 0
    for pair in pairs:
        u, v = tuple(pair)
        if u in G and v in G:
            counts.append(sum(1 for _ in nx.common_neighbors(G, u, v)))
        else:
            skipped += 1
    return counts, skipped


def cn_effect_size(cn_differing: Sequence[int], cn_random: Sequence[int]) -> dict:
    """One-sided Mann-Whitney test and common-language effect size.

    `prob_superiority` is the probability that a randomly chosen differing edge
    has more common neighbours than a randomly chosen non-edge, with ties
    counted as one half. A value of 0.5 means no effect. SciPy assigns average
    ranks to ties, so U / (n1 * n2) is the tie-corrected form of this quantity.
    """
    differing = np.asarray(cn_differing)
    baseline = np.asarray(cn_random)
    if differing.size == 0 or baseline.size == 0:
        nan = float("nan")
        return {
            "prop_CN_ge1_differing": nan,
            "prop_CN_ge1_random": nan,
            "prob_superiority": nan,
            "mannwhitney_U": nan,
            "p_greater": nan,
        }
    U, p = mannwhitneyu(differing, baseline, alternative="greater")
    return {
        "prop_CN_ge1_differing": float(np.mean(differing >= 1)),
        "prop_CN_ge1_random": float(np.mean(baseline >= 1)),
        "prob_superiority": float(U / (differing.size * baseline.size)),
        "mannwhitney_U": float(U),
        "p_greater": float(p),
    }


def common_neighbor_enrichment(
    G1: nx.Graph,
    G2: nx.Graph,
    n_baseline: int = 100_000,
    seed: int = 0,
    comparison: dict | None = None,
    directions: Sequence[str] = DIRECTIONS,
) -> dict:
    """Do the differing edges have more common neighbours than random non-edges?

    Each direction draws its own baseline from a generator seeded with `seed`,
    so the result for one direction does not depend on whether the other was
    also requested.
    """
    _, targets = _resolve(G1, G2, comparison, directions)
    nan = float("nan")
    report = {}
    for label, G_absent, missing in targets:
        differing, skipped = common_neighbor_counts(G_absent, missing)
        baseline, _ = common_neighbor_counts(
            G_absent, sample_non_edges(G_absent, n_baseline, seed=seed)
        )
        report[label] = {
            "median_CN_differing": float(np.median(differing)) if differing else nan,
            "median_CN_random": float(np.median(baseline)) if baseline else nan,
            "CN_differing": differing,
            "CN_random": baseline,
            "n_skipped": skipped,
            **cn_effect_size(differing, baseline),
        }
    return report


def plot_cn_comparison(
    cn_differing: Sequence[int],
    cn_random: Sequence[int],
    max_bin: int = 5,
    k_max_ccdf: int = 12,
) -> Figure:
    """Compare common-neighbour counts of differing edges and random non-edges.

    The counts are zero-inflated integers, so the median is uninformative.
    Left panel: proportion of pairs at each common-neighbour count, with
    everything at or above `max_bin` collapsed into one category.
    Right panel: complementary CDF, the proportion of pairs with at least k
    common neighbours.
    """
    differing = np.asarray(cn_differing)
    baseline = np.asarray(cn_random)
    if differing.size == 0 or baseline.size == 0:
        raise ValueError("both common-neighbour count arrays must be non-empty")

    def proportions(x: np.ndarray) -> np.ndarray:
        clipped = np.minimum(x, max_bin).astype(int)
        counts = np.bincount(clipped, minlength=max_bin + 1)
        return counts / counts.sum()

    def survival(x: np.ndarray, ks: np.ndarray) -> list[float]:
        return [float(np.mean(x >= k)) for k in ks]

    categories = [str(k) for k in range(max_bin)] + [f"{max_bin}+"]
    positions = np.arange(len(categories))
    width = 0.4

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    ax1.bar(positions - width / 2, proportions(differing), width,
            label="differing edges", color="#3b6ea5")
    ax1.bar(positions + width / 2, proportions(baseline), width,
            label="random non-edges", color="#bdbdbd")
    ax1.set_xticks(positions)
    ax1.set_xticklabels(categories)
    ax1.set_xlabel("common neighbours (triangles closed)")
    ax1.set_ylabel("proportion of pairs")
    ax1.legend(frameon=False)

    ks = np.arange(0, k_max_ccdf + 1)
    ax2.plot(ks, survival(differing, ks), marker="o",
             label="differing edges", color="#3b6ea5")
    ax2.plot(ks, survival(baseline, ks), marker="o",
             label="random non-edges", color="#999999")
    ax2.set_xlabel("k")
    ax2.set_ylabel("proportion with CN \u2265 k")
    ax2.legend(frameon=False)
    # A log y-axis sharpens the tail but hides exact zeros:
    # ax2.set_yscale("log")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 3. Do the differing edges simply sit near the score cutoff?
# ---------------------------------------------------------------------------

def _abs_edge_weights(
    G: nx.Graph, edges: Iterable, attr: str
) -> tuple[list[float], int]:
    """Absolute edge weights of `edges` in `G`, plus the number without `attr`."""
    weights, missing = [], 0
    for edge in edges:
        u, v = tuple(edge)
        if not G.has_edge(u, v):
            continue
        value = G[u][v].get(attr)
        if value is None or not np.isfinite(value):
            missing += 1
        else:
            weights.append(abs(value))
    return weights, missing


def threshold_proximity(
    G1: nx.Graph,
    G2: nx.Graph,
    weight_attr: str = DEFAULT_WEIGHT_ATTR,
    comparison: dict | None = None,
) -> dict:
    """Median absolute weight of the differing edges versus the shared edges.

    Differing edges sitting closer to the cutoff than shared edges indicate that
    the difference is threshold flipping rather than structural change.

    Raises `KeyError` if no edge carries `weight_attr`, which would otherwise
    return NaN for every field.
    """
    if comparison is None:
        comparison = compare_sets(G1, G2)

    result, total_missing, total_found = {}, 0, 0
    for name, G, edges in [
        ("G1_shared", G1, comparison["shared_edges"]),
        ("G1_only", G1, comparison["edges_only_in_G1"]),
        ("G2_shared", G2, comparison["shared_edges"]),
        ("G2_only", G2, comparison["edges_only_in_G2"]),
    ]:
        weights, missing = _abs_edge_weights(G, edges, weight_attr)
        total_missing += missing
        total_found += len(weights)
        result[f"{name}_median"] = float(np.median(weights)) if weights else float("nan")

    if total_found == 0:
        raise KeyError(
            f"no edge carries the attribute {weight_attr!r}; "
            "check the value passed as weight_attr"
        )
    if total_missing:
        result["n_edges_without_weight"] = total_missing
    return result


# ---------------------------------------------------------------------------
# 4. Does the community structure survive, even where the edges do not?
# ---------------------------------------------------------------------------

def partition_stability(
    G1: nx.Graph,
    G2: nx.Graph,
    resolution: float = 0.5,
    seed: int = 0,
    n_iterations: int = 2,
    weight_attr: str = DEFAULT_WEIGHT_ATTR,
) -> dict:
    """Agreement between the Leiden partitions of two networks.

    Each network is restricted to the shared node set, partitioned, and the two
    partitions are compared. `RBConfigurationVertexPartition` is modularity with
    a resolution parameter. Weights are taken as absolute values, because
    modularity is not defined for signed weights.

    Note that this restricts the graphs before clustering, whereas clustering
    each full network and then aligning the labels on the shared nodes is a
    different measurement. The two are not expected to agree exactly, and
    neither is the Leiden implementation used here identical to Cytoscape's
    Leiden app at the same resolution value. Treat the result as a relative
    stability check.
    """
    common = sorted(set(G1.nodes()) & set(G2.nodes()))
    index = {node: i for i, node in enumerate(common)}

    def leiden_labels(G: nx.Graph) -> np.ndarray:
        edges, weights = [], []
        for u, v, data in G.edges(data=True):
            if u in index and v in index:
                edges.append((index[u], index[v]))
                weights.append(abs(data.get(weight_attr, 1.0)))
        h = ig.Graph(n=len(common), edges=edges)
        partition = la.find_partition(
            h,
            la.RBConfigurationVertexPartition,
            weights=weights,
            resolution_parameter=resolution,
            n_iterations=n_iterations,
            seed=seed,
        )
        labels = np.empty(len(common), dtype=int)
        for cluster_id, members in enumerate(partition):
            for node in members:
                labels[node] = cluster_id
        return labels

    labels1, labels2 = leiden_labels(G1), leiden_labels(G2)
    return {
        "n_common_nodes": len(common),
        "ARI": adjusted_rand_score(labels1, labels2),
        "NMI": normalized_mutual_info_score(labels1, labels2),
        "AMI": adjusted_mutual_info_score(labels1, labels2),
    }


# ---------------------------------------------------------------------------
# Self test
#
# Every check below asserts the property it is meant to hold, so a change in
# behaviour fails here instead of passing unnoticed. Run with
# `python topology_compare.py`; no model output is needed. Note that `python -O`
# strips assertions and therefore disables the whole test.
# ---------------------------------------------------------------------------

def _planted_pair() -> tuple[nx.Graph, nx.Graph]:
    """Two weighted planted-partition graphs over the same four communities."""
    rng = np.random.default_rng(0)
    G1 = nx.planted_partition_graph(4, 60, 0.20, 0.02, seed=1)
    G2 = nx.planted_partition_graph(4, 60, 0.18, 0.03, seed=2)
    for G in (G1, G2):
        for u, v in G.edges():
            G[u][v][DEFAULT_WEIGHT_ATTR] = float(rng.uniform(0.75, 1.0))
    return G1, G2


def _assert_raises(exception: type[Exception], message: str, function, *args, **kwargs):
    """Assert that calling `function` raises `exception`."""
    try:
        function(*args, **kwargs)
    except exception:
        return
    raise AssertionError(message)


def _test_jaccard() -> None:
    assert jaccard({1, 2}, {2, 3}) == 1 / 3
    assert jaccard({1, 2}, {1, 2}) == 1.0
    assert jaccard({1}, {2}) == 0.0
    assert np.isnan(jaccard([], []))


def _test_compare_sets(G1: nx.Graph, G2: nx.Graph) -> None:
    comparison = compare_sets(G1, G2)
    e1, e2 = _edge_set(G1), _edge_set(G2)
    only_1 = comparison["edges_only_in_G1"]
    only_2 = comparison["edges_only_in_G2"]
    shared = comparison["shared_edges"]

    # The three sets partition the union of the two edge sets.
    assert only_1 | shared == e1
    assert only_2 | shared == e2
    assert not (only_1 & only_2) and not (only_1 & shared) and not (only_2 & shared)
    assert comparison["edge_jaccard"] == len(shared) / len(e1 | e2)
    assert 0.0 < comparison["edge_jaccard"] < 1.0
    assert comparison["node_jaccard"] == 1.0  # the planted graphs share every node

    # Swapping the arguments swaps the two directions and leaves the indices.
    swapped = compare_sets(G2, G1)
    assert swapped["edges_only_in_G1"] == only_2
    assert swapped["edges_only_in_G2"] == only_1
    assert swapped["edge_jaccard"] == comparison["edge_jaccard"]

    identical = compare_sets(G1, G1)
    assert identical["node_jaccard"] == 1.0 and identical["edge_jaccard"] == 1.0
    assert not identical["edges_only_in_G1"] and not identical["edges_only_in_G2"]


def _test_sample_non_edges(G1: nx.Graph) -> None:
    pairs = sample_non_edges(G1, 500, seed=7)
    assert len(pairs) == 500
    assert all(u != v and not G1.has_edge(u, v) for u, v in pairs)

    # The seed fixes the draw.
    assert sample_non_edges(G1, 500, seed=7) == pairs
    assert sample_non_edges(G1, 500, seed=8) != pairs

    # A complete graph has no non-adjacent pair, and the attempt limit ends the
    # loop rather than blocking.
    assert sample_non_edges(nx.complete_graph(10), 10, seed=0) == []
    assert sample_non_edges(nx.empty_graph(1), 10, seed=0) == []


def _test_endpoint_distances() -> None:
    path = nx.path_graph(6)  # 0-1-2-3-4-5

    distances, skipped = endpoint_distances(path, [(0, 2), (0, 3), (0, 5)], cutoff=3)
    assert skipped == 0
    assert distances[0] == 2 and distances[1] == 3  # triangle and square closure
    assert not np.isfinite(distances[2])            # farther than the cutoff

    # Pairs with an endpoint outside the graph are skipped, not counted.
    distances, skipped = endpoint_distances(path, [(0, 2), (0, 99)], cutoff=3)
    assert distances == [2] and skipped == 1


def _test_summarise_distances() -> None:
    summary = summarise_distances([2, 2, 3, np.inf], n_skipped=4)
    assert summary["n"] == 4 and summary["n_skipped"] == 4
    assert summary["frac_triangle_d2"] == 0.5
    assert summary["frac_square_d3"] == 0.25
    assert summary["frac_unreachable_or_far"] == 0.25
    assert summary["distance_counts"] == {2: 2, 3: 1}

    empty = summarise_distances([])
    assert empty["n"] == 0 and np.isnan(empty["frac_triangle_d2"])


def _test_common_neighbors() -> None:
    # In a 4-cycle the two opposite pairs each share both remaining nodes.
    square = nx.cycle_graph(4)
    counts, skipped = common_neighbor_counts(square, [(0, 2), (1, 3)])
    assert counts == [2, 2] and skipped == 0
    counts, skipped = common_neighbor_counts(square, [(0, 2), (0, 99)])
    assert counts == [2] and skipped == 1

    # Two identical samples give no effect; a strictly larger sample gives a
    # probability of superiority of one.
    assert cn_effect_size([1, 2, 3], [1, 2, 3])["prob_superiority"] == 0.5
    assert cn_effect_size([5, 6, 7], [1, 2, 3])["prob_superiority"] == 1.0
    assert cn_effect_size([0, 0, 1], [1, 1, 1])["prop_CN_ge1_differing"] == 1 / 3
    assert np.isnan(cn_effect_size([], [1, 2])["prob_superiority"])


def _test_directions(G1: nx.Graph, G2: nx.Graph) -> None:
    comparison = compare_sets(G1, G2)
    both = indirectness_report(G1, G2, n_baseline=500, comparison=comparison)
    one = indirectness_report(
        G1, G2, n_baseline=500, comparison=comparison, directions=("G2_minus_G1",)
    )
    assert set(both) == set(DIRECTIONS)
    assert set(one) == {"G2_minus_G1"}

    # Requesting one direction gives exactly the numbers of that direction in
    # the two-direction report.
    assert one["G2_minus_G1"] == both["G2_minus_G1"]

    # G2_minus_G1 evaluates the edges of G2 that are absent from G1, inside G1,
    # so it covers every edge in edges_only_in_G2.
    assert (
        one["G2_minus_G1"]["differing_edges"]["n"]
        == len(comparison["edges_only_in_G2"])
    )
    assert both["G1_minus_G2"]["differing_edges"]["n"] == len(
        comparison["edges_only_in_G1"]
    )

    _assert_raises(
        ValueError,
        "an unknown direction must raise ValueError",
        indirectness_report,
        G1,
        G2,
        n_baseline=10,
        directions=("G3_minus_G1",),
    )


def _test_threshold_proximity(G1: nx.Graph, G2: nx.Graph) -> None:
    result = threshold_proximity(G1, G2)
    for key in ("G1_shared_median", "G1_only_median", "G2_shared_median", "G2_only_median"):
        # The weights were drawn from the interval [0.75, 1.0].
        assert 0.75 <= result[key] <= 1.0
    assert "n_edges_without_weight" not in result

    # An unweighted graph has no edge carrying the attribute at all.
    unweighted = nx.cycle_graph(4)
    _assert_raises(
        KeyError,
        "a missing weight attribute must raise KeyError",
        threshold_proximity,
        unweighted,
        unweighted,
    )


def _test_partition_stability(G1: nx.Graph, G2: nx.Graph) -> None:
    identical = partition_stability(G1, G1)
    assert identical["n_common_nodes"] == G1.number_of_nodes()
    assert identical["ARI"] == 1.0 and identical["NMI"] == 1.0

    # Both graphs carry the same four planted communities, so their partitions
    # must agree far better than chance.
    planted = partition_stability(G1, G2)
    assert planted["n_common_nodes"] == G1.number_of_nodes()
    assert planted["ARI"] > 0.5 and planted["NMI"] > 0.5


def _test_plot(enrichment: dict) -> None:
    fig = plot_cn_comparison(enrichment["CN_differing"], enrichment["CN_random"])
    assert len(fig.axes) == 2
    plt.close(fig)

    _assert_raises(
        ValueError,
        "an empty count array must raise ValueError",
        plot_cn_comparison,
        [],
        [1, 2],
    )


def _self_test() -> None:
    """End-to-end check on two synthetic planted-partition graphs."""
    G1, G2 = _planted_pair()
    comparison = compare_sets(G1, G2)

    _test_jaccard()
    _test_compare_sets(G1, G2)
    _test_sample_non_edges(G1)
    _test_endpoint_distances()
    _test_summarise_distances()
    _test_common_neighbors()
    _test_directions(G1, G2)
    _test_threshold_proximity(G1, G2)
    _test_partition_stability(G1, G2)

    report = indirectness_report(
        G1, G2, n_baseline=2000, comparison=comparison, directions=("G2_minus_G1",)
    )["G2_minus_G1"]
    differing, random_pairs = report["differing_edges"], report["random_nonedges"]

    # An edge of G2 that is absent from G1 closes a triangle in G1 more often
    # than a random non-adjacent pair does.
    assert differing["frac_triangle_d2"] > random_pairs["frac_triangle_d2"]

    # Non-adjacent pairs are at distance two or more, so the three fractions
    # account for every pair.
    for summary in (differing, random_pairs):
        total = (
            summary["frac_triangle_d2"]
            + summary["frac_square_d3"]
            + summary["frac_unreachable_or_far"]
        )
        assert abs(total - 1.0) < 1e-9

    enrichment = common_neighbor_enrichment(
        G1, G2, n_baseline=2000, comparison=comparison, directions=("G2_minus_G1",)
    )["G2_minus_G1"]

    # The same enrichment measured through common neighbours.
    assert enrichment["n_skipped"] == 0
    assert enrichment["prob_superiority"] > 0.5
    assert enrichment["p_greater"] < 0.01
    assert enrichment["prop_CN_ge1_differing"] > enrichment["prop_CN_ge1_random"]
    assert len(enrichment["CN_differing"]) == len(comparison["edges_only_in_G2"])

    _test_plot(enrichment)

    print(f"edge jaccard          {comparison['edge_jaccard']:.3f}")
    print(f"threshold proximity   {threshold_proximity(G1, G2, comparison=comparison)}")
    print(f"differing edges       {differing}")
    print(f"random non-edges      {random_pairs}")
    print(f"prob superiority      {enrichment['prob_superiority']:.3f}")
    print(f"partition stability   {partition_stability(G1, G2)}")
    print("self test passed")


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    _self_test()
