"""Tests for `2-VAE-code/Sensitivity-analysis/topology_compare.py`.

Converted from the `_self_test()` block at the bottom of that module, which only runs
when the script is executed directly and therefore never runs in a test session. The
checks are the same; they are split into named tests so that one failure does not hide
the rest, and the synthetic graphs are built once per session rather than per check.

The module needs igraph, leidenalg and scikit-learn, which are not installed in every
environment used with this repository, so the whole file skips when they are absent.
"""
import numpy as np
import pytest

pytest.importorskip("igraph", reason="topology_compare needs python-igraph")
pytest.importorskip("leidenalg", reason="topology_compare needs leidenalg")
pytest.importorskip("sklearn", reason="topology_compare needs scikit-learn")

import matplotlib.pyplot as plt
import networkx as nx

import topology_compare as tc


@pytest.fixture(scope="session")
def graphs():
    """Two weighted planted-partition graphs over the same four communities."""
    return tc._planted_pair()


@pytest.fixture(scope="session")
def comparison(graphs):
    return tc.compare_sets(*graphs)


@pytest.fixture(scope="session")
def enrichment(graphs, comparison):
    G1, G2 = graphs
    return tc.common_neighbor_enrichment(
        G1, G2, n_baseline=2000, comparison=comparison, directions=("G2_minus_G1",)
    )["G2_minus_G1"]


class TestJaccard:
    def test_partial_overlap(self):
        assert tc.jaccard({1, 2}, {2, 3}) == 1 / 3

    def test_identical_sets(self):
        assert tc.jaccard({1, 2}, {1, 2}) == 1.0

    def test_disjoint_sets(self):
        assert tc.jaccard({1}, {2}) == 0.0

    def test_two_empty_sets_are_undefined(self):
        assert np.isnan(tc.jaccard([], []))


class TestCompareSets:
    def test_the_three_edge_sets_partition_the_union(self, graphs, comparison):
        G1, G2 = graphs
        e1, e2 = tc._edge_set(G1), tc._edge_set(G2)
        only_1 = comparison["edges_only_in_G1"]
        only_2 = comparison["edges_only_in_G2"]
        shared = comparison["shared_edges"]

        assert only_1 | shared == e1
        assert only_2 | shared == e2
        assert not (only_1 & only_2)
        assert not (only_1 & shared)
        assert not (only_2 & shared)

    def test_edge_jaccard_matches_the_edge_sets(self, graphs, comparison):
        G1, G2 = graphs
        e1, e2 = tc._edge_set(G1), tc._edge_set(G2)
        assert comparison["edge_jaccard"] == len(comparison["shared_edges"]) / len(e1 | e2)
        assert 0.0 < comparison["edge_jaccard"] < 1.0

    def test_the_planted_graphs_share_every_node(self, comparison):
        assert comparison["node_jaccard"] == 1.0

    def test_swapping_the_arguments_swaps_the_directions(self, graphs, comparison):
        G1, G2 = graphs
        swapped = tc.compare_sets(G2, G1)
        assert swapped["edges_only_in_G1"] == comparison["edges_only_in_G2"]
        assert swapped["edges_only_in_G2"] == comparison["edges_only_in_G1"]
        assert swapped["edge_jaccard"] == comparison["edge_jaccard"]

    def test_a_graph_compared_with_itself_has_no_difference(self, graphs):
        G1, _ = graphs
        identical = tc.compare_sets(G1, G1)
        assert identical["node_jaccard"] == 1.0
        assert identical["edge_jaccard"] == 1.0
        assert not identical["edges_only_in_G1"]
        assert not identical["edges_only_in_G2"]


class TestSampleNonEdges:
    def test_returns_the_requested_number_of_non_adjacent_pairs(self, graphs):
        G1, _ = graphs
        pairs = tc.sample_non_edges(G1, 500, seed=7)
        assert len(pairs) == 500
        assert all(u != v and not G1.has_edge(u, v) for u, v in pairs)

    def test_the_seed_fixes_the_draw(self, graphs):
        G1, _ = graphs
        assert tc.sample_non_edges(G1, 500, seed=7) == tc.sample_non_edges(G1, 500, seed=7)
        assert tc.sample_non_edges(G1, 500, seed=8) != tc.sample_non_edges(G1, 500, seed=7)

    def test_a_complete_graph_yields_nothing_and_terminates(self):
        # The attempt limit ends the loop instead of blocking.
        assert tc.sample_non_edges(nx.complete_graph(10), 10, seed=0) == []

    def test_a_single_node_graph_yields_nothing(self):
        assert tc.sample_non_edges(nx.empty_graph(1), 10, seed=0) == []


class TestEndpointDistances:
    def test_distances_along_a_path(self):
        path = nx.path_graph(6)  # 0-1-2-3-4-5
        distances, skipped = tc.endpoint_distances(path, [(0, 2), (0, 3), (0, 5)], cutoff=3)
        assert skipped == 0
        assert distances[0] == 2          # closing this pair would make a triangle
        assert distances[1] == 3          # closing this pair would make a square
        assert not np.isfinite(distances[2])  # farther apart than the cutoff

    def test_pairs_with_an_absent_endpoint_are_skipped(self):
        path = nx.path_graph(6)
        distances, skipped = tc.endpoint_distances(path, [(0, 2), (0, 99)], cutoff=3)
        assert distances == [2]
        assert skipped == 1


class TestSummariseDistances:
    def test_fractions_per_distance_class(self):
        summary = tc.summarise_distances([2, 2, 3, np.inf], n_skipped=4)
        assert summary["n"] == 4
        assert summary["n_skipped"] == 4
        assert summary["frac_triangle_d2"] == 0.5
        assert summary["frac_square_d3"] == 0.25
        assert summary["frac_unreachable_or_far"] == 0.25
        assert summary["distance_counts"] == {2: 2, 3: 1}

    def test_empty_input_is_undefined_rather_than_zero(self):
        empty = tc.summarise_distances([])
        assert empty["n"] == 0
        assert np.isnan(empty["frac_triangle_d2"])


class TestCommonNeighbours:
    def test_counts_in_a_four_cycle(self):
        # Each pair of opposite corners shares both remaining nodes.
        square = nx.cycle_graph(4)
        counts, skipped = tc.common_neighbor_counts(square, [(0, 2), (1, 3)])
        assert counts == [2, 2]
        assert skipped == 0

    def test_pairs_with_an_absent_endpoint_are_skipped(self):
        square = nx.cycle_graph(4)
        counts, skipped = tc.common_neighbor_counts(square, [(0, 2), (0, 99)])
        assert counts == [2]
        assert skipped == 1

    def test_identical_samples_show_no_effect(self):
        assert tc.cn_effect_size([1, 2, 3], [1, 2, 3])["prob_superiority"] == 0.5

    def test_a_strictly_larger_sample_is_always_superior(self):
        assert tc.cn_effect_size([5, 6, 7], [1, 2, 3])["prob_superiority"] == 1.0

    def test_proportion_with_at_least_one_common_neighbour(self):
        assert tc.cn_effect_size([0, 0, 1], [1, 1, 1])["prop_CN_ge1_differing"] == 1 / 3

    def test_an_empty_sample_is_undefined(self):
        assert np.isnan(tc.cn_effect_size([], [1, 2])["prob_superiority"])


class TestDirections:
    @pytest.fixture(scope="class")
    def reports(self, graphs, comparison):
        G1, G2 = graphs
        both = tc.indirectness_report(G1, G2, n_baseline=500, comparison=comparison)
        one = tc.indirectness_report(
            G1, G2, n_baseline=500, comparison=comparison, directions=("G2_minus_G1",)
        )
        return both, one

    def test_both_directions_are_reported_by_default(self, reports):
        both, _ = reports
        assert set(both) == set(tc.DIRECTIONS)

    def test_one_direction_can_be_requested(self, reports):
        _, one = reports
        assert set(one) == {"G2_minus_G1"}

    def test_one_direction_gives_the_same_numbers_as_both(self, reports):
        both, one = reports
        assert one["G2_minus_G1"] == both["G2_minus_G1"]

    def test_every_differing_edge_is_covered(self, reports, comparison):
        both, one = reports
        assert one["G2_minus_G1"]["differing_edges"]["n"] == len(
            comparison["edges_only_in_G2"]
        )
        assert both["G1_minus_G2"]["differing_edges"]["n"] == len(
            comparison["edges_only_in_G1"]
        )

    def test_an_unknown_direction_raises(self, graphs):
        G1, G2 = graphs
        with pytest.raises(ValueError):
            tc.indirectness_report(G1, G2, n_baseline=10, directions=("G3_minus_G1",))


class TestIndirectness:
    @pytest.fixture(scope="class")
    def report(self, graphs, comparison):
        G1, G2 = graphs
        return tc.indirectness_report(
            G1, G2, n_baseline=2000, comparison=comparison, directions=("G2_minus_G1",)
        )["G2_minus_G1"]

    def test_differing_edges_close_triangles_more_often_than_random_non_edges(self, report):
        assert (
            report["differing_edges"]["frac_triangle_d2"]
            > report["random_nonedges"]["frac_triangle_d2"]
        )

    @pytest.mark.parametrize("group", ["differing_edges", "random_nonedges"])
    def test_the_three_distance_fractions_account_for_every_pair(self, report, group):
        # Non-adjacent pairs are at distance two or more, so nothing falls outside.
        summary = report[group]
        total = (
            summary["frac_triangle_d2"]
            + summary["frac_square_d3"]
            + summary["frac_unreachable_or_far"]
        )
        assert abs(total - 1.0) < 1e-9


class TestCommonNeighbourEnrichment:
    def test_no_pair_is_skipped_on_graphs_with_the_same_nodes(self, enrichment):
        assert enrichment["n_skipped"] == 0

    def test_differing_edges_have_more_common_neighbours(self, enrichment):
        assert enrichment["prob_superiority"] > 0.5
        assert enrichment["prop_CN_ge1_differing"] > enrichment["prop_CN_ge1_random"]

    def test_the_effect_is_significant(self, enrichment):
        assert enrichment["p_greater"] < 0.01

    def test_every_differing_edge_is_counted(self, enrichment, comparison):
        assert len(enrichment["CN_differing"]) == len(comparison["edges_only_in_G2"])


class TestThresholdProximity:
    @pytest.mark.parametrize(
        "key", ["G1_shared_median", "G1_only_median", "G2_shared_median", "G2_only_median"]
    )
    def test_medians_lie_inside_the_planted_weight_range(self, graphs, key):
        # The fixture draws weights from the interval [0.75, 1.0].
        result = tc.threshold_proximity(*graphs)
        assert 0.75 <= result[key] <= 1.0

    def test_no_weight_is_reported_missing(self, graphs):
        assert "n_edges_without_weight" not in tc.threshold_proximity(*graphs)

    def test_a_graph_without_the_weight_attribute_raises(self):
        unweighted = nx.cycle_graph(4)
        with pytest.raises(KeyError):
            tc.threshold_proximity(unweighted, unweighted)


class TestPartitionStability:
    def test_a_graph_agrees_perfectly_with_itself(self, graphs):
        G1, _ = graphs
        identical = tc.partition_stability(G1, G1)
        assert identical["n_common_nodes"] == G1.number_of_nodes()
        assert identical["ARI"] == 1.0
        assert identical["NMI"] == 1.0

    def test_the_planted_communities_survive_the_edge_differences(self, graphs):
        G1, G2 = graphs
        planted = tc.partition_stability(G1, G2)
        assert planted["n_common_nodes"] == G1.number_of_nodes()
        assert planted["ARI"] > 0.5
        assert planted["NMI"] > 0.5


class TestPlot:
    def test_two_panels_are_drawn(self, enrichment):
        fig = tc.plot_cn_comparison(enrichment["CN_differing"], enrichment["CN_random"])
        assert len(fig.axes) == 2
        plt.close(fig)

    def test_an_empty_count_array_raises(self):
        with pytest.raises(ValueError):
            tc.plot_cn_comparison([], [1, 2])
