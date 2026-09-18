"""Tests for `parse_network` in `modpa_network_utilities.py`.

`parse_network` decides which associations become edges of the MoDPA network, so every
reported edge and node count rests on it. The filtering rules and the artefact
annotation are covered here against small hand-written association lists.

The first test is a regression guard. An earlier revision of the library ended the
query with `.filter(~pl.col("potential_artefact"))`, which silently removed part of the
network and made the published edge count impossible to reproduce. The README states
that these edges are annotated, not removed, and that is what is asserted.

The module imports igraph, leidenalg and networkx at module level, which are not
installed in every environment used with this repository, so the whole file skips when
they are absent.
"""
import pytest

pytest.importorskip("igraph", reason="modpa_network_utilities needs python-igraph")
pytest.importorskip("leidenalg", reason="modpa_network_utilities needs leidenalg")
pytest.importorskip("networkx", reason="modpa_network_utilities needs networkx")

import polars as pl

from modpa_network_utilities import count_edges, count_nodes, parse_network

COLUMNS = ["nodeA", "nodeB", "Score", "pval", "distance", "PCC", "qvalue"]


def write_associations(tmp_path, rows, name="associations.csv"):
    """Write an association list and return its path.

    `rows` holds one `(nodeA, nodeB, Score, qvalue)` tuple per association. The columns
    `parse_network` does not read are filled in so that the file has the same shape as
    the output of `calculate_sdcorr.py`.
    """
    path = tmp_path / name
    pl.DataFrame(
        [
            {
                "nodeA": node_a,
                "nodeB": node_b,
                "Score": score,
                "pval": 1e-8,
                "distance": abs(score),
                "PCC": score,
                "qvalue": qvalue,
            }
            for node_a, node_b, score, qvalue in rows
        ],
        schema={
            "nodeA": pl.String,
            "nodeB": pl.String,
            "Score": pl.Float64,
            "pval": pl.Float64,
            "distance": pl.Float64,
            "PCC": pl.Float64,
            "qvalue": pl.Float64,
        },
    ).write_csv(path)
    return path


def parse(tmp_path, rows, **kwargs):
    """Parse an association list built from `rows` and collect the result."""
    kwargs.setdefault("min_score", 0.6)
    path = write_associations(tmp_path, rows)
    return parse_network(path, **kwargs).collect()


def test_potential_artefact_edges_are_annotated_not_removed(tmp_path):
    # Two phosphosites three residues apart on one protein: same modification, inside
    # the default adjacency window, so the pair is flagged. It must still be an edge.
    edges = parse(
        tmp_path,
        [("P00001|10|S|21", "P00001|13|T|21", 0.9, 0.001)],
    )

    assert edges.height == 1
    assert edges["potential_artefact"].to_list() == [True]
    assert count_edges(edges) == 1
    assert count_nodes(edges) == 2


def test_artefact_flag_does_not_change_the_edge_count(tmp_path):
    rows = [
        ("P00001|10|S|21", "P00001|13|T|21", 0.9, 0.001),  # flagged
        ("P00002|10|S|21", "P00003|40|K|1", 0.9, 0.001),  # not flagged
    ]
    edges = parse(tmp_path, rows)

    assert count_edges(edges) == 2
    assert sorted(edges["potential_artefact"].to_list()) == [False, True]


def test_score_filter_uses_the_signed_score(tmp_path):
    rows = [
        ("P00001|10|S|21", "P00002|20|S|21", 0.7, 0.001),  # kept
        ("P00003|30|S|21", "P00004|40|S|21", -0.95, 0.001),  # anti-correlated, dropped
        ("P00005|50|S|21", "P00006|60|S|21", 0.59, 0.001),  # below threshold, dropped
    ]
    edges = parse(tmp_path, rows)

    assert edges["pair_key"].to_list() == ["P00001|10|S|21__P00002|20|S|21"]


def test_min_score_boundary_is_inclusive(tmp_path):
    edges = parse(
        tmp_path,
        [("P00001|10|S|21", "P00002|20|S|21", 0.6, 0.001)],
        min_score=0.6,
    )

    assert edges.height == 1


def test_qvalue_filter_is_strict(tmp_path):
    rows = [
        ("P00001|10|S|21", "P00002|20|S|21", 0.9, 0.049),  # kept
        ("P00003|30|S|21", "P00004|40|S|21", 0.9, 0.05),  # dropped, not < 0.05
    ]
    edges = parse(tmp_path, rows)

    assert edges["pair_key"].to_list() == ["P00001|10|S|21__P00002|20|S|21"]


def test_annotation_columns(tmp_path):
    rows = [
        # same protein, same modification, three residues apart
        ("P00001|10|S|21", "P00001|13|T|21", 0.9, 0.001),
        # same protein and position, different modification
        ("P00002|20|K|1", "P00002|20|K|34", 0.9, 0.001),
        # different proteins
        ("P00003|30|S|21", "P00004|40|S|21", 0.9, 0.001),
    ]
    edges = parse(tmp_path, rows).sort("pair_key")

    assert edges["same_protein"].to_list() == [True, True, False]
    assert edges["position_gap"].to_list() == [3, 0, None]
    assert edges["same_site"].to_list() == [False, True, False]
    assert edges["same_mod"].to_list() == [True, False, True]
    assert edges["shared_peptide"].to_list() == [True, True, False]
    assert edges["potential_artefact"].to_list() == [True, False, False]


def test_adjacency_window_boundary(tmp_path):
    rows = [
        ("P00001|10|S|21", "P00001|15|T|21", 0.9, 0.001),  # gap 5, inside the window
        ("P00002|10|S|21", "P00002|16|T|21", 0.9, 0.001),  # gap 6, outside it
    ]
    edges = parse(tmp_path, rows, adjacency_window=5).sort("pair_key")

    assert edges["position_gap"].to_list() == [5, 6]
    assert edges["shared_peptide"].to_list() == [True, False]
    assert edges["potential_artefact"].to_list() == [True, False]


def test_adjacency_window_is_configurable(tmp_path):
    rows = [("P00001|10|S|21", "P00001|20|T|21", 0.9, 0.001)]

    narrow = parse(tmp_path, rows, adjacency_window=5)
    wide = parse(tmp_path, rows, adjacency_window=10)

    assert narrow["potential_artefact"].to_list() == [False]
    assert wide["potential_artefact"].to_list() == [True]


def test_pair_key_is_independent_of_node_order(tmp_path):
    rows = [
        ("P00001|10|S|21", "P00002|20|S|21", 0.9, 0.001),
        ("P00002|20|S|21", "P00001|10|S|21", 0.8, 0.001),
    ]
    edges = parse(tmp_path, rows)

    # Both rows survive the filter; they describe one undirected edge, which is what
    # count_edges reports and what build_graph collapses.
    assert edges.height == 2
    assert set(edges["pair_key"]) == {"P00001|10|S|21__P00002|20|S|21"}
    assert count_edges(edges) == 1
    assert count_nodes(edges) == 2


@pytest.mark.parametrize(
    "node",
    [
        "P00001|10|S",  # too few fields
        "P00001|10|S|21|extra",  # too many fields
        "P00001|ten|S|21",  # non-numeric position
        "P00001||S|21",  # empty position
        "|10|S|21",  # empty accession
    ],
)
def test_malformed_node_identifier_raises(tmp_path, node):
    rows = [(node, "P00002|20|S|21", 0.9, 0.001)]

    with pytest.raises(ValueError, match="Malformed node identifier"):
        parse(tmp_path, rows)


def test_malformed_node_is_only_reported_after_filtering(tmp_path):
    # A malformed identifier on a row that the score filter removes never reaches the
    # network, so it must not fail the run.
    rows = [
        ("P00001|ten|S|21", "P00002|20|S|21", 0.1, 0.001),
        ("P00003|30|S|21", "P00004|40|S|21", 0.9, 0.001),
    ]
    edges = parse(tmp_path, rows)

    assert edges["pair_key"].to_list() == ["P00003|30|S|21__P00004|40|S|21"]
