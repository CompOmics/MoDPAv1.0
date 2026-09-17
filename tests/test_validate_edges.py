"""Tests for `3-pulse-silac-validation/validate_edges.py`.

These four functions implement the rule-based edge labelling that the pulsed SILAC
validation reports. The rules, taken from the repository README:

1. a positive association between two heavy or two light residues is valid;
2. a negative association between one heavy and one light residue is valid;
3. every other association is invalid.

Unlabelled residues, meaning any modification that is neither the light (Unimod 0) nor a
heavy (Unimod 259, 267) label, are excluded rather than counted as invalid.
"""
import numpy as np
import pandas as pd
import pytest

import validate_edges as ve


class TestGetLabel:
    @pytest.mark.parametrize(
        "node,expected",
        [
            ("P12345|10|K|0", "L"),    # light, Unimod 0
            ("P12345|10|R|0", "L"),
            ("P12345|10|K|259", "H"),  # heavy lysine, K8
            ("P12345|10|R|267", "H"),  # heavy arginine, R10
            ("P12345|10|M|35", "u"),   # oxidation, unlabelled
            ("P12345|10|C|4", "u"),
            ("P12345|10|S|21", "u"),
        ],
    )
    def test_unimod_suffix_determines_the_label(self, node, expected):
        assert ve.get_label(node) == expected

    def test_position_is_not_read_as_the_modification(self):
        # The suffix test must match the modification field, not a position that
        # happens to end in the same digits.
        assert ve.get_label("P12345|259|K|35") == "u"
        assert ve.get_label("P12345|267|R|35") == "u"
        assert ve.get_label("P12345|10|K|350") == "u"


class TestLabelToNum:
    @pytest.mark.parametrize("label,expected", [("H", 1), ("L", -1), ("u", 0)])
    def test_mapping(self, label, expected):
        assert ve.label_to_num(label) == expected

    def test_unknown_label_is_treated_as_unlabelled(self):
        assert ve.label_to_num("X") == 0


class TestValidateEdge:
    @staticmethod
    def row(label_a, label_b, score):
        return pd.Series({"labelA": label_a, "labelB": label_b, "Score": score})

    @pytest.mark.parametrize(
        "label_a,label_b,score,expected",
        [
            ("H", "H", 0.8, 1),    # rule 1, two heavy
            ("L", "L", 0.8, 1),    # rule 1, two light
            ("H", "L", -0.8, 1),   # rule 2, mixed and negative
            ("L", "H", -0.8, 1),
            ("H", "H", -0.8, 0),   # rule 3, same label but negative
            ("L", "L", -0.8, 0),
            ("H", "L", 0.8, 0),    # rule 3, mixed but positive
            ("H", "u", 0.8, 0),    # unlabelled endpoint
            ("u", "L", -0.8, 0),
            ("u", "u", 0.8, 0),
        ],
    )
    def test_rule_table(self, label_a, label_b, score, expected):
        assert ve.validate_edge(self.row(label_a, label_b, score)) == expected

    def test_zero_score_is_not_validated(self):
        # A score of exactly zero satisfies neither the > 0 nor the < 0 branch.
        assert ve.validate_edge(self.row("H", "H", 0.0)) == 0
        assert ve.validate_edge(self.row("H", "L", 0.0)) == 0

    def test_is_symmetric_in_the_two_endpoints(self):
        for label_a, label_b, score in [("H", "L", -0.5), ("H", "H", 0.5), ("u", "H", 0.5)]:
            assert ve.validate_edge(self.row(label_a, label_b, score)) == ve.validate_edge(
                self.row(label_b, label_a, score)
            )


class TestCountValidated:
    @pytest.fixture
    def edges(self):
        return pd.DataFrame(
            {
                "abs_corr": [0.9, 0.9, 0.4, 0.9, 0.8],
                "network_type": ["MoDPA", "MoDPA", "MoDPA", "random", "random"],
                "validated": [1, 0, 1, 1, 1],
            }
        )

    def test_filters_by_threshold_and_network(self, edges):
        # Only the two MoDPA edges at or above 0.5 are counted: the mean of [1, 0].
        assert ve.count_validated(edges, 0.5, "MoDPA") == 0.5

    def test_threshold_is_inclusive(self, edges):
        # The 0.4 edge enters once the threshold drops to 0.4: the mean of [1, 0, 1].
        assert ve.count_validated(edges, 0.4, "MoDPA") == pytest.approx(2 / 3)

    def test_selects_the_requested_network(self, edges):
        assert ve.count_validated(edges, 0.5, "random") == 1.0

    def test_empty_selection_returns_nan(self, edges):
        assert np.isnan(ve.count_validated(edges, 0.99, "MoDPA"))
        assert np.isnan(ve.count_validated(edges, 0.5, "degree_preserved"))

    def test_all_nan_validated_returns_nan(self):
        df = pd.DataFrame(
            {"abs_corr": [0.9], "network_type": ["MoDPA"], "validated": [np.nan]}
        )
        assert np.isnan(ve.count_validated(df, 0.5, "MoDPA"))

    def test_nan_rows_are_dropped_not_counted_as_zero(self):
        df = pd.DataFrame(
            {
                "abs_corr": [0.9, 0.9],
                "network_type": ["MoDPA", "MoDPA"],
                "validated": [1.0, np.nan],
            }
        )
        assert ve.count_validated(df, 0.5, "MoDPA") == 1.0

    def test_input_is_not_mutated(self, edges):
        before = edges.copy()
        ve.count_validated(edges, 0.5, "MoDPA")
        pd.testing.assert_frame_equal(edges, before)


class TestEndToEndOnASmallNetwork:
    """The three README rules applied to a small edge table, as the script does it."""

    def test_labelled_edges_are_scored_by_the_rules(self):
        edges = pd.DataFrame(
            {
                "nodeA": ["P1|10|K|259", "P1|20|K|0", "P1|30|K|259", "P1|40|M|35"],
                "nodeB": ["P2|11|R|267", "P2|21|R|0", "P2|31|R|0", "P2|41|S|21"],
                "Score": [0.7, 0.6, -0.8, 0.9],
            }
        )
        edges["labelA"] = edges.nodeA.map(ve.get_label)
        edges["labelB"] = edges.nodeB.map(ve.get_label)
        edges["validated"] = edges.apply(ve.validate_edge, axis=1)

        assert list(edges.labelA) == ["H", "L", "H", "u"]
        assert list(edges.labelB) == ["H", "L", "L", "u"]
        # heavy/heavy positive, light/light positive, heavy/light negative, unlabelled.
        assert list(edges.validated) == [1, 1, 1, 0]
