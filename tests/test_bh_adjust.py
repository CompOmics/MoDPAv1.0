"""Tests for `bh_adjust` in `5-pathway-ORA/reactome_network_background_ora.py`.

Every FDR reported by the network-background sensitivity analysis comes out of this
function. It is checked against `scipy.stats.false_discovery_control`, which implements
the same Benjamini-Hochberg procedure, and against hand-computed values.

The NaN handling is specific to this analysis: a pathway that could not be tested carries
a NaN p-value, and those rows must pass through without entering the correction or
inflating the number of tests.
"""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

import reactome_network_background_ora as ora

bh_adjust = ora.bh_adjust


class TestAgainstScipy:
    @pytest.mark.parametrize("seed", [0, 1, 2, 3])
    def test_matches_scipy_on_random_p_values(self, seed):
        rng = np.random.default_rng(seed)
        pvals = pd.Series(rng.uniform(0, 1, 200))
        expected = false_discovery_control(pvals.to_numpy(), method="bh")
        assert np.allclose(bh_adjust(pvals).to_numpy(), expected)

    def test_matches_scipy_with_ties(self, seed=0):
        rng = np.random.default_rng(seed)
        pvals = pd.Series(np.round(rng.uniform(0, 1, 200), 2))
        expected = false_discovery_control(pvals.to_numpy(), method="bh")
        assert np.allclose(bh_adjust(pvals).to_numpy(), expected)

    def test_matches_scipy_on_already_sorted_input(self):
        pvals = pd.Series(np.linspace(0.001, 0.9, 50))
        expected = false_discovery_control(pvals.to_numpy(), method="bh")
        assert np.allclose(bh_adjust(pvals).to_numpy(), expected)


class TestHandComputed:
    def test_five_p_values(self):
        # ranked p * m / rank, then the running minimum from the largest rank down:
        # 0.01*5/1=0.050, 0.02*5/2=0.050, 0.03*5/3=0.050, 0.04*5/4=0.050, 0.05*5/5=0.050
        result = bh_adjust(pd.Series([0.01, 0.02, 0.03, 0.04, 0.05]))
        assert np.allclose(result.to_numpy(), [0.05] * 5)

    def test_step_up_enforces_monotonicity(self):
        # Raw values are 0.004*3/1=0.012, 0.5*3/2=0.75, 0.6*3/3=0.6; the running minimum
        # from the top pulls the middle value down to 0.6.
        result = bh_adjust(pd.Series([0.004, 0.5, 0.6]))
        assert np.allclose(result.to_numpy(), [0.012, 0.6, 0.6])

    def test_single_p_value_is_unchanged(self):
        assert bh_adjust(pd.Series([0.5])).to_numpy() == pytest.approx(0.5)

    def test_large_p_values_are_pulled_down_by_the_step_up(self):
        # Raw values are 0.9*3/1=2.7, 0.95*3/2=1.425, 0.99*3/3=0.99. The running minimum
        # from the largest rank down sets all three to 0.99, so nothing reaches the clip.
        result = bh_adjust(pd.Series([0.9, 0.95, 0.99]))
        assert np.allclose(result.to_numpy(), [0.99, 0.99, 0.99])

    def test_result_never_exceeds_one(self):
        result = bh_adjust(pd.Series([0.9, 0.95, 1.0]))
        assert (result <= 1.0).all()
        assert np.allclose(result.to_numpy(), [1.0, 1.0, 1.0])


class TestNaNHandling:
    def test_nan_passes_through(self):
        result = bh_adjust(pd.Series([0.01, np.nan, 0.02]))
        assert np.isnan(result.iloc[1])
        assert result.notna().sum() == 2

    def test_nan_rows_do_not_count_towards_m(self):
        # Two testable p-values, so m is 2 and not 4.
        with_nan = bh_adjust(pd.Series([0.01, np.nan, 0.02, np.nan])).dropna()
        without_nan = bh_adjust(pd.Series([0.01, 0.02]))
        assert np.allclose(with_nan.to_numpy(), without_nan.to_numpy())

    def test_all_nan_returns_all_nan(self):
        result = bh_adjust(pd.Series([np.nan, np.nan, np.nan]))
        assert result.isna().all()
        assert len(result) == 3

    def test_input_is_not_mutated(self):
        pvals = pd.Series([0.01, 0.5, np.nan])
        before = pvals.copy()
        bh_adjust(pvals)
        pd.testing.assert_series_equal(pvals, before)


class TestStructure:
    def test_index_is_preserved(self):
        pvals = pd.Series([0.3, 0.01, 0.2], index=["R-HSA-1", "R-HSA-2", "R-HSA-3"])
        result = bh_adjust(pvals)
        assert list(result.index) == list(pvals.index)

    def test_non_unique_index_is_preserved(self):
        pvals = pd.Series([0.3, 0.01, 0.2], index=["a", "a", "b"])
        assert list(bh_adjust(pvals).index) == ["a", "a", "b"]

    def test_adjusted_values_are_never_below_the_raw_value(self):
        rng = np.random.default_rng(7)
        pvals = pd.Series(rng.uniform(0, 1, 100))
        assert (bh_adjust(pvals) >= pvals - 1e-12).all()

    def test_ordering_is_preserved(self):
        # The adjustment is monotone, so the rank order of the p-values survives it.
        rng = np.random.default_rng(3)
        pvals = pd.Series(rng.uniform(0, 1, 100))
        adjusted = bh_adjust(pvals)
        assert np.all(np.diff(adjusted.to_numpy()[np.argsort(pvals.to_numpy())]) >= -1e-12)
