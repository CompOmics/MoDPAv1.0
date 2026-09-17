"""Tests for the pair-encoding and binning helpers in `4-PTM-pairs-annotation/scripts/s15_background.py`.

These functions carry the matched permutation background. `pack` encodes an unordered
pair of PTM events as a single int64 so that membership of an evidence set becomes one
sorted-array lookup; a collision there would silently move pairs between tiers. `tier_of`
is a vectorised restatement of the SQL tier assignment in s13, and its weakest-to-strongest
ordering is what reproduces the first-match-wins semantics of that CASE expression.

The functions close over module-level globals in the original script, so each test injects
the globals it needs into the extracted module.
"""
import numpy as np
import pandas as pd
import pytest

from conftest import load_functions

FUNCTIONS = ["pack", "packp", "isin_sorted", "dec", "tier_of"]


@pytest.fixture
def s15():
    """A fresh copy of the extracted functions, so injected globals do not leak."""
    return load_functions("4-PTM-pairs-annotation/scripts/s15_background.py", FUNCTIONS)


class TestPack:
    def test_encoding_is_collision_free(self, s15):
        s15.N = 40
        lo, hi = np.triu_indices(s15.N, k=1)
        keys = s15.pack(lo, hi)
        assert len(np.unique(keys)) == len(keys)

    def test_encoding_is_order_independent(self, s15):
        s15.N = 40
        lo, hi = np.triu_indices(s15.N, k=1)
        assert np.array_equal(s15.pack(lo, hi), s15.pack(hi, lo))

    def test_self_pairs_do_not_collide_with_distinct_pairs(self, s15):
        s15.N = 10
        i, j = np.meshgrid(np.arange(10), np.arange(10))
        keys = s15.pack(i.ravel(), j.ravel())
        # min*N + max is unique over all ordered pairs once the order is normalised.
        assert len(np.unique(keys)) == 10 * 11 // 2

    def test_result_is_int64(self, s15):
        s15.N = 10
        assert s15.pack(np.array([1]), np.array([2])).dtype == np.int64

    def test_packp_uses_the_protein_count(self, s15):
        s15.N, s15.P = 1000, 7
        lo, hi = np.triu_indices(s15.P, k=1)
        keys = s15.packp(lo, hi)
        assert len(np.unique(keys)) == len(keys)
        # packp must not silently reuse N.
        assert not np.array_equal(keys, s15.pack(lo, hi))


class TestIsInSorted:
    def test_empty_reference_matches_nothing(self, s15):
        result = s15.isin_sorted(np.array([1, 2, 3]), np.array([], dtype=np.int64))
        assert result.dtype == bool
        assert not result.any()

    def test_membership(self, s15):
        reference = np.array([10, 20, 30])
        assert list(s15.isin_sorted(np.array([10, 20, 30]), reference)) == [True] * 3
        assert list(s15.isin_sorted(np.array([5, 15, 25]), reference)) == [False] * 3

    def test_value_above_the_reference_range_is_not_matched(self, s15):
        # searchsorted returns an insertion point past the end; the clamp to 0 must not
        # make such a value compare equal to the first element.
        reference = np.array([10, 20, 30])
        assert list(s15.isin_sorted(np.array([35, 10]), reference)) == [False, True]

    def test_agrees_with_numpy_isin(self, s15):
        rng = np.random.default_rng(0)
        reference = np.unique(rng.integers(0, 500, 120))
        values = rng.integers(0, 600, 400)
        assert np.array_equal(s15.isin_sorted(values, reference), np.isin(values, reference))


class TestDec:
    def test_range_is_zero_to_nine(self, s15):
        rng = np.random.default_rng(0)
        deciles = s15.dec(rng.normal(size=1000))
        assert deciles.min() == 0 and deciles.max() == 9

    def test_bins_are_nearly_equally_populated(self, s15):
        # Ranks are 1-based, so the first decile holds nine values and the last eleven.
        # Pinning this keeps a change in the binning visible.
        counts = np.bincount(s15.dec(np.arange(100)), minlength=10)
        assert list(counts) == [9] + [10] * 8 + [11]

    def test_ties_share_a_decile(self, s15):
        assert set(s15.dec(np.full(20, 7.0))) == {5}

    def test_monotonic_in_the_input(self, s15):
        deciles = s15.dec(np.arange(100))
        assert np.all(np.diff(deciles) >= 0)

    def test_shape_independent_of_distribution(self, s15):
        # Both inputs have the same ranks, so they must produce the same deciles.
        linear = s15.dec(np.arange(50))
        exponential = s15.dec(np.exp(np.arange(50)))
        assert np.array_equal(linear, exponential)


class TestTierOf:
    """Tier precedence, on a six-event fixture spanning three proteins.

    Events 0 and 1 sit on protein 0, events 2 and 3 on protein 1, events 4 and 5 on
    protein 2. Events 0-3 carry a modification annotation, events 4 and 5 do not.
    """

    @pytest.fixture
    def tiers(self, s15):
        s15.N, s15.P = 6, 3
        s15.acc_codes = np.array([0, 0, 1, 1, 2, 2])
        s15.ann = np.array([True, True, True, True, False, False])
        s15.expv = np.array([True, False, True, False, False, False])
        s15.htpv = np.array([False, True, False, False, False, False])
        s15.S_T1 = np.array([], dtype=np.int64)
        s15.S_T2 = np.array([], dtype=np.int64)
        s15.S_T4 = np.array([], dtype=np.int64)
        return s15

    @staticmethod
    def tier(module, i, j):
        return int(module.tier_of(np.array([i]), np.array([j]))[0])

    def test_unannotated_pair_is_unsupported(self, tiers):
        assert self.tier(tiers, 4, 5) == 9

    def test_both_annotated_and_evidenced_across_proteins_is_tier_5(self, tiers):
        # Events 0 and 2 are on different proteins, both annotated, both with evidence.
        assert self.tier(tiers, 0, 2) == 5

    def test_annotation_without_evidence_is_not_tier_5(self, tiers):
        # Event 3 is annotated but carries neither experimental nor high-throughput support.
        assert self.tier(tiers, 0, 3) == 9

    def test_protein_pair_evidence_gives_tier_4(self, tiers):
        tiers.S_T4 = np.unique(tiers.packp(np.array([0]), np.array([1])))
        assert self.tier(tiers, 0, 2) == 4

    def test_tier_4_does_not_apply_within_one_protein(self, tiers):
        tiers.S_T4 = np.unique(tiers.packp(np.array([0]), np.array([0])))
        # Events 0 and 1 share protein 0, so the complex/interaction tier is skipped and
        # the same-protein tier applies instead.
        assert self.tier(tiers, 0, 1) == 3

    def test_same_protein_both_annotated_is_tier_3(self, tiers):
        assert self.tier(tiers, 0, 1) == 3
        assert self.tier(tiers, 2, 3) == 3

    def test_shared_regulator_set_gives_tier_2(self, tiers):
        tiers.S_T2 = np.unique(tiers.pack(np.array([0]), np.array([1])))
        assert self.tier(tiers, 0, 1) == 2

    def test_explicit_crosstalk_set_gives_tier_1(self, tiers):
        tiers.S_T1 = np.unique(tiers.pack(np.array([0]), np.array([1])))
        assert self.tier(tiers, 0, 1) == 1

    def test_stronger_tier_wins(self, tiers):
        pair = np.unique(tiers.pack(np.array([0]), np.array([1])))
        tiers.S_T1 = pair
        tiers.S_T2 = pair
        tiers.S_T4 = np.unique(tiers.packp(np.array([0]), np.array([0])))
        assert self.tier(tiers, 0, 1) == 1

    def test_order_of_endpoints_does_not_matter(self, tiers):
        tiers.S_T2 = np.unique(tiers.pack(np.array([0]), np.array([1])))
        assert self.tier(tiers, 0, 1) == self.tier(tiers, 1, 0)

    def test_returns_int8(self, tiers):
        assert tiers.tier_of(np.array([0]), np.array([1])).dtype == np.int8

    def test_vectorises_over_many_pairs(self, tiers):
        i = np.array([0, 0, 4, 2])
        j = np.array([1, 2, 5, 3])
        assert list(tiers.tier_of(i, j)) == [3, 5, 9, 3]
