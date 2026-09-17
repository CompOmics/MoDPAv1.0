"""Tests for the tryptic co-peptide logic in `4-PTM-pairs-annotation/scripts/s10_tryptic.py`.

`cuts_of` and `shares` decide the `same_peptide_std` and `same_peptide_modblocked` flags,
which report how many associated PTM pairs could be quantified from overlapping PSM sets.
Both are pure functions and are tested here against hand-computed digests.

Digestion rules taken from the script: trypsin, cleavage after K or R, no cleavage before
proline, up to two missed cleavages.
"""
import pytest

from conftest import load_functions

s10 = load_functions("4-PTM-pairs-annotation/scripts/s10_tryptic.py", ["cuts_of", "shares"])
cuts_of = s10.cuts_of
shares = s10.shares


class TestCutsOf:
    def test_cuts_after_k_and_r(self):
        # AAK|AAR|AAA, so cleavage sites at 3 and 6, bracketed by 0 and 9.
        assert cuts_of("AAKAARAAA") == [0, 3, 6, 9]

    def test_no_cleavage_before_proline(self):
        # The K at position 3 is followed by P, so only the R at 7 cuts.
        assert cuts_of("AAKPAARAAA") == [0, 7, 10]

    def test_blocked_position_is_not_cleaved(self):
        # Blocking position 3 models a modified lysine that trypsin cannot process.
        assert cuts_of("AAKAARAAA", block={3}) == [0, 6, 9]

    def test_blocking_every_site_leaves_the_whole_protein(self):
        assert cuts_of("AAKAARAAA", block={3, 6}) == [0, 9]

    def test_c_terminal_k_is_not_duplicated(self):
        # The cut after the final K coincides with the closing sentinel.
        assert cuts_of("AAK") == [0, 3]

    def test_c_terminal_proline_does_not_suppress_the_sentinel(self):
        # K at 3 is followed by P, so it does not cut, but the sentinel still closes.
        assert cuts_of("AAKP") == [0, 4]

    def test_protein_without_k_or_r(self):
        assert cuts_of("AAAA") == [0, 4]

    def test_empty_sequence(self):
        assert cuts_of("") == [0]

    def test_block_only_affects_listed_positions(self):
        unblocked = cuts_of("AAKAARAAA")
        blocked = cuts_of("AAKAARAAA", block={3})
        assert set(blocked) < set(unblocked)


class TestShares:
    # Cuts of a protein digested into four peptides: 1-3, 4-6, 7-9, 10-12.
    CUTS = [0, 3, 6, 9, 12]

    @pytest.mark.parametrize(
        "p1,p2,expected_missed_cleavages",
        [
            (1, 2, 0),   # both inside peptide 1-3
            (1, 3, 0),   # first and last residue of one peptide
            (2, 5, 1),   # peptides 1-3 and 4-6
            (3, 4, 1),   # residues either side of a cleavage site
            (2, 8, 2),   # peptides 1-3 through 7-9, at the default limit
            (2, 11, 3),  # peptides 1-3 through 10-12, above the default limit
        ],
    )
    def test_default_limit_is_two_missed_cleavages(self, p1, p2, expected_missed_cleavages):
        assert shares(self.CUTS, p1, p2) is (expected_missed_cleavages <= 2)

    def test_max_mc_argument_is_respected(self):
        # The same pair crosses three cleavage sites, so it needs max_mc >= 3.
        assert not shares(self.CUTS, 2, 11, max_mc=2)
        assert shares(self.CUTS, 2, 11, max_mc=3)

    def test_zero_missed_cleavages_requires_one_peptide(self):
        assert shares(self.CUTS, 1, 3, max_mc=0)
        assert not shares(self.CUTS, 3, 4, max_mc=0)

    def test_identical_positions_always_share(self):
        for position in (1, 5, 9, 12):
            assert shares(self.CUTS, position, position)

    def test_position_on_a_cleavage_site_belongs_to_the_preceding_peptide(self):
        # Residue 3 is the last residue of peptide 1-3, not the first of 4-6.
        assert shares(self.CUTS, 1, 3, max_mc=0)
        assert not shares(self.CUTS, 3, 6, max_mc=0)

    def test_requires_p1_not_greater_than_p2(self):
        # Documented precondition: the caller sorts the two positions before calling.
        # With the arguments reversed the two binary searches cross and the result is
        # meaningless, so this pins the precondition rather than the reversed answer.
        assert not shares(self.CUTS, 2, 11)
        assert shares(self.CUTS, 11, 2)  # wrong answer, reached only by misuse


class TestBlockingIsAnUpperBound:
    """Blocking modified K and R can only merge peptides, never split them.

    `same_peptide_modblocked` is therefore always at least `same_peptide_std`, which is
    what makes it usable as the upper bound the script describes.
    """

    SEQUENCE = "AAKAARAAKAARAAKAARAAA"

    @pytest.mark.parametrize("p1,p2", [(1, 5), (1, 12), (2, 18), (5, 20), (1, 21)])
    def test_blocked_digest_shares_at_least_as_often(self, p1, p2):
        standard = cuts_of(self.SEQUENCE)
        blocked = cuts_of(self.SEQUENCE, block={3, 6, 9})
        assert shares(blocked, p1, p2) >= shares(standard, p1, p2)

    def test_blocking_turns_a_separated_pair_into_a_shared_one(self):
        standard = cuts_of(self.SEQUENCE)       # cuts at 3, 6, 9, 12, 15, 18
        blocked = cuts_of(self.SEQUENCE, block={3, 6, 9, 12})
        assert not shares(standard, 1, 14)
        assert shares(blocked, 1, 14)
