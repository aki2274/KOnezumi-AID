from __future__ import annotations
import pytest
from konezumiaid.nominate_candidate_stopcodon.search_candidate_codon import (
    search_candidate_index,
)

test_candidate_cds_seq = [
    "ATGCAANNN",
    "ATGCAGNNNNN",
    "ATGCGANNNNN",
    "ATGTGGNNNNN",
    "ATGCAANNNNNNTGG",  # multiple candidate stopcodon
]


expected = [[3], [3], [3], [3], [3, 12]]


@pytest.mark.parametrize(
    "test_candidate_cds_seq,expected",
    zip(test_candidate_cds_seq, expected),
)
def test_get_candidate_stopcodon_index(test_candidate_cds_seq, expected):
    assert search_candidate_index(test_candidate_cds_seq) == expected


test_not_candidate_cds_seq = [
    "ATGNCAANNN",  # CAA is not a codon (mod3 = 1)
    "ATGNNCAANNN",  # CAA is not a codon (mod3 != 2)
    "ATGNNNCGG",  # CGG is not candidate stopcodon
]


expected = [[], [], []]


@pytest.mark.parametrize(
    "test_not_candidate_cds_seq,expected",
    zip(test_not_candidate_cds_seq, expected),
)
def test_get_not_candidate_stopcodon_index(test_not_candidate_cds_seq, expected):
    assert search_candidate_index(test_not_candidate_cds_seq) == expected
