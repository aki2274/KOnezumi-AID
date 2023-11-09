from __future__ import annotations
import pytest
from src.get_candidate_stopcodon_index import get_candidate_stopcodon_index

test_candidate_cds_seq = [
    "ATGCAANNN",
    "ATGCAGNNNNN",
    "ATGCGANNNNN",
    "ATGTGGNNNNN",
    "ATGCAANNNNNNTGG",  #
]


expected = [[3], [3], [3], [3], [3, 12]]


@pytest.mark.parametrize(
    "test_candidate_cds_seq,expected",
    zip(test_candidate_cds_seq, expected),
)
def test_get_candidate_stopcodon_index(test_candidate_cds_seq, expected):
    assert get_candidate_stopcodon_index(test_candidate_cds_seq) == expected


test_not_candidate_cds_seq = ["ATGNCAANNN", "ATGNNCAANNN", "ATGNNNCGG"]


expected = [[], [], []]


@pytest.mark.parametrize(
    "test_not_candidate_cds_seq,expected",
    zip(test_not_candidate_cds_seq, expected),
)
def test_get_not_candidate_stopcodon_index(test_not_candidate_cds_seq, expected):
    assert get_candidate_stopcodon_index(test_not_candidate_cds_seq) == expected
