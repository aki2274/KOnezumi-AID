from __future__ import annotations
import pytest
from src.get_primer_isin_num_of_exon import get_primer_exon_num

primer_index_list = [
    [2, 7, 12, 17],
    [3, 8, 13, 18],
    [0, 5, 10, 15],  #  edge case
    [4, 9, 14, 19],  #  edge case
]
exon_range = [[0, 4], [5, 9], [10, 14], [15, 19]]
expected = [[0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3], [0, 1, 2, 3]]


@pytest.mark.parametrize(
    "primer_index_list,exon_range,expected",
    zip(primer_index_list, [exon_range] * len(primer_index_list), expected),
)
def test_get_primer_exon_num_param(primer_index_list, exon_range, expected):
    assert get_primer_exon_num(primer_index_list, exon_range) == expected
