import pytest
from src.get_guiderna_seq import get_revcomp

# small letter means low complexity region
input_seq = ["ATGCAAATTTGGGCCC", "NNNNNNNNNN", "ATCGN", "atgc", "nnn", "ATGCatgc"]

expected = ["GGGCCCAAATTTGCAT", "NNNNNNNNNN", "NCGAT", "GCAT", "NNN", "GCATGCAT"]


@pytest.mark.parametrize("seq,expected", zip(input_seq, expected))
def test_get_revcomp(seq, expected):
    assert get_revcomp(seq) == expected


# test the case that the unexpected character is included in the sequence.
# このエラーのテストが必要かどうかについて検討
"""
error_input_seq = ["ATGCQ", "ATGC>", "ATGCw", "ATGC1"]

expecte_error = ["KeyError: 'Q'", "KeyError: '>'", "KeyError: 'W'", "KeyError: '1'"]


@pytest.mark.parametrize("seq,expected", zip(error_input_seq, expecte_error))
def test_get_revcomp_error(seq, expected):
    assert get_revcomp(seq) == expected
"""
