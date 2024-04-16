import pytest
from konezumiaid.get_reverse_complement import get_revcomp

# small letter means low complexity region
input_seq = [
    "ATGCAAATTTGGGCCC",
    "atgcaaatttgggccc",
    "NNNNNNNNNN",
    "ATCGN",
    "atgc",
    "nnnnn",
    "ATGCatgc",
]

expected = [
    "GGGCCCAAATTTGCAT",
    "gggcccaaatttgcat",
    "NNNNNNNNNN",
    "NCGAT",
    "gcat",
    "nnnnn",
    "gcatGCAT",
]


@pytest.mark.parametrize("seq,expected", zip(input_seq, expected))
def test_get_revcomp(seq, expected):
    assert get_revcomp(seq) == expected
