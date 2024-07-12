from __future__ import annotations
import pytest
from src.konezumiaid.nominate_ptc_guide.find_candidate_seq import (
    search_c_to_t_guide_seq,
    search_g_to_a_guide_seq,
)


#################################
# find candidate gRNA sequences
#################################

# C-to-T conversion

seq_candidates = [
    "NCAANNNNNNNNNNNNNNNNNGG",
    "NCAGNNNNNNNNNNNNNNNNNGG",
    "NCGANNNNNNNNNNNNNNNNNGG",
    "NNCAANNNNNNNNNNNNNNNNGG",
    "NNCAGNNNNNNNNNNNNNNNNGG",
    "NNCGANNNNNNNNNNNNNNNNGG",
    "NNNCAANNNNNNNNNNNNNNNGG",
    "NNNCAGNNNNNNNNNNNNNNNGG",
    "NNNCGANNNNNNNNNNNNNNNGG",
    "NNCAGNNNNNNNNNNNNNNNNGGCGANNNNNNNNNNNNNNNNGG",  # multipme candidates
    "NNCAGCAANNNNNNNNNNNNNGGTGG",  # multipme candidates
]
expected = [
    ["NCAANNNNNNNNNNNNNNNNNGG"],
    ["NCAGNNNNNNNNNNNNNNNNNGG"],
    ["NCGANNNNNNNNNNNNNNNNNGG"],
    ["NNCAANNNNNNNNNNNNNNNNGG"],
    ["NNCAGNNNNNNNNNNNNNNNNGG"],
    ["NNCGANNNNNNNNNNNNNNNNGG"],
    ["NNNCAANNNNNNNNNNNNNNNGG"],
    ["NNNCAGNNNNNNNNNNNNNNNGG"],
    ["NNNCGANNNNNNNNNNNNNNNGG"],
    ["NNCAGNNNNNNNNNNNNNNNNGG", "GGCGANNNNNNNNNNNNNNNNGG"],
    ["NNCAGCAANNNNNNNNNNNNNGG", "AGCAANNNNNNNNNNNNNGGTGG"],
]


@pytest.mark.parametrize("seq, expected", zip(seq_candidates, expected))
def test_search_c_to_t_guide_seq(seq, expected):
    assert search_c_to_t_guide_seq(seq) == expected


seq_not_candidates = [
    "NNNNNNNNNNNNNNNNNNNNNGG",
    "NNAGNNNNNNNNNNNNNNNNNGG",
    "NCTANNNNNNNNNNNNNNNNNGG",
    "NNNNCGANNNNNNNNNNNNNNGG",  # no CGA in target window (-17~-19 from PAM)
    "NNNNCGANNNNNNNNNNNNNNAG",  # no PAM
]
expected = [[], [], [], [], []]


@pytest.mark.parametrize("seq, expected", zip(seq_not_candidates, expected))
def test_non_search_c_to_t_guide_seq(seq, expected):
    assert search_c_to_t_guide_seq(seq) == expected


# G-to-A conversion

seq_candidates = [
    "CCNNNNNNNNNNNNNNNTGGNNN",
    "CCNNNNNNNNNNNNNNNNTGGNN",
    "CCNNNNNNNNNNNNNNNNNTGGN",
    "CCNNNNNNNNNNNNNNNNNNTGG",
    "CCNNNNNNNNNNNNNNNNNNTGGCCANNNNNNNNNNNNNNNNNTGG",  # multipme candidates
    "CCNCCCNNNNNNNNNNNNNNTGGTGGNNNN",  # multipme candidates
]
expected = [
    ["CCNNNNNNNNNNNNNNNTGGNNN"],
    ["CCNNNNNNNNNNNNNNNNTGGNN"],
    ["CCNNNNNNNNNNNNNNNNNTGGN"],
    ["CCNNNNNNNNNNNNNNNNNNTGG"],
    ["CCNNNNNNNNNNNNNNNNNNTGG", "CCANNNNNNNNNNNNNNNNNTGG"],
    ["CCNCCCNNNNNNNNNNNNNNTGG", "CCCNNNNNNNNNNNNNNTGGTGG", "CCNNNNNNNNNNNNNNTGGTGGN"],
]


@pytest.mark.parametrize("seq, expected", zip(seq_candidates, expected))
def test_search_g_to_a_guide_seq(seq, expected):
    assert search_g_to_a_guide_seq(seq) == expected


seq_not_candidates = [
    "CCNNNNNNNNNNNNNNNTTTNNN",  # no TGG
    "CCNNNNNNNNNNNNNNTGGNNNN",  # no TGG in target window (+17~+19 from PAM)
    "CANNNNNNNNNNNNNNNTGGNNN",  # no PAM
]
expected = [[], [], []]


@pytest.mark.parametrize("seq, expected", zip(seq_not_candidates, expected))
def test_non_search_g_to_a_guide_seq(seq, expected):
    assert search_g_to_a_guide_seq(seq) == expected
