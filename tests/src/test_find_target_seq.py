import pytest
from src.get_guiderna_seq import find_ct_target_seq, find_ag_target


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
    "NNCAGNNNNNNNNNNNNNNNNGGCGANNNNNNNNNNNNNNNNGG", #multipme candidates
    "NNCAGCAANNNNNNNNNNNNNGGTGG", #multipme candidates
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
def test_get_CtoT_target(seq, expected):
    assert find_ct_target_seq(seq) == expected


seq_not_candidates = [
    "NNNNNNNNNNNNNNNNNNNNNGG",
    "NNAGNNNNNNNNNNNNNNNNNGG",
    "NCTANNNNNNNNNNNNNNNNNGG",
    "NNNNCGANNNNNNNNNNNNNNGG", # no CGA in target window (-17~-19 from PAM)
    "NNNNCGANNNNNNNNNNNNNNAG", # no PAM
]
expected = [[],[],[],[],[]]
@pytest.mark.parametrize("seq, expected", zip(seq_not_candidates, expected))
def test_get_CtoT_target_not_candidate(seq, expected):
    assert find_ct_target_seq(seq) == expected


# G-to-A conversion

seq_candidates = [
    "CCNNNNNNNNNNNNNNNTGGNNN",
    "CCNNNNNNNNNNNNNNNNTGGNN",
    "CCNNNNNNNNNNNNNNNNNTGGN",
    "CCNNNNNNNNNNNNNNNNNNTGG",
    "CCNNNNNNNNNNNNNNNNNNTGGCCANNNNNNNNNNNNNNNNNTGG", #multipme candidates
    "CCNCCCNNNNNNNNNNNNNNTGGTGGNNNN", #multipme candidates
]
expected = [
    ["CCNNNNNNNNNNNNNNNTGGNNN"],
    ["CCNNNNNNNNNNNNNNNNTGGNN"],
    ["CCNNNNNNNNNNNNNNNNNTGGN"],
    ["CCNNNNNNNNNNNNNNNNNNTGG"],
    ["CCNNNNNNNNNNNNNNNNNNTGG", "CCANNNNNNNNNNNNNNNNNTGG"],["CCNCCCNNNNNNNNNNNNNNTGG", "CCCNNNNNNNNNNNNNNTGGTGG", "CCNNNNNNNNNNNNNNTGGTGGN"],
]
@pytest.mark.parametrize("seq, expected", zip(seq_candidates, expected))
def test_get_AtoG_target(seq, expected):
    assert find_ag_target(seq) == expected


seq_not_candidates = [
    "CCNNNNNNNNNNNNNNNTTTNNN", # no TGG
    "CCNNNNNNNNNNNNNNTGGNNNN", # no TGG in target window (+17~+19 from PAM)
    "CANNNNNNNNNNNNNNNTGGNNN", # no PAM
]
expected = [[],[],[]]
@pytest.mark.parametrize("seq, expected", zip(seq_not_candidates, expected))
def test_get_AtoG_target_not_candidate(seq, expected):
    assert find_ag_target(seq) == expected
