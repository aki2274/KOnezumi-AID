from __future__ import annotations
import pytest
from konezumiaid.nominate_candidate_stopcodon.transform_candidateseq_to_index import (
    transform_ct_guideseq_to_index,
    transform_ga_guideseq_to_index,
)


#################################
# Get start index of target base
#################################

# C-to-T conversion

input_seq = [
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
    "NNCAGCAANNNNNNNNNNNNNGGGGG" # multipme candidates and G rich around PAM 
]
seq_candidates = [
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
    ["NNCAGCAANNNNNNNNNNNNNGG", "CAGCAANNNNNNNNNNNNNGGGG"],
]

expected = [[1], [1], [1], [2], [2], [2], [3], [3], [3], [2, 23], [2, 5],[2, 5]]


@pytest.mark.parametrize(
    "seq,candidates,expected", zip(input_seq, seq_candidates, expected)
)
def test_get_CtoT_target(seq, candidates, expected):
    assert transform_ct_guideseq_to_index(seq, candidates) == expected


input_seq = [
    "NNNNNNNNNNNNNNNNNNNNNGG",
    "NNAGNNNNNNNNNNNNNNNNNGG",
    "NCTANNNNNNNNNNNNNNNNNGG",
    "NNNNCGANNNNNNNNNNNNNNGG",  # no CGA in target window (-17~-19 from PAM)
    "NNNNCGANNNNNNNNNNNNNNAG",  # no PAM
]
seq_not_candidates = [[], [], [], [], []]
expected = [[], [], [], [], []]


@pytest.mark.parametrize(
    "seq,candidates,expected", zip(input_seq, seq_not_candidates, expected)
)
def test_get_CtoT_target_not_candidate(seq, candidates, expected):
    assert transform_ct_guideseq_to_index(seq, candidates) == expected


# G-to-A conversion

input_seq = [
    "CCNNNNNNNNNNNNNNNTGGNNN",
    "CCNNNNNNNNNNNNNNNNTGGNN",
    "CCNNNNNNNNNNNNNNNNNTGGN",
    "CCNNNNNNNNNNNNNNNNNNTGG",
    "CCNNNNNNNNNNNNNNNNNNTGGCCANNNNNNNNNNNNNNNNNTGG",  # multipme candidates
    "CCNCCCNNNNNNNNNNNNNNTGGTGGNNNN",  # multipme candidates
    "CCNNNNNNNNNNNNNNNTGGTGGNN" ]# multipme candidates in the same candidate
seq_candidates = [
    ["CCNNNNNNNNNNNNNNNTGGNNN"],
    ["CCNNNNNNNNNNNNNNNNTGGNN"],
    ["CCNNNNNNNNNNNNNNNNNTGGN"],
    ["CCNNNNNNNNNNNNNNNNNNTGG"],
    ["CCNNNNNNNNNNNNNNNNNNTGG", "CCANNNNNNNNNNNNNNNNNTGG"],
    ["CCNCCCNNNNNNNNNNNNNNTGG", "CCCNNNNNNNNNNNNNNTGGTGG", "CCNNNNNNNNNNNNNNTGGTGGN"],
    ["CCNNNNNNNNNNNNNNNTGGTGG"]
]
expected = [[17], [18], [19], [20], [20, 43], [20, 23],[20, 23]]


@pytest.mark.parametrize(
    "seq,candidates,expected", zip(input_seq, seq_candidates, expected)
)
def test_get_AtoG_target(seq, candidates, expected):
    assert transform_ga_guideseq_to_index(seq, candidates) == expected


input_seq = [
    "CCNNNNNNNNNNNNNNNTTTNNN",  # no TGG
    "CCNNNNNNNNNNNNNNTGGNNNN",  # no TGG in target window (+17~+19 from PAM)
    "CANNNNNNNNNNNNNNNTGGNNN",  # no PAM
]
seq_not_candidates = [[], [], []]
expected = [[], [], []]


@pytest.mark.parametrize(
    "seq,candidates,expected", zip(input_seq, seq_not_candidates, expected)
)
def test_get_AtoG_target_not_candidates(seq, candidates, expected):
    assert transform_ga_guideseq_to_index(seq, candidates) == expected
