from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.apply_nmd_rules.evaluate_primers import (
    eliminate_in_back_half,
    eliminate_in_last_exon,
    label_in_start_150bp,
    label_in_50bp_from_LEJ,
    create_candidates_list_dict
)


def test_create_candidate_list():
    ct_cand = [
        "CTGCAGAACTCGGACAATTCGGG",
        "GCAGCGGCGTTACTCGCTGTGGG",
        "AGCAGCGGCGTTACTCGCTGTGG",
        "TGCAAGTGTTCAGCTTCCGCTGG",
        "TGCAGCTTGGGCAAATCTGGAGG",
        "GGCAGAGTGGGGAGAGCGGCAGG",
    ]
    ga_cand = [
        "CCTGCGTGGCCAGCGCTGGTGGT",
        "CCTGCTCCTTTTGCATCTGGCTC",
        "CCTCCCTTGTGTCCTTGGCTTGG",
        "CCCTTGTGTCCTTGGCTTGGGCC",
    ]

    expected = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ct_seq": "AGCAGCGGCGTTACTCGCTGTGG"},
        {"ct_seq": "TGCAAGTGTTCAGCTTCCGCTGG"},
        {"ct_seq": "TGCAGCTTGGGCAAATCTGGAGG"},
        {"ct_seq": "GGCAGAGTGGGGAGAGCGGCAGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
        {"ga_seq": "CCTGCTCCTTTTGCATCTGGCTC"},
        {"ga_seq": "CCTCCCTTGTGTCCTTGGCTTGG"},
        {"ga_seq": "CCCTTGTGTCCTTGGCTTGGGCC"},
    ]

    assert create_candidates_list_dict(ct_cand, ga_cand) == expected

def test_label_in_start_150bp():
    cand_grna = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
    ]
    ds = GeneData(
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGCAGAACTCGGACAATTCGGGNNNNNNNNNNNNNNNNNNNNNNNGCAGCGGCGTTACTCGCTGTGGGCCTGCGTGGCCAGCGCTGGTGGT",
        0,
        229,
        0,
        229,
        1,
        [0],
        [229],
    )
    expected = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG", "in_start_150bp": True},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG", "in_start_150bp": True},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT", "in_start_150bp": False},
    ]
    assert label_in_start_150bp(cand_grna, ds) == expected



def test_eliminate_in_front_half():
    cand_grna = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
    ]
    ds = GeneData(
        "CTGCAGAACTCGGACAATTCGGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAGCGGCGTTACTCGCTGTGGGCCTGCGTGGCCAGCGCTGGTGGT",
        0,
        111,
        0,
        111,
        1,
        [0],
        [111],
    )
    expected = [{"ct_seq": "CTGCAGAACTCGGACAATTCGGG"}]
    assert eliminate_in_back_half(cand_grna, ds) == expected



def test_eliminate_in_last_exon():
    cand_grna = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
    ]
    ds = GeneData(
        "CTGCAGAACTCGGACAATTCGGGNNNNNNNNNNNNNNNNNNNNNNNGCAGCGGCGTTACTCGCTGTGGGCCTGCGTGGCCAGCGCTGGTGGT",
        0,
        91,
        0,
        91,
        4,
        [0, 23, 46, 69],
        [22, 45, 68, 91],
    )
    excepted = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
    ]
    assert eliminate_in_last_exon(cand_grna, ds) == excepted


def test_label_in_50bp_from_LEJ():
    cand_grna = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
    ]
    ds = GeneData(
        "CTGCAGAACTCGGACAATTCGGGNNNNNNNNNNNNNNNNNNNNNNNGCAGCGGCGTTACTCGCTGTGGGCCTGCGTGGCCAGCGCTGGTGGT",
        0,
        91,
        0,
        91,
        4,
        [0, 23, 46, 69],
        [22, 45, 68, 91],
    )
    expected = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG", "in_50bp_from_LEJ": False},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG", "in_50bp_from_LEJ": True},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT", "in_50bp_from_LEJ": True},
    ]
    assert label_in_50bp_from_LEJ(cand_grna, ds) == expected
