from __future__ import annotations
from src.create_gene_dataclass import GeneData
from src.adjust_nmd_rules.rete_position import (
    # label_in_front_half,
    eliminate_in_last_exon,
    label_in_start_150bp,
    label_in_50bp_from_LEJ,
)


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
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG", "start_150": True},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG", "start_150": True},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT", "start_150": False},
    ]
    assert label_in_start_150bp(cand_grna, ds) == expected


"""
def test_label_in_front_half():
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
        1,
        [0],
        [91],
    )
    expected = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG", "back_half": True},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG", "back_half": False},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT", "back_half": False},
    ]
    assert label_in_front_half(cand_grna, ds) == expected
"""


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
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG", "50bp_from_LEJ": True},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG", "50bp_from_LEJ": False},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT", "50bp_from_LEJ": False},
    ]
    assert label_in_50bp_from_LEJ(cand_grna, ds) == expected
