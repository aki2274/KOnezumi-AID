from __future__ import annotations
from src.adjust_nmd_rules.create_grna_lsit import create_list


def test_create_list():
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

    assert create_list(ct_cand, ga_cand) == expected
