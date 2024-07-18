from __future__ import annotations
from src.konezumiaid.create_gene_dataclass import TranscriptRecord
from src.konezumiaid.evaluate_grna.apply_nmd_rules import (
    eliminate_in_back_half,
    eliminate_in_last_exon,
    label_in_start_150bp,
    label_in_50bp_from_LEJ,
    label_in_more_than_400bp_exon,
    create_candidates_list_dict,
)


def test_eliminate_in_back_half():
    cand_grna = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
    ]
    ds = TranscriptRecord(
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
    ds = TranscriptRecord(
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


def test_label_in_start_150bp():
    cand_grna = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
    ]
    ds = TranscriptRecord(
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


def test_label_in_50bp_from_LEJ():
    cand_grna = [
        {"ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
    ]
    ds = TranscriptRecord(
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


def test_label_in_more_than_400bp_exon():
    cand_grna = [
        {"ct_seq": "CAAACTATTTACTCTTTTCTTCG"},
        {"ct_seq": "TATTATAGGAAACTACAATGCAT"},
        {"ga_seq": "ACTACAATGCATTAAAGAACCTG"},
    ]
    transcript_record = TranscriptRecord(
        "AGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCACACAAACTATTTACTCTTTTCTTCGTAAGGAAAGGTTCAACTTCTGGTAAGTACAAAACTATTTGAAATATTTCTCTTAAAAATATGGTATGCTCATTTTATTTTTGGCAATATCTTTCATATATAGAACACAAACCTCTGTTTTATCTTTAACCATTAGAACAATACATGCTATTAATAGTTATTGTCTGTTACTAAAATCTAACAGCAATCTTTGTAGCATGCTGATATTGATATAGGGAAAAGCTTTAGTGTTGGAAAATTGTGATTCCAGTTTAGTGAAATATAAGGAGAGTTTACTGTAGTTTTGCTTATAATAGCATTATATACTTATTATAGGAAACTACAATGCATTAAAGAACCTGCTTATAATTAAAAGGTTATGGCTACCACAACAACTAGTATATTTTTCTGCTGAGTATCTTTTCCCTCC",
        0,
        491,
        0,
        491,
        4,
        [0, 3, 406, 480],
        [1, 404, 450, 491],
    )
    expected = [
        {"ct_seq": "CAAACTATTTACTCTTTTCTTCG", "in_more_than_400bp_exon": True},
        {"ct_seq": "TATTATAGGAAACTACAATGCAT", "in_more_than_400bp_exon": True},
        {"ga_seq": "ACTACAATGCATTAAAGAACCTG", "in_more_than_400bp_exon": False},
    ]
    assert label_in_more_than_400bp_exon(cand_grna, transcript_record) == expected


def test_create_candidates_list_dict():
    ct_cand = [
        {"aminoacid": "Q", "seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"aminoacid": "R", "seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"aminoacid": "S", "seq": "AGCAGCGGCGTTACTCGCTGTGG"},
        {"aminoacid": "C", "seq": "TGCAAGTGTTCAGCTTCCGCTGG"},
        {"aminoacid": "L", "seq": "TGCAGCTTGGGCAAATCTGGAGG"},
        {"aminoacid": "E", "seq": "GGCAGAGTGGGGAGAGCGGCAGG"},
    ]
    ga_cand = [
        {"aminoacid": "R", "seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
        {"aminoacid": "S", "seq": "CCTGCTCCTTTTGCATCTGGCTC"},
        {"aminoacid": "P", "seq": "CCTCCCTTGTGTCCTTGGCTTGG"},
        {"aminoacid": "L", "seq": "CCCTTGTGTCCTTGGCTTGGGCC"},
    ]

    expected = [
        {"aminoacid": "Q", "ct_seq": "CTGCAGAACTCGGACAATTCGGG"},
        {"aminoacid": "R", "ct_seq": "GCAGCGGCGTTACTCGCTGTGGG"},
        {"aminoacid": "S", "ct_seq": "AGCAGCGGCGTTACTCGCTGTGG"},
        {"aminoacid": "C", "ct_seq": "TGCAAGTGTTCAGCTTCCGCTGG"},
        {"aminoacid": "L", "ct_seq": "TGCAGCTTGGGCAAATCTGGAGG"},
        {"aminoacid": "E", "ct_seq": "GGCAGAGTGGGGAGAGCGGCAGG"},
        {"aminoacid": "R", "ga_seq": "CCTGCGTGGCCAGCGCTGGTGGT"},
        {"aminoacid": "S", "ga_seq": "CCTGCTCCTTTTGCATCTGGCTC"},
        {"aminoacid": "P", "ga_seq": "CCTCCCTTGTGTCCTTGGCTTGG"},
        {"aminoacid": "L", "ga_seq": "CCCTTGTGTCCTTGGCTTGGGCC"},
    ]

    assert create_candidates_list_dict(ct_cand, ga_cand) == expected
