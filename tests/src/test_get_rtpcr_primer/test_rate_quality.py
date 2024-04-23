from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.get_rtpcr_primer.rate_quality import (
    is_crossing_juncion,
    add_intron_len,
)


def test_is_cross_junction():
    True_info = [
        # left primer is in exon junction
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 29,
            "right_start": 60,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
        # right primer is in exon junction
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 45,
            "right_start": 49,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
    ]
    exon_range = [[0, 10], [10, 50], [50, 100]]
    expected = [
        {
            "left_cross_junction": True,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 29,
            "right_start": 60,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": True,
            "intron_len": 0,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 45,
            "right_start": 49,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
    ]
    assert is_crossing_juncion(True_info, exon_range) == expected


def test_not_corss_verify_cross_junction():
    False_info = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 45,
            "right_start": 60,
            "left_exon_num": 0,
            "right_exon_num": 0,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 1000,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 30,
            "right_start": 50,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
    ]

    exon_range = [[0, 10], [10, 50], [50, 100]]

    expected = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 45,
            "right_start": 60,
            "left_exon_num": 0,
            "right_exon_num": 0,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 1000,
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_end": 30,
            "right_start": 50,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
    ]
    assert is_crossing_juncion(False_info, exon_range) == expected


# the inputdata is considered that only ~_exon_num.
def test_add_intron_len():
    candidate_info = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 5,
            "right_start": 60,
            "left_exon_num": 0,
            "right_exon_num": 1,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 5,
            "right_start": 60,
            "left_exon_num": 1,
            "right_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 15,
            "right_start": 49,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
        # 0 intron length
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 15,
            "right_start": 49,
            "left_exon_num": 0,
            "right_exon_num": 0,
        },
    ]
    set_data = GeneData("sample_base", 0, 100, 10, 90, 3, [0, 25, 60], [15, 50, 95])
    expected = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 11,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 5,
            "right_start": 60,
            "left_exon_num": 0,
            "right_exon_num": 1,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 11,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 5,
            "right_start": 60,
            "left_exon_num": 1,
            "right_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 46,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 15,
            "right_start": 49,
            "left_exon_num": 0,
            "right_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 15,
            "right_start": 49,
            "left_exon_num": 0,
            "right_exon_num": 0,
        },
    ]
    result = add_intron_len(candidate_info, set_data)
    assert result == expected
