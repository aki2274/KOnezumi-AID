from __future__ import annotations
import pytest
from dataclasses import dataclass
from get_rtpcr_primer.rate_quality import (
    check_exon_junction,
    rewrite_primer_pair_intron_len,
)


def test_cross_check_exon_junction():
    True_primer_info = [
        # left primer is in exon junction
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 29,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
        # right primer is in exon junction
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 45,
            "right_primer_start": 49,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]
    exon_range = [[0, 10], [10, 50], [50, 100]]
    expected = [
        {
            "left_cross_junction": True,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 29,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": True,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 45,
            "right_primer_start": 49,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]
    assert check_exon_junction(True_primer_info, exon_range) == expected


def test_not_cross_check_exon_junction():
    False_primer_info = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 45,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 1000,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 10,
            "right_primer_start": 50,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]

    exon_range = [[0, 10], [10, 50], [50, 100]]

    expected = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 45,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 1000,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 10,
            "right_primer_start": 50,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]
    assert check_exon_junction(False_primer_info, exon_range) == expected


@dataclass
class GeneData:
    orf_seq: str
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


# the inputdata is considered that only ~_exon_num.
def test_rewrite_primer_pair_intron_len():
    candidate_primer_info = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 5,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 1,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 5,
            "right_primer_start": 60,
            "left_primer_exon_num": 1,
            "right_primer_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 15,
            "right_primer_start": 49,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
        # 0 intron length
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 15,
            "right_primer_start": 49,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
    ]
    set_data = GeneData("sample_base", 0, 100, 10, 90, 3, [0, 25, 60], [15, 50, 95])
    expected = [
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 11,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 5,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 1,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 11,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 5,
            "right_primer_start": 60,
            "left_primer_exon_num": 1,
            "right_primer_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 46,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 15,
            "right_primer_start": 49,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
        {
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 15,
            "right_primer_start": 49,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
    ]
    result = rewrite_primer_pair_intron_len(candidate_primer_info, set_data)
    assert result == expected
