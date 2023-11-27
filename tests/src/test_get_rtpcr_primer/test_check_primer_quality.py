from __future__ import annotations
import pytest
from src.get_rtpcr_primer.check_primer_quality import (
    check_exon_junction,
    rewrite_primer_pair_intron_len,
)


def test_true_check_exon_junction():
    True_primer_info = [
        # left primer is in exon junction
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 5,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
        # right primer is in exon junction
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 15,
            "right_primer_start": 49,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]
    exon_range = [[0, 10], [10, 50], [50, 100]]
    expected = [True, True]
    assert check_exon_junction(True_primer_info, exon_range) == expected


def test_false_check_exon_junction():
    False_primer_info = [
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 27,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
        {
            "primer_score": 999,
            "intron_len": 1000,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 10,
            "right_primer_start": 50,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]

    exon_range = [[0, 10], [10, 50], [50, 100]]

    expected = [False, False]
    assert check_exon_junction(False_primer_info, exon_range) == expected
