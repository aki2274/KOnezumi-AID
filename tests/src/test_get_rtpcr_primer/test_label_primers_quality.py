from __future__ import annotations
import pytest
from src.get_rtpcr_primer.label_primers_quality import label_primers_quality

candidate_primer_info = [
    [
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 27,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 10,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 27,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 1,
        },
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 27,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]
]

exon_junction_bools = [[True, False, False, False]]

expected_output = [
    [
        {
            "primer_score": 1,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 27,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 10,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 0,
        },
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 27,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 1,
        },
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_start": 27,
            "ringt_primer_start": 114,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },
    ]
]


@pytest.mark.parametrize(
    "candidate_primer_info,exon_junction_bool,expected",
    zip(candidate_primer_info, exon_junction_bools, expected_output),
)
def test_label_primers_quality(
    self, candidate_primer_info, exon_junction_bool, expected
):
    assert label_primers_quality(candidate_primer_info, exon_junction_bool) == expected
