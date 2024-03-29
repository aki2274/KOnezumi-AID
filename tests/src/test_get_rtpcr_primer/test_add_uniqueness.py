from __future__ import annotations
from src.get_rtpcr_primer.add_uniqueness import read_uniqueness, add_uniqueness


def test_read_uniqueness():
    expected = [
        ["23", "ATGC"],
        ["7", "TCGA"],
        ["7", "TTTT"],
        ["3", "AAAA"],
        ["12", "CCCC"],
    ]

    assert read_uniqueness("tests/data/uniq/test_counts.txt") == expected


def test_add_uniqueness():
    candidate = [
        {"left_primer": "ATGC", "right_primer": "TCGA"},
        {"left_primer": "TTTT", "right_primer": "AAAA"},
        {"left_primer": "CCCC", "right_primer": "GGGG"},
    ]
    expected = [
        {
            "left_primer": "ATGC",
            "right_primer": "TCGA",
            "left_0_uniq": "23",
            "right_0_uniq": "7",
            "left_1_uniq": "23",
            "right_1_uniq": "7",
            "left_2_uniq": "23",
            "right_2_uniq": "7",
        },
        {
            "left_primer": "TTTT",
            "right_primer": "AAAA",
            "left_0_uniq": "7",
            "right_0_uniq": "3",
            "left_1_uniq": "7",
            "right_1_uniq": "3",
            "left_2_uniq": "7",
            "right_2_uniq": "3",
        },
        {
            "left_primer": "CCCC",
            "right_primer": "GGGG",
            "left_0_uniq": "12",
            "right_0_uniq": 0,
            "left_1_uniq": "12",
            "right_1_uniq": 0,
            "left_2_uniq": "12",
            "right_2_uniq": 0,
        },
    ]
    assert (
        add_uniqueness(
            candidate,
            "tests/data/uniq/test_counts.txt",
            "tests/data/uniq/test_counts.txt",
            "tests/data/uniq/test_counts.txt",
        )
        == expected
    )
