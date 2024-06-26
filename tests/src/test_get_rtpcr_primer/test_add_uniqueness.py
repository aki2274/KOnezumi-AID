from __future__ import annotations
from konezumiaid.get_rtpcr_primer.add_uniqueness import read_uniqueness, add_uniqueness


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
        {"left": "ATGC", "right": "TCGA"},
        {"left": "TTTT", "right": "AAAA"},
        {"left": "CCCC", "right": "GGGG"},
    ]

    expected = [
        {
            "left": "ATGC",
            "right": "TCGA",
            "left_0_mismatch": "23",
            "right_0_mismatch": "7",
            "left_1_mismatch": "23",
            "right_1_mismatch": "7",
            "left_2_mismatch": "23",
            "right_2_mismatch": "7",
        },
        {
            "left": "TTTT",
            "right": "AAAA",
            "left_0_mismatch": "7",
            "right_0_mismatch": "3",
            "left_1_mismatch": "7",
            "right_1_mismatch": "3",
            "left_2_mismatch": "7",
            "right_2_mismatch": "3",
        },
        {
            "left": "CCCC",
            "right": "GGGG",
            "left_0_mismatch": "12",
            "right_0_mismatch": 0,
            "left_1_mismatch": "12",
            "right_1_mismatch": 0,
            "left_2_mismatch": "12",
            "right_2_mismatch": 0,
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
