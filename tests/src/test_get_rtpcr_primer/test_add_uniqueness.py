import pytest
from src.get_rtpcr_primer.add_uniqueness import read_uniqueness, add_uniqueness


def test_read_uniqueness():
    expected = [
        ["23", "ATGC"],
        ["7", "TCGA"],
        ["7", "TTTT"],
        ["3", "AAAA"],
        ["12", "CCCC"],
    ]

    assert read_uniqueness(r"tests\data\uniq\test_counts.txt") == expected
