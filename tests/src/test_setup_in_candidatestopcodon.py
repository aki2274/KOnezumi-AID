import pytest
import pandas as pd
from dataclasses import dataclass
from src.get_candidate_stopcodon_index import setup

# Make test data
# txStart とexonStart[0]は一致している必要がある
# test_df = pd.read_csv("test/data/test_genedata.csv")
gene_data = [
    [
        "A",
        "t1",
        "1",
        "+",
        "0",
        "30",
        "10",
        "25",
        "3",
        "0,5,15",
        "3,13,29",
    ]
]

col_names = [
    "geneName",
    "name",
    "chrom",
    "strand",
    "txStart",
    "txEnd",
    "cdsStart",
    "cdsEnd",
    "exonCount",
    "exonStarts",
    "exonEnds",
]
test_df = pd.DataFrame(gene_data, columns=col_names)

test_name = ["t1"]  # test_df["name"]

test_seq = [
    {"t1::1:0-100": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN"}
]  # dict key is "{transcripts_name}::{chrom}:{txStart}-{txEnd}". Value is orf_seq.


# setup module return DataClass
@dataclass
class expected_dataclass:
    orf_seq: str
    data: pd.DataFrame
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


expected_return = expected_dataclass(
    "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN",
    test_df[test_df["name"] == "t1"],
    0,
    30,
    10,
    25,
    3,
    [0, 5, 15],
    [3, 13, 29],
)
expected = [expected_return]


@pytest.mark.parametrize(
    "transcripts_name",
    "gene_df",
    "gene_seq_data",
    "expected",
    zip(test_name, [test_df], test_seq, expected),
)
def test_setup_data(transcripts_name, gene_df, gene_seq_data, expected):
    assert setup(transcripts_name, gene_df, gene_seq_data) == expected
