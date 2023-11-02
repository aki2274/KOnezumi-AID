import pytest
import pandas as pd
from src.get_candidate_stopcodon_index import setup

# Make test data
# txStart とexonStart[0]は一致している必要がある

gene_data = [
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

test_name = ["t1"]

test_seq = {"t1::1:0-100": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN"}

expected = [
    {
        "orf_seq": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN",
        "data": test_df[test_df["name"] == "t1"],
        "txStart": 0,
        "txend": 30,
        "cdsStart": 10,
        "cdsEnd": 25,
        "exonCount": 3,
        "exon_start_list": [0, 5, 15],
        "exon_end_list": [3, 13, 29],
    }
]


@pytest.mark.parametrize(
    "transcripts_name,gene_df,gene_seq_data,expected",
    zip(test_name, test_df, test_seq, expected),
)
def test_get_exon_seq(transcripts_name, gene_df, gene_seq_data, expected):
    assert setup(transcripts_name, gene_df, gene_seq_data) == expected
