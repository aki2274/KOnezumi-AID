import os
import pytest
import pandas as pd
from dataclasses import dataclass
from src import set_gene_dataframe as sgd

# Make test data
# txStart とexonStart[0]は一致している必要がある
test_data_path = os.path.join(
    os.path.dirname(__file__), "..", "data", "test_genedata.csv"
)
test_df = pd.read_csv(test_data_path)


test_name = ["t1"]  # test_df["name"]

test_seq = [
    {"t1::chr1:0-30": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN"}
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
    test_df[test_df["name"] == "t1"].reset_index(),
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
    "transcripts_name,gene_df,gene_seq_data,expected",
    zip(test_name, [test_df], test_seq, expected),
)
def test_setup_data(transcripts_name, gene_df, gene_seq_data, expected):
    assert sgd.setup(transcripts_name, gene_df, gene_seq_data) == expected
