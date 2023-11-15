from __future__ import annotations
import os
import pytest
from dataclasses import dataclass, asdict
from src.import_dataset import read_csv
from src.set_gene_dataclass import set_dataclass

# Make test data
# txStart とexonStart[0]は一致している必要がある
test_data_path = os.path.join(
    os.path.dirname(__file__), "..", "data", "test_genedata.csv"
)
input_genedata = read_csv(test_data_path)


input_transcript_name = ["t1"]  # test_df["name"]

input_seq = [
    {"t1::chr1:0-30": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN"}
]  # dict key is "{transcripts_name}::{chrom}:{txStart}-{txEnd}". Value is orf_seq.


# setup module return DataClass
@dataclass
class DataClass:
    orf_seq: str
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


expected_return = DataClass(
    "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN",
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
    "transcripts_name,gene_data,gene_seq_data,expected",
    zip(input_transcript_name, [input_genedata], input_seq, expected),
)
def test_setup_data(transcripts_name, gene_data, gene_seq_data, expected):
    assert asdict(set_dataclass(transcripts_name, gene_data, gene_seq_data)) == asdict(
        expected
    )