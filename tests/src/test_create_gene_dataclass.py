from __future__ import annotations
import pytest
from dataclasses import dataclass, asdict
from src.konezumiaid.create_gene_dataclass import create_dataclass

# Make test data
# txStart とexonStart[0]は一致している必要がある

input_genedata = [
    {
        "geneName": "A",
        "name": "t1",
        "chrom": "chr1",
        "strand": "+",
        "txStart": "0",
        "txEnd": "30",
        "cdsStart": "10",
        "cdsEnd": "25",
        "exonCount": "3",
        "exonStarts": "0,5,15",
        "exonEnds": "3,13,29",
    }
]
input_transcript_name = ["t1"]  # test_df["name"]

input_seq = [
    {"t1": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN"}
]  # dict key is "{transcripts_name}::{chrom}:{txStart}-{txEnd}". Value is orf_seq.


# setup module return DataClass
@dataclass
class DataClass:
    transcript_seq: str
    transcript_start: int
    transcript_end: int
    cds_start: int
    cds_end: int
    exon_count: int
    exon_start_positions: list[int]
    exon_end_positions: list[int]


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
def test_create_dataclass(transcripts_name, gene_data, gene_seq_data, expected):
    assert asdict(create_dataclass(transcripts_name, gene_data, gene_seq_data)) == asdict(expected)
