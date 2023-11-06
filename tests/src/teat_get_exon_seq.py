import os
import pytest
import numpy as np
import pandas as pd
from dataclasses import dataclass
from src.set_gene_dataframe import setup
from src.get_candidate_stopcodon_index import (
    get_exon_seq,
    get_startcodon_exonindex,
    get_stopcodon_exonindex,
    get_cdsseq,
)

#####
# test dataset
#####
test_data_path = os.path.join(
    os.path.dirname(__file__), "..", "data", "test_genedata.csv"
)
test_df = pd.read_csv(test_data_path)


@dataclass
class test_dataclass:
    orf_seq: str
    data: pd.DataFrame
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


test_seq = {"t1::chr1:0-30": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN"}

test_dataset = setup("t1", test_df, test_seq)

expected_return = ["ATGTNNNNNNNNNNN"]


def test_get_exon_seq(test_dataset):
    assert get_exon_seq(test_dataset) == expected_return
