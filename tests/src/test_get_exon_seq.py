import os
import pytest
import numpy as np
import pandas as pd
from dataclasses import dataclass
from src.set_gene_dataclass import set_dataclass
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
test_seq = {"t1::chr1:0-30": "NNNNNNNNNNATGTNNNNNNNNNNNNNNNN"}
test_gene_dataclass = set_dataclass("t1", test_df, test_seq)

#####
expected_return = ["NNNNNNNNATGNNNNNNNNNNNNNN"]  # NNN + NNNNNATG + NNNNNNNNNNNNNN


@pytest.mark.parametrize(
    "test_dataclass,expected",
    zip([test_gene_dataclass], expected_return),
)
def test_get_exon_seq(test_dataclass, expected):
    assert get_exon_seq(test_dataclass) == expected


#####
expected_return = [1]


@pytest.mark.parametrize(
    "test_dataclass,expected",
    zip([test_gene_dataclass], expected_return),
)
def test_get_startcodon_exonindex(test_dataclass, expected):
    assert get_startcodon_exonindex(test_dataclass) == expected


#####
expected_return = [2]


@pytest.mark.parametrize(
    "test_dataclass,expected", zip([test_gene_dataclass], expected_return)
)
def test_get_stopcodon_exonindex(test_dataclass, expected):
    assert get_stopcodon_exonindex(test_dataclass) == expected


#####
exon_seq = ["NNNNNNNNATGNNNNNNNNNNNNNN"]
startcodon_exonindex = [1]
endcodon_exonindex = [2]
expected_return = ["ATGNNNNNNNNNNNN"]

"""
@pytest.mark.parametrize(
    "test_dataclass,exon_seq,start_index,end_index,expected",
    zip(
        [test_gene_dataclass],
        exon_seq,
        startcodon_exonindex,
        endcodon_exonindex,
        expected_return,
    ),
)
def test_get_cdsseq(test_dataclass, exon_seq, start_index, end_index, expected):
    assert get_cdsseq(test_dataclass, exon_seq, start_index, end_index) == expected
"""
