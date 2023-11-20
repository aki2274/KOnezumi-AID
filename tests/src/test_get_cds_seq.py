from __future__ import annotations
import os
import pytest
import numpy as np
from dataclasses import dataclass
from src.import_dataset import read_pkl
from src.set_gene_dataclass import set_dataclass
from src.get_candidate_stopcodon_index import (
    get_exon_seq,
    get_startcodon_exon_num,
    get_stopcodon_exon_num,
    get_cdsseq,
)

#####
# test dataset
#####
test_data_path = os.path.join(
    os.path.dirname(__file__), "..", "data", "test_genedata.pkl"
)
test_list_dict = read_pkl(test_data_path)
orf_seq_dict = {"t1::chr1:0-30": "NNNNNNNNNNATGTNANNNNNNNNNNNNNN"}
test_gene_dataclass = set_dataclass("t1", test_list_dict, orf_seq_dict)

#####
expected_return = ["NNNNNNNNATGANNNNNNNNNNNNN"]  # NNN + NNNNNATG + NNNNNNNNNNNNNNN


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
def test_get_startcodon_exon_num(test_dataclass, expected):
    assert get_startcodon_exon_num(test_dataclass) == expected


#####
expected_return = [2]


@pytest.mark.parametrize(
    "test_dataclass,expected", zip([test_gene_dataclass], expected_return)
)
def test_get_stopcodon_exon_num(test_dataclass, expected):
    assert get_stopcodon_exon_num(test_dataclass) == expected


#####
exon_seq = ["NNNNNNNNATGANNNNNNNNNNNNN"]
startcodon_exonindex = [1]
endcodon_exonindex = [2]
expected_return = ["ATGANNNNNNNNN"]  # ATG +ANNNNNNNNN


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
