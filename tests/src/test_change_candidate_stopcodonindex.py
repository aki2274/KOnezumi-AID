from __future__ import annotations
import os
import pytest
from src.set_gene_dataclass import set_dataclass
from src.import_dataset import read_csv
from src.get_candidate_stopcodon_index import get_startcodon_exonindex
from src.change_exon_index_to_gene import (
    get_range_of_exon,
    get_candidate_codon_exon_index,
    add_num_to_orf_index,
)

test_data_path = os.path.join(
    os.path.dirname(__file__), "..", "data", "test_genedata.csv"
)
test_df = read_csv(test_data_path)
test_seq = {"t1::chr1:0-30": "NNNNNNNNNNATGTNANNNNNNNNNNNNNN"}
test_gene_dataclass = set_dataclass("t1", test_df, test_seq)

test_candidate_codon_exon_index = [get_startcodon_exonindex(test_gene_dataclass)]

expected = [[[0, 8], [8, 22]]]


@pytest.mark.parametrize(
    "test_dataclass,test_cdsStart_index,expected",
    zip([test_gene_dataclass], test_candidate_codon_exon_index, expected),
)
def test_get_range_of_exon(test_dataclass, test_cdsStart_index, expected):
    assert get_range_of_exon(test_dataclass, test_cdsStart_index) == expected
