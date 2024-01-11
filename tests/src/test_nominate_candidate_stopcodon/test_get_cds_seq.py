from __future__ import annotations
import pytest
from dataclasses import dataclass
from create_gene_dataclass import set_dataclass
from src.nominate_candidate_stopcodon.get_cds_seq import (
    get_exon_seq,
    get_startcodon_exon_num,
    get_stopcodon_exon_num,
    get_cdsseq,
)

#####
# test dataset
#####
input_genedata = [
    {
        "geneName": "A",
        "name": "t1",
        "chrom": "chr1",
        "strand": "+",
        "txStart": "0",
        "txEnd": "100",
        "cdsStart": "10",
        "cdsEnd": "90",
        "exonCount": "3",
        "exonStarts": "0,25,60",
        "exonEnds": "15,50,95",
    },
    {
        "geneName": "B",
        "name": "t2",
        "chrom": "chr1",
        "strand": "+",
        "txStart": "0",
        "txEnd": "30",
        "cdsStart": "10",
        "cdsEnd": "25",
        "exonCount": "3",
        "exonStarts": "0,5,15",
        "exonEnds": "3,13,29",
    },
]
orf_seq_dict = {
    "t1::chr1:0-100": "1NNNNNNNNNATGNNNNNNNNNNNN2CAGNNNNNNNNNNNNNNNNNGGNNNNNNNNNNNN3NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNENNNNN",
    "t2::chr1:0-30": "NNNNNNNNNNATGTNANNNNNNNNNNNNNN",
}

test_name = ["t1", "t2"]

expected_return = [
    "1NNNNNNNNNATGNN2CAGNNNNNNNNNNNNNNNNNGGNN3NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNE",
    "NNNNNNNNATGANNNNNNNNNNNNN",  # NNN + NNNNNATG + NNNNNNNNNNNNNNN
]


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        expected_return,
    ),
)
def test_get_exon_seq(test_name, input_genedata, orf_seq_dict, expected):
    assert (
        get_exon_seq(set_dataclass(test_name, input_genedata, orf_seq_dict)) == expected
    )


#####
expected_return = [0, 1]


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        expected_return,
    ),
)
def test_get_startcodon_exon_num(test_name, input_genedata, orf_seq_dict, expected):
    assert (
        get_startcodon_exon_num(set_dataclass(test_name, input_genedata, orf_seq_dict))
        == expected
    )


#####
expected_return = [2, 2]


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        expected_return,
    ),
)
def test_get_stopcodon_exon_num(test_name, input_genedata, orf_seq_dict, expected):
    assert (
        get_stopcodon_exon_num(set_dataclass(test_name, input_genedata, orf_seq_dict))
        == expected
    )


#####
exon_seq = [
    "1NNNNNNNNNATGNN2CAGNNNNNNNNNNNNNNNNNGGNN3NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNE",
    "NNNNNNNNATGANNNNNNNNNNNNN",
]
startcodon_exonindex = [0, 1]
endcodon_exonindex = [2, 2]
expected_return = [
    "ATGNN2CAGNNNNNNNNNNNNNNNNNGGNN3NNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
    "ATGANNNNNNNNN",
]  # ATG +ANNNNNNNNN


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,exon_seq,start_index,end_index,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        exon_seq,
        startcodon_exonindex,
        endcodon_exonindex,
        expected_return,
    ),
)
def test_get_cdsseq(
    test_name, input_genedata, orf_seq_dict, exon_seq, start_index, end_index, expected
):
    assert (
        get_cdsseq(
            set_dataclass(test_name, input_genedata, orf_seq_dict),
            exon_seq,
            start_index,
            end_index,
        )
        == expected
    )
