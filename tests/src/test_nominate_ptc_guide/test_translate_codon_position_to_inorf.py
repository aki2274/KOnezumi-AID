from __future__ import annotations
import pytest
from src.konezumiaid.create_gene_dataclass import create_dataclass
from src.konezumiaid.get_exon_range import get_exon_range
from src.konezumiaid.nominate_ptc_guide.translate_codon_position_to_inorf import (
    translate_cds_position_to_exon,
    get_exonindex_in_cand_codon,
    translate_position_in_splicedexon_to_orf,
)

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
    "t1": "1NNNNNNNNNATGNNNNNNNNNNNN2CAGNNNNNNNNNNNNNNNNNGGNNNNNNNNNNNN3NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNENNNNN",
    "t2": "NNNNNNNNNNATGTNANNNNNNNNNNNNNN",
}
test_name = ["t1", "t2"]
expected = [[0, 14], [15, 39], [40, 74]], [[0, 2], [3, 10], [11, 24]]


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        expected,
    ),
)
def test_get_exon_range(test_name, input_genedata, orf_seq_dict, expected):
    assert get_exon_range(create_dataclass(test_name, input_genedata, orf_seq_dict)) == expected


# there is the candidate stopcodon
test_name = ["t1"]
candidate_stopcodon_index_incds = [[6]]  # len(ATGNN2CAG)-3. length of C in cds seq
exon_range = [[[0, 14], [15, 39], [40, 74]]]
exon_num = [0]
expected = [[16]]  # len(NNNNNNNNNNATGNN2CAG)-3. length of C in exon seq


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,candidate_stopcodon_index_incds,exon_range,exon_num,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        candidate_stopcodon_index_incds,
        exon_range,
        exon_num,
        expected,
    ),
)
def test_translate_cds_position_to_exon(
    test_name,
    input_genedata,
    orf_seq_dict,
    candidate_stopcodon_index_incds,
    exon_range,
    exon_num,
    expected,
):
    assert (
        translate_cds_position_to_exon(
            create_dataclass(test_name, input_genedata, orf_seq_dict),
            candidate_stopcodon_index_incds,
            exon_num,
            exon_range,
        )
        == expected
    )


candidate_stopcodon_index_inexon = [[16]]
expected = [[1]]


@pytest.mark.parametrize(
    "candidate_stopcodon,exon_range,expected",
    zip(
        candidate_stopcodon_index_inexon,
        exon_range,
        expected,
    ),
)
def test_get_exonindex_in_cand_codon(candidate_stopcodon, exon_range, expected):
    assert get_exonindex_in_cand_codon(candidate_stopcodon, exon_range) == expected


exon_index = [[1]]
exprcted = [[26]]


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,candidate_stopcodon_index_incds,exon_num,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        candidate_stopcodon_index_incds,
        exon_num,
        exprcted,
    ),
)
def test_translate_position_in_splicedexon_to_orf(
    test_name,
    input_genedata,
    orf_seq_dict,
    candidate_stopcodon_index_incds,
    exon_num,
    expected,
):
    assert (
        translate_position_in_splicedexon_to_orf(
            create_dataclass(test_name, input_genedata, orf_seq_dict),
            candidate_stopcodon_index_incds,
            exon_num,
        )
        == expected
    )


# there is no candidate stopcodon
test_name = ["t2"]
candidate_stopcodon_index_incds = [[]]
exon_range = [[[0, 2], [3, 10], [11, 24]]]
exon_num = [1]
expected = [[]]


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,candidate_stopcodon,exon_range,exon_num,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        candidate_stopcodon_index_incds,
        exon_range,
        exon_num,
        expected,
    ),
)
def test_non_translate_cds_position_to_exon(
    test_name,
    input_genedata,
    orf_seq_dict,
    candidate_stopcodon,
    exon_range,
    exon_num,
    expected,
):
    assert (
        translate_cds_position_to_exon(
            create_dataclass(test_name, input_genedata, orf_seq_dict),
            candidate_stopcodon,
            exon_num,
            exon_range,
        )
        == expected
    )


candidate_stopcodon_index_inexon = [[]]
expected = [[]]


@pytest.mark.parametrize(
    "candidate_stopcodon,exon_range,expected",
    zip(
        candidate_stopcodon_index_inexon,
        exon_range,
        expected,
    ),
)
def test_non_get_exonindex_in_cand_codon(candidate_stopcodon, exon_range, expected):
    assert get_exonindex_in_cand_codon(candidate_stopcodon, exon_range) == expected


exon_index = [[]]
exprcted = [[]]


@pytest.mark.parametrize(
    "test_name,input_genedata,orf_seq_dict,candidate_stopcodon_index_incds,exon_num,expected",
    zip(
        test_name,
        [input_genedata] * len(test_name),
        [orf_seq_dict] * len(test_name),
        candidate_stopcodon_index_incds,
        exon_num,
        exprcted,
    ),
)
def test_non_translate_position_in_splicedexon_to_orf(
    test_name,
    input_genedata,
    orf_seq_dict,
    candidate_stopcodon_index_incds,
    exon_num,
    expected,
):
    assert (
        translate_position_in_splicedexon_to_orf(
            create_dataclass(test_name, input_genedata, orf_seq_dict),
            candidate_stopcodon_index_incds,
            exon_num,
        )
        == expected
    )
