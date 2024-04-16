from konezumiaid.nominate_candidate_stopcodon.make_grna_from_index import (
    GeneData,
    convert_ct_grna,
    convert_ga_grna,
)


def test_convert_ct_grna():
    ds = GeneData(
        orf_seq="ATGNCAGNNNNNNNNNNNNNNNNPGGATCGANNNNCAGNNNNNNNNNNNNNNPGGCATATAA",
        txStart=0,
        txend=80,
        cdsStart=0,
        cdsEnd=80,
        exonCount=1,
        exon_start_list=[0],
        exon_end_list=[80],
    )
    indices = [4, 35]
    expected_grna = [
        "NCAGNNNNNNNNNNNNNNNNPGG",
        "NNNCAGNNNNNNNNNNNNNNPGG",
    ]
    assert convert_ct_grna(ds, indices) == expected_grna


def test_convert_ga_grna():
    ds = GeneData(
        orf_seq="ATGCCPNNNNNNNNNNNNNNTGGNNNATCGAGCCPNNNNNNNNNNNNNNNNNTGGCATATAA",
        txStart=0,
        txend=80,
        cdsStart=0,
        cdsEnd=80,
        exonCount=1,
        exon_start_list=[0],
        exon_end_list=[80],
    )
    indices = [20, 52]
    expected_grna = [
        "CCPNNNNNNNNNNNNNNTGGNNN",
        "CCPNNNNNNNNNNNNNNNNNTGG",
    ]
    assert convert_ga_grna(ds, indices) == expected_grna


def test_nocandidate_convert_ct_grna():
    ds = GeneData(
        orf_seq="ATGNCAGNNNNNNNNNNNNNNNNPGGATCGANNNNCAGNNNNNNNNNNNNNNPGGCATATAANNNNCAGNNNNNNNTTTTNNNNNPGGNNN",
        txStart=0,
        txend=91,
        cdsStart=0,
        cdsEnd=91,
        exonCount=1,
        exon_start_list=[0],
        exon_end_list=[91],
    )
    indices = [66]
    expected_grna = []
    assert convert_ct_grna(ds, indices) == expected_grna


def test_nocandidate_convert_ga_grna():
    ds = GeneData(
        orf_seq="ATGCCPNNNNNNNNNNNNNNTGGNNNATCGAGCCPNNNNNNNNNNNNNNNNNTGGCATATCCPNNNNNNNNNNNNNNTTTTGGAA",
        txStart=0,
        txend=85,
        cdsStart=0,
        cdsEnd=85,
        exonCount=1,
        exon_start_list=[0],
        exon_end_list=[85],
    )
    indices = [80]
    expected_grna = []
    assert convert_ga_grna(ds, indices) == expected_grna
