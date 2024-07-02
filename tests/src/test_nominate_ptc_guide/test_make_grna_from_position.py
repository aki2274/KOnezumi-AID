from konezumiaid.nominate_ptc_guide.make_grna_from_position import (
    TranscriptRecord,
    extract_c_to_t_grna_from_position,
    extract_g_to_a_grna_from_position,
)


def test_convert_ct_grna():
    ds = TranscriptRecord(
        transcript_seq="ATGNCAGNNNNNNNNNNNNNNNNPGGATCGANNNNCAGNNNNNNNNNNNNNNPGGCATATAA",
        transcript_start=0,
        transcript_end=80,
        cds_start=0,
        cds_end=80,
        exon_count=1,
        exon_start_positions=[0],
        exon_end_positions=[80],
    )
    indices = [4, 35]
    expected_grna = [
        "NCAGNNNNNNNNNNNNNNNNPGG",
        "NNNCAGNNNNNNNNNNNNNNPGG",
    ]
    assert extract_c_to_t_grna_from_position(ds, indices) == expected_grna


def test_convert_ga_grna():
    ds = TranscriptRecord(
        transcript_seq="ATGCCPNNNNNNNNNNNNNNTGGNNNATCGAGCCPNNNNNNNNNNNNNNNNNTGGCATATAA",
        transcript_start=0,
        transcript_end=80,
        cds_start=0,
        cds_end=80,
        exon_count=1,
        exon_start_positions=[0],
        exon_end_positions=[80],
    )
    indices = [20, 52]
    expected_grna = [
        "CCPNNNNNNNNNNNNNNTGGNNN",
        "CCPNNNNNNNNNNNNNNNNNTGG",
    ]
    assert extract_g_to_a_grna_from_position(ds, indices) == expected_grna


def test_nocandidate_convert_ct_grna():
    ds = TranscriptRecord(
        transcript_seq="ATGNCAGNNNNNNNNNNNNNNNNPGGATCGANNNNCAGNNNNNNNNNNNNNNPGGCATATAANNNNCAGNNNNNNNTTTTNNNNNPGGNNN",
        transcript_start=0,
        transcript_end=91,
        cds_start=0,
        cds_end=91,
        exon_count=1,
        exon_start_positions=[0],
        exon_end_positions=[91],
    )
    indices = [66]
    expected_grna = []
    assert extract_c_to_t_grna_from_position(ds, indices) == expected_grna


def test_nocandidate_convert_ga_grna():
    ds = TranscriptRecord(
        transcript_seq="ATGCCPNNNNNNNNNNNNNNTGGNNNATCGAGCCPNNNNNNNNNNNNNNNNNTGGCATATCCPNNNNNNNNNNNNNNAAAAGGAA",
        transcript_start=0,
        transcript_end=85,
        cds_start=0,
        cds_end=85,
        exon_count=1,
        exon_start_positions=[0],
        exon_end_positions=[85],
    )
    indices = [80]
    expected_grna = []
    assert extract_g_to_a_grna_from_position(ds, indices) == expected_grna
