from src.konezumiaid.format_output import format_single_transcript_result, extract_multiple_transcripts_match
from src.konezumiaid.create_gene_dataclass import TranscriptRecord


def test_format_single_transcript_result():
    candidate_multi_exon = [
        {
            "seq": "AAAA",
            "in_start_150bp": False,
            "in_50bp_from_LEJ": False,
            "aminoacid": "1M",
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
        {
            "seq": "TTTT",
            "in_start_150bp": True,
            "in_50bp_from_LEJ": False,
            "aminoacid": "12L",
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
        {
            "seq": "CCCC",
            "in_start_150bp": False,
            "in_50bp_from_LEJ": True,
            "aminoacid": "32P",
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
        {
            "seq": "GGGG",
            "in_start_150bp": True,
            "in_50bp_from_LEJ": True,
            "aminoacid": "451A",
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
    ]
    candidate_single_exon = [
        {
            "seq": "AAAA",
            "in_start_150bp": False,
            "aminoacid": "123M",
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
        {
            "seq": "TTTT",
            "in_start_150bp": True,
            "aminoacid": "12L",
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
    ]
    candidate_splice_site = [
        {
            "seq": "AAAA",
            "exon_index": 1,
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
        {
            "seq": "TTTT",
            "exon_index": 2,
            "link_to_crisprdirect": "https://crisprdirect.org",
        },
    ]
    transcript_record_multi_exon = TranscriptRecord(
        "ATGCGTACG",
        100,
        200,
        120,
        180,
        2,
        [100, 150],
        [120, 200],
    )
    transcript_record_single_sxon = TranscriptRecord(
        "ATGCGTACG",
        100,
        200,
        120,
        180,
        1,
        [100, 150],
        [120, 200],
    )
    expect_multi_exon = [
        {
            "Target sequence (20mer + PAM)": "AAAA",
            "Recommended": True,
            "Target amino acid": "1M",
            "link to CRISPRdirect": "https://crisprdirect.org",
        },
        {
            "Target sequence (20mer + PAM)": "TTTT",
            "Recommended": False,
            "Target amino acid": "12L",
            "link to CRISPRdirect": "https://crisprdirect.org",
        },
        {
            "Target sequence (20mer + PAM)": "CCCC",
            "Recommended": False,
            "Target amino acid": "32P",
            "link to CRISPRdirect": "https://crisprdirect.org",
        },
        {
            "Target sequence (20mer + PAM)": "GGGG",
            "Recommended": False,
            "Target amino acid": "451A",
            "link to CRISPRdirect": "https://crisprdirect.org",
        },
    ]

    except_single_exon = [
        {
            "Target sequence (20mer + PAM)": "AAAA",
            "Recommended": True,
            "Target amino acid": "123M",
            "link to CRISPRdirect": "https://crisprdirect.org",
        },
        {
            "Target sequence (20mer + PAM)": "TTTT",
            "Recommended": False,
            "Target amino acid": "12L",
            "link to CRISPRdirect": "https://crisprdirect.org",
        },
    ]

    expect_splice_site = [
        {"Target sequence (20mer + PAM)": "AAAA", "Exon index": 1, "link to CRISPRdirect": "https://crisprdirect.org"},
        {"Target sequence (20mer + PAM)": "TTTT", "Exon index": 2, "link to CRISPRdirect": "https://crisprdirect.org"},
    ]

    print(format_single_transcript_result(candidate_multi_exon, transcript_record_multi_exon, flag_ptc=True))

    assert (
        format_single_transcript_result(candidate_multi_exon, transcript_record_multi_exon, flag_ptc=True).to_dict(
            orient="records"
        )
        == expect_multi_exon
    )
    assert (
        format_single_transcript_result(candidate_single_exon, transcript_record_single_sxon, flag_ptc=True).to_dict(
            orient="records"
        )
        == except_single_exon
    )
    assert (
        format_single_transcript_result(candidate_splice_site, transcript_record_multi_exon).to_dict(orient="records")
        == expect_splice_site
    )
