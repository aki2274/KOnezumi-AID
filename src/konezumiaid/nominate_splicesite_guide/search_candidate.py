from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def search_site_candidate(
    transcript_record: GeneData,
) -> tuple[list[dict[int, str], list[dict[int, str]]]]:
    acceptor_cands = [
        {
            "seq": transcript_record.orf_seq[start - 22 : start + 3][
                cc_idx : cc_idx + 23
            ],
            "exon_index": i + 2,
        }
        for i, start in enumerate(transcript_record.exon_start_list[1:])
        if "CC" in transcript_record.orf_seq[start - 22 : start + 3][:4]
        and "AG" in transcript_record.orf_seq[start - 2 : start]
        for cc_idx in [
            idx
            for idx in range(3)
            if transcript_record.orf_seq[start - 22 : start + 3][idx : idx + 2] == "CC"
        ]
    ]

    donor_cands = [
        {
            "seq": transcript_record.orf_seq[end - 21 : end + 4][cc_idx : cc_idx + 23],
            "exon_index": i + 1,
        }
        for i, end in enumerate(transcript_record.exon_end_list[:-1])
        if "CC" in transcript_record.orf_seq[end - 21 : end + 4][:4]
        and "GT" in transcript_record.orf_seq[end : end + 2]
        for cc_idx in [
            idx
            for idx in range(3)
            if transcript_record.orf_seq[end - 21 : end + 4][idx : idx + 2] == "CC"
        ]
    ]

    acceptor_candidates = [
        candidate for candidate in acceptor_cands if "TTTT" not in candidate["seq"]
    ]
    donor_candidates = [
        candidate for candidate in donor_cands if "TTTT" not in candidate["seq"]
    ]

    return acceptor_candidates, donor_candidates
