from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def search_site_candidate(
    transcript_record: GeneData,
) -> tuple[list[dict[int, str], list[dict[int, str]]]]:
    acceptor_cands = [
        {
            "seq": transcript_record.orf_seq[start - 22 : start + 3][
                cc_idx : cc_idx + 23
            ],  # get 25bp sequence ,then extract 23bp sequence(PAM + 20bp) from the 25bp sequence
            "exon_index": i
            + 2,  # exon index is i+2 because the first exon is not included in the list and the index starts from 0
        }
        for i, start in enumerate(
            transcript_record.exon_start_list[1:]
        )  # eliminate the first exon
        if "CC"
        in transcript_record.orf_seq[start - 22 : start + 3][
            :4
        ]  # Check if the site is a candidate(wheather it has PAM site or not)
        and "AG"
        in transcript_record.orf_seq[
            start - 2 : start
        ]  # Check if the acceptor site consensus seq is present
        for cc_idx in [
            idx
            for idx in range(3)
            if transcript_record.orf_seq[start - 22 : start + 3][idx : idx + 2] == "CC"
        ]
    ]

    donor_cands = [
        {
            "seq": transcript_record.orf_seq[end - 21 : end + 4][
                cc_idx : cc_idx + 23
            ],  # get 25bp sequence ,then extract 23bp sequence(PAM + 20bp) from the 25bp sequence
            "exon_index": i + 1,  # exon index is i+1 because the index starts from 0
        }
        for i, end in enumerate(transcript_record.exon_end_list[:-1])
        if "CC"
        in transcript_record.orf_seq[end - 21 : end + 4][
            :4
        ]  # Check if the site is a candidate(wheather it has PAM site or not)
        and "GT"
        in transcript_record.orf_seq[
            end : end + 2
        ]  # Check if the donor site consensus seq is present
        for cc_idx in [
            idx
            for idx in range(3)
            if transcript_record.orf_seq[end - 21 : end + 4][idx : idx + 2] == "CC"
        ]
    ]

    acceptor_candidates = [
        candidate for candidate in acceptor_cands if "AAAA" not in candidate["seq"][3:]
    ]
    donor_candidates = [
        candidate for candidate in donor_cands if "AAAA" not in candidate["seq"][3:]
    ]

    return acceptor_candidates, donor_candidates
