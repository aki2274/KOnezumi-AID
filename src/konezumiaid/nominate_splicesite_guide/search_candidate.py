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

    index_exon_has_3_utr = next(
        (
            i
            for i, (start, end) in enumerate(
                zip(transcript_record.exon_start_list, transcript_record.exon_end_list)
            )
            if start < transcript_record.cdsEnd <= end
        )
    )

    # remove candidates that have "AAAA" or the exon length is a multiple of 3 or the exon has only 3'UTR
    acceptor_candidates = [
        cand
        for cand in acceptor_cands
        if "AAAA" not in cand["seq"][3:]  # exclude PAM
        and (
            transcript_record.exon_end_list[cand["exon_index"] - 1]
            - transcript_record.exon_start_list[cand["exon_index"] - 1]
        )
        % 3
        != 0
        and cand["exon_index"] <= index_exon_has_3_utr
    ]

    donor_candidates = [
        cand
        for cand in donor_cands
        if "AAAA" not in cand["seq"][3:]  # exclude PAM
        and (
            transcript_record.exon_end_list[cand["exon_index"] - 1]
            - transcript_record.exon_start_list[cand["exon_index"] - 1]
        )
        % 3
        != 0
        and cand["exon_index"] < index_exon_has_3_utr
    ]

    return acceptor_candidates, donor_candidates
