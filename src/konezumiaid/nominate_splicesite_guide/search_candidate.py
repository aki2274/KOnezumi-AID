from __future__ import annotations
from konezumiaid.get_reverse_complement import get_revcomp
from konezumiaid.create_gene_dataclass import GeneData


def search_site_candidate(
    transcript_record: GeneData,
) -> tuple[list[dict[int, str], list[dict[int, str]]]]:
    """
    Search gRNA candidate for targeting the splice site. The candidate must have PAM site and the consensus sequence of the splice site.
    Exclude candidates that have "AAAA" or the exon length is a multiple of 3 or the exon has 3'UTR.
    ( If the target site is acceptor site, It is permissible to set the target exon as having a 3'UTR. )
    """
    orf = transcript_record.orf_seq
    acceptor_cands = [
        {
            "seq": get_revcomp(orf[start - 22 : start + 3][cc_idx : cc_idx + 23]),
            # get 25bp sequence ,then extract 23bp sequence(PAM + 20bp) from the 25bp sequence
            "exon_index": i + 2,
            # exon index is i+2 because the first exon is not included in the list and the index starts from 0
        }
        for i, start in enumerate(transcript_record.exon_start_list[1:])
        # eliminate the first exon
        if "CC" in orf[start - 22 : start + 3][:4] and "AG" in orf[start - 2 : start]
        # Check if the site is a candidate(wheather it has PAM site or not) the acceptor site consensus seq is present
        for cc_idx in [
            idx
            for idx in range(3)
            if orf[start - 22 : start + 3][idx : idx + 2] == "CC"
        ]
    ]

    donor_cands = [
        {
            "seq": get_revcomp(orf[end - 21 : end + 4][cc_idx : cc_idx + 23]),
            # get 25bp sequence ,then extract 23bp sequence(PAM + 20bp) from the 25bp sequence
            "exon_index": i + 1,
            # exon index is i+1 because the index starts from 0
        }
        for i, end in enumerate(transcript_record.exon_end_list[:-1])
        if "CC" in orf[end - 21 : end + 4][:4] and "GT" in orf[end : end + 2]
        # Check if the site is a candidate(wheather it has PAM site or not) the donor site consensus seq is present
        for cc_idx in [
            idx for idx in range(3) if orf[end - 21 : end + 4][idx : idx + 2] == "CC"
        ]
    ]
    exon_start_pos = transcript_record.exon_start_list
    exon_end_pos = transcript_record.exon_end_list

    index_exon_has_3_utr = next(
        (
            i
            for i, (start, end) in enumerate(zip(exon_start_pos, exon_end_pos))
            if start < transcript_record.cdsEnd <= end
        )
    )

    # remove candidates that have "AAAA" or the exon length is a multiple of 3 or the exon has only 3'UTR
    acceptor_candidates = [
        cand
        for cand in acceptor_cands
        if "TTTT" not in cand["seq"][:20]  # exclude PAM
        and (
            exon_end_pos[cand["exon_index"] - 1]
            - exon_start_pos[cand["exon_index"] - 1]
        )
        % 3
        != 0
        and (cand["exon_index"] - 1) <= index_exon_has_3_utr
        # It is permissible to set the target exon as having a 3'UTR.
    ]

    donor_candidates = [
        cand
        for cand in donor_cands
        if "TTTT" not in cand["seq"][:20]  # exclude PAM
        and (
            transcript_record.exon_end_list[cand["exon_index"] - 1]
            - exon_start_pos[cand["exon_index"] - 1]
        )
        % 3
        != 0
        and (cand["exon_index"] - 1) < index_exon_has_3_utr
        # It is NOT acceptable.exon as having a 3'UTR.
    ]

    return acceptor_candidates, donor_candidates
