from __future__ import annotations
from konezumiaid.get_reverse_complement import get_revcomp
from konezumiaid.create_gene_dataclass import TranscriptRecord
from konezumiaid.evaluate_grna.add_grna_info import link_to_crisperdirect


def find_splice_site_candidate(
    orf: str, positions: list[int], offset: int, target_seq: str, index_adjustment: int, acc_flag: bool = True
) -> list[dict[int, str]]:
    search_length = 25
    seq_length = 23
    candidates = []

    for i, pos in enumerate(positions):
        orf_segment = orf[pos - offset : pos - offset + search_length]
        target_segment = orf[pos - 2 : pos] if acc_flag else orf[pos : pos + 2]

        if "CC" in orf_segment[:4] and target_seq in target_segment:
            for cc_idx in [idx for idx in range(3) if orf_segment[idx : idx + 2] == "CC"]:
                candidate_seq = get_revcomp(orf_segment[cc_idx : cc_idx + seq_length])
                candidates.append(
                    {
                        "seq": candidate_seq,
                        "exon_index": i + index_adjustment,
                    }
                )

    return candidates


def filter_candidate(
    candidates: list[dict[int, str]],
    index_exon_has_3_utr: int,
    start_pos: list[int],
    end_pos: list[int],
    acc_flag: bool = True,
) -> list[dict[int, str]]:
    return [
        cand
        for cand in candidates
        if "TTTT" not in cand["seq"][:20]  # exclude PAM
        and (end_pos[cand["exon_index"] - 1] - start_pos[cand["exon_index"] - 1]) % 3 != 0
        and (
            (cand["exon_index"] - 1) <= index_exon_has_3_utr
            if acc_flag
            else (cand["exon_index"] - 1) < index_exon_has_3_utr
        )
    ]


def search_site_candidate(
    transcript_record: TranscriptRecord,
) -> tuple[list[dict[int, str], list[dict[int, str]]]]:
    """
    Search gRNA candidate for targeting the splice site. The candidate must have PAM site and the consensus sequence of the splice site.
    Exclude candidates that have "AAAA" or the exon length is a multiple of 3 or the exon has 3'UTR.
    ( If the target site is acceptor site, It is permissible to set the target exon as having a 3'UTR. )
    """
    orf = transcript_record.transcript_seq
    exon_start_pos = transcript_record.exon_start_positions
    exon_end_pos = transcript_record.exon_end_positions

    acceptor_cands = find_splice_site_candidate(orf, exon_start_pos[1:], 22, "AG", 2)
    donor_cands = find_splice_site_candidate(orf, exon_end_pos[:-1], 21, "GT", 1, False)
    index_exon_has_3_utr = next(
        (
            i
            for i, (start, end) in enumerate(zip(exon_start_pos, exon_end_pos))
            if start < transcript_record.cds_end <= end
        )
    )

    acceptor_candidates = filter_candidate(acceptor_cands, index_exon_has_3_utr, exon_start_pos, exon_end_pos)
    donor_candidates = filter_candidate(donor_cands, index_exon_has_3_utr, exon_start_pos, exon_end_pos, False)

    acceptor_candidates = link_to_crisperdirect(acceptor_candidates)
    donor_candidates = link_to_crisperdirect(donor_candidates)

    return acceptor_candidates, donor_candidates
