from __future__ import annotations
import primer3
from konezumiaid.get_reverse_complement import get_revcomp


def export_candidate(exon_seq: str) -> dict:
    # the result is a dict of primer3 output.
    # those parameters are referred from TAKARA bio.
    result = primer3.bindings.design_primers(
        seq_args={
            "SEQUENCE_TEMPLATE": exon_seq,
        },
        global_args={
            "PRIMER_NUM_RETURN": 30,
            "PRIMER_PRODUCT_SIZE_RANGE": [[80, 150], [150, 300]],
            "PRIMER_MIN_SIZE": 17,
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MAX_SIZE": 25,
            "PRIMER_PAIR_MAX_DIFF_TM": 4,
            "PRIMER_OPT_TM": 62.0,
            "PRIMER_MIN_TM": 60.0,
            "PRIMER_MAX_TM": 65.0,
            "PRIMER_MIN_GC": 40.0,
            "PRIMER_OPT_GC_PERCENT": 50.0,
            "PRIMER_MAX_GC": 60.0,
            "PRIMER_MAX_SELF_ANY": 8,
            "PRIMER_MAX_SELF_END": 3,
            "PRIMER_MAX_POLY_X": 4,
            "PRIMER_WT_TM_GT": 0.5,
            "PRIMER_WT_TM_LT": 0.5,
            "PRIMER_WT_SELF_ANY": 0.5,
            "PRIMER_WT_SELF_END": 1,
            "PRIMER_PAIR_WT_DIFF_TM": 0.5,
            "PRIMER_PAIR_WT_COMPL_ANY": 0.5,
            "PRIMER_PAIR_WT_COMPL_END": 1,
        },
    )
    return result


def generate_candidate_info(
    exon_seq: str,
    exon_range: list[list[int, int]],
) -> list[dict]:
    """
    Get dict of primer quality info and the primer pairs
    _cross_junction : the primer cross the exon junction or not
    intron_len : the length of intron between the two primers
    left_primer_end : the end position of left primer
    right_primer_start : the start position of right primer
    _primer_exon_num : the exon number of left primer
    """
    primer_result = export_candidate(exon_seq)
    candidate_primer_info = [
        {
            "left": left_primer_data["SEQUENCE"],
            "right": right_primer_data["SEQUENCE"],
            "left_tm": left_primer_data["TM"],
            "right_tm": right_primer_data["TM"],
            "left_end": exon_seq.find(left_primer_data["SEQUENCE"])
            + len(left_primer_data["SEQUENCE"]),
            "right_start": exon_seq.find(get_revcomp(right_primer_data["SEQUENCE"])),
            "product_size": pair_data["PRODUCT_SIZE"],
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
        }
        for (pair_data, left_primer_data, right_primer_data) in zip(
            primer_result["PRIMER_PAIR"],
            primer_result["PRIMER_LEFT"],
            primer_result["PRIMER_RIGHT"],
        )
    ]
    for primer in candidate_primer_info:
        exon_start = 0
        exon_end = 1
        for s in range(len(exon_range)):
            if (
                exon_range[s][exon_start]
                <= primer["left_end"]
                <= exon_range[s][exon_end]
            ):
                primer["left_exon_num"] = s
            if (
                exon_range[s][exon_start]
                <= primer["right_start"]
                <= exon_range[s][exon_end]
            ):
                primer["right_exon_num"] = s
    return candidate_primer_info
