from __future__ import annotations
import primer3


def export_candidate_primer(exon_seq: str) -> dict:
    # the result is a dict of primer3 output.
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


def get_candidate_primer_pairs(
    primer_result: dict,
) -> list[dict]:
    # get dict primer quality score and the primer pairs
    # the primer quality score is will change in src/select_good_primers.py
    return [
        {
            "primer_score": 999,
            "intron_len": 0,
            "left_primer": left_primer_data["SEQUENCE"],
            "right_primer": right_primer_data["SEQUENCE"],
        }
        for (left_primer_data, right_primer_data) in zip(
            primer_result["PRIMER_LEFT"], primer_result["PRIMER_RIGHT"]
        )
    ]
