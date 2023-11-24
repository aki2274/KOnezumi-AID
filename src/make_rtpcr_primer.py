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


# select primers from the result of export_candidate_primer


def get_left_candidate_primer(
    primer_result: dict,
) -> list[str]:
    # get the left primer sequence from primer_result
    return [primer_data["SEQUENCE"] for primer_data in primer_result["PRIMER_LEFT"]]


def get_right_candidate_primer(
    primer_result: dict,
) -> list[str]:
    # get the right primer sequence from primer_result
    return [primer_data["SEQUENCE"] for primer_data in primer_result["PRIMER_RIGHT"]]


def get_candidate_primer_pairs(
    left_candidate_primer: list[str],
    right_candidate_primer: list[str],
) -> list[str]:
    # get the primer pairs from primer_result
    return list(zip(left_candidate_primer, right_candidate_primer))
