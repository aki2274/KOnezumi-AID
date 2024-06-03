from __future__ import annotations
import primer3
from konezumiaid.get_reverse_complement import get_revcomp


def export_candidate(exon_seq: str) -> dict:
    """
    Export candidate primers for a given exon sequence.
    those parameters are referred from TAKARA bio.
    https://www.takara-bio.co.jp/research/prt/pdfs/prt3-1.pdf

    Args:
        exon_seq (str): The sequence of the exon.

    Returns:
        dict: A dictionary containing the primer3 output.

    """
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
    spliced_exon_seq: str,
    exon_range: list[list[int, int]],
) -> list[dict]:
    """
    Get dict of primer quality info and the primer pairs

    Args:
        exon_seq (str): The sequence of the exon.
        exon_range (list[list[int, int]]): The range of each exon.

    Returns:
        list[dict]: A list of dictionaries containing primers information.

    Each dictionary in the list contains the following keys:
        - left_seq (str): The sequence of the left primer.
        - right_seq (str): The sequence of the right primer.
        - left_tm (float): The melting temperature of the left primer.
        - right_tm (float): The melting temperature of the right primer.
        - left_end (int): The end position of the left primer in the exon sequence.
        - right_start (int): The start position of the right primer in the exon sequence.
        - product_size (int): The size of the PCR product amplified by the primer pair.
        - left_cross_junction (bool): Indicates whether the left primer crosses the exon junction.
        - right_cross_junction (bool): Indicates whether the right primer crosses the exon junction.
        - intron_len (int): The length of the intron between the two primers.
        - left_exon_num (int): The exon number of the left primer.
        - right_exon_num (int): The exon number of the right primer.
    """
    primer3_return = export_candidate(spliced_exon_seq)
    candidate_primers = [
        {
            "left": primer_left["SEQUENCE"],
            "right": primer_right["SEQUENCE"],
            "left_tm": primer_left["TM"],
            "right_tm": primer_right["TM"],
            "left_end": spliced_exon_seq.find(primer_left["SEQUENCE"]) + len(primer_left["SEQUENCE"]),
            "right_start": spliced_exon_seq.find(get_revcomp(primer_right["SEQUENCE"])),
            "product_size": primer_pair["PRODUCT_SIZE"],
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
        }
        for (primer_pair, primer_left, primer_right) in zip(
            primer3_return["PRIMER_PAIR"],
            primer3_return["PRIMER_LEFT"],
            primer3_return["PRIMER_RIGHT"],
        )
    ]
    result = []
    for primer in candidate_primers:
        for s, (exon_start, exon_end) in enumerate(exon_range):
            left_end_within_exon = exon_start <= primer["left_end"] <= exon_end
            right_start_within_exon = exon_start <= primer["right_start"] <= exon_end
            if left_end_within_exon:
                primer["left_exon_num"] = s + 1
            if right_start_within_exon:
                primer["right_exon_num"] = s + 1
        if primer["left_exon_num"] == primer["right_exon_num"]:
            continue
        result.append(primer)
    return result
