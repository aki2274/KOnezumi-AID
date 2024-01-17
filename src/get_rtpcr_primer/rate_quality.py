from __future__ import annotations
from dataclasses import dataclass


@dataclass
class GeneData:
    orf_seq: str
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


def verify_crossing_exonjunction(
    candidate_primer_info: list[dict],
    exon_range: list[int],
) -> list[dict]:
    # primer_index_list is the index in exon. if return True, then the primer is in exon junction.
    return_info = candidate_primer_info.copy()
    for primer_data in return_info:
        left_primer_length = len(primer_data["left_primer"])
        right_primer_length = len(primer_data["right_primer"])
        # check if the primer is in any exon junction
        for exon_num in range(
            len(exon_range) - 1
        ):  # -1 because the last exon has no junction in the right side.
            if (exon_range[exon_num][1]) in range(
                primer_data["left_primer_end"] - +left_primer_length,
                primer_data["left_primer_end"],
            ):
                primer_data["left_cross_junction"] = True

            if (exon_range[exon_num][1] - 1) in range(
                primer_data["right_primer_start"],
                primer_data["right_primer_start"] + right_primer_length,
            ):
                primer_data["right_cross_junction"] = True
    return return_info


def autocorrect_intron_len(
    candidate_primer_info: list[dict],
    ds: GeneData,
    # dsは、単一遺伝子のGeneDataクラスのインスタンス。転写産物名だけだと複数マッチしてる可能性ありなので修正いる？
) -> list[dict]:
    # get intron length between primer pairs
    intron_len_list = []
    # the length is not considered the length of the exons in primer pairs.
    for primer_data in candidate_primer_info:
        intron_len = (
            ds.exon_start_list[primer_data["right_primer_exon_num"]]
            - ds.exon_end_list[primer_data["left_primer_exon_num"]]
            + 1
        )
        if intron_len > 0:
            primer_data["intron_len"] = intron_len
    return candidate_primer_info
