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


def check_exon_junction(
    candidate_primer_info: list[dict],
    exon_range: list[int],
) -> bool:
    # primer_index_list is the index in exon. if return True, then the primer is in exon junction.
    result = []
    for primer_data in candidate_primer_info:
        left_primer_longth = len(primer_data["left_primer"])
        right_primer_longth = len(primer_data["right_primer"])
        bool_list = []
        # check if the primer is in any exon junction
        for exon_num in range(
            len(exon_range) - 1
        ):  # -1 because the last exon has no junction in the right side.
            if (exon_range[exon_num][1] - 1) in range(
                primer_data["left_primer_start"],
                primer_data["left_primer_start"] + left_primer_longth + 1,
            ) or (exon_range[exon_num][1] - 1) in range(
                primer_data["right_primer_start"],
                primer_data["right_primer_start"] + right_primer_longth + 1,
            ):
                bool_list.append(True)
            else:
                bool_list.append(False)
        result.append(any(bool_list))

    return result


def rewrite_primer_pair_intron_len(
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
        )
        intron_len_list.append(intron_len)
    # if the primer pair is in the same exon, the intron length is 0.
    # (the showed number of intron length is negative if primer_num is the same.)
    for intron_len, primer_data in zip(intron_len_list, candidate_primer_info):
        if intron_len > 0:
            primer_data["intron_len"] = intron_len
    return candidate_primer_info
