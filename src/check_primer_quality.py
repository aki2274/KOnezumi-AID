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
    candidate_primer_pairs: list[dict],
    primer_index_list: list[int],  # マッパーで得られるプライマーの位置なので型が違うかも。
    exon_range: list[int],
) -> bool:
    # primer_index_list is the index in exon. if return True, then the primer is in exon junction.
    result = []
    for left_primer, right_primer in zip(
        candidate_primer_pairs["left_primer"], candidate_primer_pairs["right_primer"]
    ):
        left_primer_longth = len(left_primer)
        right_primer_longth = len(right_primer)
        bool_list = []
        # check if the left primer is in any exon junction
        for exon_num in range(len(exon_range)):
            if exon_range[exon_num][1] in range(
                primer_index_list[exon_num],
                primer_index_list[exon_num] + left_primer_longth + 1,
            ) or exon_range[exon_num][1] in range(
                primer_index_list[exon_num],
                primer_index_list[exon_num] + right_primer_longth + 1,
            ):
                bool_list.append(True)
            else:
                bool_list.append(False)
        result.append(any(bool_list))

    return result


def get_primer_pair_intron_len(
    candidate_primer_pairs: list[dict],
    left_primer_exon_num: list[int],
    right_primer_exon_num: list[int],
    ds: GeneData,
) -> list[int]:
    # get intron length between primer pairs
    intron_len_list = []
    # the length is not considered the length of the exons in primer pairs.
    for primer_num in range(len(left_primer_exon_num)):
        intron_len = (
            ds.exon_start_list[right_primer_exon_num[primer_num]]
            - ds.exon_end_list[left_primer_exon_num[primer_num]]
        )
        intron_len_list.append(intron_len)
    # if the primer pair is in the same exon, the intron length is 0.
    # (the showed number of intron length is negative if primer_num is the same.)
    for intron_len, primer_data in zip(intron_len_list, candidate_primer_pairs):
        if intron_len > 0:
            primer_data["intron_len"] = intron_len
    return candidate_primer_pairs
