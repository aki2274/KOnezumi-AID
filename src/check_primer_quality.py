import __future__
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
    primer_list: list[str], primer_index_list: list[int], exon_range: list[int]
) -> bool:
    # primer_index_list is the index in exon. if return True, then the primer is in exon junction.
    result = []
    for primer_num in range(len(primer_index_list)):
        primer_longth = len(primer_list[primer_num])
        bool_list = []
        # check if the primer is in any exon junction
        for exon_num in range(len(exon_range)):
            if exon_range[exon_num][1] in range(
                primer_index_list[primer_num],
                primer_index_list[primer_num] + primer_longth + 1,
            ):
                bool_list.append(True)
            else:
                bool_list.append(False)
        result.append(any(bool_list))

    return result


def get_primer_pair_intron_len(
    left_primer_num_list: list[int], right_primer_num_list: list[int], ds: GeneData
) -> list[int]:
    # get intron length between primer pairs
    intron_len_list = []
    # the length is not considered the length of the exons in primer pairs.
    for primer_num in range(len(left_primer_num_list)):
        intron_len = (
            ds.exon_start_list[right_primer_num_list[primer_num]]
            - ds.exon_end_list[left_primer_num_list[primer_num]]
        )
        intron_len_list.append(intron_len)
    # if the primer pair is in the same exon, the intron length is 0.
    # (the showed number of intron length is negative if primer_num is the same.)
    intron_natural_len = [max(0, i) for i in intron_len_list]
    return intron_natural_len
