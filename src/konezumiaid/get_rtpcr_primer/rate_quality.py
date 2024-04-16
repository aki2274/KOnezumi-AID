from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def verify_cross_junction(
    primer3_result: list[dict],
    exon_range: list[int],
) -> list[dict]:
    # primer_index_list means the indexes in exon. if return True, then the primer is in exon junction.
    cross_valided = primer3_result.copy()
    for primer_pair in cross_valided:
        left_length = len(primer_pair["left_seq"])
        right_length = len(primer_pair["right_seq"])
        # check if the primer is in any exon junction
        for exon_num in range(
            len(exon_range) - 1
        ):  # -1 because the last exon has no junction in the right side.
            if (exon_range[exon_num][1]) in range(
                primer_pair["left_end"] - +left_length,
                primer_pair["left_end"],
            ):
                primer_pair["left_cross_junction"] = True

            if (exon_range[exon_num][1] - 1) in range(
                primer_pair["right_start"],
                primer_pair["right_start"] + right_length,
            ):
                primer_pair["right_cross_junction"] = True
    return cross_valided


def add_intron_len(
    candidate_primers: list[dict],
    ds: GeneData,
) -> list[dict]:
    # the length is not considered the length of the exons that primer pairs are in.
    intron_len_added = candidate_primers.copy()
    for primer_pair in intron_len_added :
        intron_len = (
            ds.exon_start_list[primer_pair["right_exon_num"]]
            - ds.exon_end_list[primer_pair["left_exon_num"]]
            + 1
        )
        if intron_len > 0:# if the intron length is negative, then the primer pair is in the same exon.
            primer_pair["intron_len"] = intron_len
    return intron_len_added
