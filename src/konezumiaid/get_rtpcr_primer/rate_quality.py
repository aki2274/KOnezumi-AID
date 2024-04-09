from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def verify_crossing_exonjunction(
    candidate_primer_info: list[dict],
    exon_range: list[int],
) -> list[dict]:
    # primer_index_list means the indexes in exon. if return True, then the primer is in exon junction.
    return_info = candidate_primer_info.copy()
    for primer_data in return_info:
        left_primer_length = len(primer_data["left"])
        right_primer_length = len(primer_data["right"])
        # check if the primer is in any exon junction
        for exon_num in range(
            len(exon_range) - 1
        ):  # -1 because the last exon has no junction in the right side.
            if (exon_range[exon_num][1]) in range(
                primer_data["left_end"] - +left_primer_length,
                primer_data["left_end"],
            ):
                primer_data["left_cross_junction"] = True

            if (exon_range[exon_num][1] - 1) in range(
                primer_data["right_start"],
                primer_data["right_start"] + right_primer_length,
            ):
                primer_data["right_cross_junction"] = True
    return return_info


def autocorrect_intron_len(
    candidate_primer_info: list[dict],
    ds: GeneData,
) -> list[dict]:
    # the length is not considered the length of the exons that primer pairs are in.
    for primer_data in candidate_primer_info:
        intron_len = (
            ds.exon_start_list[primer_data["right_exon_num"]]
            - ds.exon_end_list[primer_data["left_exon_num"]]
            + 1
        )
        if intron_len > 0:
            primer_data["intron_len"] = intron_len
    return candidate_primer_info
