from __future__ import annotations
from konezumiaid.create_gene_dataclass import GeneData


def is_crossing_juncion(
    candidate_primers: list[dict],
    exon_range: list[tuple[int, int]],
) -> list[dict]:
    validated_primers = candidate_primers.copy()

    for primer_pair in validated_primers:
        left_length = len(primer_pair["left_seq"])
        right_length = len(primer_pair["right_seq"])
        left_end = primer_pair["left_end"]
        right_start = primer_pair["right_start"]

        is_crossing_left = any((left_end - left_length < exon_end <= left_end) for _, exon_end in exon_range)

        is_crossing_right = any((right_start < exon_end <= right_start + right_length) for _, exon_end in exon_range)

        primer_pair["left_cross_junction"] = is_crossing_left
        primer_pair["right_cross_junction"] = is_crossing_right

    return validated_primers


def add_intron_len(
    candidate_primers: list[dict],
    transcript_record: GeneData,
) -> list[dict]:
    # The additional length does not account for the size of the exons housing the primer pairs."
    added_len = candidate_primers.copy()
    for primer_pair in added_len:
        intron_len = (
            transcript_record.exon_start_list[primer_pair["right_exon_num"] - 1]
            - transcript_record.exon_end_list[primer_pair["left_exon_num"] - 1]
            + 1
        )
        if intron_len > 0:
            # If the intron length is negative, then the primer pair is in the same exon.
            primer_pair["intron_len"] = intron_len
    return added_len
