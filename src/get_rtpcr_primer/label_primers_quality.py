from __future__ import annotations


# before the function, we need to declare primer_quality_score. that list's length is the same number of primers.
def label_primers_quality(
    candidate_primer_info: list[dict],
    exon_junction_bool: list[bool],
) -> list[int]:
    # select primers that are in exon junction. label "1"
    # exon_junction_bool is the list of bool of whether the primer is in exon junction
    primer_num = 0
    for primer_data in candidate_primer_info:
        if exon_junction_bool[primer_num]:
            primer_data["primer_score"] = 1
        primer_num += 1
    # select primers that have intron length >1000, label "2", and intron length >500, label "3".
    for primer_data in candidate_primer_info:
        if primer_data["primer_score"] == 1:
            continue

        elif primer_data["intron_len"] >= 1000:
            primer_data["primer_score"] = 2

        elif primer_data["intron_len"] >= 500:
            primer_data["primer_score"] = 3

        else:
            primer_data["primer_score"] = 4

    return candidate_primer_info
