from __future__ import annotations


# before the function, we need to declare primer_quality_score. that list's length is the same number of primers.
def select_exon_juncton_primers(
    primer_quality_score: list[int],
    left_junction_bool: list[bool],
    right_junction_bool: list[bool],
) -> list[int]:
    # select primers that are in exon junction. label "1"
    # left_junction_bool and right_junction_bool are the list of bool of whether the primer is in exon junction
    for primer_num in range(len(primer_quality_score)):
        if any(left_junction_bool[primer_num], right_junction_bool[primer_num]):
            primer_quality_score.append(1)
    return primer_quality_score


def label_primers_quality(
    primer_quality_score: list[int],
    primer_pair_intron_len: list[int],
) -> list[str]:
    # select primers that have intron length >1000, label "2", and intron length >500, label "3".
    for primer_num in range(len(primer_quality_score)):
        if primer_quality_score[primer_num] == 1:
            continue
        elif primer_pair_intron_len[primer_num] > 1000:
            primer_quality_score.append(2)
        elif primer_pair_intron_len[primer_num] > 500:
            primer_quality_score.append(3)
        else:
            primer_quality_score.append(4)
    return primer_quality_score
