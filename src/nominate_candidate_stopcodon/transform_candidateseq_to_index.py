from __future__ import annotations
import re


def transform_ct_guideseq_to_index(orf_seq: str, targets: list[str]) -> list[int]:
    # Get the index of "C" in "CAA", "CAG", or "CGA" in candidate gRNA (the index on the ORF)
    positions = []
    for target in targets:
        # Add the index of the candidate gRNA start position in the ORF and the index of "C" form gRNA start position
        for match in re.finditer(target, orf_seq):
            position = (
                re.search(r"(CAA|CAG|CGA)", match.group()).start() + match.start()
            )
            positions.append(position)
    positions = list(dict.fromkeys(positions))
    return positions


def transform_ga_guideseq_to_index(orf_seq: str, targets: list[str]) -> list[int]:
    # Get the index of "T" in "TGG" in candidate gRNA (the index on the ORF)
    # To reduce the work of comparison with candidate PTC (TGG), obtain the index of "T".
    positions = []
    for target in targets:
        for match in re.finditer(target, orf_seq):
            rev_match = match.group()[::-1]
            add_num = re.search("GGT", rev_match).start()
            # Get the index of the candidate gRNA end position in the ORF and back to "T" in "TGG"
            positions.append(match.start() + len(match.group()) - 3 - add_num)
    positions = list(dict.fromkeys(positions))
    return positions
