from __future__ import annotations
import re


def get_candidate_stopcodon_index(cds_seq: str) -> list:
    # search the position of candidate stopcodon (the index of C in CAA, CAG, CGA, or T in TGG)
    matches = re.finditer(r"(?=(CAA)|(?=(CAG))|(?=(CGA))|(?=(TGG)))", cds_seq)
    candidate_codon_index_list = [
        match.start() for match in matches if (match.start() % 3) == 0
    ]
    return candidate_codon_index_list
