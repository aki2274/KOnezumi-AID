from __future__ import annotations
import re


def search_candidate_codon_position(cds_seq: str) -> tuple[list[int], list[str]]:
    # search the position of candidate PTC (the index of C in CAA, CAG, CGA, or T in TGG)
    candidate_codons = re.finditer(r"(?=(CAA)|(?=(CAG))|(?=(CGA))|(?=(TGG)))", cds_seq)
    candidate_codon_positions = []
    candidate_codon = []
    for codon in candidate_codons:
        if codon.start() % 3 == 0:
            candidate_codon_positions.append(codon.start())
            candidate_codon.append(cds_seq[codon.start() : codon.start() + 3])
    return candidate_codon_positions, candidate_codon
