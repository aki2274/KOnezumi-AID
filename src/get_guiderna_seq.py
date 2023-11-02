from __future__ import annotations
import re


def get_revcomp(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement[nt] for nt in seq[::-1])


def find_ct_target_seq(seq: str) -> list[str]:
    matches = re.findall(r"(?=(\w{21}GG))", seq)
    # Check for the presence of "CAA", "CAG", or "CGA" in the next 3 positions.
    targets = [
        match
        for match in matches
        if any(match[i : i + 3] in {"CAA", "CAG", "CGA"} for i in range(1, 4))
    ]
    # Remove duplicates while preserving order
    return list(dict.fromkeys(targets))


def find_ga_target_seq(seq: str) -> list[str]:
    matches = re.findall(r"(?=(CC\w{18,21}))", seq)
    # Check for the presence of "TGG" in the next 4 positions.
    targets = [
        match
        for match in matches
        if any(match[i : i + 3] == "TGG" for i in range(17, 21))
    ]
    # Remove duplicates while preserving order
    targets = list(dict.fromkeys(targets))
    return targets
