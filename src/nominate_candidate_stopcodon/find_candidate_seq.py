from __future__ import annotations
import re


def find_ct_target_seq(orf_seq: str) -> list[str]:
    matches = re.findall(r"(?=(\w{21}GG))", orf_seq)
    # Check for the presence of "CAA", "CAG", or "CGA" in the next 3 positions.
    targets = [
        match
        for match in matches
        if any(match[i : i + 3] in {"CAA", "CAG", "CGA"} for i in range(1, 4))
    ]
    # Remove duplicates while preserving order
    return list(dict.fromkeys(targets))


def find_ga_target_seq(orf_seq: str) -> list[str]:
    matches = re.findall(r"(?=(CC\w{18,21}))", orf_seq)
    # Check for the presence of "TGG" in the next 4 positions.
    targets = [
        match
        for match in matches
        if any(match[i : i + 3] == "TGG" for i in range(17, 21))
    ]
    # Remove duplicates while preserving order
    return list(dict.fromkeys(targets))
