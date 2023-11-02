from __future__ import annotations

import re
import numpy as np
import pandas as pd


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


def get_index_of_ct_target_seq(seq: str, targets: list[str]) -> list[int]:
    # Get the index of "C" in "CAA", "CAG", or "CGA" in candidate gRNA
    positions = []
    for target in targets:
        # Add the index of the candidate gRNA start position in the ORF and the index of "C" form gRNA start position
        for match in re.finditer(target, seq):
            position = (
                re.search(r"(CAA|CAG|CGA)", match.group()).start() + match.start()
            )
            positions.append(position)
    positions = list(dict.fromkeys(positions))
    return positions


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


def get_index_of_ga_target_seq(seq: str, targets: list[str]) -> list[int]:
    # Get the index of "T" in "TGG" in candidate gRNA
    positions = []
    for target in targets:
        for match in re.finditer(target, seq):
            rev_match = match.group()[::-1]
            add_num = re.search("GGT", rev_match).start()
            # Get the index of the candidate gRNA end position in the ORF and back to "T" in "TGG"
            positions.append(match.start() + len(rev_match) - 3 - add_num)
    positions = list(dict.fromkeys(positions))
    return positions


def get_targetstart_num(
    name: str, target: list[int], gene_df: pd.DataFrame
) -> list[int]:
    data = gene_df[gene_df["name"] == name]
    num = data["txStart"].astype(int).to_numpy()
    targetable_seq_idex = np.array(target) + num
    return targetable_seq_idex  # list かnp.arrayかは後の処理にあわせる
