from __future__ import annotations


def get_revcomp(seq: str) -> str:
    """Reverse complement a sequence."""
    complement = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "N": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "n": "n",
    }
    return "".join(complement[nt] for nt in seq[::-1])
