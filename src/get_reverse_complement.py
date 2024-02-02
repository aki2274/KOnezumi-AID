from __future__ import annotations


def get_revcomp(orf_seq: str) -> str:
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
    return "".join(complement[nt] for nt in orf_seq[::-1])
