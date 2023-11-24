from __future__ import annotations


def get_revcomp(orf_seq: str) -> str:
    """Reverse complement a DNA sequence."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(complement[nt.upper()] for nt in orf_seq[::-1])
