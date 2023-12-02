from __future__ import annotations
import pickle


def read_fasta(path_fasta: str) -> dict[str, str]:
    fasta = {}
    identifier = None
    sequence = ""

    with open(path_fasta, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    fasta[identifier] = sequence
                identifier = line.lstrip(">")
                sequence = ""
            else:
                sequence += line
        fasta[identifier] = sequence

    return fasta


def read_pkl(path_pkl: str) -> list[dict[str, str]]:
    with open(path_pkl, mode="rb") as pklfile:
        return pickle.load(pklfile)
