from __future__ import annotations
from pathlib import Path


def read_fasta(file_path: str | Path) -> dict:
    sequences = {}
    seq_id = None
    current_seq = ""

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                if seq_id is not None:
                    sequences[seq_id] = current_seq
                    current_seq = ""
                seq_id = line[1:]
            else:
                current_seq += line

        if seq_id is not None and current_seq != "":
            sequences[seq_id] = current_seq

    return sequences
