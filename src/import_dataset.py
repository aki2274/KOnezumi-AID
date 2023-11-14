from __future__ import annotations
import csv


# import ORF data each gene.
def read_fasta(path: str) -> dict[str, str]:
    sequences = {}
    seq_id = None
    current_seq = ""
    with open(path, "r") as file:
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


# import anotation data each transcription
def read_csv(path: str) -> list[dict[str, str]]:
    result = []
    with open(path, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=",")
        for row in reader:
            result.append(row)
    return result
