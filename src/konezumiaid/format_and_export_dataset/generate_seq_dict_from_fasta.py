from __future__ import annotations
import pandas as pd
from konezumiaid.get_reverse_complement import get_revcomp


def read_fasta(path_fasta: str) -> dict[str, str]:
    fasta = {}
    transcript_key = None
    sequence = ""
    with open(path_fasta, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if transcript_key is not None:
                    fasta[transcript_key] = sequence
                transcript_key = line.lstrip(">")
                sequence = ""
            else:
                sequence += line
        fasta[transcript_key] = sequence
    return fasta


def create_strand_plus_seq_dict(
    df_refflat: pd.DataFrame,
    transcript_seq_dict: dict[str, str],
) -> dict[str, str]:
    """convert the sequence dictionary from fasta to the sequence dictionary.
    Reverse the sequence and change key if the transcript is negative strand."""
    fasta_keys = df_refflat["name"].tolist()
    result_dict = {}

    for fasta_key in fasta_keys:
        strand = df_refflat.loc[df_refflat["name"] == fasta_key, "strand"].iloc[0]
        if strand == "-":
            converted_seq = get_revcomp(transcript_seq_dict[fasta_key])
        else:
            converted_seq = transcript_seq_dict[fasta_key]

        result_dict[fasta_key] = converted_seq.upper()
    return result_dict
