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


def create_dict_keys(df_refflat: pd.DataFrame) -> dict[str, str]:
    """Create a unique key for each transcript for the sequence dictionary."""
    transcript_names = df_refflat["name"].tolist()
    keys = []
    for transcript_name in transcript_names:
        df_transcript = df_refflat[df_refflat["name"] == transcript_name]
        chrom = str(df_transcript["chrom"].iloc[0])
        txStart = str(df_transcript["txStart"].iloc[0])
        txEnd = str(df_transcript["txEnd"].iloc[0])
        query = f"{transcript_name}::{chrom}:{txStart}-{txEnd}"
        keys.append(query)
    return keys


def create_strand_plus_seq_dict(
    df_refflat: pd.DataFrame,
    df_refflat_sorted: pd.DataFrame,
    transcript_seq_dict: dict[str, str],
) -> dict[str, str]:
    """convert the sequence dictionary from fasta to the sequence dictionary.
    Reverse the sequence and change key if the transcript is negative strand."""
    transcripts_fasta_keys = create_dict_keys(df_refflat)
    sorted_keys = create_dict_keys(df_refflat_sorted)
    name_list = df_refflat["name"].tolist()
    result_dict = {}

    for transcript_name, fasta_key, sorted_key in zip(
        name_list, transcripts_fasta_keys, sorted_keys
    ):
        strand = df_refflat.loc[df_refflat["name"] == transcript_name, "strand"].iloc[0]
        if strand == "-":
            converted_seq = get_revcomp(transcript_seq_dict[fasta_key])
        else:
            converted_seq = transcript_seq_dict[fasta_key]

        result_dict[sorted_key] = converted_seq.upper()
    return result_dict
