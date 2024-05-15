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
    gene_data: pd.DataFrame, sort_gene_dataframe: pd.DataFrame, gene_seq: dict[str, str]
) -> dict[str, str]:
    """convert the sequence dictionary from fasta to the sequence dictionary.
    reverse the sequence if the transcript is negative strand."""
    normal_query = create_dict_keys(gene_data)
    sorted_query = create_dict_keys(sort_gene_dataframe)
    name_list = gene_data["name"].tolist()
    copy_dict = gene_seq.copy()
    for transcript_name, n_que in zip(name_list, normal_query):
        fil_data = gene_data[gene_data["name"] == transcript_name]
        if (fil_data["strand"] == "-").any():
            copy_dict[n_que] = get_revcomp(gene_seq[n_que])
    result = {}
    for n_que, s_que in zip(normal_query, sorted_query):
        result[s_que] = copy_dict[n_que].upper()
    return result
