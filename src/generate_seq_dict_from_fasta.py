from __future__ import annotations
from pandas import DataFrame
from pathlib import Path
from src.get_reverse_complement import get_revcomp
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


def create_dict_keys(genedata: DataFrame) -> dict[str, str]:
    # get query list, which is the key of the sequence dictionary
    name_list = genedata["name"].tolist()
    query_list = []
    for transcript_name in name_list:
        data_filtered_transcript = genedata[genedata["name"] == transcript_name]
        chrom = str(data_filtered_transcript["chrom"].iloc[0])
        txStart = str(data_filtered_transcript["txStart"].iloc[0])
        txEnd = str(data_filtered_transcript["txEnd"].iloc[0])
        query = f"{transcript_name}::{chrom}:{txStart}-{txEnd}"
        query_list.append(query)
    return query_list


def create_sorted_seq_dict(
    gene_data: DataFrame, sort_gene_dataframe: DataFrame, gene_seq: dict[str, str]
) -> dict[str, str]:
    # convert the sequence dictionary from fasta to the sorted sequence dictionary.
    # reverse the sequence if the transcript is on the negative strand
    normal_query = create_dict_keys(gene_data)
    sorted_query = create_dict_keys(sort_gene_dataframe)
    name_list = gene_data["name"].tolist()
    for transcript_name, n_que in zip(name_list, normal_query):
        fil_data = gene_data[gene_data["name"] == transcript_name]
        if (fil_data["strand"] == "-").any():
            gene_seq[n_que] = get_revcomp(gene_seq[n_que])
    result = {}
    for n_que, s_que in zip(normal_query, sorted_query):
        result[s_que] = gene_seq[n_que]
    return result
