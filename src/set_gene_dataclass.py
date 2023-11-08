from __future__ import annotations
import pandas as pd
from dataclasses import dataclass


@dataclass
class DataClass:
    orf_seq: str
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


def set_dataclass(
    transcripts_name: str, gene_df: pd.DataFrame, gene_seq_data: dict
) -> DataClass:
    # make dataset from refFlat.txt(gene_df) data
    data = gene_df[gene_df["name"] == transcripts_name].reset_index()
    chrom = str(data.loc[0, "chrom"])
    txStart = str(data.loc[0, "txStart"])
    txEnd = str(data.loc[0, "txEnd"])
    query = f"{transcripts_name}::{chrom}:{txStart}-{txEnd}"
    txStart = int(data["txStart"].iloc[0])
    txEnd = int(data["txEnd"].iloc[0])
    orf_seq = gene_seq_data[query]
    cdsStart = int(data["cdsStart"].iloc[0])
    cdsEnd = int(data["cdsEnd"].iloc[0])
    exonCount = int(data["exonCount"].iloc[0])
    start = data.at[0, "exonStarts"]
    end = data.at[0, "exonEnds"]
    exon_start_list = [int(x) for x in list(start.split(","))]
    exon_end_list = [int(x) for x in list(end.split(","))]
    set_data = DataClass(
        orf_seq,
        txStart,
        txEnd,
        cdsStart,
        cdsEnd,
        exonCount,
        exon_start_list,
        exon_end_list,
    )
    return set_data
