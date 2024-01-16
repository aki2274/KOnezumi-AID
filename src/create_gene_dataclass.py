from __future__ import annotations
from dataclasses import dataclass


@dataclass
class GeneData:
    orf_seq: str
    txStart: int
    txend: int
    cdsStart: int
    cdsEnd: int
    exonCount: int
    exon_start_list: list[int]
    exon_end_list: list[int]


def create_dataclass(
    transcript_name: str, refflat_data: list[dict], gene_seq_data: dict
) -> GeneData:
    # make dataset from refFlat.txt(exported csv) data
    data_filtered_transcript = [
        gene_data for gene_data in refflat_data if gene_data["name"] == transcript_name
    ][
        0
    ]  # 転写産物名に重複がないことを確認する
    chrom = str(data_filtered_transcript["chrom"])
    txStart = str(data_filtered_transcript["txStart"])
    txEnd = str(data_filtered_transcript["txEnd"])
    query = f"{transcript_name}::{chrom}:{txStart}-{txEnd}"  # 重複がないことを確認する必要がある。
    orf_seq = gene_seq_data[query]
    txStart = int(txStart)
    txEnd = int(txEnd)
    cdsStart = int(data_filtered_transcript["cdsStart"])
    cdsEnd = int(data_filtered_transcript["cdsEnd"])
    exonCount = int(data_filtered_transcript["exonCount"])
    start = data_filtered_transcript["exonStarts"]
    end = data_filtered_transcript["exonEnds"]
    exon_start_list = [int(x) for x in list(start.split(","))]
    exon_end_list = [int(x) for x in list(end.split(","))]
    set_data = GeneData(
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
