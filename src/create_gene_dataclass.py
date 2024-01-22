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
    data_filtered_transcript = next(
        (
            gene_data
            for gene_data in refflat_data
            if gene_data["name"] == transcript_name
        ),
        None,
    )
    if data_filtered_transcript is None:
        raise ValueError("Transcript name not found in refflat data")

    chrom = str(data_filtered_transcript["chrom"])
    txStart = int(data_filtered_transcript["txStart"])
    txEnd = int(data_filtered_transcript["txEnd"])
    cdsStart = int(data_filtered_transcript["cdsStart"])
    cdsEnd = int(data_filtered_transcript["cdsEnd"])
    exonCount = int(data_filtered_transcript["exonCount"])
    exon_start_list = [
        int(x) for x in data_filtered_transcript["exonStarts"].split(",")
    ]
    exon_end_list = [int(x) for x in data_filtered_transcript["exonEnds"].split(",")]

    query = f"{transcript_name}::{chrom}:{txStart}-{txEnd}"
    orf_seq = gene_seq_data.get(query)
    if orf_seq is None:
        raise ValueError("Query not found in gene sequence data")

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
