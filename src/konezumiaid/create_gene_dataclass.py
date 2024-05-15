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
    transcript_name: str, refflat: list[dict], transcript_seq_dict: dict
) -> GeneData:
    # Create dataclass from the transcript name.
    transcript_filtered = next(
        (
            transcript_data
            for transcript_data in refflat
            if transcript_data["name"] == transcript_name
        ),
        None,
    )
    if transcript_filtered is None:
        raise ValueError("Transcript name doesn't exist in the refflat")

    chrom = str(transcript_filtered["chrom"])
    txStart = int(transcript_filtered["txStart"])
    txEnd = int(transcript_filtered["txEnd"])
    cdsStart = int(transcript_filtered["cdsStart"])
    cdsEnd = int(transcript_filtered["cdsEnd"])
    exonCount = int(transcript_filtered["exonCount"])
    exon_start_list = [int(x) for x in transcript_filtered["exonStarts"].split(",")]
    exon_end_list = [int(x) for x in transcript_filtered["exonEnds"].split(",")]

    key = f"{transcript_name}::{chrom}:{txStart}-{txEnd}"
    orf_seq = transcript_seq_dict.get(key)
    if orf_seq is None:
        raise ValueError("Transcript sequence doesn't exist in the transcript_seq_dict")

    transcript_record = GeneData(
        orf_seq,
        txStart,
        txEnd,
        cdsStart,
        cdsEnd,
        exonCount,
        exon_start_list,
        exon_end_list,
    )
    return transcript_record
