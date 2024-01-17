from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import subprocess

from src.get_rtpcr_primer.make_rtpcr_primer import (
    export_candidate,
    generate_candidate_info,
)
from get_rtpcr_primer.rate_quality import (
    verify_crossing_exonjunction,
    autocorrect_intron_len,
)
from src.get_rtpcr_primer.export_fasta import export_fasta
from src.get_rtpcr_primer.add_uniqueness import add_uniqueness


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


Path("data", "uniq").mkdir(parents=True, exist_ok=True)

miss_0_path = Path("data", "uniq", "0_miss_counts.txt")
miss_1_path = Path("data", "uniq", "1_miss_counts.txt")
miss_2_path = Path("data", "uniq", "2_miss_counts.txt")


def main(exon_seq: str, exon_range: tuple[int, int], ds: GeneData) -> list[dict]:
    """
    Export candidate rt-qPCR primers, based on exon seq.

    Args:
        exon_seq (str): a seq include cds and UTR.
        exon_range(list[int,int]): the range of the exon.
        ds(dataclass): a dataclass of one transcription.

    Returns:
        list[dict]: the candidate candidate rt-qPCR primer pairs.

    Example:
        >>> exon_seq = "ATGCAT...ATGCAT"
        >>> exon_range = [[1, 10], [10, 30], [30, 40], [40, 80]]
        >>> ds = GeneData(...)
        >>> main(exon_seq, exon_range, ds)
        {
            "left_cross_junction": 0,
            "right_cross_junction": 0,
            "intron_len": 0,
            "left_primer": "agcaaaagtgtgaagcgccc",
            "right_primer": "atctcgatcaccacgggctg",
            "left_primer_end": 29,
            "right_primer_start": 60,
            "left_primer_exon_num": 0,
            "right_primer_exon_num": 2,
        },....
    """
    primer3_result = export_candidate(exon_seq)
    candidate_pairs = generate_candidate_info(exon_seq, primer3_result, exon_range)
    candidate_pairs = verify_crossing_exonjunction(candidate_pairs, exon_range)
    candidate_pairs = autocorrect_intron_len(candidate_pairs, ds)
    export_fasta(candidate_pairs)
    path_get_uniqueness = Path("src", "get_rtpcr_primer", "get_uniqueness.sh")
    subprocess.run([path_get_uniqueness])
    candidate_pairs = add_uniqueness(
        candidate_pairs, miss_0_path, miss_1_path, miss_2_path
    )
    return candidate_pairs
