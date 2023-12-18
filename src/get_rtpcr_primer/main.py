from __future__ import annotations
from dataclasses import dataclass
import subprocess

from src.get_rtpcr_primer.make_rtpcr_primer import (
    export_candidate_primer,
    get_candidate_primer_pairs,
)
from src.get_rtpcr_primer.check_primer_quality import (
    check_exon_junction,
    rewrite_primer_pair_intron_len,
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
        >>> import cstag
        >>> exon_seq = "ATGCAT...ATGCAT"
        >>> exon_range = [[1, 10], [10, 30], [30, 40], [40, 80]]
        >>> ds = GeneData(...)
        >>> cstag.call(exon_seq, exon_range, ds)
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
    primer3_result = export_candidate_primer(exon_seq)
    candidate_pairs = get_candidate_primer_pairs(exon_seq, primer3_result, exon_range)
    candidate_pairs = check_exon_junction(candidate_pairs, exon_range)
    candidate_pairs = rewrite_primer_pair_intron_len(candidate_pairs, ds)
    export_fasta(candidate_pairs)
    subprocess.run([r"src\get_rtpcr_primer\get_uniqueness.sh"])
    candidate_pairs = add_uniqueness(candidate_pairs)
    return candidate_pairs
