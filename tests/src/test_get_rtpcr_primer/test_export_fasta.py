from __future__ import annotations
import pytest
from src.get_rtpcr_primer.export_fasta import export_fasta


def test_exprt_fasta():
    data = [
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
        }
    ]
    export_fasta(data)
    with open("data/reads/candidateprimer.fa", "r") as file:
        assert (
            file.read()
            == ">agcaaaagtgtgaagcgccc\nagcaaaagtgtgaagcgccc\n>atctcgatcaccacgggctg\natctcgatcaccacgggctg\n"
        )
