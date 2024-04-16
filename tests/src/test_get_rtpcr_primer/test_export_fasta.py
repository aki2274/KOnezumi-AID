from __future__ import annotations
from konezumiaid.get_rtpcr_primer.export_fasta import export_fasta


def test_exprt_fasta():
    data = [
        {
            "left_cross_junction": 0,
            "right_cross_junction": 0,
            "intron_len": 0,
            "left": "agcaaaagtgtgaagcgccc",
            "right": "atctcgatcaccacgggctg",
            "left_end": 29,
            "right_start": 60,
            "left_exon_num": 0,
            "right_exon_num": 2,
        }
    ]
    export_fasta(data, "tests/data/reads/candidateprimer.fa")
    with open("tests/data/reads/candidateprimer.fa", "r") as file:
        assert (
            file.read()
            == ">agcaaaagtgtgaagcgccc\nagcaaaagtgtgaagcgccc\n>atctcgatcaccacgggctg\natctcgatcaccacgggctg\n"
        )
