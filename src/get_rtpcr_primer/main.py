from __future__ import annotations
import pickle
from pathlib import Path
import subprocess
from src.create_gene_dataclass import GeneData
from src.get_range_of_exon import get_exon_range
from src.nominate_candidate_stopcodon.generate_cds_seq import generate_exon_seq
from src.get_rtpcr_primer.make_rtpcr_primer import generate_candidate_info
from src.get_rtpcr_primer.rate_quality import (
    verify_crossing_exonjunction,
    autocorrect_intron_len,
)
from src.get_rtpcr_primer.export_fasta import export_fasta
from src.get_rtpcr_primer.add_uniqueness import add_uniqueness


def export_primers(ds: GeneData) -> list[dict]:
    """
    Export candidate rt-qPCR primers, based on exon seq.
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
    Path("data", "uniq").mkdir(parents=True, exist_ok=True)

    miss_0_path = Path("data", "uniq", "0_miss_counts.txt")
    miss_1_path = Path("data", "uniq", "1_miss_counts.txt")
    miss_2_path = Path("data", "uniq", "2_miss_counts.txt")

    # 1. get exon range
    exon_range = get_exon_range(ds)
    # 2. get exon seq
    exon_seq = generate_exon_seq(ds)
    # 4. get candidate primer info
    candidate_pairs = generate_candidate_info(exon_seq, exon_range)
    # 5. rate quality of candidate primer
    candidate_pairs = verify_crossing_exonjunction(candidate_pairs, exon_range)

    candidate_pairs = autocorrect_intron_len(candidate_pairs, ds)
    export_fasta(candidate_pairs)
    path_get_uniqueness = Path("src", "get_rtpcr_primer", "get_uniqueness.sh")
    process = subprocess.Popen(
        ["bash", path_get_uniqueness], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    _, stderr = process.communicate()
    candidate_pairs = add_uniqueness(
        candidate_pairs, miss_0_path, miss_1_path, miss_2_path
    )
    return candidate_pairs
