from __future__ import annotations
from pathlib import Path
import subprocess
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.get_range_of_exon import get_exon_range
from konezumiaid.nominate_candidate_stopcodon.generate_cds_seq import generate_exon_seq
from konezumiaid.get_rtpcr_primer.make_rtpcr_primer import generate_candidate_info
from konezumiaid.get_rtpcr_primer.rate_quality import (
    verify_cross_junction,
    add_intron_len,
)
from konezumiaid.get_rtpcr_primer.export_fasta import export_fasta
from konezumiaid.get_rtpcr_primer.add_uniqueness import add_uniqueness


def export_primers(ds: GeneData) -> list[dict]:
    """
    Export candidate rt-qPCR primers, based on exon seq.
    Args:
        ds(GeneData): The Dataclass named GeneData.
    Returns:
        list[dict]: A list of dictionaries containing the candidate primer information.
    Example output:
        [{
            "left_seq": "agcaaaagtgtgaagcgccc",
            "right_seq": "atctcgatcaccacgggctg",
            "left_tm": 61.996786,
            "right_tm": 61.996786,
            "left_end": 2121,
            "right_start": 2162,
            "product_size": 81,
            "left_cross_junction": False,
            "right_cross_junction": False,
            "intron_len": 0,
            "left_exon_num": 2,
            "right_exon_num": 2,
            "left_0_mismatch": 1
            "left_1_mismatch": 10
            "left_2_mismatch": 25
            "right_0_mismatch": 1
            "right_1_mismatch": 16
            "right_2_mismatch": 56
        }...]
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
    primer3_result = generate_candidate_info(exon_seq, exon_range)
    # 5. rate quality of candidate primer
    cross_validated = verify_cross_junction(primer3_result, exon_range)

    intron_len_added = add_intron_len(cross_validated, ds)
    fasta_path = Path("data", "uniq", "candidateprimer.fa")
    export_fasta(intron_len_added, fasta_path)
    path_get_uniqueness = Path(
        "src", "konezumiaid", "get_rtpcr_primer", "get_uniqueness.sh"
    )
    subprocess.run(["bash", path_get_uniqueness], check=True)
    candidate_pairs = add_uniqueness(
        intron_len_added, miss_0_path, miss_1_path, miss_2_path
    )
    return candidate_pairs

