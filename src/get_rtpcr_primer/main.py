from __future__ import annotations
import pickle
from pathlib import Path
import subprocess
from src.create_gene_dataclass import GeneData, create_dataclass
from src.get_range_of_exon import get_exon_range
from src.nominate_candidate_stopcodon.generate_cds_seq import generate_exon_seq

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


def export_primers(
    transcript_name: str,
    refflat_path: str,
    seq_path: str,
) -> list[dict]:
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
    # 0. get gene data
    refflat = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    ds = create_dataclass(transcript_name, refflat, seq_dict)

    Path("data", "uniq").mkdir(parents=True, exist_ok=True)

    miss_0_path = Path("data", "uniq", "0_miss_counts.txt")
    miss_1_path = Path("data", "uniq", "1_miss_counts.txt")
    miss_2_path = Path("data", "uniq", "2_miss_counts.txt")

    # 1. get exon range
    exon_range = get_exon_range(ds)
    # 2. get exon seq
    exon_seq = generate_exon_seq(ds)
    # 3. get candidate primer
    primer3_result = export_candidate(exon_seq)
    # 4. get candidate primer info
    candidate_pairs = generate_candidate_info(exon_seq, primer3_result, exon_range)
    # 5. rate quality of candidate primer
    candidate_pairs = verify_crossing_exonjunction(candidate_pairs, exon_range)

    candidate_pairs = autocorrect_intron_len(candidate_pairs, ds)
    export_fasta(candidate_pairs)
    path_get_uniqueness = Path("src", "get_rtpcr_primer", "get_uniqueness.sh")
    subprocess.run([path_get_uniqueness])
    candidate_pairs = add_uniqueness(
        candidate_pairs, miss_0_path, miss_1_path, miss_2_path
    )
    return candidate_pairs
