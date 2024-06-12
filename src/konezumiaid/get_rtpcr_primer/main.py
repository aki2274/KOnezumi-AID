from __future__ import annotations
from pathlib import Path
import subprocess
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.get_exon_range import get_exon_range
from konezumiaid.nominate_ptc_guide.generate_cds_seq import generate_exon_seq
from konezumiaid.get_rtpcr_primer.make_rtpcr_primer import generate_candidate_info
from konezumiaid.get_rtpcr_primer.rate_quality import (
    is_crossing_juncion,
    add_intron_len,
)
from konezumiaid.get_rtpcr_primer.export_fasta import export_fasta
from konezumiaid.get_rtpcr_primer.add_uniqueness import add_uniqueness


def export_primers(transcript_record: GeneData) -> list[dict]:
    """
    This function exports candidate rt-qPCR primers based on the exon sequence of a given transcript record.
    It performs the following steps:
    1. Creates a directory named "uniq" inside the "data" directory if it doesn't exist.
    2. Defines file paths for storing mismatch counts.
    3. Retrieves the exon range of the transcript record.
    4. Generates the exon sequence.
    5. Generates candidate primer information based on the exon sequence and exon range.
    6. Determines if the candidate primers cross the junction.
    7. Adds the intron length to the candidate primers.
    8. Exports the candidate primers to a FASTA file.
    9. Executes a shell script to calculate the uniqueness of the candidate primers.
    10. Adds the uniqueness information to the candidate primers.
    11. Returns the list of candidate primer dictionaries.

    Args:
        transcript_record (GeneData): The Dataclass named GeneData.

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
            "left_0_mismatch": 1,
            "left_1_mismatch": 10,
            "left_2_mismatch": 25,
            "right_0_mismatch": 1,
            "right_1_mismatch": 16,
            "right_2_mismatch": 56
        }...]
    """
    # 1 Create a directory named "uniq" inside the "data" directory if it doesn't exist.
    Path("data", "uniq").mkdir(parents=True, exist_ok=True)
    # 2 Define file paths for storing mismatch counts.
    miss_0_path = Path("data", "uniq", "missmatch_0_counts.txt")
    miss_1_path = Path("data", "uniq", "missmatch_1_counts.txt")
    miss_2_path = Path("data", "uniq", "missmatch_2_counts.txt")
    # 3 Retrieve the exon range of the transcript record.
    exon_range = get_exon_range(transcript_record)
    # 4 Generate the exon sequence.
    exon_seq = generate_exon_seq(transcript_record)
    # 5 Generate candidate primer information based on the exon sequence and exon range.
    primer3_result = generate_candidate_info(exon_seq, exon_range)
    # 6 Determine if the candidate primers cross the junction.
    cross_validated = is_crossing_juncion(primer3_result, exon_range)
    # 7 Add the intron length to the candidate primers.
    intron_len_added = add_intron_len(cross_validated, transcript_record)
    # 8 Export the candidate primers to a FASTA file.
    fasta_path = Path("data", "uniq", "candidateprimer.fa")
    export_fasta(intron_len_added, fasta_path)
    # 9 Execute a shell script to calculate the uniqueness of the candidate primers.
    path_get_uniqueness = Path("src", "konezumiaid", "get_rtpcr_primer", "get_uniqueness.sh")
    subprocess.run(["bash", path_get_uniqueness], check=True)
    # 10 Add the uniqueness information to the candidate primers.
    candidate_pairs = add_uniqueness(intron_len_added, miss_0_path, miss_1_path, miss_2_path)
    return candidate_pairs
