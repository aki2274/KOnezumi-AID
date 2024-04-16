from __future__ import annotations
from pathlib import Path


def read_uniqueness(file_path) -> list[list]:
    with open(file_path, "r") as file:
        uniq_counts = [line.strip().split() for line in file]
    return uniq_counts


def add_uniqueness(
    candidate_primers: list[dict],
    miss_0_path: Path,
    miss_1_path: Path,
    miss_2_path: Path,
) -> list[dict]:
    # import primer uniqueness from files
    # data format: [[count of matched in the genome: int, primer_seq: str],...]
    miss_0 = read_uniqueness(miss_0_path)
    miss_1 = read_uniqueness(miss_1_path)
    miss_2 = read_uniqueness(miss_2_path)
    # add uniqueness info to candidate
    return_candidates = candidate_primers.copy()
    for primer_pair in return_candidates:
        left_primer = primer_pair["left_seq"]
        right_primer = primer_pair["right_seq"]

        left_0_mismatch = next((x[0] for x in miss_0 if x[1] == left_primer), 0)
        right_0_mismatch = next((x[0] for x in miss_0 if x[1] == right_primer), 0)
        left_1_mismatch = next((x[0] for x in miss_1 if x[1] == left_primer), 0)
        right_1_mismatch = next((x[0] for x in miss_1 if x[1] == right_primer), 0)
        left_2_mismatch = next((x[0] for x in miss_2 if x[1] == left_primer), 0)
        right_2_mismatch = next((x[0] for x in miss_2 if x[1] == right_primer), 0)

        primer_pair["left_0_mismatch"] = left_0_mismatch    
        primer_pair["left_1_mismatch"] = left_1_mismatch
        primer_pair["left_2_mismatch"] = left_2_mismatch
        primer_pair["right_0_mismatch"] = right_0_mismatch
        primer_pair["right_1_mismatch"] = right_1_mismatch
        primer_pair["right_2_mismatch"] = right_2_mismatch

    return return_candidates
