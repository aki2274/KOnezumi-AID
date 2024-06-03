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
        left_primer = primer_pair["left"]
        right_primer = primer_pair["right"]

        for allowed_miss, matched_count in [(0, miss_0), (1, miss_1), (2, miss_2)]:
            left_mismatch = next((x[0] for x in matched_count if x[1] == left_primer), 0)
            right_mismatch = next((x[0] for x in matched_count if x[1] == right_primer), 0)

            primer_pair[f"left_{allowed_miss}_mismatch"] = left_mismatch
            primer_pair[f"right_{allowed_miss}_mismatch"] = right_mismatch

    return return_candidates
