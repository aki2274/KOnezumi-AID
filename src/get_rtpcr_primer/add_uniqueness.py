from __future__ import annotations
from pathlib import Path


def read_uniqueness(file_path) -> list[list]:
    with open(file_path, "r") as file:
        uniq_counts = [line.strip().split() for line in file]
    return uniq_counts


def add_uniqueness(
    candidate: list[dict],
    miss_0_path: Path,
    miss_1_path: Path,
    miss_2_path: Path,
) -> list[dict]:
    # import primer uniqueness from files
    miss_0 = read_uniqueness(miss_0_path)
    miss_1 = read_uniqueness(miss_1_path)
    miss_2 = read_uniqueness(miss_2_path)
    # add uniqueness info to candidate
    return_candidate = candidate.copy()
    for primers in return_candidate:
        left_primer = primers["left"]
        right_primer = primers["right"]

        left_0_uniq = next((x[0] for x in miss_0 if x[1] == left_primer), 0)
        right_0_uniq = next((x[0] for x in miss_0 if x[1] == right_primer), 0)
        left_1_uniq = next((x[0] for x in miss_1 if x[1] == left_primer), 0)
        right_1_uniq = next((x[0] for x in miss_1 if x[1] == right_primer), 0)
        left_2_uniq = next((x[0] for x in miss_2 if x[1] == left_primer), 0)
        right_2_uniq = next((x[0] for x in miss_2 if x[1] == right_primer), 0)

        primers["left_0_uniq"] = left_0_uniq
        primers["right_0_uniq"] = right_0_uniq
        primers["left_1_uniq"] = left_1_uniq
        primers["right_1_uniq"] = right_1_uniq
        primers["left_2_uniq"] = left_2_uniq
        primers["right_2_uniq"] = right_2_uniq

    return return_candidate
