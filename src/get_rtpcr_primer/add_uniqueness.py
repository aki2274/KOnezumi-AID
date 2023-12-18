def read_uniqueness(file_path) -> list[list]:
    with open(file_path, "r") as file:
        uniq_counts = [line.strip().split() for line in file]
    return uniq_counts


def add_uniqueness(
    candidate: list[dict],
    miss_0_path: str = "data/uniq/0_miss_counts.txt",
    miss_1_path: str = "data/uniq/1_miss_counts.txt",
    miss_2_path: str = "data/uniq/2_miss_counts.txt",
) -> list[dict]:
    miss_0 = read_uniqueness(miss_0_path)
    miss_1 = read_uniqueness(miss_1_path)
    miss_2 = read_uniqueness(miss_2_path)

    for primer in candidate:
        left_primer = primer["left_primer"]
        right_primer = primer["right_primer"]

        left_0_uniq = next((x[0] for x in miss_0 if x[1] == left_primer), 0)
        right_0_uniq = next((x[0] for x in miss_0 if x[1] == right_primer), 0)
        left_1_uniq = next((x[0] for x in miss_1 if x[1] == left_primer), 0)
        right_1_uniq = next((x[0] for x in miss_1 if x[1] == right_primer), 0)
        left_2_uniq = next((x[0] for x in miss_2 if x[1] == left_primer), 0)
        right_2_uniq = next((x[0] for x in miss_2 if x[1] == right_primer), 0)

        primer["left_0_uniq"] = left_0_uniq
        primer["right_0_uniq"] = right_0_uniq
        primer["left_1_uniq"] = left_1_uniq
        primer["right_1_uniq"] = right_1_uniq
        primer["left_2_uniq"] = left_2_uniq
        primer["right_2_uniq"] = right_2_uniq

    return candidate
