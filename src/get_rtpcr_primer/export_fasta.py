from __future__ import annotations


def export_fasta(data: list[dict]) -> None:
    file_path = "data/reads/candidateprimer.fa"
    processed_names = set()

    with open(file_path, "w") as file:
        for item in data:
            left_primer = item["left_primer"]
            right_primer = item["right_primer"]

            left_primer_name = left_primer
            if left_primer_name not in processed_names:
                processed_names.add(left_primer_name)
                file.write(f">{left_primer_name}\n{left_primer}\n")
            right_primer_name = right_primer
            if right_primer_name not in processed_names:
                processed_names.add(right_primer_name)
                file.write(f">{right_primer_name}\n{right_primer}\n")
