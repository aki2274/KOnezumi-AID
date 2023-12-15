from __future__ import annotations


def export_fasta(data: list[dict]) -> None:
    fasta_file = "data/reads/candidateprimer.fa"

    with open(fasta_file, "w") as file:
        for num, item in enumerate(data, start=1):
            left_primer = item["left_primer"]
            right_primer = item["right_primer"]

            file.write(f">L{num}\n{left_primer}\n")
            file.write(f">R{num}\n{right_primer}\n")
