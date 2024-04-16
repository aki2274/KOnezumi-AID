from __future__ import annotations
from pathlib import Path

def export_fasta(
    candidate_primers: list[dict], output_path: Path
) -> None:
    """
    export candidate primer seq to fasta file, for exceuting bowtie
    example output file: index and seqs are equal.
        >ATGCATGC
        ATCGATCG
        >AAAAAAAAAAA
        AAAAAAAAAAA
    """
    processed_names = set()

    with open(output_path, "w") as file:
        for primer_pair in candidate_primers:
            left_primer = primer_pair["left_seq"]
            right_primer = primer_pair["right_seq"]

            left_primer_name = left_primer
            if left_primer_name not in processed_names:
                processed_names.add(left_primer_name)
                file.write(f">{left_primer_name}\n{left_primer}\n")
            right_primer_name = right_primer
            if right_primer_name not in processed_names:
                processed_names.add(right_primer_name)
                file.write(f">{right_primer_name}\n{right_primer}\n")
