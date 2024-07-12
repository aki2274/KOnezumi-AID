from __future__ import annotations
import pandas as pd
from pathlib import Path
from konezumiaid.create_gene_dataclass import TranscriptRecord


def format_output(
    candidate_primers: list[dict],
    transcript_record: TranscriptRecord,
    flag_ptc: bool = False,
) -> pd.DataFrame:
    if flag_ptc:
        target_aminoacids = [d["aminoacid"] for d in candidate_primers]
        if transcript_record.exon_count == 1:
            is_recomend = [not d["in_start_150bp"] for d in candidate_primers]
        else:
            is_recomend = [not any([d["in_start_150bp"], d["in_50bp_from_LEJ"]]) for d in candidate_primers]
        output_primers = [
            {
                "Target sequence (20mer + PAM)": d["seq"],
                "Recommended": r,
                "Target amino acid": aa,
                "link to CRISPRdirect": d["link_to_crisprdirect"],
            }
            for d, aa, r in zip(candidate_primers, target_aminoacids, is_recomend)
        ]
    else:
        exon_index = [d["exon_index"] for d in candidate_primers]
        output_primers = [
            {
                "Target sequence (20mer + PAM)": d["seq"],
                "Exon index": ei,
                "link to CRISPRdirect": d["link_to_crisprdirect"],
            }
            for d, ei in zip(candidate_primers, exon_index)
        ]
    return pd.DataFrame(output_primers)


def extract_matching_seqs(candidate_gene_primers: list[list], flag_ptc: bool = False) -> list[dict]:
    seq_sets = [set(d["seq"] for d in lst) for lst in candidate_gene_primers]
    common_seq = set.intersection(*seq_sets)
    common_data = [d for d in candidate_gene_primers[0] if d["seq"] in common_seq]
    if flag_ptc:
        result = [
            {
                "Target sequence (20mer + PAM)": d["seq"],
                "Target amino acid": d["aminoacid"][-1:],
                "link to CRISPRdirect": d["link_to_crisprdirect"],
            }
            for d in common_data
        ]
    else:
        result = [
            {
                "Target sequence (20mer + PAM)": d["seq"],
                "link to CRISPRdirect": d["link_to_crisprdirect"],
            }
            for d in common_data
        ]

    return pd.DataFrame(result)


def show_table(
    df_ptc_gRNA: pd.DataFrame,
    df_acceptor_cand: pd.DataFrame,
    df_donor_cand: pd.DataFrame,
) -> None:
    print("List of gRNAs to generate PTC (premature termination codon)")
    if df_ptc_gRNA.empty:
        print("No gRNA found.")
    else:
        print(df_ptc_gRNA.to_string())
    print("List of gRNAs to disrupt splice acceptor site")
    if df_acceptor_cand.empty:
        print("No gRNA found.")
    else:
        print(df_acceptor_cand.to_string())
    print("List of gRNAs to disrupt splice donor site")
    if df_donor_cand.empty:
        print("No gRNA found.")
    else:
        print(df_donor_cand.to_string())


def export_csv(
    name: str,
    df_ptc_gRNA: pd.DataFrame,
    df_acceptor_cand: pd.DataFrame,
    df_donor_cand: pd.DataFrame,
) -> None:
    output_folder = Path("data", "output")
    output_folder.mkdir(parents=True, exist_ok=True)
    df_splice_gRNA = pd.concat([df_acceptor_cand, df_donor_cand])
    if not df_ptc_gRNA.empty:
        df_ptc_gRNA.to_csv(output_folder / f"{name}_ptc_gRNA.csv", index=False)
    if not df_splice_gRNA.empty:
        df_splice_gRNA.to_csv(output_folder / f"{name}_splice_gRNA.csv", index=False)
