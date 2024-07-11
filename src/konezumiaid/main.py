from __future__ import annotations
import pandas as pd
import pickle
import argparse
import sys
import shutil
from pathlib import Path
from konezumiaid.create_gene_dataclass import TranscriptRecord
from konezumiaid.format_and_export_dataset.main import execute_export
from konezumiaid.create_gene_dataclass import create_dataclass
from konezumiaid.nominate_ptc_guide.main import nominate_candidate_stopcodon
from konezumiaid.nominate_splicesite_guide.search_candidate import search_site_candidate
from konezumiaid.evaluate_grna.main import apply_nmd_rules

# from konezumiaid.get_rtpcr_primer.main import export_primers

parser = argparse.ArgumentParser(
    description="This is KonezumiAID. A software to automate the design of gRNA for multiplex KO mouse using Target-AID"
)
parser.add_argument("-n", "--name", type=str, help="Gene name or transcript name (Refseq ID) you want to.")

subparsers = parser.add_subparsers(dest="subcommand")

parser_preprocess = subparsers.add_parser("preprocess", help="Format and export the dataset as pickle files.")
parser_preprocess.add_argument("refflat_path", type=Path, help="Path to the refFlat text file.")
parser_preprocess.add_argument(
    "chromosome_fasta_path", type=Path, help="Path to the chromosome fasta file.(e.g. mm39.fa)"
)


args = parser.parse_args()


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


def konezumiaid_main(
    transcript_record: TranscriptRecord,
) -> tuple[list[dict], list[dict], list[dict]]:
    ct_cand, ga_cand = nominate_candidate_stopcodon(transcript_record)
    applied_nmd_rules_gRNA = apply_nmd_rules(transcript_record, ct_cand, ga_cand)
    acceptor_cand, donor_cand = search_site_candidate(transcript_record)
    return applied_nmd_rules_gRNA, acceptor_cand, donor_cand


def execute(input_name: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    name = input_name
    refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    seq_path = Path("data", "sorted_seq_dict.pkl")
    refflat_dic = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    if not name.startswith("NM_"):
        try:
            df_ref = pd.DataFrame(refflat_dic)
            df_symbol = df_ref[df_ref["geneName"] == name]
        except Exception as e:
            raise Exception(f"Error: {e}")
        if df_symbol.empty:
            raise Exception(f"Gene name {name} not found in the dataset.")
        symbol_transcript_names = df_symbol["name"]
        if len(symbol_transcript_names) == 1:
            name = symbol_transcript_names.values[0]
        else:
            tmp_ptc_cand = []
            tmp_acceptor_cand = []
            tmp_donor_cand = []
            for transcript_name in symbol_transcript_names:
                print(f"Processing {transcript_name}...")
                transcript_record = create_dataclass(transcript_name, refflat_dic, seq_dict)
                if transcript_record is None:
                    continue
                ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
                tmp_ptc_cand.append(ptc_cand)
                tmp_acceptor_cand.append(acceptor_cand)
                tmp_donor_cand.append(donor_cand)
            df_ptcp_cand = extract_matching_seqs(tmp_ptc_cand, flag_ptc=True)
            df_acceptor_cand = extract_matching_seqs(tmp_acceptor_cand)
            df_donor_cand = extract_matching_seqs(tmp_donor_cand)
            show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand)
            export_csv(name, df_ptcp_cand, df_acceptor_cand, df_donor_cand)
            return None

    transcript_record = create_dataclass(name, refflat_dic, seq_dict)
    ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
    df_ptcp_cand = format_output(ptc_cand, transcript_record, flag_ptc=True)
    df_acceptor_cand = format_output(acceptor_cand, transcript_record)
    df_donor_cand = format_output(donor_cand, transcript_record)

    show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand)
    export_csv(input_name, df_ptcp_cand, df_acceptor_cand, df_donor_cand)
    return None


def main():
    path_refFlat = Path("data", "refFlat_genedata_sorted.pkl")
    path_seq = Path("data", "sorted_seq_dict.pkl")
    if args.subcommand == "preprocess":
        if not shutil.which("bedtools"):
            raise Exception("bedtools is not installed. Please install bedtools before running this script.")
        if path_refFlat.exists() or path_seq.exists():
            print("Pre-processing is already done.")
            sys.exit(0)
        execute_export(args.refflat_path, args.chromosome_fasta_path)
    else:
        if not path_refFlat.exists() or not path_seq.exists():
            raise FileNotFoundError(
                "The dataset is not found. Please preprocess the dataset by running 'konezumiaid preprocess'."
            )
        gene_name = args.name
        execute(gene_name)
