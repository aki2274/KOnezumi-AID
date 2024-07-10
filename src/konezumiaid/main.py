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
    list_of_dict: list[dict],
    transcript_record: TranscriptRecord,
    flag_ptc: bool = False,
) -> pd.DataFrame:
    if flag_ptc:
        target_aminoacids = [d["aminoacid"] for d in list_of_dict]
        if transcript_record.exon_count == 1:
            recomend = [not d["in_start_150bp"] for d in list_of_dict]
        else:
            recomend = [not any([d["in_start_150bp"], d["in_50bp_from_LEJ"]]) for d in list_of_dict]
    else:
        target_aminoacids = [None for _ in list_of_dict]
        recomend = [None for _ in list_of_dict]
    result_dict = [
        {
            "Target sequence (20mer + PAM)": d["seq"],
            "Recommended": r,
            "Target amino acid": aa,
            "link to CRISPRdirect": d["link_to_crisprdirect"],
        }
        for d, aa, r in zip(list_of_dict, target_aminoacids, recomend)
    ]
    return pd.DataFrame(result_dict)


def extract_matching_seqs(lists, flag_ptc=False) -> list[dict]:
    seq_sets = [set(d["seq"] for d in lst) for lst in lists]
    common_seqs = set.intersection(*seq_sets)
    common_seqs = [d["seq"] for d in lists[0] if d["seq"] in common_seqs]
    link_to_crisprdirect = [d["link_to_crisprdirect"] for d in lists[0] if d["seq"] in common_seqs]
    if flag_ptc:
        target_aminoacids = [d["aminoacid"] for d in lists[0] if d["seq"] in common_seqs]
        target_aminoacids = [aa[-1:] for aa in target_aminoacids]
    else:
        target_aminoacids = [None for d in lists[0] if d["seq"] in common_seqs]
    result = [
        {"Target sequence (20mer + PAM)": seq, "Target amino acid": aa, "link to CRISPRdirect": link}
        for seq, link, aa in zip(common_seqs, link_to_crisprdirect, target_aminoacids)
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
    df_ptc_gRNA.to_csv(output_folder / f"{name}_ptc.csv", index=False)
    df_acceptor_cand.to_csv(output_folder / f"{name}_acceptor_site.csv", index=False)
    df_donor_cand.to_csv(output_folder / f"{name}_donor_site.csv", index=False)


def konezumiaid_main(
    transcript_record: TranscriptRecord,
) -> tuple[list[dict], list[dict], list[dict]]:
    ct_cand, ga_cand = nominate_candidate_stopcodon(transcript_record)
    applied_nmd_rules_gRNA = apply_nmd_rules(transcript_record, ct_cand, ga_cand)
    acceptor_cand, donor_cand = search_site_candidate(transcript_record)
    return applied_nmd_rules_gRNA, acceptor_cand, donor_cand


def execute(name: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    seq_path = Path("data", "sorted_seq_dict.pkl")
    refflat_dic = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    if name.startswith("NM_"):
        transcript_record = create_dataclass(name, refflat_dic, seq_dict)
        ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
        df_ptcp_cand = format_output(ptc_cand, transcript_record, flag_ptc=True)
        df_acceptor_cand = format_output(acceptor_cand, transcript_record)
        df_donor_cand = format_output(donor_cand, transcript_record)

        show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand)
        export_csv(name, df_ptcp_cand, df_acceptor_cand, df_donor_cand)
    else:
        try:
            df_ref = pd.DataFrame(refflat_dic)
            df_symbol = df_ref[df_ref["geneName"] == name]
        except Exception as e:
            raise Exception(f"Error: {e}")
        if df_symbol.empty:
            raise Exception(f"Gene name {name} not found in the dataset.")
        symbol_transcript_names = df_symbol["name"]
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
