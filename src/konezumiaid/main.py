from __future__ import annotations
import pandas as pd
import pickle
import argparse
import sys
import shutil
from pathlib import Path
import importlib.metadata
from konezumiaid.create_gene_dataclass import TranscriptRecord
from konezumiaid.format_and_export_dataset.main import execute_export
from konezumiaid.create_gene_dataclass import create_dataclass
from konezumiaid.nominate_ptc_guide.main import nominate_candidate_stopcodon
from konezumiaid.nominate_splicesite_guide.search_candidate import search_site_candidate
from konezumiaid.evaluate_grna.main import apply_nmd_rules
from konezumiaid.format_output import (
    flag_single_exon,
    format_single_transcript_result,
    extract_multiple_transcripts_match,
    show_table,
    export_csv,
)

# from konezumiaid.get_rtpcr_primer.main import export_primers

try:
    version = importlib.metadata.version("konezumiaid")
except importlib.metadata.PackageNotFoundError:
    version = "could not be retrieved"

parser = argparse.ArgumentParser(
    description="This is KonezumiAID. A software to automate the design of gRNA for multiplex KO mouse using Target-AID"
)
parser.add_argument("-n", "--name", type=str, help="Gene name or transcript name (Refseq ID) you want to.")
parser.add_argument("-v", "--version", action="version", version=f"KOnezumi-AID version {version}")

subparsers = parser.add_subparsers(dest="subcommand")

parser_preprocess = subparsers.add_parser("preprocess", help="Format and export the dataset as pickle files.")
parser_preprocess.add_argument("refflat_path", type=Path, help="Path to the refFlat text file.")
parser_preprocess.add_argument(
    "chromosome_fasta_path", type=Path, help="Path to the chromosome fasta file.(e.g. mm39.fa)"
)

parser_batch = subparsers.add_parser("batch", help="Batch processing for multiple genes.")
parser_batch.add_argument("-f", "--file", type=Path, help="Path to the gene CSV or Excel file", required=True)

args = parser.parse_args()


def konezumiaid_main(
    transcript_record: TranscriptRecord,
) -> tuple[list[dict], list[dict], list[dict]]:
    ct_cand, ga_cand = nominate_candidate_stopcodon(transcript_record)
    applied_nmd_rules_gRNA = apply_nmd_rules(transcript_record, ct_cand, ga_cand)
    acceptor_cand, donor_cand = search_site_candidate(transcript_record)
    return applied_nmd_rules_gRNA, acceptor_cand, donor_cand


def execute(input_name: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    name = input_name
    refflat_path = Path("konezumiaid_data", "unique_coding_refflat.pkl")
    seq_path = Path("konezumiaid_data", "refseq_to_transcribed_region.pkl")
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
            df_ptcp_cand = extract_multiple_transcripts_match(tmp_ptc_cand, flag_ptc=True)
            df_acceptor_cand = extract_multiple_transcripts_match(tmp_acceptor_cand)
            df_donor_cand = extract_multiple_transcripts_match(tmp_donor_cand)
            show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand)
            export_csv(name, df_ptcp_cand, df_acceptor_cand, df_donor_cand)
            return None

    transcript_record = create_dataclass(name, refflat_dic, seq_dict)
    ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
    df_ptcp_cand = format_single_transcript_result(ptc_cand, transcript_record, flag_ptc=True)
    df_acceptor_cand = format_single_transcript_result(acceptor_cand, transcript_record)
    df_donor_cand = format_single_transcript_result(donor_cand, transcript_record)
    flag_sinle_exon = flag_single_exon(transcript_record)

    show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand, flag_sinle_exon)
    export_csv(input_name, df_ptcp_cand, df_acceptor_cand, df_donor_cand, flag_sinle_exon)
    return None


def main():
    konezumiaid_dir = Path("konezumiaid_data")
    konezumiaid_dir.mkdir(parents=True, exist_ok=True)
    path_refFlat = Path("konezumiaid_data", "unique_coding_refflat.pkl")
    path_seq = Path("konezumiaid_data", "refseq_to_transcribed_region.pkl")
    if args.subcommand == "preprocess":
        if not shutil.which("bedtools"):
            raise Exception("bedtools is not installed. Please install bedtools before running this script.")
        if path_refFlat.exists() or path_seq.exists():
            print("Pre-processing is already done.")
            sys.exit(0)
        execute_export(args.refflat_path, args.chromosome_fasta_path)
    elif args.subcommand == "batch":
        if not args.file.exists():
            raise FileNotFoundError(f"{args.file} is not found.")
        if args.file.suffix not in [".csv", ".xlsx"]:
            raise ValueError("The file must be either CSV or Excel file.")
        df = pd.read_csv(args.file, header=None) if args.file.suffix == ".csv" else pd.read_excel(args.file, header=None)
        for gene_name in df[0]:
            print(f"Processing {gene_name}...")
            execute(gene_name)
    else:
        if not path_refFlat.exists() or not path_seq.exists():
            raise FileNotFoundError(
                "The dataset is not found. Please preprocess the dataset by running 'konezumiaid preprocess'."
            )
        gene_name = args.name
        if gene_name is None:
            raise ValueError("Please specify the gene name.")
        execute(gene_name)
