from __future__ import annotations
import pandas as pd
import pickle
import argparse
from pathlib import Path
from konezumiaid.create_gene_dataclass import GeneData

# from konezumiaid.format_and_export_dataset.main import execute_export
from konezumiaid.create_gene_dataclass import create_dataclass
from konezumiaid.nominate_ptc_guide.main import nominate_candidate_stopcodon
from konezumiaid.nominate_splicesite_guide.search_candidate import search_site_candidate
from konezumiaid.apply_nmd_rules.main import apply_nmd_rules

# from konezumiaid.get_rtpcr_primer.main import export_primers

parser = argparse.ArgumentParser(
    description="This is KonezumiAID. A software to automate the design of gRNA for multiplex KO mouse using Target-AID"
)
parser.add_argument("-n", "--gene_name", type=str, help="Gene name or transcript name you want to.")

subparsers = parser.add_subparsers(dest="subcommand")

parser_create = subparsers.add_parser("create", help="Format and export the dataset as pickle files.")
parser_create.add_argument("refflat_path", type=Path, help="Path to the refFlat text file.")
parser_create.add_argument("chromosome_fasta_path", type=Path, help="Path to the chromosome fasta file.(ex. mm39.fa)")


args = parser.parse_args()


def extract_matching_seqs(*lists):
    seq_sets = [set(d["seq"] for d in lst) for lst in lists]
    common_seqs = set.intersection(*seq_sets)
    result = [{"seq": seq} for seq in common_seqs]
    return result


def show_table(
    df_ptc_gRNA: pd.DataFrame,
    df_acceptor_cand: pd.DataFrame,
    df_donor_cand: pd.DataFrame,
) -> None:
    print("PTC gRNA")
    if df_ptc_gRNA.empty:
        print("No PTC gRNA found.")
    else:
        print(df_ptc_gRNA.to_string())
    print("Acceptor gRNA")
    if df_acceptor_cand.empty:
        print("No Acceptor gRNA found.")
    else:
        print(df_acceptor_cand.to_string())
    print("Donor gRNA")
    if df_donor_cand.empty:
        print("No Donor gRNA found.")
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
    df_ptc_gRNA.to_csv(output_folder / f"{name}_ptc_gRNA.csv", index=False)
    df_acceptor_cand.to_csv(output_folder / f"{name}_acceptor_cand.csv", index=False)
    df_donor_cand.to_csv(output_folder / f"{name}_donor_cand.csv", index=False)


def konezumiaid_main(
    transcript_record: GeneData,
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
        df_ptcp_cand = pd.DataFrame(ptc_cand)
        df_acceptor_cand = pd.DataFrame(acceptor_cand)
        df_donor_cand = pd.DataFrame(donor_cand)
        show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand)
        export_csv(name, df_ptcp_cand, df_acceptor_cand, df_donor_cand)
    else:
        try:
            df_ref = pd.DataFrame(refflat_dic)
            df_symbol = df_ref[df_ref["geneName"] == name]
        except KeyError:
            raise ValueError(f"Gene name {name} not found in the dataset.")
        symbol_transcript_names = df_symbol["name"]
        for i, transcript_name in enumerate(symbol_transcript_names):
            print(f"Processing {transcript_name}...")
            transcript_record = create_dataclass(transcript_name, refflat_dic, seq_dict)
            ptc_cand, acceptor_cand, donor_cand = konezumiaid_main(transcript_record)
            if i == 0:
                symbol_ptc_cand = ptc_cand
                symbol_acceptor_cand = acceptor_cand
                symbol_donor_cand = donor_cand
            else:
                symbol_ptc_cand = extract_matching_seqs(symbol_ptc_cand, ptc_cand)
                symbol_acceptor_cand = extract_matching_seqs(symbol_acceptor_cand, acceptor_cand)
                symbol_donor_cand = extract_matching_seqs(symbol_donor_cand, donor_cand)
        df_ptcp_cand = pd.DataFrame(symbol_ptc_cand)
        df_acceptor_cand = pd.DataFrame(symbol_acceptor_cand)
        df_donor_cand = pd.DataFrame(symbol_donor_cand)
        show_table(df_ptcp_cand, df_acceptor_cand, df_donor_cand)
        export_csv(name, df_ptcp_cand, df_acceptor_cand, df_donor_cand)


def main():
    if args.subcommand == "create":
        execute_export(args.refflat_path, args.chromosome_fasta_path)
    else:
        gene_name = args.gene_name
        if gene_name is None:
            raise ValueError("Please provide a gene name.")
        execute(gene_name)
