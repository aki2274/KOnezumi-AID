from __future__ import annotations
import pandas as pd
import pickle
import sys
from pathlib import Path
from konezumiaid.create_gene_dataclass import GeneData
from konezumiaid.format_and_export_dataset.main import export_pkl
from konezumiaid.create_gene_dataclass import create_dataclass
from konezumiaid.nominate_candidate_stopcodon.main import nominate_candidate_stopcodon
from konezumiaid.nominate_spliceside_guide.search_candidate import search_candidate
from konezumiaid.apply_nmd_rules.main import apply_nmd_rules
from konezumiaid.get_rtpcr_primer.main import export_primers


def show_table(
    adjusted_gRNA_df: pd.DataFrame,
    acceptor_cand_df: pd.DataFrame,
    donor_cand_df: pd.DataFrame,
    candidate_primers_df: pd.DataFrame,
) -> None:
    print("PTC gRNA")
    print(adjusted_gRNA_df.to_string())
    print("Acceptor gRNA")
    print(acceptor_cand_df.to_string())
    print("Donor gRNA")
    print(donor_cand_df.to_string())
    print("Primer")
    print(candidate_primers_df.to_string())


def konezumiaid_main(
    ds: GeneData,
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    ct_acand, ga_acand = nominate_candidate_stopcodon(ds)
    adjusted_gRNA = apply_nmd_rules(ds, ct_acand, ga_acand)
    acceptor_cand, donor_cand = search_candidate(ds)
    candidate_primers = export_primers(ds)
    return adjusted_gRNA, acceptor_cand, donor_cand, candidate_primers


def find_common_dicts_by_key(key, lists):
    value_sets = [{d[key] for d in lst if key in d} for lst in lists]
    common_values = set.intersection(*value_sets)
    result = [{key: value} for value in common_values]
    return result


def excecute(name: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    seq_path = Path("data", "sorted_seq_dict.pkl")
    refflat = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    if name.startswith("NM_") or name.startswith("NR_"):
        ds = create_dataclass(name, refflat, seq_dict)
        stop_cand, acceptor_cand, donor_cand, candidate_primers = konezumiaid_main(ds)
        stop_cand_df = pd.DataFrame(stop_cand)
        acceptor_cand_df = pd.DataFrame(acceptor_cand)
        donor_cand_df = pd.DataFrame(donor_cand)
        candidate_primers_df = pd.DataFrame(candidate_primers)
        show_table(stop_cand_df, acceptor_cand_df, donor_cand_df, candidate_primers_df)
        # return konezumiaid_main(ds)

    else:
        transcript_names = [d["name"] for d in refflat if d["geneName"] == name]
        gene_cand = {}
        for transcript in transcript_names:
            ds = create_dataclass(transcript, refflat, seq_dict)
            adjusted_gRNA, acceptor_cand, donor_cand, candidate_primers = (
                konezumiaid_main(ds)
            )
            stop_cand_df = pd.DataFrame(adjusted_gRNA)
            acceptor_cand_df = pd.DataFrame(acceptor_cand)
            donor_cand_df = pd.DataFrame(donor_cand)
            candidate_primers_df = pd.DataFrame(candidate_primers)
            gene_cand[transcript] = {
                "stop_cand": stop_cand_df,
                "acceptor_cand": acceptor_cand_df,
                "donor_cand": donor_cand_df,
                "candidate_primers": candidate_primers_df,
            }
        for k, v in gene_cand.items():
            print(k)
            show_table(
                v["stop_cand"],
                v["acceptor_cand"],
                v["donor_cand"],
                v["candidate_primers"],
            )


def export():
    if len(sys.argv) != 3:
        raise ValueError("Please provide a refflat file and a fasta file as arguments.")
    refflat_path = Path(sys.argv[1])
    fasta_path = Path(sys.argv[2])
    export_pkl(refflat_path, fasta_path)


def main():
    if len(sys.argv) != 2:
        raise ValueError("Please provide a gene name as an argument.")
    gene_name = sys.argv[1]
    excecute(gene_name)
