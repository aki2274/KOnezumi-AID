from __future__ import annotations
import sys

sys.path.append("..")
import pandas as pd
import ast
import pickle
from pathlib import Path
from src.create_gene_dataclass import GeneData
from src.create_gene_dataclass import create_dataclass
from src.nominate_candidate_stopcodon.main import nominate_candidate_stopcodon
from src.nominate_spliceside_guide.search_candidate import search_splice_candidate
from src.adjust_nmd_rules.main import adjust_nmd_rules
from src.get_rtpcr_primer.main import export_primers


def show_table(
    adj_gRNA: list[dict], acand: list[dict], dcand: list[dict], cprimer: list[dict]
) -> None:
    adjusted_gRNA_df = pd.DataFrame(adj_gRNA)
    acceptor_cand_df = pd.DataFrame(acand)
    donor_cand_df = pd.DataFrame(dcand)
    candidate_primers_df = pd.DataFrame(cprimer)
    print("PTC gRNA")
    print(adjusted_gRNA_df.to_string())
    print("Acceptor gRNA")
    print(acceptor_cand_df.to_string())
    print("Donor gRNA")
    print(donor_cand_df.to_string())
    print("Primer")
    print(candidate_primers_df.to_string())


def konezumi_main(
    ds: GeneData,
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    ct_acand, ga_acand = nominate_candidate_stopcodon(ds)
    adjusted_gRNA = adjust_nmd_rules(ds, ct_acand, ga_acand)
    acceptor_cand, donor_cand = search_splice_candidate(ds)
    candidate_primers = export_primers(ds)
    return adjusted_gRNA, acceptor_cand, donor_cand, candidate_primers


def find_common_dicts_by_key(key, lists):
    value_sets = [{d[key] for d in lst if key in d} for lst in lists]
    common_values = set.intersection(*value_sets)
    result = [{key: value} for value in common_values]
    return result


def excecute_main(name: str) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    refflat_path = Path("data", "refFlat_genedata_sorted.pkl")
    seq_path = Path("data", "sorted_seq_dict.pkl")
    refflat = pickle.load(open(refflat_path, "rb"))
    seq_dict = pickle.load(open(seq_path, "rb"))
    if name.startswith("NM_") or name.startswith("NR_"):
        ds = create_dataclass(name, refflat, seq_dict)
        stop_cand, acceptor_cand, donor_cand, candidate_primers = konezumi_main(ds)
        show_table(stop_cand, acceptor_cand, donor_cand, candidate_primers)
        # return konezumi_main(ds)
    """
    else:
        transcript_names = [entry.name for entry in refflat if entry.geneName == name]
        stop_temp = []
        acceptor_temp = []
        donor_temp = []
        primer_temp = []
        for gene_name in transcript_names:
            ds = create_dataclass(gene_name, refflat, seq_dict)
            adjusted_gRNA, acceptor_cand, donor_cand, candidate_primers = konezumi_main(
                ds
            )
            stop_temp.append(
                ast.literal_eval(adjusted_gRNA["ct_seq"].replace("'", '"'))
            )
            stop_temp.append(
                ast.literal_eval(adjusted_gRNA["ga_seq"].replace("'", '"'))
            )
            acceptor_temp.append(
                ast.literal_eval(acceptor_cand["seq"].replace("'", '"'))
            )
            donor_temp.append(ast.literal_eval(donor_cand["seq"].replace("'", '"')))
            primer_temp.append(
                ast.literal_eval(candidate_primers["primer"].replace("'", '"'))
            )
        common_ct_stop = find_common_dicts_by_key("ct_seq", stop_temp)
        common_ga_stop = find_common_dicts_by_key("ga_seq", stop_temp)
        common_stop = common_ct_stop + common_ga_stop
        common_acceptor = find_common_dicts_by_key("seq", acceptor_temp)
        common_donor = find_common_dicts_by_key("seq", donor_temp)
        common_primer = find_common_dicts_by_key("primer", primer_temp)
        show_table(common_stop, common_acceptor, common_donor, common_primer)
        # return common_stop, common_acceptor, common_donor, common_primer
"""


excecute_main("NM_001011874")
