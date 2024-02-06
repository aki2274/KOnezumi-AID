from __future__ import annotations
import re
from src.create_gene_dataclass import GeneData
from src.nominate_candidate_stopcodon.generate_cds_seq import generate_cdsseq

#####
# check grna position
#####


def label_in_start_150bp(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    result = cand_grna.copy()
    for grna in result:
        if "ct_seq" in grna:
            if re.search(grna["ct_seq"], ds.orf_seq).start() + 1 <= ds.cdsStart + 150:
                grna["start_150"] = True
            else:
                grna["start_150"] = False
        elif "ga_seq" in grna:
            if re.search(grna["ga_seq"], ds.orf_seq).start() + 16 <= ds.cdsStart + 150:
                grna["start_150"] = True
            else:
                grna["start_150"] = False
    return result


"""
# the case of 1 exon
# Targeting exons whose sequences differ by more than half of the coding region
def label_in_front_half(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    result = cand_grna.copy()
    cds = generate_cdsseq(ds)
    for grna in result:
        if "ct_seq" in grna:
            if re.search(grna["ct_seq"], cds).start() + 3 > len(cds) / 2:
                grna["back_half"] = True
            else:
                grna["back_half"] = False
        elif "ga_seq" in grna:
            if re.search(grna["ga_seq"], cds).start() + 19 > len(cds) / 2:
                grna["back_half"] = True
            else:
                grna["back_half"] = False
    return result
"""


# the case of multi exon
def eliminate_in_last_exon(cand_grna: list[str], ds: GeneData) -> list[str]:
    result = []
    for grna in cand_grna:
        if "ct_seq" in grna:
            if (
                re.search(grna["ct_seq"], ds.orf_seq).start() + 1
                <= ds.exon_start_list[-1]
            ):
                result.append(grna)
        elif "ga_seq" in grna:
            if (
                re.search(grna["ga_seq"], ds.orf_seq).start() + 16
                <= ds.exon_start_list[-1]
            ):
                result.append(grna)
    return result


def label_in_50bp_from_LEJ(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    result = cand_grna.copy()
    for grna in result:
        if "ct_seq" in grna:
            if (
                re.search(grna["ct_seq"], ds.orf_seq).start() + 1
                <= ds.exon_end_list[-2] - 50
            ):
                grna["50bp_from_LEJ"] = True
            else:
                grna["50bp_from_LEJ"] = False
        elif "ga_seq" in grna:
            if (
                re.search(grna["ga_seq"], ds.orf_seq).start() + 16
                <= ds.exon_end_list[-2] - 50
            ):
                grna["50bp_from_LEJ"] = True
            else:
                grna["50bp_from_LEJ"] = False
    return result
