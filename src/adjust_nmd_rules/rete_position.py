from __future__ import annotations
import re
from src.create_gene_dataclass import GeneData
from src.nominate_candidate_stopcodon.generate_cds_seq import generate_cdsseq

#####
# check grna position
#####


# the case of 1 exon


# Targeting exons whose sequences differ by more than half of the coding region
def label_in_back_half(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    result = cand_grna.copy()
    cds = generate_cdsseq(ds)
    for grna in result:
        if re.search(grna["seq"], cds).start() > len(cds) / 2:
            grna["back_half"] = True
        else:
            grna["back_half"] = False
    return result


# the case of multi exon
def eliminate_in_last_exon(cand_grna: list[str], ds: GeneData) -> list[str]:
    cand_grna = [
        grna
        for grna in cand_grna
        if re.search(grna, ds.orf_seq).start() <= ds.exon_end_list[-1]
    ]
    return cand_grna


def label_in_start_150bp(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    result = cand_grna.copy()
    for grna in result:
        if re.search(grna["seq"], ds.orf_seq).start() <= ds.cdsStart + 150:
            grna["start_150"] = True
        else:
            grna["start_150"] = False
    return result


def label_in_50bp_from_LEJ(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    result = cand_grna.copy()
    for grna in result:
        if re.search(grna["seq"], ds.orf_seq).start() >= ds.exon_end_list[-2] - 50:
            grna["end_50"] = True
        else:
            grna["end_50"] = False
    return result
