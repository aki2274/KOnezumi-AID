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


# the case of single exon
def label_in_front_half(cand_grna: list[dict], ds: GeneData) -> list[dict]:
    # check if the PTC is not in the front half of the CDS
    result = []
    cds = generate_cdsseq(ds)
    for grna in cand_grna:
        if "ct_seq" in grna:
            ptc_index = re.search(grna["ct_seq"], cds)
            if ptc_index is not None and ptc_index.start() + 3 < len(cds) / 2:
                result.append(grna)
        elif "ga_seq" in grna:
            ptc_index = re.search(grna["ga_seq"], cds)
            if ptc_index is not None and ptc_index.start() + 19 < len(cds) / 2:
                result.append(grna)
    return result


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
                >= ds.exon_end_list[-2] - 50
            ):
                grna["50bp_from_LEJ"] = True
            else:
                grna["50bp_from_LEJ"] = False
        elif "ga_seq" in grna:
            if (
                re.search(grna["ga_seq"], ds.orf_seq).start() + 16
                >= ds.exon_end_list[-2] - 50
            ):
                grna["50bp_from_LEJ"] = True
            else:
                grna["50bp_from_LEJ"] = False
    return result
