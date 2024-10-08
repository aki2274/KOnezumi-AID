from pathlib import Path
import pandas as pd
from src.konezumiaid.format_and_export_dataset.generate_sorted_genedata_from_refflat import (
    read_refflat,
    built_gene_dataframe,
    clean_refflat,
    remove_transcript_duplicates,
    extract_autosomes_and_sex_chr,
)


def test_read_refFlat():
    path_refFlat = Path("tests/data/example_refFlat.txt")
    test = read_refflat(path_refFlat)
    expected = [
        [
            "Xkr4",
            "NM_001011874",
            "chr1",
            "-",
            "3284704",
            "3741721",
            "3286244",
            "3741571",
            "3",
            "3284704,3491924,3740774",
            "3287191,3492124,3741721",
        ],
        [
            "Rp1",
            "NM_001370921",
            "chr1",
            "-",
            "4190088",
            "4430526",
            "4190226",
            "4423048",
            "23",
            "4190088,4212834,4218034,4218834,4234077,4240427,4267756,4276882,4296833,4298665,4301275,4313355,4313639,4313765,4315253,4331749,4337691,4354988,4363148,4381492,4422132,4422424,4430422",
            "4190296,4212989,4218186,4218967,4234164,4240627,4267864,4277060,4297046,4298842,4301367,4313485,4313671,4313842,4315329,4331828,4337843,4355121,4363235,4381656,4422304,4423060,4430526",
        ],
    ]
    assert test == expected


def test_built_gene_dataframe():
    path_refFlat = Path("tests/data/example_refFlat.txt")
    test_dataframe = built_gene_dataframe(path_refFlat)
    expected = [
        {
            "geneName": "Xkr4",
            "name": "NM_001011874",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "3284704",
            "txEnd": "3741721",
            "cdsStart": "3286244",
            "cdsEnd": "3741571",
            "exonCount": "3",
            "exonStarts": "3284704,3491924,3740774",
            "exonEnds": "3287191,3492124,3741721",
        },
        {
            "geneName": "Rp1",
            "name": "NM_001370921",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "4190088",
            "txEnd": "4430526",
            "cdsStart": "4190226",
            "cdsEnd": "4423048",
            "exonCount": "23",
            "exonStarts": "4190088,4212834,4218034,4218834,4234077,4240427,4267756,4276882,4296833,4298665,4301275,4313355,4313639,4313765,4315253,4331749,4337691,4354988,4363148,4381492,4422132,4422424,4430422",
            "exonEnds": "4190296,4212989,4218186,4218967,4234164,4240627,4267864,4277060,4297046,4298842,4301367,4313485,4313671,4313842,4315329,4331828,4337843,4355121,4363235,4381656,4422304,4423060,4430526",
        },
    ]
    assert test_dataframe.to_dict(orient="records") == expected


def test_clean_refflat():
    test = [
        {
            "geneName": "Xkr4",
            "name": "NM_001011874",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "3284704",
            "txEnd": "3741721",
            "cdsStart": "3286244",
            "cdsEnd": "3741571",
            "exonCount": "3",
            "exonStarts": "3284704,3491924,3740774",
            "exonEnds": "3287191,3492124,3741721",
        },
        {
            "geneName": "Rp1",
            "name": "NM_001370921",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "4190088",
            "txEnd": "4430526",
            "cdsStart": "4190226",
            "cdsEnd": "4423048",
            "exonCount": "23",
            "exonStarts": "4190088,4212834,4218034,4218834,4234077,4240427,4267756,4276882,4296833,4298665,4301275,4313355,4313639,4313765,4315253,4331749,4337691,4354988,4363148,4381492,4422132,4422424,4430422",
            "exonEnds": "4190296,4212989,4218186,4218967,4234164,4240627,4267864,4277060,4297046,4298842,4301367,4313485,4313671,4313842,4315329,4331828,4337843,4355121,4363235,4381656,4422304,4423060,4430526",
        },
        {
            "geneName": "Lypla1",
            "name": "NM_001355712",
            "chrom": "chr1",
            "strand": "+",
            "txStart": "4878045",
            "txEnd": "4916958",
            "cdsStart": "4878136",
            "cdsEnd": "4915239",
            "exonCount": "8",
            "exonStarts": "4878045,4878677,4898806,4900490,4902533,4907223,4911178,4915185",
            "exonEnds": "4878205,4878709,4898872,4900538,4902604,4907297,4911355,4916958",
        },
    ]
    inpt_dataframe = pd.DataFrame(test)
    result_dataframe = clean_refflat(inpt_dataframe)
    excepted = [
        {
            "geneName": "Xkr4",
            "name": "NM_001011874",
            "chrom": "chr1",
            "strand": "-",
            "txStart": 0,
            "txEnd": 457017,
            "cdsStart": 150,
            "cdsEnd": 455477,
            "exonCount": "3",
            "exonStarts": "0,249597,454530",
            "exonEnds": "947,249797,457017",
        },
        {
            "geneName": "Rp1",
            "name": "NM_001370921",
            "chrom": "chr1",
            "strand": "-",
            "txStart": 0,
            "txEnd": 240438,
            "cdsStart": 7478,
            "cdsEnd": 240300,
            "exonCount": "23",
            "exonStarts": "0,7466,8222,48870,67291,75405,92683,98698,115197,116684,116855,117041,129159,131684,133480,153466,162662,189899,196362,211559,212340,217537,240230",
            "exonEnds": "104,8102,8394,49034,67378,75538,92835,98777,115273,116761,116887,117171,129251,131861,133693,153644,162770,190099,196449,211692,212492,217692,240438",
        },
        {
            "geneName": "Lypla1",
            "name": "NM_001355712",
            "chrom": "chr1",
            "strand": "+",
            "txStart": 0,
            "txEnd": 38913,
            "cdsStart": 91,
            "cdsEnd": 37194,
            "exonCount": "8",
            "exonStarts": "0,632,20761,22445,24488,29178,33133,37140",
            "exonEnds": "160,664,20827,22493,24559,29252,33310,38913",
        },
    ]
    assert result_dataframe.to_dict(orient="records") == excepted


def test_remove_transcript_duplicates():
    test = [
        {
            "geneName": "Xkr4",
            "name": "NM_001011874",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "3284704",
            "txEnd": "3741721",
            "cdsStart": "3286244",
            "cdsEnd": "3741571",
            "exonCount": "3",
            "exonStarts": "3284704,3491924,3740774",
            "exonEnds": "3287191,3492124,3741721",
        },
        {
            "geneName": "Rp1",
            "name": "NM_001370921",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "4190088",
            "txEnd": "4430526",
            "cdsStart": "4190226",
            "cdsEnd": "4423048",
            "exonCount": "23",
            "exonStarts": "4190088,4212834,4218034,4218834,4234077,4240427,4267756,4276882,4296833,4298665,4301275,4313355,4313639,4313765,4315253,4331749,4337691,4354988,4363148,4381492,4422132,4422424,4430422",
            "exonEnds": "4190296,4212989,4218186,4218967,4234164,4240627,4267864,4277060,4297046,4298842,4301367,4313485,4313671,4313842,4315329,4331828,4337843,4355121,4363235,4381656,4422304,4423060,4430526",
        },
        {
            "geneName": "Xkr4",
            "name": "NM_001011874",
            "chrom": "chr1",
            "strand": "+",
            "txStart": "4878045",
            "txEnd": "4916958",
            "cdsStart": "4878136",
            "cdsEnd": "4915239",
            "exonCount": "8",
            "exonStarts": "4878045,4878677,4898806,4900490,4902533,4907223,4911178,4915185",
            "exonEnds": "4878205,4878709,4898872,4900538,4902604,4907297,4911355,4916958",
        },
    ]
    inpt_dataframe = pd.DataFrame(test)
    result_dataframe = remove_transcript_duplicates(inpt_dataframe)
    excepted = [
        {
            "geneName": "Rp1",
            "name": "NM_001370921",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "4190088",
            "txEnd": "4430526",
            "cdsStart": "4190226",
            "cdsEnd": "4423048",
            "exonCount": "23",
            "exonStarts": "4190088,4212834,4218034,4218834,4234077,4240427,4267756,4276882,4296833,4298665,4301275,4313355,4313639,4313765,4315253,4331749,4337691,4354988,4363148,4381492,4422132,4422424,4430422",
            "exonEnds": "4190296,4212989,4218186,4218967,4234164,4240627,4267864,4277060,4297046,4298842,4301367,4313485,4313671,4313842,4315329,4331828,4337843,4355121,4363235,4381656,4422304,4423060,4430526",
        },
    ]
    assert result_dataframe.to_dict(orient="records") == excepted


def test_extract_autosomes_and_sex_chr():
    test_chr = ["chr1", "chr2", "chrX", "chrY", "chrM", "chr1_KI270706v1_random"]
    pd_df = pd.DataFrame({"chrom": test_chr})
    excepted = pd.DataFrame({"chrom": ["chr1", "chr2", "chrX", "chrY"]})
    result = extract_autosomes_and_sex_chr(pd_df)
    assert result.to_dict(orient="records") == excepted.to_dict(orient="records")
