import pandas as pd
from konezumiaid.format_and_export_dataset.convert_refflat_to_bed6 import convert_refFlat_to_bed6


def test_convert_refFlat_to_bed6():
    column_names = [
        "geneName",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
    ]
    sample_refFlat = [
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
            "3284704,3491924,3740774,",
            "3287191,3492124,3741721,",
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
            "4190088,4212834,4218034,4218834,4234077,4240427,4267756,4276882,4296833,4298665,4301275,4313355,4313639,4313765,4315253,4331749,4337691,4354988,4363148,4381492,4422132,4422424,4430422,",
            "4190296,4212989,4218186,4218967,4234164,4240627,4267864,4277060,4297046,4298842,4301367,4313485,4313671,4313842,4315329,4331828,4337843,4355121,4363235,4381656,4422304,4423060,4430526,",
        ],
        [
            "Rp1",
            "NM_001195662",
            "chr1",
            "-",
            "4361068",
            "4479464",
            "4363203",
            "4479410",
            "4",
            "4361068,4422132,4422424,4479392,",
            "4363235,4422304,4423060,4479464,",
        ],
    ]
    df_refflat = pd.DataFrame(sample_refFlat, columns=column_names)
    test_path = "tests/data/test_convert_refFlat_to_bed6.bed"
    expected_output = "chr1\t3284704\t3741721\tNM_001011874\t0\t-\nchr1\t4190088\t4430526\tNM_001370921\t0\t-\nchr1\t4361068\t4479464\tNM_001195662\t0\t-\n"

    convert_refFlat_to_bed6(df_refflat, test_path)

    with open(test_path, "r") as f:
        assert f.read() == expected_output
