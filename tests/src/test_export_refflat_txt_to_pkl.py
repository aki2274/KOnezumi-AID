from src import export_refflat_txt_to_pkl
import pickle


def test_read_refFlat():
    path_refFlat = "tests/data/example_refFlat.txt"
    test = export_refflat_txt_to_pkl.read_refFlat(path_refFlat)
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


def test_make_genedata_dataframe():
    path_refFlat = "tests/data/example_refFlat.txt"
    test_dataframe = export_refflat_txt_to_pkl.make_genedata_dataframe(path_refFlat)
    expected_df = expected = [
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


def test_export_genedata_pkl():
    path_refFlat = "tests/data/example_refFlat.txt"
    test_dataframe = export_refflat_txt_to_pkl.make_genedata_dataframe(path_refFlat)
    export_refflat_txt_to_pkl.export_genedata_pkl(test_dataframe, "/tmp/test.pkl")
    with open("/tmp/test.pkl", mode="rb") as pklfile:
        test = pickle.load(pklfile)
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
    assert test == expected
