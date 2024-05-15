from src.konezumiaid.export_dataset_as_pkl.generate_seq_dict_from_fasta import (
    read_fasta,
    create_dict_keys,
    create_strand_plus_seq_dict,
)
from src.konezumiaid.export_dataset_as_pkl.generate_sorted_genedata_from_refflat import (
    built_gene_dataframe,
    clean_refflat,
)
from tempfile import NamedTemporaryFile
from pathlib import Path
import pandas as pd


def test_read_fasta():
    # Create a temporary FASTA file for testing
    fasta_content = ">gene1\nATGC\n>gene2\nTTGGCC\n>gene3\nAAATTT\nCCCGGG"
    with NamedTemporaryFile(mode="w", delete=False) as temp_file:
        temp_file.write(fasta_content)
        temp_file_name = temp_file.name

    # Expected result
    expected = {"gene1": "ATGC", "gene2": "TTGGCC", "gene3": "AAATTTCCCGGG"}

    # Test the function
    assert read_fasta(temp_file_name) == expected


def test_create_dict_keys():
    path_refFlat = Path("tests/data/example_refFlat.txt")
    genedata = built_gene_dataframe(path_refFlat)
    expected = [
        "NM_001011874::chr1:3284704-3741721",
        "NM_001370921::chr1:4190088-4430526",
    ]
    assert create_dict_keys(genedata) == expected


def test_create_sorted_seq_dict():
    gene_data = [
        {
            "geneName": "Xkr4",
            "name": "NM_001011874",
            "chrom": "chr1",
            "strand": "-",
            "txStart": "10",
            "txEnd": "110",
            "cdsStart": "10",
            "cdsEnd": "110",
            "exonCount": "3",
            "exonStarts": "60,80,100",
            "exonEnds": "70,90,110",
        },
        {
            "geneName": "Lypla1",
            "name": "NM_001355712",
            "chrom": "chr1",
            "strand": "+",
            "txStart": "10",
            "txEnd": "110",
            "cdsStart": "10",
            "cdsEnd": "110",
            "exonCount": "4",
            "exonStarts": "10,20,30,40",
            "exonEnds": "15,25,35,45",
        },
    ]
    sort_gene_data = [
        {
            "geneName": "Xkr4",
            "name": "NM_001011874",
            "chrom": "chr1",
            "strand": "-",
            "txStart": 0,
            "txEnd": 100,
            "cdsStart": 0,
            "cdsEnd": 100,
            "exonCount": "3",
            "exonStarts": "0,20,40",
            "exonEnds": "10,30,50",
        },
        {
            "geneName": "Lypla1",
            "name": "NM_001355712",
            "chrom": "chr1",
            "strand": "+",
            "txStart": 0,
            "txEnd": 100,
            "cdsStart": 0,
            "cdsEnd": 100,
            "exonCount": "4",
            "exonStarts": "0,10,20,30",
            "exonEnds": "5,15,25,35",
        },
    ]
    seq_dict = {
        "NM_001011874::chr1:10-110": "ATGC",
        "NM_001355712::chr1:10-110": "TTGGCC",
    }
    expected = {
        "NM_001011874::chr1:0-100": "GCAT",
        "NM_001355712::chr1:0-100": "TTGGCC",
    }
    assert (
        create_strand_plus_seq_dict(
            pd.DataFrame(gene_data), pd.DataFrame(sort_gene_data), seq_dict
        )
        == expected
    )
