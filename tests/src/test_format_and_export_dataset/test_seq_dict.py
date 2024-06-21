from src.konezumiaid.format_and_export_dataset.generate_seq_dict_from_fasta import (
    read_fasta,
    create_strand_plus_seq_dict,
)

from tempfile import NamedTemporaryFile
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


def test_create_strand_plus_seq_dict():
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
    seq_dict = {
        "NM_001011874": "ATGC",
        "NM_001355712": "TTGGCC",
    }
    expected = {
        "NM_001011874": "GCAT",
        "NM_001355712": "TTGGCC",
    }
    assert create_strand_plus_seq_dict(pd.DataFrame(gene_data), seq_dict) == expected
