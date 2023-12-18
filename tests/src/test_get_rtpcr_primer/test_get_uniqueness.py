import pytest
import subprocess

# test bowtie's working
test_shell = """
    bowtie tests/data/indexes/e_coli tests/data/reads/e_coli_1000.fq > /tmp/bowtie_test.sam
"""
subprocess.run(["bash", "-c", test_shell])


def read_sam(file_path: str):
    with open(file_path, "r") as f:
        sam = f.readlines()
    return sam


def test_get_uniqueness():
    expected_sam = read_sam("tests/data/uniq/e_coli_1000.sam")
    assert read_sam("/tmp/bowtie_test.sam") == expected_sam
