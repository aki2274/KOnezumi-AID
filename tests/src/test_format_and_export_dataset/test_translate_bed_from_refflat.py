import subprocess
import pytest
from pathlib import Path


def test_bedtools_getfasta():
    tmp_path = Path("tests/data/tmp")
    # Setup test files
    expected_output = """>feature1
    TAGC
    >feature2
    CTAC"""

    tmp_fasta = Path(tmp_path, "test.fa")
    tmp_bed = Path(tmp_path, "test.bed")
    output_file = Path("tests", "data", "bedtools.fa")
    # Run the script
    test_shell = "bedtools getfasta -nameOnly -fi $1 -bed $2 -fo $3"
    process = subprocess.Popen(
        ["bash", "-c", test_shell, "_", str(tmp_fasta), str(tmp_bed), str(output_file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _, stderr = process.communicate()

    # 終了コードが成功(0)でない場合は例外を発生させ、エラーメッセージを表示する
    exit_code = process.returncode
    if exit_code != 0:
        raise Exception(stderr.decode())

    # Check the output
    with open(output_file, "r") as f:
        output = f.read()
    assert output == expected_output
