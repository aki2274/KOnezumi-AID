from __future__ import annotations
import subprocess


# test bowtie's working
def read_sam(file_path: str):
    with open(file_path, "r") as f:
        sam = f.readlines()
    return sam


def test_get_uniqueness():
    test_shell = "bowtie tests/data/indexes/e_coli tests/data/reads/e_coli_1000.fq > /tmp/bowtie_test.sam"

    process = subprocess.Popen(
        ["bash", "-c", test_shell], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    _, stderr = process.communicate()

    # 終了コードが成功(0)でない場合は例外を発生させ、エラーメッセージを表示する
    exit_code = process.returncode
    if exit_code != 0:
        raise Exception(stderr.decode())

    expected_sam = read_sam("tests/data/uniq/e_coli_1000.sam")
    assert read_sam("/tmp/bowtie_test.sam") == expected_sam
