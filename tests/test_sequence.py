"""核心功能的单元测试。运行: pytest 或 python -m pytest"""

import pytest

from dna_toolkit import (
    complement,
    gc_content,
    kmer_counts,
    nucleotide_counts,
    reverse_complement,
    transcribe,
    translate,
    validate_dna,
)


def test_validate_dna_ok():
    assert validate_dna("atgc") == "ATGC"


def test_validate_dna_rejects_invalid():
    with pytest.raises(ValueError):
        validate_dna("ATGX")


def test_nucleotide_counts():
    assert nucleotide_counts("AATGC") == {"A": 2, "C": 1, "G": 1, "T": 1}


def test_gc_content():
    assert gc_content("GGCC") == 100.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("ATGC") == 50.0


def test_gc_content_empty():
    assert gc_content("") == 0.0


def test_complement():
    assert complement("ATGC") == "TACG"


def test_reverse_complement():
    # 经典例子：reverse complement of AAAACCCGGT 是 ACCGGGTTTT
    assert reverse_complement("AAAACCCGGT") == "ACCGGGTTTT"


def test_transcribe():
    assert transcribe("ATGC") == "AUGC"


def test_translate_basic():
    # ATG=M, GCC=A, ATT=I, TGA=终止
    assert translate("ATGGCCATTTGA") == "MAI"


def test_translate_ignores_trailing():
    # 末尾不足 3 个碱基被忽略
    assert translate("ATGGC") == "M"


def test_kmer_counts():
    counts = kmer_counts("AAATTT", 2)
    assert counts["AA"] == 2
    assert counts["TT"] == 2
    assert counts["AT"] == 1
