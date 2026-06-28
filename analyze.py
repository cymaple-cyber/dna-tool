"""命令行演示：分析一个 FASTA 文件中的每条序列。

用法:
    python analyze.py data/example.fasta
    python analyze.py data/example.fasta --kmer 3
"""

import argparse

from dna_toolkit import (
    gc_content,
    most_common_kmers,
    nucleotide_counts,
    read_fasta,
    reverse_complement,
    translate,
)


def main():
    parser = argparse.ArgumentParser(description="分析 FASTA 文件中的 DNA 序列")
    parser.add_argument("fasta", help="输入的 FASTA 文件路径")
    parser.add_argument("--kmer", type=int, default=3, help="统计的 k-mer 长度（默认 3）")
    args = parser.parse_args()

    for header, seq in read_fasta(args.fasta):
        counts = nucleotide_counts(seq)
        print(f"\n>>> {header}")
        print(f"  长度        : {len(seq)} bp")
        print(f"  碱基组成    : {counts}")
        print(f"  GC 含量     : {gc_content(seq):.1f}%")
        print(f"  反向互补    : {reverse_complement(seq)}")
        print(f"  翻译产物    : {translate(seq) or '(空)'}")
        top = most_common_kmers(seq, args.kmer, n=3)
        top_str = ", ".join(f"{kmer}×{n}" for kmer, n in top)
        print(f"  高频{args.kmer}-mer : {top_str}")


if __name__ == "__main__":
    main()
