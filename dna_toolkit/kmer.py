"""k-mer 频率统计。

k-mer 指长度为 k 的子序列，是基因组组装、序列比对、
物种鉴定等很多算法的基础特征。
"""

from collections import Counter

from .sequence import validate_dna


def kmer_counts(seq, k):
    """统计序列中所有长度为 k 的子串出现次数，返回 Counter。"""
    if k <= 0:
        raise ValueError("k 必须为正整数")
    seq = validate_dna(seq)
    if k > len(seq):
        return Counter()
    return Counter(seq[i:i + k] for i in range(len(seq) - k + 1))


def most_common_kmers(seq, k, n=5):
    """返回出现频率最高的 n 个 k-mer，列表形式 [(kmer, count), ...]。"""
    return kmer_counts(seq, k).most_common(n)
