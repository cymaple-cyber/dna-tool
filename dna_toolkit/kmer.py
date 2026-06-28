"""k-mer 频率统计。"""

from collections import Counter

from .sequence import validate_dna


def kmer_counts(seq, k):
    if k <= 0:
        raise ValueError("k 必须为正整数")
    seq = validate_dna(seq)
    if k > len(seq):
        return Counter()
    return Counter(seq[i:i + k] for i in range(len(seq) - k + 1))


def most_common_kmers(seq, k, n=5):
    return kmer_counts(seq, k).most_common(n)
