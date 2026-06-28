"""序列基本操作：互补、转录、翻译、GC 含量等。"""

from collections import Counter

_DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# 标准遗传密码子表（DNA 版，用 T）。* 为终止密码子。
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def clean(seq):
    """转大写并去掉空白。"""
    return "".join(seq.split()).upper()


def validate_dna(seq):
    """校验只含 A/T/G/C，否则抛 ValueError，返回规范化后的序列。"""
    seq = clean(seq)
    invalid = set(seq) - set(_DNA_COMPLEMENT)
    if invalid:
        raise ValueError(f"序列含有非法碱基: {sorted(invalid)}")
    return seq


def nucleotide_counts(seq):
    seq = validate_dna(seq)
    counts = Counter(seq)
    return {base: counts.get(base, 0) for base in "ACGT"}


def gc_content(seq):
    """GC 含量，百分比。"""
    seq = validate_dna(seq)
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100


def complement(seq):
    seq = validate_dna(seq)
    return "".join(_DNA_COMPLEMENT[base] for base in seq)


def reverse_complement(seq):
    return complement(seq)[::-1]


def transcribe(seq):
    """DNA -> mRNA，T 换成 U。"""
    seq = validate_dna(seq)
    return seq.replace("T", "U")


def translate(seq, to_stop=True):
    """按密码子翻译成氨基酸序列。末尾不足 3 个碱基忽略。"""
    seq = validate_dna(seq)
    protein = []
    for i in range(0, len(seq) - 2, 3):
        amino_acid = CODON_TABLE[seq[i:i + 3]]
        if amino_acid == "*" and to_stop:
            break
        protein.append(amino_acid)
    return "".join(protein)
