"""核心序列分析功能：互补、转录、翻译、GC 含量等。

本模块只处理 DNA/RNA 的字符串表示，不依赖任何第三方库，
方便理解每一步生物学操作背后的计算逻辑。
"""

from collections import Counter

# DNA 碱基的互补配对关系：A-T, G-C
_DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# 标准遗传密码子表（DNA 版本，T 而非 U）。
# 每 3 个碱基（一个密码子）对应一个氨基酸，"*" 表示终止密码子。
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
    """规范化序列：转大写并去除空白字符。"""
    return "".join(seq.split()).upper()


def validate_dna(seq):
    """校验是否为合法 DNA 序列（只含 A/T/G/C），非法时抛出 ValueError。"""
    seq = clean(seq)
    invalid = set(seq) - set(_DNA_COMPLEMENT)
    if invalid:
        raise ValueError(f"序列含有非法碱基: {sorted(invalid)}")
    return seq


def nucleotide_counts(seq):
    """统计四种碱基的数量，返回有序字典 {A, C, G, T}。"""
    seq = validate_dna(seq)
    counts = Counter(seq)
    return {base: counts.get(base, 0) for base in "ACGT"}


def gc_content(seq):
    """计算 GC 含量（百分比）。GC 含量是衡量序列稳定性的常用指标。"""
    seq = validate_dna(seq)
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100


def complement(seq):
    """返回互补链（不反向）。"""
    seq = validate_dna(seq)
    return "".join(_DNA_COMPLEMENT[base] for base in seq)


def reverse_complement(seq):
    """返回反向互补链——双链 DNA 中与给定链配对的另一条链（5'->3' 方向）。"""
    return complement(seq)[::-1]


def transcribe(seq):
    """DNA 转录为 mRNA：将 T 替换为 U。"""
    seq = validate_dna(seq)
    return seq.replace("T", "U")


def translate(seq, to_stop=True):
    """将 DNA 编码链翻译为氨基酸序列（单字母表示）。

    以 3 个碱基为一组读取密码子；末尾不足 3 个的碱基被忽略。
    to_stop=True 时遇到终止密码子停止翻译。
    """
    seq = validate_dna(seq)
    protein = []
    for i in range(0, len(seq) - 2, 3):
        amino_acid = CODON_TABLE[seq[i:i + 3]]
        if amino_acid == "*" and to_stop:
            break
        protein.append(amino_acid)
    return "".join(protein)
