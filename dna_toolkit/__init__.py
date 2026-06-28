"""dna-toolkit: 一个轻量的 DNA 序列分析工具包。"""

from .io import read_fasta, read_fasta_dict
from .kmer import kmer_counts, most_common_kmers
from .sequence import (
    complement,
    gc_content,
    nucleotide_counts,
    reverse_complement,
    transcribe,
    translate,
    validate_dna,
)

__version__ = "0.1.0"

__all__ = [
    "read_fasta",
    "read_fasta_dict",
    "kmer_counts",
    "most_common_kmers",
    "complement",
    "gc_content",
    "nucleotide_counts",
    "reverse_complement",
    "transcribe",
    "translate",
    "validate_dna",
]
