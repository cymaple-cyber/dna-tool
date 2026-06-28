"""FASTA 文件读写。

FASTA 是生物信息学中最常见的序列格式：
每条序列以 '>' 开头的描述行起始，后跟一行或多行序列数据。
"""


def read_fasta(path):
    """读取 FASTA 文件，逐条产出 (header, sequence) 元组。

    header 不含开头的 '>'；多行序列会被拼接为一个字符串。
    使用生成器，可处理较大的文件而不必全部载入内存。
    """
    header = None
    chunks = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def read_fasta_dict(path):
    """读取 FASTA 文件为字典 {header: sequence}。"""
    return dict(read_fasta(path))
