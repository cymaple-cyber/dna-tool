"""FASTA 读取。"""


def read_fasta(path):
    """逐条产出 (header, sequence)，header 不含 '>'。"""
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
    return dict(read_fasta(path))
