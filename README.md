# dna-toolkit

DNA 序列分析的一些常用操作，纯 Python，没有第三方依赖。学习计算生物学时写的练习项目。

包含：FASTA 读取、碱基组成与 GC 含量、互补 / 反向互补、转录、翻译、k-mer 统计。

## 用法

命令行：

```bash
python analyze.py data/example.fasta
python analyze.py data/example.fasta --kmer 3
```

当作库用：

```python
from dna_toolkit import gc_content, reverse_complement, translate

seq = "ATGGCCATTGTAATGGGCCGCTGA"
gc_content(seq)          # 54.2
reverse_complement(seq)  # 'TCAGCGGCCCATTACAATGGCCAT'
translate(seq)           # 'MAIVMGR'
```

## 测试

```bash
pytest
```
