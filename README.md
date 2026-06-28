# dna-toolkit

一个轻量的 **DNA 序列分析工具包**，纯 Python 实现、零第三方依赖，作为计算生物学入门项目。

## 功能

- 读取 FASTA 格式文件
- 碱基组成统计与 GC 含量计算
- 互补链 / 反向互补链
- DNA → mRNA 转录
- DNA → 蛋白质翻译（标准遗传密码子表）
- k-mer 频率统计

## 安装

```bash
git clone <你的仓库地址>
cd dna-toolkit
pip install -e ".[dev]"   # 可编辑安装，并附带测试依赖
```

> 也可以不安装，直接在项目根目录运行脚本——代码本身不依赖任何第三方库。

## 快速上手

命令行分析一个 FASTA 文件：

```bash
python analyze.py data/example.fasta
python analyze.py data/example.fasta --kmer 3
```

在 Python 中作为库使用：

```python
from dna_toolkit import gc_content, reverse_complement, translate

seq = "ATGGCCATTGTAATGGGCCGCTGA"
print(gc_content(seq))         # GC 含量百分比
print(reverse_complement(seq)) # 反向互补链
print(translate(seq))          # 翻译为氨基酸序列
```

## 项目结构

```
dna-toolkit/
├── dna_toolkit/        # 核心代码包
│   ├── sequence.py     # 互补、转录、翻译、GC 含量
│   ├── kmer.py         # k-mer 频率统计
│   └── io.py           # FASTA 读取
├── data/example.fasta  # 示例数据
├── tests/              # 单元测试
└── analyze.py          # 命令行演示脚本
```

## 运行测试

```bash
pytest
```

## 许可证

MIT
