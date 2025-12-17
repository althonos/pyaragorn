# 👑 PyARAGORN [![Stars](https://img.shields.io/github/stars/althonos/pyaragorn.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyaragorn/stargazers)

*Cython bindings and Python interface to [ARAGORN](https://www.trna.se/), a (t|mt|tm)RNA gene finder*.

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyaragorn/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyaragorn/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyaragorn?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyaragorn/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pyaragorn.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyaragorn)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyaragorn?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyaragorn)
[![AUR](https://img.shields.io/aur/version/python-pyaragorn?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyaragorn)
[![Wheel](https://img.shields.io/pypi/wheel/pyaragorn.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyaragorn/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyaragorn.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyaragorn/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyaragorn.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyaragorn/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyaragorn/)
[![Mirror](https://img.shields.io/badge/mirror-LUMC-003EAA.svg?maxAge=2678400&style=flat-square)](https://git.lumc.nl/mflarralde/pyaragorn/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/pyaragorn.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyaragorn/issues)
[![Docs](https://img.shields.io/readthedocs/pyaragorn/latest?style=flat-square&maxAge=600)](https://pyaragorn.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyaragorn/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pyaragorn?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pyaragorn)


## 🗺️ Overview

[ARAGORN](https://trna.se) is a fast method developed
by Dean Laslett & Björn Canback[\[1\]](#ref1) to identify tRNA and tmRNA
genes in genomic sequences using heuristics to detect potential high-scoring
stem-loop structures. The complementary method ARWEN, developed by the same
authors[\[2\]](#ref2) to support the detection of metazoan mitochondrial
RNA (mtRNA) genes, was later integrated into ARAGORN.

`pyaragorn` is a Python module that provides bindings to ARAGORN and ARWEN
using [Cython](https://cython.org/). It directly interacts with the
ARAGORN internals, which has the following advantages:

- **single dependency**: PyARAGORN is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  ARAGORN binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you fully control, so you don't have to invoke the ARAGORN CLI using a
  sub-process and temporary files. Sequences can be passed directly as
  strings, bytes, or any buffer objects, which avoids the overhead of
  formatting your input to FASTA for ARAGORN.
- **no output parsing**: The detected RNA genes are returned as Python
  objects with transparent attributes, which facilitate handling the output
  of ARAGORN compared to parsing the output tables.
- **same results**: PyARAGORN is tested to ensure it produces the same results
  as ARAGORN `v1.2.41`, the latest release.


### 📋 Features

PyARAGORN currently supports the following features from the ARAGORN
command line:

- [x] tRNA gene detection (`aragorn -t`).
- [x] tmRNA gene detection (`aragorn -m`).
- [ ] mtRNA gene detection (`aragorn -mt`).
- [x] Reporting of batch mode metadata (`aragorn -w`).
- [x] Alternative genetic code (`aragorn -gc`).
- [ ] Custom genetic code (`aragorn -gc<n>,BBB=<aa>`).
- [x] Circular and linear topologies (`aragorn -c` | `aragorn -l`).
- [ ] Intron length configuration (`aragorn -i`).
- [x] Scoring threshold configuration (`aragorn -ps`).
- [x] Sequence extraction from RNA gene (`aragorn -seq`).
- [ ] Secondary structure extraction from each gene (`aragorn -br`).

### 🧶 Thread-safety

`pyaragorn.RNAFinder` instances are thread-safe. In addition, the `find_rna`
method is re-entrant. This means you can parameterize a `RNAFinder` instance
once, and then use a pool to process sequences in parallel:

```python
import multiprocessing.pool
import pyaragorn

rna_finder = pyaragorn.RNAFinder()

with multiprocessing.pool.ThreadPool() as pool:
    predictions = pool.map(rna_finder.find_rna, sequences)
```

## 🔧 Installing

This project is supported on Python 3.7 and later.

PyARAGORN can be installed directly from [PyPI](https://pypi.org/project/pyaragorn/),
which hosts some pre-built wheels for the x86-64 architecture (Linux/MacOS/Windows)
and the Aarch64 architecture (Linux/MacOS), as well as the code required to compile
from source with Cython:
```console
$ pip install pyaragorn
```

## 💡 Example

Let's load a sequence from a
[GenBank](http://www.insdc.org/files/feature_table.html) file,
use a `RNAFinder` to find all the tRNA genes it contains,
and print the anticodon and corresponding amino-acids of the detected
tRNAs.


### 🔬 [Biopython](https://github.com/biopython/biopython)

To use the `RNAFinder` to detect tRNA and tmRNA genes, the default operation
mode, but using the bacterial genetic code (translation table 11):

```python
import Bio.SeqIO
import pyaragorn

record = Bio.SeqIO.read("sequence.gbk", "genbank")

rna_finder = pyaragorn.RNAFinder(translation_table=11)
genes = rna_finder.find_rna(bytes(record.seq))

for gene in genes:
    if gene.type == "tRNA":
        print(
            gene.amino_acid,   # 3-letter code
            gene.begin,        # 1-based, inclusive
            gene.end,
            gene.strand,       # +1 or -1 for direct and reverse strand
            gene.energy,
            gene.anticodon
        )
```

*On older versions of Biopython (before 1.79) you will need to use
`record.seq.encode()` instead of `bytes(record.seq)`*.


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue tracker](https://github.com/althonos/pyaragorn/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


<!-- ### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyaragorn/blob/main/CONTRIBUTING.md)
for more details. -->


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyaragorn/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0 or later](https://choosealicense.com/licenses/gpl-3.0/).
ARAGORN and ARWEN were developed by Dean Laslett and are distributed under the
terms of the GPLv3 or later as well. See `vendor/aragorn` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the ARAGORN authors. It was developed by
[Martin Larralde](https://github.com/althonos/) during his PhD project
at the [Leiden University Medical Center](https://www.lumc.nl/en/) in
the [Zeller Lab](https://zellerlab.org).*


## 📚 References

- <a id="ref1">\[1\]</a> Laslett, Dean, and Bjorn Canback. “ARAGORN, a program to detect tRNA genes and tmRNA genes in nucleotide sequences.” Nucleic acids research vol. 32,1 11-6. 2 Jan. 2004, [doi:10.1093/nar/gkh152](https://doi.org/10.1093/nar/gkh152)
- <a id="ref2">\[2\]</a> Laslett, Dean, and Björn Canbäck. “ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences.” Bioinformatics (Oxford, England) vol. 24,2 (2008): 172-5. [doi:10.1093/bioinformatics/btm573](https://doi.org/10.1093/bioinformatics/btm573)
