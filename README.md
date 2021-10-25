# `peptides.py` [![Stars](https://img.shields.io/github/stars/althonos/peptides.py.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/peptides.py/stargazers)

*Physicochemical properties, indices and descriptors for amino-acid sequences.*

[![Actions](https://img.shields.io/github/workflow/status/althonos/peptides.py/Test/main?logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/peptides.py/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/peptides.py?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/peptides.py/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/peptides.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/peptides)
[![Wheel](https://img.shields.io/pypi/wheel/peptides.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/peptides/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/peptides.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/peptides/#files)
[![Python Implementations](https://img.shields.io/badge/impl-universal-success.svg?style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/peptides/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/peptides.py/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/peptides.py/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/peptides.py.svg?style=flat-square&maxAge=600)](https://github.com/althonos/peptides.py/issues)
[![Docs](https://img.shields.io/readthedocs/peptides/latest?style=flat-square&maxAge=600)](https://peptides.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/peptides.py/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpeptides)](https://pepy.tech/project/peptides)

## 🗺️ Overview

`peptides.py` is a pure-Python package to compute common descriptors for
protein sequences. It started as a port of [`Peptides`](https://cran.r-project.org/web/packages/Peptides/index.html), the R package written by
[Daniel Osorio](https://orcid.org/0000-0003-4424-8422), but now also provides
some additional features from [EMBOSS](http://emboss.bioinformatics.nl/cgi-bin/emboss/),
[ExPASy Protein Identification and Analysis Tools](https://web.expasy.org/protparam/), and [Rcpi](https://bioconductor.org/packages/release/bioc/html/Rcpi.html).

This library has no external dependency and is available for all modern Python
versions (3.6+).

### 📋 Features

A non-exhaustive list of available features:

- Peptide statistics: amino acid counts and frequencies.
- **QSAR** descriptors: BLOSUM indices, Cruciani properties, FASGAI vectors, Kidera factors, MS-WHIM scores, PCP descriptors, ProtFP descriptors, Sneath vectors, ST-scales, T-scales, VHSE-scales, Z-scales.
- Sequence profiles: hydrophobicity, hydrophobic moment, membrane position.
- Physicochemical properties: aliphatic index, instability index, theoretical net charge, isoelectric point, molecular weight (with isotope labelling support).
- Biological properties: structural class prediction.

*If this library is missing a useful statistic or descriptor, feel free to
reach out and open a feature request on the [issue tracker](https://github.com/althonos/peptides.py/issues)
of the [project repository](https://github.com/althonos/peptides.py).*

### 🧊 Vectorization

Most of the descriptors for a protein sequence are simple averages of values
taken in a lookup table, so computing them can be done in a vectorized manner.
If [`numpy`](https://numpy.org/) can be imported, relevant functions
(like `numpy.sum` or `numpy.take`) will be used, otherwise a fallback
implementation using [`array.array`](https://docs.python.org/3/library/array.html#array.array)
from the standard library is available.

## 🔧 Installing

Install the `peptides` package directly from [PyPi](https://pypi.org/project/peptides)
which hosts universal wheels that can be installed with `pip`:
```console
$ pip install peptides
```

<!--
Otherwise, `peptides.py` is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda peptides-py
``` -->

## 📖 Documentation

A complete [API reference](https://peptides.readthedocs.io/en/stable/api.html)
can be found in the [online documentation](https://peptides.readthedocs.io/),
or directly from the command line using
[`pydoc`](https://docs.python.org/3/library/pydoc.html):
```console
$ pydoc peptides.Peptide
```

## 💡 Example

Start by creating a `Peptide` object from a protein sequence:
```python
>>> import peptides
>>> peptide = peptides.Peptide("MLKKRFLGALAVATLLTLSFGTPVMAQSGSAVFTNEGVTPFAISYPGGGT")
```

Then use the appropriate methods to compute the descriptors you want:
```python
>>> peptide.aliphatic_index()
89.8...
>>> peptide.boman()
-0.2097...
>>> peptide.charge(pH=7.4)
1.99199...
>>> peptide.isoelectric_point()
10.2436...
```

Methods that return more than one scalar value (for instance, `Peptide.blosum_indices`)
will return a dedicated named tuple:
```python
>>> peptide.ms_whim_scores()
MSWHIMScores(mswhim1=-0.436399..., mswhim2=0.4916..., mswhim3=-0.49200...)
```

Use the `Peptide.descriptors` method to get a dictionary with every available
descriptor. This makes it very easy to create a `pandas.DataFrame` with
descriptors for several protein sequences:
```python
>>> seqs = ["SDKEVDEVDAALSDLEITLE", "ARQQNLFINFCLILIFLLLI", "EGVNDNECEGFFSAR"]
>>> df = pandas.DataFrame([ peptides.Peptide(s).descriptors() for s in seqs ])
>>> df
    BLOSUM1   BLOSUM2  BLOSUM3   BLOSUM4  ...        Z2        Z3        Z4        Z5
0  0.367000 -0.436000   -0.239  0.014500  ... -0.711000 -0.104500 -1.486500  0.429500
1 -0.697500 -0.372500   -0.493  0.157000  ... -0.307500 -0.627500 -0.450500  0.362000
2  0.479333 -0.001333    0.138  0.228667  ... -0.299333  0.465333 -0.976667  0.023333

[3 rows x 66 columns]
```

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue
tracker](https://github.com/althonos/peptides.py/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.

### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/peptides.py/blob/main/CONTRIBUTING.md)
for more details.

## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
The original R `Peptides` package was written by [Daniel Osorio](https://orcid.org/0000-0003-4424-8422),
[Paola Rondón-Villarreal](https://orcid.org/0000-0001-8209-3885) and
[Rodrigo Torres](https://orcid.org/0000-0003-1113-3020), and is licensed under
the terms of the [GNU General Public License v2.0](https://choosealicense.com/licenses/gpl-2.0/).
The [EMBOSS](http://emboss.bioinformatics.nl/cgi-bin/emboss/) applications are
released under the [GNU General Public License v1.0](https://www.gnu.org/licenses/old-licenses/gpl-1.0.html).

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original `Peptides` authors](https://github.com/dosorio). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
