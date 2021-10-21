# `peptides.py` [![Stars](https://img.shields.io/github/stars/althonos/peptides.py.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/peptides.py/stargazers)

*Physicochemical properties and indices for amino-acid sequences.*

[![Actions](https://img.shields.io/github/workflow/status/althonos/peptides.py/Test/main?logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/peptides.py/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/peptides.py?style=flat-square&maxAge=3600)](https://codecov.io/gh/althonos/peptides.py/)
[![PyPI](https://img.shields.io/pypi/v/peptides.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/peptides)
[![Wheel](https://img.shields.io/pypi/wheel/peptides.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/peptides/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/peptides.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/peptides/#files)
[![Python Implementations](https://img.shields.io/badge/impl-universal-success.svg?style=flat-square&maxAge=3600&label=impl)](https://pypi.org/project/peptides/#files)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/peptides.py/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/peptides.py/)
[![GitHub issues](https://img.shields.io/github/issues/althonos/peptides.py.svg?style=flat-square&maxAge=600)](https://github.com/althonos/peptides.py/issues)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/peptides.py/blob/master/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpeptides)](https://pepy.tech/project/peptides)

## 🗺️ Overview

`peptides.py` is a pure-Python package to compute common descriptors for
protein sequences. It is a port of `Peptides`, the R package written by
[Daniel Osorio](https://orcid.org/0000-0003-4424-8422) for the same purpose.
This library has no external dependency and is available for all modern Python
versions (3.6+).

## 🔧 Installing

Install the `peptides` package directly from [PyPi](https://pypi.org/project/peptides)
which hosts universal wheels that can be installed with `pip`:
```console
$ pip install peptides
```

<!--
Otherwise, Pyrodigal is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda peptides-py
``` -->

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
the terms of the GPLv2.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [original `Peptides` authors](https://github.com/dosorio). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
