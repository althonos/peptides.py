.. peptides documentation master file, created by
   sphinx-quickstart on Sat Oct 23 18:17:35 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

`peptides.py` |stars|
=====================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/peptides.py.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/peptides.py/stargazers

*Physicochemical properties, indices and descriptors for amino-acid sequences.*

|Actions| |Coverage| |License| |PyPI| |Bioconda| |Wheel| |Versions| |Implementations| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/peptides.py/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/peptides.py/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/peptides.py?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/peptides.py/

.. |PyPI| image:: https://img.shields.io/pypi/v/peptides.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/peptides

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/peptides?style=flat-square&maxAge=3600&logo=anaconda
   :target: https://anaconda.org/bioconda/peptides

.. |Wheel| image:: https://img.shields.io/pypi/wheel/peptides?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/peptides/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/peptides.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/peptides/#files

.. |Implementations| image:: https://img.shields.io/badge/impl-universal-success.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/peptides/#files

.. |License| image:: https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/gpl-3.0/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/peptides.py/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/peptides.py/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/peptides.py.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/peptides.py/issues

.. |Docs| image:: https://img.shields.io/readthedocs/peptides?style=flat-square&maxAge=3600
   :target: http://peptides.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/peptides.py/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/pypi/dm/peptides?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/peptides


Overview
--------

``peptides.py`` is a pure-Python package to compute common descriptors for
protein sequences. It started as a port of `Peptides <https://cran.r-project.org/web/packages/Peptides/index.html>`_,
the R package written by `Daniel Osorio <https://orcid.org/0000-0003-4424-8422>`_
for the same purpose, but now also provides some more features from
`EMBOSS <http://emboss.bioinformatics.nl/cgi-bin/emboss/>`_,
`ExPASy Protein Identification and Analysis Tools <https://web.expasy.org/protparam/>`_,
and `Rcpi <https://bioconductor.org/packages/release/bioc/html/Rcpi.html>`_. This library has no external dependency and is
available for all modern Python versions (3.6+).

A non-exhaustive list of available features:

- Amino-acid Statistics:

  - Number of occurrences in the peptide sequence
  - Frequency in the peptide sequence

- `QSAR <https://en.wikipedia.org/wiki/Quantitative_structure%E2%80%93activity_relationship>`_ descriptors:

  - `BLOSUM indices <https://doi.org/10.1089/cmb.2008.0173>`_
  - `Cruciani properties <https://doi.org/10.1002/cem.856>`_
  - `FASGAI vectors <https://doi.org/10.1111/j.1747-0285.2008.00641.x>`_
  - `Kidera factors <https://doi.org/10.1007/BF01025492>`_
  - `MS-WHIM scores <https://doi.org/10.1021/ci980211b>`_
  - `PCP descriptors <https://doi.org/10.1007/s00894-001-0058-5>`_
  - `ProtFP descriptors <https://doi.org/10.1186/1758-2946-5-41>`_
  - `Sneath vectors <https://doi.org/10.1016/0022-5193(66)90112-3>`_
  - `ST-scales <https://doi.org/10.1007/s00726-009-0287-y>`_
  - `SVGER descriptors <https://doi.org/10.1002/minf.201501023>`_
  - `T-scales <https://doi.org/10.1016/j.molstruc.2006.07.004>`_
  - `VHSE-scales <https://doi.org/10.1002/bip.20296>`_
  - `Z-scales <https://doi.org/10.1021/jm9700575>`_

- Sequence profiles:

  - Hydrophobicity profile using one of 39 proposed scales.
  - Hydrophobic moment profile based on `Eisenberg, Weiss and Terwilliger (1984) <https://doi.org/10.1073/pnas.81.1.140>`_.
  - Membrane position based on `Eisenberg (1984) <https://doi.org/10.1146/annurev.bi.53.070184.003115>`_.

- Physical-chemical properties:

  - Aliphatic index proposed in `Ikai (1980) <https://pubmed.ncbi.nlm.nih.gov/7462208/>`_.
  - Instability index proposed in `Boman (2003) <https://doi.org/10.1046/j.1365-2796.2003.01228.x>`_.
  - Theoretical net charge based on the `Henderson-Hasselbach equation <https://en.wikipedia.org/wiki/Henderson%E2%80%93Hasselbalch_equation>`_.
  - Isoelectric point using one of 8 pKa scales.
  - Molecular weight, taking into account isotope labelling, using one of 3 average weight tables.

- Biological properties:

  - Structural class using methods and reference data from either
    `Nakashima, Nishikawa & Ooi (1985) <https://doi.org/10.1093/oxfordjournals.jbchem.a135454>`_,
    `Chou (1989) <https://doi.org/10.1007/978-1-4613-1571-1>`_,
    `Chou & Zhang (1992) <https://doi.org/10.1111/j.1432-1033.1992.tb17067.x>`_,
    or `Chou, Liu, Maggiora & Zhang (1998) <https://pubmed.ncbi.nlm.nih.gov/9552161/>`_.


Setup
-----

Run ``pip install peptides`` in a shell to download the latest release, or have
a look at the :doc:`Installation page <install>` to find other ways to install
``peptides.py``.


Library
-------

.. toctree::
  :maxdepth: 2

  Installation <install>
  Contributing <contributing>
  API Reference <api>
  Changelog <changes>


License
-------

This library is provided under the `GNU General Public License v3.0 <https://choosealicense.com/licenses/gpl-3.0/>`_.
The original R `Peptides` package was written by `Daniel Osorio <https://orcid.org/0000-0003-4424-8422>`_,
`Paola Rondón-Villarreal <https://orcid.org/0000-0001-8209-3885>`_ and
`Rodrigo Torres <https://orcid.org/0000-0003-1113-3020>`_, and is licensed under
the terms of the `GNU General Public License v2.0 <https://choosealicense.com/licenses/gpl-2.0/>`_.
The `EMBOSS <http://emboss.bioinformatics.nl/cgi-bin/emboss/>`_ applications are
released under the `GNU General Public License v1.0 <https://www.gnu.org/licenses/old-licenses/gpl-1.0.html>`_.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the original* `Peptides`_ *authors. It was developed by*
`Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
