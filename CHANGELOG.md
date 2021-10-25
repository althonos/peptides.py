# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/peptides.py/compare/v0.2.0...HEAD


## [v0.2.0] - 2021-10-21
[v0.2.0]: https://github.com/althonos/peptides.py/compare/v0.1.0...v0.2.0

### Added
- `Peptide.counts` method to get the number of occurences of each amino acid in the peptide.
- `Peptide.frequencies` to get the frequencies of each amino acid in the peptide.
- `Peptide.pcp_descriptors` to compute the PCP descriptors from Mathura & Braun (2001).
- `Peptide.sneath_vectors` to compute the descriptors from Sneath (1966).
- Hydrophilicity descriptors from Barley (2018).
- `Peptide.structural_class` to predict the structural class of a protein using one of three reference datasets and one of four distance metrics.

### Changed
- `Peptide.aliphatic_index` now supports unknown Leu/Ile residue (code *J*).
- Swap order of `Peptide.hydrophobic_moment` arguments for consistency with profile methods.
- Some `Peptide` functions now support vectorized code using `numpy` if available.


## [v0.1.0] - 2021-10-21
[v0.1.0]: https://github.com/althonos/peptides.py/compare/14f254e9...v0.1.0

Initial release.
