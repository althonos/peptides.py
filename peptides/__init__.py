"""Physicochemical properties and indices for amino-acid sequences.
"""

import array
import collections
import functools
import math
import operator
import random
import statistics
import typing

from . import tables, datasets

try:
    import numpy
    _sum = numpy.sum
except ImportError:
    numpy = None
    _sum = sum

__all__ = ["Peptide", "tables", "datasets"]
__version__ = "0.4.0"
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPL-3.0-or-later"
__credits__ = """
Daniel Osorio, Paola Rondón-Villarreal and Rodrigo Torres for ``Peptides``.
Alan Bleasby for the ``hmoment`` binary of the EMBOSS.
Dong-Sheng Cao, Nan Xiao, Qing-Song Xu, and Alex F. Chen for ``Rcpi``.
""".strip()

# --- Constants --------------------------------------------------------------

# fmt: off
_CODE1 = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
    "O", "U", "B", "Z", "J", "X"
]

# fmt: off
_CODE3 = [
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
    "Pyl", "Sec", "Asx", "Glx", "Xle", "Xaa"
]

_TRANS = bytes(
    _CODE1.index(chr(x).upper()) if chr(x).upper() in _CODE1 else _CODE1.index("X")
    for x in range(256)
)


# --- Helper classes ---------------------------------------------------------

class BLOSUMIndices(typing.NamedTuple):
    """The BLOSUM62-derived indices of a peptide.

    BLOSUM indices were derived of physicochemical properties that have
    been subjected to a VARIMAX analysis and an alignment matrix of the
    20 natural AAs using the BLOSUM62 matrix.

    References:
        - Georgiev, A. G.
          *Interpretable Numerical Descriptors of Amino Acid Space*.
          Journal of Computational Biology. May 2009;16(5):703–23.
          :doi:`10.1089/cmb.2008.0173`. :pmid:`19432540`.

    """
    blosum1: float
    blosum2: float
    blosum3: float
    blosum4: float
    blosum5: float
    blosum6: float
    blosum7: float
    blosum8: float
    blosum9: float
    blosum10: float


class CrucianiProperties(typing.NamedTuple):
    """The Cruciani properties of a peptide.

    The Cruciani properties are a collection of scaled principal
    component scores that summarize a broad set of descriptors
    calculated based on the interaction of each amino acid residue with
    several chemical groups (or "probes"), such as charged ions, methyl,
    hydroxyl groups, and so forth.

    References:
        - Cruciani, G., M. Baroni, E. Carosati, M. Clementi, R. Valigi,
          and S. Clementi.
          *Peptide Studies by Means of Principal Properties of Amino
          Acids Derived from MIF Descriptors*. Journal of Chemometrics.
          2004;18(3-4):146–55. :doi:`10.1002/cem.856`.

    """
    pp1: float
    pp2: float
    pp3: float


class FasgaiVectors(typing.NamedTuple):
    """FASGAI vectors of a peptide.

    The FASGAI vectors (Factor Analysis Scales of Generalized Amino
    Acid Information) are a set of amino acid descriptors, that reflect
    hydrophobicity, alpha and turn propensities, bulky properties,
    compositional characteristics, local flexibility, and electronic
    properties, that can be utilized to represent the sequence
    structural features of peptides or protein motifs.

    References:
        - Liang, G., G. Chen, W. Niu, and Z. Li.
          *Factor Analysis Scales of Generalized Amino Acid Information
          as Applied in Predicting Interactions between the Human
          Amphiphysin-1 SH3 Domains and Their Peptide Ligands*.
          Chemical Biology & Drug Design. Apr 2008;71(4):345–51.
          :doi:`10.1111/j.1747-0285.2008.00641.x`. :pmid:`18318694`.

    """
    f1: float
    f2: float
    f3: float
    f4: float
    f5: float
    f6: float


class KideraFactors(typing.NamedTuple):
    """The Kidera factors of a peptide.

    The Kidera Factors were originally derived by applying multivariate
    analysis to 188 physical properties of the 20 amino acids and using
    dimension reduction techniques.

    Attributes:
        kf1 (`float`): A factor modeling the helix / bend preference.
        kf2 (`float`): A factor modeling the side-chain size of each
            residue (:aaindex:`KIDA850101`).
        kf3 (`float`): A factor modeling the extended structured preference.
        kf4 (`float`): A factor representing the hydrophobicity.
        kf5 (`float`): A factor modeling the double-bend preference.
        kf6 (`float`): A factor modeling the partial specific volume.
        kf7 (`float`): A factor modeling the flat extended preference.
        kf8 (`float`): A factor modeling the occurence in alpha regions.
        kf9 (`float`): A factor encoding the pK-C.
        kf10 (`float`): A factor representing the surrounding hydrophobicity.

    References:
        - Kidera, A., Y. Konishi, M. Oka, T. Ooi, and H. A. Scheraga.
          *Statistical Analysis of the Physical Properties of the 20
          Naturally Occurring Amino Acids*. Journal of Protein Chemistry.
          Feb 1985;4(1):23–55. :doi:`10.1007/BF01025492`.

    """
    kf1: float
    kf2: float
    kf3: float
    kf4: float
    kf5: float
    kf6: float
    kf7: float
    kf8: float
    kf9: float
    kf10: float


class AtchleyFactors(typing.NamedTuple):
    """The Atchley factors of a peptide.

    The Atchley Factors were originally derived by applying multivariate
    analysis to 494 physical properties of the 20 amino acids and using
    dimension reduction techniques.

    Attributes:
        af1 (`float`): A factor modeling polarity, accessibility, and hydrophobicity.
        af2 (`float`): A factor modeling propensity for secondary structure.
        af3 (`float`): A factor modeling molecular size.
        af4 (`float`): A factor modeling codon composition.
        af5 (`float`): A factor modeling electrostatic charge.

    References:
        - Atchley, W. R., Zhao, J., Fernandes, A. D., Drüke, T.
          *Solving the protein sequence metric problem*.
          Proceedings of the National Academy of Sciences.
          Apr 2005;102(18):6395-6400. :doi:`10.1073/pnas.040867710`.

    """
    af1: float
    af2: float
    af3: float
    af4: float
    af5: float


class MSWHIMScores(typing.NamedTuple):
    """The MS-WHIM scores of a peptide.

    MS-WHIM scores were derived from 36 electrostatic potential
    properties derived from the three-dimensional structure of the
    20 natural amino acids.

    References:
        - Bravi, G., E. Gancia, P. Mascagni, M. Pegna, R. Todeschini,
          and A. Zaliani.
          *MS-WHIM, New 3D Theoretical Descriptors Derived from
          Molecular Surface Properties: A Comparative 3D QSAR Study in a
          Series of Steroids*. Journal of Computer-Aided Molecular
          Design. Jan 1997;11(1):79-92.
          :doi:`10.1023/a:1008079512289`. :pmid:`9139115`
        - Gancia, E., G. Bravi, P. Mascagni, and A. Zaliani.
          *Global 3D-QSAR Methods: MS-WHIM and Autocorrelation*. Journal
          of Computer-Aided Molecular Design. Mar 2000;14(3):293–306.
          :doi:`10.1023/a:1008142124682`. :pmid:`10756483`.
        - Zaliani, A., and E. Gancia.
          *MS-WHIM Scores for Amino Acids: A New 3D-Description for
          Peptide QSAR and QSPR Studies*. Journal of Chemical
          Information and Computer Sciences. May 1999;39(3):525–33.
          :doi:`10.1021/ci980211b`.

    """
    mswhim1: float
    mswhim2: float
    mswhim3: float


class PhysicalDescriptors(typing.NamedTuple):
    """The Physical Descriptors of a peptide.

    The PP descriptors were constructed by improving on existing
    PCA-derived descriptors (Z-scales, MS-WHIM and T-scales) after
    correcting for the hydrophilicity of Methionine, Asparagine and
    Tryptophan based on Feng *et al*.

    Attributes:
        pd1 (`float`): A descriptor related to residue volume.
        pd2 (`float`): A descriptor related to hydrophilicity.

    Note:
        Barley *et al* insisted on maintaining a minimal number of
        descriptors as a way to reduce the chances of finding spurious
        QSAM models that would be affected by mutation between
        interaction sites.

    References:
        - Barley, M. H., N. J. Turner, and R. Goodacre.
          *Improved Descriptors for the Quantitative Structure–Activity
          Relationship Modeling of Peptides and Proteins*. Journal of
          Chemical Information and Modeling. Feb 2018;58(2):234–43.
          :doi:`10.1021/acs.jcim.7b00488`. :pmid:`29338232`.
        - Feng, X., J. Sanchis, M. T. Reetz, and H. Rabitz.
          *Enhancing the Efficiency of Directed Evolution in Focused
          Enzyme Libraries by the Adaptive Substituent Reordering
          Algorithm*. Chemistry. Apr 2012;18(18):5646–54.
          :doi:`10.1002/chem.201103811`. :pmid:`22434591`.

    """
    pd1: float
    pd2: float


class PCPDescriptors(typing.NamedTuple):
    """The Physical-Chemical Properties descriptors of a peptide.

    The PCP descriptors were constructed by performing multidimensional
    scaling of 237 physical-chemical properties.

    References:
        - Mathura, V. S., and W. Braun.
          *New Quantitative Descriptors of Amino Acids Based on
          Multidimensional Scaling of a Large Number of
          Physical–Chemical Properties*.
          Molecular Modeling Annual. Dec 2001;7(12):445–53.
          :doi:`10.1007/s00894-001-0058-5`.
        - Mathura, V. S., D. Paris, and M. J. Mullan.
          *A Novel Physico-Chemical Property Based Model for Studying
          the Effects of Mutation on the Aggregation of Peptides*.
          Protein and Peptide Letters. 2009;16(8):991–98.
          :doi:`10.2174/092986609788923220`. :pmid:`19689427`.

    """
    e1: float
    e2: float
    e3: float
    e4: float
    e5: float


class PRINComponents(typing.NamedTuple):
    """The PRIN components of a peptide.

    The PRIN components were constructed in Vicator *et al* (2005)
    by running PCA on different amino-acid properties.

    Attributes:
        prin1 (`float`): The PRIN1 component, representing hydrophobicity
            properties.
        prin2 (`float`): The PRIN2 component, representing the residue
            size.
        prin3 (`float`): The PRIN3 component, representing the pkN values
            of the amino-acids.

    References:
        - Vicatos, S., Reddy, B.V., and Y. Kaznessis.
          *Prediction of distant residue contacts with the use of
          evolutionary information*. Proteins. 2005 Mar 1;58(4):935-49.
          :doi:`10.1002/prot.20370` :pmid:`15645442`.

    """
    prin1: float
    prin2: float
    prin3: float


class ProtFPDescriptors(typing.NamedTuple):
    """The ProtFP descriptors of a peptide.

    The ProtFP set was constructed from a large initial selection of
    indices obtained from the `AAindex <https://www.genome.jp/aaindex/>`_
    database for all 20 naturally occurring amino acids.

    References:
        - van Westen, G. J., R. F. Swier, J. K. Wegner, A. P. Ijzerman,
          H. W. van Vlijmen, and A. Bender.
          *Benchmarking of Protein Descriptor Sets in Proteochemometric
          Modeling (Part 1): Comparative Study of 13 Amino Acid
          Descriptor Sets*. Journal of Cheminformatics. Sep 2013;5(1):41.
          :doi:`10.1186/1758-2946-5-41`. :pmid:`24059694`.
        - van Westen, G. J., R. F. Swier, I. Cortes-Ciriano,
          J. K. Wegner, J. P. Overington, A. P. Ijzerman,
          H. W. van Vlijmen, and A. Bender.
          *Benchmarking of Protein Descriptor Sets in Proteochemometric
          Modeling (Part 2): Modeling Performance of 13 Amino Acid
          Descriptor Sets*. Journal of Cheminformatics. Sep 2013;5(1):42.
          :doi:`10.1186/1758-2946-5-42`. :pmid:`24059743`.

    """
    protfp1: float
    protfp2: float
    protfp3: float
    protfp4: float
    protfp5: float
    protfp6: float
    protfp7: float
    protfp8: float


class SneathVectors(typing.NamedTuple):
    """The Sneath vectors of a peptide.

    These vectors were obtained in Sneath (1996) by running PCA on the
    `ϕ coefficient <https://en.wikipedia.org/wiki/Phi_coefficient>`_
    to explain the dissimilarity between the 20 natural amino acids
    based on binary state encoding of 134 physical and chemical
    properties (such as presence/absence of a —CH₃ group, step-wise
    optical rotation, etc.).

    Attributes:
        sv1 (`float`): A descriptor representing mainly aliphatic properties
            of each residue (:aaindex:`SNEP660101`).
        sv2 (`float`): A descriptor putatively modeling the number of
            reactive groups (:aaindex:`SNEP660102`).
        sv3 (`float`): A descriptor representing the aromatic properties
            of each residue (:aaindex:`SNEP660103`).
        sv4 (`float`): A descriptor with uncertain interpretation
            (:aaindex:`SNEP660104`).

    References:
        - Sneath, P. H. A.
          *Relations between Chemical Structure and Biological Activity
          in Peptides*.
          Journal of Theoretical Biology. Nov 1996;12(2):157–95.
          :doi:`10.1016/0022-5193(66)90112-3`. :pmid:`4291386`.

    """
    sv1: float
    sv2: float
    sv3: float
    sv4: float


class STScales(typing.NamedTuple):
    """The ST-scales of a peptide.

    The ST-scales were proposed in Yang *et al* (2010), taking 827
    properties into account which are mainly constitutional,
    topological, geometrical, hydrophobic, electronic, and steric
    properties of a total set of 167 amino acids.

    References:
        - Yang, L., M. Shu, K. Ma, H. Mei, Y. Jiang, and Z. Li.
          *ST-Scale as a Novel Amino Acid Descriptor and Its Application
          in QSAM of Peptides and Analogues*.
          Amino Acids. Mar 2010;38(3):805–16.
          :doi:`10.1007/s00726-009-0287-y`. :pmid:`19373543`.

    """
    st1: float
    st2: float
    st3: float
    st4: float
    st5: float
    st6: float
    st7: float
    st8: float


class SVGERDescriptors(typing.NamedTuple):
    """The SVGER descriptors of a peptide.

    SVGER descriptors were constructed by Principal Component Analysis
    of 74 geometrical descriptors (`svger1` to `svger6`), 44 eigenvalue
    descriptors (`svger7`, `svger8` and `svger9`), and 41 Randić
    descriptors (`svger10` and `svger11`) computed for the 20 proteinogenic
    amino acids.

    References:
        - Tong, J., L. Li, M. Bai, and K. Li.
          *A New Descriptor of Amino Acids-SVGER and Its Applications in
          Peptide QSAR*. Molecular Informatics 36, no. 5–6 (2017): 1501023.
          :doi:`10.1002/minf.201501023`.
        - Randic, M.
          *Molecular Shape Profiles*. Journal of Chemical Information and
          Computer Sciences 35, no. 3 (1 May 1995): 373–82.
          :doi:`10.1021/ci00025a005`.

    """
    svger1: float
    svger2: float
    svger3: float
    svger4: float
    svger5: float
    svger6: float
    svger7: float
    svger8: float
    svger9: float
    svger10: float
    svger11: float


class TScales(typing.NamedTuple):
    """The T-scales of a peptide.

    The T-scales are based on 67 common topological descriptors of 135
    amino acids. These topological descriptors are based on the
    connectivity table of amino acids alone, and to not explicitly
    consider 3D properties of each structure.

    References:
        - Tian, F., P. Zhou, and Z. Li.
          *T-Scale as a Novel Vector of Topological Descriptors for
          Amino Acids and Its Application in QSARs of Peptides*.
          Journal of Molecular Structure. Mar 2007;830(1):106–15.
          :doi:`10.1016/j.molstruc.2006.07.004`.

    """
    t1: float
    t2: float
    t3: float
    t4: float
    t5: float


class VHSEScales(typing.NamedTuple):
    """The VHSE-scales of a peptide.

    The VHSE-scales (principal components score Vectors of Hydrophobic,
    Steric, and Electronic properties), are derived from principal
    components analysis (PCA) on independent families of 18 hydrophobic
    properties, 17 steric properties, and 15 electronic properties,
    respectively, which are included in total 50 physicochemical
    variables of 20 coded amino acids.

    Attributes:
        vhse1 (`float`): A descriptor representing hydrophobic
            properties.
        vhse2 (`float`): Another descriptor representing hydrophobic
            properties.
        vhse3 (`float`): A descriptor representing steric
            properties.
        vhse4 (`float`): Another descriptor representing steric
            properties.
        vhse5 (`float`): A descriptor representing electronic
            properties.
        vhse6 (`float`): A second descriptor representing electronic
            properties.
        vhse7 (`float`): A third descriptor representing electronic
            properties.
        vhse8 (`float`): A fourth descriptor representing electronic
            properties.

    References:
        - Mei, H., Z. H. Liao, Y. Zhou, and S. Z. Li. *A New Set of
          Amino Acid Descriptors and Its Application in Peptide QSARs*.
          Biopolymers. 2005;80(6):775-86.
          :doi:`10.1002/bip.20296`. :pmid:`15895431`.

    """
    vhse1: float
    vhse2: float
    vhse3: float
    vhse4: float
    vhse5: float
    vhse6: float
    vhse7: float
    vhse8: float


class VSTPVDescriptors(typing.NamedTuple):
    """The VSTPV descriptors of a peptide.

    The VSTPV (vector of structural and topological variables) are derived
    from principal component analysis (PCA) of 85 structural variables of
    166 amino acids.

    References:
        - Shu, M., Yu, R., Zhang, Y., Wang, J., Wang, L., Lin, Z.
          *Predicting the Activity of Antimicrobial Peptides with Amino Acid
          Topological Information*.
          Medicinal Chemistry. Feb 2013;9(1):32-44.
          :doi:`10.2174/157340613804488350`.

    """
    vstpv1: float
    vstpv2: float
    vstpv3: float
    vstpv4: float
    vstpv5: float
    vstpv6: float


class ZScales(typing.NamedTuple):
    """The Z-scales of a peptide.

    The Z-scales were proposed in Sandberg *et al* (1998) based on
    physicochemical properties of proteogenic and non-proteogenic
    amino acids, including NMR data and thin-layer chromatography
    (TLC) data.

    Attributes:
        z1 (`float`): A descriptor quantifying lipophilicity.
        z2 (`float`): A descriptor modeling steric properties like steric
            bulk and polarizability.
        z3 (`float`): A descriptor quantifying electronic properties like
            polarity and charge.
        z4 (`float`): A descriptor relating to electronegativity, heat of
            formation, electrophilicity and hardness.
        z5 (`float`): Another descriptor relating to electronegativity,
            heat of formation, electrophilicity and hardness.

    References:
        - Sandberg, M., L. Eriksson, J. Jonsson, M. Sjöström, and
          S. Wold. *New Chemical Descriptors Relevant for the Design of
          Biologically Active Peptides. A Multivariate Characterization
          of 87 Amino Acids*.
          Journal of Medicinal Chemistry. Jul 1998;41(14):2481–91.
          :doi:`10.1021/jm9700575`. :pmid:`9651153`.

    """
    z1: float
    z2: float
    z3: float
    z4: float
    z5: float


# --- Peptide class ----------------------------------------------------------

class Peptide(typing.Sequence[str]):
    """A sequence of amino acids.

    Attributes:
        sequence (`str`): The peptide primary sequence, encoded in the IUPAC
            one-letter code.

    """

    # --- Class methods ------------------------------------------------------

    @classmethod
    def sample(
        cls,
        length: int,
        frequencies: str = "SwissProt2021",
    ) -> "Peptide":
        """Generate a peptide with the given amino-acid frequencies.

        This method is useful for testing, but using amino-acid frequencies
        to generate a peptide is not a biologically accurate method, instead
        consider sampling based on dipeptide frequencies in a particular
        organism, or using k-mer shuffling.

        Arguments:
            length (`int`): The desired length for the generated peptide.
            frequencies (`str`): The name of the amino-acid frequency table
                to use: either *KingJukes* to use the amino-acid frequencies
                for vertebrate organisms reported in King & Jukes (1969),
                or *SwissProt2021* to use the amino-acid frequencies in all
                the proteins from the January 2021 release of SwissProt.

        Returns:
            `~peptides.Peptide`: A new peptide. The first amino-acid will
            always be a Methionine for biological accuracy.

        References:
            - King, J. L., and T. H. Jukes.
              *Non-Darwinian Evolution*.
              Science. May 1969;164(3881):788–98.
              :doi:`10.1126/science.164.3881.788`. :pmid:`5767777`.
            - The UniProt Consortium.
              *UniProt: The Universal Protein Knowledgebase in 2021*.
              Nucleic Acids Research. Jan 2021;49(D1):D480–89.
              :doi:`10.1093/nar/gkaa1100`. :pmid:`33237286`.

        """
        table = tables.AA_FREQUENCIES.get(frequencies)
        if table is None:
            raise ValueError(f"Invalid amino acid frequencies: {frequencies!r}")

        if length == 0:
            return cls("")

        cumfreq = 0
        cumulative_frequencies = {}
        for k,v in table.items():
            cumfreq += v
            cumulative_frequencies[k] = cumfreq

        residues = ["M"]
        for i in range(1, length):
            x = random.random()
            r = next((k for k,v in cumulative_frequencies.items() if x <= v), "X")
            residues.append(r)

        return cls("".join(residues))

    # --- Magic methods ------------------------------------------------------

    def __init__(self, sequence: str) -> None:
        """Create a new peptide object with the given sequence.

        Arguments:
            sequence (`str`): A sequence of amino acids encoded with the
                IUPAC one-letter code. Non-standard (*O*, *U*), ambiguous
                (*B*, *Z*, *J*) and unknown (*X*) residues are supported
                in some methods, but not all of them.

        """
        # store the sequence in text format
        self.sequence: str = sequence
        # encode the sequence as ASCII bytes
        binary = sequence.encode('ascii')
        # encode the peptide with the given encoding
        self.encoded = array.array('B', binary.translate(_TRANS))

    def __len__(self) -> int:
        return len(self.sequence)

    @typing.overload
    def __getitem__(self, index: slice) -> "Peptide":
        pass

    @typing.overload
    def __getitem__(self, index: int) -> str:
        pass

    def __getitem__(
        self, index: typing.Union[int, slice]
    ) -> typing.Union[str, "Peptide"]:
        if isinstance(index, slice):
            return Peptide(self.sequence[index])
        return self.sequence[index]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.sequence!r})"

    # --- Generic methods ----------------------------------------------------

    def descriptors(self) -> typing.Dict[str, float]:
        """Create a dictionary containing every protein descriptor available.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLE")
            >>> sorted(peptide.descriptors().keys())
            ['AF1', ..., 'F1', ..., 'KF1', ..., 'MSWHIM1', ..., 'PP1', ...]

        Hint:
            Use this method to create a `~pandas.DataFrame` containing the
            descriptors for several sequences.

        """
        d = {}
        for prefix, method in self.__DESCRIPTORS.items():
            for i, x in enumerate(method(self)):
                d[f"{prefix}{i+1}"] = x
        return d

    def auto_correlation(
        self, table: typing.Dict[str, float], lag: int = 1, center: bool = True
    ) -> float:
        """Compute the auto-correlation index of a peptide sequence.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLE")
            >>> table = peptides.tables.HYDROPHOBICITY["KyteDoolittle"]
            >>> peptide.auto_correlation(table=table)
            -0.3519908...
            >>> peptide.auto_correlation(table=table, lag=5)
            0.00113355...

        """
        # center the table if requested
        if center:
            mu = statistics.mean(table.values())
            sigma = statistics.stdev(table.values())
            table = {k: (v - mu) / sigma for k, v in table.items()}
        # build look up table
        lut = [table.get(aa, 0.0) for aa in _CODE1]
        # compute using Cruciani formula
        if numpy is None:
            s1 = s2 = 0.0
            for aa1, aa2 in zip(self.encoded[:-lag], self.encoded[lag:]):
                s1 += lut[aa1] * lut[aa2]
                s2 += lut[aa1] ** 2
        else:
            v1 = numpy.take(lut, self.encoded[:-lag])
            v2 = numpy.take(lut, self.encoded[lag:])
            s1 = numpy.sum(v1*v2)
            s2 = numpy.sum(v1**2)
        return s1 / s2

    def auto_covariance(
        self, table: typing.Dict[str, float], lag: int = 1, center: bool = True
    ) -> float:
        """Compute the auto-covariance index of a peptide sequence.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLE")
            >>> table = peptides.tables.HYDROPHOBICITY["KyteDoolittle"]
            >>> peptide.auto_covariance(table)
            -0.414005...
            >>> peptide.auto_covariance(table, lag=5)
            0.0010003...

        """
        # center the table if requested
        if center:
            mu = statistics.mean(table.values())
            sigma = statistics.stdev(table.values())
            table = {k: (v - mu) / sigma for k, v in table.items()}
        # build the lookup table
        lut = [table.get(aa, 0.0) for aa in _CODE1]
        # compute correlation using Cruciani formula
        if numpy is None:
            s = 0.0
            for aa1, aa2 in zip(self.encoded[:-lag], self.encoded[lag:]):
                s += lut[aa1] * lut[aa2]
        else:
            v1 = numpy.take(lut, self.encoded[:-lag])
            v2 = numpy.take(lut, self.encoded[lag:])
            s = numpy.sum(v1*v2)
        return s / len(self)

    def cross_covariance(
        self,
        table1: typing.Dict[str, float],
        table2: typing.Dict[str, float],
        lag: int = 1,
        center: bool = True,
    ) -> float:
        """Compute the cross-covariance index of a peptide sequence.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLE")
            >>> table1 = peptides.tables.HYDROPHOBICITY["KyteDoolittle"]
            >>> table2 = peptides.tables.HYDROPHOBICITY["Eisenberg"]
            >>> peptide.cross_covariance(table1, table2)
            -0.3026609...
            >>> peptide.cross_covariance(table1, table2, lag=5)
            0.0259803...

        """
        # center the tables if requested
        if center:
            mu1 = statistics.mean(table1.values())
            sigma1 = statistics.stdev(table1.values())
            table1 = {k: (v - mu1) / sigma1 for k, v in table1.items()}
            mu2 = statistics.mean(table2.values())
            sigma2 = statistics.stdev(table2.values())
            table2 = {k: (v - mu2) / sigma2 for k, v in table2.items()}

        # build the lookup table
        lut1 = [table1.get(aa, 0.0) for aa in _CODE1]
        lut2 = [table2.get(aa, 0.0) for aa in _CODE1]

        # compute using Cruciani formula
        if numpy is None:
            s = 0.0
            for aa1, aa2 in zip(self.encoded[:-lag], self.encoded[lag:]):
                s += lut1[aa1] * lut2[aa2]
        else:
            v1 = numpy.take(lut1, self.encoded[:-lag])
            v2 = numpy.take(lut2, self.encoded[lag:])
            s = numpy.sum(v1*v2)
        return s / len(self)

    def profile(
        self,
        table: typing.Dict[str, float],
        window: int = 1,
        default: float = 0.0,
    ) -> typing.Sequence[float]:
        """Compute a generic per-residue profile from per-residue indices.

        Arguments:
            table (`dict`): The values per residue to apply to the whole
                protein sequence.
            window (`int`): The window size for computing the profile.
                Leave as *1* to return per-residue values.
            default (`float`): The default value to use for amino-acids
                that are not present in the given table.

        Returns:
            `collections.abc.Sequence` of `float`: The per-residue profile
            values, averaged in the given window size. When ``window`` is
            larger than the available number of resiudes, an empty sequence
            is returned.

        Example:
            >>> peptide = Peptide("PKLVCLKKC")
            >>> peptide.profile(peptides.tables.CHARGE['sign'])
            [0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 1.0, -1.0]
            >>> peptide.profile(peptides.tables.MOLECULAR_WEIGHT['expasy'], 5)
            [108..., 111..., 111..., 114..., 115...]

        .. versionadded:: 0.3.0

        """
        if window < 1:
            raise ValueError("Window must be strictly positive")

        # skip computing profile is window is larger than the available
        # number of residues in the peptide sequence
        if len(self) >= window:
            # build a look-up table and index values
            lut = [table.get(aa, default) for aa in _CODE1]
            if numpy is None:
                values = [lut[i] for i in self.encoded]
            else:
                values = numpy.take(lut, self.encoded)  # type: ignore
            # don't perform window averaging if window is 1
            if window <= 1:
                return list(values)
            elif window > 1:
                p = []
                # use a rolling sum over the window
                s = 0.0
                for i in range(window):
                    s += values[i]
                for j in range(window, len(self)):
                    p.append(s / window)
                    s -= values[j-window]
                    s += values[j]
                p.append(s / window)
        else:
            p = []

        return p

    # --- Sequence properties -----------------------------------------------

    def counts(self) -> typing.Dict[str, int]:
        """Return a table of amino-acid counts in the peptide.

        Returns:
            `dict`: A dictionary mapping each amino-acid code to the number
            of times it occurs in the peptide sequence.

        Example:
            >>> p = Peptide("SDKEVDEVDAALS")
            >>> {k:v for k,v in p.counts().items() if v != 0}
            {'A': 2, 'D': 3, 'E': 2, 'L': 1, 'K': 1, 'S': 2, 'V': 2}

        """
        # NB: This is consistently faster than using a `collections.Counter`
        #     because `str.count` has a fast C implementation based on
        #     `memchr` while `collections.Counter` has to iterate, which has
        #     some overhead.
        return {
            aa: self.sequence.count(aa)
            for aa in _CODE1
        }

    def frequencies(self) -> typing.Dict[str, float]:
        """Return a table of amino-acid frequencies in the peptide.

        Returns:
            `dict`: A dictionary mapping each amino-acid code to its
            frequency in the peptide sequence.

        Example:
            >>> p = Peptide("AALS")
            >>> {k:v for k,v in p.frequencies().items() if v != 0}
            {'A': 0.5, 'L': 0.25, 'S': 0.25}

        """
        return {
            aa:count/len(self)
            for aa,count in self.counts().items()
        }

    # --- Physico-chemical properties ----------------------------------------

    def aliphatic_index(self) -> float:
        r"""Compute the aliphatic index of the peptide.

        The aliphatic index of a protein was proposed in Ikai (1980). It is
        defined as the relative volume occupied by aliphatic side chains
        (Alanine, Valine, Isoleucine, and Leucine):

        .. math::

            \text{aliphatic index} = A + 2.9 V + 3.9 (I + L + J)


        It may be regarded as a positive factor for the increase of
        thermostability of globular proteins.

        Returns:
            `float`: The computed aliphatic index for the peptide
            sequence, between *0.0* and *390.0*.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLE")
            >>> peptide.aliphatic_index()
            117.0

        References:
            - Ikai, A.
              *Thermostability and Aliphatic Index of Globular Proteins*.
              Journal of Biochemistry. Dec 1980;88(6):1895–98.
              :pmid:`7462208`.

        """
        # NOTE(@althonos): Counting 5 times is faster than using a loop-based
        #                  counter or a `collections.Counter` as it will use
        #                  the C-level `memchr` if available, which is much
        #                  faster than iterating in Python space.
        # count aliphatic residues
        ala = self.sequence.count("A") / len(self.sequence)
        val = self.sequence.count("V") / len(self.sequence)
        leu = self.sequence.count("L") / len(self.sequence)
        ile = self.sequence.count("I") / len(self.sequence)
        # support unknown Leu/Ile residues
        xle = self.sequence.count("J") / len(self.sequence)
        # return aliphatic index
        return (ala + 2.9 * val + 3.9 * (leu + ile + xle)) * 100

    def boman(self) -> float:
        """Compute the Boman (potential peptide interaction) index.

        The potential interaction index proposed by Boman (2003) is an index
        computed by averaging the solubility values for all residues in a
        sequence. It can be used to give an overall estimate of the
        potential of a peptide to bind to membranes or other proteins.

        Returns:
            `float`: The Boman index for the peptide. A value greater than
            *2.48* indicates that a protein has high binding potential.

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.boman()
            -1.2358...

        Note:
            The potential protein interaction index was originally proposed
            as an easy way to differentiate between the action mechanism of
            hormones (protein/protein) and antimicrobial peptides
            (protein/membrane).

        References:
            - Boman, H. G.
              *Antibacterial Peptides: Basic Facts and Emerging Concepts*.
              Journal of Internal Medicine. 2003 Sep;254(3):197–215.
              :doi:`10.1046/j.1365-2796.2003.01228.x`. :pmid:`12930229`.

        """
        return -_sum(self.profile(tables.BOMAN["Boman"])) / len(self)

    def charge(self, pH: float = 7, pKscale: str = "Lehninger") -> float:
        """Compute the theoretical net charge of a peptide sequence.

        This function computes the theoretical net charge of a peptide
        sequence, based on the Henderson-Hasselbach equation described by
        Dexter S. Moore (1985). The net charge can be computed at a given pH
        using one of the 9 pKa scales available.

        Arguments:
            pH (`float`): The pH value for which to compute the charge.
            pKscale (`str`): The name of the pKa scale to be used. A list of
                all the allowed values can be retrieved from the keys of the
                `peptides.tables.PK` dictionary.

        Returns:
            `float`: The net charge of the peptide.

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.charge(pKscale="Bjellqvist")
            2.7373...
            >>> peptide.charge(pKscale="EMBOSS")
            2.9141...
            >>> peptide.charge(pKscale="Murray")
            2.9075...
            >>> peptide.charge(pKscale="Sillero")
            2.9198...
            >>> peptide.charge(pKscale="Solomon")
            2.8444...
            >>> peptide.charge(pKscale="Stryer")
            2.8765...
            >>> peptide.charge(pKscale="Lehninger")
            2.8731...
            >>> peptide.charge(pKscale="Dawson")
            2.8444...
            >>> peptide.charge(pKscale="Rodwell")
            2.8197...

        References:
            - Bjellqvist, B., G. J. Hughes, C. Pasquali, N. Paquet,
              F. Ravier, J. C. Sanchez, S. Frutiger, and D. Hochstrasser.
              *The Focusing Positions of Polypeptides in Immobilized pH
              Gradients Can Be Predicted from Their Amino Acid Sequences.*
              Electrophoresis. 1993 Oct;14(10):1023–31.
              :doi:`10.1002/elps.11501401163`. :pmid:`8125050`.
            - Dawson, R. M. C. and D. C. Elliott.
              *Data for Biochemical Research.*
              Oxford: Clarendon Press. 2002;3:592.
              ISBN:978-0-19-855299-4.
            - Kiraga, J.
              *Analysis and computer simulations of variability of
              isoelectric  point of proteins in the proteomes.*
              PhD thesis, University of Wroclaw, Poland. 2008.
            - Lehninger, A. L., D. L. Nelson, and M. M. Cox.
              *Lehninger Principles of Biochemistry*. 4th ed.
              New York: W.H. Freeman. 2005;4:1100.
              ISBN:978-0-7167-4339-2.
            - Murray, R. K.
              *Harper’s Illustrated Biochemistry.*
              New York: Lange Medical Books/McGraw-Hill. 2006;27.
              ISBN:978-0-07-146197-9.
            - Rodwell, J.D.
              *Heterogeneity of Component Bands in Isoelectric Focusing
              Patterns*. Analytical Biochemistry. 1982 Jan;119(2):440-49.
              :doi:`10.1016/0003-2697(82)90611-x`. :pmid:`7072964`.
            - Sillero, A., and A. Maldonado.
              *Isoelectric Point Determination of Proteins and Other
              Macromolecules: Oscillating Method*.
              Computers in Biology and Medicine. 2006 Feb;36(2): 157–66.
              :doi:`10.1016/j.compbiomed.2004.09.006`. :pmid:`16389075`.
            - Solomons, T. W. G.
              *Fundamentals of Organic Chemistry*. New York: Wiley. 1997;5.
              ISBN:978-0-471-28298-3.
            - Stryer, L., J. Augustyniak, and J. Michejda.
              *Biochemia*. Warszawa: Wydawnictwo Naukowe PWN. 2000.
              ISBN:978-83-01-12044-3.

        """
        # get chosen the pKa scale
        scale_pKa = tables.PK.get(pKscale)
        if scale_pKa is None:
            raise ValueError(f"Invalid pK scale: {scale!r}")
        scale_sign = tables.CHARGE["sign"]

        # build a look-up table for the pKa scale and the charge sign
        lut_pKa = [scale_pKa.get(aa, 0.0) for aa in _CODE1]
        lut_sign = [scale_sign.get(aa, 0.0) for aa in _CODE1]

        # compute charge of each amino-acid, and sum
        if numpy is None:
            charge = 0.0
            for aa in self.encoded:
                pKa = lut_pKa[aa]
                sign = lut_sign[aa]
                charge += sign / (1.0 + 10 ** (sign * (pH - pKa)))
        else:
            pKa = numpy.take(lut_pKa, self.encoded)
            sign = numpy.take(lut_sign, self.encoded)
            charge = numpy.sum(sign / (1.0 + 10**(sign * (pH - pKa))))

        # add charge for C-terminal and N-terminal ends of the peptide
        if "nTer" in scale_pKa:
            charge += 1.0 / (1.0 + 10 ** (pH - scale_pKa["nTer"]))
        if "cTer" in scale_pKa:
            charge += -1.0 / (1.0 + 10 ** (scale_pKa["cTer"] - pH))

        # return the net protein charge
        return charge

    def hydrophobic_moment(self, window: int = 11, angle: int = 100) -> float:
        """Compute the maximal hydrophobic moment of a protein sequence.

        This function computes the hydrophobic moment based on Eisenberg
        *et al* (1984). Hydrophobic moment is a quantitative measure of
        the amphiphilicity perpendicular to the axis of any periodic peptide
        structure, such as the α-helix or β-sheet.

        Arguments:
            angle (`int`): A protein rotational angle, in **degrees**.
                Usual values are *100* for α-helix, and *160* for β-sheet.
            window (`int`): The size of the sliding window for which to
                compute the local hydrophobic moment.

        Returns:
            `float`: The maximal hydrophobic moment of the peptide.

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.hydrophobic_moment(angle=100)
            0.519922...
            >>> peptide.hydrophobic_moment(angle=160)
            0.270590...

        See Also:
            The `~Peptide.hydrophobic_moment_profile` method, which builds
            a profile for each amino acid position instead of simply
            extracting the global maximum.

        References:
            - Eisenberg, D., R. M. Weiss, and T. C. Terwilliger.
              *The Hydrophobic Moment Detects Periodicity in Protein
              Hydrophobicity*. Proceedings of the National Academy of
              Sciences of the United States of America. 1984 Jan;81(1):140–44.
              :doi:`10.1073/pnas.81.1.140`. :pmid:`6582470`.

        """
        window = min(window, len(self))
        scale = tables.HYDROPHOBICITY["Eisenberg"]
        lut = [scale.get(aa, 0.0) for aa in _CODE1]
        angles = [(angle * i) % 360 for i in range(window)]

        if numpy is None:
            angsin = [math.sin(math.radians(theta)) for theta in angles]
            angcos = [math.cos(math.radians(theta)) for theta in angles]
        else:
            angsin = numpy.sin(numpy.radians(angles))
            angcos = numpy.cos(numpy.radians(angles))

        maxnorm = 0.0
        for i in range(len(self.sequence) - window + 1):
            # compute sin and cos of angles
            if numpy is None:
                sumsin = sumcos = 0
                for aa, s, c in zip(self.encoded[i:i+window], angsin, angcos):
                    sumsin += lut[aa]*s
                    sumcos += lut[aa]*c
            else:
                hvec = numpy.take(lut, self.encoded[i:i+window])
                sumsin = numpy.sum(hvec * angsin)
                sumcos = numpy.sum(hvec * angcos)
            # compute only the distance component (this way we can avoid
            # computing the square root in each iteration)
            norm = sumsin**2 + sumcos**2
            if norm > maxnorm:
                maxnorm = norm

        # compute the angular moment from the norm
        return math.sqrt(maxnorm) / window

    def hydrophobicity(self, scale: str = "KyteDoolittle") -> float:
        """Compute the hydrophobicity index of a protein sequence.

        This function calculates the hydrophobicity index of an amino
        acid sequence by averaging the hydrophobicity values of each residue
        using one of the 39 scales from different sources.

        Arguments:
            scale (`str`): The name of the hydrophobicity scale to be used.
                A list of all the allowed values can be retrieved from the
                keys of the `peptides.tables.HYDROPHOBICITY` dictionary.

        Returns:
            `float`: The hydrophobicity index of the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> peptide.hydrophobicity(scale="Aboderin")
            3.84...
            >>> peptide.hydrophobicity(scale="AbrahamLeo")
            0.092...

        Note:
             The hydrophobicity is an important stabilization force in
             protein folding; this force changes depending on the solvent
             in which the protein is found.

        References:
            - Aboderin, A. A.
              *An Empirical Hydrophobicity Scale for α-Amino-Acids and Some
              of Its Applications*.
              International Journal of Biochemistry. 1971 Oct;2(11):537–44.
              :doi:`10.1016/0020-711X(71)90023-1`.
            - Abraham, D.J., and A. J. Leo.
              *Extension of the Fragment Method to Calculate Amino Acid
              Zwitterion and Side Chain Partition Coefficients*.
              Proteins: Structure, Function, and Genetics. 1987;2(2):130–52.
              :doi:`10.1002/prot.340020207`.
            - Argos, P., J. K. Rao, and P. A. Hargrave.
              *Structural Prediction of Membrane-Bound Proteins*.
              European Journal of Biochemistry. Nov 1982;128(2–3):565–75.
              :doi:`10.1111/j.1432-1033.1982.tb07002.x`. :pmid:`7151796`.
            - Barley, M. H., N. J. Turner, and R. Goodacre.
              *Improved Descriptors for the Quantitative Structure–Activity
              Relationship Modeling of Peptides and Proteins*. Journal of
              Chemical Information and Modeling. Feb 2018;58(2):234–43.
              :doi:`10.1021/acs.jcim.7b00488`. :pmid:`29338232`.
            - Black, S. D., and D. R. Mould.
              *Development of Hydrophobicity Parameters to Analyze Proteins
              Which Bear Post- or Cotranslational Modifications*.
              Analytical Biochemistry. Feb 1991;193(1):72–82.
              :doi:`10.1016/0003-2697(91)90045-u`. :pmid:`2042744`.
            - Bull, H. B., and K. Breese.
              *Surface Tension of Amino Acid Solutions: A Hydrophobicity
              Scale of the Amino Acid Residues*. Archives of Biochemistry
              and Biophysics. Apr 1974;161(2):665–70.
              :doi:`10.1016/0003-9861(74)90352-x`. :pmid:`4839053`.
            - Casari, G., and M. J. Sippl.
              *Structure-Derived Hydrophobic Potential. Hydrophobic
              Potential Derived from X-Ray Structures of Globular Proteins
              Is Able to Identify Native Folds*.
              Journal of Molecular Biology. Apr 1992;224(3):725–32.
              :doi:`10.1016/0022-2836(92)90556-y`. :pmid:`1569551`.
            - Chothia, C.
              *The Nature of the Accessible and Buried Surfaces in Proteins*.
              Journal of Molecular Biology. Jul 2917;105(1):1–12.
              :doi:`10.1016/0022-2836(76)90191-1`. :pmid:`994183`.
            - Cid, H., M. Bunster, M. Canales, and F. Gazitúa.
              *Hydrophobicity and Structural Classes in Proteins*.
              Protein Engineering. Jul 1992;5(5):373–75.
              :doi:`10.1093/protein/5.5.373`. :pmid:`1518784`.
            - Cowan, R., and R. G. Whittaker.
              *Hydrophobicity Indices for Amino Acid Residues as Determined
              by High-Performance Liquid Chromatography*.
              Peptide Research. Apr 1990;3(2):75–80.
              :pmid:`2134053`.
            - Eisenberg, D., E. Schwarz, M. Komaromy, and R. Wall.
              *Analysis of Membrane and Surface Protein Sequences with
              the Hydrophobic Moment Plot*.
              Journal of Molecular Biology. Oct 1984;179(1):125–42.
              :doi:`10.1016/0022-2836(84)90309-7`. :pmid:`6502707`.
            - Engelman, D. M., T. A. Steitz, and A. Goldman.
              *Identifying Nonpolar Transbilayer Helices in Amino Acid
              Sequences of Membrane Proteins*. Annual Review of Biophysics
              and Biophysical Chemistry. 1986;15:321–53.
              :doi:`10.1146/annurev.bb.15.060186.001541`. :pmid:`3521657`.
            - Fasman, G. D.
              *Prediction of Protein Structure and the Principles of Protein
              Conformation*. Springer US. 1989.
              :doi:`10.1007/978-1-4613-1571-1`. ISBN:978-0-306-43131-9.
            - Fauchère, J-L., and Pliska V.
              *Hydrophobic Parameters π of Amino-Acid Side Chains from the
              Partitioning of N-Acetyl-Amino-Acid Amides*. European Journal
              of Medicinal Chemistry. 1983;18(4):369–75.
            - Goldsack, D. E., and R. C. Chalifoux.
              *Contribution of the Free Energy of Mixing of Hydrophobic Side
              Chains to the Stability of the Tertiary Structure of Proteins*.
              Journal of Theoretical Biology. Jun 1973;39(3):645–51.
              :doi:`10.1016/0022-5193(73)90075-1`. :pmid:`4354159`.
            - Guy, H. R.
              *Amino Acid Side-Chain Partition Energies and Distribution
              of Residues in Soluble Proteins*.
              Biophysical Journal. Jan 1985;47(1):61–70.
              :doi:`10.1016/S0006-3495(85)83877-7`. :pmid:`3978191`.
            - Hopp, T. P., and K. R. Woods.
              *Prediction of Protein Antigenic Determinants from Amino Acid
              Sequences*. Proceedings of the National Academy of Sciences
              of the United States of America. Jun 1981;78(6):3824–28.
              :doi:`10.1073/pnas.78.6.3824`. :pmid:`6167991`.
            - Janin, J.
              *Surface and inside Volumes in Globular Proteins*.
              Nature. Feb 1979;277(5696):491–92.
              :doi:`10.1038/277491a0`. :pmid:`763335`.
            - Jones, D. D.
              *Amino Acid Properties and Side-Chain Orientation in Proteins:
              A Cross Correlation Approach*.
              Journal of Theoretical Biology. Mar 1975;50(1):167–83.
              :doi:`10.1016/0022-5193(75)90031-4`. :pmid:`1127956`.
            - Juretić, D., D. Zucić, B. Lucić, and N. Trinajstić.
              *Preference Functions for Prediction of Membrane-Buried
              Helices in Integral Membrane Proteins*.
              Computers & Chemistry. Jun 1998;22(4):279–94.
              :doi:`10.1016/s0097-8485(97)00070-3`. :pmid:`9680689`.
            - Kawashima, S., H. Ogata, and M. Kanehisa.
              *AAindex: Amino Acid Index Database*.
              Nucleic Acids Research. Jan 1999;27(1):368–69.
              :doi:`10.1093/nar/27.1.368`. :pmid:`9847231`.
            - Kawashima, S., and M. Kanehisa.
              *AAindex: Amino Acid Index Database*.
              Nucleic Acids Research. Jan 2000;28(1):374.
              :doi:`10.1093/nar/28.1.374`. :pmid:`10592278`.
            - Kawashima, S., P. Pokarowski, M. Pokarowska, A. Kolinski,
              T. Katayama, and M. Kanehisa.
              *AAindex: Amino Acid Index Database, Progress Report 2008*.
              Nucleic Acids Research. Jan 2008;36:D202-205.
              :doi:`10.1093/nar/gkm998`. :pmid:`17998252`.
            - Kidera, A., Y. Konishi, M. Oka, T. Ooi, and H. A. Scheraga.
              *Statistical Analysis of the Physical Properties of the
              20 Naturally Occurring Amino Acids*.
              Journal of Protein Chemistry. Feb 1985;4(1):23-55.
              :doi:`10.1007/BF01025492`.
            - Kuhn, L. A., C. A. Swanson, M. E. Pique, J. A. Tainer, and
              E. D. Getzoff.
              *Atomic and Residue Hydrophilicity in the Context of Folded
              Protein Structures*. Proteins. Dec 1995;23(4):536–47.
              :doi:`10.1002/prot.340230408`. :pmid:`8749849`.
            - Kyte, J., and R. F. Doolittle.
              *A Simple Method for Displaying the Hydropathic Character of a
              Protein*. Journal of Molecular Biology. May 1982;157(1):105–32.
              :doi:`10.1016/0022-2836(82)90515-0`. :pmid:`7108955`.
            - Levitt, M.
              *A Simplified Representation of Protein Conformations for
              Rapid Simulation of Protein Folding*.
              Journal of Molecular Biology. Jun 1976;104(1):59–107.
              :doi:`10.1016/0022-2836(76)90004-8`. :pmid:`957439`.
            - Manavalan, P., and P. K. Ponnuswamy.
              *Hydrophobic Character of Amino Acid Residues in Globular
              Proteins*. Nature. Oct 1978;275(5681):673–74.
              :doi:`10.1038/275673a0`. :pmid:`703834`.
            - Miyazawa, S., and R. L. Jernigan.
              *Estimation of Effective Interresidue Contact Energies from
              Protein Crystal Structures: Quasi-Chemical Approximation*.
              Macromolecules. Mar 1985;18(3):534–52.
              :doi:`10.1021/ma00145a039`.
            - Nakai, K., A. Kidera, and M. Kanehisa.
              *Cluster Analysis of Amino Acid Indices for Prediction of
              Protein Structure and Function*.
              Protein Engineering. Jul 1988;2(2):93–100.
              :doi:`10.1093/protein/2.2.93`. :pmid:`3244698`.
            - Nozaki, Y., and C. Tanford.
              *The Solubility of Amino Acids and Two Glycine Peptides in
              Aqueous Ethanol and Dioxane Solutions. Establishment of a
              Hydrophobicity Scale*.
              The Journal of Biological Chemistry. Apr 1971;246(7):2211–17.
              :pmid:`5555568`.
            - Parker, J. M., D. Guo, and R. S. Hodges.
              *New Hydrophilicity Scale Derived from High-Performance Liquid
              Chromatography Peptide Retention Data: Correlation of
              Predicted Surface Residues with Antigenicity and X-Ray-Derived
              Accessible Sites*. Biochemistry. 1986;25(19):5425–32.
              :doi:`10.1021/bi00367a013`. :pmid:`2430611`.
            - Ponnuswamy, P. K.
              *Hydrophobic Characteristics of Folded Proteins*. Progress in
              Biophysics and Molecular Biology. 1993;59(1):57–103.
              :doi:`10.1016/0079-6107(93)90007-7`. :pmid:`8419986`.
            - Prabhakaran, M.
              *The Distribution of Physical, Chemical and Conformational
              Properties in Signal and Nascent Peptides*.
              The Biochemical Journal. Aug 1990;269(3):691–96.
              :doi:`10.1042/bj2690691`. :pmid:`2390062`.
            - Rao, J. K. M., and P. Argos.
              *A Conformational Preference Parameter to Predict Helices in
              Integral Membrane Proteins*.
              Biochimica Et Biophysica Acta. Jan 1986;869(2):197–214.
              :doi:`10.1016/0167-4838(86)90295-5`. :pmid:`2935194`.
            - Rose, G. D., A. R. Geselowitz, G. J. Lesser, R. H. Lee, and
              M. H. Zehfus.
              *Hydrophobicity of Amino Acid Residues in Globular Proteins*.
              Science (New York, N.Y.). Aug 1985;229(4716):834–38.
              :doi:`10.1126/science.4023714`. :pmid:`4023714`.
            - Roseman, M. A.
              *Hydrophilicity of Polar Amino Acid Side-Chains Is Markedly
              Reduced by Flanking Peptide Bonds*.
              Journal of Molecular Biology. Apr 1988;200(3):513–22.
              :doi:`10.1016/0022-2836(88)90540-2`. :pmid:`3398047`.
            - Sweet, R. M., and D. Eisenberg.
              *Correlation of Sequence Hydrophobicities Measures Similarity
              in Three-Dimensional Protein Structure*.
              Journal of Molecular Biology. Dec 1983;171(4):479-88.
              :doi:`10.1016/0022-2836(83)90041-4`. :pmid:`6663622`.
            - Tomii, K., and M. Kanehisa.
              *Analysis of Amino Acid Indices and Mutation Matrices for
              Sequence Comparison and Structure Prediction of Proteins*.
              Protein Engineering. Jan 1996;9(1):27–36.
              :doi:`10.1093/protein/9.1.27`. :pmid:`9053899`.
            - Welling, G. W., W. J. Weijer, R. van der Zee R, and
              S. Welling-Wester.
              *Prediction of Sequential Antigenic Regions in Proteins*.
              FEBS Letters. Feb 1985;188(2):215-8.
              :doi:`10.1016/0014-5793(85)80374-4`. :pmid:`2411595`.
            - White, S. H., and W. C. Wimley.
              *Membrane Protein Folding and Stability: Physical Principles*.
              Annual Review of Biophysics and Biomolecular Structure.
              1999;28:319–65.
              :doi:`10.1146/annurev.biophys.28.1.319`. :pmid:`10410805`
            - White, S. H., and W. C. Wimley.
              *Hydrophobic Interactions of Peptides with Membrane Interfaces*.
              Biochimica Et Biophysica Acta. Nov 1998;1376(3):339-52.
              :doi:`10.1016/s0304-4157(98)00021-5`. :pmid:`9804985`.
            - Wilson, K. J., A. Honegger, R. P. Stötzel, and G. J. Hughes.
              *The Behaviour of Peptides on Reverse-Phase Supports during
              High-Pressure Liquid Chromatography*.
              The Biochemical Journal. Oct 1981;199(1):31-41.
              :doi:`10.1042/bj1990031`. :pmid:`7337711`.
            - Wimley, W. C., and S. H. White.
              *Experimentally Determined Hydrophobicity Scale for Proteins
              at Membrane Interfaces*.
              Nature Structural Biology. Oct 1996;3(10):842–48.
              :doi:`10.1038/nsb1096-842`. :pmid:`8836100`.
            - Wimley, W. C., T. P. Creamer, and S. H. White.
              *Solvation Energies of Amino Acid Side Chains and Backbone in
              a Family of Host-Guest Pentapeptides*.
              Biochemistry. Apr 1996;35(16):5109–24.
              :doi:`10.1021/bi9600153`. :pmid:`8611495`.
            - Wolfenden, R., L. Andersson, P. M. Cullis, and C. C. Southgate.
              *Affinities of Amino Acid Side Chains for Solvent Water*.
              Biochemistry. Feb 1981;20(4):849–55.
              :doi:`10.1021/bi00507a030`. :pmid:`7213619`.
            - Zimmerman, J. M., N. Eliezer, and R. Simha.
              *The Characterization of Amino Acid Sequences in Proteins by
              Statistical Methods*.
              Journal of Theoretical Biology. Nov 1968;21(2):170–201.
              :doi:`10.1016/0022-5193(68)90069-6`. :pmid:`5700434`.

        """
        table = tables.HYDROPHOBICITY.get(scale)
        if table is None:
            raise ValueError(f"Invalid hydrophobicity scale: {scale!r}")
        return _sum(self.profile(table)) / len(self)

    def instability_index(self) -> float:
        """Compute the instability index of a protein sequence.

        This function calculates the instability index proposed by
        Guruprasad *et al* (1990). This index predicts the stability of a
        protein based on its dipeptide composition.

        Returns:
            `float`: The instability index of the peptide. A protein whose
            instability index is smaller than 40 is predicted as stable, a
            value above 40 predicts that the protein may be unstable.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> round(peptide.instability_index(), 2)
            83.68

        References:
            - Guruprasad, K., B.V. Bhasker Reddy, and M. W. Pandit.
              *Correlation between Stability of a Protein and Its Dipeptide
              Composition: A Novel Approach for Predicting in Vivo
              Stability of a Protein from Its Primary Sequence*. Protein
              Engineering, Design and Selection. 1990 Dec;4(2):155–61.
              :doi:`10.1093/protein/4.2.155`. :pmid:`2075190`.

        """
        scale = tables.INSTABILITY["Guruprasad"]
        gp = sum(scale.get(self.sequence[i : i + 2], 1.0) for i in range(len(self.sequence) - 1))
        return gp * 10 / (len(self.sequence))

    def isoelectric_point(self, pKscale: str = "EMBOSS") -> float:
        """Compute the isoelectric point of a protein sequence.

        The isoelectric point (*pI*), is the *pH* at which a particular
        molecule or surface carries no net electrical charge.

        Arguments:
            pKscale (`str`): The name of the pKa scale to be used. A list of
                all the allowed values can be retrieved from the keys of the
                `peptides.tables.PK` dictionary.

        Returns:
            `float`: The pH at which the peptide has a neutral net charge.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> peptide.isoelectric_point(pKscale="EMBOSS")
            9.71...
            >>> peptide.isoelectric_point(pKscale="Murray")
            9.81...
            >>> peptide.isoelectric_point(pKscale="Sillero")
            9.89...
            >>> peptide.isoelectric_point(pKscale="Solomon")
            9.58...
            >>> peptide.isoelectric_point(pKscale="Stryer")
            9.62...
            >>> peptide.isoelectric_point(pKscale="Lehninger")
            9.93...
            >>> peptide.isoelectric_point(pKscale="Dawson")
            9.56...
            >>> peptide.isoelectric_point(pKscale="Rodwell")
            9.71...

        Note:
            The pI is a variable that affects the solubility of the peptides
            under certain conditions of pH. When the pH of the solvent is
            equal to the pI of the protein, it tends to precipitate and lose
            its biological function.

        References:
            - Rice, P., I. Longden, and A. Bleasby.
              *EMBOSS: The European Molecular Biology Open Software Suite*.
              Trends in Genetics. June 2000;16(6):276–77.
              :doi:`10.1016/s0168-9525(00)02024-2`. :pmid:`10827456`

        """
        # use a simple bissecting loop to minimize the charge function
        top, bottom, x = 0.0, 14.0, 7.0
        while not math.isclose(top, bottom):
            x = (top + bottom) / 2
            c = self.charge(pH=x, pKscale=pKscale)
            if c >= 0:
                top = x
            if c <= 0:
                bottom = x
        return x

    def mass_shift(
        self,
        aa_shift: typing.Union[str, typing.Dict[str, float], None] = "silac_13c",
        monoisotopic: bool = True,
    ) -> float:
        """Compute the mass difference of modified peptides.

        This function calculates the mass difference of peptides introduced
        by chemical modifications or heavy isotope labelling.

        Arguments:
            aa_shift (`str` or `dict`): Either the key to a pre-defined
                isotope label (see `peptides.tables.MASS_SHIFT`), or a
                dictionary mapping each amino acid to it mass difference
                in Dalton (use ``nTer`` and ``cTer`` keys for N-terminal
                and C-terminal modifications).
            monoisotopic (`bool`): Flag whether monoisotopic weights of
                amino-acids should be used.

        Return:
            `float`: The mass difference of the modified peptide.

        Example:
            >>> peptide = Peptide("EGVNDNECEGFFSAR")
            >>> peptide.mass_shift(aa_shift="silac_13c")
            6.020129...
            >>> peptide.mass_shift(aa_shift=dict(R=10.00827))
            10.00827...

        References:
            - Ong, S-E., I. Kratchmarova, and M. Mann.
              *Properties of 13C-Substituted Arginine in Stable Isotope
              Labeling by Amino Acids in Cell Culture (SILAC)*.
              Journal of Proteome Research. Apr 2003;2(2):173–81.
              :doi:`10.1021/pr0255708`. :pmid:`12716131`.
            - Picotti, P., B. Bodenmiller, L. N. Mueller, B. Domon,
              and R. Aebersold.
              *Full Dynamic Range Proteome Analysis of S. Cerevisiae by
              Targeted Proteomics*. Cell. Aug 2009;138(4):795–806.
              :doi:`10.1016/j.cell.2009.05.051`. :pmid:`19664813`.

        """
        if isinstance(aa_shift, str):
            table = tables.MASS_SHIFT.get(aa_shift)
            if table is None:
                raise ValueError(f"Invalid mass shift scale: {aa_shift!r}")
            scale = {}
            if aa_shift == "silac_13c":
                scale["K"] = table["K"] - 0.064229 * (not monoisotopic)
                scale["R"] = table["R"] - 0.064229 * (not monoisotopic)
            elif aa_shift == "silac_13c15n":
                scale["K"] = table["K"] - 0.071499 * (not monoisotopic)
                scale["R"] = table["R"] - 0.078669 * (not monoisotopic)
            elif aa_shift == "15n":
                for k, v in table.items():
                    scale[k] = v * 0.997035 - 0.003635 * (not monoisotopic)
        elif isinstance(aa_shift, dict):
            scale = aa_shift
        else:
            raise TypeError(
                f"Expected str or dict, found {aa_shift.__class__.__name__}"
            )

        # sum the mass-shift of each amino acid
        shift = _sum(self.profile(scale))
        # return the shift with C-terminal and N-terminal ends
        return shift + scale.get("nTer", 0.0) + scale.get("cTer", 0.0)

    def molecular_weight(
        self,
        average: str = "expasy",
        aa_shift: typing.Union[str, typing.Dict[str, float], None] = None,
    ) -> float:
        """Compute the molecular weight of a protein sequence.

        This function calculates the molecular weight of a protein sequence.
        It is calculated as the sum of the mass of each amino acid using
        one of the 3 available scales. It also supports mass calculation
        of proteins with predefined or custom stable isotope mass labels.

        Arguments:
            average (`str`): The name of the average amino acid average
                weight scale. See `peptides.tables.MOLECULAR_WEIGHT` for a
                list of appropriate values.
            aa_shift (`str`, `dict` or `None`): Either an appropriate shift
                value to pass to `Peptide.mass_shift`, or `None` to get the
                unmodified weight.

        Return:
            `float`: The molecular weight of the peptide, in Dalton.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> peptide.molecular_weight()
            2485.91...
            >>> peptide.molecular_weight(average="mascot")
            2485.89...
            >>> peptide.molecular_weight(average="monoisotopic")
            2484.11...

        References:
            - Wilkins, M. R., E. Gasteiger, A. Bairoch, J. C. Sanchez,
              K. L. Williams, R. D. Appel, and D. F. Hochstrasser.
              *Protein Identification and Analysis Tools in the ExPASy
              Server*. Methods in Molecular Biology. 1992;112: 531–52.
              :doi:`10.1385/1-59259-584-7:531`. :pmid:`10027275`

        """
        scale = tables.MOLECULAR_WEIGHT.get(average)
        if scale is None:
            raise ValueError(f"Invalid average weight scale: {average!r}")

        # sum the weight of each amino acid and add weight of water molecules
        mass = _sum(self.profile(scale)) + scale["H2O"]
        # add mass shift for labeled proteins
        if aa_shift is not None:
            mass += self.mass_shift(
                aa_shift=aa_shift, monoisotopic=average == "monoisotopic"
            )

        return mass

    def mz(
        self,
        charge: int = 2,
        aa_shift: typing.Union[str, typing.Dict[str, float], None] = None,
        cysteins: float = 57.021464,
    ) -> float:
        """Compute the m/z (mass over charge) ratio for a peptide.

        This function calculates the (monoisotopic) mass over charge ratio
        (m/z) for peptides, as measured in mass spectrometry.

        Arguments:
            charge (`int`): The net charge for which the m/z should be
                computed.
            aa_shift (`str`, `dict` or `None`): Either an appropriate shift
                value to pass to `Peptide.mass_shift`, or `None` to get the
                unmodified weight.
            cysteins (`float`): The mass shift (in Dalton) of blocked
                cysteins. Default corresponds to cysteins blocked by
                iodoacetamide.

        Returns:
            `float`: The m/z ratio of the peptide.

        Example:
            >>> peptide = Peptide("EGVNDNECEGFFSAR")
            >>> peptide.mz()
            865.857...
            >>> peptide.mz(aa_shift=dict(K=6.020129, R=6.020129))
            868.867...
            >>> peptide.mz(aa_shift="silac_13c", cysteins=58.005479)
            869.359...

        """
        if not isinstance(charge, int):
            raise TypeError(f"Expected int, found {charge.__class__.__name__!r}")

        # compute the mass of the uncharged peptide
        mass = self.molecular_weight(average="monoisotopic", aa_shift=aa_shift)
        # add modification at cysteins
        mass += self.sequence.count("C") * cysteins
        # modify for charged peptides
        if charge >= 0:
            mass += charge * 1.007276  # weights of H+1 ions
            mass /= charge  # divide by charge state

        return mass

    # --- Profiles -----------------------------------------------------------

    def hydrophobicity_profile(
        self, window: int = 11, scale: str = "KyteDoolittle"
    ) -> typing.Sequence[float]:
        """Build a hydrophobicity profile of a sliding window.

        Example:
            >>> peptide = Peptide("ARQQNLFINFCLILIFLLLI")
            >>> h = peptide.hydrophobicity_profile(window=12, scale="Eisenberg")
            >>> [round(x, 3) for x in h]
            [0.083, 0.147, 0.446, 0.632, 0.802, 0.955, 0.955, 0.944, 0.944]

        """
        if scale not in tables.HYDROPHOBICITY:
            raise ValueError(f"Invalid hydrophobicity scale: {scale!r}")

        profile = array.array('f')
        for i in range(len(self.sequence) - window + 1):
            profile.append(self[i : i + window].hydrophobicity(scale=scale))

        return profile

    def hydrophobic_moment_profile(
        self, window: int = 11, angle: int = 100
    ) -> typing.Sequence[float]:
        """Build a hydrophobic moment profile of a sliding window.

        This function builds a profile computing the hydrophobic moment of
        a section of the peptide based on the primary sequecne.

        Arguments:
            window (`int`): The size of the sliding window for which to
                compute the local hydrophobic moment.
            angle (`int`): A protein rotational angle, in **degrees**.
                Usual values are *100* for α-helix, and *160* for β-sheet.

        Example:
            >>> peptide = Peptide("ARQQNLFINFCLILIFLLLI")
            >>> uH = peptide.hydrophobic_moment_profile(window=12, angle=100)
            >>> [round(x, 3) for x in uH]
            [0.353, 0.317, 0.274, 0.274, 0.253, 0.113, 0.113, 0.108, 0.132]

        See Also:
            The `~Peptide.hydrophobic_moment` method, which computes the
            maximal hydrophobic moment instead of building a profile.

        """
        profile = array.array("f")
        for i in range(len(self.sequence) - window + 1):
            profile.append(
                self[i : i + window].hydrophobic_moment(window=window - 1, angle=angle)
            )

        return profile

    def membrane_position_profile(
        self, window: int = 11, angle: int = 100
    ) -> typing.List[str]:
        """Compute the theoretical class of a protein sequence.

        This function builds a profile predicting the theoretical class
        of a section of the peptide based on the relationship between the
        hydrophobic moment and hydrophobicity scale as proposed by
        Eisenberg (1984).

        Arguments:
            window (`int`): The window size to consider when building the
                profile.
            angle (`int`): The protein rotational angle, in **degrees**, for
                which to compute the hydrophobic moment profile. Usual
                values are *100* for α-helix, and *160* for β-sheet.

        Returns:
            `list` of `str`: A list containing a one-character code for
            each window starting position: either `'G'` for globular,
            `'T'` for transmembrane, or `'S'` for surface.

        Example:
            >>> peptide = Peptide("ARQQNLFINFCLILIFLLLI")
            >>> peptide.membrane_position_profile(window=12, angle=100)
            ['G', 'G', 'G', 'T', 'S', 'T', 'T', 'T', 'T']
            >>> peptide.membrane_position_profile(window=12, angle=160)
            ['G', 'G', 'G', 'S', 'S', 'S', 'S', 'S', 'S']

        References:
            - Eisenberg, D.
              *Three-Dimensional Structure of Membrane and Surface Proteins*.
              Annual Review of Biochemistry. July 1984;53:595–623.
              :doi:`10.1146/annurev.bi.53.070184.003115`. :pmid:`6383201`.
            - Eisenberg, D., E. Schwarz, M. Komaromy, and R. Wall.
              *Analysis of Membrane and Surface Protein Sequences with
              the Hydrophobic Moment Plot*.
              Journal of Molecular Biology. Oct 1984;179(1):125–42.
              :doi:`10.1016/0022-2836(84)90309-7`. :pmid:`6502707`.
            - Eisenberg, D., R. M. Weiss, and T. C. Terwilliger.
              *The Helical Hydrophobic Moment: A Measure of the
              Amphiphilicity of a Helix*. Nature. Sep 1982;299(5881):371–74.
              :doi:`10.1038/299371a0`. :pmid:`7110359`

        """
        profile_H = self.hydrophobicity_profile(window=window, scale="Eisenberg")
        profile_uH = self.hydrophobic_moment_profile(window=window, angle=angle)

        profile = []
        for h, uh in zip(profile_H, profile_uH):
            m = h * -0.421 + 0.579
            if uh <= m:
                if h >= 0.5:
                    profile.append("T")
                else:
                    profile.append("G")
            else:
                profile.append("S")

        return profile

    # --- Structural properties ----------------------------------------------

    def structural_class(
        self,
        frequencies: str = "Nakashima",
        distance: str = "mahalanobis",
    ) -> str:
        """Predict the structural class of the peptide from its sequence.

        The structural class of a protein, as defined in Levitt and Chothia
        (1976), can be either α, β, α+β, or α/β, with ζ being later defined
        for irregular proteins. It depends on the secondary structure of the
        protein. Several methods have been proposed to elucidate the
        structural class from the amino acid sequence, all based on
        similarity with proteins which structures have been elucidated.

        Arguments:
            frequencies (`str`): The frequencies of the amino acids in
                proteins of different structural classes to use as reference
                centroids. Use `"Chou"` to load the frequencies of the 64
                proteins analyzed in Chou (1989), `"Nakashima"` to use
                the normalized frequencies of the 135 proteins analyzed in
                Nakashima *et al.* (1986) and Zhang & Chou (1995), or
                `"ChouZhang"` to load the frequencies of 120 proteins used
                in Chou & Zhang (1995).
            distance (`str`): The distance metric to use in the 20-D space
                formed by the 20 usual amino acid to find the nearest
                structural class for the peptide. Use `"cityblock"` to use
                the Manhattan distance like in Chou (1989), `"euclidean"`
                to use the Euclidean distance like in Nakashima *et al*
                (1986), `"correlation"` to use the correlation distance
                like in Chou & Zhang (1992), `"mahalanobis"` to use
                the Mahalanobis distance like in Chou & Zhang (1995),
                or `"discriminant"` to use the Bayes discriminant like in
                Chou *et al.* (1998).

        Returns:
            `str`: The structural class the protein most likely belongs to.
            Note that some classes may not be predictable, depending on the
            reference frequencies being used (at the moment, the ζ class
            can only be predicted from the *Nakashima* frequencies with
            *euclidean* or *manhattan* distances).

        Example:
            Predict the structural class of the skipjack tuna Cytochrome C,
            (`P0025 <https://www.uniprot.org/uniprot/P00025>`_), an α
            protein ::

                >>> p = Peptide(
                ...     "MGDVAKGKKTFVQKCAQCHTVENGGKHKVGPNLWGLFGRKTGQAEGYSYT"
                ...     "DANKSKGIVWNENTLMEYLENPKKYIPGTKMIFAGIKKKGERQDLVAYLK"
                ...     "SATS"
                ... )
                >>> p.structural_class("Nakashima", distance="mahalanobis")
                'alpha'
                >>> p.structural_class("ChouZhang", distance="mahalanobis")
                'beta'
                >>> p.structural_class("Chou", distance="correlation")
                'alpha'
                >>> p.structural_class("Nakashima", distance="euclidean")
                'alpha'
                >>> p.structural_class("Chou", distance="cityblock")
                'alpha+beta'

            Predict the structural class of the sea krait Erabutoxin B
            (`Q90VW1 <https://www.uniprot.org/uniprot/Q90VW1>`_), a β
            protein::

                >>> p = Peptide(
                ...     "MKTLLLTLVVVTIVCLDLGYTRICFNHQSSQPQTTKTCSPGESSCYHKQW"
                ...     "SDFRGTIIERGCGCPTVKPGIKLSCCESEVCNN"
                ... )
                >>> p.structural_class("Nakashima", distance="mahalanobis")
                'beta'
                >>> p.structural_class("ChouZhang", distance="mahalanobis")
                'alpha+beta'
                >>> p.structural_class("Chou", distance="correlation")
                'beta'
                >>> p.structural_class("Nakashima", distance="euclidean")
                'zeta'
                >>> p.structural_class("Chou", distance="cityblock")
                'beta'

            Predict the structural class of the *Arthrospira platensis*
            Ferredoxin (`P00246 <https://www.uniprot.org/uniprot/P00246>`_),
            a ζ protein::

                >>> p = Peptide(
                ...     "MATYKVTLINEAEGINETIDCDDDTYILDAAEEAGLDLPYSCRAGACSTC"
                ...     "AGTITSGTIDQSDQSFLDDDQIEAGYVLTCVAYPTSDCTIKTHQEEGLY"
                ... )
                >>> p.structural_class("Nakashima", distance="euclidean")
                'zeta'

        References:
            - Chou, K-C., W-M. Liu, G. M. Maggiora, and C-T. Zhang.
              *Prediction and Classification of Domain Structural Classes*.
              Proteins: Structure, Function, and Genetics.
              Apr 1998;31(1):97–103. :pmid:`9552161`.
            - Chou, K-C., and C-T. Zhang.
              *Prediction of Protein Structural Classes*. Critical Reviews
              in Biochemistry and Molecular Biology. Feb 1995;30:275–349.
              :doi:`10.3109/10409239509083488`. :pmid:`7587280`.
            - Chou, K-C., and C-T. Zhang.
              *A Correlation-Coefficient Method to Predicting
              Protein-Structural Classes from Amino Acid Compositions*.
              European Journal of Biochemistry. 1992;207(2):429–33.
              :doi:`10.1111/j.1432-1033.1992.tb17067.x`. :pmid:`1633801`.
            - Chou, P. Y.
              *Prediction of Protein Structural Classes from Amino Acid
              Compositions*. In Prediction of Protein Structure and the
              Principles of Protein Conformation, edited by G. D. Fasman.
              Springer US. 1989:549–86.
              :doi:`10.1007/978-1-4613-1571-1`. ISBN:978-0-306-43131-9.
            - Nakashima, H., K. Nishikawa, and T. Ooi.
              *The Folding Type of a Protein Is Relevant to the Amino Acid
              Composition*. Journal of Biochemistry. Jan 1986;99(1):153–62.
              :doi:`10.1093/oxfordjournals.jbchem.a135454`. :pmid:`3957893`.
            - Zhang, Chun-Ting, and Kuo-Chen Chou.
              *An Eigenvalue-Eigenvector Approach to Predicting Protein
              Folding Types*.
              Journal of Protein Chemistry. Jul 1995;14(5):309–26.
              :doi:`10.1007/BF01886788`. :pmid:`8590599`.
            - Zhou, G.P., and N. Assa-Munt.
              *Some Insights into Protein Structural Class Prediction*.
              Proteins: Structure, Function, and Bioinformatics.
              2001;44(1):57–59. :doi:`10.1002/prot.1071`. :pmid:`11354006`.

        """
        # get peptide frequencies
        pep_frequencies = self.frequencies()

        # get reference frequencies
        if frequencies == "Chou":
            dataset = datasets.CHOU1989
        elif frequencies == "Nakashima":
            dataset = datasets.NAKASHIMA1986
            # Nakashima frequencies are normalized, so we must normalize
            # the peptide frequencies too in that case
            mean = datasets.NAKASHIMA1986["all"]["mean"]
            sd = datasets.NAKASHIMA1986["all"]["sd"]
            pep_frequencies = {
                aa:(pep_frequencies[aa]-mean[aa])/sd[aa] for aa in mean
            }
        elif frequencies == "ChouZhang":
            dataset = datasets.CHOUZHANG1995
        else:
            raise ValueError(f"Invalid amino acid frequencies: {frequencies!r}")

        # get distances to each structural class centroids
        distances = {}
        if distance == "correlation":
            for name,tables in dataset.items():
                if name == "all":
                    continue
                mean = tables["mean"]
                s1 = sum(pep_frequencies[x]*mean[x] for x in mean)
                s2 = sum(pep_frequencies[x]**2 for x in mean)
                s3 = sum(mean[x]**2 for x in mean)
                distances[name] = 1 - s1 / math.sqrt(s2*s3)
        elif distance == "cityblock":
            for name,tables in dataset.items():
                if name == "all":
                    continue
                table = tables["mean"]
                s = sum(abs(pep_frequencies[x]-table[x]) for x in table)
                distances[name] = s
        elif distance == "euclidean":
            for name,tables in dataset.items():
                if name == "all":
                    continue
                table = tables["mean"]
                s = sum((pep_frequencies[x]-table[x])**2 for x in table)
                distances[name] = math.sqrt(s)
        elif distance == "mahalanobis" or distance == "discriminant":
            x = [pep_frequencies[aa] for aa in sorted(_CODE1[:20])]
            for name,tables in dataset.items():
                if name == "all":
                    continue
                # load parameters
                mean = tables["mean"]
                eivals = tables.get("eigenvalues")
                eivecs = tables.get("eigenvectors")
                if eivals is None or eivecs is None:
                    continue
                # compute Mahalanobis distance using the components of
                # the eigendecomposition
                xm = [mean[aa] for aa in sorted(_CODE1[:20])]
                y = [
                    sum(eivecs[i][j] * (x[j] - xm[j]) for j in range(20))
                    for i in range(20)
                ]
                distances[name] = sum(y[i]**2 / eivals[i] for i in range(1, 20))
                # if Bayes discriminant is requested, add the logarithm of
                # the product of all non-null eigenvalues
                if distance == "discriminant":
                    product = functools.reduce(operator.prod, einvals[1:])
                    distances[name] += math.log(product)
            if not distances:
                raise ValueError(
                    f"Cannot use {frequencies!r} frequencies with "
                    f"{distance!r} distance (no eigendecomposition available)"
                )
        else:
            raise ValueError(f"Invalid distance: {distance!r}")

        # find the most likely structural class based on the distance
        return min(distances, key=distances.get)

    def linker_preference_profile(
        self, window: int = 15
    ) -> typing.Sequence[float]:
        """Compute the linker preference profile of a protein sequence.

        The linker preference profile is a measure used as a basis for the
        DomCut method in Suyama & Ohara (2002). The resulting profile can
        then be used to identify putative domain boundaries in the input
        protein, either:

        - Using prior knowledge of the estimated domain count :math:`D`, in
          which case the :math:`D-1` global minimums in the sequence can
          be used as cutting points
        - Without prior knowledge of the domain count, using a fixed
          threshold to estimate the number of domains and linkers. A
          cutoff value of :math:`-0.09` was selected by the authors
          optimizing on the specificity / selectivity tradeoff.

        References:
            - Suyama, M., O. Ohara O.
              *DomCut: prediction of inter-domain linker regions in amino
              acid sequences.* Bioinformatics. Mar 2003;19(5):673-4.
              :doi:`10.1093/bioinformatics/btg031`. :pmid:`12651735`.

        .. versionadded:: 0.3.0

        """
        return self.profile(tables.LINKER_INDEX["Suyama"], window=window)

    # --- Biosynthetic costs -------------------------------------------------

    def nutrient_cost(
        self,
        nutrient: str = "glucose",
        organism: str = "yeast",
        relative: bool = False,
    ):
        """Estimate the nutrient cost to biosynthesize a peptide.

        The nutrient cost was proposed by Barton *et al.* to estimate
        the energy cost required to biosynthesize a peptide based on
        genome-scale metabolic modeling in *Saccharomyces cerevisiae* and
        *Escherichia coli*. This approach offers advantages to estimate
        costs in nutrient-limited environments compared to energy-based
        methods.

        Arguments:
            nutrient (`str`): The name of the nutrient, one of ``glucose``,
                ``sulphate`` or ``ammonia``.
            organism (`str`): The reference organism for which the values
                were computed, either ``yeast`` or ``ecoli``.
            relative (`bool`): Whether to use the absolute or relative
                values when computing costs.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLE")
            >>> peptide.nutrient_cost("glucose", "ecoli")
            15.763...
            >>> peptide.nutrient_cost("ammonia", "yeast", relative=True)
            10.839...

        References:
            - Barton, M. D., Delneri, D., Oliver, S. G., Rattray, M.,
              & Bergman, C. M. (2010). *Evolutionary systems biology of amino
              acid biosynthetic cost in yeast*. PloS one, 5(8), e11935.
              :pmid:`20808905` :doi:`10.1371/journal.pone.0011935`.

        """
        if nutrient == "glucose":
            n = "car"
        elif nutrient == "ammonia":
            n = "nit"
        elif nutrient == "sulphate":
            n = "sul"
        else:
            raise ValueError(f"invalid nutrient: {nutrient!r}")

        if organism not in ("ecoli", "yeast"):
            raise ValueError(f"invalid organism: {organism!r}")

        r = "rel" if relative else "abs"
        table = f"{organism}_{n}_{r}"

        p = self.profile(tables.NUTRIENT_COST[table])
        return _sum(p)

    def energy_cost(
        self,
        scale: str = "Akashi",
        *,
        mode: typing.Optional[str] = None,
    ):
        """Estimate the energy cost required to biosynthesize a peptide.

        Arguments:
            scale (`str`): The name of the energy estimation scale to use.
                Supports the following values:

                Akashi
                    The energetic cost computed by Akashi & Gojobori (2002) 
                    based on major codon usage values in *Escherichia coli* 
                    and *Bacillus subtilis*.
                Craig
                    The energetic cost computed by Craig & Weber ()
                    from amino-acid substitution probabilities in 
                    *Escherichia coli*.
                Heizer
                    The energetic cost computed by Heizer *et al.* (2006), 
                    derived from Akashi & Gojobori (2002) for photoautotrophs 
                    (capable of the Calvin cycle reactions).
                Wagner
                    The energetic cost computed by Wagner (2005)
                    based on expression data in *Saccharomyces cerevisiae*.

        Keyword Arguments:
            mode (`str`): For the ``Wagner`` scale, the mode of growth of 
                the source organism, either ``respiration`` (the default) or 
                ``fermentation``.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLEYLKW")
            >>> peptide.energy_cost()
            550.5...
            >>> peptide.energy_cost("Heizer")
            554.5...

        References:
            - Akashi, H., & Gojobori, T. (2002). *Metabolic efficiency and
              amino acid composition in the proteomes of Escherichia coli and
              Bacillus subtilis*. Proceedings of the National Academy of
              Sciences of the United States of America, 99(6), 3695–3700.
              :pmid:`11904428`. :doi:`10.1073/pnas.062526999`.
            - Heizer, E. M., Jr, Raiford, D. W., Raymer, M. L., Doom, T. E.,
              Miller, R. V., & Krane, D. E. (2006). *Amino acid cost and
              codon-usage biases in 6 prokaryotic genomes: a whole-genome
              analysis*. Molecular biology and evolution, 23(9), 1670–1680.
              :pmid:`16754641`. :doi:`10.1093/molbev/msl029`.
            - Wagner A. (2005). *Energy constraints on the evolution of gene
              expression*. Molecular biology and evolution, 22(6), 1365–1374.
              :pmid:`15758206`. :doi:`10.1093/molbev/msi126`.

        """
        if scale == "Wagner":
            if mode == "respiration" or mode is None:
                table = tables.ENERGY_COST["Wagner_resp"]
            elif mode == "fermentation":
                table = tables.ENERGY_COST["Wagner_ferm"]
            else:
                raise ValueError(f"invalid lifestyle for {scale!r} scale: {mode!r}")
        elif scale in tables.ENERGY_COST:
            table = tables.ENERGY_COST[scale]
        else:
            raise ValueError(f"invalid scale: {scale!r}")

        p = self.profile(table)
        return _sum(p)


    # --- Descriptors --------------------------------------------------------

    def blosum_indices(self) -> BLOSUMIndices:
        """Compute the BLOSUM62-derived indices of the peptide.

        See `~peptides.BLOSUMIndices` for more information.

        Returns:
            `peptides.BLOSUMIndices`: The computed average BLOSUM indices
            for all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("KLKLLLLLKLK")
            >>> for i, b in enumerate(peptide.blosum_indices()):
            ...     print(f"BLOSUM{i+1:<3} {b: .4f}")
            BLOSUM1   -0.4827
            BLOSUM2   -0.5618
            BLOSUM3   -0.8509
            BLOSUM4   -0.4173
            BLOSUM5    0.3173
            BLOSUM6    0.2527
            BLOSUM7    0.1464
            BLOSUM8    0.1427
            BLOSUM9   -0.2145
            BLOSUM10  -0.3218

        """
        out = []
        for i in range(len(tables.BLOSUM)):
            p = self.profile(tables.BLOSUM[f"BLOSUM{i+1}"])
            out.append(_sum(p) / len(self))
        return BLOSUMIndices(*out)

    def cruciani_properties(self) -> CrucianiProperties:
        """Compute the Cruciani properties of the peptide.

        See `~peptides.CrucianiProperties` for more information.

        Returns:
            `peptides.CrucianiProperties`: The computed average Cruciani
            properties of all the amino acids in the corresponding peptide
            sequence.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, b in enumerate(peptide.cruciani_properties()):
            ...     print(f"PP{i+1:<3} {b: .4f}")
            PP1   -0.1130
            PP2   -0.0220
            PP3    0.2735

        """
        out = []
        for i in range(len(tables.CRUCIANI)):
            p = self.profile(tables.CRUCIANI[f"PP{i+1}"])
            out.append(_sum(p) / len(self))
        return CrucianiProperties(*out)

    def fasgai_vectors(self) -> FasgaiVectors:
        """Compute the FASGAI vectors of the peptide.

        See `~peptides.FasgaiVectors` for more information.

        Returns:
            `peptides.FasgaiVectors`: The computed average FASGAI vectors
            for all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, b in enumerate(peptide.fasgai_vectors()):
            ...     print(f"F{i+1:<3} {b: .5f}")
            F1   -0.13675
            F2   -0.45485
            F3   -0.11695
            F4   -0.45800
            F5   -0.38015
            F6    0.52740

        """
        out = []
        for i in range(len(tables.FASGAI)):
            p = self.profile(tables.FASGAI[f"F{i+1}"])
            out.append(_sum(p) / len(self))
        return FasgaiVectors(*out)

    def kidera_factors(self) -> KideraFactors:
        """Compute the Kidera factors of the peptide.

        See `~peptides.KideraFactors` for more information.

        Returns:
            `peptides.KideraFactors`: The computed average Kidera factors
            for all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("KLKLLLLLKLK")
            >>> for i, kf in enumerate(peptide.kidera_factors()):
            ...     print(f"KF{i+1:<3} {kf: .4f}")
            KF1   -0.7855
            KF2    0.2982
            KF3   -0.2364
            KF4   -0.0818
            KF5    0.2100
            KF6   -1.8936
            KF7    1.0291
            KF8   -0.5127
            KF9    0.1118
            KF10   0.8100

        """
        out = []
        for i in range(len(tables.KIDERA)):
            p = self.profile(tables.KIDERA[f"KF{i+1}"])
            out.append(_sum(p) / len(self))
        return KideraFactors(*out)

    def atchley_factors(self) -> AtchleyFactors:
        """Compute the Atchley factors of the peptide.

        See `~peptides.AtchleyFactors` for more information.

        Returns:
            `peptides.AtchleyFactors`: The computed average Atchley factors
            for all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("KLKLLLLLKLK")
            >>> for i, kf in enumerate(peptide.atchley_factors()):
            ...     print(f"AF{i+1:<3} {kf: .4f}")
            AF1    0.0176
            AF2   -0.8321
            AF3   -0.7636
            AF4    0.7048
            AF5    0.0189

        """
        out = []
        for i in range(len(tables.ATCHLEY)):
            p = self.profile(tables.ATCHLEY[f"AF{i+1}"])
            out.append(_sum(p) / len(self))
        return AtchleyFactors(*out)

    def ms_whim_scores(self) -> MSWHIMScores:
        """Compute the MS-WHIM scores of the peptide.

        See `~peptides.MSWHIMScores` for more information.

        Returns:
            `peptides.MSWHIMScores`: The compute average of MS-WHIM scores
            of all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("KLKLLLLLKLK")
            >>> for i, mw in enumerate(peptide.ms_whim_scores()):
            ...     print(f"MSWHIM{i+1:<3} {mw: .4f}")
            MSWHIM1   -0.6564
            MSWHIM2    0.4873
            MSWHIM3    0.1164

        """
        out = []
        for i in range(len(tables.MSWHIM)):
            p = self.profile(tables.MSWHIM[f"MSWHIM{i+1}"])
            out.append(_sum(p) / len(self))
        return MSWHIMScores(*out)

    def pcp_descriptors(self) -> PCPDescriptors:
        """Compute the Physical-Chemical Properties descriptors of the peptide.

        See `~peptides.PCPDescriptors` for more information.

        Returns:
            `peptides.PCPDescriptors`: The computed average of PCP
            descriptors of all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, pcp in enumerate(peptide.pcp_descriptors()):
            ...     print(f"E{i+1:<3} {pcp: .5f}")
            E1    0.01090
            E2    0.03810
            E3    0.12505
            E4    0.04095
            E5   -0.10595

        """
        out = []
        for i in range(len(tables.PCP_DESCRIPTORS)):
            p = self.profile(tables.PCP_DESCRIPTORS[f"E{i+1}"])
            out.append(_sum(p) / len(self))
        return PCPDescriptors(*out)

    def physical_descriptors(self) -> PhysicalDescriptors:
        """Compute the Physical Descriptors of the peptide.

        See `~peptides.PhysicalDescriptors` for more information.

        Returns:
            `peptides.PhyiscalDescriptors`: The computed average of Physical
            Descriptors of all the amino acids in the peptide. *PD1* is
            related to volume while *PD2* is related to hydrophilicity.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, pd in enumerate(peptide.physical_descriptors()):
            ...     print(f"PD{i+1:<3} {pd: .4f}")
            PD1    0.1190
            PD2    0.2825

        """
        out = []
        for i in range(len(tables.PHYSICAL_DESCRIPTORS)):
            p = self.profile(tables.PHYSICAL_DESCRIPTORS[f"PD{i+1}"])
            out.append(_sum(p) / len(self))
        return PhysicalDescriptors(*out)

    def prin_components(self) -> PRINComponents:
        """Compute the PRIN components of the peptide.

        See `~peptides.PRINComponents` for more information.

        Returns:
            `peptides.PRINComponents`: The computed average of PRIN
            components of all the amino acids in the peptide.

        """
        out = []
        for i in range(len(tables.PRIN)):
            p = self.profile(tables.PRIN[f"PRIN{i+1}"])
            out.append(_sum(p) / len(self))
        return PRINComponents(*out)

    def protfp_descriptors(self) -> ProtFPDescriptors:
        """Compute the ProtFP descriptors of the peptide.

        See `~peptides.ProtFPDescriptors` for more information.

        Returns:
            `peptides.ProtFPDescriptors`: The computed average of ProtFP
            descriptors of all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, fp in enumerate(peptide.protfp_descriptors()):
            ...     print(f"ProtFP{i+1:<3} {fp: .4f}")
            ProtFP1    0.2065
            ProtFP2   -0.0565
            ProtFP3    1.9930
            ProtFP4   -0.2845
            ProtFP5    0.7315
            ProtFP6    0.7000
            ProtFP7    0.1715
            ProtFP8    0.1135

        """
        out = []
        for i in range(len(tables.PROTFP)):
            p = self.profile(tables.PROTFP[f"ProtFP{i+1}"])
            out.append(_sum(p) / len(self))
        return ProtFPDescriptors(*out)

    def sneath_vectors(self) -> SneathVectors:
        """Compute the Sneath vectors for the peptide.

        See `~peptides.SneathVectors` for more information.

        Returns:
            `peptides.SneathVectors`: The computed average of Sneath vectors
            of all the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, fp in enumerate(peptide.sneath_vectors()):
            ...     print(f"SV{i+1:<3} {fp: .5f}")
            SV1    0.19620
            SV2    0.04655
            SV3    0.04050
            SV4    0.02775

        """
        out = []
        for i in range(len(tables.SNEATH)):
            p = self.profile(tables.SNEATH[f"SV{i+1}"])
            out.append(_sum(p) / len(self))
        return SneathVectors(*out)

    def st_scales(self) -> STScales:
        """Compute the ST-scales of the peptide.

        See `~peptides.STScales` for more information.

        Returns:
            `peptides.STScales`: The computed average of ST-scales of all
            the amino acids in the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, st in enumerate(peptide.st_scales()):
            ...     print(f"ST{i+1:<3} {st: .5f}")
            ST1   -0.63760
            ST2    0.07965
            ST3    0.05150
            ST4    0.07135
            ST5   -0.27905
            ST6   -0.80995
            ST7    0.58020
            ST8    0.54400

        """
        out = []
        for i in range(len(tables.ST_SCALES)):
            p = self.profile(tables.ST_SCALES[f"ST{i+1}"])
            out.append(_sum(p) / len(self))
        return STScales(*out)

    def svger_descriptors(self) -> SVGERDescriptors:
        """Compute the SVGER descriptors of the peptide.

        See `~peptides.SVGERDescriptors` for more information.

        Returns:
            `peptides.SVGERDescriptors`: The computed average of SVGER
            descriptors of all the amino acids in the peptide.

        .. versionadded:: 0.3.2

        """
        out = []
        for i in range(len(tables.SVGER)):
            p = self.profile(tables.SVGER[f"SVGER{i+1}"])
            out.append(_sum(p) / len(self))
        return SVGERDescriptors(*out)

    def t_scales(self) -> TScales:
        """Compute the T-scales of the peptide.

        See `~peptides.TScales` for more information.

        Returns:
            `peptides.TScales`: The computed average of T-scales of all the
            amino acids in the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, t in enumerate(peptide.t_scales()):
            ...     print(f"T{i+1:<3} {t: .4f}")
            T1   -3.2700
            T2   -0.0035
            T3   -0.3855
            T4   -0.1475
            T5    0.7585

        """
        out = []
        for i in range(len(tables.T_SCALES)):
            p = self.profile(tables.T_SCALES[f"T{i+1}"])
            out.append(_sum(p) / len(self))
        return TScales(*out)

    def vhse_scales(self) -> VHSEScales:
        """Compute the VHSE-scales of the peptide.

        See `~peptides.VHSEScales` for more information.

        Returns:
            `peptides.VHSEScales`: The computed average of VHSE-scales of
            the amino acids in the peptide. *VHSE1* and *VHSE2* represent
            hydrophobic properties, *VHSE3* and *VHSE4* represent steric
            properties, while *VHSE5*, *VHSE6*, *VHSE7* and *VHSE8*
            represent electronic properties.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, vhse in enumerate(peptide.vhse_scales()):
            ...     print(f"VHSE{i+1:<3} {vhse: .4f}")
            VHSE1   -0.1150
            VHSE2    0.0630
            VHSE3   -0.0055
            VHSE4    0.7955
            VHSE5    0.4355
            VHSE6    0.2485
            VHSE7    0.1740
            VHSE8   -0.0960

        """
        out = []
        for i in range(len(tables.VHSE)):
            p = self.profile(tables.VHSE[f"VHSE{i+1}"])
            out.append(_sum(p) / len(self))
        return VHSEScales(*out)

    def vstpv_descriptors(self) -> VSTPVDescriptors:
        """Compute the VSTPV descriptors of the peptide.

        See `~peptides.VSTPVDescriptors` for more information.

        Returns:
            `peptides.VSTPVDescriptors`: The computed VSTPV descriptors for
            the peptide.

        """
        out = []
        for i in range(len(tables.VSTPV)):
            p = self.profile(tables.VSTPV[f"VSTPV{i+1}"])
            out.append(_sum(p) / len(self))
        return VSTPVDescriptors(*out)

    def z_scales(self) -> ZScales:
        """Compute the Z-scales of the peptide.

        See `~peptides.ZScales` for more information.

        Returns:
            `peptides.ZScales`: The computed average of Z-scales of all
            the amino acid in the peptide.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, z in enumerate(peptide.z_scales()):
            ...     print(f"Z{i+1:<3} {0.0+round(z,5): .4f}")
            Z1    0.5520
            Z2    0.0985
            Z3    0.0000
            Z4    0.8130
            Z5   -0.8285

        """
        out = []
        for i in range(len(tables.Z_SCALES)):
            p = self.profile(tables.Z_SCALES[f"Z{i+1}"])
            out.append(_sum(p) / len(self))
        return ZScales(*out)

    __DESCRIPTORS = {
        "AF": atchley_factors,
        "BLOSUM": blosum_indices,
        "PP": cruciani_properties,
        "F": fasgai_vectors,
        "KF": kidera_factors,
        "MSWHIM": ms_whim_scores,
        "E": pcp_descriptors,
        "PD": physical_descriptors,
        "PRIN": prin_components,
        "ProtFP": protfp_descriptors,
        "SV": sneath_vectors,
        "ST": st_scales,
        "SVGER": svger_descriptors,
        "T": t_scales,
        "VHSE": vhse_scales,
        "VSTPV": vstpv_descriptors,
        "Z": z_scales,
    }
