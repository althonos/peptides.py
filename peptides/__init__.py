import array
import math
import statistics
import typing

from .data import tables

__all__ = ["Peptide", "tables"]

__version__ = "0.1.0"
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
__credits__ = """
Daniel Osorio, Paola Rondón-Villarreal and Rodrigo Torres for ``Peptides``.
Alan Bleasby for the ``hmoment`` binary of the EMBOSS.
""".strip()


class BLOSUMIndices(typing.NamedTuple):
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
    pp1: float
    pp2: float
    pp3: float


class FasgaiVectors(typing.NamedTuple):
    f1: float
    f2: float
    f3: float
    f4: float
    f5: float
    f6: float


class KideraFactors(typing.NamedTuple):
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


class MSWHIMScores(typing.NamedTuple):
    mswhim1: float
    mswhim2: float
    mswhim3: float


class ProtFPDescriptors(typing.NamedTuple):
    protfp1: float
    protfp2: float
    protfp3: float
    protfp4: float
    protfp5: float
    protfp6: float
    protfp7: float
    protfp8: float


class STScales(typing.NamedTuple):
    st1: float
    st2: float
    st3: float
    st4: float
    st5: float
    st6: float
    st7: float
    st8: float


class TScales(typing.NamedTuple):
    t1: float
    t2: float
    t3: float
    t4: float
    t5: float


class VHSEScales(typing.NamedTuple):
    vhse1: float
    vhse2: float
    vhse3: float
    vhse4: float
    vhse5: float
    vhse6: float
    vhse7: float
    vhse8: float


class ZScales(typing.NamedTuple):
    z1: float
    z2: float
    z3: float
    z4: float
    z5: float


class Peptide(object):

    # --- Magic methods ------------------------------------------------------

    def __init__(self, sequence: str) -> None:
        self.sequence: str = sequence

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
            ['BLOSUM1', ..., 'F1', ..., 'KF1', ..., 'MSWHIM1', ..., 'PP1', ...]

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
        # compute using Cruciani formula
        s1 = s2 = 0.0
        for aa1, aa2 in zip(self.sequence, self.sequence[lag:]):
            s1 += table.get(aa1, 0.0) * table.get(aa2, 0.0)
            s2 += table.get(aa1, 0.0) ** 2
        # return correlation
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
        # compute using Cruciani formula
        s = 0.0
        for aa1, aa2 in zip(self.sequence, self.sequence[lag:]):
            s += table.get(aa1, 0.0) * table.get(aa2, 0.0)
        # return correlation
        return s / len(self.sequence)

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
        # center the table if requested
        if center:
            mu1 = statistics.mean(table1.values())
            sigma1 = statistics.stdev(table1.values())
            table1 = {k: (v - mu1) / sigma1 for k, v in table1.items()}
            mu2 = statistics.mean(table2.values())
            sigma2 = statistics.stdev(table2.values())
            table2 = {k: (v - mu2) / sigma2 for k, v in table2.items()}
        # compute using Cruciani formula
        s = 0.0
        for aa1, aa2 in zip(self.sequence, self.sequence[lag:]):
            s += table1.get(aa1, 0.0) * table2.get(aa2, 0.0)
        # return correlation
        return s / len(self.sequence)

    # --- Physico-chemical properties ----------------------------------------

    def aliphatic_index(self) -> float:
        """Compute the aliphatic index of the peptide.

        Example:
            >>> peptide = Peptide("SDKEVDEVDAALSDLEITLE")
            >>> peptide.aliphatic_index()
            117.0

        """
        ala = self.sequence.count("A") / len(self.sequence)
        val = self.sequence.count("V") / len(self.sequence)
        leu = self.sequence.count("L") / len(self.sequence)
        ile = self.sequence.count("I") / len(self.sequence)
        return (ala + 2.9 * val + 3.9 * (leu + ile)) * 100

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
              doi:10.1046/j.1365-2796.2003.01228.x. PMID:12930229.

        """
        scale = tables.BOMAN["Boman"]
        return -sum(scale.get(aa, 0.0) for aa in self.sequence) / len(self.sequence)

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
              doi:10.1002/elps.11501401163. PMID:8125050.
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
              doi:10.1016/0003-2697(82)90611-x. PMID:7072964.
            - Sillero, A., and A. Maldonado.
              *Isoelectric Point Determination of Proteins and Other
              Macromolecules: Oscillating Method*.
              Computers in Biology and Medicine. 2006 Feb;36(2): 157–66.
              doi:10.1016/j.compbiomed.2004.09.006. PMID:16389075.
            - Solomons, T. W. G.
              *Fundamentals of Organic Chemistry*. New York: Wiley. 1997;5.
              ISBN:978-0-471-28298-3.
            - Stryer, L., J. Augustyniak, and J. Michejda.
              *Biochemia*. Warszawa: Wydawnictwo Naukowe PWN. 2000.
              ISBN:978-83-01-12044-3.

        """
        sign_scale = tables.CHARGE["sign"]
        scale = tables.PK.get(pKscale)
        if scale is None:
            raise ValueError(f"Invalid pK scale: {scale!r}")

        # nterm
        charge = 1.0 / (1.0 + 10 ** (1.0 * (pH - scale["nTer"])))
        # aa
        for aa in self.sequence:
            sign = sign_scale.get(aa, 0)
            charge += sign / (1 + 10 ** (sign * (pH - scale.get(aa, 0))))
        # cterm
        charge += -1.0 / (1.0 + 10 ** (-1.0 * (pH - scale["cTer"])))

        return charge

    def hydrophobic_moment(self, angle: int = 100, window: int = 11) -> float:
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

        References:
            - Eisenberg, D., R. M. Weiss, and T. C. Terwilliger.
              *The Hydrophobic Moment Detects Periodicity in Protein
              Hydrophobicity*. Proceedings of the National Academy of
              Sciences of the United States of America. 1984 Jan;81(1):140–44.
              doi:10.1073/pnas.81.1.140. PMID:6582470.

        """
        scale = tables.HYDROPHOBICITY["Eisenberg"]
        angles = [(angle * i) % 360 for i in range(window)]
        moment = 0.0

        for i in range(len(self.sequence) - window + 1):
            # compute sin and cos of angles
            sumsin = sumcos = 0.0
            for aa, theta in zip(self.sequence[i : i + window], angles):
                sumsin += scale[aa] * math.sin(math.radians(theta))
                sumcos += scale[aa] * math.cos(math.radians(theta))
            # compute hydrophobic moment of window
            hm = math.sqrt(sumsin ** 2 + sumcos ** 2) / window
            if hm > moment:
                moment = hm

        return moment

    def hydrophobicity(self, scale: str = "KyteDoolittle") -> float:
        """Compute the hydrophobicity index of a protein sequence.

        This function calculates the hydrophobicity index of an amino
        acid sequence by averaging the hydrophobicity values of each residue
        using one of the 38 scales from different sources.

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
              doi:10.1016/0020-711X(71)90023-1.
            - Abraham, D.J., and A. J. Leo.
              *Extension of the Fragment Method to Calculate Amino Acid
              Zwitterion and Side Chain Partition Coefficients*.
              Proteins: Structure, Function, and Genetics. 1987;2(2):130–52.
              doi:10.1002/prot.340020207.
            - Argos, P., J. K. Rao, and P. A. Hargrave.
              *Structural Prediction of Membrane-Bound Proteins*.
              European Journal of Biochemistry. Nov 1982;128(2–3):565–75.
              doi:10.1111/j.1432-1033.1982.tb07002.x. PMID:7151796.
            - Black, S. D., and D. R. Mould.
              *Development of Hydrophobicity Parameters to Analyze Proteins
              Which Bear Post- or Cotranslational Modifications*.
              Analytical Biochemistry. Feb 1991;193(1):72-82.
              doi:10.1016/0003-2697(91)90045-u. PMID:2042744.
            - Bull, H. B., and K. Breese.
              *Surface Tension of Amino Acid Solutions: A Hydrophobicity
              Scale of the Amino Acid Residues*. Archives of Biochemistry
              and Biophysics. Apr 1974;161(2):665–70.
              doi:10.1016/0003-9861(74)90352-x. PMID: 4839053
            - Casari, G., and M. J. Sippl.
              *Structure-Derived Hydrophobic Potential. Hydrophobic
              Potential Derived from X-Ray Structures of Globular Proteins
              Is Able to Identify Native Folds*.
              Journal of Molecular Biology. Apr 1992;224(3):725–32.
              doi:10.1016/0022-2836(92)90556-y. PMID:1569551.
            - Chothia, C.
              *The Nature of the Accessible and Buried Surfaces in Proteins*.
              Journal of Molecular Biology. Jul 2917;105(1):1–12.
              doi:10.1016/0022-2836(76)90191-1. PMID:994183.
            - Cid, H., M. Bunster, M. Canales, and F. Gazitúa.
              *Hydrophobicity and Structural Classes in Proteins*.
              Protein Engineering. Jul 1992;5(5):373–75.
              doi:10.1093/protein/5.5.373. PMID:1518784.
            - Cowan, R., and R. G. Whittaker.
              *Hydrophobicity Indices for Amino Acid Residues as Determined
              by High-Performance Liquid Chromatography*.
              Peptide Research. Apr 1990;3(2):75–80.
              PMID:2134053.
            - Eisenberg, D., E. Schwarz, M. Komaromy, and R. Wall.
              *Analysis of Membrane and Surface Protein Sequences with
              the Hydrophobic Moment Plot*.
              Journal of Molecular Biology. Oct 1984;179(1):125–42.
              doi:10.1016/0022-2836(84)90309-7. PMID:6502707.
            - Engelman, D. M., T. A. Steitz, and A. Goldman.
              *Identifying Nonpolar Transbilayer Helices in Amino Acid
              Sequences of Membrane Proteins*. Annual Review of Biophysics
              and Biophysical Chemistry. 1986;15:321–53.
              doi:10.1146/annurev.bb.15.060186.001541. PMID:3521657.
            - Fasman, G. D.
              *Prediction of Protein Structure and the Principles of Protein
              Conformation*. Springer US. 1989.
              doi:10.1007/978-1-4613-1571-1. ISBN:978-0-306-43131-9.
            - Fauchère, J-L., and Pliska V.
              *Hydrophobic Parameters π of Amino-Acid Side Chains from the
              Partitioning of N-Acetyl-Amino-Acid Amides*. European Journal
              of Medicinal Chemistry. 1983;18(4):369–75.
            - Goldsack, D. E., and R. C. Chalifoux.
              *Contribution of the Free Energy of Mixing of Hydrophobic Side
              Chains to the Stability of the Tertiary Structure of Proteins*.
              Journal of Theoretical Biology. Jun 1973;39(3):645–51.
              doi:10.1016/0022-5193(73)90075-1. PMID:4354159.
            - Guy, H. R.
              *Amino Acid Side-Chain Partition Energies and Distribution
              of Residues in Soluble Proteins*.
              Biophysical Journal. Jan 1985;47(1):61–70.
              doi:10.1016/S0006-3495(85)83877-7. PMID:3978191.
            - Hopp, T. P., and K. R. Woods.
              *Prediction of Protein Antigenic Determinants from Amino Acid
              Sequences*. Proceedings of the National Academy of Sciences
              of the United States of America. Jun 1981;78(6):3824–28.
              doi:10.1073/pnas.78.6.3824. PMID:6167991.
            - Janin, J.
              *Surface and inside Volumes in Globular Proteins*.
              Nature. Feb 1979;277(5696):491–92.
              doi:10.1038/277491a0. PMID:763335.
            - Jones, D. D.
              *Amino Acid Properties and Side-Chain Orientation in Proteins:
              A Cross Correlation Approach*.
              Journal of Theoretical Biology. Mar 1975;50(1):167–83.
              doi:10.1016/0022-5193(75)90031-4. PMID:1127956.
            - Juretić, D., D. Zucić, B. Lucić, and N. Trinajstić.
              *Preference Functions for Prediction of Membrane-Buried
              Helices in Integral Membrane Proteins*.
              Computers & Chemistry. Jun 1998;22(4):279–94.
              doi:10.1016/s0097-8485(97)00070-3. PMID:9680689.
            - Kawashima, S., H. Ogata, and M. Kanehisa.
              *AAindex: Amino Acid Index Database*.
              Nucleic Acids Research. Jan 1999;27(1):368–69.
              doi:10.1093/nar/27.1.368. PMID:9847231.
            - Kawashima, S., and M. Kanehisa.
              *AAindex: Amino Acid Index Database*.
              Nucleic Acids Research. Jan 2000;28(1):374.
              doi:10.1093/nar/28.1.374. PMID:10592278.
            - Kawashima, S., P. Pokarowski, M. Pokarowska, A. Kolinski,
              T. Katayama, and M. Kanehisa.
              *AAindex: Amino Acid Index Database, Progress Report 2008*.
              Nucleic Acids Research. Jan 2008;36:D202-205.
              doi:10.1093/nar/gkm998. PMID:17998252.
            - Kidera, A., Y. Konishi, M. Oka, T. Ooi, and H. A. Scheraga.
              *Statistical Analysis of the Physical Properties of the
              20 Naturally Occurring Amino Acids*.
              Journal of Protein Chemistry. Feb 1985;4(1):23-55.
              doi:10.1007/BF01025492.
            - Kuhn, L. A., C. A. Swanson, M. E. Pique, J. A. Tainer, and
              E. D. Getzoff.
              *Atomic and Residue Hydrophilicity in the Context of Folded
              Protein Structures*. Proteins. Dec 1995;23(4):536–47.
              doi:10.1002/prot.340230408. PMID:8749849.
            - Kyte, J., and R. F. Doolittle.
              *A Simple Method for Displaying the Hydropathic Character of a
              Protein*. Journal of Molecular Biology. May 1982;157(1):105–32.
              doi:10.1016/0022-2836(82)90515-0. PMID:7108955.
            - Levitt, M.
              *A Simplified Representation of Protein Conformations for
              Rapid Simulation of Protein Folding*.
              Journal of Molecular Biology. Jun 1976;104(1):59–107.
              doi:10.1016/0022-2836(76)90004-8. PMID:957439.
            - Manavalan, P., and P. K. Ponnuswamy.
              *Hydrophobic Character of Amino Acid Residues in Globular
              Proteins*. Nature. Oct 1978;275(5681):673–74.
              doi:10.1038/275673a0. PMID:703834.
            - Miyazawa, S., and R. L. Jernigan.
              *Estimation of Effective Interresidue Contact Energies from
              Protein Crystal Structures: Quasi-Chemical Approximation*.
              Macromolecules. Mar 1985;18(3):534–52.
              doi:10.1021/ma00145a039.
            - Nakai, K., A. Kidera, and M. Kanehisa.
              *Cluster Analysis of Amino Acid Indices for Prediction of
              Protein Structure and Function*.
              Protein Engineering. Jul 1988;2(2):93–100.
              doi:10.1093/protein/2.2.93. PMID:3244698.
            - Nozaki, Y., and C. Tanford.
              *The Solubility of Amino Acids and Two Glycine Peptides in
              Aqueous Ethanol and Dioxane Solutions. Establishment of a
              Hydrophobicity Scale*.
              The Journal of Biological Chemistry. Apr 1971;246(7):2211–17.
              PMID:5555568.
            - Parker, J. M., D. Guo, and R. S. Hodges.
              *New Hydrophilicity Scale Derived from High-Performance Liquid
              Chromatography Peptide Retention Data: Correlation of
              Predicted Surface Residues with Antigenicity and X-Ray-Derived
              Accessible Sites*. Biochemistry. 1986;25(19):5425–32.
              doi:10.1021/bi00367a013. PMID:2430611.
            - Ponnuswamy, P. K.
              *Hydrophobic Characteristics of Folded Proteins*. Progress in
              Biophysics and Molecular Biology. 1993;59(1):57–103.
              doi:10.1016/0079-6107(93)90007-7. PMID:8419986.
            - Prabhakaran, M.
              *The Distribution of Physical, Chemical and Conformational
              Properties in Signal and Nascent Peptides*.
              The Biochemical Journal. Aug 1990;269(3):691–96.
              doi:10.1042/bj2690691. PMID:2390062.
            - Rao, J. K. M., and P. Argos.
              *A Conformational Preference Parameter to Predict Helices in
              Integral Membrane Proteins*.
              Biochimica Et Biophysica Acta. Jan 1986;869(2):197–214.
              doi:10.1016/0167-4838(86)90295-5. PMID:2935194.
            - Rose, G. D., A. R. Geselowitz, G. J. Lesser, R. H. Lee, and
              M. H. Zehfus.
              *Hydrophobicity of Amino Acid Residues in Globular Proteins*.
              Science (New York, N.Y.). Aug 1985;229(4716):834–38.
              doi:10.1126/science.4023714. PMID:4023714.
            - Roseman, M. A.
              *Hydrophilicity of Polar Amino Acid Side-Chains Is Markedly
              Reduced by Flanking Peptide Bonds*.
              Journal of Molecular Biology. Apr 1988;200(3):513–22.
              doi:10.1016/0022-2836(88)90540-2. PMID:3398047.
            - Sweet, R. M., and D. Eisenberg.
              *Correlation of Sequence Hydrophobicities Measures Similarity
              in Three-Dimensional Protein Structure*.
              Journal of Molecular Biology. Dec 1983;171(4):479-88.
              doi:10.1016/0022-2836(83)90041-4. PMID:6663622.
            - Tomii, K., and M. Kanehisa.
              *Analysis of Amino Acid Indices and Mutation Matrices for
              Sequence Comparison and Structure Prediction of Proteins*.
              Protein Engineering. Jan 1996;9(1):27–36.
              doi:10.1093/protein/9.1.27. PMID:9053899.
            - Welling, G. W., W. J. Weijer, R. van der Zee R, and
              S. Welling-Wester.
              *Prediction of Sequential Antigenic Regions in Proteins*.
              FEBS Letters. Feb 1985;188(2):215-8.
              doi:10.1016/0014-5793(85)80374-4. PMID:2411595.
            - White, S. H., and W. C. Wimley.
              *Membrane Protein Folding and Stability: Physical Principles*.
              Annual Review of Biophysics and Biomolecular Structure.
              1999;28:319–65.
              doi:10.1146/annurev.biophys.28.1.319. PMID:10410805
            - White, S. H., and W. C. Wimley.
              *Hydrophobic Interactions of Peptides with Membrane Interfaces*.
              Biochimica Et Biophysica Acta. Nov 1998;1376(3):339-52.
              doi:10.1016/s0304-4157(98)00021-5. PMID:9804985.
            - Wilson, K. J., A. Honegger, R. P. Stötzel, and G. J. Hughes.
              *The Behaviour of Peptides on Reverse-Phase Supports during
              High-Pressure Liquid Chromatography*.
              The Biochemical Journal. Oct 1981;199(1):31-41.
              doi:10.1042/bj1990031. PMID:7337711.
            - Wimley, W. C., and S. H. White.
              *Experimentally Determined Hydrophobicity Scale for Proteins
              at Membrane Interfaces*.
              Nature Structural Biology. Oct 1996;3(10):842–48.
              doi:10.1038/nsb1096-842. PMID:8836100.
            - Wimley, W. C., T. P. Creamer, and S. H. White.
              *Solvation Energies of Amino Acid Side Chains and Backbone in
              a Family of Host-Guest Pentapeptides*.
              Biochemistry. Apr 1996;35(16):5109–24.
              doi:10.1021/bi9600153. PMID:8611495.
            - Wolfenden, R., L. Andersson, P. M. Cullis, and C. C. Southgate.
              *Affinities of Amino Acid Side Chains for Solvent Water*.
              Biochemistry. Feb 1981;20(4):849–55.
              doi:10.1021/bi00507a030. PMID:7213619.
            - Zimmerman, J. M., N. Eliezer, and R. Simha.
              *The Characterization of Amino Acid Sequences in Proteins by
              Statistical Methods*.
              Journal of Theoretical Biology. Nov 1968;21(2):170–201.
              doi:10.1016/0022-5193(68)90069-6. PMID:5700434.

        """
        table = tables.HYDROPHOBICITY.get(scale)
        if table is None:
            raise ValueError(f"Invalid hydrophobicity scale: {scale!r}")
        return sum(table[aa] for aa in self.sequence) / len(self.sequence)

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
              doi:10.1093/protein/4.2.155. PMID:2075190.

        """
        scale = tables.INSTABILITY["Guruprasad"]
        gp = sum(scale[self.sequence[i : i + 2]] for i in range(len(self.sequence) - 1))
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
            6.020129
            >>> peptide.mass_shift(aa_shift=dict(R=10.00827))
            10.00827

        References:
            - Ong, S-E., I. Kratchmarova, and M. Mann.
              *Properties of 13C-Substituted Arginine in Stable Isotope
              Labeling by Amino Acids in Cell Culture (SILAC)*.
              Journal of Proteome Research. Apr 2003;2(2):173–81.
              doi:10.1021/pr0255708. PMID:12716131.
            - Picotti, P., B. Bodenmiller, L. N. Mueller, B. Domon,
              and R. Aebersold.
              *Full Dynamic Range Proteome Analysis of S. Cerevisiae by
              Targeted Proteomics*. Cell. Aug 2009;138(4):795–806.
              doi:10.1016/j.cell.2009.05.051. PMID:19664813.

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

        s = scale.get("nTer", 0.0) + scale.get("cTer", 0.0)
        s += sum(scale.get(aa, 0.0) for aa in self.sequence)
        return s

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
              doi:10.1385/1-59259-584-7:531. PMID:10027275

        """
        scale = tables.MOLECULAR_WEIGHT.get(average)
        if scale is None:
            raise ValueError(f"Invalid average weight scale: {average!r}")

        # sum the weight of each amino acid
        mass = sum(scale.get(aa, 0.0) for aa in self.sequence)
        # add weight of water molecules
        mass += scale["H2O"]
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

        profile = array.array("d")
        for i in range(len(self.sequence) - window + 1):
            profile.append(self[i : i + window].hydrophobicity(scale=scale))

        return profile

    def hydrophobic_moment_profile(
        self, window: int = 11, angle: int = 100
    ) -> typing.Sequence[float]:
        """Build a hydrophobic moment profile of a sliding window.

        Example:
            >>> peptide = Peptide("ARQQNLFINFCLILIFLLLI")
            >>> uH = peptide.hydrophobic_moment_profile(window=12, angle=100)
            >>> [round(x, 3) for x in uH]
            [0.353, 0.317, 0.274, 0.274, 0.253, 0.113, 0.113, 0.108, 0.132]

        """
        profile = array.array("d")
        for i in range(len(self.sequence) - window + 1):
            profile.append(
                self[i : i + window].hydrophobic_moment(window=window - 1, angle=angle)
            )

        return profile

    def membrane_position_profile(
        self, window: int = 11, angle: int = 100
    ) -> typing.Sequence[str]:
        """Compute the theoretical class of a protein sequence.

        Example:
            >>> peptide = Peptide("ARQQNLFINFCLILIFLLLI")
            >>> peptide.membrane_position_profile(window=12, angle=100)
            ['G', 'G', 'G', 'T', 'S', 'T', 'T', 'T', 'T']
            >>> peptide.membrane_position_profile(window=12, angle=160)
            ['G', 'G', 'G', 'S', 'S', 'S', 'S', 'S', 'S']

        """
        profile_H = self.hydrophobicity_profile(window=window, scale="Eisenberg")
        profile_uH = self.hydrophobic_moment_profile(window=window, angle=angle)

        profile = []
        for h, uh in zip(profile_H, profile_uH):
            m = h * -0.421 + 0.579
            if uh <= m and h >= 0.5:
                profile.append("T")
            elif uh <= m and h <= 0.5:
                profile.append("G")
            elif uh >= m:
                profile.append("S")
            else:
                profile.append("?")

        return profile

    # --- Descriptors --------------------------------------------------------

    def blosum_indices(self) -> BLOSUMIndices:
        """Compute the BLOSUM62-derived indices of a peptide sequence.

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
        out = array.array("d")
        for i in range(len(tables.BLOSUM)):
            scale = tables.BLOSUM[f"BLOSUM{i+1}"]
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return BLOSUMIndices(*out)

    def cruciani_properties(self) -> CrucianiProperties:
        """Compute the Cruciani properties of protein sequence.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, b in enumerate(peptide.cruciani_properties()):
            ...     print(f"PP{i+1:<3} {b: .4f}")
            PP1   -0.1130
            PP2   -0.0220
            PP3    0.2735

        """
        out = array.array("d")
        for i in range(len(tables.CRUCIANI)):
            scale = tables.CRUCIANI[f"PP{i+1}"]
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return CrucianiProperties(*out)

    def fasgai_vectors(self) -> FasgaiVectors:
        """Compute the FASGAI vectors of a protein sequence.

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
        out = array.array("d")
        for i in range(len(tables.FASGAI)):
            scale = tables.FASGAI[f"F{i+1}"]
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return FasgaiVectors(*out)

    def kidera_factors(self) -> KideraFactors:
        """Compute the Kidera factors of a protein sequence.

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
        out = array.array("d")
        for i in range(len(tables.KIDERA)):
            scale = tables.KIDERA[f"KF{i+1}"]
            out.append(
                sum(scale.get(aa, 0.0) for aa in self.sequence) / len(self.sequence)
            )
        return KideraFactors(*out)

    def ms_whim_scores(self) -> MSWHIMScores:
        """Compute the MS-WHIM scores of a protein sequence.

        Example:
            >>> peptide = Peptide("KLKLLLLLKLK")
            >>> for i, mw in enumerate(peptide.ms_whim_scores()):
            ...     print(f"MSWHIM{i+1:<3} {mw: .4f}")
            MSWHIM1   -0.6564
            MSWHIM2    0.4873
            MSWHIM3    0.1164

        """
        out = array.array("d")
        for i in range(len(tables.MSWHIM)):
            scale = tables.MSWHIM[f"MSWHIM{i+1}"]
            out.append(
                sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence)
            )
        return MSWHIMScores(*out)

    def protfp_descriptors(self) -> ProtFPDescriptors:
        """Compute the protFP descriptors of a protein sequence.

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
        out = array.array("d")
        for i in range(len(tables.PROTFP)):
            scale = tables.PROTFP[f"ProtFP{i+1}"]
            out.append(
                sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence)
            )
        return ProtFPDescriptors(*out)

    def st_scales(self) -> STScales:
        """Compute the ST-scales of a protein sequence.

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
        out = array.array("d")
        for i in range(len(tables.ST_SCALES)):
            scale = tables.ST_SCALES[f"ST{i+1}"]
            out.append(
                sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence)
            )
        return STScales(*out)

    def t_scales(self) -> TScales:
        """Compute the T-scales of a protein sequence.

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
        out = array.array("d")
        for i in range(len(tables.T_SCALES)):
            scale = tables.T_SCALES[f"T{i+1}"]
            out.append(
                sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence)
            )
        return TScales(*out)

    def vhse_scales(self) -> VHSEScales:
        """Compute the VHSE-scales of a protein sequence.

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
        out = array.array("d")
        for i in range(len(tables.VHSE)):
            scale = tables.VHSE[f"VHSE{i+1}"]
            out.append(
                sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence)
            )
        return VHSEScales(*out)

    def z_scales(self) -> ZScales:
        """Compute the Z-scales of a protein sequence.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> for i, z in enumerate(peptide.z_scales()):
            ...     print(f"Z{i+1:<3} {z: .4f}")
            Z1    0.5520
            Z2    0.0985
            Z3    0.0000
            Z4    0.8130
            Z5   -0.8285

        """
        out = array.array("d")
        for i in range(len(tables.Z_SCALES)):
            scale = tables.Z_SCALES[f"Z{i+1}"]
            out.append(
                sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence)
            )
        return ZScales(*out)

    __DESCRIPTORS = {
        "BLOSUM": blosum_indices,
        "PP": cruciani_properties,
        "F": fasgai_vectors,
        "KF": kidera_factors,
        "MSWHIM": ms_whim_scores,
        "ProtFP": protfp_descriptors,
        "ST": st_scales,
        "T": t_scales,
        "VHSE": vhse_scales,
        "Z": z_scales,
    }
