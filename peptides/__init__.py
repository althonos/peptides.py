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

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.boman()
            -1.2358...

        """
        scale = tables.BOMAN["Boman"]
        return -sum(scale.get(aa, 0.0) for aa in self.sequence) / len(self.sequence)

    def charge(self, pH: float = 7, pKscale: str = "Lehninger") -> float:
        """Compute the theoretical net charge of a peptide sequence.

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

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.hydrophobic_moment(angle=100)
            0.519922...
            >>> peptide.hydrophobic_moment(angle=160)
            0.270590...

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

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> peptide.hydrophobicity(scale="Aboderin")
            3.84...
            >>> peptide.hydrophobicity(scale="AbrahamLeo")
            0.092...

        """
        table = tables.HYDROPHOBICITY.get(scale)
        if table is None:
            raise ValueError(f"Invalid hydrophobicity scale: {scale!r}")
        return sum(table[aa] for aa in self.sequence) / len(self.sequence)

    def instability_index(self) -> float:
        """Compute the instability index of a protein sequence.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> round(peptide.instability_index(), 2)
            83.68

        """
        scale = tables.INSTABILITY["Guruprasad"]
        gp = sum(scale[self.sequence[i : i + 2]] for i in range(len(self.sequence) - 1))
        return gp * 10 / (len(self.sequence))

    def isoelectric_point(self, pKscale: str = "EMBOSS") -> float:
        """Compute the isoelectric point of a protein sequence.

        The isoelectric point (*pI*), is the *pH* at which a particular
        molecule or surface carries no net electrical charge.

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

        Example:
            >>> peptide = Peptide("EGVNDNECEGFFSAR")
            >>> peptide.mass_shift(aa_shift="silac_13c")
            6.020129
            >>> peptide.mass_shift(aa_shift=dict(R=10.00827))
            10.00827

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

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> peptide.molecular_weight()
            2485.91...
            >>> peptide.molecular_weight(average="mascot")
            2485.89...
            >>> peptide.molecular_weight(average="monoisotopic")
            2484.11...

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
        """Compute the m/z (mass/charge) ratio for a peptide.

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
