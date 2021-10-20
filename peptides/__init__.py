# import numpy

import array
import math

from .data import tables

__version__ = "0.1.0"
__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
__credits__ = """
Daniel Osorio *et al.* for the ``Peptides`` R package.
""".strip()


class Peptide(object):

    def __init__(self, sequence: str):
        self.sequence = sequence

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

    def auto_correlation(self, lag, property, center=True):
        """Compute the auto-correlation index of a peptide sequence.
        """
        raise NotImplementedError("auto_correlation")

    def auto_covariance(self, lag, property, center=True):
        """Compute the auto-covariance index of a peptide sequence.
        """
        raise NotImplementedError("auto_covariance")

    def blosum_indices(self):
        """Compute the BLOSUM62-derived indices of a peptide sequence.
        """
        out = array.array("d")
        for i in range(10):
            scale = tables.BLOSUM[f"BLOSUM{i+1}"]
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def boman(self):
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
            raise ValueError("Invalid pK scale name: {!r}".format(pKscale))

        # nterm
        charge = 1.0 / (1.0 + 10**(1.0 * (pH - scale['nTer'])))
        # aa
        for aa in self.sequence:
            sign = sign_scale.get(aa, 0)
            charge += sign / (1 + 10**(sign * (pH - scale.get(aa, 0))))
        # cterm
        charge += -1.0 / (1.0 + 10**(-1.0 * (pH - scale['cTer'])))

        return charge

    def cross_covariance(self, lag, property1, property2, center=True):
        """Compute the cross-covariance index of a peptide sequence.
        """
        raise NotImplementedError("cross_covariance")

    def cruciani_properties(self):
        """Compute the Cruciani properties of protein sequence.
        """
        out = array.array("d")
        for i in range(3):
            scale = tables.CRUCIANI[f"scale{i+1}"]
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def fasgai_vectors(self):
        """Compute the FASGAI vectors of a protein sequence.
        """
        out = array.array("d")
        for i in range(6):
            scale = tables.FASGAI[f"scale{i+1}"]
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def hydrophobic_moment(self, angle: int = 100, window: int = 11):
        """Compute the maximal hydrophobic moment of a protein sequence.

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.hydrophobic_moment(angle=100)
            0.519922...
            >>> peptide.hydrophobic_moment(angle=160)
            0.270590...

        """
        moment = 0.0
        # load the hydrophobicity scale
        scale = tables.HYDROPHOBICITY["Eisenberg"]
        # create a vector of angles
        angles = [(angle*i)%360 for i in range(window)]

        # compute the moment using a sliding window
        for i in range(len(self.sequence) - window):
            # compute sin and cos of angles
            sumsin = sumcos = 0.0
            for aa, theta in zip(self.sequence[i:i+window], angles):
                sumsin += scale[aa] * math.sin(math.radians(theta))
                sumcos += scale[aa] * math.cos(math.radians(theta))
            # compute hydrophobic moment of window
            hm = math.sqrt(sumsin**2 + sumcos**2) / window
            # record if maximum
            if hm > moment:
                moment = hm

        # return the maximal hydrophobic moment
        return moment

    def hydrophobicity(self, scale: str = "KyteDoolittle"):
        """Compute the hydrophobicity index of a protein sequence.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> peptide.hydrophobicity(scale="Aboderin")
            3.84...
            >>> peptide.hydrophobicity(scale="AbrahamLeo")
            0.092...

        """
        scale = _aadata.HYDROPHOBICITY.get(scale)
        if scale is None:
            raise ValueError(f"Invalid hydrophobicity scale: {scale!r}")
        return sum(scale[aa] for aa in self.sequence) / len(self.sequence)

    def instability_index(self):
        """Compute the instability index of a protein sequence.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> round(peptide.instability_index(), 2)
            83.68

        """
        scale = tables.INSTABILITY["Guruprasad"]
        gp = sum(scale[self.sequence[i:i+2]] for i in range(len(self.sequence) - 1))
        return gp * 10 / (len(self.sequence))

    def kidera_factors(self):
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
        for i in range(10):
            scale = tables.KIDERA[f"KF{i+1}"]
            out.append(sum(scale.get(aa, 0.0) for aa in self.sequence) / len(self.sequence))
        return out

    def mass_shift(self):
        """Compute the mass difference of modified peptides.
        """
        raise NotImplementedError("mass_shift")

    def membrane_position(self):
        """Compute the theoretical class of a protein sequence.
        """
        raise NotImplementedError("membrane_position")

    def ms_whim_scores(self):
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
        for i in range(3):
            scale = tables.MSWHIM[f"MSWHIM{i+1}"]
            out.append(sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence))
        return out

    def molecular_weight(self, monoisotopic=False, average_scale="expasy", label=None, aa_shift=None):
        """Compute the molecular weight of a protein sequence.
        """
        raise NotImplementedError("molecular_weight")

    def mz(self, charge=2, label=None, aa_shift=None, cysteins=57.021464):
        """Compute the m/z (mass/charge) ratio for a peptide.
        """
        raise NotImplementedError("mz")

    def isoelectric_point(self):
        """Compute the isoelectric point of a protein sequence.

        The isoelectric point (*pI*), is the *pH* at which a particular
        molecule or surface carries no net electrical charge.

        """
        raise NotImplementedError("isoelectric_point")

    def protfp_descriptors(self):
        """Compute the protFP descriptors of a protein sequence.
        """
        out = array.array("d")
        for i in range(8):
            scale = tables.PROTFP[f"ProtFP{i+1}"]
            out.append(sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence))
        return out

    def st_scales(self):
        """Compute the ST-scales of a protein sequence.
        """
        out = array.array("d")
        for i in range(8):
            scale = tables.ST_SCALES[f"ST{i+1}"]
            out.append(sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence))
        return out

    def t_scales(self):
        """Compute the T-scales of a protein sequence.
        """
        out = array.array("d")
        for i in range(5):
            scale = tables.T_SCALES[f"T{i+1}"]
            out.append(sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence))
        return out

    def vhse_scales(self):
        """Compute the VHSE-scales of a protein sequence.
        """
        out = array.array("d")
        for i in range(8):
            scale = tables.VHSE[f"VHSE{i+1}"]
            out.append(sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence))
        return out

    def z_scales(self):
        """Compute the Z-scales of a protein sequence.
        """
        out = array.array("d")
        for i in range(8):
            scale = tables.Z_SCALES[f"Z{i+1}"]
            out.append(sum(scale.get(aa, 0) for aa in self.sequence) / len(self.sequence))
        return out
