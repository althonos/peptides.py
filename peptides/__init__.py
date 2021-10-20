# import numpy

import array

from . import _aadata

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
        for scale in _aadata.BLOSUM.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    # fmt: off
    _BOMAN_SCALE = dict(
        L = 4.92, I = 4.92, V = 4.04, F = 2.98, M = 2.35, W = 2.33,
        A = 1.81, C = 1.28, G = 0.94, Y = -0.14, T = -2.57, S = -3.40,
        H = -4.66, Q = -5.54, K = -5.55, N = -6.64, E = -6.81, D = -8.72,
        R = -14.92,
    )

    def boman(self):
        """Compute the Boman (potential peptide interaction) index.

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.boman()
            -1.2358...

        """
        return -sum(self._BOMAN_SCALE.get(aa, 0.0) for aa in self.sequence) / len(self.sequence)

    # fmt: off
    _PK_SCALES = {
        "Bjellqvist": dict(
            C = 9.00, D = 4.05, E = 4.45, H = 5.98, K = 10.0, R = 12.0,
            Y = 10.0, c = 3.55, n = 7.50,
        ),
        "Dawson": dict(
            C = 8.3, D = 3.9, E = 4.3, H = 6.0, K = 10.5, R = 12.0, Y = 10.1,
            c = 3.2, n = 8.2,
        ),
        "EMBOSS": dict(
            C = 8.5, D = 3.9, E = 4.1, H = 6.5, K = 10.8, R = 12.5, Y = 10.1,
            c = 3.6, n = 8.6,
        ),
        "Murray": dict(
            C = 8.33, D = 3.68, E = 4.25, H = 6.0, K = 11.50, R = 11.50,
            Y = 10.07, c = 2.15, n = 9.52,
        ),
        "Sillero": dict(
            C = 9, D = 4, E = 4.5, H = 6.4, K = 10.4, R = 12.0, Y = 10.0,
            c = 3.2, n = 8.2,
        ),
        "Solomon": dict(
            C = 8.3, D = 3.9, E = 4.3, H = 6.0, K = 10.5, R = 12.5, Y = 10.1,
            c = 3.2, n = 8.2,
        ),
        "Stryer": dict(
            C = 8.5, D = 4.4, E = 4.4, H = 6.5, K = 10.0, R = 12.0, Y = 10.0,
            c = 3.2, n = 8.2,
        ),
        "Lehninger": dict(
            C = 8.18, D = 3.65, E = 4.25, H = 6.0, K = 10.53, R = 12.48,
            Y = 10.07, c = 2.34, n = 9.69,
        ),
        "Rodwell": dict(
            C = 8.33, D = 3.86, E = 4.25, H = 6, K = 11.50, R = 11.50,
            Y = 10.70, c = 3.10, n = 8,
        )
    }

    # fmt: off
    _SIGN = dict(
        R = 1.0,
        H = 1.0,
        K = 1.0,
        D = -1.0,
        E = -1.0,
        C = -1.0,
        Y = -1.0,
    )

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
        scale = self._PK_SCALES.get(pKscale)
        if scale is None:
            raise ValueError("Invalid pK scale name: {!r}".format(pKscale))

        # nterm
        charge = 1.0 / (1.0 + 10**(1.0 * (pH - scale['n'])))
        # aa
        for aa in self.sequence:
            sign = self._SIGN.get(aa, 0)
            charge += sign / (1 + 10**(sign * (pH - scale.get(aa, 0))))
        # cterm
        charge += -1.0 / (1.0 + 10**(-1.0 * (pH - scale['c'])))

        return charge

    def cross_covariance(self, lag, property1, property2, center=True):
        """Compute the cross-covariance index of a peptide sequence.
        """
        raise NotImplementedError("cross_covariance")

    def cruciani_properties(self):
        """Compute the Cruciani properties of protein sequence.
        """
        out = array.array("d")
        for scale in _aadata.CRUCIANIPROPERTIES.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def fasgai_vectors(self):
        """Compute the FASGAI vectors of a protein sequence.
        """
        out = array.array("d")
        for scale in _aadata.FASGAI.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def hydrophobic_moment(self, angle: int = 100, window: int = 11):
        """Compute the maximal hydrophobic moment of a protein sequence.

        Example:
            >>> peptide = Peptide("FLPVLAGLTPSIVPKLVCLLTKKC")
            >>> peptide.hydrophobic_moment(angle=100)
            0.5199226...
            >>> peptide.hydrophobic_moment(angle=160)
            0.2705906...

        """
        raise NotImplementedError("hydrophobic_moment")

    def hydrophobicity(self, scale: str = "KyteDoolittle"):
        raise NotImplementedError("hydrophobicity")

    # fmt: off
    _GURUPRASAD_SCALE = dict(
        WW = 1, WC = 1, WM = 24.68, WH = 24.68, WY = 1, WF = 1, WQ = 1,
        WN = 13.34, WI = 1, WR = 1, WD = 1, WP = 1, WT = -14.03, WK = 1,
        WE = 1, WV = -7.49, WS = 1, WG = -9.37, WA = -14.03, WL = 13.34,
        CW = 24.68, CC = 1, CM = 33.6, CH = 33.6, CY = 1, CF = 1, CQ = -6.54,
        CN = 1, CI = 1, CR = 1, CD = 20.26, CP = 20.26, CT = 33.6, CK = 1,
        CE = 1, CV = -6.54, CS = 1, CG = 1, CA = 1, CL = 20.26, MW = 1,
        MC = 1, MM = -1.88, MH = 58.28, MY = 24.68, MF = 1, MQ = -6.54,
        MN = 1, MI = 1, MR = -6.54, MD = 1, MP = 44.94, MT = -1.88, MK = 1,
        ME = 1, MV = 1, MS = 44.94, MG = 1, MA = 13.34, ML = 1, HW = -1.88,
        HC = 1, HM = 1, HH = 1, HY = 44.94, HF = -9.37, HQ = 1, HN = 24.68,
        HI = 44.94, HR = 1, HD = 1, HP = -1.88, HT = -6.54, HK = 24.68,
        HE = 1, HV = 1, HS = 1, HG = -9.37, HA = 1, HL = 1, YW = -9.37,
        YC = 1, YM = 44.94, YH = 13.34, YY = 13.34, YF = 1, YQ = 1, YN = 1,
        YI = 1, YR = -15.91, YD = 24.68, YP = 13.34, YT = -7.49, YK = 1,
        YE = -6.54, YV = 1, YS = 1, YG = -7.49, YA = 24.68, YL = 1,
        FW = 1, FC = 1, FM = 1, FH = 1, FY = 33.6, FF = 1, FQ = 1, FN = 1,
        FI = 1, FR = 1, FD = 13.34, FP = 20.26, FT = 1, FK = -14.03, FE = 1,
        FV = 1, FS = 1, FG = 1, FA = 1, FL = 1, QW = 1, QC = -6.54, QM = 1,
        QH = 1, QY = -6.54, QF = -6.54, QQ = 20.26, QN = 1, QI = 1, QR = 1,
        QD = 20.26, QP = 20.26, QT = 1, QK = 1, QE = 20.26, QV = -6.54,
        QS = 44.94, QG = 1, QA = 1, QL = 1, NW = -9.37, NC = -1.88, NM = 1,
        NH = 1, NY = 1, NF = -14.03, NQ = -6.54, NN = 1, NI = 44.94, NR = 1,
        ND = 1, NP = -1.88, NT = -7.49, NK = 24.68, NE = 1, NV = 1, NS = 1,
        NG = -14.03, NL = 1, IW = 1, IC = 1, IM = 1, IH = 13.34, IY = 1,
        IF = 1, IQ = 1, IN = 1, II = 1, IR = 1, ID = 1, IP = -1.88, IT = 1,
        IK = -7.49, IE = 44.94, IV = -7.49, IS = 1, IG = 1, IA = 1, IL = 20.26,
        RW = 58.28, RC = 1, RM = 1, RH = 20.26, RY = -6.54, RF = 1, RQ = 20.26,
        RN = 13.34, RI = 1, RR = 58.28, RD = 1, RP = 20.26, RT = 1, RK = 1,
        RE = 1, RV = 1, RS = 44.94, RG = -7.49, RA = 1, RL = 1, DW = 1, DC = 1,
        DM = 1, DH = 1, DY = 1, DF = -6.54, DQ = 1, DN = 1, DI = 1, DR = -6.54,
        DD = 1, DP = 1, DT = -14.03, DK = -7.49, DE = 1, DV = 1, DS = 20.26,
        DG = 1, DA = 1, DL = 1, PW = -1.88, PC = -6.54, PM = -6.54, PH = 1,
        PY = 1, PF = 20.26, PQ = 20.26, PN = 1, PI = 1, PR = -6.54, PD = -6.54,
        PP = 20.26, PT = 1, PK = 1, PE = 18.38, PV = 20.26, PS = 20.26, PG = 1,
        PA = 20.26, PL = 1, TW = -14.03, TC = 1, TM = 1, TH = 1, TY = 1, TF = 13.34,
        TQ = -6.54, TN = -14.03, TI = 1, TR = 1, TD = 1, TP = 1, TT = 1, TK = 1,
        TE = 20.26, TV = 1, TS = 1, TG = -7.49, TA = 1, TL = 1, KW = 1, KC = 1,
        KM = 33.6, KH = 1, KY = 1, KF = 1, KQ = 24.68, KN = 1, KI = -7.49,
        KR = 33.6, KD = 1, KP = -6.54, KT = 1, KK = 1, KE = 1, KV = -7.49,
        KS = 1, KG = -7.49, KA = 1, KL = -7.49, EW = -14.03, EC = 44.94,
        EM = 1, EH = -6.54, EY = 1, EF = 1, EQ = 20.26, EN = 1, EI = 20.26,
        ER = 1, ED = 20.26, EP = 20.26, ET = 1, EK = 1, EE = 33.6, EV = 1,
        ES = 20.26, EG = 1, EA = 1, EL = 1, VW = 1, VC = 1, VM = 1, VH = 1,
        VY = -6.54, VF = 1, VQ = 1, VN = 1, VI = 1, VR = 1, VD = -14.03,
        VP = 20.26, VT = -7.49, VK = -1.88, VE = 1, VV = 1, VS = 1, VG = -7.49,
        VA = 1, VL = 1, SW = 1, SC = 33.6, SM = 1, SH = 1, SY = 1, SF = 1,
        SQ = 20.26, SN = 1, SI = 1, SR = 20.26, SD = 1, SP = 44.94, ST = 1,
        SK = 1, SE = 20.26, SV = 1, SS = 20.26, SG = 1, SA = 1, SL = 1,
        GW = 13.34, GC = 1, GM = 1, GH = 1, GY = -7.49, GF = 1, GQ = 1,
        GN = -7.49, GI = -7.49, GR = 1, GD = 1, GP = 1, GT = -7.49, GK = -7.49,
        GE = -6.54, GV = 1, GS = 1, GG = 13.34, GA = -7.49, GL = 1, AW = 1,
        AC = 44.94, AM = 1, AH = -7.49, AY = 1, AF = 1, AQ = 1, AN = 1,
        AI = 1, AR = 1, AD = -7.49, AP = 20.26, AT = 1, AK = 1, AE = 1,
        AV = 1, AS = 1, AG = 1, AA = 1, AL = 1, LW = 24.68, LC = 1, LM = 1,
        LH = 1, LY = 1, LF = 1, LQ = 33.6, LN = 1, LI = 1, LR = 20.26, LD = 1,
        LP = 20.26, LT = 1, LK = -7.49, LE = 1, LV = 1, LS = 1, LG = 1, LA = 1,
        LL = 1,
    )

    def instability_index(self):
        """Compute the instability index of a protein sequence.

        Example:
            >>> peptide = Peptide("QWGRRCCGWGPGRRYCVRWC")
            >>> round(peptide.instability_index(), 2)
            83.68

        """
        gp = sum(
            self._GURUPRASAD_SCALE[self.sequence[i:i+2]]
            for i in range(len(self.sequence) - 1)
        )
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
        for scale in _aadata.KIDERAFACTORS.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
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
        for scale in _aadata.MSWHIM.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
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
        for scale in _aadata.PROTFP.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def st_scales(self):
        """Compute the ST-scales of a protein sequence.
        """
        out = array.array("d")
        for scale in _aadata.STSCALES.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def t_scales(self):
        """Compute the T-scales of a protein sequence.
        """
        out = array.array("d")
        for scale in _aadata.TSCALES.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def vhse_scales(self):
        """Compute the VHSE-scales of a protein sequence.
        """
        out = array.array("d")
        for scale in _aadata.VHSE.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out

    def z_scales(self):
        """Compute the Z-scales of a protein sequence.
        """
        # FIXME!
        out = array.array("d")
        for scale in _aadata.ZSCALES.values():
            out.append(sum(scale[aa] for aa in self.sequence) / len(self.sequence))
        return out
