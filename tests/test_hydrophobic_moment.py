import unittest

import peptides


class TestHydrophobicMoment(unittest.TestCase):

    def test_small_peptides(self):
        # https://github.com/althonos/peptides.py/issues/1
        p1 = peptides.Peptide("MLK")
        self.assertAlmostEqual(p1.hydrophobic_moment(window=5, angle=100), 0.8099386)
        p2 = peptides.Peptide("AACQ")
        self.assertAlmostEqual(p2.hydrophobic_moment(window=5, angle=100), 0.3152961)
        p3 = peptides.Peptide("FGGIQ")
        self.assertAlmostEqual(p3.hydrophobic_moment(window=5, angle=100), 0.3184719)
