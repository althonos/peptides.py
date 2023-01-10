import unittest

import peptides


class TestProfile(unittest.TestCase):

    def test_large_window(self):
        p1 = peptides.Peptide("PKLV")
        self.assertEqual(p1.profile(peptides.tables.CHARGE['sign'], window=5), [])

    def test_empty(self):
        p1 = peptides.Peptide("")
        self.assertEqual(p1.profile(peptides.tables.CHARGE['sign']), [])
