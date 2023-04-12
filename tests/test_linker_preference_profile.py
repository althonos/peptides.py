import unittest
import urllib.request
from urllib.parse import urlencode

import peptides


class TestLinkerPreferenceProfile(unittest.TestCase):

    DOMCUT_URL = "http://www.bork.embl-heidelberg.de/Docu/mikita/domplot.cgi"

    def assertAsExpected(self, peptide):
        data = urlencode({"seqnam": b"", "sequence": peptide.sequence.encode("ascii"), "outform": b"dat"})
        req = urllib.request.Request(self.DOMCUT_URL, data=data.encode())
        with urllib.request.urlopen(req) as res:
            expected = [float(l.split()[2]) for l in res.readlines() if not l.startswith(b"#")]
        profile = peptide.linker_preference_profile()
        for i, (x, y) in enumerate(zip(profile, expected[7:])):
            self.assertAlmostEqual(x, y, places=2)

    def test_colicin_A(self):
        peptide = peptides.Peptide("MPGFNYGGKGDGTGWSSERGSGPEPGGGSHGNSGGHDRGDSSNVGNESVTVMKPGDSYNTPWGKVIINAAGQPTMNGTVMTADNSSMVPYGRGFTRVLNSLVNNPVSPAGQNGGKSPVQTAVENYLMVQSGNLPPGYWLSNGKVMTEVREERTSGGGGKNGNERTWTVKVPREVPQLTASYNEGMRIRQEAADRARAEANARALAEEEARAIASGKSKAEFDAGKRVEAAQAAINTAQLNVNNLSGAVSAANQVITQKQAEMTPLKNELAAANQRVQETLKFINDPIRSRIHFNMRSGLIRAQHNVDTKQNEINAAVANRDALNSQLSQANNILQNARNEKSAADAALSAATAQRLQAEAALRAAAEAAEKARQRQAEEAERQRQAMEVAEKAKDERELLEKTSELIAGMGDKIGEHLGDKYKAIAKDIADNIKNFQGKTIRSFDDAMASLNKITANPAMKINKADRDALVNAWKHVDAQDMANKLGNLSKAFKVADVVMKVEKVREKSIEGYETGNWGPLMLEVESWVLSGIASSVALGIFSATLGAYALSLGVPAIAVGIAGILLAAVVGALIDDKFADALNNEIIRPAH")
        self.assertAsExpected(peptide)

    def test_trpc(self):
        peptide = peptides.Peptide("MADSGLVDHSPHHPTKAAQLSTASNVILIDNYDSFTWNVYQYLVLEGATVNVFRNDQITLEELIAKKPTQLVISPGPGHPETDAGISSAAIQYFSGKIPIFGVCMGQQCIITCFGGKVDVTGEILHGKTSPLKHDGKGAYEGLPGSLAVTRYHSLAGTHATIPDCLEVSSSVQLADDSNKDVIMGVRHKKLAVEGVQFHPESILTEYGRIMFRNFLKLTAGTWEGNGKHFGEQSSTTKATVPSNPPPKTDKKLSILERIYDHRRAAVAVQKTIPSQRPADLQAAYDLNLAPPQIPFPARLRQSPYPLSLMAEIKRASPSKGMIAENACAPAQARQYAKAGASVISVLTEPEWFKGSIDDLRAVRQSLEGMTNRPAILRKEFVFDEYQILEARLAGADTVLLIVKMLSVELLTRLYHYSRSLGMEPLVEVNTPEEMKIAVDLGAEVIGVNNRDLTSFEVDLGTTSRLMDQVPSSTIVCALSGISGPKDVEAYKKEGVKAILVGEALMRAADTATFIAELLGGSSQTVSSESRRSPLVKICGTRSEEAARAAIEAGADLIGIIMVQGRTGCVPDDVALPISQVVRSTPKPASQALHTSQEPPAATSVEYFDHSAKILRHPSRALLVGVFQNQPLDYILSQQQKLGLDVVQLHGSEPLEWAKLIPVPVIRKFGLDEPAIARRAYHSLPLLDSGVGGTGELLDQSRVQNVLDKDCGLRVILAGGLDPTNVAGIVQKLGESGRKVVGVDVSSGVESDGAQDLNKIRAFVQAVRGL")
        self.assertAsExpected(peptide)