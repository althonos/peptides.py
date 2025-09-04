import unittest
import math
import peptides


class TestSequenceVetting(unittest.TestCase):
    """Test cases for the new sequence vetting functions: entropy, max_frequency, and longest_run."""

    def test_entropy_basic(self):
        """Test basic entropy calculations."""
        # Single amino acid sequence (minimum entropy)
        peptide1 = peptides.Peptide("AAAA")
        self.assertEqual(peptide1.entropy(), 0.0)
        
        # Two amino acids with equal frequency
        peptide2 = peptides.Peptide("AALS")
        expected2 = -0.5 * math.log2(0.5) - 0.25 * math.log2(0.25) - 0.25 * math.log2(0.25)
        self.assertAlmostEqual(peptide2.entropy(), expected2, places=10)
        
        # All 20 standard amino acids (maximum entropy for standard set)
        peptide3 = peptides.Peptide("ACDEFGHIKLMNPQRSTVWY")
        expected3 = math.log2(20)  # Maximum entropy for 20 equally frequent amino acids
        self.assertAlmostEqual(peptide3.entropy(), expected3, places=10)

    def test_entropy_edge_cases(self):
        """Test entropy function with edge cases."""
        # Empty sequence
        peptide_empty = peptides.Peptide("")
        self.assertEqual(peptide_empty.entropy(), 0.0)
        
        # Single amino acid
        peptide_single = peptides.Peptide("A")
        self.assertEqual(peptide_single.entropy(), 0.0)
        
        # Two identical amino acids
        peptide_two = peptides.Peptide("AA")
        self.assertEqual(peptide_two.entropy(), 0.0)

    def test_entropy_diverse_sequences(self):
        """Test entropy with diverse sequences."""
        # Mixed sequence with known frequencies
        peptide = peptides.Peptide("AALLLSSS")
        # A: 2/8 = 0.25, L: 3/8 = 0.375, S: 3/8 = 0.375
        expected = -0.25 * math.log2(0.25) - 0.375 * math.log2(0.375) - 0.375 * math.log2(0.375)
        self.assertAlmostEqual(peptide.entropy(), expected, places=10)
        
        # Sequence with all 26 possible amino acid codes
        peptide_full = peptides.Peptide("ACDEFGHIKLMNPQRSTVWYOUBZXJ")
        # Each amino acid appears once, so frequency = 1/26
        expected_full = math.log2(26)
        self.assertAlmostEqual(peptide_full.entropy(), expected_full, places=10)

    def test_max_frequency_basic(self):
        """Test basic max_frequency calculations."""
        # Single amino acid sequence
        peptide1 = peptides.Peptide("AAAA")
        self.assertEqual(peptide1.max_frequency(), 1.0)
        
        # Two amino acids with equal frequency
        peptide2 = peptides.Peptide("AALS")
        self.assertEqual(peptide2.max_frequency(), 0.5)
        
        # All 20 standard amino acids (equal frequency)
        peptide3 = peptides.Peptide("ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(peptide3.max_frequency(), 0.05)  # 1/20

    def test_max_frequency_edge_cases(self):
        """Test max_frequency function with edge cases."""
        # Empty sequence
        peptide_empty = peptides.Peptide("")
        self.assertEqual(peptide_empty.max_frequency(), 0.0)
        
        # Single amino acid
        peptide_single = peptides.Peptide("A")
        self.assertEqual(peptide_single.max_frequency(), 1.0)
        
        # Two identical amino acids
        peptide_two = peptides.Peptide("AA")
        self.assertEqual(peptide_two.max_frequency(), 1.0)

    def test_max_frequency_mixed_sequences(self):
        """Test max_frequency with mixed sequences."""
        # Mixed sequence with known frequencies
        peptide = peptides.Peptide("AALLLSSS")
        # A: 2/8 = 0.25, L: 3/8 = 0.375, S: 3/8 = 0.375
        self.assertEqual(peptide.max_frequency(), 0.375)
        
        # Complex sequence
        peptide_complex = peptides.Peptide("AAALLLSSSS")
        # A: 3/10 = 0.3, L: 3/10 = 0.3, S: 4/10 = 0.4
        self.assertEqual(peptide_complex.max_frequency(), 0.4)

    def test_longest_run_basic(self):
        """Test basic longest_run calculations."""
        # Single amino acid sequence
        peptide1 = peptides.Peptide("AAAA")
        self.assertEqual(peptide1.longest_run(), 4)
        
        # Two amino acids with equal frequency
        peptide2 = peptides.Peptide("AALS")
        self.assertEqual(peptide2.longest_run(), 2)
        
        # All 20 standard amino acids (no runs)
        peptide3 = peptides.Peptide("ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(peptide3.longest_run(), 1)

    def test_longest_run_edge_cases(self):
        """Test longest_run function with edge cases."""
        # Empty sequence
        peptide_empty = peptides.Peptide("")
        self.assertEqual(peptide_empty.longest_run(), 0)
        
        # Single amino acid
        peptide_single = peptides.Peptide("A")
        self.assertEqual(peptide_single.longest_run(), 1)
        
        # Two identical amino acids
        peptide_two = peptides.Peptide("AA")
        self.assertEqual(peptide_two.longest_run(), 2)

    def test_longest_run_mixed_sequences(self):
        """Test longest_run with mixed sequences."""
        # Mixed sequence with known runs
        peptide = peptides.Peptide("AALLLSSS")
        self.assertEqual(peptide.longest_run(), 3)  # LLL
        
        # Complex sequence with multiple runs
        peptide_complex = peptides.Peptide("AAALLLSSSS")
        self.assertEqual(peptide_complex.longest_run(), 4)  # SSSS
        
        # Sequence with runs at different positions
        peptide_positions = peptides.Peptide("LLAAASSLL")
        self.assertEqual(peptide_positions.longest_run(), 3)  # AAA

    def test_integration_all_functions(self):
        """Test that all three functions work together correctly."""
        peptide = peptides.Peptide("AALLLSSS")
        
        # Get all three metrics
        entropy_val = peptide.entropy()
        max_freq_val = peptide.max_frequency()
        longest_run_val = peptide.longest_run()
        
        # Verify the values make sense together
        self.assertGreaterEqual(entropy_val, 0.0)
        self.assertLessEqual(entropy_val, math.log2(26))  # Maximum possible entropy
        
        self.assertGreaterEqual(max_freq_val, 0.125)  # 1/8 for this sequence
        self.assertLessEqual(max_freq_val, 1.0)
        
        self.assertGreaterEqual(longest_run_val, 1)
        self.assertLessEqual(longest_run_val, len(peptide.sequence))
        
        # Verify specific values for this sequence
        self.assertAlmostEqual(entropy_val, 1.5613, places=3)
        self.assertEqual(max_freq_val, 0.375)
        self.assertEqual(longest_run_val, 3)

    def test_entropy_max_frequency_relationship(self):
        """Test the mathematical relationship between entropy and max_frequency."""
        # Sequences with high max_frequency should have low entropy
        peptide_high_freq = peptides.Peptide("AAAA")
        self.assertEqual(peptide_high_freq.max_frequency(), 1.0)
        self.assertEqual(peptide_high_freq.entropy(), 0.0)
        
        # Sequences with low max_frequency should have high entropy
        peptide_low_freq = peptides.Peptide("ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(peptide_low_freq.max_frequency(), 0.05)
        self.assertAlmostEqual(peptide_low_freq.entropy(), math.log2(20), places=10)

    def test_longest_run_independence(self):
        """Test that longest_run is independent of overall frequencies."""
        # These sequences have the same frequencies but different run patterns
        peptide1 = peptides.Peptide("AALLLSSS")
        peptide2 = peptides.Peptide("ALSLASLS")
        
        # Same frequencies and entropy
        self.assertEqual(peptide1.frequencies(), peptide2.frequencies())
        self.assertAlmostEqual(peptide1.entropy(), peptide2.entropy(), places=10)
        self.assertEqual(peptide1.max_frequency(), peptide2.max_frequency())
        
        # But different longest runs
        self.assertEqual(peptide1.longest_run(), 3)  # LLL
        self.assertEqual(peptide2.longest_run(), 1)  # No runs

    def test_real_protein_sequences(self):
        """Test with real protein sequences."""
        # Real protein sequence
        real_protein = peptides.Peptide("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDFEWHHFYTLIASRLAWNFKLGGK")
        
        # Verify all functions work with long sequences
        entropy_val = real_protein.entropy()
        max_freq_val = real_protein.max_frequency()
        longest_run_val = real_protein.longest_run()
        
        # Check reasonable ranges
        self.assertGreater(entropy_val, 3.0)  # Should be high for diverse protein
        self.assertLess(entropy_val, 4.7)     # Less than maximum possible
        self.assertGreater(max_freq_val, 0.01) # Some amino acids should be common
        self.assertLess(max_freq_val, 0.2)     # But not dominant
        self.assertGreaterEqual(longest_run_val, 1)  # At least 1
        self.assertLessEqual(longest_run_val, 5)     # Reasonable upper bound for real proteins

    def test_ambiguous_amino_acids(self):
        """Test that functions work with ambiguous amino acid codes."""
        # Sequence with ambiguous codes (B, Z, X, J, O, U)
        peptide_ambig = peptides.Peptide("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        
        # All functions should work without errors
        entropy_val = peptide_ambig.entropy()
        max_freq_val = peptide_ambig.max_frequency()
        longest_run_val = peptide_ambig.longest_run()
        
        # Verify reasonable values
        self.assertGreater(entropy_val, 0.0)
        self.assertLessEqual(entropy_val, math.log2(26))
        self.assertEqual(max_freq_val, 1.0 / 26)  # Each amino acid appears once
        self.assertEqual(longest_run_val, 1)       # No consecutive identical amino acids

    def test_detect_outlier(self):
        """Test the detect_outlier method."""
        # Test with a normal sequence
        normal_peptide = peptides.Peptide("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIGLNKALELVKQKLKELN")
        
        # This should work if distribution data exists
        try:
            result = normal_peptide.detect_outlier()
            # Verify the result structure
            self.assertIsInstance(result.is_outlier, bool)
            self.assertIsInstance(result.issues, list)
            self.assertIsInstance(result.metrics, dict)
            
            # Verify metrics are present
            expected_metrics = ['entropy', 'max_frequency', 'longest_run', 'sequence_length']
            for metric in expected_metrics:
                self.assertIn(metric, result.metrics)
                
        except FileNotFoundError:
            # Skip test if distribution data not available
            self.skipTest("SwissProt distribution data not available. Run generate_swissprot_distributions_for_vetting_functions.py first.")
        
        # Test with a clearly problematic sequence
        problematic_peptide = peptides.Peptide("AAAA")
        try:
            result = problematic_peptide.detect_outlier()
            # This should definitely be an outlier
            self.assertTrue(result.is_outlier)
            self.assertGreater(len(result.issues), 0)
            
        except FileNotFoundError:
            self.skipTest("SwissProt distribution data not available.")


if __name__ == "__main__":
    unittest.main()
