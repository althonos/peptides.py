Sequence Vetting Functions
==========================

.. currentmodule:: peptides

The sequence vetting functions provide tools for analyzing the complexity of protein sequences. These functions are particularly useful for quality control.

Entropy
-------

.. autofunction:: peptides.Peptide.entropy

The Shannon entropy measures the diversity of amino acids in a peptide sequence. It is calculated using the formula:

.. math::

    H = -\\sum_{i=1}^{n} p_i \\log_2(p_i)

where :math:`p_i` is the frequency of amino acid :math:`i` and :math:`n` is the number of possible amino acids (26, including ambiguous codes).

**Examples:**

.. code-block:: python

    from peptides import Peptide
    
    # Single amino acid sequence (minimum entropy)
    peptide = Peptide("AAAA")
    print(peptide.entropy())  # Output: 0.0
    
    # Diverse sequence (high entropy)
    peptide = Peptide("ACDEFGHIKLMNPQRSTVWY")
    print(peptide.entropy())  # Output: ~4.32
    
    # Mixed sequence
    peptide = Peptide("AALS")
    print(peptide.entropy())  # Output: 1.5

**Interpretation:** Most values fall between 3.71 and 4.18 bits for SwissProt data.
**Range:** 0.0 to log₂(26) ≈ 4.70 bits

Max Frequency
-------------

.. autofunction:: peptides.Peptide.max_frequency

The maximum frequency identifies the amino acid that appears most frequently in the sequence and returns its frequency. This metric is useful for identifying dominant amino acids and assessing sequence diversity.

**Examples:**

.. code-block:: python

    from peptides import Peptide
    
    # Single amino acid sequence
    peptide = Peptide("AAAA")
    print(peptide.max_frequency())  # Output: 1.0
    
    # Two amino acids with equal frequency
    peptide = Peptide("AALS")
    print(peptide.max_frequency())  # Output: 0.5
    
    # Diverse sequence
    peptide = Peptide("ACDEFGHIKLMNPQRSTVWY")
    print(peptide.max_frequency())  # Output: 0.05

**Interpretation:** Most values fall between 0.085 and 0.172 for SwissProt data.
**Range:** 1/sequence_length to 1.0

Longest Run
-----------

.. autofunction:: peptides.Peptide.longest_run

The longest run identifies the length of the longest stretch of consecutive identical amino acids. This metric is useful for detecting repetitive regions, low complexity sequences, and potential sequencing artifacts.

**Examples:**

.. code-block:: python

    from peptides import Peptide
    
    # Single amino acid sequence
    peptide = Peptide("AAAA")
    print(peptide.longest_run())  # Output: 4
    
    # Mixed sequence with runs
    peptide = Peptide("AALLLSSS")
    print(peptide.longest_run())  # Output: 3
    
    # No consecutive identical amino acids
    peptide = Peptide("ACDEFGHIKLMNPQRSTVWY")
    print(peptide.longest_run())  # Output: 1

**Interpretation:** Most values fall between 2 and 5 for SwissProt data.
**Range:** 1 to sequence_length

Outlier Detection
----------------

.. autofunction:: peptides.Peptide.detect_outlier

The outlier detection function analyzes a peptide sequence using composition-based vetting metrics and compares them against established SwissProt protein distributions to identify potential outliers, artifacts, or unusual sequences. This function provides an automated way to flag sequences that may require further investigation.

**How it works:**

The function evaluates the following metrics against 5th and 95th percentile thresholds derived from SwissProt protein analysis:

1. **Entropy**: Shannon entropy of amino acid composition
2. **Max Frequency**: Maximum frequency of any single amino acid
3. **Longest Run**: Length of longest consecutive identical amino acids

Note: The result also returns the sequence length in the `metrics` field for reference, but it is not used to determine outliers.

**Examples:**

.. code-block:: python

    from peptides import Peptide
    
    # Normal protein sequence
    normal_peptide = Peptide("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAIGLNKALELVKQKLKELN")
    result = normal_peptide.detect_outlier()
    print(f"Is outlier: {result['is_outlier']}")  # Output: False
    print(f"Issues: {result['issues']}")          # Output: []
    
    # Problematic sequence
    problematic_peptide = Peptide("AAAA")
    result = problematic_peptide.detect_outlier()
    print(f"Is outlier: {result['is_outlier']}")  # Output: True
    print(f"Issues: {result['issues']}")          # Output: ['Entropy 0.000 below 5th percentile (3.714)', 'Max frequency 1.000 above 95th percentile (0.172)']

**Return Value:**

The function returns a dictionary with three keys:

- **is_outlier**: Boolean indicating if the sequence is flagged as an outlier
- **issues**: List of specific problems found (empty if no issues)
- **metrics**: Dictionary containing the calculated values for all four metrics

**Thresholds (based on SwissProt analysis):**

- **Entropy**: 5th percentile = 3.714, 95th percentile = 4.185
- **Max Frequency**: 5th percentile = 0.085, 95th percentile = 0.172
- **Longest Run**: 5th percentile = 2.0, 95th percentile = 5.0

**Use Cases:**

1. **Quality Control**: Automatically flag problematic sequences in large datasets
2. **Sequence Validation**: Ensure sequences meet expected biological parameters
3. **Artifact Detection**: Identify potential sequencing errors or contaminants
4. **Database Filtering**: Remove or flag unusual sequences before analysis
5. **Research Validation**: Verify sequence quality in experimental workflows

**Best Practices:**

- Use as part of automated quality control pipelines
- Review flagged sequences manually to confirm issues
- Consider context when interpreting results (some unusual sequences may be biologically valid)
- Combine with other validation methods for comprehensive quality assessment

**References:**

- Shannon, C. E. (1948). *A Mathematical Theory of Communication*. Bell System Technical Journal, 27(3), 379-423.