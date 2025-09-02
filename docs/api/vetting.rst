Sequence Vetting
================

.. currentmodule:: peptides

The sequence vetting functions provide tools for analyzing the complexity of 
protein sequences. These functions are particularly useful for quality control.


Outlier Detection
-----------------

The outlier detection function (`Peptide.detect_outlier`) analyzes a peptide 
sequence using composition-based vetting metrics and compares them against 
established SwissProt protein distributions to identify potential outliers, 
artifacts, or unusual sequences. This function provides an automated way to 
flag sequences that may require further investigation.

How it works
^^^^^^^^^^^^

The function evaluates the following metrics against 5th and 95th percentile thresholds derived from SwissProt protein analysis:

1. **Entropy**: Shannon entropy of amino acid composition
2. **Max Frequency**: Maximum frequency of any single amino acid
3. **Longest Run**: Length of longest consecutive identical amino acids

Use Cases
^^^^^^^^^

1. **Quality Control**: Automatically flag problematic sequences in large datasets.
2. **Sequence Validation**: Ensure sequences meet expected biological parameters.
3. **Artifact Detection**: Identify potential sequencing errors or contaminants.
4. **Database Filtering**: Remove or flag unusual sequences before analysis.
5. **Research Validation**: Verify sequence quality in experimental workflows.

Best Practices
^^^^^^^^^^^^^^

- Use as part of automated quality control pipelines.
- Review flagged sequences manually to confirm issues.
- Consider context when interpreting results (some unusual sequences may be biologically valid).
- Combine with other validation methods for comprehensive quality assessment.

Result
^^^^^^

.. autoclass:: peptides.OutlierResult(typing.NamedTuple)