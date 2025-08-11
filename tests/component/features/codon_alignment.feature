Feature: Codon alignment on real sequences
  Scenario: Align sample sequences to reference
    Given sample sequences "tests/component/data/samplesmall.fas" and reference "tests/component/data/hiv1_reference.fasta"
    And minimap2 from "https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2" is available
    When they are codon-aligned with options "-k 6 -w 3 --score-N 0 --secondary no" from 790 to 2085
    Then the number of aligned pairs is 5
