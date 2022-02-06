# Power Analysis

## Exome sequencing depth

normally a simple binomial given pr success (VAF) and n draws (coverage), 
however more realistic to also model the coverage as a random variable too (poisson)
ignoring the binomial sampling of cells from tissue - but large n of cells should
make this variance small and ignorable

[](power_analysis.R)
## Mutational signature identification and quantification

Followed method as outlined in S1 text from SparseSignatures paper, "De novo mutational signature discovery in tumor genomes using SparseSignatures", by Avantika Lal,Keli Liu,Robert Tibshirani,Arend Sidow,Daniele Ramazzotti 

code in mutation_signature_simulation.R

simulated 6 high + low UV exposed eye quadrants using COSMIC 7a and 7b signatures on Lal defined clock like muational backgrounds.

SparseSignatures correctly identified signatures, 'S3' with highest correlation with COSMIC 7a/7b.
Correctly identified samples from higher UV region.
