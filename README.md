# Power Analysis

## Exome sequencing depth required for sensitivity for low Variant Allele Frequency

normally a simple binomial given pr success (VAF) and n draws (coverage), 
however more realistic to also model the coverage as a random variable too (poisson)
ignoring the binomial sampling of cells from tissue - but large n of cells should
make this variance small and ignorable

[code is here](/power_analysis.R)

Justifies 500X versus 100X: with variant allele frequency of 1%, 87% power to find variants at 500X; drops to only 8% for 100X ;(!

## Mutational signature identification and quantification

Followed method as outlined in S1 text from SparseSignatures paper, "De novo mutational signature discovery in tumor genomes using SparseSignatures", by Avantika Lal,Keli Liu,Robert Tibshirani,Arend Sidow,Daniele Ramazzotti.
 
[paper is here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009119#abstract0)

simulated 6 high + low UV exposed eye quadrants using COSMIC 7a and 7b signatures on Lal defined clock like muational backgrounds.

SparseSignatures correctly identified signatures, 'S3' with highest correlation with COSMIC 7a/7b.
Correctly identified samples from higher UV region.
 See results below (confusingly 'S3' is labelled Siganture 2 in the plot - must be a minor bug in the R package):

![identified_signatures](/mut_signature.png)

[code is here](/mutation_signature_simulation.R)
