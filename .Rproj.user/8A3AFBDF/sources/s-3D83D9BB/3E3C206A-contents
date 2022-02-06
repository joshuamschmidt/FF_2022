# simulating mutational signatures
# based off of "De novo mutational signature discovery in tumor genomes using 
# SparseSignatures"
# Avantika Lal,Keli Liu,Robert Tibshirani,Arend Sidow,Daniele Ramazzotti 
# No code! So working from written description in Supp S1.text
library("SparseSignatures")
library("BSgenome.Hsapiens.1000genomes.hs37d5")

# need to know data structure to simulate, so use vignette example
bsg = BSgenome.Hsapiens.1000genomes.hs37d5
data(mutation_categories)
head(mutation_categories)

imported_data = import.trinucleotides.counts(data=ssm560_reduced,reference=bsg)
head(imported_data)
is.matrix(imported_data)
str(imported_data)
# so data is a n sample * 96 mut class matrix

# get mutation class from COSMIC: 
#cancer.sanger.ac.uk/signatures/documents/453/COSMIC_v3.2_SBS_GRCh38.txt
cosmic <- read.table(file='COSMIC_v3.2_SBS_GRCh38.txt',header = T,sep="\t",
                     row.names = 1 )

norm_cosmic <- apply(cosmic,MARGIN = 2, function(x) x/sum(x))
# sSignatures
# Lal et. al. find that SBS5 is analogous to their independently determined
# background/clock like mutation signature
# SBS 1 is another clock like signature (deamination of methylcytosine)
# UV are SBS 7a, 7b.
# following Lal et. al.
# SBS background rnbinom(1, mu = 2000, size = 1.5)
# SBS1 rnbinom(1, mu = 400, size = 1.5)
# UV signature rnbinom(1, mu = 400, size = 1.5), or 10% increase rnbinom(1, mu = 220, size = 1.5)

generate_sample_signature <- function(uv_increase=20/100){
  uv_increase <- round(uv_increase*200)
  signature <- round(rnbinom(1, mu = 2500, size = 1.5) * norm_cosmic[,"SBS5"]) +
  round(rnbinom(1, mu = 100, size = 1.5) * norm_cosmic[,"SBS1"]) +
  round(rnbinom(1, mu = 400+uv_increase, size = 1.5) * norm_cosmic[,"SBS7a"]) +
  round(rnbinom(1, mu = 400+uv_increase, size = 1.5) * norm_cosmic[,"SBS7b"])
  # noise
  noise <- round(runif(96,0,25))
  return(signature + noise)
}

# generate simulated data set. Assume 50% increase in UV exposure between quadrants
set.seed(42)
n_samples_per_quadrant <- 6
affected <- t(sapply(1:n_samples_per_quadrant, function(x) generate_sample_signature(uv_increase=50/100),USE.NAMES = T))
rownames(affected) <- paste0("A",1:n_samples_per_quadrant)
non_affected <- t(sapply(1:n_samples_per_quadrant, function(x) generate_sample_signature(uv_increase=0/100),USE.NAMES = T))
rownames(non_affected) <- paste0("NA",1:n_samples_per_quadrant)
total <- rbind(affected, non_affected)
colnames(total) <- mutation_categories$cat

# initial beta
starting_betas <- startingBetaEstimation(x=total,
                                         K=3:5,
                                         background_signature=norm_cosmic[,"SBS5"])

beta = starting_betas_example[["3_signatures","Value"]]
res = nmfLasso(x = total, K = 2, beta = beta, background_signature = norm_cosmic[,"SBS5"], seed = 97)
signatures = res$beta
signatures.plot(beta=signatures, xlabels=FALSE)


sort(apply(norm_cosmic,MARGIN = 2, function(x) cor.test(signatures['Background',],x)$estimate))

sort(apply(norm_cosmic,MARGIN = 2, function(x) cor.test(signatures['S1',],x)$estimate))
sort(apply(norm_cosmic,MARGIN = 2, function(x) cor.test(signatures['S2',],x)$estimate))
sort(apply(norm_cosmic,MARGIN = 2, function(x) cor.test(signatures['S3',],x)$estimate))
# UV identified as signature 3. NB incorrectly plotted as Signature 2!
# check with apply(signatures,MARGIN=1,function(x) which(x==max(x)))
# rho
sort(apply(norm_cosmic,MARGIN = 2, function(x) cor.test(signatures['S4',],x)$estimate))


c <- cor.test(signatures['Background',],norm_cosmic[,"SBS1"])


