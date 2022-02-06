# power to detect different variant allele frequency VAFs
# normally modelled as a simple binomial distribution given pr. success (VAF) and n draws (coverage), 
# however more realistic to also model the coverage as a random variable too (poisson)
# ignoring the binomial sampling of cells from tissue - but large n of cells should
# make this variance small and ignorable
set.seed(97)
n_sites <- 5000
coverage <- 500
min_supporting_reads <- 3
vafs <- c(1,5,10)
power <- vector()
for (v in vafs){
  realised_coverage_per_site <- rpois(n_sites, coverage)
  alt_reads_per_site <- sapply(realised_coverage_per_site, function(x) rbinom(1,size = x,prob = v/100))
  power <- c(power,length(which(alt_reads_per_site >= min_supporting_reads))/n_sites*100)
}
names(power) <- paste0(vafs,"%")
#> power
#1%     5%    10% 
#87.42 100.00 100.00 

# only 100X????
coverage <- 100
power <- vector()
for (v in vafs){
  realised_coverage_per_site <- rpois(n_sites, coverage)
  alt_reads_per_site <- sapply(realised_coverage_per_site, function(x) rbinom(1,size = x,prob = v/100))
  power <- c(power,length(which(alt_reads_per_site >= min_supporting_reads))/n_sites*100)
}
names(power) <- paste0(vafs,"%")
#1%    5%   10% 
#8.24 87.48 99.70 

# 4 genes, each with 25% expected frequency. 25 patients
n_genes <- 4
n_patients <- 6
n_sims <- 10000
test_dist <- rmultinom(n=n_sims, size=n_patients, prob=rep(1/n_genes,n_genes))
1-(mean(rowSums(test_dist == 0))/n_sims)

n_genes <- 10
test_dist <- rmultinom(n=n_sims, size=n_patients, prob=rep(1/n_genes,n_genes))
1-(mean(rowSums(test_dist == 0))/n_sims)




rnbinom(1, mu = 200, size = 1.5)
