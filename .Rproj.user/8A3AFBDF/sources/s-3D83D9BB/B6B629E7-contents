# simulate multinomial draws, becuase I dont know yet how to do it with density functions.....


# 4 genes, each with 25% expected frequency. 25 patients
n_genes <- 4
n_patients <- 25
n_sims <- 10000
test_dist <- rmultinom(n=n_sims, size=n_patients, prob=rep(1/n_genes,n_genes))
1-(mean(rowSums(test_dist == 0))/n_sims)

n_genes <- 10
test_dist <- rmultinom(n=n_sims, size=n_patients, prob=rep(1/n_genes,n_genes))
1-(mean(rowSums(test_dist == 0))/n_sims)




rnbinom(1, mu = 200, size = 1.5)
