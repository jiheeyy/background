# Some explorations of the "background" topic model using the
# newsgroups data.
library(Matrix)
library(fastTopics)
k <- 10
set.seed(1)
load("../data/newsgroups.RData")
t0 <- proc.time()
fit0 <- fit_poisson_nmf(counts,k = k,numiter = 40,method = "em",
                        init.method = "random",verbose = "detailed",
                        control = list(nc = 4))
fit_scd <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 100,method = "scd",
                           control = list(nc = 4,extrapolate = TRUE),
                           verbose = "detailed")
t1 <- proc.time()
print(t1 - t0)
