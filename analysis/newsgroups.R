# Some explorations of the "background" topic model using the
# newsgroups data.
library(Matrix)
library(fastTopics)
source("../docs/code/bg_implement.R")
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
n <- nrow(fit_scd$progress)
loglik <- max(fit_scd$progress$loglik)
plot(1:n,log10(loglik - fit_scd$progress$loglik + 1e-4),type = "l",
     xlab = "iteration",ylab = "dist. to best solution")

# Visualize the fastTopics output using a Structure plot.
set.seed(1)
topic_colors <- c(k1 = "darkmagenta",k2 = "limegreen",k3 = "darkorange",
                  k4 = "gold",k5 = "magenta",k6 = "forestgreen",
                  k7 = "darkblue",k8 = "skyblue",k9 = "red",
                  k10 = "white")
p <- structure_plot(fit_scd,grouping = topics,topics = c(10,8,1:7,9),
                    gap = 10,colors = topic_colors)
print(p)

# Fit a "background" topic model.
set.seed(1)
a <- c(0.1,1,10,100)
state0 <- initialize_state(counts,k,a,L_init = "gamma")
t0 <- proc.time()
fit_btm <- run_ebpmf_bg_cavi(state0,counts,a,max_niter = 100,
                             elbo_every = 1,early_stop = FALSE)
t1 <- proc.time()
print(t1 - t0)
