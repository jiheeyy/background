rank1 <-function(X_rs, X_cs, l0, f0, pm_func,
                 ql, gl, kl_l,
                 qf, gf, kl_f){
  # Update f first
  s = sum(l0 * ql$mean) * f0
  fit = do.call(pm_func$f, c(list(x = X_cs, s = s, shape = "estimate", scale = "max", g_init = gf, fix_g = FALSE)))
  qf = fit$posterior
  gf = fit$fited_g
  kl_f = compute_kl_ebpm(y = X_cs, s = s, posterior = qf, ll = fit$log_likelihood)
  rm(fit)
  rm(s)
  qg = list(ql = ql, gl = gl, kl_l = kl_l, qf = qf, gf = gf, kl_f = kl_f)
  return(qg)}


# PASTED FROM what g_init is
scale2gammamix_init <- function(shape, scale){
  n = length(shape)
  return(gammamix(pi = replicate(n, 1)/n, shape = shape, scale =  scale))
}

# PASTED FROM EBPM_GAMMA_MIXTURE
b = 1/g_init$scale ##  from here use gamma(shape = a, rate = b)  where E = a/b
a = g_init$shape
tmp <-  compute_L(x,s,a, b)
L =  tmp$L
l_rowmax = tmp$l_rowmax
fit <- try(mixsqp(L, x0 = g_init$pi, control = control))
if (class(fit) == "try-error"){
  g_init$pi = g_init$pi + 1e-8
  fit <- try(mixsqp(L, x0 = g_init$pi, control = control))
}
pi = fit$x
fitted_g = gammamix(pi = pi, shape = a,  scale  = 1/b)
log_likelihood = sum(log(exp(l_rowmax) * L %*%  pi))



### Previously installed packages ####
# install.packages("matrixStats")
# install_github("stephenslab/ebpm")
# install_github('linxihui/NNLM')
# install_github("stephenslab/ebpmf.alpha")
# install.packages("pheatmap")

#### BG - Load Packages ####

# rm(list = ls())
# library(devtools)
library(ebpm)
library(Matrix)
library(fastTopics)
library(pheatmap)
# library(NNLM)
#library(ebpmf.alpha)
#library(flashier)
#source("/Users/jiheeyou/ebpmf.alpha/R/ebpmf.R", echo = TRUE)

#### BG - Generate Synthetic Data, Initialize PNMF, Visualization ####

set.seed(123)
n = 12
p = 30
K = 3
LL = matrix(0, nrow=n, ncol=K)
FF = matrix(0, nrow=K, ncol=p)
LL[1:(n/3),1] = 1
LL[((n/3)+1):(2*n/3),2] = 1
LL[((2*n/3)+1):n,3] = 1
LL = LL + matrix(runif(n*K,0,0.5),nrow=n)
FF[1,1:(p/3)] = 1+10
FF[2,((p/3)+1):(2*p/3)] = 1+10
FF[3,((2*p/3)+1):p] = 1+10
FF = FF + matrix(rnorm(p*K,0,1),ncol=p)
FF = pmax(FF,0)
lambda = LL %*% FF
X = matrix(rpois(n=length(lambda),lambda),nrow=n)

# Plot the true loadings matrix LL,FF
pheatmap(LL,
         main = "True Loadings (LL)",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         xlab = "Factor", ylab = "Sample")
pheatmap(FF,
         main = "True Factors (FF)",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         xlab = "Feature", ylab = "Factor")

bg_prior <- function(){
  aL = c(seq(0.01, 0.10, 0.01), seq(0.2, 0.9, 0.1), seq(1,15,2), 20, 50, 75, 100, 200, 1e3, 1e-8, 1e-16)
  D = length(aL)
  g = ebpm::gammamix(pi = replicate(D, 1/D), shape = aL, scale = 1/aL)
  return(g)
}

compute_kl_ebpm <- function(y,s, posterior, ll){
  if(length(s) == 1){s = replicate(length(y),s)}
  mask <- (y != 0)
  E_loglik = - sum(s * posterior$mean) + sum(y[mask] * log(s[mask])) + sum(y[mask]*posterior$mean_log[mask]) - sum(lgamma(y[mask] + 1))
  return(E_loglik - ll)
}

update_qg <- function(tmp, qg, k){
  qg$qls_mean[,k] = tmp$ql$mean
  qg$qls_mean_log[,k] = tmp$ql$mean_log
  qg$qfs_mean[,k] = tmp$qf$mean
  qg$qfs_mean_log[,k] = tmp$qf$mean_log
  qg$gls[k] = list(tmp$gl)
  qg$gfs[k] = list(tmp$gf)
  qg$kl_l[k] = tmp$kl_l
  qg$kl_f[k] = tmp$kl_f
  return(qg)
}

ebpmf_init = function(X, K, X_dim, init='fastTopics', maxiter = 50, method = 'topicscore',
                      pm_func = list(f = ebpm::ebpm_gamma_mixture, l = ebpm::ebpm_gamma_mixture)){
  if(init == 'fastTopics'){

    # Initial Fitting: Poisson NMF.
    # Note: NMF also results in negative b_k_max
    init_fit = fastTopics::fit_poisson_nmf(X,K,numiter=maxiter, init.method=method)
    L = init_fit$L
    F = init_fit$F
  }

  # Normalize, Add eps to L, F
  L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
  l0 = apply(L, 1, mean)
  f0 = apply(F, 1, mean)
  L = L/l0
  F = F/f0
  L[L == 0] = 1e-10
  F[F == 0] = 1e-10



  qg = list(qls_mean = L, qls_mean_log = log(L), kl_l = replicate(K,0),
            qfs_mean = F, qfs_mean_log = log(F), kl_f = replicate(K,0),
            gls = replicate(ncol(L0), list(bg_prior())), gfs = replicate(ncol(L0), list(bg_prior())))

  # Compute a, the rolling maximum of factor k's contribution to lambda
  a = replicate(length(X_dim$x), -Inf) # CHANGE from zero

  for (k in 1:K){
    b_k_tmp <- qg$qls_mean_log[X_dim$i, k] + qg$qfs_mean_log[X_dim$j, k]
    a <- pmax(a, b_k_tmp) # b_k_tmp is usually negative
  }

  # Compute b: at the end, lambda_ij = b + a
  b = qg$qls_mean_log[X_dim$i, 1] + qg$qfs_mean_log[X_dim$j, 1] - a
  for (k in 2:K){
    b_k = qg$qls_mean_log[X_dim$i, k] + qg$qfs_mean_log[X_dim$j, k] - a
    b <- log(exp(b) + exp(b_k))
  }

  return(list(l0=l0, f0=f0, qg=qg, b=b, a=a))
}

rank1 <-function(X_rs, X_cs, l0, f0, pm_func,
                 ql, gl, kl_l,
                 qf, gf, kl_f){
  # Update f first
  s = sum(l0 * ql$mean) * f0
  fit = do.call(pm_func$f, c(list(x = X_cs, s = s, shape = "estimate", scale = "max", g_init = gf, fix_g = FALSE)))
  qf = fit$posterior
  gf = fit$fited_g
  kl_f = compute_kl_ebpm(y = X_cs, s = s, posterior = qf, ll = fit$log_likelihood)
  rm(fit)
  rm(s)

  s = sum(f0 * qf$mean) * l0
  fit = do.call(pm_func$l, c(list(x = X_rs, s = s, shape = "estimate", scale = "max",g_init = gl, fix_g = FALSE)))
  ql = fit$posterior
  gl = fit$fited_g
  kl_l = compute_kl_ebpm(y = X_rs, s = s, posterior = ql, ll = fit$log_likelihood)
  rm(fit)
  rm(s)
  qg = list(ql = ql, gl = gl, kl_l = kl_l, qf = qf, gf = gf, kl_f = kl_f)

  return(qg)}

#### BG -  ####

ebpmf_bg = function(X, K){
  X = Matrix(X, sparse=TRUE)
  X_dim = Matrix::summary(X)
  X_rs = Matrix::rowSums(X)
  X_cs = Matrix::colSums(X)

  ll_const = sum(apply.nonzeros(X = X, f = function(x) lgamma(x + 1))) # Poisson log likelihood constant term
  # nonzero_idx = cbind(x$i,x$j)
  # nonzero indices - 338 out of 360 print(nrow(nonzero_idx))

  init_tmp = ebpmf_init(X,K,X_dim)
  qg = init_tmp$qg
  b = init_tmp$b
  a = init_tmp$b
  l0 = init_tmp$l0
  f0 = init_tmp$f0
  rm(init_tmp)

  # Main loop
  ELBOs = c()
  KLs = c()

  # point mass 0 L, 1 F
  # math part: z sum to x, does that ressult in conv cancel otu restrict var
  for (i in 1:maxiter){
    b_k_max = replicate(length(X_dim$x),-Inf)
    for (k in 1:K){
      b_k = qg$qls_mean_log[X_dim$i,k] + qg$qfs_mean_log[X_dim$j, k] - a
      Ez_matrix = sparseMatrix(i = X_dim$i, j = X_dim$j, x = (X_dim$x * exp(b_k - b)))
      Ez = list(rs = Matrix::rowSums(Ez_matrix), cs = Matrix::colSums(Ez_matrix))

      # qg update step, assuming neither L nor F are fixed.
      rank1_qg = rank1(X_rs = Ez$rs, X_cs = Ez$cs, l0 = l0, f0 = f0, pm_func = list(f = ebpm::ebpm_gamma_mixture, l = ebpm::ebpm_gamma_mixture),
                       ql = list(mean = qg$qls_mean[,k],
                                 mean_log = qg$qls_mean_log[,k]),
                       qf = list(mean = qg$qfs_mean[,k],
                                 mean_log = qg$qfs_mean_log[,k]),
                       gl = qg$gls[[k]],
                       gf = qg$gfs[[k]],
                       kl_l = qg$kl_l[k],
                       kl_f = qg$kl_f[k])
      qg = update_qg(tmp = rank1_qg, qg=qg, k = k)
      rm(rank1_qg)

      # Update b at the end of each k loop
      b_k0 = b_k
      b_k = qg$qls_mean_log[X_dim$i,k] + qg$qfs_mean_log[X_dim$j, k] - a
      b = log( exp(b) - exp(b_k0) + exp(b_k)  )
      b_k_max = pmax(b_k, b_k_max)
    }
    # Update l0, f0. Calculate ELBO at the end of each iteration
    denom <- colSums( t(qg$qls_mean) * colSums(f0 * qg$qfs_mean))
    l0 <- X_rs/denom
    denom <- colSums( t(qg$qfs_mean) * colSums(l0 * qg$qls_mean))
    f0 <- X_cs/denom

    KL = sum(qg$kl_l) + sum(qg$kl_f)
    ELBO = - sum( colSums(l0 * qg$qls_mean) * colSums(f0 * qg$qfs_mean) ) + sum(X_dim$x * (log(l0[X_dim$i]) + log(f0[X_dim$j]) + b + a) ) - KL - ll_const
    ELBOs <- c(ELBOs, ELBO)
    KLs <- c(KLs, KL)
    print("iter         ELBO")
    print(sprintf("%d:    %f", i, ELBO))
  }
  return(list(l0 = l0, f0 = f0, qg = qg, ELBO = ELBOs, KL = KLs))
}


start_time = Sys.time()
ebpmf_bg(X=X, K=3)
end_time = Sys.time()
print(end_time - start_time)

