# Create gamma mixture, each component mean is 1
init_prior <- function(){
  len_grid = length(shape_grid)
  g = gammamix()
  return(g)
}

initialize <- function(X, K,
                       # 3 parameters used in initial fitting
                       init='fastTopics', init_iter = 50){

  X_dim = Matrix::summary(X)

  if(init == 'fastTopics'){
    # Initial Fitting: Poisson NMF.
    # Note: NMF also results in negative b_k_max
    init_fit = fastTopics::fit_poisson_nmf(X, K, numiter=init_iter)
    L = init_fit$L
    F = init_fit$F
  }

  # Normalize L, F
  # Then add eps
  L[L <  1e-8] <- 1e-8
  F[F <  1e-8] <- 1e-8
  l_mean = apply(L, 1, mean) # get mean of each row (res: length n vector)
  f_mean = apply(F, 1, mean) # get mean of each row (res: length p vector)
  L = L/l_mean
  F = F/f_mean
  L[L == 0] = 1e-10
  F[F == 0] = 1e-10

  qg = list(ql_norm = L, qllog_norm = log(L), kl_l = replicate(K,NA),
            qf_norm = F, qflog_norm = log(F), kl_f = replicate(K,NA),
            gl = replicate(K, list(gammamix())), gl = replicate(K, list(gammamix())))

  # Compute a, the rolling maximum of factor k's contribution to lambda_ij
  # length of a is number of nonzero elements in X
  a = replicate(nrow(X_dim), -Inf)

  for (k in 1:K){
    b_k_tmp <- qg$qllog_norm[X_dim$i, k] + qg$qflog_norm[X_dim$j, k] # log(l_ik * f_jk) = log(lambda_ijk)
    a <- pmax(a, b_k_tmp) # b_k_tmp is usually negative
  }

  # Update b over k = 2 to K
  b = qg$qllog_norm[X_dim$i, 1] + qg$qflog_norm[X_dim$j, 1] - a # b = log(lambda_ij1) - a
  for (k in 2:K){
    b_k = qg$qllog_norm[X_dim$i, k] + qg$qflog_norm[X_dim$j, k] - a
    b <- log(exp(b) + exp(b_k))} # b <- b + b_k
  # At the end, log(lambda_ij) = b + a

  return(list(qg=qg, b=b, a=a))
}
