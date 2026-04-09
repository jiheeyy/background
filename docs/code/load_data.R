#### BG - Generate Synthetic Data, Visualization ####

normalize_w_matrix <- function(w) {
  M <- nrow(w)
  cs <- colSums(w)
  # Identify indices where sum is 0
  zero_cols <- which(cs == 0)

  if (length(zero_cols) > 0) {
    w[, zero_cols] <- 1/M
    # Update column sums for the next step so we don't divide by zero
    cs[zero_cols] <- 1
  }

  # 3. Normalize: Divide every value by its column sum
  # scale() divides columns by the 'scale' argument
  w_norm <- scale(w, center = FALSE, scale = cs)

  # Remove the "scaled:scale" attribute that R adds, just for cleanliness
  attr(w_norm, "scaled:scale") <- NULL

  return(as.matrix(w_norm))
}

###########################
# Define mixture weights for each L_k and F_k
# 21st mixture component (shape, rate) value is 1
l1 <- c(1, rep(0, 33)) # shape = 1e-16
l2 <- c(rep(0, 20),1, rep(0,13)) # shape = 1
l3 <- c(rep(0, 33), 1) # shape = 1000
f1 <- c(1, rep(0, 33))
f2 <- c(rep(0, 20),1, rep(0,13))
f3 <- c(rep(0, 33), 1)
# Weights will be normalized by column inside visualize_mixtures or generate_synthetic
default_w_matrix <- cbind(l1,l2,l3,f1,f2,f3)
###########################

general_load <- function(data_type=NA, w_matrix = default_w_matrix, n=12, p=30, K=3){
  if (data_type == "mix"){
    # w_matrix is M by K
    w_matrix <- normalize_w_matrix(w_matrix)
    colnames(w_matrix) <- c(paste0("L", 1:K), paste0("F", 1:K))

    LL = matrix(0, nrow=n, ncol=K)
    FF = matrix(0, nrow=p, ncol=K)

    grid_info = gammamix()
    shapes <- grid_info$shape
    rates  <- grid_info$rate
    len_grid <- length(shapes)

    for (k in 1:K){
      w_L <- w_matrix[, k]
      # Pick mixture component for each of n values in L_k
      idx_L <- sample(1:len_grid, size = n, replace = TRUE, prob = w_L)
      LL[, k] <- rgamma(n, shape = shapes[idx_L], rate = rates[idx_L])

      w_F <- w_matrix[, K + k]
      idx_F <- sample(1:len_grid, size = p, replace = TRUE, prob = w_F)
      FF[, k] <- rgamma(p, shape = shapes[idx_F], rate = rates[idx_F])
    }

    lambda = LL %*% t(FF)
    X = matrix(rpois(n=length(lambda),lambda),nrow=n)
    return(list(L = LL, F = FF, X = X))
  }

  if (data_type == "unif"){
    F <- rand(p,K) # fasttopics uses 4 args, but current rand requires 2
    L <- rand(n,K)
    X <- matrix(as.double(rpois(n*p,tcrossprod(L,F))),n,p)
    X <- as(X,"CsparseMatrix")
    return(list(X=X, L=L, F=F))
  }

  # n, p must be multiples of 3 for indexing to not break
  if (data_type == "clear"){
    LL = matrix(0, nrow=n, ncol=K)
    FF = matrix(0, nrow=p, ncol=K)
    LL[1:(n/3),1] = 1
    LL[((n/3)+1):(2*n/3),2] = 1
    LL[((2*n/3)+1):n,3] = 1
    LL = LL + matrix(runif(n*K,0,0.5),nrow=n)
    FF[1:(p/3), 1] = 1+10
    FF[((p/3)+1):(2*p/3), 2] = 1+10
    FF[((2*p/3)+1):p, 3] = 1+10
    FF = FF + matrix(rnorm(p*K,0,1),ncol=K)
    FF = pmax(FF,0)
    lambda = LL %*% t(FF)
    X = matrix(rpois(n=length(lambda),lambda),nrow=n)
    return(list(L = LL, F = FF, X = X))
  }

  if (data_type == "news"){
    load(here("docs","data","newsgroups.RData"))
    X <- counts
    y <- topics
    return(list(X = X,y=y))

  }

  if (data_type == "gtex"){

  }
}

# library(ggplot2)
visualize_mixtures <- function(w_mat) {
  # 1. Get Grid Parameters
  w_mat <- normalize_w_matrix((w_mat))
  grid_params <- gammamix() # shape = grid, rate = grid
  shapes <- grid_params$shape
  rates  <- grid_params$rate

  # 2. Define X range for plotting
  # Since mean is 1, range 0 to 4 covers most relevant density
  x_vals <- seq(1e-4, 4, length.out = 500)

  # 3. Calculate Densities
  plot_data <- data.frame()

  # Loop through each column (L1, L2... F3)
  for (col_name in colnames(w_mat)) {
    weights <- w_mat[, col_name]

    # Calculate mixture density: Sum(weight_i * dgamma(x, shape_i, rate_i))
    y_vals <- sapply(x_vals, function(x) {
      densities <- dgamma(x, shape = shapes, rate = rates)
      sum(weights * densities)
    })

    # Store
    tmp <- data.frame(x = x_vals, density = y_vals, Component = col_name)
    plot_data <- rbind(plot_data, tmp)
  }

  # 4. Plot with ggplot2
  ggplot(plot_data, aes(x = x, y = density, color = Component, fill = Component)) +
    geom_line(linewidth = 1) +
    geom_area(alpha = 0.2, position = "identity") +
    facet_wrap(~Component, scales = "free_y") +
    theme_light() +
    labs(title = "Theoretical Gamma Mixture Densities",
         subtitle = "Visualizing the 'Spikiness' of each defined vector",
         x = "Value (Mean centered at 1)",
         y = "Density") +
    theme(legend.position = "none")
}

### TODO KUSHAL
# Downloaded from https://github.com/stephenslab/count-clustering/blob/master/project/rdas/gtexv6brain.k6fit.rda
# k6fit selected since this is the version visualized in final paper

# Fitting results, not the raw data
# Raw data could be https://github.com/kkdey/GTExV6Brain/tree/master/data, but not sure

# ---

gom_model_fit <- get(load(here("docs","data","gtexv6brain.k6fit.rda")))
omega <- gom_model_fit$omega
colnames(omega) <- c(1:NCOL(omega))

### Paused since I can't find tissue name index txt (samples_id.txt in Kushal's code)
