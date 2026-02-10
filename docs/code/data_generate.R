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


generate_synthetic <- function(data_type, w_matrix = NA, n=12, p=30, K=3){
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

  # n, p must be multiples of 3 for indexing to not break
  if (data_type == "clear_membership"){
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

# Optional: Plot the true loadings matrix LL,FF
#library(pheatmap)
visualize_lf <- function(L, F, lf_type){
  pheatmap(L,
           main = paste(lf_type, " Loadings"),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           xlab = "Factor", ylab = "Sample")
  pheatmap(F,
           main = paste(lf_type," Factors"),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           xlab = "Feature", ylab = "Factor")
}
