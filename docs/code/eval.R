# Evaluate L, F using mean absolute error from
# ground truth L, F (if known from simulation)

library(lsa)
library(ggplot2)

eval_mae <- function(est_L, est_F, true_L, true_F, X){
  est_loglam = compute_log_Lambda_mat(log_L = log(est_L),
                                      log_F = log(est_F),
                                      X = X)
  loglam = compute_log_Lambda_mat(log_L = log(true_L),
                                  log_F = log(true_F),
                                  X = X)
  mae = mean(abs(exp(est_loglam) - exp(loglam)))
  return(mae)
}

eval_csim <- function(est_L, est_F, true_L, true_F){

  k = ncol(true_L)
  estk = ncol(est_L)
  col_names <- c(paste0("true_", 1:k), paste0("est_", 1:estk))

  combined_L <- cbind(true_L, est_L)
  colnames(combined_L) <- col_names
  cos_L <- lsa::cosine(combined_L)

  combined_F <- cbind(true_F, est_F)
  colnames(combined_F) <- col_names
  cos_F <- lsa::cosine(combined_F)

  return(list(cos_L = cos_L, cos_F = cos_F))
}

# Input 1 matrix, text title
plot_upper_heatmap <- function(sim_mat, plot_title) {
  rows <- row(sim_mat)
  cols <- col(sim_mat)
  upper_tri_mask <- rows <= cols

  df <- data.frame(
    Row = rownames(sim_mat)[rows[upper_tri_mask]],
    Col = colnames(sim_mat)[cols[upper_tri_mask]],
    Value = sim_mat[upper_tri_mask]
  )

  df$Row <- factor(df$Row, levels = rev(rownames(sim_mat)))
  df$Col <- factor(df$Col, levels = colnames(sim_mat))

  # --- NEW: Dynamically find the quadrant boundaries ---
  # Count how many "true" and "est" factors exist
  num_true <- sum(grepl("^true_", colnames(sim_mat)))
  num_est  <- sum(grepl("^est_", colnames(sim_mat)))

  # X-axis goes left-to-right. The line goes after the last "true_"
  x_intercept <- num_true + 0.5

  # Y-axis goes bottom-to-top (because we reversed the levels).
  # The line goes above the highest "est_"
  y_intercept <- num_est + 0.5

  # Generate the plot
  p <- ggplot(df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", mid = "white",
                         midpoint = 0, limit = c(0, 1),
                         name = "Cosine\nSimilarity",
                         na.value = "grey85") +
    geom_text(aes(label = ifelse(is.na(Value), "NaN", sprintf("%.2f", Value))),
              color = "black", size = 3.5) +

    # --- NEW: Draw the quadrant lines ---
    geom_vline(xintercept = x_intercept, color = "black", linewidth = 1) +
    geom_hline(yintercept = y_intercept, color = "black", linewidth = 1) +

    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    coord_fixed() +
    ggtitle(plot_title)

  return(p)
}
