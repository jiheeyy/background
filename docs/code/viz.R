library(pheatmap)
library(gridExtra)

# Plot matrix L,F
plot_lf <- function(L, F, lf_str, breaks=seq(0, 1, length.out = 101)){
  pl <- pheatmap(L,
           main = paste(lf_str, " Loadings"),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           silent = TRUE,
           breaks = breaks,
           xlab = "Factor", ylab = "Sample")
  pf <- pheatmap(F,
           main = paste(lf_str," Factors"),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           silent= TRUE,
           breaks = breaks,
           xlab = "Feature", ylab = "Factor")
  grid.arrange(pl$gtable, pf$gtable, ncol=2)
}

plot_em_ll <- function(init_ll, ll_list, gt_ll){
  # Plot Log likelihood over EM iterations
  y_bounds <- range(c(init_ll, ll_list, gt_ll), na.rm = TRUE)
  plot(1:length(ll_list), ll_list,
       xlim = c(0, length(ll_list)),
       ylim = y_bounds,
       main = sprintf("Log Likelihood over EM Iterations\nTrue LL (blue line): %.2f", gt_ll),
       ylab = "Log-Likelihood",
       xlab = "EM Iteration",
       pch = 16)
  points(0, init_ll, col = "red", pch = 16, cex = 1.5) # Initial LL point at x = 0
  abline(h = gt_ll, col = "blue", lty = 2, lwd = 2) # True LL
}
