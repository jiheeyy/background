# Functions useful after running ebpmf-bg-cavi

# Plot the initial parameter values.
initial_heatmap <- function(
    value,
    main,
    row_label,
    row_axis_at = NULL,
    row_axis_labels = NULL
) {
  color_limit <- range(value, finite = TRUE)
  if (diff(color_limit) == 0) {
    padding <- max(abs(color_limit[1L]), 1) * 1e-8
    color_limit <- color_limit + c(-padding, padding)
  }

  image(
    x = seq_len(ncol(value)),
    y = seq_len(nrow(value)),
    z = t(value),
    col = hcl.colors(80L, "Viridis"),
    zlim = color_limit,
    xlab = "Factor",
    ylab = row_label,
    main = "",
    axes = FALSE,
    useRaster = TRUE
  )
  axis(1, at = seq_len(ncol(value)), labels = seq_len(ncol(value)))

  if (is.null(row_axis_at)) {
    row_axis_at <- unique(c(1L, nrow(value)))
    row_axis_labels <- row_axis_at
  }
  axis(
    2,
    at = row_axis_at,
    labels = row_axis_labels,
    las = 1,
    cex.axis = 0.8
  )
  box()
  # Reserve separate margin lines for the heading and value range.
  title(main = main, line = 2)
  mtext(
    sprintf("range: %.3g to %.3g", min(value), max(value)),
    side = 3,
    line = 0.2,
    cex = 0.75
  )
}

column_cosine_similarity <- function(A, B) {
  denominator <- outer(
    sqrt(colSums(A^2)),
    sqrt(colSums(B^2))
  )
  crossprod(A, B) / pmax(denominator, .Machine$double.eps)
}

solve_maximum_assignment <- function(score) {
  K <- nrow(score)
  stopifnot(ncol(score) == K, K >= 1L, K <= 20L, all(is.finite(score)))

  if (K == 1L) return(1L)

  n_state <- bitwShiftL(1L, K)
  best <- rep(-Inf, n_state)
  previous <- chosen <- rep(NA_integer_, n_state)
  used_count <- vapply(
    0:(n_state - 1L),
    function(mask) sum(as.integer(intToBits(mask)[seq_len(K)])),
    integer(1)
  )
  best[1] <- 0

  for (mask in 0:(n_state - 2L)) {
    topic <- used_count[mask + 1L] + 1L
    for (candidate in seq_len(K)) {
      bit <- bitwShiftL(1L, candidate - 1L)
      if (bitwAnd(mask, bit) != 0L) {
        next
      }

      next_mask <- bitwOr(mask, bit)
      value <- best[mask + 1L] + score[topic, candidate]
      if (value > best[next_mask + 1L]) {
        best[next_mask + 1L] <- value
        previous[next_mask + 1L] <- mask
        chosen[next_mask + 1L] <- candidate
      }
    }
  }

  assignment <- integer(K)
  mask <- n_state - 1L
  for (topic in K:1L) {
    assignment[topic] <- chosen[mask + 1L]
    mask <- previous[mask + 1L]
  }
  assignment
}

# Shared matching core for simulated truth and between-model comparisons.
# "rate" compares rank-one rate contributions; "profiles" averages two cosines.
match_factor_order <- function(L_reference, F_reference, L_fitted, F_fitted,
                               criterion = c("rate", "profiles")) {
  criterion <- match.arg(criterion)
  stopifnot(
    ncol(L_reference) == ncol(F_reference),
    ncol(L_fitted) == ncol(F_fitted),
    ncol(L_reference) == ncol(L_fitted)
  )
  document_similarity <- column_cosine_similarity(L_reference, L_fitted)
  word_similarity <- column_cosine_similarity(F_reference, F_fitted)
  score <- if (criterion == "rate") {
    document_similarity * word_similarity
  } else {
    (document_similarity + word_similarity) / 2
  }
  list(
    permutation = solve_maximum_assignment(score),
    document_similarity = document_similarity,
    word_similarity = word_similarity
  )
}

# Align factor order to the simulated truth without rescaling.
align_k <- function(synthetic, state) {
  L_true <- synthetic$L
  B_true <- sweep(synthetic$F, 1L, synthetic$f0, "*")
  L_hat <- state$L
  F_hat <- state$qF$mean
  B_hat <- sweep(F_hat, 1L, state$f0, "*")
  K <- ncol(L_true)

  stopifnot(
    ncol(L_hat) == K,
    ncol(B_true) == K,
    ncol(B_hat) == K,
    ncol(state$pi) == K
  )

  order <- match_factor_order(L_true, B_true, L_hat, B_hat)$permutation

  list(
    L = L_hat[, order, drop = FALSE],
    F = F_hat[, order, drop = FALSE],
    B = B_hat[, order, drop = FALSE],
    pi = state$pi[, order, drop = FALSE],
    permutation = order
  )
}

# Rescale L columns to mean one while preserving L %*% t(B).
canonicalize <- function(L, B) {
  component_scale <- colMeans(L)

  list(
    L = sweep(L, 2, component_scale, "/"),
    B = sweep(B, 2, component_scale, "*")
  )
}

# Compare true and estimated values against the identity line.
plot_comparison <- function(truth, estimate, main, xlab, ylab) {
  limits <- range(c(truth, estimate), finite = TRUE)

  plot(
    truth, estimate,
    pch = 19,
    col = rgb(0.1, 0.35, 0.75, 0.55),
    xlim = limits,
    ylim = limits,
    xlab = xlab,
    ylab = ylab,
    main = main
  )

  abline(0, 1, col = "firebrick", lwd = 2, lty = 2)

  legend(
    "topleft",
    legend = sprintf("cor = %.3f", cor(truth, estimate)),
    bty = "n"
  )
}

# Draw a wide heatmap on the log1p scale.
draw_heatmap_wide <- function(
    x,
    main,
    zlim,
    row_labels = NULL,
    xlab = "Index"
) {
  x_plot <- log1p(x)

  image(
    x = seq_len(ncol(x_plot)),
    y = seq_len(nrow(x_plot)),
    z = t(x_plot),
    zlim = zlim,
    col = hcl.colors(100, "Viridis"),
    useRaster = TRUE,
    axes = FALSE,
    xlab = xlab,
    ylab = "",
    main = main
  )

  # Show a few horizontal index labels.
  x_ticks <- unique(round(seq(
    1,
    ncol(x_plot),
    length.out = min(6, ncol(x_plot))
  )))

  axis(1, at = x_ticks, labels = x_ticks)

  if (is.null(row_labels)) {
    row_labels <- seq_len(nrow(x_plot))
  }

  axis(
    2,
    at = seq_len(nrow(x_plot)),
    labels = row_labels,
    las = 1
  )

  box()
}

# Draw true and fitted heatmaps with a shared color scale.
draw_pair_wide <- function(
    truth,
    estimate,
    label,
    row_labels = NULL,
    xlab = "Index"
) {
  pair_range <- range(
    c(log1p(truth), log1p(estimate)),
    finite = TRUE
  )

  if (diff(pair_range) == 0) {
    pair_range <- pair_range + c(-1, 1) * 1e-8
  }

  draw_heatmap_wide(
    truth,
    main = paste("True", label, ", Max value ", round(max(truth),1)),
    zlim = pair_range,
    row_labels = row_labels,
    xlab = xlab
  )

  draw_heatmap_wide(
    estimate,
    main = paste("Estimated", label, ", Max value ", round(max(estimate),1)),
    zlim = pair_range,
    row_labels = row_labels,
    xlab = xlab
  )
}

# Display mixture proportions with numeric labels.
draw_mixture_proportions <- function(pi, main, a = NULL) {
  M <- nrow(pi)
  K <- ncol(pi)

  image(
    x = seq_len(K),
    y = seq_len(M),
    z = t(pi),
    zlim = c(0, 1),
    col = hcl.colors(100, "Blues 3", rev = TRUE),
    axes = FALSE,
    xlab = "Factor",
    ylab = "",
    main = main,
    useRaster = TRUE
  )

  axis(
    1,
    at = seq_len(K),
    labels = paste0("k", seq_len(K))
  )

  mixture_labels <- paste0("m", seq_len(M))
  if (!is.null(a)) {
    mixture_labels <- sprintf(
      "m%d (a = %g)",
      seq_len(M),
      a
    )
  }

  axis(
    2,
    at = seq_len(M),
    labels = mixture_labels,
    las = 1
  )

  cells <- expand.grid(
    factor = seq_len(K),
    mixture = seq_len(M)
  )
  values <- as.vector(t(pi))

  text(
    x = cells$factor,
    y = cells$mixture,
    labels = sprintf("%.3f", values),
    col = ifelse(values >= 0.55, "white", "black"),
    font = 2,
    cex = 1.15
  )

  box()
}

# Plot the state entering CAVI.
plot_initial_parameters <- function(state, a) {
  L_initial <- state$L
  F_initial <- state$qF$mean
  pi_initial <- state$pi
  f0_initial <- state$f0
  old_par <- par(mfrow = c(2, 2), mar = c(2.5, 3, 3.5, 0.5))
  on.exit(par(old_par), add = TRUE)

  initial_heatmap(
    L_initial,
    main = "Initial L",
    row_label = "Document"
  )

  initial_heatmap(
    F_initial,
    main = "Initial E_q[F]",
    row_label = "Feature"
  )

  initial_heatmap(
    pi_initial,
    main = "Initial pi",
    row_label = "Gamma component",
    row_axis_at = seq_len(nrow(pi_initial)),
    row_axis_labels = paste0("a=", a)
  )

  plot(
    seq_along(f0_initial),
    f0_initial,
    type = "l",
    lwd = 1.5,
    col = "steelblue",
    xlab = "Feature",
    ylab = "f0",
    main = "Initial f0"
  )
  abline(h = 1, col = "grey50", lty = 2)
  grid(col = "grey90")

}

# Compare aligned, canonicalized loadings and effective feature factors.
plot_factor_recovery <- function(truth, estimate) {
  K <- ncol(truth$L)
  L_true_plot <- truth$L
  L_hat_plot <- estimate$L
  B_true_plot <- truth$B
  B_hat_plot <- estimate$B
  old_par <- par(mfrow = c(2, K), mar = c(4, 4, 3, 1))
  on.exit(par(old_par), add = TRUE)

  for (k in seq_len(K)) {
    plot_comparison(
      L_true_plot[, k],
      L_hat_plot[, k],
      main = sprintf("L, component %d", k),
      xlab = "True normalized L",
      ylab = "Estimated normalized L"
    )
  }

  for (k in seq_len(K)) {
    plot_comparison(
      B_true_plot[, k],
      B_hat_plot[, k],
      main = sprintf("f0 * F, component %d", k),
      xlab = "True effective F",
      ylab = "Estimated effective F"
    )
  }

}

# Plot mixture trajectories with generating values and stabilization times.
plot_mixture_trace <- function(trace, pi, a, stabilization) {
  K <- dim(trace)[3L]
  iterations <- seq_len(dim(trace)[1L]) - 1L
  mixture_colors <- setNames(
    hcl.colors(length(a), "Dark 3"),
    paste0("a = ", a)
  )

  old_par <- par(
    mfrow = c(1, K),
    mar = c(4, 4, 3, 1)
  )
  on.exit(par(old_par), add = TRUE)

  for (k in seq_len(K)) {
    matplot(
      iterations,
      trace[, , k],
      type = "l",
      lty = 1,
      lwd = 2,
      col = mixture_colors,
      ylim = c(0, 1),
      xlab = "Iteration",
      ylab = "Mixture proportion",
      main = paste("True factor", k)
    )

    # Generating mixture proportions.
    for (m in seq_len(length(a))) {
      abline(
        h = pi[m, k],
        col = mixture_colors[m],
        lty = 2
      )
    }

    # Empirical stabilization iteration.
    if (!is.na(stabilization[k])) {
      abline(
        v = stabilization[k],
        col = "grey40",
        lty = 3
      )
    }

    legend(
      "right",
      legend = names(mixture_colors),
      col = mixture_colors,
      lty = 1,
      lwd = 2,
      cex = 0.85,
      bty = "n"
    )
  }

}

# Plot the log10 gap from the best ELBO, omitting the final checkpoint.
plot_elbo_gap <- function(elbo_trace) {
  best_elbo <- max(elbo_trace$elbo)
  elbo_by_iteration <- head(elbo_trace[, c("iteration", "elbo")], -1L)
  elbo_by_iteration$log10_gap <- log10(
    best_elbo - elbo_by_iteration$elbo
  )

  plot(
    log10_gap ~ iteration,
    data = elbo_by_iteration,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "steelblue",
    xlab = "CAVI iteration",
    ylab = expression(log[10](max(ELBO) - ELBO)),
    main = "ELBO convergence"
  )
  grid(col = "grey90")

  # Label the final visible checkpoint.
  last_point <- tail(elbo_by_iteration, 1L)
  text(
    last_point$iteration,
    last_point$log10_gap,
    labels = paste("Iteration", last_point$iteration),
    pos = 2
  )
}

# First iteration after which each factor stays within tolerance of its final pi.
# trace has dimensions iteration x mixture x factor; iterations start at zero.
mixture_stabilization <- function(trace, tolerance = 0.01) {
  vapply(seq_len(dim(trace)[3L]), function(k) {
    trajectory <- matrix(trace[, , k], nrow = dim(trace)[1L])
    final <- trajectory[nrow(trajectory), ]
    distance <- apply(abs(sweep(trajectory, 2L, final, "-")), 1L, max)
    future_max <- rev(cummax(rev(distance)))
    which(future_max <= tolerance)[1L] - 1L
  }, integer(1))
}

# Plot aligned mixture trajectories and return alignment and convergence summaries.
plot_mixprop_convergence <- function(pi_trace, pi_true, pi_hat, permutation, a, tolerance = 0.01) {
  K <- ncol(pi_true)
  niter <- dim(pi_trace)[1L] - 1L

  # Apply the final factor order to every iteration.
  pi_trace_aligned <- pi_trace[, , permutation, drop = FALSE]

  dimnames(pi_trace_aligned)$estimated_factor <-
    paste0("k", seq_len(K))

  # Confirm that the final trace values equal the aligned pi used
  # by the mixture-proportion heatmap.
  stopifnot(
    isTRUE(all.equal(
      unname(pi_trace_aligned[niter + 1L, , ]),
      unname(pi_hat),
      tolerance = 1e-10
    ))
  )

  # First iteration after which all proportions remain within the
  # tolerance of their final values.
  stabilization_iteration <- mixture_stabilization(
    pi_trace_aligned, tolerance = tolerance
  )

  plot_mixture_trace(
    pi_trace_aligned, pi_true, a, stabilization_iteration
  )

  final_mixprop <- expand.grid(
    mixture = paste0("a=", a),
    factor = paste0("k", seq_len(K))
  )

  final_mixprop$generating <- as.vector(pi_true)
  final_mixprop$estimated <- as.vector(pi_hat)

  list(
    alignment = data.frame(
      true_factor = paste0("k", seq_len(K)),
      estimated_column = permutation
    ),
    stabilization = data.frame(
      factor = paste0("k", seq_len(K)),
      iteration = stabilization_iteration
    ),
    final_proportions = final_mixprop
  )
}

# Plot reconstructed Poisson means and return recovery metrics.
plot_lambda <- function(L_true, B_true, L_hat, B_hat) {
  lambda_true <- L_true %*% t(B_true)
  lambda_hat  <- L_hat  %*% t(B_hat)

  lambda_limits <- range(c(lambda_true, lambda_hat))

  plot(
    as.vector(lambda_true),
    as.vector(lambda_hat),
    pch = 19,
    cex = 0.7,
    col = rgb(0, 0, 0, 0.25),
    xlim = lambda_limits,
    ylim = lambda_limits,
    xlab = "True Poisson mean",
    ylab = "Estimated Poisson mean",
    main = "Reconstructed mean: true versus fitted"
  )
  abline(0, 1, col = "firebrick", lwd = 2)

  lambda_metrics <- data.frame(
    correlation = cor(as.vector(lambda_true), as.vector(lambda_hat)),
    relative_error = sqrt(sum((lambda_hat - lambda_true)^2)) /
      sqrt(sum(lambda_true^2)),
    rmse = sqrt(mean((lambda_hat - lambda_true)^2))
  )

  lambda_metrics
}

# Plot raw background recovery and return error metrics.
plot_bgquality <- function(L_true, L_hat, f0_true, f0_hat) {
  loading_scale <- c(
    true = mean(rowSums(L_true)),
    estimated = mean(rowSums(L_hat))
  )

  if (any(abs(loading_scale - 1) > 1e-2)) {
    warning(sprintf(
      "Loading scales are not both within 1e-2 of 1 (true = %.4f, estimated = %.4f); raw f0 values may not be directly comparable.",
      loading_scale["true"], loading_scale["estimated"]
    ))
  }

  log_f0_true <- log(f0_true)
  log_f0_hat <- log(f0_hat)
  plot_limits <- range(log_f0_true, log_f0_hat)
  feature <- seq_along(f0_true)

  old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  on.exit(par(old_par), add = TRUE)

  matplot(
    feature,
    cbind(log_f0_true, log_f0_hat),
    type = "l",
    lty = c(1, 2),
    lwd = 2,
    col = c("black", "steelblue"),
    xlab = "Feature",
    ylab = "log(f0)",
    main = "Background profile"
  )
  legend(
    "topright",
    legend = c("True", "Estimated"),
    col = c("black", "steelblue"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )

  plot(
    log_f0_true,
    log_f0_hat,
    pch = 19,
    col = rgb(0.1, 0.35, 0.75, 0.55),
    xlim = plot_limits,
    ylim = plot_limits,
    xlab = "True log(f0)",
    ylab = "Estimated log(f0)",
    main = "Background recovery"
  )
  abline(0, 1, col = "firebrick", lwd = 2, lty = 2)


  log_error <- log_f0_hat - log_f0_true
  background_metrics <- data.frame(
    correlation = cor(f0_true, f0_hat),
    log_correlation = cor(log_f0_true, log_f0_hat),
    log_rmse = sqrt(mean(log_error^2)),
    median_fold_error = exp(median(abs(log_error)))
  )

  background_metrics
}

# Plot canonicalized loadings and effective feature factors.
# Heatmaps use log1p values because L and f0*F likely contain most entries
# close to zero and a few much larger values. The effective factor f0*F is
# emphasized because it is the feature-side quantity that enters the Poisson
# rate. However, separate differences in f0 and F can compensate for one another.
plot_lf <- function(truth, estimate) {
  K <- ncol(truth$L)
  L_true_plot <- truth$L
  L_hat_plot <- estimate$L
  B_true_plot <- truth$B
  B_hat_plot <- estimate$B

  factor_labels <- paste0("k", seq_len(K))

  old_par <- par(
    mfrow = c(3, 2),
    mar = c(2.5, 3, 2, 0.5),
    mgp = c(2.4, 0.7, 0),
    tcl = -0.25,
    cex.axis = 0.9,
    cex.lab = 1,
    cex.main = 1.05
  )
  on.exit(par(old_par), add = TRUE)

  # K rows × n columns
  draw_pair_wide(
    t(L_true_plot),
    t(L_hat_plot),
    label = "L",
    row_labels = factor_labels,
    xlab = "Sample"
  )

  # K rows × p columns
  draw_pair_wide(
    t(B_true_plot),
    t(B_hat_plot),
    label = "f0 * F",
    row_labels = factor_labels,
    xlab = "Feature"
  )
}

# Compare true and aligned estimated mixture proportions.
plot_mixprop <- function(pi_true, pi_hat, a) {
  old_par <- par(
    mfrow = c(1, 2),
    mar = c(3, 4, 2, 0.5),
    mgp = c(2.4, 0.7, 0),
    tcl = -0.25
  )
  on.exit(par(old_par), add = TRUE)

  draw_mixture_proportions(
    pi_true,
    main = "True mixture proportions",
    a = a
  )

  draw_mixture_proportions(
    pi_hat,
    main = "Estimated mixture proportions",
    a = a
  )
}

# Plot background/factor error compensation and return its correlation.
plot_compensate <- function(B_true, B_hat, f0_true, f0_hat) {
  # If correlation high, error in f0 is offsetting error in F
  F_true_canonical <- sweep(B_true, 1, f0_true, "/")
  F_hat_canonical  <- sweep(B_hat, 1, f0_hat, "/")

  log_f0_error <- log(f0_hat / f0_true)
  mean_log_F_error <- rowMeans(log(F_hat_canonical / F_true_canonical))

  plot(
    log_f0_error,
    -mean_log_F_error,
    pch = 19,
    xlab = "log(f0_hat / f0_true)",
    ylab = "-mean log(F_hat / F_true)",
    main = "Compensation between f0 and F"
  )
  abline(0, 1, col = "firebrick", lwd = 2)

  cor(log_f0_error, -mean_log_F_error)
}

# ---- Newsgroups data, fitting, and comparison diagnostics ----

# pi is stored as mixture components by topics; display topics as rows.
plot_news_mixprop <- function(pi, a) {
  stopifnot(is.matrix(pi), nrow(pi) == length(a), ncol(pi) > 0L,
            all(is.finite(pi)), all(pi >= 0 & pi <= 1))
  cells <- expand.grid(component = seq_along(a), topic = seq_len(ncol(pi)))
  cells$pi <- as.vector(pi)
  ggplot2::ggplot(cells, ggplot2::aes(component, topic, fill = pi)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_x_continuous(breaks = seq_along(a), labels = a,
                              expand = c(0, 0)) +
    ggplot2::scale_y_reverse(breaks = seq_len(ncol(pi)), expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), name = expression(pi)) +
    ggplot2::labs(x = "Mixture component a", y = "Topic k",
                  title = "Final EBPMF-BG mixture proportions") +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# Load the complete sparse count matrix and observed newsgroup labels.
load_news_data <- function(data_path, map_path) {
  # Load only the required objects without adding them to the working environment.
  newsgroups_data <- local({
    data_env <- new.env(parent = emptyenv())
    load(data_path, envir = data_env)
    mget(c("counts", "topics"), envir = data_env, inherits = FALSE)
  })
  X <- newsgroups_data$counts
  y_numeric <- newsgroups_data$topics

  # Match numeric topic IDs to the original 20-newsgroups names.
  # Note train.map is already ordered 1 through 20
  news_map <- read.table(
    map_path,
    col.names = c("label", "id"),
    stringsAsFactors = FALSE
  )
  stopifnot(identical(news_map$id, seq_len(nrow(news_map))))
  news_order <- c(
    "rec.sport.hockey",
    "rec.sport.baseball",
    "sci.med",
    "comp.graphics",
    "comp.windows.x",
    "comp.os.ms-windows.misc",
    "comp.sys.ibm.pc.hardware",
    "comp.sys.mac.hardware",
    "misc.forsale",
    "sci.electronics",
    "sci.space",
    "alt.atheism",
    "soc.religion.christian",
    "talk.religion.misc",
    "rec.autos",
    "rec.motorcycles",
    "sci.crypt",
    "talk.politics.misc",
    "talk.politics.guns",
    "talk.politics.mideast"
  )
  y_str <- factor(
    news_map$label[y_numeric],
    levels = news_order
  )
  stopifnot(nrow(X) == length(y_str), !anyNA(y_str))
  list(X = X, labels = y_str, terms = colnames(X))
}

# Optional reproducible document sampling and selection of the most frequent terms.
sample_news_data <- function(data, max_documents = 1000L, max_terms = 2000L, seed = 1L) {
  X <- data$X
  set.seed(seed)
  rows <- sample(seq_len(nrow(X)), min(max_documents, nrow(X)))
  totals <- as.numeric(Matrix::colSums(X))
  columns <- head(order(-totals, seq_along(totals), method = "radix"),
                  min(max_terms, ncol(X)))
  X <- X[rows, columns, drop = FALSE]
  X <- X[, Matrix::colSums(X) > 0, drop = FALSE]
  list(X = X, labels = droplevels(data$labels[rows]), terms = colnames(X))
}

normalize_rows <- function(X) {
  sweep(X, 1, pmax(rowSums(X), .Machine$double.eps), "/")
}

normalize_columns <- function(X) {
  sweep(X, 2, pmax(colSums(X), .Machine$double.eps), "/")
}

# Align fastTopics to EBPMF using document and word profiles; retain one shared mapping.
align_news_models <- function(fits) {
  K_news <- ncol(fits$ebpmf$state$L)
  # Raw factors reconstruct the Poisson rates.
  L_ebpmf_rate <- fits$ebpmf$state$L
  F_ebpmf_rate <- sweep(fits$ebpmf$state$qF$mean, 1L, fits$ebpmf$state$f0, "*")
  L_fast_rate <- fits$fasttopics$L
  F_fast_rate <- fits$fasttopics$F

  # Multinomial forms are used for document proportions and word frequencies.
  # Weight each loading by its topic's total word mass before row normalization.
  # This preserves the fitted rates when converting to multinomial proportions.
  L_ebpmf_topics <- normalize_rows(
    sweep(L_ebpmf_rate, 2L, colSums(F_ebpmf_rate), "*")
  )
  F_ebpmf_topics <- normalize_columns(F_ebpmf_rate)
  fasttopics_topics <- poisson2multinom(fits$fasttopics)
  L_fast_topics <- normalize_rows(fasttopics_topics$L)
  F_fast_topics <- normalize_columns(fasttopics_topics$F)

  matching <- match_factor_order(
    L_ebpmf_topics, F_ebpmf_topics, L_fast_topics, F_fast_topics,
    criterion = "profiles"
  )
  document_similarity <- matching$document_similarity
  word_similarity <- matching$word_similarity
  fast_topic_order <- matching$permutation

  L_fast_rate <- L_fast_rate[, fast_topic_order, drop = FALSE]
  F_fast_rate <- F_fast_rate[, fast_topic_order, drop = FALSE]
  L_fast_topics <- L_fast_topics[, fast_topic_order, drop = FALSE]
  F_fast_topics <- F_fast_topics[, fast_topic_order, drop = FALSE]

  factor_names <- paste0("k", seq_len(K_news))
  colnames(L_ebpmf_topics) <- colnames(F_ebpmf_topics) <- factor_names
  colnames(L_fast_topics) <- colnames(F_fast_topics) <- factor_names

  alignment_table <- data.frame(
    EBPMF_topic = factor_names,
    original_fastTopics_topic = paste0("k", fast_topic_order),
    document_cosine = document_similarity[
      cbind(seq_len(K_news), fast_topic_order)
    ],
    word_cosine = word_similarity[
      cbind(seq_len(K_news), fast_topic_order)
    ]
  )

  list(
    L_ebpmf_rate = L_ebpmf_rate,
    F_ebpmf_rate = F_ebpmf_rate,
    L_fast_rate = L_fast_rate,
    F_fast_rate = F_fast_rate,
    L_ebpmf_topics = L_ebpmf_topics,
    F_ebpmf_topics = F_ebpmf_topics,
    L_fast_topics = L_fast_topics,
    F_fast_topics = F_fast_topics,
    factor_names = factor_names,
    alignment_table = alignment_table,
    fast_topic_order = fast_topic_order
  )
}

# Evaluate all counts without materializing full dense count or rate matrices.
news_fit_metrics <- function(X, aligned, block_size = 64L) {
  # Only small document-by-word blocks are dense, even for the full corpus.
  stopifnot(block_size >= 1L)
  X <- methods::as(X, "dgCMatrix")
  total <- sum(X)
  squared_total <- sum(X@x^2)
  log_factorial <- sum(lgamma(X@x + 1))
  saturated <- sum(X@x * log(X@x) - X@x - lgamma(X@x + 1), na.rm = TRUE)
  row_totals <- Matrix::rowSums(X)
  col_totals <- Matrix::colSums(X)
  positive_rows <- row_totals > 0
  positive_cols <- col_totals > 0
  null_loglik <- sum(row_totals[positive_rows] * log(row_totals[positive_rows])) +
    sum(col_totals[positive_cols] * log(col_totals[positive_cols])) -
    total * log(total) - total - log_factorial
  null_deviance <- 2 * (saturated - null_loglik)
  summarize <- function(method, L, F) {
    loglik <- -log_factorial
    squared_error <- 0
    for (start in seq.int(1L, ncol(X), by = block_size)) {
      columns <- seq.int(start, min(start + block_size - 1L, ncol(X)))
      observed <- as.matrix(X[, columns, drop = FALSE])
      rate <- tcrossprod(L, F[columns, , drop = FALSE])
      squared_error <- squared_error + sum((observed - rate)^2)
      rate <- pmax(rate, 1e-12)
      loglik <- loglik + sum(observed * log(rate) - rate)
    }
    deviance <- 2 * (saturated - loglik)
    data.frame(
      method = method, log_likelihood = loglik,
      log_likelihood_per_token = loglik / total,
      mean_deviance = deviance / (as.double(nrow(X)) * ncol(X)),
      explained_deviance = 1 - deviance / null_deviance,
      relative_rmse = sqrt(squared_error / squared_total)
    )
  }
  rbind(
    summarize("EBPMF", aligned$L_ebpmf_rate, aligned$F_ebpmf_rate),
    summarize("fastTopics", aligned$L_fast_rate, aligned$F_fast_rate)
  )
}

# Round the displayed fit metrics without altering the underlying results.
format_news_metrics <- function(metrics) {
  transform(
    metrics, log_likelihood = round(log_likelihood),
    log_likelihood_per_token = round(log_likelihood_per_token, 4),
    mean_deviance = round(mean_deviance, 4),
    explained_deviance = round(explained_deviance, 4),
    relative_rmse = round(relative_rmse, 4)
  )
}

# Newsgroups diagnostic: fit-comparison.
plot_news_fit <- function(fit_metrics, X) {
  method_order <- c(
    "EBPMF",
    "fastTopics"
  )
  fit_comparison_metrics <- fit_metrics
  fit_comparison_metrics$method <- factor(
    fit_comparison_metrics$method,
    levels = method_order
  )

  fit_plot_data <- rbind(
    data.frame(
      method = fit_comparison_metrics$method,
      metric = "Explained deviance\n(higher is better)",
      value = fit_comparison_metrics$explained_deviance
    ),
    data.frame(
      method = fit_comparison_metrics$method,
      metric = "Relative RMSE\n(lower is better)",
      value = fit_comparison_metrics$relative_rmse
    )
  )

  print(ggplot(fit_plot_data, aes(method, value, fill = method)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(
      aes(label = sprintf("%.3f", value)),
      vjust = -0.4,
      size = 3.5
    ) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values = c("steelblue", "darkorange")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(
      x = NULL,
      y = NULL,
      title = sprintf(
        "Poisson fit on the same %s documents and %s terms",
        format(nrow(X), big.mark = ","),
        format(ncol(X), big.mark = ",")
      )
    ) +
    theme_cowplot(font_size = 11) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain")
    ))

  transform(
    fit_comparison_metrics,
    log_likelihood = round(log_likelihood),
    log_likelihood_per_token = round(log_likelihood_per_token, 4),
    mean_deviance = round(mean_deviance, 4),
    explained_deviance = round(explained_deviance, 4),
    relative_rmse = round(relative_rmse, 4)
  )
}

# Newsgroups diagnostic: convergence.
plot_news_convergence <- function(fits) {
  fasttopics_fit <- fits$fasttopics
  ebpmf_elbo <- fits$ebpmf$elbo_trace
  names(ebpmf_elbo)[1:2] <- c("iter", "total")

  fast_progress <- fasttopics_fit$progress
  fast_progress$loglik_gap <- max(fast_progress$loglik) -
    fast_progress$loglik
  positive_gap <- fast_progress$loglik_gap[fast_progress$loglik_gap > 0]
  gap_floor <- if (length(positive_gap)) min(positive_gap) / 10 else 1e-8
  fast_progress$loglik_gap <- pmax(fast_progress$loglik_gap, gap_floor)

  p_fast_loglik <- ggplot(
    fast_progress,
    aes(iter, loglik_gap)
  ) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 1.2) +
    scale_y_log10() +
    labs(
      x = "Iteration",
      y = "Gap to best log-likelihood",
      title = "fastTopics likelihood convergence"
    ) +
    theme_cowplot(font_size = 11)

  p_fast_residual <- ggplot(fast_progress, aes(iter, res)) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    geom_point(color = "firebrick", size = 1.2) +
    scale_y_log10() +
    labs(
      x = "Iteration",
      y = "KKT residual",
      title = "fastTopics stationarity diagnostic"
    ) +
    theme_cowplot(font_size = 11)

  p_ebpmf_elbo <- ggplot(ebpmf_elbo, aes(iter, total)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 2) +
    labs(
      x = "Iteration",
      y = "ELBO",
      title = "EBPMF objective convergence"
    ) +
    theme_cowplot(font_size = 11)

  ebpmf_change <- ebpmf_elbo[
    is.finite(ebpmf_elbo$relative_change) & ebpmf_elbo$relative_change > 0,
  ]

  p_ebpmf_change <- ggplot(
    ebpmf_change,
    aes(iter, relative_change)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_y_log10() +
    labs(
      x = "Iteration",
      y = "Relative ELBO change",
      title = "EBPMF checkpoint convergence"
    ) +
    theme_cowplot(font_size = 11) +
    theme(legend.position = "bottom")

  plot_grid(
    p_fast_loglik,
    p_fast_residual,
    p_ebpmf_elbo,
    p_ebpmf_change,
    ncol = 2
  )
}

# Plot a shared document sample or all documents, with independent t-SNE ordering.
plot_news_structure <- function(aligned, labels, n = 1000L, seed = 1L,
                                embedding_seed = 1L, plot_sample = TRUE) {
  profiles <- if (is.matrix(aligned) || inherits(aligned, "Matrix")) {
    list("EBPMF initialization" = normalize_rows(as.matrix(aligned)))
  } else {
    list(EBPMF = aligned$L_ebpmf_topics, fastTopics = aligned$L_fast_topics)
  }
  n_documents <- nrow(profiles[[1L]])
  stopifnot(is.logical(plot_sample), length(plot_sample) == 1L,
            !is.na(plot_sample))
  stopifnot(length(n) == 1L, n >= 1, length(labels) == n_documents,
            all(vapply(profiles, nrow, integer(1)) == n_documents))
  # Restore the caller's RNG state after selecting the same rows for each panel.
  rows <- if (plot_sample && n < n_documents) {
    withr::with_seed(seed, sample.int(n_documents, n))
  } else {
    seq_len(n_documents)
  }
  profiles <- lapply(profiles, function(L) L[rows, , drop = FALSE])
  K_news <- ncol(profiles[[1L]])
  comparison_labels <- droplevels(labels[rows])
  # structure_plot embeds each newsgroup separately for each model.
  embed_group <- function(fit, ...) {
    # t-SNE needs at least four documents; retain input order for smaller groups.
    if (nrow(fit$L) < 4L) return(seq_len(nrow(fit$L)))
    drop(fastTopics::tsne_from_topics(fit, dims = 1L, verbose = FALSE))
  }

  topic_colors <- c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
    "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"
  )

  plots <- lapply(names(profiles), function(name) {
    # Reset the embedding seed for each panel without changing the caller's RNG.
    withr::with_seed(embedding_seed, structure_plot(
      profiles[[name]], topics = seq_len(K_news),
      grouping = comparison_labels, loadings_order = "embed",
      n = nrow(profiles[[name]]), embed_method = embed_group,
      colors = topic_colors, gap = 6
    )) +
      ggtitle(name) +
      theme_cowplot(font_size = 10) +
      theme(
        plot.title = element_text(face = "plain"),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 7)
      )
  })
  if (length(plots) == 1L) return(plots[[1L]])
  plot_grid(plotlist = plots, ncol = 1, align = "v")
}

# Newsgroups diagnostic: factor-agreement.
plot_news_agreement <- function(aligned, terms) {
  L_ebpmf_topics <- aligned$L_ebpmf_topics
  F_ebpmf_topics <- aligned$F_ebpmf_topics
  L_fast_topics <- aligned$L_fast_topics
  F_fast_topics <- aligned$F_fast_topics
  factor_names <- aligned$factor_names
  K_news <- ncol(aligned$L_ebpmf_topics)
  comparison_terms <- terms
  document_plot_data <- expand.grid(
    document = seq_len(nrow(L_ebpmf_topics)),
    factor = factor_names
  )
  document_plot_data$EBPMF <- as.vector(L_ebpmf_topics)
  document_plot_data$fastTopics <- as.vector(L_fast_topics)

  word_plot_data <- expand.grid(
    term = comparison_terms,
    factor = factor_names
  )
  word_plot_floor <- 1e-6
  word_plot_data$EBPMF <- pmax(
    as.vector(F_ebpmf_topics),
    word_plot_floor
  )
  word_plot_data$fastTopics <- pmax(
    as.vector(F_fast_topics),
    word_plot_floor
  )
  word_frequency_limits <- range(
    unlist(
      word_plot_data[c("EBPMF", "fastTopics")],
      use.names = FALSE
    )
  )

  p_document_agreement <- ggplot(
    document_plot_data,
    aes(EBPMF, fastTopics)
  ) +
    geom_point(color = "steelblue", alpha = 0.25, size = 0.6) +
    geom_abline(linetype = 2, color = "firebrick") +
    facet_wrap(~ factor, ncol = 5) +
    coord_equal() +
    labs(
      x = "EBPMF document-topic proportion",
      y = "fastTopics document-topic proportion",
      title = "Aligned document loadings"
    ) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  p_word_agreement <- ggplot(
    word_plot_data,
    aes(EBPMF, fastTopics)
  ) +
    geom_point(color = "steelblue", alpha = 0.25, size = 0.6) +
    geom_abline(linetype = 2, color = "firebrick") +
    scale_x_log10(limits = word_frequency_limits) +
    scale_y_log10(limits = word_frequency_limits) +
    facet_wrap(~ factor, ncol = 5) +
    coord_equal() +
    labs(
      x = "EBPMF word frequency",
      y = "fastTopics word frequency",
      title = "Aligned topic-word frequencies"
    ) +
    theme_cowplot(font_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  print(plot_grid(
    p_document_agreement,
    p_word_agreement,
    ncol = 1
  ))

  factor_agreement <- data.frame(
    factor = factor_names,
    document_correlation = vapply(
      seq_len(K_news),
      function(k) cor(L_ebpmf_topics[, k], L_fast_topics[, k]),
      numeric(1)
    ),
    word_correlation = vapply(
      seq_len(K_news),
      function(k) cor(F_ebpmf_topics[, k], F_fast_topics[, k]),
      numeric(1)
    )
  )

  transform(
    factor_agreement,
    document_correlation = round(document_correlation, 3),
    word_correlation = round(word_correlation, 3)
  )
}

# Newsgroups diagnostic: keyword-topic.
plot_news_keyword <- function(aligned, X, terms) {
  F_ebpmf_topics <- aligned$F_ebpmf_topics
  F_fast_topics <- aligned$F_fast_topics
  factor_names <- aligned$factor_names
  comparison_counts <- X
  comparison_terms <- terms
  word_plot_floor <- 1e-6
  term_document_frequency <- as.numeric(
    Matrix::colSums(comparison_counts > 0)
  ) / nrow(comparison_counts)

  keyword_score <- data.frame(
    factor = factor_names,
    EBPMF = colSums(F_ebpmf_topics * term_document_frequency),
    fastTopics = colSums(F_fast_topics * term_document_frequency)
  )
  keyword_score$mean_score <- rowMeans(
    keyword_score[c("EBPMF", "fastTopics")]
  )
  keyword_factor <- which.max(keyword_score$mean_score)

  make_keyword_data <- function(F, method, k, label_count = 12L) {
    topic_frequency <- pmax(F[, k], word_plot_floor)
    other_frequency <- pmax(
      apply(F[, -k, drop = FALSE], 1, max),
      word_plot_floor
    )
    label_score <- topic_frequency * pmax(
      log2(topic_frequency / other_frequency),
      0
    )
    label_index <- head(order(label_score, decreasing = TRUE), label_count)

    data.frame(
      word = comparison_terms,
      other_frequency = other_frequency,
      topic_frequency = topic_frequency,
      label = ifelse(seq_along(comparison_terms) %in% label_index,
                     comparison_terms, ""),
      method = method
    )
  }

  keyword_plot_data <- rbind(
    make_keyword_data(F_ebpmf_topics, "EBPMF", keyword_factor),
    make_keyword_data(F_fast_topics, "fastTopics", keyword_factor)
  )
  keyword_frequency_limits <- range(c(
    keyword_plot_data$other_frequency,
    keyword_plot_data$topic_frequency
  ))

  print(ggplot(
    keyword_plot_data,
    aes(other_frequency, topic_frequency, label = label)
  ) +
    geom_point(color = "steelblue", alpha = 0.35, size = 1) +
    geom_abline(linetype = 2, color = "firebrick") +
    ggrepel::geom_text_repel(
      color = "grey35",
      size = 2.5,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 1
    ) +
    scale_x_log10(limits = keyword_frequency_limits) +
    scale_y_log10(limits = keyword_frequency_limits) +
    coord_equal() +
    facet_wrap(~ method, nrow = 1) +
    labs(
      x = "Highest frequency in another topic",
      y = paste("Frequency in", factor_names[keyword_factor]),
      title = paste("Keyword-like aligned topic:",
                    factor_names[keyword_factor])
    ) +
    theme_cowplot(font_size = 11))

  transform(
    keyword_score,
    EBPMF = round(EBPMF, 4),
    fastTopics = round(fastTopics, 4),
    mean_score = round(mean_score, 4)
  )
}

# Newsgroups diagnostic: keywords.
# Rank EBPMF posterior mean F within each topic.
# Include shared background f0 and effective f0 * F without changing the ranking.
news_top_words <- function(fits, terms, n = 5L) {
  F_ebpmf <- fits$ebpmf$state$qF$mean
  f0 <- fits$ebpmf$state$f0
  stopifnot(nrow(F_ebpmf) == length(terms),
            length(f0) == length(terms),
            length(n) == 1L, is.finite(n), n >= 1L, n == floor(n))
  result <- do.call(rbind, lapply(seq_len(ncol(F_ebpmf)), function(k) {
    ebpmf_rows <- head(order(-F_ebpmf[, k], seq_along(terms)), n)
    data.frame(
      topic = paste0("k", k), rank = seq_along(ebpmf_rows),
      EBPMF_word = terms[ebpmf_rows], EBPMF_F = F_ebpmf[ebpmf_rows, k],
      f0 = f0[ebpmf_rows],
      effective_F = f0[ebpmf_rows] * F_ebpmf[ebpmf_rows, k]
    )
  }))
  # Report words selected in more than one distinct topic, without markup.
  pairs <- unique(result[c("topic", "EBPMF_word")])
  topics_by_word <- split(pairs$topic, pairs$EBPMF_word)
  recurring <- topics_by_word[lengths(topics_by_word) > 1L]
  duplicate_keywords <- data.frame(
    word = as.character(names(recurring)),
    topics = unname(vapply(recurring, paste, collapse = ", ", character(1)))
  )
  list(top_words = result, duplicate_keywords = duplicate_keywords)
}

# List topic-neutral words: posterior mean F closest to its prior mean of one.
# Include f0 and f0 * F so neutral multipliers can be distinguished from common
# background words, and report words selected for more than one topic.
news_background_words <- function(fits, terms, n = 5L) {
  F_ebpmf <- fits$ebpmf$state$qF$mean
  f0 <- fits$ebpmf$state$f0
  stopifnot(nrow(F_ebpmf) == length(terms),
            length(f0) == length(terms),
            length(n) == 1L, is.finite(n), n >= 1L, n == floor(n))

  result <- do.call(rbind, lapply(seq_len(ncol(F_ebpmf)), function(k) {
    distance_from_one <- abs(F_ebpmf[, k] - 1)
    rows <- head(order(distance_from_one, seq_along(terms)), n)
    data.frame(
      topic = paste0("k", k), rank = seq_along(rows),
      EBPMF_word = terms[rows], EBPMF_F = F_ebpmf[rows, k],
      distance_from_1 = distance_from_one[rows], f0 = f0[rows],
      effective_F = f0[rows] * F_ebpmf[rows, k]
    )
  }))

  pairs <- unique(result[c("topic", "EBPMF_word")])
  topics_by_word <- split(pairs$topic, pairs$EBPMF_word)
  recurring <- topics_by_word[lengths(topics_by_word) > 1L]
  duplicate_background_words <- data.frame(
    word = as.character(names(recurring)),
    topics = unname(vapply(recurring, paste, collapse = ", ", character(1)))
  )

  list(
    background_words = result,
    duplicate_background_words = duplicate_background_words
  )
}

# Inspect documents containing a word, ordered by EBPMF topic membership.
# text, if supplied, must be a character vector aligned exactly with rows of X.
news_word_documents <- function(news, theta, word, k, n = 2L,
                                lowest = NULL, text = NULL, preview_words = 30L) {
  j <- match(word, news$terms)
  if (is.na(j)) stop("Word is not in the vocabulary: ", word)
  stopifnot(nrow(theta) == nrow(news$X), length(k) == 1L,
            k >= 1L, k <= ncol(theta), k == floor(k), n >= 1L,
            is.null(text) || length(text) == nrow(news$X))
  counts <- as.numeric(news$X[, j])
  containing <- which(counts > 0)
  # By default select n documents at each extreme; groups may overlap.
  extremes <- if (is.null(lowest)) c(highest = FALSE, lowest = TRUE) else
    setNames(lowest, if (lowest) "lowest" else "highest")
  selected <- lapply(extremes, function(low) {
    head(containing[order(if (low) theta[containing, k] else
      -theta[containing, k], containing)], n)
  })
  rows <- as.integer(unlist(selected, use.names = FALSE))
  documents <- data.frame(
    membership = rep(names(selected), lengths(selected)),
    document_row = rows, newsgroup = as.character(news$labels[rows]),
    topic_proportion = unname(theta[rows, k]), word_count = counts[rows],
    total_tokens = as.numeric(Matrix::rowSums(news$X[rows, , drop = FALSE]))
  )
  documents$word_counts <- vapply(rows, function(i) {
    x <- as.numeric(news$X[i, ])
    index <- head(order(-x, seq_along(x)), preview_words)
    index <- index[x[index] > 0]
    paste(paste0(news$terms[index], " (", x[index], ")"), collapse = ", ")
  }, character(1))
  if (!is.null(text)) documents$text <- text[rows]
  list(
    summary = data.frame(word = word, topic = paste0("k", k),
      documents_containing = length(containing),
      fraction_of_documents = length(containing) / nrow(news$X),
      total_occurrences = sum(counts)),
    documents = documents
  )
}

# Document frequency and raw F for each topic's selected top words.
news_top_word_prevalence <- function(news, fits, n = 5L) {
  words <- news_top_words(fits, news$terms, n)$top_words
  j <- match(words$EBPMF_word, news$terms)
  frequency <- Matrix::colSums(news$X > 0)
  totals <- Matrix::colSums(news$X)
  words$documents_containing <- as.numeric(frequency[j])
  words$fraction_of_documents <- as.numeric(frequency[j]) / nrow(news$X)
  words$total_occurrences <- as.numeric(totals[j])
  words
}

# Soft size sums membership; dominant size assigns each document to one topic.
news_topic_sizes <- function(news, theta) {
  stopifnot(nrow(theta) == nrow(news$X))
  mass <- colSums(theta)
  data.frame(
    effective_documents = mass,
    dominant_documents = tabulate(max.col(theta, ties.method = "first"),
                                  nbins = ncol(theta)),
    row.names = paste0("k", seq_len(ncol(theta)))
  )
}

# Compare global raw-F keywords with empirical counts in the strongest documents.
# This is descriptive: selecting documents using the fitted model is not a
# statistical test, and local word counts are not re-estimated F parameters.
news_strong_document_words <- function(news, fits, theta, k,
                                       n_documents = 10L, n_words = 5L) {
  F <- fits$ebpmf$state$qF$mean
  stopifnot(nrow(theta) == nrow(news$X), nrow(F) == ncol(news$X),
            k >= 1L, k <= ncol(theta), k == floor(k),
            n_documents >= 1L, n_words >= 1L)
  rows <- head(order(-theta[, k], seq_len(nrow(theta))), n_documents)
  counts <- as.numeric(Matrix::colSums(news$X[rows, , drop = FALSE]))
  global <- head(order(-F[, k], seq_len(nrow(F))), n_words)
  local <- which(counts > 0)
  local <- head(local[order(-counts[local], local)], n_words)
  intersection <- intersect(global, local)
  list(
    documents = data.frame(document_row = rows,
      newsgroup = as.character(news$labels[rows]), topic_proportion = theta[rows, k]),
    global_top_words = data.frame(word = news$terms[global], F = F[global, k],
      count_in_selected_documents = counts[global]),
    selected_document_top_words = data.frame(word = news$terms[local],
      count = counts[local], frequency = counts[local] / sum(counts),
      global_F = F[local, k]),
    comparison = data.frame(topic = paste0("k", k),
      documents_selected = length(rows), words_shared = length(intersection),
      jaccard = length(intersection) / length(union(global, local))),
    shared_words = news$terms[intersection]
  )
}

# Compare EBPMF states fit to the SAME document and vocabulary rows.
# Match each fit to the first using document-topic and topic-word profiles.
news_word_stability <- function(states, terms, n = 5L) {
  if (length(states) < 2L) stop("Provide at least two EBPMF states fitted to the same data.")
  stopifnot(n >= 1L, n == floor(n))
  if (is.null(names(states))) names(states) <- paste0("fit", seq_along(states))
  reference <- states[[1L]]
  profiles <- function(state) {
    B <- sweep(state$qF$mean, 1L, state$f0, "*")
    list(L = normalize_rows(sweep(state$L, 2L, colSums(B), "*")),
         F = normalize_columns(B))
  }
  ref <- profiles(reference)
  K <- ncol(reference$L)
  words <- comparisons <- vector("list", length(states))
  for (i in seq_along(states)) {
    state <- states[[i]]
    stopifnot(identical(dim(state$L), dim(reference$L)),
              identical(dim(state$qF$mean), dim(reference$qF$mean)),
              nrow(state$qF$mean) == length(terms))
    current <- profiles(state)
    matching <- match_factor_order(ref$L, ref$F, current$L, current$F,
                                   criterion = "profiles")
    permutation <- if (i == 1L) seq_len(K) else matching$permutation
    F <- state$qF$mean[, permutation, drop = FALSE]
    words[[i]] <- do.call(rbind, lapply(seq_len(K), function(k) {
      j <- head(order(-F[, k], seq_along(terms)), n)
      data.frame(fit = names(states)[i], topic = paste0("k", k),
                 rank = seq_along(j), word = terms[j], F = F[j, k])
    }))
    comparisons[[i]] <- do.call(rbind, lapply(seq_len(K), function(k) {
      original <- head(order(-reference$qF$mean[, k], seq_along(terms)), n)
      fitted <- head(order(-F[, k], seq_along(terms)), n)
      shared <- length(intersect(original, fitted))
      data.frame(fit = names(states)[i], topic = paste0("k", k),
        fitted_column = permutation[k], words_shared = shared,
        fraction_retained = shared / length(original),
        jaccard = shared / length(union(original, fitted)),
        document_cosine = matching$document_similarity[k, permutation[k]],
        word_cosine = matching$word_similarity[k, permutation[k]])
    }))
  }
  words <- do.call(rbind, words)
  presence <- unique(words[c("fit", "topic", "word")])
  recurrence <- aggregate(fit ~ topic + word, presence, length)
  names(recurrence)[3L] <- "fits_containing"
  recurrence$fraction_of_fits <- recurrence$fits_containing / length(states)
  list(comparison_to_reference = do.call(rbind, comparisons),
       top_words = words, recurrence = recurrence)
}

news_distinctive_words <- function(aligned, terms) {
  F_ebpmf_topics <- aligned$F_ebpmf_topics
  F_fast_topics <- aligned$F_fast_topics
  factor_names <- aligned$factor_names
  comparison_terms <- terms
  top_distinctive_words <- function(F, n = 8L) {
    vapply(seq_len(ncol(F)), function(k) {
      topic_frequency <- F[, k]
      other_frequency <- apply(F[, -k, drop = FALSE], 1, max)
      score <- topic_frequency * pmax(
        log2((topic_frequency + 1e-12) / (other_frequency + 1e-12)),
        0
      )
      paste(comparison_terms[head(order(score, decreasing = TRUE), n)],
            collapse = ", ")
    }, character(1))
  }

  keyword_comparison <- data.frame(
    factor = factor_names,
    EBPMF = top_distinctive_words(F_ebpmf_topics),
    fastTopics = top_distinctive_words(F_fast_topics),
    check.names = FALSE
  )

  keyword_comparison
}
