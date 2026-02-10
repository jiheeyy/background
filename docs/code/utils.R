# Follows mostly from github.com/stephenslab/ebpm/gammamix.R with little changes
# Takes 3 parameters, then creates data frame with 3 corresponding columns
# Assigns custom class "gammamix"
gammamix <- function(pi=NULL) {
  grid = c(1e-16, 1e-8, seq(0.01, 0.10, 0.01), seq(0.2, 0.9, 0.1), seq(1,15,2), 20, 50, 75, 100, 200, 1000)
  len_grid = length(grid)

  # Default: equal weight
  if (is.null(pi)) {pi = rep(1/len_grid, len_grid)}

  # Safety check
  if (length(pi) != len_grid) stop("Length of 'pi' must match the grid length.")

  return(structure(data.frame(pi = pi, shape = grid, rate = grid),
                       class = "gammamix"))
  }

dnbinom_cts_log <- function(x, size, prob) {
  # x:    Data vector (length N)
  # size: Shape parameter vector (length M) (Your 'rates' or 'shapes')
  # prob: Probability vector (length M) (Your 'nbprob')

  # 1. Expand inputs to matrices (N x M)
  # This handles the 'recycling' correctly so every data point
  # is evaluated against every mixture component.

  # Term 1: x * log(1-p)
  # Use outer product to match every x with every p
  log_1_p_term <- outer(x, log(1 - prob), "*")

  # Term 2: size * log(p)
  # Create a matrix where every row is size * log(p)
  size_log_p   <- size * log(prob)
  # Expand to N rows
  size_log_p_term <- matrix(size_log_p, nrow = length(x),
                            ncol = length(size), byrow = TRUE)

  # Term 3: lgamma(size)
  lgamma_size_term <- matrix(lgamma(size), nrow = length(x),
                             ncol = length(size), byrow = TRUE)

  # Term 4: lgamma(x + size)
  lgamma_x_size_term <- lgamma(outer(x, size, "+"))

  # Term 5: lgamma(x + 1)
  # This only varies by x (rows), repeated across columns
  lgamma_x_1_term <- matrix(lgamma(x + 1), nrow = length(x),
                            ncol = length(size), byrow = FALSE)

  # Combine
  out <- lgamma_x_size_term - lgamma_size_term - lgamma_x_1_term +
    size_log_p_term + log_1_p_term

  # Handle the 0 * -Inf case for x=0
  # If x is 0, the x*log(1-p) term should be 0, even if log(1-p) is -Inf
  if (any(x == 0)) {
    # Find indices where x is 0
    zero_indices <- which(x == 0)
    # Set the relevant rows of the log_1_p_term to 0 (since 0 * log(1-p) = 0)
    # We must correct 'out' by removing the NaN/Inf from that term and adding 0
    # Easier strategy: Re-calculate 'out' for these rows carefully or
    # rely on R's finite check.
    # Better: fix the input term directly.
    log_1_p_term[zero_indices, ] <- 0

    # Re-sum with the corrected term
    out <- lgamma_x_size_term - lgamma_size_term - lgamma_x_1_term +
      size_log_p_term + log_1_p_term
  }

  return(out)
}

