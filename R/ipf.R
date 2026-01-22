#' Estimate Leunbach Score Parameters using Iterative Proportional Fitting
#'
#' @description
#' Estimates score parameters for Leunbach test equating based on the power series
#' distribution framework.  The total score distribution in Rasch models follows a
#' power series distribution where score parameters are defined by generalized
#' symmetric functions of the item parameters.
#'
#' @param data A data. frame or matrix with two columns:
#'             Column 1: Sum scores from Test 1
#'             Column 2: Sum scores from Test 2
#' @param max_score1 Maximum possible score for Test 1 (default: max observed)
#' @param max_score2 Maximum possible score for Test 2 (default: max observed)
#' @param max_iter Maximum number of iterations for IPF
#' @param tol Convergence tolerance
#' @param verbose Print iteration progress
#'
#' @return A list of class "leunbach_ipf" containing:
#'   - gamma: Score parameters for Test1 (generalized symmetric functions)
#'   - delta: Score parameters for Test2 (generalized symmetric functions)
#'   - sigma: Score parameters for total score (Test1 + Test2)
#'   - log_gamma: Log of gamma parameters
#'   - log_delta: Log of delta parameters
#'   - log_sigma: Log of sigma parameters
#'   - fitted: Fitted frequency table
#'   - observed:  Observed frequency table
#'   - iterations: Number of iterations to convergence
#'   - converged:  Logical indicating convergence
#'   - test1_scores: Score values for Test1
#'   - test2_scores: Score values for Test2
#'   - total_scores: Score values for total score
#'   - g_sq: Likelihood ratio test statistic
#'   - chi_sq: Pearson chi-square statistic
#'   - df:  Degrees of freedom
#'   - p_value: P-value for likelihood ratio test
#'
#' @details
#' The model is based on the power series distribution (PSD) for sum scores in
#' Rasch models.  For a test with k polytomous items, the probability of obtaining
#' sum score r is: 
#'
#' \deqn{P(R = r | \theta) = \frac{\gamma_r \theta^r}{\sum_s \gamma_s \theta^s}}
#'
#' where \eqn{\gamma_r} is the generalized symmetric function of order r of the
#' item parameters, and \eqn{\theta = exp(\beta)} is the person parameter.
#'
#' The model is over-parameterized, so we apply identification constraints: 
#' - \eqn{\gamma_{min}} = 1 (lowest observed score parameter fixed at 1)
#' - \eqn{\delta_{min}} = 1 (lowest observed score parameter fixed at 1)
#' - The highest score parameters are constrained so that their product
#'   equals the product of all score parameters (Rasch-style constraint)
#'
#' The sigma parameters for the total score (S = X + Y) are computed as:
#' \deqn{\sigma_s = \sum_{x+y=s} \gamma_x \delta_y}
#'
#' @references
#' Leunbach, G. (1976). A probabilistic measurement model for assessing whether
#' two tests measure the same personal factor. Copenhagen: Danish Institute for
#' Educational Research.
#' 
#' Adroher, N. D., Kreiner, S., Young, C., Mills, R., & Tennant, A. (2019). 
#' Test equating sleep scales: Applying the Leunbach’s model. 
#' BMC Medical Research Methodology, 19(1), 141. 
#' <https://doi.org/10.1186/s12874-019-0768-y>
#'
#' Kreiner, S. (2007). Validity and objectivity:  Reflections on the role and
#' nature of Rasch models. Nordic Psychology, 59(3), 268-298. 
#' <https://doi.org/10.1027/1901-2276.59.3.268>
#'
#' @examples
#' # Simulate test score data
#' set.seed(123)
#' n <- 500
#' theta <- rnorm(n)
#' test1 <- pmin(pmax(round(3 + 1.5 * theta + rnorm(n, sd = 0.8)), 0), 6)
#' test2 <- pmin(pmax(round(2.5 + 1.3 * theta + rnorm(n, sd = 0.7)), 0), 5)
#' score_data <- data.frame(test1 = test1, test2 = test2)
#'
#' # Estimate parameters (specifying max possible scores)
#' fit <- leunbach_ipf(score_data, max_score1 = 6, max_score2 = 5, verbose = TRUE)
#' print(fit)
#' summary(fit)
#'
#' @export
leunbach_ipf <- function(data, max_score1 = NULL, max_score2 = NULL,
                         max_iter = 1000, tol = 1e-10, verbose = FALSE) {
  
  # Validate input
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data.frame or matrix")
  }
  
  data <- as.matrix(data)
  
  if (ncol(data) != 2) {
    stop("Input must have exactly two columns (Test1 scores, Test2 scores)")
  }
  
  # Extract scores
  test1 <- data[, 1]
  test2 <- data[, 2]
  
  # Remove any rows with NA
  complete_cases <- complete.cases(test1, test2)
  if (sum(! complete_cases) > 0) {
    warning(sprintf("Removed %d rows with missing values", sum(!complete_cases)))
    test1 <- test1[complete_cases]
    test2 <- test2[complete_cases]
  }
  
  # Determine score ranges
  # Use provided max scores or default to max observed
  if (is.null(max_score1)) {
    max_score1 <- max(test1)
  }
  
  if (is.null(max_score2)) {
    max_score2 <- max(test2)
  }
  
  # Full possible score ranges (for parameter vectors)
  test1_scores <- 0:max_score1
  test2_scores <- 0:max_score2
  
  n_test1 <- length(test1_scores)
  n_test2 <- length(test2_scores)
  
  # Total score range
  smin <- 0
  smax <- max_score1 + max_score2
  total_scores <- smin: smax
  
  # Create observed frequency table
  obs <- table(factor(test1, levels = test1_scores),
               factor(test2, levels = test2_scores))
  obs <- as.matrix(obs)
  rownames(obs) <- test1_scores
  colnames(obs) <- test2_scores
  
  # Get observed marginals (sufficient statistics)
  row_margins <- rowSums(obs)
  col_margins <- colSums(obs)
  total_n <- sum(obs)
  
  # Find observed score ranges (XMIN, XMAX, YMIN, YMAX in DIGRAM)
  # These are the first/last scores with observations > 0
  xmin <- test1_scores[min(which(row_margins > 0))]
  xmax <- test1_scores[max(which(row_margins > 0))]
  ymin <- test2_scores[min(which(col_margins > 0))]
  ymax <- test2_scores[max(which(col_margins > 0))]
  
  # Observed total score range
  sfra <- xmin + ymin
  stil <- xmax + ymax
  
  if (verbose) {
    cat("Leunbach Test Equating - Score Parameter Estimation\n")
    cat("====================================================\n\n")
    cat("Power Series Distribution Model with Generalized Symmetric Functions\n\n")
    cat(sprintf("Data summary:\n"))
    cat(sprintf("  N = %d observations\n", total_n))
    cat(sprintf("  Test1: scores %d to %d (observed:  %d to %d)\n",
                min(test1_scores), max(test1_scores), xmin, xmax))
    cat(sprintf("  Test2: scores %d to %d (observed:  %d to %d)\n",
                min(test2_scores), max(test2_scores), ymin, ymax))
    cat(sprintf("  Total:  scores %d to %d (observed: %d to %d)\n\n",
                smin, smax, sfra, stil))
    cat("Estimating score parameters via IPF...\n\n")
  }
  
  # Initialize score parameters
  # gamma[x] = 1 if x has observations, 0 otherwise
  # delta[y] = 1 if y has observations, 0 otherwise
  gamma <- rep(0, n_test1)
  names(gamma) <- test1_scores
  gamma[row_margins > 0] <- 1
  
  delta <- rep(0, n_test2)
  names(delta) <- test2_scores
  delta[col_margins > 0] <- 1
  
  # Create matrix of total scores for each cell
  score_matrix <- outer(test1_scores, test2_scores, "+")
  rownames(score_matrix) <- test1_scores
  colnames(score_matrix) <- test2_scores
  
  # Compute observed frequency for each total score (SMARGIN)
  total_score_freq <- numeric(length(total_scores))
  names(total_score_freq) <- total_scores
  for (s in total_scores) {
    cells <- which(score_matrix == s, arr.ind = TRUE)
    if (nrow(cells) > 0) {
      total_score_freq[as.character(s)] <- sum(obs[cells])
    }
  }
  
  # Initial adjustment of gamma parameters
  result <- adjust_gamma(gamma, delta, test1_scores, test2_scores,
                         xmin, xmax, ymin, ymax, sfra, stil,
                         smin, smax, total_scores)
  gamma <- result$gamma
  delta <- result$delta
  sigma <- result$sigma
  
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    old_gamma <- gamma
    old_delta <- delta
    
    # Compute expected frequencies given current parameters
    expected <- matrix(0, nrow = n_test1, ncol = n_test2)
    rownames(expected) <- test1_scores
    colnames(expected) <- test2_scores
    
    for (s in total_scores) {
      s_char <- as.character(s)
      if (total_score_freq[s_char] > 0) {
        cells <- which(score_matrix == s, arr.ind = TRUE)
        if (nrow(cells) > 0) {
          cell_weights <- sapply(1:nrow(cells), function(i) {
            x <- test1_scores[cells[i, 1]]
            y <- test2_scores[cells[i, 2]]
            gamma[as.character(x)] * delta[as.character(y)]
          })
          if (sum(cell_weights) > 0) {
            cell_probs <- cell_weights / sum(cell_weights)
            for (i in 1:nrow(cells)) {
              expected[cells[i, 1], cells[i, 2]] <-
                total_score_freq[s_char] * cell_probs[i]
            }
          }
        }
      }
    }
    
    # Update gamma to match row marginals
    expected_row_margins <- rowSums(expected)
    for (x in test1_scores) {
      x_char <- as.character(x)
      x_idx <- which(test1_scores == x)
      if (expected_row_margins[x_idx] > 0 && row_margins[x_idx] > 0) {
        gamma[x_char] <- gamma[x_char] * (row_margins[x_idx] / expected_row_margins[x_idx])
      }
    }
    
    # Recompute expected with updated gamma
    expected <- matrix(0, nrow = n_test1, ncol = n_test2)
    for (s in total_scores) {
      s_char <- as.character(s)
      if (total_score_freq[s_char] > 0) {
        cells <- which(score_matrix == s, arr.ind = TRUE)
        if (nrow(cells) > 0) {
          cell_weights <- sapply(1:nrow(cells), function(i) {
            x <- test1_scores[cells[i, 1]]
            y <- test2_scores[cells[i, 2]]
            gamma[as.character(x)] * delta[as.character(y)]
          })
          if (sum(cell_weights) > 0) {
            cell_probs <- cell_weights / sum(cell_weights)
            for (i in 1:nrow(cells)) {
              expected[cells[i, 1], cells[i, 2]] <-
                total_score_freq[s_char] * cell_probs[i]
            }
          }
        }
      }
    }
    
    # Update delta to match column marginals
    expected_col_margins <- colSums(expected)
    for (y in test2_scores) {
      y_char <- as.character(y)
      y_idx <- which(test2_scores == y)
      if (expected_col_margins[y_idx] > 0 && col_margins[y_idx] > 0) {
        delta[y_char] <- delta[y_char] * (col_margins[y_idx] / expected_col_margins[y_idx])
      }
    }
    
    # Apply ADJUST_GAMMA normalization
    result <- adjust_gamma(gamma, delta, test1_scores, test2_scores,
                           xmin, xmax, ymin, ymax, sfra, stil,
                           smin, smax, total_scores)
    gamma <- result$gamma
    delta <- result$delta
    sigma <- result$sigma
    
    # Check convergence
    max_diff_gamma <- max(abs(gamma - old_gamma))
    max_diff_delta <- max(abs(delta - old_delta))
    max_diff <- max(max_diff_gamma, max_diff_delta)
    
    if (verbose) {
      cat(sprintf("  Iteration %d: max param change = %.2e (gamma:  %.2e, delta: %.2e)\n",
                  iter, max_diff, max_diff_gamma, max_diff_delta))
    }
    
    if (max_diff < tol) {
      converged <- TRUE
      break
    }
  }
  
  # Final normalization and compute sigma
  result <- adjust_gamma(gamma, delta, test1_scores, test2_scores,
                         xmin, xmax, ymin, ymax, sfra, stil,
                         smin, smax, total_scores)
  gamma <- result$gamma
  delta <- result$delta
  sigma <- result$sigma
  
  # Final expected frequencies
  expected <- matrix(0, nrow = n_test1, ncol = n_test2)
  rownames(expected) <- test1_scores
  colnames(expected) <- test2_scores
  
  for (s in total_scores) {
    s_char <- as.character(s)
    if (total_score_freq[s_char] > 0) {
      cells <- which(score_matrix == s, arr.ind = TRUE)
      if (nrow(cells) > 0) {
        cell_weights <- sapply(1:nrow(cells), function(i) {
          x <- test1_scores[cells[i, 1]]
          y <- test2_scores[cells[i, 2]]
          gamma[as.character(x)] * delta[as.character(y)]
        })
        if (sum(cell_weights) > 0) {
          cell_probs <- cell_weights / sum(cell_weights)
          for (i in 1:nrow(cells)) {
            expected[cells[i, 1], cells[i, 2]] <-
              total_score_freq[s_char] * cell_probs[i]
          }
        }
      }
    }
  }
  
  # Compute log parameters
  log_gamma <- log(gamma)
  log_gamma[! is.finite(log_gamma)] <- NA
  log_delta <- log(delta)
  log_delta[!is.finite(log_delta)] <- NA
  log_sigma <- log(sigma)
  log_sigma[!is.finite(log_sigma)] <- NA
  
  if (verbose) {
    cat("\n")
    if (converged) {
      cat(sprintf("Converged after %d iterations\n", iter))
    } else {
      cat(sprintf("Did not converge after %d iterations\n", max_iter))
    }
  }
  
  # Build result object
  result <- list(
    gamma = gamma,
    delta = delta,
    sigma = sigma,
    log_gamma = log_gamma,
    log_delta = log_delta,
    log_sigma = log_sigma,
    fitted = expected,
    observed = obs,
    iterations = iter,
    converged = converged,
    test1_scores = test1_scores,
    test2_scores = test2_scores,
    total_scores = total_scores,
    n = total_n,
    row_margins = row_margins,
    col_margins = col_margins,
    total_score_freq = total_score_freq,
    # Observed score ranges
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    sfra = sfra, stil = stil
  )
  
  class(result) <- "leunbach_ipf"
  
  # Calculate fit statistics
  fit_stats <- calculate_statistics(result)
  result$g_sq <- fit_stats$g_sq
  result$chi_sq <- fit_stats$chi_sq
  result$df <- fit_stats$df
  result$p_value <- fit_stats$p_g_sq
  result$std_residuals <- fit_stats$std_residuals
  result$loglike <- fit_stats$loglike
  result$loglike1 <- fit_stats$loglike1
  result$orbit_df <- fit_stats$orbit_df
  result$gk_gamma_observed <- fit_stats$gk_gamma_observed
  result$gk_gamma_expected <- fit_stats$gk_gamma_expected
  result$gk_gamma_se <- fit_stats$gk_gamma_se
  result$gk_gamma_z <- fit_stats$gk_gamma_z
  result$gk_gamma_p <- fit_stats$gk_gamma_p
  
  return(result)
}


#' Adjust gamma and delta parameters (ADJUST_GAMMA from Pascal)
#'
#' @description
#' Applies the identification constraints to the score parameters:
#' 1. Normalize so lowest observed score parameter = 1
#' 2. Apply log-linear adjustment so product of highest parameters is constrained
#' 3. Compute sigma (total score parameters)
#'
#' @keywords internal
adjust_gamma <- function(gamma, delta, test1_scores, test2_scores,
                         xmin, xmax, ymin, ymax, sfra, stil,
                         smin, smax, total_scores) {
  
  # Step 1: Normalize gamma so gamma[xmin] = 1
  alpha <- gamma[as.character(xmin)]
  if (! is.na(alpha) && alpha != 1.0 && alpha > 0) {
    gamma <- gamma / alpha
  }
  
  # Step 2: Normalize delta so delta[ymin] = 1
  alpha <- delta[as.character(ymin)]
  if (!is.na(alpha) && alpha != 1.0 && alpha > 0) {
    delta <- delta / alpha
  }
  
  # Step 3: Log-linear adjustment
  # Compute alpha = -ln(gamma[xmax] * delta[ymax]) / (stil - sfra)
  gamma_xmax <- gamma[as.character(xmax)]
  delta_ymax <- delta[as.character(ymax)]
  
  if (! is.na(gamma_xmax) && !is.na(delta_ymax) &&
      gamma_xmax > 0 && delta_ymax > 0 && stil != sfra) {
    
    alpha <- -log(gamma_xmax * delta_ymax) / (stil - sfra)
    
    if (alpha != 0.0) {
      # Adjust gamma:  gamma[x] = exp(alpha * (x - xmin) + ln(gamma[x]))
      for (x in test1_scores) {
        x_char <- as.character(x)
        if (! is.na(gamma[x_char]) && gamma[x_char] > 0) {
          gamma[x_char] <- exp(alpha * (x - xmin) + log(gamma[x_char]))
        }
      }
      
      # Adjust delta: delta[y] = exp(alpha * (y - ymin) + ln(delta[y]))
      for (y in test2_scores) {
        y_char <- as.character(y)
        if (!is.na(delta[y_char]) && delta[y_char] > 0) {
          delta[y_char] <- exp(alpha * (y - ymin) + log(delta[y_char]))
        }
      }
    }
  }
  
  # Step 4: Compute sigma (total score parameters)
  # sigma[s] = sum over all x,y where x+y=s of gamma[x] * delta[y]
  sigma <- numeric(length(total_scores))
  names(sigma) <- total_scores
  
  for (s in total_scores) {
    s_char <- as.character(s)
    sigma[s_char] <- 0
    for (x in test1_scores) {
      y <- s - x
      if (y >= min(test2_scores) && y <= max(test2_scores)) {
        x_char <- as.character(x)
        y_char <- as.character(y)
        if (!is.na(gamma[x_char]) && !is.na(delta[y_char])) {
          sigma[s_char] <- sigma[s_char] + gamma[x_char] * delta[y_char]
        }
      }
    }
  }
  
  list(gamma = gamma, delta = delta, sigma = sigma)
}


#' Calculate fit statistics (CALCULATE_STATISTICS from Pascal)
#'
#' @description
#' Computes the likelihood ratio test statistic, Goodman-Kruskal Gamma test,
#' and degrees of freedom following the DIGRAM implementation.
#'
#' @param object A leunbach_ipf object
#' @return A list with fit statistics
#'
#' @keywords internal
calculate_statistics <- function(object) {
  
  gamma <- object$gamma
  delta <- object$delta
  sigma <- object$sigma
  obs <- object$observed
  expected <- object$fitted
  smargin <- object$total_score_freq
  
  # Get observed score ranges
  xmin <- object$xmin
  xmax <- object$xmax
  ymin <- object$ymin
  ymax <- object$ymax
  
  # Initialize log-likelihoods
  loglike <- 0.0
  loglike1 <- 0.0
  chi_sq <- 0.0
  
  # Matrix for standardized residuals
  std_residuals <- matrix(0, nrow = length(object$test1_scores),
                          ncol = length(object$test2_scores))
  rownames(std_residuals) <- object$test1_scores
  colnames(std_residuals) <- object$test2_scores
  
  for (x in object$test1_scores) {
    x_char <- as.character(x)
    x_idx <- which(object$test1_scores == x)
    
    if (! is.na(gamma[x_char]) && gamma[x_char] > 0) {
      lngamx <- log(gamma[x_char])
    } else {
      lngamx <- 0
    }
    
    for (y in object$test2_scores) {
      y_char <- as.character(y)
      y_idx <- which(object$test2_scores == y)
      s <- x + y
      s_char <- as.character(s)
      
      if (!is.na(delta[y_char]) && delta[y_char] > 0) {
        lngamy <- log(delta[y_char])
      } else {
        lngamy <- 0
      }
      
      if (! is.na(sigma[s_char]) && sigma[s_char] > 0) {
        lngams <- log(sigma[s_char])
      } else {
        lngams <- 0
      }
      
      n_xy <- obs[x_idx, y_idx]
      
      if (n_xy > 0) {
        loglike <- loglike + n_xy * (lngamx + lngamy - lngams)
        loglike1 <- loglike1 + n_xy * (log(n_xy) - log(smargin[s_char]))
      }
      
      exp_xy <- expected[x_idx, y_idx]
      if (! is.na(exp_xy) && exp_xy > 0) {
        d <- n_xy - exp_xy
        chi_sq_cell <- d * d / exp_xy
        chi_sq <- chi_sq + chi_sq_cell
        
        if (d > 0) {
          std_residuals[x_idx, y_idx] <- sqrt(chi_sq_cell)
        } else {
          std_residuals[x_idx, y_idx] <- -sqrt(chi_sq_cell)
        }
      }
    }
  }
  
  # Likelihood ratio statistic
  g_sq <- 2.0 * (loglike1 - loglike)
  
  # Degrees of freedom calculation following Pascal code
  orbit_df <- numeric(length(object$total_scores))
  names(orbit_df) <- object$total_scores
  
  for (s in object$total_scores) {
    s_char <- as.character(s)
    
    xfra_orbit <- max(xmin, s - ymax)
    xtil_orbit <- min(xmax, s - ymin)
    
    if (xtil_orbit >= xfra_orbit) {
      orbit_df[s_char] <- xtil_orbit - xfra_orbit
    } else {
      orbit_df[s_char] <- 0
    }
  }
  
  # Sum ORBIT_DF for each total score with observations
  df <- 0
  for (s in object$total_scores) {
    s_char <- as.character(s)
    if (! is.na(smargin[s_char]) && smargin[s_char] > 0) {
      df <- df + orbit_df[s_char]
    }
  }
  
  # Subtract parameter count using OBSERVED score ranges
  df <- df - (xmax - xmin) - (ymax - ymin) + 1
  
  # Add back 1 for each score where gamma = 0
  for (x in object$test1_scores) {
    x_char <- as.character(x)
    if (is.na(gamma[x_char]) || gamma[x_char] == 0) {
      df <- df + 1
    }
  }
  
  # Add back 1 for each score where delta = 0
  for (y in object$test2_scores) {
    y_char <- as.character(y)
    if (is.na(delta[y_char]) || delta[y_char] == 0) {
      df <- df + 1
    }
  }
  
  df <- max(0, df)
  
  if (df > 0) {
    p_g_sq <- pchisq(g_sq, df, lower.tail = FALSE)
    p_chi_sq <- pchisq(chi_sq, df, lower.tail = FALSE)
  } else {
    p_g_sq <- NA
    p_chi_sq <- NA
  }
  
  # Calculate Goodman-Kruskal Gamma test
  gk_gamma <- calculate_gamma_test(obs, expected)
  
  list(
    g_sq = g_sq,
    chi_sq = chi_sq,
    df = df,
    p_g_sq = p_g_sq,
    p_chi_sq = p_chi_sq,
    loglike = loglike,
    loglike1 = loglike1,
    std_residuals = std_residuals,
    orbit_df = orbit_df,
    # Goodman-Kruskal Gamma results
    gk_gamma_observed = gk_gamma$gamma_observed,
    gk_gamma_expected = gk_gamma$gamma_expected,
    gk_gamma_se = gk_gamma$se_gamma,
    gk_gamma_z = gk_gamma$z_statistic,
    gk_gamma_p = gk_gamma$p_value
  )
}


#' Print method for leunbach_ipf objects
#' @export
print.leunbach_ipf <- function(x, ...) {
  cat("Leunbach Score Parameter Estimation\n")
  cat("====================================\n")
  cat("Power Series Distribution with Generalized Symmetric Functions\n\n")
  
  cat(sprintf("N = %d observations\n\n", x$n))
  
  cat("Test 1 score parameters (gamma):\n\n")
  gamma_df <- data.frame(
    Score = x$test1_scores,
    Frequency = as.numeric(x$row_margins),
    Gamma = round(x$gamma, 6),
    Log_Gamma = round(x$log_gamma, 6)
  )
  print(gamma_df, row.names = FALSE)
  
  cat("\nTest 2 score parameters (delta):\n\n")
  delta_df <- data.frame(
    Score = x$test2_scores,
    Frequency = as.numeric(x$col_margins),
    Delta = round(x$delta, 6),
    Log_Delta = round(x$log_delta, 6)
  )
  print(delta_df, row.names = FALSE)
  
  cat("\nTotal score parameters (sigma = Test1 + Test2):\n\n")
  sigma_df <- data.frame(
    Score = x$total_scores,
    Frequency = as.numeric(x$total_score_freq),
    Sigma = round(x$sigma, 6),
    Log_Sigma = round(x$log_sigma, 6)
  )
  print(sigma_df, row.names = FALSE)
  
  cat(sprintf("\nConverged:  %s (after %d iterations)\n", x$converged, x$iterations))
  
  cat("\n--- Goodness of Fit ---\n\n")
  
  # 1. Likelihood Ratio Test
  cat("1. Likelihood Ratio Test:\n")
  cat(sprintf("   LR = %.2f  DF = %d  p = %.4f\n\n", x$g_sq, x$df, x$p_value))
  
  # 2. Goodman-Kruskal Gamma Test (one-sided)
  cat("2. Goodman-Kruskal Gamma Test (one-sided):\n")
  if (! is.na(x$gk_gamma_observed)) {
    cat(sprintf("   Gamma (observed) = %.4f\n", x$gk_gamma_observed))
    cat(sprintf("   Gamma (expected) = %.4f\n", x$gk_gamma_expected))
    cat(sprintf("   SE = %.4f\n", x$gk_gamma_se))
    cat(sprintf("   Z = %.2f  p = %.4f\n\n", x$gk_gamma_z, x$gk_gamma_p))
  } else {
    cat("   Could not be calculated\n\n")
  }
  
  # 3. Note about orbit analysis (done separately)
  cat("3. Orbit Analysis:\n")
  cat("   Use analyze_orbits() to assess person fit within total score strata\n")
  
  invisible(x)
}


#' Summary method for leunbach_ipf objects
#' @export
summary.leunbach_ipf <- function(object, ...) {
  cat("Leunbach Score Parameter Estimation - Summary\n")
  cat("==============================================\n")
  cat("Power Series Distribution with Generalized Symmetric Functions\n\n")
  
  # Basic info
  cat(sprintf("N = %d observations\n", object$n))
  cat(sprintf("Test 1: scores %d to %d (observed:  %d to %d)\n",
              min(object$test1_scores), max(object$test1_scores),
              object$xmin, object$xmax))
  cat(sprintf("Test 2: scores %d to %d (observed: %d to %d)\n",
              min(object$test2_scores), max(object$test2_scores),
              object$ymin, object$ymax))
  cat(sprintf("Total:  scores %d to %d (observed:  %d to %d)\n\n",
              min(object$total_scores), max(object$total_scores),
              object$sfra, object$stil))
  
  cat("=== Goodness of Fit ===\n\n")
  
  # 1. Likelihood Ratio Test
  cat("1. Likelihood Ratio Test (observed vs expected counts):\n")
  cat(sprintf("   Likelihood ratio G² = %8.2f (df = %d, p = %.4f)\n",
              object$g_sq, object$df, object$p_value))
  cat(sprintf("   Pearson chi-square  = %8.2f\n\n", object$chi_sq))
  
  # 2. Goodman-Kruskal Gamma Test (one-sided)
  cat("2. Goodman-Kruskal Gamma Test (one-sided):\n")
  cat("   Tests if observed correlation exceeds expected under the model.\n")
  if (!is.na(object$gk_gamma_observed)) {
    cat(sprintf("   Gamma (observed)    = %8.4f\n", object$gk_gamma_observed))
    cat(sprintf("   Gamma (expected)    = %8.4f\n", object$gk_gamma_expected))
    cat(sprintf("   Standard error      = %8.4f\n", object$gk_gamma_se))
    cat(sprintf("   Z statistic         = %8.2f (p = %.4f, one-sided)\n\n", 
                object$gk_gamma_z, object$gk_gamma_p))
  } else {
    cat("   Could not be calculated\n\n")
  }
  
  # 3. Orbit analysis note
  cat("3. Orbit Analysis (person fit):\n")
  cat("   Run analyze_orbits() separately to assess the number of cases\n")
  cat("   outside 95% confidence regions of orbit distributions.\n\n")
  
  cat(sprintf("Converged: %s (after %d iterations)\n",
              object$converged, object$iterations))
  
  invisible(list(
    g_sq = object$g_sq,
    chi_sq = object$chi_sq,
    df = object$df,
    p_value = object$p_value,
    gk_gamma_observed = object$gk_gamma_observed,
    gk_gamma_expected = object$gk_gamma_expected,
    gk_gamma_z = object$gk_gamma_z,
    gk_gamma_p = object$gk_gamma_p
  ))
}

#' Plot method for leunbach_ipf objects
#' @export
plot.leunbach_ipf <- function(x, type = c("parameters", "residuals", "observed", "sigma"), ...) {
  type <- match.arg(type)
  
  if (type == "parameters") {
    oldpar <- par(mfrow = c(1, 2))
    on.exit(par(oldpar))
    
    plot(x$test1_scores, x$log_gamma,
         type = "b", pch = 19, col = "blue",
         xlab = "Score", ylab = "Log(gamma)",
         main = "Test 1: Log Score Parameters")
    grid()
    
    plot(x$test2_scores, x$log_delta,
         type = "b", pch = 19, col = "red",
         xlab = "Score", ylab = "Log(delta)",
         main = "Test 2: Log Score Parameters")
    grid()
    
  } else if (type == "sigma") {
    plot(x$total_scores, x$log_sigma,
         type = "b", pch = 19, col = "darkgreen",
         xlab = "Total Score (Test1 + Test2)", ylab = "Log(sigma)",
         main = "Total Score:  Log Score Parameters")
    grid()
    
  } else if (type == "residuals") {
    obs <- x$observed
    exp <- x$fitted
    std_resid <- (obs - exp) / sqrt(exp)
    std_resid[! is.finite(std_resid)] <- 0
    
    image(x$test1_scores, x$test2_scores, std_resid,
          col = colorRampPalette(c("blue", "white", "red"))(21),
          xlab = "Test 1 Score", ylab = "Test 2 Score",
          main = "Standardized Residuals")
    
  } else if (type == "observed") {
    image(x$test1_scores, x$test2_scores, log1p(x$observed),
          col = colorRampPalette(c("white", "steelblue"))(20),
          xlab = "Test 1 Score", ylab = "Test 2 Score",
          main = "Observed Frequencies (log scale)")
  }
  
  invisible(x)
}