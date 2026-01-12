#' @title Helper Functions for Leunbach Model
#' @description Utility functions for bootstrap and equating

#' Generate parametric bootstrap table
#'
#' @description
#' Generates a bootstrap contingency table by sampling from the fitted model.
#' For each total score, samples are drawn from the multinomial distribution
#' defined by the conditional probabilities P(X=x | S=s).
#'
#' @param test1_scores Score values for Test 1
#' @param test2_scores Score values for Test 2
#' @param total_scores Total score values
#' @param total_score_freq Observed frequency for each total score
#' @param gamma Score parameters for Test 1
#' @param delta Score parameters for Test 2
#' @param sigma Total score parameters
#' @param xmin Minimum observed Test 1 score
#' @param xmax Maximum observed Test 1 score
#' @param ymin Minimum observed Test 2 score
#' @param ymax Maximum observed Test 2 score
#'
#' @return A matrix (contingency table) of bootstrap counts
#' @keywords internal
parametric_eq_bootstrap <- function(test1_scores, test2_scores, total_scores,
                                    total_score_freq, gamma, delta, sigma,
                                    xmin, xmax, ymin, ymax) {
  
  n_test1 <- length(test1_scores)
  n_test2 <- length(test2_scores)
  
  boot_table <- matrix(0, nrow = n_test1, ncol = n_test2)
  rownames(boot_table) <- test1_scores
  colnames(boot_table) <- test2_scores
  
  for (s in total_scores) {
    s_char <- as.character(s)
    n_s <- total_score_freq[s_char]
    
    if (! is.na(n_s) && n_s > 0) {
      valid_cells <- list()
      probs <- numeric()
      
      for (x in test1_scores) {
        y <- s - x
        if (y >= min(test2_scores) && y <= max(test2_scores)) {
          x_char <- as.character(x)
          y_char <- as.character(y)
          
          gamma_x <- gamma[x_char]
          delta_y <- delta[y_char]
          
          if (! is.na(gamma_x) && ! is.na(delta_y) && gamma_x > 0 && delta_y > 0) {
            valid_cells[[length(valid_cells) + 1]] <- c(x, y)
            probs <- c(probs, gamma_x * delta_y)
          }
        }
      }
      
      if (length(probs) > 0 && sum(probs) > 0) {
        probs <- probs / sum(probs)
        
        if (length(probs) == 1) {
          counts <- n_s
        } else {
          counts <- as.vector(rmultinom(1, size = n_s, prob = probs))
        }
        
        for (i in seq_along(valid_cells)) {
          x <- valid_cells[[i]][1]
          y <- valid_cells[[i]][2]
          x_idx <- which(test1_scores == x)
          y_idx <- which(test2_scores == y)
          boot_table[x_idx, y_idx] <- counts[i]
        }
      }
    }
  }
  
  return(boot_table)
}


#' Convert contingency table to data frame
#'
#' @description
#' Converts a contingency table (matrix of counts) to a data frame with
#' one row per observation, suitable for input to leunbach_ipf().
#'
#' @param tab Contingency table (matrix)
#' @param test1_scores Score values for Test 1
#' @param test2_scores Score values for Test 2
#'
#' @return A data frame with columns test1 and test2
#' @keywords internal
table_to_data <- function(tab, test1_scores, test2_scores) {
  
  test1 <- numeric()
  test2 <- numeric()
  
  for (i in seq_along(test1_scores)) {
    for (j in seq_along(test2_scores)) {
      count <- tab[i, j]
      if (count > 0) {
        test1 <- c(test1, rep(test1_scores[i], count))
        test2 <- c(test2, rep(test2_scores[j], count))
      }
    }
  }
  
  data.frame(test1 = test1, test2 = test2)
}


#' Calculate Bootstrap Standard Error of Equating
#'
#' @description
#' Computes the standard error of equating from bootstrap samples
#' using the standard deviation of the bootstrap distribution.
#'
#' @param boot_matrix Matrix of bootstrap equated scores (nsim x n_scores)
#'
#' @return Named vector of SEE values for each score
#' @keywords internal
calculate_bootstrap_see <- function(boot_matrix) {
  
  n_scores <- ncol(boot_matrix)
  see <- numeric(n_scores)
  names(see) <- colnames(boot_matrix)
  
  for (i in 1:n_scores) {
    boot_values <- boot_matrix[, i]
    valid <- ! is.na(boot_values)
    S <- sum(valid)
    
    if (S > 1) {
      boot_valid <- boot_values[valid]
      boot_mean <- sum(boot_valid) / S
      see[i] <- sqrt(sum((boot_valid - boot_mean)^2) / (S - 1))
    } else {
      see[i] <- NA
    }
  }
  
  return(see)
}


#' Compute frequency of bootstrap equating errors
#'
#' @description
#' Computes the frequency distribution of equating errors (differences between
#' bootstrap equated scores and observed equated scores) for each source score.
#'
#' @param boot_rounded Matrix of bootstrap rounded scores (nsim x n_scores)
#' @param observed_rounded Vector of observed rounded scores
#' @param scores Vector of source scores
#'
#' @return Matrix of error frequencies (n_scores x 5) for errors -2, -1, 0, +1, +2
#' @keywords internal
compute_error_frequencies <- function(boot_rounded, observed_rounded, scores) {
  
  n_scores <- length(scores)
  error_levels <- c(-2, -1, 0, 1, 2)
  n_levels <- length(error_levels)
  
  freq_matrix <- matrix(0, nrow = n_scores, ncol = n_levels)
  rownames(freq_matrix) <- scores
  colnames(freq_matrix) <- error_levels
  
  for (i in 1:n_scores) {
    obs <- observed_rounded[i]
    if (is.na(obs)) next
    
    boot_vals <- boot_rounded[, i]
    boot_vals <- boot_vals[! is.na(boot_vals)]
    
    if (length(boot_vals) == 0) next
    
    errors <- boot_vals - obs
    
    for (j in 1:n_levels) {
      e <- error_levels[j]
      if (e == -2) {
        freq_matrix[i, j] <- sum(errors <= e)
      } else if (e == 2) {
        freq_matrix[i, j] <- sum(errors >= e)
      } else {
        freq_matrix[i, j] <- sum(errors == e)
      }
    }
    
    total <- length(boot_vals)
    freq_matrix[i, ] <- round(100 * freq_matrix[i, ] / total, 1)
  }
  
  return(freq_matrix)
}


#' Print SEE table for one direction
#'
#' @description
#' Prints a formatted table showing equating results with bootstrap
#' confidence intervals, standard errors, and error frequencies.
#'
#' @param x A leunbach_bootstrap object
#' @param direction "1to2" or "2to1"
#'
#' @keywords internal
print_see_table <- function(x, direction = "1to2") {
  
  conf_pct <- round(x$conf_level * 100)
  
  if (direction == "1to2") {
    cat(sprintf("Equating Test1 to Test2 (with %d%% CI)\n", conf_pct))
    cat("=====================================================\n\n")
    
    eq_table <- x$eq_1to2$equating_table
    see <- x$see_1to2
    ci_lower <- x$ci_lower_1to2
    ci_upper <- x$ci_upper_1to2
    error_freq <- x$error_freq_1to2
    source_min <- x$fit$xmin
    source_max <- x$fit$xmax
    scores <- x$fit$test1_scores
  } else {
    cat(sprintf("Equating Test2 to Test1 (with %d%% CI)\n", conf_pct))
    cat("=====================================================\n\n")
    
    eq_table <- x$eq_2to1$equating_table
    see <- x$see_2to1
    ci_lower <- x$ci_lower_2to1
    ci_upper <- x$ci_upper_2to1
    error_freq <- x$error_freq_2to1
    source_min <- x$fit$ymin
    source_max <- x$fit$ymax
    scores <- x$fit$test2_scores
  }
  
  valid <- scores >= source_min & scores <= source_max
  
  cat("                                                         Frequency of bootstrap errors\n")
  cat(sprintf("Score  Rounded  Expected    %d%% CI          SEE      -2    -1     0    +1    +2\n", conf_pct))
  cat("------------------------------------------------------------------------------------\n")
  
  for (i in which(valid)) {
    ci_str <- sprintf("[%5.2f, %5.2f]", ci_lower[i], ci_upper[i])
    cat(sprintf("%5d  %7d    %5.2f   %15s  %5.2f   %5.1f %5.1f %5.1f %5.1f %5.1f\n",
                eq_table[i, 1], eq_table[i, 3], eq_table[i, 2],
                ci_str, see[i],
                error_freq[i, 1], error_freq[i, 2], error_freq[i, 3],
                error_freq[i, 4], error_freq[i, 5]))
  }
  
  cat("------------------------------------------------------------------------------------\n")
  
  avg_see <- mean(see[valid], na.rm = TRUE)
  cat(sprintf("Average SEE:  %.2f\n", avg_see))
}


#' Diagnostic function for equating
#'
#' @description
#' Diagnoses potential issues with person parameter estimation by displaying
#' theta values and expected scores for each source score.
#'
#' @param fit A leunbach_ipf object
#' @param direction Direction of equating:  "1to2" or "2to1"
#' @param method Optimization method: "optimize" or "newton"
#'
#' @export
diagnose_equating <- function(fit, direction = c("1to2", "2to1"),
                              method = c("optimize", "newton")) {
  
  direction <- match.arg(direction)
  method <- match.arg(method)
  
  if (direction == "1to2") {
    source_scores <- fit$test1_scores
    source_gamma <- fit$gamma
    source_min <- fit$xmin
    source_max <- fit$xmax
    target_gamma <- fit$delta
    target_min <- fit$ymin
    target_max <- fit$ymax
  } else {
    source_scores <- fit$test2_scores
    source_gamma <- fit$delta
    source_min <- fit$ymin
    source_max <- fit$ymax
    target_gamma <- fit$gamma
    target_min <- fit$xmin
    target_max <- fit$xmax
  }
  
  cat("Diagnostic for equating\n")
  cat("=======================\n\n")
  cat(sprintf("Method: %s\n", method))
  cat(sprintf("Source range: %d to %d\n", source_min, source_max))
  cat(sprintf("Target range:  %d to %d\n\n", target_min, target_max))
  
  cat("Score  Theta (raw)    Theta (log)    Expected Target\n")
  cat("-----------------------------------------------------\n")
  
  for (x in source_scores[source_scores >= source_min & source_scores <= source_max]) {
    theta <- estimate_person_parameter(x, source_gamma, source_min, source_max,
                                       method = method)
    
    if (! is.na(theta) && theta > 0) {
      true_score <- calculate_true_score(theta, target_gamma, target_min, target_max)
      cat(sprintf("%5d  %12.6f  %12.6f  %12.2f\n",
                  x, theta, log(theta), true_score))
    } else {
      cat(sprintf("%5d  %12s  %12s  %12s\n", x, "NA", "NA", "NA"))
    }
  }
}


#' Extract equating table with bootstrap confidence intervals
#'
#' @description
#' Extracts equating results from a bootstrap object as a clean data frame
#' with confidence intervals and standard errors.
#'
#' @param boot A leunbach_bootstrap object
#' @param direction "1to2" or "2to1"
#'
#' @return A data frame with equating results and CIs
#' @export
get_equating_table <- function(boot, direction = c("1to2", "2to1")) {
  
  if (! inherits(boot, "leunbach_bootstrap")) {
    stop("Input must be a leunbach_bootstrap object")
  }
  
  direction <- match.arg(direction)
  
  if (direction == "1to2") {
    eq_table <- boot$eq_1to2$equating_table
    see <- boot$see_1to2
    ci_lower <- boot$ci_lower_1to2
    ci_upper <- boot$ci_upper_1to2
    source_min <- boot$fit$xmin
    source_max <- boot$fit$xmax
    scores <- boot$fit$test1_scores
  } else {
    eq_table <- boot$eq_2to1$equating_table
    see <- boot$see_2to1
    ci_lower <- boot$ci_lower_2to1
    ci_upper <- boot$ci_upper_2to1
    source_min <- boot$fit$ymin
    source_max <- boot$fit$ymax
    scores <- boot$fit$test2_scores
  }
  
  valid <- scores >= source_min & scores <= source_max
  
  result <- data.frame(
    source_score = eq_table[valid, 1],
    rounded = eq_table[valid, 3],
    expected = round(eq_table[valid, 2], 2),
    ci_lower = round(ci_lower[valid], 2),
    ci_upper = round(ci_upper[valid], 2),
    see = round(see[valid], 2)
  )
  
  colnames(result)[1] <- colnames(eq_table)[1]
  
  return(result)
}

#' Calculate Goodman-Kruskal Gamma and test statistic
#'
#' @description
#' Computes Goodman & Kruskal's Gamma coefficient for both observed and
#' expected tables, and tests whether they differ significantly.
#' Gamma measures the strength of association between two ordinal variables. 
#' Uses a one-sided test following the DIGRAM implementation.
#'
#' @param observed Observed contingency table
#' @param expected Expected (fitted) contingency table
#'
#' @return A list containing:
#'   - gamma_observed: Gamma coefficient for observed data
#'   - gamma_expected: Gamma coefficient for expected data
#'   - se_gamma:  Standard error of gamma
#'   - z_statistic: Z statistic (expected - observed) / SE
#'   - p_value: One-sided p-value from standard normal
#'
#' @keywords internal
calculate_gamma_test <- function(observed, expected) {
  
  # Calculate Gamma for a contingency table
  # Gamma = (P - Q) / (P + Q)
  # where P = concordant pairs, Q = discordant pairs
  
  calc_gamma <- function(tab) {
    n_row <- nrow(tab)
    n_col <- ncol(tab)
    
    P <- 0  # Concordant pairs
    Q <- 0  # Discordant pairs
    
    for (i in 1:(n_row - 1)) {
      for (j in 1:(n_col - 1)) {
        n_ij <- tab[i, j]
        if (n_ij > 0) {
          # Concordant:  cells below and to the right
          for (k in (i + 1):n_row) {
            for (l in (j + 1):n_col) {
              P <- P + n_ij * tab[k, l]
            }
          }
          # Discordant: cells below and to the left
          for (k in (i + 1):n_row) {
            for (l in 1:(j - 1)) {
              if (l >= 1) {
                Q <- Q + n_ij * tab[k, l]
              }
            }
          }
        }
      }
    }
    
    if ((P + Q) == 0) {
      return(list(gamma = NA, P = P, Q = Q))
    }
    
    gamma_val <- (P - Q) / (P + Q)
    return(list(gamma = gamma_val, P = P, Q = Q))
  }
  
  # Calculate gamma for observed and expected
  obs_result <- calc_gamma(observed)
  exp_result <- calc_gamma(expected)
  
  gamma_observed <- obs_result$gamma
  gamma_expected <- exp_result$gamma
  
  if (is.na(gamma_observed) || is.na(gamma_expected)) {
    return(list(
      gamma_observed = gamma_observed,
      gamma_expected = gamma_expected,
      se_gamma = NA,
      z_statistic = NA,
      p_value = NA
    ))
  }
  
  # Calculate standard error of gamma using ASE1
  # ASE1 = (2 / (P + Q)) * sqrt(sum_ij n_ij * (Q * C_ij - P * D_ij)^2)
  # where C_ij = sum of cells concordant with (i,j)
  #       D_ij = sum of cells discordant with (i,j)
  
  n_row <- nrow(observed)
  n_col <- ncol(observed)
  P <- obs_result$P
  Q <- obs_result$Q
  
  # Pre-compute C_ij and D_ij for each cell
  C_matrix <- matrix(0, nrow = n_row, ncol = n_col)
  D_matrix <- matrix(0, nrow = n_row, ncol = n_col)
  
  for (i in 1:n_row) {
    for (j in 1:n_col) {
      # C_ij:  cells that would be concordant with (i,j)
      # Cells above-left and below-right
      c_sum <- 0
      d_sum <- 0
      
      # Above and to the left
      if (i > 1 && j > 1) {
        for (k in 1:(i - 1)) {
          for (l in 1:(j - 1)) {
            c_sum <- c_sum + observed[k, l]
          }
        }
      }
      
      # Below and to the right
      if (i < n_row && j < n_col) {
        for (k in (i + 1):n_row) {
          for (l in (j + 1):n_col) {
            c_sum <- c_sum + observed[k, l]
          }
        }
      }
      
      # Above and to the right (discordant)
      if (i > 1 && j < n_col) {
        for (k in 1:(i - 1)) {
          for (l in (j + 1):n_col) {
            d_sum <- d_sum + observed[k, l]
          }
        }
      }
      
      # Below and to the left (discordant)
      if (i < n_row && j > 1) {
        for (k in (i + 1):n_row) {
          for (l in 1:(j - 1)) {
            d_sum <- d_sum + observed[k, l]
          }
        }
      }
      
      C_matrix[i, j] <- c_sum
      D_matrix[i, j] <- d_sum
    }
  }
  
  # ASE1 calculation (variance, not SD yet)
  sum_sq <- 0
  for (i in 1:n_row) {
    for (j in 1:n_col) {
      n_ij <- observed[i, j]
      if (n_ij > 0) {
        term <- Q * C_matrix[i, j] - P * D_matrix[i, j]
        sum_sq <- sum_sq + n_ij * term^2
      }
    }
  }
  
  # SE calculation
  if ((P + Q) > 0 && sum_sq > 0) {
    se_gamma <- (2 / (P + Q)) * sqrt(sum_sq)
  } else {
    se_gamma <- NA
  }
  
  # Z test:  (expected - observed) / SE
  # One-sided p-value following DIGRAM
  if (! is.na(se_gamma) && se_gamma > 0) {
    z_statistic <- (gamma_expected - gamma_observed) / se_gamma
    p_value <- pnorm(z_statistic)  # One-sided p-value
  } else {
    z_statistic <- NA
    p_value <- NA
  }
  
  list(
    gamma_observed = gamma_observed,
    gamma_expected = gamma_expected,
    se_gamma = se_gamma,
    z_statistic = z_statistic,
    p_value = p_value
  )
}