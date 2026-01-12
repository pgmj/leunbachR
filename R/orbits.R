#' Analyze Orbits for Leunbach Model
#'
#' @description
#' Performs orbit analysis for the Leunbach model.  For each total score (orbit),
#' computes the expected distribution of (Test1, Test2) score pairs and 
#' cumulative probabilities for assessing person fit.
#'
#' @param fit A leunbach_ipf object from leunbach_ipf()
#' @param alpha Significance level for critical values (default 0.05)
#' @param verbose Print detailed output
#'
#' @return A list of class "leunbach_orbits" containing:
#'   - orbits: Matrix of expected percentages within each total score
#'   - left_right:  Cumulative probabilities P(X ≤ x | S = s)
#'   - right_left: Cumulative probabilities P(X ≥ x | S = s)
#'   - crit_left: Critical values for left tail by total score
#'   - crit_right: Critical values for right tail by total score
#'   - crit_values: Combined critical values (two-tailed)
#'   - expected_critical: Expected number of persons with significant differences
#'   - variance_expected: Variance of expected critical count
#'   - n_significant: Number of observed persons with significant differences
#'   - orbit_df:  Degrees of freedom for each orbit
#'
#' @details
#' For each total score S = X + Y, the orbit consists of all valid (X, Y) pairs.
#' The expected probability of each cell within an orbit is:
#' P(X = x | S = s) = (gamma_x * delta_y) / sigma_s
#'
#' Cumulative probabilities are computed in both directions:  
#' - Left-to-right: P(X ≤ x | S = s) - tests if Test1 score is unusually low
#' - Right-to-left: P(X ≥ x | S = s) - tests if Test1 score is unusually high
#'
#' @export
analyze_orbits <- function(fit, alpha = 0.05, verbose = FALSE) {
  
  if (!inherits(fit, "leunbach_ipf")) {
    stop("Input must be a leunbach_ipf object")
  }
  
  xmin <- min(fit$test1_scores)
  xmax <- max(fit$test1_scores)
  ymin <- min(fit$test2_scores)
  ymax <- max(fit$test2_scores)
  smin <- min(fit$total_scores)
  smax <- max(fit$total_scores)
  
  n_test1 <- length(fit$test1_scores)
  n_test2 <- length(fit$test2_scores)
  
  expected <- fit$fitted
  smargin <- fit$total_score_freq
  obs <- fit$observed
  
  # Initialize orbit matrices
  orbits <- matrix(0, nrow = n_test1, ncol = n_test2)
  rownames(orbits) <- fit$test1_scores
  colnames(orbits) <- fit$test2_scores
  
  left_right_orbits <- orbits
  right_left_orbits <- orbits
  
  # Compute orbit percentages:  100 * expected / smargin
  for (x in fit$test1_scores) {
    x_idx <- which(fit$test1_scores == x)
    for (y in fit$test2_scores) {
      y_idx <- which(fit$test2_scores == y)
      s <- x + y
      s_char <- as.character(s)
      
      if (! is.na(smargin[s_char]) && smargin[s_char] > 0) {
        orbits[x_idx, y_idx] <- 100.0 * expected[x_idx, y_idx] / smargin[s_char]
      } else {
        orbits[x_idx, y_idx] <- 0
      }
    }
  }
  
  # Compute cumulative probabilities:  Left-to-Right
  # For each orbit s, accumulate from (xmin, s-xmin) towards (xmax, s-xmax)
  for (s in fit$total_scores) {
    for (x in fit$test1_scores) {
      y <- s - x
      if (y >= ymin && y <= ymax) {
        x_idx <- which(fit$test1_scores == x)
        y_idx <- which(fit$test2_scores == y)
        
        if (x == xmin || y == ymax) {
          # First cell in orbit (from left)
          left_right_orbits[x_idx, y_idx] <- orbits[x_idx, y_idx]
        } else {
          # Accumulate from previous cell
          prev_x_idx <- which(fit$test1_scores == (x - 1))
          prev_y_idx <- which(fit$test2_scores == (y + 1))
          if (length(prev_x_idx) > 0 && length(prev_y_idx) > 0) {
            left_right_orbits[x_idx, y_idx] <- 
              left_right_orbits[prev_x_idx, prev_y_idx] + orbits[x_idx, y_idx]
          } else {
            left_right_orbits[x_idx, y_idx] <- orbits[x_idx, y_idx]
          }
        }
      }
    }
  }
  
  # Compute cumulative probabilities: Right-to-Left
  # For each orbit s, accumulate from (xmax, s-xmax) towards (xmin, s-xmin)
  for (s in fit$total_scores) {
    for (y in fit$test2_scores) {
      x <- s - y
      if (x >= xmin && x <= xmax) {
        x_idx <- which(fit$test1_scores == x)
        y_idx <- which(fit$test2_scores == y)
        
        if (x == xmax || y == ymin) {
          # First cell in orbit (from right)
          right_left_orbits[x_idx, y_idx] <- orbits[x_idx, y_idx]
        } else {
          # Accumulate from previous cell
          prev_x_idx <- which(fit$test1_scores == (x + 1))
          prev_y_idx <- which(fit$test2_scores == (y - 1))
          if (length(prev_x_idx) > 0 && length(prev_y_idx) > 0) {
            right_left_orbits[x_idx, y_idx] <- 
              right_left_orbits[prev_x_idx, prev_y_idx] + orbits[x_idx, y_idx]
          } else {
            right_left_orbits[x_idx, y_idx] <- orbits[x_idx, y_idx]
          }
        }
      }
    }
  }
  
  # Identify critical levels for each total score
  crit_left <- numeric(length(fit$total_scores))
  names(crit_left) <- fit$total_scores
  crit_right <- numeric(length(fit$total_scores))
  names(crit_right) <- fit$total_scores
  crit_values <- numeric(length(fit$total_scores))
  names(crit_values) <- fit$total_scores
  
  # Degrees of freedom for each orbit
  orbit_df <- numeric(length(fit$total_scores))
  names(orbit_df) <- fit$total_scores
  
  crit_threshold <- alpha * 100  # Convert to percentage (default 5%)
  
  for (s in fit$total_scores) {
    s_char <- as.character(s)
    crit_left[s_char] <- 0
    crit_right[s_char] <- 0
    
    orbit_size <- 0
    
    # Find critical left value
    for (x in fit$test1_scores) {
      y <- s - x
      if (y >= ymin && y <= ymax) {
        x_idx <- which(fit$test1_scores == x)
        y_idx <- which(fit$test2_scores == y)
        orbit_size <- orbit_size + 1
        
        lr_val <- left_right_orbits[x_idx, y_idx]
        if (lr_val <= crit_threshold && lr_val >= crit_left[s_char]) {
          crit_left[s_char] <- lr_val
        }
      }
    }
    
    # Find critical right value
    for (y in fit$test2_scores) {
      x <- s - y
      if (x >= xmin && x <= xmax) {
        x_idx <- which(fit$test1_scores == x)
        y_idx <- which(fit$test2_scores == y)
        
        rl_val <- right_left_orbits[x_idx, y_idx]
        if (rl_val <= crit_threshold && rl_val >= crit_right[s_char]) {
          crit_right[s_char] <- rl_val
        }
      }
    }
    
    # Combined critical value (two-tailed)
    crit_values[s_char] <- (crit_left[s_char] + crit_right[s_char]) / 100
    
    # Orbit df = orbit_size - 1
    orbit_df[s_char] <- max(0, orbit_size - 1)
  }
  
  # Calculate expected number of persons with significant differences
  expected_critical <- 0
  variance_expected <- 0
  
  for (s in fit$total_scores) {
    s_char <- as.character(s)
    if (!is.na(smargin[s_char]) && smargin[s_char] > 0) {
      expected_critical <- expected_critical + smargin[s_char] * crit_values[s_char]
      variance_expected <- variance_expected + 
        smargin[s_char] * crit_values[s_char] * (1 - crit_values[s_char])
    }
  }
  
  # Count observed persons with significant differences
  sign_table <- obs
  n_significant <- 0
  significant_cells <- list()
  
  for (x in fit$test1_scores) {
    x_idx <- which(fit$test1_scores == x)
    for (y in fit$test2_scores) {
      y_idx <- which(fit$test2_scores == y)
      
      lr_val <- left_right_orbits[x_idx, y_idx]
      rl_val <- right_left_orbits[x_idx, y_idx]
      
      if (lr_val > crit_threshold && rl_val > crit_threshold) {
        sign_table[x_idx, y_idx] <- 0
      } else {
        cell_count <- obs[x_idx, y_idx]
        if (cell_count > 0) {
          n_significant <- n_significant + cell_count
          significant_cells[[length(significant_cells) + 1]] <- list(
            x = x, y = y, n = cell_count,
            p_left = if (lr_val < 50) lr_val / 100 else NA,
            p_right = if (rl_val < 50) rl_val / 100 else NA
          )
        }
      }
    }
  }
  
  # Compute chi-square test for significant differences
  n_total <- fit$n
  p_sign_observed <- n_significant / n_total
  p_sign_expected <- expected_critical / n_total
  sd_expected <- sqrt(variance_expected) / n_total
  
  if (sd_expected > 0) {
    chi2_sign <- ((p_sign_expected - p_sign_observed) / sd_expected)^2
    p_chi2_sign <- pchisq(chi2_sign, df = 1, lower.tail = FALSE)
  } else {
    chi2_sign <- NA
    p_chi2_sign <- NA
  }
  
  # Confidence interval for expected proportion
  ci_low <- 100 * (p_sign_expected - 1.96 * sd_expected)
  ci_high <- 100 * (p_sign_expected + 1.96 * sd_expected)
  
  result <- list(
    orbits = orbits,
    left_right = left_right_orbits,
    right_left = right_left_orbits,
    crit_left = crit_left,
    crit_right = crit_right,
    crit_values = crit_values,
    orbit_df = orbit_df,
    expected_critical = expected_critical,
    variance_expected = variance_expected,
    n_significant = n_significant,
    pct_significant = 100 * n_significant / n_total,
    pct_expected = 100 * expected_critical / n_total,
    ci_low = ci_low,
    ci_high = ci_high,
    chi2_sign = chi2_sign,
    p_chi2_sign = p_chi2_sign,
    sign_table = sign_table,
    significant_cells = significant_cells,
    n_total = n_total,
    alpha = alpha,
    fit = fit
  )
  
  class(result) <- "leunbach_orbits"
  
  if (verbose) {
    print(result)
  }
  
  return(result)
}


#' Print method for leunbach_orbits objects
#' @export
print.leunbach_orbits <- function(x, ...) {
  cat("Leunbach Orbit Analysis\n")
  cat("=======================\n\n")
  
  cat(sprintf("N = %d observations\n", x$n_total))
  cat(sprintf("Significance level:  %.1f%%\n\n", x$alpha * 100))
  
  cat("Critical levels for person fit assessment:\n\n")
  crit_df <- data.frame(
    Score = as.numeric(names(x$crit_left)),
    N = as.numeric(x$fit$total_score_freq),
    Crit_Left = round(x$crit_left, 3),
    Crit_Right = round(x$crit_right, 3),
    Crit_Combined = round(x$crit_values * 100, 3),
    DF = x$orbit_df
  )
  crit_df <- crit_df[crit_df$N > 0, ]
  print(crit_df, row.names = FALSE)
  
  cat("\n")
  cat(sprintf("%d (%.1f%%) persons with significant differences between measurements\n",
              x$n_significant, x$pct_significant))
  cat(sprintf("%.1f (%.1f%%) expected\n\n", 
              x$expected_critical, x$pct_expected))
  
  cat(sprintf("95%% Confidence interval: [%.1f%%, %.1f%%]\n", x$ci_low, x$ci_high))
  
  if (! is.na(x$chi2_sign)) {
    cat(sprintf("Chi-square = %.2f, df = 1, p = %.4f\n", x$chi2_sign, x$p_chi2_sign))
  }
  
  invisible(x)
}


#' Summary method for leunbach_orbits objects
#' @export
summary.leunbach_orbits <- function(object, ...) {
  cat("Leunbach Orbit Analysis - Summary\n")
  cat("==================================\n\n")
  
  # Show significant cells
  if (length(object$significant_cells) > 0) {
    cat("Significant differences between Test1 and Test2:\n")
    cat("  Test1  Test2     N    P(T1<T2)   P(T1>T2)\n")
    cat("  --------------------------------------------\n")
    
    for (cell in object$significant_cells) {
      p_left_str <- if (! is.na(cell$p_left)) sprintf("%8.4f", cell$p_left) else "        "
      p_right_str <- if (!is.na(cell$p_right)) sprintf("%8.4f", cell$p_right) else "        "
      cat(sprintf("  %5d  %5d  %5d  %s  %s\n", 
                  cell$x, cell$y, cell$n, p_left_str, p_right_str))
    }
    cat("\n")
  }
  
  cat(sprintf("%d (%.1f%%) persons with significant differences\n",
              object$n_significant, object$pct_significant))
  cat(sprintf("%.1f (%.1f%%) expected\n\n", 
              object$expected_critical, object$pct_expected))
  
  cat(sprintf("95%% CI: [%.1f%%, %.1f%%]\n", object$ci_low, object$ci_high))
  if (!is.na(object$chi2_sign)) {
    cat(sprintf("Chi-square = %.2f, df = 1, p = %.4f\n", 
                object$chi2_sign, object$p_chi2_sign))
  }
  
  invisible(object)
}


#' Plot method for leunbach_orbits objects
#' @export
plot.leunbach_orbits <- function(x, type = c("orbits", "cumulative", "significant"), ...) {
  type <- match.arg(type)
  
  fit <- x$fit
  
  if (type == "orbits") {
    # Heatmap of orbit probabilities
    image(fit$test1_scores, fit$test2_scores, x$orbits,
          col = colorRampPalette(c("white", "steelblue", "darkblue"))(20),
          xlab = "Test 1 Score", ylab = "Test 2 Score",
          main = "Orbit Distributions (Expected %)")
    
    # Add contour lines for total scores
    contour(fit$test1_scores, fit$test2_scores, 
            outer(fit$test1_scores, fit$test2_scores, "+"),
            add = TRUE, col = "gray50", lty = 2)
    
  } else if (type == "cumulative") {
    oldpar <- par(mfrow = c(1, 2))
    on.exit(par(oldpar))
    
    image(fit$test1_scores, fit$test2_scores, x$left_right,
          col = colorRampPalette(c("white", "orange", "red"))(20),
          xlab = "Test 1 Score", ylab = "Test 2 Score",
          main = "P(Test1 ≤ x | Total) %")
    
    image(fit$test1_scores, fit$test2_scores, x$right_left,
          col = colorRampPalette(c("white", "lightblue", "blue"))(20),
          xlab = "Test 1 Score", ylab = "Test 2 Score",
          main = "P(Test1 ≥ x | Total) %")
    
  } else if (type == "significant") {
    # Show cells with significant differences
    sig_matrix <- x$sign_table
    sig_matrix[sig_matrix == 0] <- NA
    
    image(fit$test1_scores, fit$test2_scores, log1p(sig_matrix),
          col = colorRampPalette(c("lightyellow", "orange", "red"))(20),
          xlab = "Test 1 Score", ylab = "Test 2 Score",
          main = "Significant Differences (log count)")
  }
  
  invisible(x)
}


#' Get orbit distribution for a specific total score
#'
#' @param orbits A leunbach_orbits object
#' @param total_score The total score to examine
#' @return A data frame with the orbit distribution
#' @export
get_orbit <- function(orbits, total_score) {
  if (! inherits(orbits, "leunbach_orbits")) {
    stop("Input must be a leunbach_orbits object")
  }
  
  fit <- orbits$fit
  s <- total_score
  
  xmin <- min(fit$test1_scores)
  xmax <- max(fit$test1_scores)
  ymin <- min(fit$test2_scores)
  ymax <- max(fit$test2_scores)
  
  result <- data.frame(
    test1 = integer(),
    test2 = integer(),
    expected_pct = numeric(),
    cum_left = numeric(),
    cum_right = numeric(),
    observed = integer()
  )
  
  for (x in fit$test1_scores) {
    y <- s - x
    if (y >= ymin && y <= ymax) {
      x_idx <- which(fit$test1_scores == x)
      y_idx <- which(fit$test2_scores == y)
      
      result <- rbind(result, data.frame(
        test1 = x,
        test2 = y,
        expected_pct = round(orbits$orbits[x_idx, y_idx], 2),
        cum_left = round(orbits$left_right[x_idx, y_idx], 2),
        cum_right = round(orbits$right_left[x_idx, y_idx], 2),
        observed = fit$observed[x_idx, y_idx]
      ))
    }
  }
  
  result
}