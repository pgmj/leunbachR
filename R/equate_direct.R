#' Estimate person parameter (theta) for a given score
#'
#' @description
#' Estimates the person parameter theta such that the expected score equals
#' the observed score. 
#'
#' @param score The observed score
#' @param gamma Score parameters (named vector)
#' @param score_min Minimum observed score
#' @param score_max Maximum observed score
#' @param method Optimization method: "optimize" (default) uses stats::optimize() 
#'        with Brent's method, "newton" uses custom Newton-Raphson with bisection fallback
#' @param tol Convergence tolerance
#' @param max_iter Maximum iterations (only for method = "newton")
#'
#' @return Estimated theta value
#' @keywords internal
estimate_person_parameter <- function(score, gamma, score_min, score_max,
                                      method = c("optimize", "newton"),
                                      tol = 1e-10, max_iter = 500) {
  
  method <- match.arg(method)
  
  # Get valid score range
  scores <- as.numeric(names(gamma)[gamma > 0])
  scores <- scores[scores >= score_min & scores <= score_max]
  
  if (length(scores) < 2) return(NA)
  
  # Boundary scores should be handled by the calling function
  # This function estimates theta for intermediate scores only
  if (score <= score_min || score >= score_max) {
    return(NA)
  }
  
  if (method == "optimize") {
    estimate_theta_optimize(score, gamma, score_min, score_max, tol)
  } else {
    estimate_theta_newton(score, gamma, score_min, score_max, tol, max_iter)
  }
}


#' Estimate theta using stats::optimize() with Brent's method
#'
#' @keywords internal
estimate_theta_optimize <- function(score, gamma, score_min, score_max, tol = 1e-10) {
  
  # Objective function: squared difference between expected and target score
  objective <- function(log_theta) {
    theta <- exp(log_theta)
    expected <- calculate_expected_score(theta, gamma, score_min, score_max)
    if (is.na(expected)) return(1e10)
    (expected - score)^2
  }
  
  # Use optimize() for 1D bounded optimization (Brent's method)
  # Search on log scale for better numerical properties
  result <- tryCatch({
    opt <- optimize(
      f = objective,
      interval = c(-25, 25),  # theta from ~1e-11 to ~1e11
      tol = tol
    )
    exp(opt$minimum)
  }, error = function(e) {
    NA
  })
  
  # Verify the result
  if (! is.na(result)) {
    expected <- calculate_expected_score(result, gamma, score_min, score_max)
    if (is.na(expected) || abs(expected - score) > 0.01) {
      # Optimization failed to find good solution
      result <- NA
    }
  }
  
  return(result)
}


#' Estimate theta using Newton-Raphson with bisection fallback
#'
#' @keywords internal
estimate_theta_newton <- function(score, gamma, score_min, score_max,
                                  tol = 1e-7, max_iter = 500) {
  
  # Better initial theta estimate using logit-like transformation
  score_range <- score_max - score_min
  score_centered <- (score - score_min) / score_range
  score_centered <- max(0.001, min(0.999, score_centered))
  theta <- score_centered / (1 - score_centered)
  
  # Try multiple starting points if needed
  starting_thetas <- c(
    theta,
    theta * 0.1,
    theta * 10,
    exp((score - score_min) / 2),
    1.0
  )
  
  best_theta <- NA
  best_error <- Inf
  
  for (start_theta in starting_thetas) {
    theta_try <- start_theta
    
    for (iter in 1:max_iter) {
      result <- expected_score_and_deriv(theta_try, gamma, score_min, score_max)
      expected <- result$expected
      deriv <- result$deriv
      
      if (is.na(expected) || is.na(deriv)) {
        break
      }
      
      diff <- expected - score
      
      if (abs(diff) < tol) {
        if (abs(diff) < best_error) {
          best_theta <- theta_try
          best_error <- abs(diff)
        }
        break
      }
      
      if (abs(deriv) < 1e-12) {
        step <- diff * 0.1
      } else {
        step <- diff / deriv
      }
      
      max_step <- theta_try * 2
      if (abs(step) > max_step) {
        step <- sign(step) * max_step
      }
      
      theta_new <- theta_try - step
      theta_new <- max(1e-10, min(1e10, theta_new))
      
      if (iter > 50 && abs(theta_new - theta_try) / theta_try < 1e-10) {
        break
      }
      
      theta_try <- theta_new
    }
    
    if (! is.na(theta_try) && theta_try > 0) {
      result <- expected_score_and_deriv(theta_try, gamma, score_min, score_max)
      if (! is.na(result$expected)) {
        error <- abs(result$expected - score)
        if (error < best_error) {
          best_theta <- theta_try
          best_error <- error
        }
      }
    }
  }
  
  # Bisection fallback if Newton-Raphson failed
  if (is.na(best_theta) || best_error > 0.01) {
    bisect_result <- estimate_theta_bisection(score, gamma, score_min, score_max)
    if (!is.na(bisect_result)) {
      expected <- calculate_expected_score(bisect_result, gamma, score_min, score_max)
      if (!is.na(expected)) {
        error <- abs(expected - score)
        if (error < best_error) {
          best_theta <- bisect_result
          best_error <- error
        }
      }
    }
  }
  
  return(best_theta)
}


#' Estimate theta using bisection method
#'
#' @keywords internal
estimate_theta_bisection <- function(score, gamma, score_min, score_max,
                                     tol = 1e-6, max_iter = 100) {
  
  theta_low <- 1e-10
  theta_high <- 1e10
  
  result_low <- calculate_expected_score(theta_low, gamma, score_min, score_max)
  result_high <- calculate_expected_score(theta_high, gamma, score_min, score_max)
  
  if (is.na(result_low) || is.na(result_high)) {
    return(NA)
  }
  
  if (result_low > score || result_high < score) {
    theta_try <- 1e-5
    while (theta_try < 1e8) {
      result <- calculate_expected_score(theta_try, gamma, score_min, score_max)
      if (!is.na(result)) {
        if (result < score) {
          theta_low <- theta_try
        } else {
          theta_high <- theta_try
          break
        }
      }
      theta_try <- theta_try * 10
    }
  }
  
  for (iter in 1:max_iter) {
    theta_mid <- sqrt(theta_low * theta_high)
    
    expected <- calculate_expected_score(theta_mid, gamma, score_min, score_max)
    
    if (is.na(expected)) {
      return(NA)
    }
    
    diff <- expected - score
    
    if (abs(diff) < tol) {
      return(theta_mid)
    }
    
    if (diff < 0) {
      theta_low <- theta_mid
    } else {
      theta_high <- theta_mid
    }
    
    if (abs(log(theta_high) - log(theta_low)) < 1e-10) {
      return(theta_mid)
    }
  }
  
  return(sqrt(theta_low * theta_high))
}


#' Calculate expected score given theta
#'
#' @keywords internal
calculate_expected_score <- function(theta, gamma, score_min, score_max) {
  
  scores <- as.numeric(names(gamma))
  valid <- scores >= score_min & scores <= score_max & gamma > 0
  scores <- scores[valid]
  gam <- as.numeric(gamma[valid])
  
  if (length(scores) < 2 || theta <= 0) {
    return(NA)
  }
  
  log_theta <- log(theta)
  log_weights <- log(gam) + scores * log_theta
  
  max_log_weight <- max(log_weights)
  weights <- exp(log_weights - max_log_weight)
  
  sum_weights <- sum(weights)
  
  if (sum_weights == 0 || !is.finite(sum_weights)) {
    return(NA)
  }
  
  sum_x_weights <- sum(scores * weights)
  expected <- sum_x_weights / sum_weights
  
  return(expected)
}


#' Calculate expected score and its derivative with respect to theta
#'
#' @keywords internal
expected_score_and_deriv <- function(theta, gamma, score_min, score_max) {
  
  scores <- as.numeric(names(gamma))
  valid <- scores >= score_min & scores <= score_max & gamma > 0
  scores <- scores[valid]
  gam <- as.numeric(gamma[valid])
  
  if (length(scores) < 2 || theta <= 0) {
    return(list(expected = NA, deriv = NA))
  }
  
  log_theta <- log(theta)
  log_weights <- log(gam) + scores * log_theta
  
  max_log_weight <- max(log_weights)
  weights <- exp(log_weights - max_log_weight)
  
  sum_weights <- sum(weights)
  
  if (sum_weights == 0 || !is.finite(sum_weights)) {
    return(list(expected = NA, deriv = NA))
  }
  
  sum_x_weights <- sum(scores * weights)
  sum_x2_weights <- sum(scores^2 * weights)
  
  expected <- sum_x_weights / sum_weights
  
  expected_x2 <- sum_x2_weights / sum_weights
  variance <- expected_x2 - expected^2
  variance <- max(0, variance)
  
  deriv <- variance / theta
  
  return(list(expected = expected, deriv = deriv))
}


#' Calculate expected (true) score on target test given theta
#'
#' @keywords internal
calculate_true_score <- function(theta, gamma, score_min, score_max) {
  calculate_expected_score(theta, gamma, score_min, score_max)
}


#' Equate Scores Between Tests using Leunbach Model
#'
#' @description
#' Creates equating tables to convert scores from one test to another using
#' the estimated score parameters from the Leunbach model.
#'
#' @param fit A leunbach_ipf object from leunbach_ipf()
#' @param direction Direction of equating:  "1to2" (Test1 to Test2) or "2to1" (Test2 to Test1)
#' @param method Optimization method for person parameter estimation:  
#'        "optimize" (default) uses stats::optimize() with Brent's method,
#'        "newton" uses custom Newton-Raphson with bisection fallback
#' @param verbose Print detailed output
#'
#' @return A list of class "leunbach_equating" containing:
#'   - equating_table: Data frame with original scores, theta, expected equated scores, and rounded scores
#'   - direction: Direction of equating
#'   - method: Optimization method used
#'   - fit: Original leunbach_ipf object
#'
#' @export
leunbach_equate <- function(fit, direction = c("1to2", "2to1"),
                            method = c("optimize", "newton"),
                            verbose = FALSE) {
  
  if (! inherits(fit, "leunbach_ipf")) {
    stop("Input must be a leunbach_ipf object")
  }
  
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
    source_name <- "Test1"
    target_name <- "Test2"
  } else {
    source_scores <- fit$test2_scores
    source_gamma <- fit$delta
    source_min <- fit$ymin
    source_max <- fit$ymax
    target_gamma <- fit$gamma
    target_min <- fit$xmin
    target_max <- fit$xmax
    source_name <- "Test2"
    target_name <- "Test1"
  }
  
  n_scores <- length(source_scores)
  theta_values <- rep(NA_real_, n_scores)
  expected_score <- rep(NA_real_, n_scores)
  rounded_score <- rep(NA_integer_, n_scores)
  names(theta_values) <- source_scores
  names(expected_score) <- source_scores
  names(rounded_score) <- source_scores
  
  # For scores below source_min, equated score is target_min (no theta)
  for (x in source_scores[source_scores < source_min]) {
    x_char <- as.character(x)
    expected_score[x_char] <- target_min
    rounded_score[x_char] <- target_min
  }
  
  # For scores above source_max, equated score is target_max (no theta)
  for (x in source_scores[source_scores > source_max]) {
    x_char <- as.character(x)
    expected_score[x_char] <- target_max
    rounded_score[x_char] <- target_max
  }
  
  # Handle minimum boundary score:  theta = exp(-5), equated = target_min
  source_min_char <- as.character(source_min)
  theta_values[source_min_char] <- exp(-5)
  expected_score[source_min_char] <- target_min
  rounded_score[source_min_char] <- target_min
  
  # Handle maximum boundary score:  theta = exp(5), equated = target_max
  source_max_char <- as.character(source_max)
  theta_values[source_max_char] <- exp(5)
  expected_score[source_max_char] <- target_max
  rounded_score[source_max_char] <- target_max
  
  # Estimate theta for intermediate scores (source_min < x < source_max)
  for (x in source_scores[source_scores > source_min & source_scores < source_max]) {
    x_char <- as.character(x)
    
    # Estimate theta for this score
    theta <- estimate_person_parameter(x, source_gamma, source_min, source_max,
                                       method = method)
    
    if (! is.na(theta) && theta > 0) {
      # Bound theta:  exp(-5) to exp(5)
      theta <- max(exp(-5), min(exp(5), theta))
      theta_values[x_char] <- theta
      
      # Calculate true score on target test
      true_score <- calculate_true_score(theta, target_gamma, target_min, target_max)
      expected_score[x_char] <- true_score
      rounded_score[x_char] <- round(true_score)
    }
  }
  
  equating_table <- data.frame(
    source_score = source_scores,
    theta = theta_values,
    expected_target = expected_score,
    rounded_target = as.integer(rounded_score)
  )
  
  colnames(equating_table) <- c(source_name, "Theta", 
                                paste0("Expected_", target_name),
                                paste0("Rounded_", target_name))
  
  if (verbose) {
    cat(sprintf("Equating %s to %s (method: %s)\n", source_name, target_name, method))
    cat("==========================================\n\n")
    print(equating_table[! is.na(equating_table[, 3]), ], row.names = FALSE)
  }
  
  result <- list(
    equating_table = equating_table,
    direction = direction,
    method = method,
    source_name = source_name,
    target_name = target_name,
    source_min = source_min,
    source_max = source_max,
    target_min = target_min,
    target_max = target_max,
    fit = fit
  )
  
  class(result) <- "leunbach_equating"
  return(result)
}


#' Print method for leunbach_equating objects
#' @export
print.leunbach_equating <- function(x, ...) {
  cat(sprintf("Leunbach Equating:  %s to %s\n", x$source_name, x$target_name))
  cat(sprintf("Method: %s\n", x$method))
  cat("==========================================\n\n")
  
  tab <- x$equating_table
  # Filter to valid range (source_min to source_max)
  valid_idx <- tab[, 1] >= x$source_min & tab[, 1] <= x$source_max
  tab <- tab[valid_idx, ]
  
  # Format for display
  display_tab <- data.frame(
    Score = tab[, 1],
    Log_Theta = ifelse(is.na(tab[, 2]) | tab[, 2] <= 0, 
                       "      NA", 
                       sprintf("%8.4f", log(tab[, 2]))),
    Expected = sprintf("%6.2f", tab[, 3]),
    Rounded = tab[, 4]
  )
  
  colnames(display_tab) <- c(x$source_name, "Theta", 
                             paste0("Expected_", x$target_name),
                             paste0("Rounded_", x$target_name))
  
  print(display_tab, row.names = FALSE)
  
  invisible(x)
}
