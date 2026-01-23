#' Indirect Equating via an Anchor Test
#'
#' @description
#' Performs indirect equating from Test A to Test C via an anchor Test B.
#' This chains two direct equatings:  A → B and B → C.
#'
#' @param fit_ab A leunbach_ipf object for the A-B equating (Tests A and B)
#' @param fit_bc A leunbach_ipf object for the B-C equating (Tests B and C)
#' @param direction_ab Direction for A-B equating:  "1to2" or "2to1"
#' @param direction_bc Direction for B-C equating: "1to2" or "2to1"
#' @param method Optimization method: "optimize" (default) or "newton"
#' @param verbose Print detailed output
#'
#' @return A list of class "leunbach_indirect" containing:
#'   - equating_table: Data frame with source scores, expected equated scores, and rounded scores
#'   - eq_ab: Direct equating object for A → B
#'   - eq_bc: Direct equating object for B → C
#'   - fit_ab: Original leunbach_ipf object for A-B
#'   - fit_bc:  Original leunbach_ipf object for B-C
#'
#' @details
#' Indirect equating works by chaining two direct equatings: 
#' 1. For a score x on Test A, find the expected score on Test B (typically non-integer)
#' 2. Find expected Test C scores for the integer B scores below and above
#' 3. Interpolate to get the expected Test C score for the non-integer B score
#' 4. Round to get the equated integer score
#'
#' @examples
#' # Fit models for A-B and B-C
#' fit_ab <- leunbach_ipf(data_ab)
#' fit_bc <- leunbach_ipf(data_bc)
#'
#' # Indirect equating:  Test1 of fit_ab → Test2 of fit_ab → Test2 of fit_bc
#' indirect <- leunbach_indirect_equate(fit_ab, fit_bc, 
#'                                       direction_ab = "1to2", 
#'                                       direction_bc = "1to2")
#' print(indirect)
#'
#' @export
leunbach_indirect_equate <- function(fit_ab, fit_bc,
                                     direction_ab = c("1to2", "2to1"),
                                     direction_bc = c("1to2", "2to1"),
                                     method = c("optimize", "newton"),
                                     verbose = FALSE) {
  
  if (! inherits(fit_ab, "leunbach_ipf")) {
    stop("fit_ab must be a leunbach_ipf object")
  }
  if (!inherits(fit_bc, "leunbach_ipf")) {
    stop("fit_bc must be a leunbach_ipf object")
  }
  
  direction_ab <- match.arg(direction_ab)
  direction_bc <- match.arg(direction_bc)
  method <- match.arg(method)
  
  # Get direct equatings
  eq_ab <- leunbach_equate(fit_ab, direction = direction_ab, method = method)
  eq_bc <- leunbach_equate(fit_bc, direction = direction_bc, method = method)
  
  # Verify that the anchor test (target of AB, source of BC) matches
  # The target of eq_ab should align with the source of eq_bc
  anchor_min_ab <- eq_ab$target_min
  anchor_max_ab <- eq_ab$target_max
  anchor_min_bc <- eq_bc$source_min
  anchor_max_bc <- eq_bc$source_max
  
  if (verbose) {
    cat("Indirect Equating\n")
    cat("=================\n\n")
    cat(sprintf("Step 1: %s → %s (anchor test)\n", "Test A", "Test B"))
    cat(sprintf("        Anchor range from eq_ab: %d to %d\n", anchor_min_ab, anchor_max_ab))
    cat(sprintf("Step 2: %s → %s\n", "Test B", "Test C"))
    cat(sprintf("        Anchor range from eq_bc: %d to %d\n", anchor_min_bc, anchor_max_bc))
    cat(sprintf("Result: %s → %s (indirect)\n\n", "Test A", "Test C"))
  }
  
  # Perform indirect equating using interpolation
  indirect_table <- compute_indirect_equating(eq_ab, eq_bc)
  
  if (verbose) {
    print(indirect_table[! is.na(indirect_table$expected), ], row.names = FALSE)
  }
  
  result <- list(
    equating_table = indirect_table,
    eq_ab = eq_ab,
    eq_bc = eq_bc,
    fit_ab = fit_ab,
    fit_bc = fit_bc,
    direction_ab = direction_ab,
    direction_bc = direction_bc,
    method = method,
    source_name = "Test A",
    anchor_name = "Test B",
    target_name = "Test C",
    source_min = eq_ab$source_min,
    source_max = eq_ab$source_max,
    anchor_min = max(anchor_min_ab, anchor_min_bc),
    anchor_max = min(anchor_max_ab, anchor_max_bc),
    target_min = eq_bc$target_min,
    target_max = eq_bc$target_max
  )
  
  class(result) <- "leunbach_indirect"
  return(result)
}


#' Compute indirect equating table using interpolation
#'
#' @keywords internal
compute_indirect_equating <- function(eq_ab, eq_bc) {
  
  # Get source scores (Test A)
  source_scores <- eq_ab$equating_table[, 1]
  theta_ab <- eq_ab$equating_table[, 2]        # Theta values from A->B
  expected_anchor <- eq_ab$equating_table[, 3] # Expected B scores (continuous)
  
  # Get anchor -> target equating table
  anchor_scores <- eq_bc$equating_table[, 1]
  theta_bc <- eq_bc$equating_table[, 2]        # Theta values from B->C
  expected_target <- eq_bc$equating_table[, 3] # Expected C scores for integer B
  
  # Create lookup for B -> C (integer B scores only)
  bc_lookup <- setNames(expected_target, anchor_scores)
  
  # Source score range
  source_min <- eq_ab$source_min
  source_max <- eq_ab$source_max
  
  # Anchor score range (use intersection of both equatings)
  anchor_min <- max(eq_ab$target_min, eq_bc$source_min)
  anchor_max <- min(eq_ab$target_max, eq_bc$source_max)
  
  # Target score range
  target_min <- eq_bc$target_min
  target_max <- eq_bc$target_max
  
  # Initialize result vectors
  n_scores <- length(source_scores)
  theta_indirect <- rep(NA_real_, n_scores)
  expected_indirect <- rep(NA_real_, n_scores)
  rounded_indirect <- rep(NA_integer_, n_scores)
  
  for (i in seq_along(source_scores)) {
    x <- source_scores[i]
    
    # Skip if outside observed source range
    if (x < source_min || x > source_max) {
      next
    }
    
    # Get theta and expected anchor score (B) for this source score (A)
    theta_a <- theta_ab[i]
    exp_b <- expected_anchor[i]
    
    if (is.na(exp_b)) {
      next
    }
    
    # Store the theta from source test
    theta_indirect[i] <- theta_a
    
    # Handle boundary cases - strictly outside the anchor range
    if (exp_b < anchor_min) {
      # Extrapolate below: use slope at anchor_min
      c_at_min <- bc_lookup[as.character(anchor_min)]
      c_at_min_plus_1 <- bc_lookup[as.character(anchor_min + 1)]
      if (!is.na(c_at_min) && !is.na(c_at_min_plus_1)) {
        slope <- c_at_min_plus_1 - c_at_min
        expected_indirect[i] <- c_at_min + slope * (exp_b - anchor_min)
        expected_indirect[i] <- max(target_min, expected_indirect[i])
      } else if (!is.na(c_at_min)) {
        expected_indirect[i] <- c_at_min
      } else {
        expected_indirect[i] <- target_min
      }
      rounded_indirect[i] <- round(expected_indirect[i])
      next
    }
    
    if (exp_b > anchor_max) {
      # Extrapolate above: use slope at anchor_max
      c_at_max <- bc_lookup[as.character(anchor_max)]
      c_at_max_minus_1 <- bc_lookup[as.character(anchor_max - 1)]
      if (!is.na(c_at_max) && !is.na(c_at_max_minus_1)) {
        slope <- c_at_max - c_at_max_minus_1
        expected_indirect[i] <- c_at_max + slope * (exp_b - anchor_max)
        expected_indirect[i] <- min(target_max, expected_indirect[i])
      } else if (!is.na(c_at_max)) {
        expected_indirect[i] <- c_at_max
      } else {
        expected_indirect[i] <- target_max
      }
      rounded_indirect[i] <- round(expected_indirect[i])
      next
    }
    
    # Interpolate between integer anchor scores
    b_low <- floor(exp_b)
    b_high <- ceiling(exp_b)
    
    # Get expected target scores for integer anchor scores
    c_low <- bc_lookup[as.character(b_low)]
    c_high <- bc_lookup[as.character(b_high)]
    
    if (is.na(c_low) || is.na(c_high)) {
      # Try to find nearest valid anchor scores
      c_low <- find_nearest_equated(b_low, bc_lookup, direction = "down")
      c_high <- find_nearest_equated(b_high, bc_lookup, direction = "up")
      
      if (is.na(c_low) || is.na(c_high)) {
        next
      }
    }
    
    # Linear interpolation
    if (b_low == b_high) {
      # exp_b is an integer
      expected_indirect[i] <- c_low
    } else {
      frac <- exp_b - b_low
      expected_indirect[i] <- c_low + frac * (c_high - c_low)
    }
    
    # Bound to target range and round
    expected_indirect[i] <- max(target_min, min(target_max, expected_indirect[i]))
    rounded_indirect[i] <- round(expected_indirect[i])
  }
  
  # Create output table with named columns
  data.frame(
    source = source_scores,
    theta = theta_indirect,
    expected = expected_indirect,
    rounded = as.integer(rounded_indirect)
  )
}

#' Find nearest valid equated score
#'
#' @keywords internal
find_nearest_equated <- function(score, lookup, direction = c("down", "up")) {
  direction <- match.arg(direction)
  
  available <- as.numeric(names(lookup)[! is.na(lookup)])
  
  if (length(available) == 0) {
    return(NA)
  }
  
  if (direction == "down") {
    valid <- available[available <= score]
    if (length(valid) == 0) return(NA)
    nearest <- max(valid)
  } else {
    valid <- available[available >= score]
    if (length(valid) == 0) return(NA)
    nearest <- min(valid)
  }
  
  return(lookup[as.character(nearest)])
}


#' Print method for leunbach_indirect objects
#' @export
print.leunbach_indirect <- function(x, ...) {
  cat("Leunbach Indirect Equating\n")
  cat("==========================\n\n")
  cat(sprintf("Path: %s -> %s -> %s\n", x$source_name, x$anchor_name, x$target_name))
  cat(sprintf("Method: %s\n\n", x$method))
  
  cat(sprintf("Source (%s) range:  %d to %d\n", x$source_name, x$source_min, x$source_max))
  cat(sprintf("Anchor (%s) range: %d to %d\n", x$anchor_name, x$anchor_min, x$anchor_max))
  cat(sprintf("Target (%s) range: %d to %d\n\n", x$target_name, x$target_min, x$target_max))
  
  tab <- x$equating_table
  # Filter to valid range
  valid_idx <- tab$source >= x$source_min & tab$source <= x$source_max
  tab <- tab[valid_idx, ]
  
  # Format for display
  display_tab <- data.frame(
    Score = tab$source,
    Log_Theta = ifelse(is.na(tab$theta) | tab$theta <= 0, 
                       "      NA", 
                       sprintf("%8.4f", log(tab$theta))),
    Expected = sprintf("%6.2f", tab$expected),
    Rounded = tab$rounded
  )
  
  colnames(display_tab) <- c(x$source_name, "Theta",
                             paste0("Expected_", x$target_name), 
                             paste0("Rounded_", x$target_name))
  
  print(display_tab, row.names = FALSE)
  
  invisible(x)
}


#' Bootstrap for Indirect Equating
#'
#' @description
#' Performs parametric bootstrapping for indirect equating to assess
#' standard errors of the indirect equated scores.
#'
#' @param fit_ab A leunbach_ipf object for the A-B equating
#' @param fit_bc A leunbach_ipf object for the B-C equating
#' @param direction_ab Direction for A-B equating:  "1to2" or "2to1"
#' @param direction_bc Direction for B-C equating: "1to2" or "2to1"
#' @param nsim Number of bootstrap samples (default: 1000)
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param see_type Type of SEE calculation:  "rounded" or "expected"
#' @param method Optimization method: "optimize" (default) or "newton"
#' @param parallel Use parallel processing if mirai is available (default: TRUE)
#' @param n_cores Number of cores for parallel processing
#' @param verbose Print progress messages
#' @param seed Random seed for reproducibility
#'
#' @return A list of class "leunbach_indirect_bootstrap" containing:
#'   - indirect_eq:  The observed indirect equating object
#'   - Bootstrap results and standard errors
#'   - Bootstrap p-values for LR and Gamma tests for both equatings
#'
#' @export
leunbach_indirect_bootstrap <- function(fit_ab, fit_bc,
                                        direction_ab = c("1to2", "2to1"),
                                        direction_bc = c("1to2", "2to1"),
                                        nsim = 1000,
                                        conf_level = 0.95,
                                        see_type = c("rounded", "expected"),
                                        method = c("optimize", "newton"),
                                        parallel = TRUE,
                                        n_cores = NULL,
                                        verbose = FALSE,
                                        seed = NULL) {
  
  if (! inherits(fit_ab, "leunbach_ipf")) {
    stop("fit_ab must be a leunbach_ipf object")
  }
  if (!inherits(fit_bc, "leunbach_ipf")) {
    stop("fit_bc must be a leunbach_ipf object")
  }
  
  direction_ab <- match.arg(direction_ab)
  direction_bc <- match.arg(direction_bc)
  see_type <- match.arg(see_type)
  method <- match.arg(method)
  
  # Check for mirai availability
  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
  
  if (parallel && ! use_parallel) {
    message("Install 'mirai' package for parallel processing:  install.packages('mirai')")
    message("Running sequentially...")
  }
  
  if (use_parallel) {
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
    n_cores <- min(n_cores, nsim)
  }
  
  if (! is.null(seed)) {
    set.seed(seed)
  }
  
  boot_seeds <- sample.int(.Machine$integer.max, nsim)
  
  # Get observed values
  lr_ab_observed <- fit_ab$g_sq
  lr_bc_observed <- fit_bc$g_sq
  gk_gamma_z_ab_observed <- fit_ab$gk_gamma_z
  gk_gamma_z_bc_observed <- fit_bc$gk_gamma_z
  
  # Get observed indirect equating
  indirect_eq <- leunbach_indirect_equate(fit_ab, fit_bc,
                                          direction_ab = direction_ab,
                                          direction_bc = direction_bc,
                                          method = method,
                                          verbose = FALSE)
  
  observed_expected <- indirect_eq$equating_table$expected
  observed_rounded <- indirect_eq$equating_table$rounded
  source_scores <- indirect_eq$equating_table$source
  n_scores <- length(source_scores)
  
  if (verbose) {
    cat("Parametric Bootstrap for Indirect Equating\n")
    cat("===========================================\n\n")
    cat(sprintf("Path: %s -> %s -> %s\n", "Test A", "Test B", "Test C"))
    cat(sprintf("Optimization method: %s\n", method))
    if (use_parallel) {
      cat(sprintf("Running %d bootstrap samples using %d cores...\n\n", nsim, n_cores))
    } else {
      cat(sprintf("Running %d bootstrap samples sequentially...\n\n", nsim))
    }
  }
  
  # Prepare bootstrap data
  boot_data_list <- list(
    # Fit AB data
    fit_ab_data = list(
      test1_scores = fit_ab$test1_scores,
      test2_scores = fit_ab$test2_scores,
      total_scores = fit_ab$total_scores,
      total_score_freq = fit_ab$total_score_freq,
      gamma = fit_ab$gamma,
      delta = fit_ab$delta,
      sigma = fit_ab$sigma,
      xmin = fit_ab$xmin,
      xmax = fit_ab$xmax,
      ymin = fit_ab$ymin,
      ymax = fit_ab$ymax
    ),
    # Fit BC data
    fit_bc_data = list(
      test1_scores = fit_bc$test1_scores,
      test2_scores = fit_bc$test2_scores,
      total_scores = fit_bc$total_scores,
      total_score_freq = fit_bc$total_score_freq,
      gamma = fit_bc$gamma,
      delta = fit_bc$delta,
      sigma = fit_bc$sigma,
      xmin = fit_bc$xmin,
      xmax = fit_bc$xmax,
      ymin = fit_bc$ymin,
      ymax = fit_bc$ymax
    ),
    direction_ab = direction_ab,
    direction_bc = direction_bc,
    method = method,
    n_scores = n_scores
  )
  
  # Run bootstrap
  if (use_parallel) {
    boot_results <- run_indirect_bootstrap_parallel(
      nsim = nsim,
      boot_seeds = boot_seeds,
      boot_data_list = boot_data_list,
      n_cores = n_cores,
      verbose = verbose
    )
  } else {
    boot_results <- run_indirect_bootstrap_sequential(
      nsim = nsim,
      boot_seeds = boot_seeds,
      boot_data_list = boot_data_list,
      verbose = verbose
    )
  }
  
  # Extract results
  lr_ab_bootstrap <- boot_results$lr_ab_bootstrap
  lr_bc_bootstrap <- boot_results$lr_bc_bootstrap
  gk_gamma_z_ab_bootstrap <- boot_results$gk_gamma_z_ab_bootstrap
  gk_gamma_z_bc_bootstrap <- boot_results$gk_gamma_z_bc_bootstrap
  boot_expected <- boot_results$boot_expected
  boot_rounded <- boot_results$boot_rounded
  boot_failed <- boot_results$boot_failed
  
  colnames(boot_expected) <- source_scores
  colnames(boot_rounded) <- source_scores
  colnames(boot_failed) <- source_scores
  
  # Compute bootstrap p-values for LR tests
  valid_lr_ab <- ! is.na(lr_ab_bootstrap)
  valid_lr_bc <- !is.na(lr_bc_bootstrap)
  
  if (sum(valid_lr_ab) > 0) {
    n_significant_lr_ab <- sum(lr_ab_bootstrap[valid_lr_ab] >= lr_ab_observed)
    p_lr_ab <- n_significant_lr_ab / sum(valid_lr_ab)
  } else {
    n_significant_lr_ab <- NA
    p_lr_ab <- NA
  }
  
  if (sum(valid_lr_bc) > 0) {
    n_significant_lr_bc <- sum(lr_bc_bootstrap[valid_lr_bc] >= lr_bc_observed)
    p_lr_bc <- n_significant_lr_bc / sum(valid_lr_bc)
  } else {
    n_significant_lr_bc <- NA
    p_lr_bc <- NA
  }
  
  # Compute bootstrap p-values for Gamma tests (one-sided)
  valid_gamma_ab <- !is.na(gk_gamma_z_ab_bootstrap)
  valid_gamma_bc <- !is.na(gk_gamma_z_bc_bootstrap)
  
  if (sum(valid_gamma_ab) > 0 && !is.na(gk_gamma_z_ab_observed)) {
    n_significant_gamma_ab <- sum(gk_gamma_z_ab_bootstrap[valid_gamma_ab] >= gk_gamma_z_ab_observed)
    p_gamma_ab <- n_significant_gamma_ab / sum(valid_gamma_ab)
  } else {
    n_significant_gamma_ab <- NA
    p_gamma_ab <- NA
  }
  
  if (sum(valid_gamma_bc) > 0 && !is.na(gk_gamma_z_bc_observed)) {
    n_significant_gamma_bc <- sum(gk_gamma_z_bc_bootstrap[valid_gamma_bc] >= gk_gamma_z_bc_observed)
    p_gamma_bc <- n_significant_gamma_bc / sum(valid_gamma_bc)
  } else {
    n_significant_gamma_bc <- NA
    p_gamma_bc <- NA
  }
  
  # Confidence intervals
  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)
  
  ci_lower <- apply(boot_expected, 2, quantile, probs = probs[1], na.rm = TRUE)
  ci_upper <- apply(boot_expected, 2, quantile, probs = probs[2], na.rm = TRUE)
  
  # Standard Error of Equating
  if (see_type == "rounded") {
    see <- calculate_bootstrap_see(boot_rounded)
  } else {
    see <- calculate_bootstrap_see(boot_expected)
  }
  
  # Error frequencies
  error_freq <- compute_error_frequencies(boot_rounded, observed_rounded, source_scores)
  
  # Proportion failed
  prop_failed <- colMeans(boot_failed, na.rm = TRUE) * 100
  
  # Valid range for average SEE
  source_min <- indirect_eq$source_min
  source_max <- indirect_eq$source_max
  valid_range <- source_scores >= source_min & source_scores <= source_max
  
  avg_see <- mean(see[valid_range], na.rm = TRUE)
  
  # Count valid bootstrap samples
  n_valid <- sum(rowSums(! is.na(boot_expected)) > 0)
  
  if (verbose) {
    cat(sprintf("\nBootstrap complete.\n"))
    cat(sprintf("  Valid samples: %d of %d\n", n_valid, nsim))
    cat(sprintf("  Bootstrap p-value for LR test (A-B): %.3f\n", p_lr_ab))
    cat(sprintf("  Bootstrap p-value for LR test (B-C): %.3f\n", p_lr_bc))
    cat(sprintf("  Bootstrap p-value for Gamma test (A-B): %.3f\n", p_gamma_ab))
    cat(sprintf("  Bootstrap p-value for Gamma test (B-C): %.3f\n", p_gamma_bc))
    cat(sprintf("  Average SEE: %.2f\n", avg_see))
  }
  
  result <- list(
    nsim = nsim,
    n_valid = n_valid,
    conf_level = conf_level,
    see_type = see_type,
    method = method,
    parallel = use_parallel,
    n_cores = if (use_parallel) n_cores else 1L,
    # Fit A-B results
    lr_ab_observed = lr_ab_observed,
    df_ab = fit_ab$df,
    lr_ab_bootstrap = lr_ab_bootstrap,
    p_lr_ab = p_lr_ab,
    gk_gamma_z_ab_observed = gk_gamma_z_ab_observed,
    gk_gamma_z_ab_bootstrap = gk_gamma_z_ab_bootstrap,
    p_gamma_ab = p_gamma_ab,
    # Fit B-C results
    lr_bc_observed = lr_bc_observed,
    df_bc = fit_bc$df,
    lr_bc_bootstrap = lr_bc_bootstrap,
    p_lr_bc = p_lr_bc,
    gk_gamma_z_bc_observed = gk_gamma_z_bc_observed,
    gk_gamma_z_bc_bootstrap = gk_gamma_z_bc_bootstrap,
    p_gamma_bc = p_gamma_bc,
    # Indirect equating results
    indirect_eq = indirect_eq,
    source_scores = source_scores,
    observed_expected = observed_expected,
    observed_rounded = observed_rounded,
    boot_expected = boot_expected,
    boot_rounded = boot_rounded,
    see = see,
    avg_see = avg_see,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    error_freq = error_freq,
    prop_failed = prop_failed,
    source_min = source_min,
    source_max = source_max,
    fit_ab = fit_ab,
    fit_bc = fit_bc
  )
  
  class(result) <- "leunbach_indirect_bootstrap"
  return(result)
}


#' Run single indirect bootstrap iteration
#' @keywords internal
run_single_indirect_bootstrap <- function(seed, data_list) {
  
  set.seed(seed)
  
  fit_ab_data <- data_list$fit_ab_data
  fit_bc_data <- data_list$fit_bc_data
  direction_ab <- data_list$direction_ab
  direction_bc <- data_list$direction_bc
  method <- data_list$method
  n_scores <- data_list$n_scores
  
  # Initialize output
  result <- list(
    lr_ab = NA,
    lr_bc = NA,
    gk_gamma_z_ab = NA,
    gk_gamma_z_bc = NA,
    expected = rep(NA, n_scores),
    rounded = rep(NA, n_scores),
    failed = rep(1, n_scores)
  )
  
  # Bootstrap fit_ab
  boot_table_ab <- parametric_eq_bootstrap(
    fit_ab_data$test1_scores, fit_ab_data$test2_scores, fit_ab_data$total_scores,
    fit_ab_data$total_score_freq, fit_ab_data$gamma, fit_ab_data$delta, 
    fit_ab_data$sigma, fit_ab_data$xmin, fit_ab_data$xmax, 
    fit_ab_data$ymin, fit_ab_data$ymax
  )
  
  boot_data_ab <- table_to_data(boot_table_ab, fit_ab_data$test1_scores, 
                                fit_ab_data$test2_scores)
  
  if (nrow(boot_data_ab) == 0) return(result)
  
  boot_fit_ab <- tryCatch({
    leunbach_ipf(boot_data_ab,
                 max_score1 = max(fit_ab_data$test1_scores),
                 max_score2 = max(fit_ab_data$test2_scores),
                 verbose = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(boot_fit_ab)) return(result)
  
  result$lr_ab <- boot_fit_ab$g_sq
  result$gk_gamma_z_ab <- boot_fit_ab$gk_gamma_z
  
  # Bootstrap fit_bc
  boot_table_bc <- parametric_eq_bootstrap(
    fit_bc_data$test1_scores, fit_bc_data$test2_scores, fit_bc_data$total_scores,
    fit_bc_data$total_score_freq, fit_bc_data$gamma, fit_bc_data$delta, 
    fit_bc_data$sigma, fit_bc_data$xmin, fit_bc_data$xmax, 
    fit_bc_data$ymin, fit_bc_data$ymax
  )
  
  boot_data_bc <- table_to_data(boot_table_bc, fit_bc_data$test1_scores, 
                                fit_bc_data$test2_scores)
  
  if (nrow(boot_data_bc) == 0) return(result)
  
  boot_fit_bc <- tryCatch({
    leunbach_ipf(boot_data_bc,
                 max_score1 = max(fit_bc_data$test1_scores),
                 max_score2 = max(fit_bc_data$test2_scores),
                 verbose = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(boot_fit_bc)) return(result)
  
  result$lr_bc <- boot_fit_bc$g_sq
  result$gk_gamma_z_bc <- boot_fit_bc$gk_gamma_z
  
  # Compute indirect equating
  boot_indirect <- tryCatch({
    leunbach_indirect_equate(boot_fit_ab, boot_fit_bc,
                             direction_ab = direction_ab,
                             direction_bc = direction_bc,
                             method = method,
                             verbose = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(boot_indirect)) return(result)
  
  # Indirect equating table has 4 columns: source, theta, expected, rounded
  result$expected <- boot_indirect$equating_table$expected  # Column 3
  result$rounded <- boot_indirect$equating_table$rounded    # Column 4
  result$failed <- as.numeric(is.na(result$expected))
  
  return(result)
}


#' Run indirect bootstrap sequentially
#' @keywords internal
run_indirect_bootstrap_sequential <- function(nsim, boot_seeds, boot_data_list, 
                                              verbose = FALSE) {
  
  n_scores <- boot_data_list$n_scores
  
  # Initialize storage - ADD lr and gamma_z for both fits
  lr_ab_bootstrap <- numeric(nsim)
  lr_bc_bootstrap <- numeric(nsim)
  gk_gamma_z_ab_bootstrap <- numeric(nsim)
  gk_gamma_z_bc_bootstrap <- numeric(nsim)
  boot_expected <- matrix(NA, nrow = nsim, ncol = n_scores)
  boot_rounded <- matrix(NA, nrow = nsim, ncol = n_scores)
  boot_failed <- matrix(1, nrow = nsim, ncol = n_scores)
  
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }
  
  for (sim in 1:nsim) {
    result <- run_single_indirect_bootstrap(boot_seeds[sim], boot_data_list)
    
    lr_ab_bootstrap[sim] <- result$lr_ab
    lr_bc_bootstrap[sim] <- result$lr_bc
    gk_gamma_z_ab_bootstrap[sim] <- result$gk_gamma_z_ab
    gk_gamma_z_bc_bootstrap[sim] <- result$gk_gamma_z_bc
    boot_expected[sim, ] <- result$expected
    boot_rounded[sim, ] <- result$rounded
    boot_failed[sim, ] <- result$failed
    
    if (verbose) {
      setTxtProgressBar(pb, sim)
    }
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  list(
    lr_ab_bootstrap = lr_ab_bootstrap,
    lr_bc_bootstrap = lr_bc_bootstrap,
    gk_gamma_z_ab_bootstrap = gk_gamma_z_ab_bootstrap,
    gk_gamma_z_bc_bootstrap = gk_gamma_z_bc_bootstrap,
    boot_expected = boot_expected,
    boot_rounded = boot_rounded,
    boot_failed = boot_failed
  )
}


#' Run indirect bootstrap in parallel
#' @keywords internal
run_indirect_bootstrap_parallel <- function(nsim, boot_seeds, boot_data_list,
                                            n_cores, verbose = FALSE) {
  
  n_scores <- boot_data_list$n_scores
  
  # Initialize storage - ADD lr and gamma_z for both fits
  lr_ab_bootstrap <- numeric(nsim)
  lr_bc_bootstrap <- numeric(nsim)
  gk_gamma_z_ab_bootstrap <- numeric(nsim)
  gk_gamma_z_bc_bootstrap <- numeric(nsim)
  boot_expected <- matrix(NA, nrow = nsim, ncol = n_scores)
  boot_rounded <- matrix(NA, nrow = nsim, ncol = n_scores)
  boot_failed <- matrix(1, nrow = nsim, ncol = n_scores)
  
  mirai:: daemons(n_cores)
  on.exit(mirai::daemons(0), add = TRUE)
  
  if (verbose) {
    cat(sprintf("Starting %d daemons...\n", n_cores))
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    completed <- 0
  }
  
  tasks <- lapply(1:nsim, function(sim) {
    mirai::mirai(
      {
        run_single_indirect_bootstrap(seed, data_list)
      },
      seed = boot_seeds[sim],
      data_list = boot_data_list,
      run_single_indirect_bootstrap = run_single_indirect_bootstrap,
      parametric_eq_bootstrap = parametric_eq_bootstrap,
      table_to_data = table_to_data,
      leunbach_ipf = leunbach_ipf,
      leunbach_equate = leunbach_equate,
      leunbach_indirect_equate = leunbach_indirect_equate,
      compute_indirect_equating = compute_indirect_equating,
      find_nearest_equated = find_nearest_equated,
      adjust_gamma = adjust_gamma,
      calculate_statistics = calculate_statistics,
      calculate_gamma_test = calculate_gamma_test,
      estimate_person_parameter = estimate_person_parameter,
      estimate_theta_optimize = estimate_theta_optimize,
      estimate_theta_newton = estimate_theta_newton,
      estimate_theta_bisection = estimate_theta_bisection,
      expected_score_and_deriv = expected_score_and_deriv,
      calculate_expected_score = calculate_expected_score,
      calculate_true_score = calculate_true_score
    )
  })
  
  for (sim in 1:nsim) {
    result <- mirai::call_mirai(tasks[[sim]])$data
    
    if (! inherits(result, "errorValue")) {
      lr_ab_bootstrap[sim] <- result$lr_ab
      lr_bc_bootstrap[sim] <- result$lr_bc
      gk_gamma_z_ab_bootstrap[sim] <- result$gk_gamma_z_ab
      gk_gamma_z_bc_bootstrap[sim] <- result$gk_gamma_z_bc
      boot_expected[sim, ] <- result$expected
      boot_rounded[sim, ] <- result$rounded
      boot_failed[sim, ] <- result$failed
    } else {
      lr_ab_bootstrap[sim] <- NA
      lr_bc_bootstrap[sim] <- NA
      gk_gamma_z_ab_bootstrap[sim] <- NA
      gk_gamma_z_bc_bootstrap[sim] <- NA
    }
    
    if (verbose) {
      completed <- completed + 1
      setTxtProgressBar(pb, completed)
    }
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  list(
    lr_ab_bootstrap = lr_ab_bootstrap,
    lr_bc_bootstrap = lr_bc_bootstrap,
    gk_gamma_z_ab_bootstrap = gk_gamma_z_ab_bootstrap,
    gk_gamma_z_bc_bootstrap = gk_gamma_z_bc_bootstrap,
    boot_expected = boot_expected,
    boot_rounded = boot_rounded,
    boot_failed = boot_failed
  )
}



#' Print method for leunbach_indirect_bootstrap objects
#' @export
print.leunbach_indirect_bootstrap <- function(x, ...) {
  cat("Leunbach Indirect Equating - Parametric Bootstrap Results\n")
  cat("==========================================================\n\n")
  
  cat(sprintf("Path: %s -> %s -> %s\n", "Test A", "Test B", "Test C"))
  cat(sprintf("Bootstrap samples: %d (%d valid)\n", x$nsim, x$n_valid))
  if (x$parallel) {
    cat(sprintf("Processing:  parallel (%d cores)\n", x$n_cores))
  } else {
    cat("Processing: sequential\n")
  }
  cat(sprintf("Optimization method: %s\n", x$method))
  cat(sprintf("SEE type: %s scores\n\n", x$see_type))
  
  cat("Assessment of significance by parametric bootstrapping:\n\n")
  
  # Equating A-B
  cat(sprintf("Equating A-B (%s -> %s):\n", 
              "Test A", "Test B"))
  cat("  1. Likelihood Ratio Test:\n")
  cat(sprintf("     Observed LR = %.2f (df = %d)\n", x$lr_ab_observed, x$df_ab))
  cat(sprintf("     Asymptotic p-value:    p = %.4f\n", x$fit_ab$p_value))
  cat(sprintf("     Bootstrap p-value:    p = %.4f\n", x$p_lr_ab))
  cat("  2. Goodman-Kruskal Gamma Test (one-sided):\n")
  if (! is.na(x$gk_gamma_z_ab_observed)) {
    cat(sprintf("     Observed Z = %.2f\n", x$gk_gamma_z_ab_observed))
    cat(sprintf("     Asymptotic p-value:   p = %.4f\n", x$fit_ab$gk_gamma_p))
    cat(sprintf("     Bootstrap p-value:    p = %.4f\n\n", x$p_gamma_ab))
  } else {
    cat("     Could not be calculated\n\n")
  }
  
  # Equating B-C
  cat(sprintf("Equating B-C (%s -> %s):\n", 
              "Test B", "Test C"))
  cat("  1. Likelihood Ratio Test:\n")
  cat(sprintf("     Observed LR = %.2f (df = %d)\n", x$lr_bc_observed, x$df_bc))
  cat(sprintf("     Asymptotic p-value:    p = %.4f\n", x$fit_bc$p_value))
  cat(sprintf("     Bootstrap p-value:    p = %.4f\n", x$p_lr_bc))
  cat("  2. Goodman-Kruskal Gamma Test (one-sided):\n")
  if (!is.na(x$gk_gamma_z_bc_observed)) {
    cat(sprintf("     Observed Z = %.2f\n", x$gk_gamma_z_bc_observed))
    cat(sprintf("     Asymptotic p-value:   p = %.4f\n", x$fit_bc$gk_gamma_p))
    cat(sprintf("     Bootstrap p-value:    p = %.4f\n\n", x$p_gamma_bc))
  } else {
    cat("     Could not be calculated\n\n")
  }
  
  print_indirect_see_table(x)
  
  invisible(x)
}


#' Summary method for leunbach_indirect_bootstrap objects
#' @export
summary.leunbach_indirect_bootstrap <- function(object, ...) {
  cat("Leunbach Indirect Equating - Bootstrap Summary\n")
  cat("===============================================\n\n")
  
  cat(sprintf("Path: %s -> %s -> %s\n", "Test A", "Test B", "Test C"))
  cat(sprintf("Bootstrap samples:  %d (%d valid)\n", object$nsim, object$n_valid))
  cat(sprintf("Confidence level: %d%%\n", round(object$conf_level * 100)))
  cat(sprintf("SEE type: %s scores\n\n", object$see_type))
  
  cat("Model Fit Summary:\n\n")
  
  # A-B fit
  cat(sprintf("Equating A-B (%s -> %s):\n", 
              "Test A", "Test B"))
  cat(sprintf("  LR test:     asymptotic p = %.4f, bootstrap p = %.4f\n", 
              object$fit_ab$p_value, object$p_lr_ab))
  if (!is.na(object$gk_gamma_z_ab_observed)) {
    cat(sprintf("  Gamma test:  asymptotic p = %.4f, bootstrap p = %.4f\n\n", 
                object$fit_ab$gk_gamma_p, object$p_gamma_ab))
  } else {
    cat("  Gamma test: could not be calculated\n\n")
  }
  
  # B-C fit
  cat(sprintf("Equating B-C (%s -> %s):\n", 
              "Test B", "Test C"))
  cat(sprintf("  LR test:     asymptotic p = %.4f, bootstrap p = %.4f\n", 
              object$fit_bc$p_value, object$p_lr_bc))
  if (!is.na(object$gk_gamma_z_bc_observed)) {
    cat(sprintf("  Gamma test: asymptotic p = %.4f, bootstrap p = %.4f\n\n", 
                object$fit_bc$gk_gamma_p, object$p_gamma_bc))
  } else {
    cat("  Gamma test: could not be calculated\n\n")
  }
  
  cat(sprintf("Average SEE: %.2f\n", object$avg_see))
  
  # Report failure rates for extreme scores
  valid <- object$source_scores >= object$source_min & 
    object$source_scores <= object$source_max
  failed_rates <- object$prop_failed[valid]
  high_fail <- which(failed_rates > 5)
  
  if (length(high_fail) > 0) {
    cat("\nScores with >5% bootstrap failures:\n")
    scores_with_failures <- object$source_scores[valid][high_fail]
    for (j in seq_along(high_fail)) {
      cat(sprintf("  Score %d: %.1f%% failed\n", 
                  scores_with_failures[j], failed_rates[high_fail[j]]))
    }
  }
  
  invisible(object)
}


#' Print SEE table for indirect equating
#' @keywords internal
print_indirect_see_table <- function(x) {
  
  conf_pct <- round(x$conf_level * 100)
  
  cat(sprintf("Indirect Equating:  %s → %s (with %d%% CI)\n",
              "Test A", "Test C", conf_pct))
  cat("==========================================================\n\n")
  
  valid <- x$source_scores >= x$source_min & x$source_scores <= x$source_max
  
  cat("                                                         Frequency of bootstrap errors\n")
  cat(sprintf("Score  Rounded  Expected    %d%% CI          SEE      -2    -1     0    +1    +2   Failed%%\n", conf_pct))
  cat("---------------------------------------------------------------------------------------------\n")
  
  for (i in which(valid)) {
    if (is.na(x$observed_expected[i])) next
    
    ci_str <- sprintf("[%5.2f, %5.2f]", x$ci_lower[i], x$ci_upper[i])
    cat(sprintf("%5d  %7d    %5.2f   %15s  %5.2f   %5.1f %5.1f %5.1f %5.1f %5.1f   %5.1f%%\n",
                x$source_scores[i], x$observed_rounded[i], x$observed_expected[i],
                ci_str, x$see[i],
                x$error_freq[i, 1], x$error_freq[i, 2], x$error_freq[i, 3],
                x$error_freq[i, 4], x$error_freq[i, 5],
                x$prop_failed[i]))
  }
  
  cat("---------------------------------------------------------------------------------------------\n")
  cat(sprintf("Average SEE:   %.2f\n", x$avg_see))
}


#' Plot method for leunbach_indirect_bootstrap objects
#' @export
plot.leunbach_indirect_bootstrap <- function(x, type = c("equating", "see"), ...) {
  type <- match.arg(type)
  
  valid <- x$source_scores >= x$source_min & x$source_scores <= x$source_max &
    ! is.na(x$observed_expected)
  
  if (type == "equating") {
    plot(x$source_scores[valid], x$observed_expected[valid], type = "b", pch = 19,
         ylim = range(c(x$ci_lower[valid], x$ci_upper[valid]), na.rm = TRUE),
         xlab = x$indirect_eq$source_name,
         ylab = paste("Expected", "Test C"),
         main = sprintf("Indirect Equating: %s → %s (%d%% CI)",
                        "Test A", "Test C",
                        round(x$conf_level * 100)))
    
    polygon(c(x$source_scores[valid], rev(x$source_scores[valid])),
            c(x$ci_lower[valid], rev(x$ci_upper[valid])),
            col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
    lines(x$source_scores[valid], x$observed_expected[valid], type = "b", pch = 19)
    grid()
    
  } else if (type == "see") {
    plot(x$source_scores[valid], x$see[valid], type = "b", pch = 19,
         xlab = x$indirect_eq$source_name,
         ylab = "SEE",
         main = sprintf("Standard Error of Equating: %s → %s",
                        "Test A", "Test C"))
    grid()
    abline(h = x$avg_see, col = "red", lty = 2)
  }
  
  invisible(x)
}


#' Get indirect equating table
#'
#' @param boot A leunbach_indirect_bootstrap object
#' @return A data frame with equating results and CIs
#' @export
get_indirect_equating_table <- function(boot) {
  
  if (!inherits(boot, "leunbach_indirect_bootstrap")) {
    stop("Input must be a leunbach_indirect_bootstrap object")
  }
  
  valid <- boot$source_scores >= boot$source_min & 
    boot$source_scores <= boot$source_max
  
  # Get theta values from indirect equating table and convert to log
  theta_values <- boot$indirect_eq$equating_table$theta[valid]
  log_theta <- ifelse(is.na(theta_values) | theta_values <= 0, NA, round(log(theta_values), 4))
  
  result <- data.frame(
    source = boot$source_scores[valid],
    log_theta = log_theta,
    rounded = boot$observed_rounded[valid],
    expected = round(boot$observed_expected[valid], 2),
    ci_lower = round(boot$ci_lower[valid], 2),
    ci_upper = round(boot$ci_upper[valid], 2),
    see = round(boot$see[valid], 2),
    pct_failed = round(boot$prop_failed[valid], 1)
  )
  
  colnames(result)[1] <- boot$indirect_eq$source_name
  
  return(result)
}