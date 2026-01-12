#' Parametric Bootstrap for Leunbach Model with Standard Error of Equating
#'
#' @description
#' Performs parametric bootstrapping to assess significance of tests and
#' compute standard errors of equating (SEE) with confidence intervals.
#' Supports parallel processing using the mirai package.
#'
#' @param fit A leunbach_ipf object from leunbach_ipf()
#' @param nsim Number of bootstrap samples (default: 1000)
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param see_type Type of SEE calculation: "rounded" uses rounded (integer)
#'        scores, "expected" uses continuous expected scores
#' @param method Optimization method for person parameter estimation: 
#'        "optimize" (default) uses stats::optimize() with Brent's method,
#'        "newton" uses custom Newton-Raphson with bisection fallback
#' @param parallel Use parallel processing if mirai package is available (default: TRUE)
#' @param n_cores Number of cores to use for parallel processing. 
#'        Default NULL uses all available cores minus one.
#' @param verbose Print progress messages
#' @param seed Random seed for reproducibility (optional)
#'
#' @return A list of class "leunbach_bootstrap" containing bootstrap results
#'   and standard errors of equating
#'
#' @export
leunbach_bootstrap <- function(fit, nsim = 1000, conf_level = 0.95,
                               see_type = c("rounded", "expected"),
                               method = c("optimize", "newton"),
                               parallel = TRUE, n_cores = NULL,
                               verbose = FALSE, seed = NULL) {
  
  if (!inherits(fit, "leunbach_ipf")) {
    stop("Input must be a leunbach_ipf object")
  }
  
  see_type <- match.arg(see_type)
  method <- match.arg(method)
  
  # Check for mirai availability
  use_parallel <- parallel && requireNamespace("mirai", quietly = TRUE)
  
  if (parallel && ! use_parallel) {
    message("Install 'mirai' package for parallel processing:  install.packages('mirai')")
    message("Running sequentially...")
  }
  
  # Set up cores for parallel processing
  if (use_parallel) {
    if (is.null(n_cores)) {
      stop(paste0("For parallel processing, you need to specify how many cores to use,
                  by setting the option `n_cores` to a useful number.
                  It seems like your computer has ",parallel::detectCores()," cores available."))
    }
    n_cores <- min(n_cores, nsim)
  }
  
  # Set seed for reproducibility
  if (! is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate seeds for each bootstrap sample
  boot_seeds <- sample.int(.Machine$integer.max, nsim)
  
  # Extract info from fitted model
  gamma <- fit$gamma
  delta <- fit$delta
  sigma <- fit$sigma
  total_score_freq <- fit$total_score_freq
  
  test1_scores <- fit$test1_scores
  test2_scores <- fit$test2_scores
  total_scores <- fit$total_scores
  
  xmin <- fit$xmin
  xmax <- fit$xmax
  ymin <- fit$ymin
  ymax <- fit$ymax
  
  n_test1 <- length(test1_scores)
  n_test2 <- length(test2_scores)
  
  lr_observed <- fit$g_sq
  
  # Get observed equating (point estimates)
  eq_1to2 <- leunbach_equate(fit, direction = "1to2", method = method)
  eq_2to1 <- leunbach_equate(fit, direction = "2to1", method = method)
  
  observed_eq_1to2 <- eq_1to2$equating_table[, 2]
  observed_eq_2to1 <- eq_2to1$equating_table[, 2]
  observed_rd_1to2 <- eq_1to2$equating_table[, 3]
  observed_rd_2to1 <- eq_2to1$equating_table[, 3]
  
  if (verbose) {
    cat("Parametric Bootstrap for Leunbach Model\n")
    cat("========================================\n\n")
    cat(sprintf("Optimization method: %s\n", method))
    if (use_parallel) {
      cat(sprintf("Running %d bootstrap samples using %d cores...\n\n", nsim, n_cores))
    } else {
      cat(sprintf("Running %d bootstrap samples sequentially...\n\n", nsim))
    }
  }
  
  # Prepare data for bootstrap function
  boot_data_list <- list(
    test1_scores = test1_scores,
    test2_scores = test2_scores,
    total_scores = total_scores,
    total_score_freq = total_score_freq,
    gamma = gamma,
    delta = delta,
    sigma = sigma,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    n_test1 = n_test1,
    n_test2 = n_test2,
    method = method
  )
  
  # Run bootstrap (parallel or sequential)
  if (use_parallel) {
    boot_results <- run_bootstrap_parallel(
      nsim = nsim,
      boot_seeds = boot_seeds,
      boot_data_list = boot_data_list,
      n_cores = n_cores,
      verbose = verbose
    )
  } else {
    boot_results <- run_bootstrap_sequential(
      nsim = nsim,
      boot_seeds = boot_seeds,
      boot_data_list = boot_data_list,
      verbose = verbose
    )
  }
  
  # Extract results
  lr_bootstrap <- boot_results$lr_bootstrap
  gk_gamma_z_bootstrap <- boot_results$gk_gamma_z_bootstrap
  boot_eq_1to2 <- boot_results$boot_eq_1to2
  boot_eq_2to1 <- boot_results$boot_eq_2to1
  boot_rd_1to2 <- boot_results$boot_rd_1to2
  boot_rd_2to1 <- boot_results$boot_rd_2to1
  failed_1to2 <- boot_results$failed_1to2
  failed_2to1 <- boot_results$failed_2to1
  
  # Count significant LR results
  valid_lr <- !is.na(lr_bootstrap)
  n_significant_lr <- sum(lr_bootstrap[valid_lr] >= lr_observed)
  p_lr <- n_significant_lr / sum(valid_lr)
  
  # Get observed Gamma Z statistic and compute bootstrap p-value
  gk_gamma_z_observed <- fit$gk_gamma_z
  
  # Compute bootstrap p-value for Gamma test
  # (proportion of bootstrap |Z| >= observed |Z|)
  valid_gamma <- !is.na(gk_gamma_z_bootstrap)
  if (sum(valid_gamma) > 0 && !is.na(gk_gamma_z_observed)) {
    n_significant_gamma <- sum(abs(gk_gamma_z_bootstrap[valid_gamma]) >= abs(gk_gamma_z_observed))
    p_gamma <- n_significant_gamma / sum(valid_gamma)
  } else {
    n_significant_gamma <- NA
    p_gamma <- NA
  }
  
  # Confidence interval quantiles
  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)
  
  ci_lower_1to2 <- apply(boot_eq_1to2, 2, quantile, probs = probs[1], na.rm = TRUE)
  ci_upper_1to2 <- apply(boot_eq_1to2, 2, quantile, probs = probs[2], na.rm = TRUE)
  ci_lower_2to1 <- apply(boot_eq_2to1, 2, quantile, probs = probs[1], na.rm = TRUE)
  ci_upper_2to1 <- apply(boot_eq_2to1, 2, quantile, probs = probs[2], na.rm = TRUE)
  
  # Standard Error of Equating
  if (see_type == "rounded") {
    see_1to2 <- calculate_bootstrap_see(boot_rd_1to2)
    see_2to1 <- calculate_bootstrap_see(boot_rd_2to1)
  } else {
    see_1to2 <- calculate_bootstrap_see(boot_eq_1to2)
    see_2to1 <- calculate_bootstrap_see(boot_eq_2to1)
  }
  
  # Average SEE
  valid_1to2 <- test1_scores >= xmin & test1_scores <= xmax
  valid_2to1 <- test2_scores >= ymin & test2_scores <= ymax
  
  avg_see_1to2 <- mean(see_1to2[valid_1to2], na.rm = TRUE)
  avg_see_2to1 <- mean(see_2to1[valid_2to1], na.rm = TRUE)
  
  # Frequency of bootstrap errors
  error_freq_1to2 <- compute_error_frequencies(boot_rd_1to2, observed_rd_1to2, test1_scores)
  error_freq_2to1 <- compute_error_frequencies(boot_rd_2to1, observed_rd_2to1, test2_scores)
  
  # Proportion of failed cases
  prop_failed_1to2 <- colMeans(failed_1to2, na.rm = TRUE) * 100
  prop_failed_2to1 <- colMeans(failed_2to1, na.rm = TRUE) * 100
  
  if (verbose) {
    cat(sprintf("\nBootstrap complete.\n"))
    cat(sprintf("  Valid samples: %d of %d\n", sum(valid_lr), nsim))
    cat(sprintf("  Bootstrap p-value for LR test: %.3f\n", p_lr))
    cat(sprintf("  Average SEE (Test1 to Test2): %.2f\n", avg_see_1to2))
    cat(sprintf("  Average SEE (Test2 to Test1): %.2f\n", avg_see_2to1))
  }
  
  result <- list(
    nsim = nsim,
    n_valid = sum(valid_lr),
    conf_level = conf_level,
    see_type = see_type,
    method = method,
    parallel = use_parallel,
    n_cores = if (use_parallel) n_cores else 1L,
    lr_observed = lr_observed,
    df = fit$df,
    lr_bootstrap = lr_bootstrap,
    p_lr = p_lr,
    n_significant_lr = n_significant_lr,
    gk_gamma_z_observed = gk_gamma_z_observed,
    gk_gamma_z_bootstrap = gk_gamma_z_bootstrap,
    p_gamma = p_gamma,
    n_significant_gamma = n_significant_gamma,
    eq_1to2 = eq_1to2,
    eq_2to1 = eq_2to1,
    boot_eq_1to2 = boot_eq_1to2,
    boot_eq_2to1 = boot_eq_2to1,
    boot_rd_1to2 = boot_rd_1to2,
    boot_rd_2to1 = boot_rd_2to1,
    see_1to2 = see_1to2,
    see_2to1 = see_2to1,
    avg_see_1to2 = avg_see_1to2,
    avg_see_2to1 = avg_see_2to1,
    ci_lower_1to2 = ci_lower_1to2,
    ci_upper_1to2 = ci_upper_1to2,
    ci_lower_2to1 = ci_lower_2to1,
    ci_upper_2to1 = ci_upper_2to1,
    error_freq_1to2 = error_freq_1to2,
    error_freq_2to1 = error_freq_2to1,
    prop_failed_1to2 = prop_failed_1to2,
    prop_failed_2to1 = prop_failed_2to1,
    fit = fit
  )
  
  class(result) <- "leunbach_bootstrap"
  return(result)
}


#' Run a single bootstrap iteration
#' @keywords internal
run_single_bootstrap <- function(seed, data_list) {
  
  set.seed(seed)
  
  test1_scores <- data_list$test1_scores
  test2_scores <- data_list$test2_scores
  total_scores <- data_list$total_scores
  total_score_freq <- data_list$total_score_freq
  gamma <- data_list$gamma
  delta <- data_list$delta
  sigma <- data_list$sigma
  xmin <- data_list$xmin
  xmax <- data_list$xmax
  ymin <- data_list$ymin
  ymax <- data_list$ymax
  n_test1 <- data_list$n_test1
  n_test2 <- data_list$n_test2
  method <- data_list$method
  
  # Initialize output
  result <- list(
    lr = NA,
    gk_gamma_z = NA,
    eq_1to2 = rep(NA, n_test1),
    eq_2to1 = rep(NA, n_test2),
    rd_1to2 = rep(NA, n_test1),
    rd_2to1 = rep(NA, n_test2),
    failed_1to2 = rep(1, n_test1),
    failed_2to1 = rep(1, n_test2)
  )
  
  # Generate bootstrap table
  boot_table <- parametric_eq_bootstrap(
    test1_scores, test2_scores, total_scores,
    total_score_freq, gamma, delta, sigma,
    xmin, xmax, ymin, ymax
  )
  
  # Create bootstrap data from table
  boot_data <- table_to_data(boot_table, test1_scores, test2_scores)
  
  if (nrow(boot_data) == 0) {
    return(result)
  }
  
  # Re-estimate model on bootstrap sample
  boot_fit <- tryCatch({
    leunbach_ipf(boot_data,
                 max_score1 = max(test1_scores),
                 max_score2 = max(test2_scores),
                 verbose = FALSE)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(boot_fit)) {
    return(result)
  }
  
  result$lr <- boot_fit$g_sq
  result$gk_gamma_z <- boot_fit$gk_gamma_z
  
  # Compute equating for this bootstrap sample
  boot_eq_1to2_result <- tryCatch({
    leunbach_equate(boot_fit, direction = "1to2", method = method)
  }, error = function(e) NULL)
  
  boot_eq_2to1_result <- tryCatch({
    leunbach_equate(boot_fit, direction = "2to1", method = method)
  }, error = function(e) NULL)
  
  if (!is.null(boot_eq_1to2_result)) {
    result$eq_1to2 <- boot_eq_1to2_result$equating_table[, 2]
    result$rd_1to2 <- boot_eq_1to2_result$equating_table[, 3]
    result$failed_1to2 <- as.numeric(is.na(result$eq_1to2))
  }
  
  if (!is.null(boot_eq_2to1_result)) {
    result$eq_2to1 <- boot_eq_2to1_result$equating_table[, 2]
    result$rd_2to1 <- boot_eq_2to1_result$equating_table[, 3]
    result$failed_2to1 <- as.numeric(is.na(result$eq_2to1))
  }
  
  return(result)
}


#' Run bootstrap samples in parallel using mirai
#' @keywords internal
run_bootstrap_parallel <- function(nsim, boot_seeds, boot_data_list, n_cores, verbose = FALSE) {
  
  n_test1 <- boot_data_list$n_test1
  n_test2 <- boot_data_list$n_test2
  test1_scores <- boot_data_list$test1_scores
  test2_scores <- boot_data_list$test2_scores
  
  # Initialize storage
  lr_bootstrap <- numeric(nsim)
  gk_gamma_z_bootstrap <- numeric(nsim)
  boot_eq_1to2 <- matrix(NA, nrow = nsim, ncol = n_test1)
  boot_eq_2to1 <- matrix(NA, nrow = nsim, ncol = n_test2)
  boot_rd_1to2 <- matrix(NA, nrow = nsim, ncol = n_test1)
  boot_rd_2to1 <- matrix(NA, nrow = nsim, ncol = n_test2)
  failed_1to2 <- matrix(0, nrow = nsim, ncol = n_test1)
  failed_2to1 <- matrix(0, nrow = nsim, ncol = n_test2)
  
  colnames(boot_eq_1to2) <- test1_scores
  colnames(boot_eq_2to1) <- test2_scores
  colnames(boot_rd_1to2) <- test1_scores
  colnames(boot_rd_2to1) <- test2_scores
  colnames(failed_1to2) <- test1_scores
  colnames(failed_2to1) <- test2_scores
  
  # Set up mirai daemons
  mirai:: daemons(n_cores)
  on.exit(mirai:: daemons(0), add = TRUE)
  
  if (verbose) {
    cat(sprintf("Starting %d daemons...\n", n_cores))
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    completed <- 0
  }
  
  # Submit all tasks
  tasks <- lapply(1:nsim, function(sim) {
    mirai::mirai(
      {
        run_single_bootstrap(seed, data_list)
      },
      seed = boot_seeds[sim],
      data_list = boot_data_list,
      run_single_bootstrap = run_single_bootstrap,
      parametric_eq_bootstrap = parametric_eq_bootstrap,
      table_to_data = table_to_data,
      leunbach_ipf = leunbach_ipf,
      leunbach_equate = leunbach_equate,
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
  
  # Collect results as they complete
  for (sim in 1:nsim) {
    result <- mirai::call_mirai(tasks[[sim]])$data
    
    if (! inherits(result, "errorValue")) {
      lr_bootstrap[sim] <- result$lr
      gk_gamma_z_bootstrap[sim] <- result$gk_gamma_z
      boot_eq_1to2[sim, ] <- result$eq_1to2
      boot_eq_2to1[sim, ] <- result$eq_2to1
      boot_rd_1to2[sim, ] <- result$rd_1to2
      boot_rd_2to1[sim, ] <- result$rd_2to1
      failed_1to2[sim, ] <- result$failed_1to2
      failed_2to1[sim, ] <- result$failed_2to1
    } else {
      lr_bootstrap[sim] <- NA
      gk_gamma_z_bootstrap[sim] <- NA
      failed_1to2[sim, ] <- 1
      failed_2to1[sim, ] <- 1
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
    lr_bootstrap = lr_bootstrap,
    gk_gamma_z_bootstrap = gk_gamma_z_bootstrap,
    boot_eq_1to2 = boot_eq_1to2,
    boot_eq_2to1 = boot_eq_2to1,
    boot_rd_1to2 = boot_rd_1to2,
    boot_rd_2to1 = boot_rd_2to1,
    failed_1to2 = failed_1to2,
    failed_2to1 = failed_2to1
  )
}


#' Run bootstrap samples sequentially
#' @keywords internal
run_bootstrap_sequential <- function(nsim, boot_seeds, boot_data_list, verbose = FALSE) {
  
  n_test1 <- boot_data_list$n_test1
  n_test2 <- boot_data_list$n_test2
  test1_scores <- boot_data_list$test1_scores
  test2_scores <- boot_data_list$test2_scores
  
  # Initialize storage
  lr_bootstrap <- numeric(nsim)
  gk_gamma_z_bootstrap <- numeric(nsim)
  boot_eq_1to2 <- matrix(NA, nrow = nsim, ncol = n_test1)
  boot_eq_2to1 <- matrix(NA, nrow = nsim, ncol = n_test2)
  boot_rd_1to2 <- matrix(NA, nrow = nsim, ncol = n_test1)
  boot_rd_2to1 <- matrix(NA, nrow = nsim, ncol = n_test2)
  failed_1to2 <- matrix(0, nrow = nsim, ncol = n_test1)
  failed_2to1 <- matrix(0, nrow = nsim, ncol = n_test2)
  
  colnames(boot_eq_1to2) <- test1_scores
  colnames(boot_eq_2to1) <- test2_scores
  colnames(boot_rd_1to2) <- test1_scores
  colnames(boot_rd_2to1) <- test2_scores
  colnames(failed_1to2) <- test1_scores
  colnames(failed_2to1) <- test2_scores
  
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  }
  
  for (sim in 1:nsim) {
    result <- run_single_bootstrap(boot_seeds[sim], boot_data_list)
    
    lr_bootstrap[sim] <- result$lr
    gk_gamma_z_bootstrap[sim] <- result$gk_gamma_z
    boot_eq_1to2[sim, ] <- result$eq_1to2
    boot_eq_2to1[sim, ] <- result$eq_2to1
    boot_rd_1to2[sim, ] <- result$rd_1to2
    boot_rd_2to1[sim, ] <- result$rd_2to1
    failed_1to2[sim, ] <- result$failed_1to2
    failed_2to1[sim, ] <- result$failed_2to1
    
    if (verbose) {
      setTxtProgressBar(pb, sim)
    }
  }
  
  if (verbose) {
    close(pb)
    cat("\n")
  }
  
  list(
    lr_bootstrap = lr_bootstrap,
    gk_gamma_z_bootstrap = gk_gamma_z_bootstrap,
    boot_eq_1to2 = boot_eq_1to2,
    boot_eq_2to1 = boot_eq_2to1,
    boot_rd_1to2 = boot_rd_1to2,
    boot_rd_2to1 = boot_rd_2to1,
    failed_1to2 = failed_1to2,
    failed_2to1 = failed_2to1
  )
}


#' Print method for leunbach_bootstrap objects
#' @export
print.leunbach_bootstrap <- function(x, ...) {
  cat("Leunbach Model - Parametric Bootstrap Results\n")
  cat("==============================================\n\n")
  
  cat(sprintf("Bootstrap samples: %d (%d valid)\n", x$nsim, x$n_valid))
  if (x$parallel) {
    cat(sprintf("Processing:  parallel (%d cores)\n", x$n_cores))
  } else {
    cat("Processing: sequential\n")
  }
  cat(sprintf("Optimization method: %s\n", x$method))
  cat(sprintf("SEE type: %s scores\n\n", x$see_type))
  
  cat("Assessment of significance by parametric bootstrapping:\n\n")
  
  # 1. Likelihood Ratio Test
  cat("1. Likelihood Ratio Test:\n")
  cat(sprintf("   Observed LR = %.2f (df = %d)\n", x$lr_observed, x$df))
  cat(sprintf("   Asymptotic p-value:  p = %.4f\n", x$fit$p_value))
  cat(sprintf("   Bootstrap p-value:   p = %.4f\n\n", x$p_lr))
  
  # 2. Goodman-Kruskal Gamma Test
  cat("2. Goodman-Kruskal Gamma Test:\n")
  if (!is.na(x$gk_gamma_z_observed)) {
    cat(sprintf("   Observed |Z| = %.2f\n", abs(x$gk_gamma_z_observed)))
    cat(sprintf("   Asymptotic p-value:  p = %.4f\n", x$fit$gk_gamma_p))
    cat(sprintf("   Bootstrap p-value:   p = %.4f\n\n", x$p_gamma))
  } else {
    cat("   Could not be calculated\n\n")
  }
  
  print_see_table(x, direction = "1to2")
  cat("\n")
  print_see_table(x, direction = "2to1")
  
  invisible(x)
}