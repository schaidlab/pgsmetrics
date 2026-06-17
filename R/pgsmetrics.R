#' Compare Polygenic Scores (PGS) using bootstrapping
#'
#' This function compares Polygenic Scores (PGS) by bootstrapping the difference
#' between them using various metrics. It supports single PGS evaluation, pairwise
#' comparisons, and ranking of multiple PGS.
#'
#' @param data A data.table containing necessary data: PGS, dependent variable, and covariates.
#' @param K Double, between 0 and 1. Prevalence of the binary trait in the population. Used by $R^2$ on the liability scale.
#' @param pgs Character vector. Metrics are computed for these pgs. Must exist in `data`.
#' @param dep Character. Name of the dependent variable in `data`. Binary variable must be 0/1 (otherwise continuous assumed).
#' @param covars Character vector. Names of covariate to adjust for. Must be in columns in `data`.
#' @param boot Boolean. Whether to perform bootstrap or not.
#' @param boot_method Character vector. Either `standard` or `blb`. `standard` is the usual non-parametric bootstrap, `blb` is the Bag of Little Bootstraps, useful for analysis on large datasets (f.x., $n>500000$) with many PGS. However, BLB and the standard bootstrap may not yield identical results. The standard bootstrap uses the observed statistics (all observations) and gives percentile CIs with resampling. The BLB bootstrap estimates the mean across all resamples (instead of the observed value) and percentile CIs.
#' @param blb_s Integer. TODO. Default: 20.
#' @param blb_r Integer. TODO. Default: 100.
#' @param blb_b Double between 0 and 1. TODO. Default: 0.8.
#' @param n_boot Integer. Number of bootstrap iterations. Default is 1000.
#' @param n_cores Integer. Number of cores for parallel processing. Default is 1.
#' @param missing Warns if missing values in `data`. Use missing="drop" to drop rows with missing values. However, we recommend imputing missing values instead. Default: "warn".
#' @param return_boot_stats Return bootstrap output before processing (for debugging). Default: FALSE.
#' @param custom_metrics_list List of custom metrics to be summarized. 
#'
#' @return An object of class "pgsmetrics" containing:
#'   \item{metrics}{A data.table with calculated metrics for each PGS. In this data.table, `type` can be partial/full/cov, where partial=full-cov, full=model with covariates+pgs, cov=model with only covariates. When using `boot_method="blb"`, the `observed` column is the median value of the metric statistic, but when using `boot_method="standard"`, it is the value from an analysis on the unmodified dataset.}
#'   \item{boot_output}{The raw output from the bootstrapping process.}
#'   \item{call}{The original function call.}
#'
#' @examples
#' # Simulate data
#' data <- simulate_data(n = 1000, n_pgs = 3)
#'
#' # Single PGS evaluation
#' result_single <- pgsmetrics(data,
#'     pgs = c("pgs1", "pgs2", "pgs3"),
#'     dep = "y_bin",
#'     covars = c("age", "sex"),
#'     n_boot = 100
#' )
#' # print(result_single)
#'
#' @export
pgsmetrics <- function(
    data,
    pgs,
    dep,
    covars = NULL,
    K = NULL,
    boot = TRUE,
    boot_method = "standard",
    blb_s = 10, blb_r = 50, blb_b = 0.9,
    n_boot = 1e3,
    n_cores = 1,
    missing = "warn",
    return_boot_stats = FALSE,
    custom_metrics_list = NULL) {

    validate_inputs(
        data,
        K,
        pgs,
        dep,
        covars,
        boot_method,
        blb_s,
        blb_r,
        blb_b,
        n_boot,
        n_cores,
        boot,
        missing
    )

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }

    data <- handle_missing_values(data, pgs, dep, covars, missing)

    quant <- determine_outcome_type(data[[dep]])

    metrics_setup <- setup_metrics(quant, custom_metrics_list)
    metrics_list <- metrics_setup$metrics_list
    n_metrics <- metrics_setup$n_metrics
    ns_metrics <- metrics_setup$ns_metrics

    data_m <- prepare_data_matrix(data, covars, pgs)

    if (boot) {
        message(paste0("bootstrapping (", boot_method, ")..."))
        bootstrap_results <- run_bootstrap(
            boot_method,
            n_boot,
            data,
            data_m,
            quant,
            K,
            pgs,
            covars,
            data[[dep]],
            dep,
            metrics_list,
            n_cores,
            blb_b,
            blb_s,
            blb_r
        )
        # -> return_bootstrap_results
        if (return_boot_stats) {
            return(bootstrap_results)
        }
        message("processing resamples...")
        bootstrap_results <- rbindlist(bootstrap_results)
        processed_results <- process_metrics(bootstrap_results, boot_method)
        metric_differences <- process_diff(bootstrap_results, pgs, boot_method)
        metrics <- processed_results$metrics
        ranks <- processed_results$ranks
    }

    if (boot_method == "standard") {
        message("computing observed values...")
        # FIXME y stuff
        observed <- compute_observed(data_m, quant, K, pgs, covars, data[[dep]], metrics_list)
        message("computing differences...")
        observed_differences <- compute_observed_differences(observed$observed, pgs)
        # Combine results
        if (boot) {
            metrics <- observed$metrics[metrics, on = c("pgs", "metric", "type")]
            ranks <- observed$ranks[ranks, on = c("pgs", "metric", "type")]
            metric_differences <- metric_differences[
                observed_differences[, .(pgs1_pgs2, metric, observed)],
                on = c("pgs1_pgs2", "metric")
            ]
            setcolorder(metric_differences, "observed", after = "metric")
        } else {
            metric_differences <- observed_differences
            metric_differences[, c("ci_lower", "ci_upper") := NA_real_]
            metrics <- observed$metrics
            metrics[, c("ci_lower", "ci_upper") := NA_real_]
            ranks <- observed$ranks
            ranks[, c("ci_lower", "ci_upper") := NA_real_]
            bootstrap_results <- NULL
        }
    }

    # Finalize result structure
    metrics[, pgs := factor(pgs, levels = unique(pgs))]
    ranks[, pgs := factor(pgs, levels = unique(pgs))]
    metric_differences[, `:=`(
        pgs1 = factor(pgs1, levels = pgs),
        pgs2 = factor(pgs2, levels = pgs),
        pgs1_pgs2 = factor(pgs1_pgs2, levels = unique(pgs1_pgs2))
    )]
    setkey(metric_differences, pgs1, pgs2, metric)
    setcolorder(metrics, "observed", before = "ci_lower")
    setcolorder(ranks, "observed", before = "ci_lower")

    structure(
        list(
            metrics = metrics,
            ranks = ranks,
            diff = metric_differences,
            boot_output = bootstrap_results,
            call = call,
            sample_size = nrow(data),
            dependent_variable = dep,
            covariates = covars,
            bootstrap_method = boot_method,
            blb_b = if (boot_method == "blb") blb_b else NA,
            blb_s = if (boot_method == "blb") blb_s else NA,
            blb_r = if (boot_method == "blb") blb_r else NA,
            metrics_list = metrics_list
        ),
        class = "pgsmetrics"
    )
}


handle_missing_values <- function(data, pgs, dep, covars, missing) {
    if (missing == "drop") {
        n_before <- nrow(data)
        na_omit_cols <- c(pgs, dep, covars)
        data <- na.omit(data, cols = na_omit_cols)
        n_after <- nrow(data)

        if (n_before > n_after) {
            message(sprintf(
                "Removed %d rows with missing values. Using remaining %d rows.",
                n_before - n_after, n_after
            ))
        }
    }
    data
}


determine_outcome_type <- function(dep_var) {
    dep_unique <- unique(dep_var)
    dep_unique_length <- length(dep_unique)
    if (dep_unique_length == 1) {
        stop("Outcome variable must have >1 values.")
    }
    if (length(dep_unique) == 2) {
        if (!all(sort(dep_unique) == c(0, 1))) {
            stop("Binary outcome variable must be coded as 0/1")
        }
        message("Binary outcome variable detected")
        return(FALSE)
    }
    message("Continuous outcome variable detected")
    return(TRUE)
}


setup_metrics <- function(quant, custom_metrics_list) {
    metrics_type <- if (quant) "continuous" else "binary"
    metrics_list <- get_default_metrics() |>
        filter_metrics(metrics_type)

    if (!is.null(custom_metrics_list)) {
        message("Adding custom metrics")
        for (metric in custom_metrics_list) {
            if (!is.null(metrics_list[[metric$name]])) {
                message(sprintf("Overriding default metric: %s", metric$name))
            }
            metrics_list[[metric$name]] <- metric
        }
    }

    list(
        metrics_list = metrics_list,
        n_metrics = length(metrics_list),
        ns_metrics = names(metrics_list)
    )
}


prepare_data_matrix <- function(data, covars, pgs) {
    data_cols <- c(covars, pgs)
    as.matrix(cbind(intercept = 1, data[, ..data_cols]))
}


run_bootstrap <- function(
    boot_method,
    n_boot,
    data,
    data_m,
    quant,
    K,
    pgs,
    covars,
    y,
    dep,
    metrics_list,
    n_cores,
    blb_b,
    blb_s,
    blb_r
) {
    if (boot_method == "standard") {
        boot_stats <- pbapply::pblapply(1:n_boot, function(i) {
            boot_standard(
                data_m = data_m,
                ind = sample(nrow(data_m), replace = TRUE),
                boot_i = i,
                quant = quant,
                K = K,
                pgs = pgs,
                covars = covars,
                y = y,
                w = NULL,
                metrics_list = metrics_list
            )
        }, cl = n_cores)
    }
    if (boot_method == "blb") {
        b <- floor(nrow(data_m)^blb_b)
        idx <- CJ(s = 1:blb_s, r = 1:blb_r)
        inds <- list()
        for (s in 1:blb_s) {
            inds[[s]] <- sample(nrow(data_m), b)
        }
        boot_stats <- pbapply::pblapply(
            1:nrow(idx),
            function(j) {
                boot_blb(
                    j = idx$s[j],
                    ind = inds[[idx$s[j]]],
                    # FIXME remove
                    data = data,
                    n = nrow(data_m),
                    b = b,
                    r = idx$r[j],
                    K = K,
                    pgs = pgs,
                    quant = quant,
                    dep = dep,
                    covars = covars
                )
            },
            cl = n_cores
        )
    }
    boot_stats
}


compute_observed <- function(data_m, quant, K, pgs, covars, y, metrics_list) {
    observed <- boot_standard(
        data_m,
        ind = 1:nrow(data_m),
        boot_i = NA,
        quant,
        K,
        pgs,
        covars,
        y,
        w = NULL,
        metrics_list
    )
    observed_metrics <- melt(observed[, boot_i := NULL], id.vars = c("pgs", "metric"), variable.name = "type", value.name = "observed")
    observed_metrics$observed <- round(observed_metrics$observed, 5)
    setkey(observed_metrics, pgs, metric, type)

    observed_ranks <- observed_metrics[, .(
        pgs = pgs,
        observed = frank(observed, ties.method = "average")
    ),
    by = c("metric", "type")
    ]
    setcolorder(observed_ranks, c("pgs", "metric", "type"))

    list(
        observed = observed,
        metrics = observed_metrics,
        ranks = observed_ranks
    )
}


compute_observed_differences <- function(observed_values, pgs) {
    temp <- copy(observed_values)
    temp$boot_i <- 1
    observed_differences <- process_diff(temp, pgs, "standard")
    observed_differences[, ci_lower := NULL]
    setnames(observed_differences, "ci_upper", "observed")
    observed_differences
}
