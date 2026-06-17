#' Metric Class
#'
#' @description
#' An R6 class representing a metric for evaluating Polygenic Scores (PGS).
#'
#' @field name Character. The name of the metric.
#' @field fn Function. The function used to calculate the metric.
#' @field type binary/continuous/both.
#' @field description Longer name.
#'
#' @importFrom R6 R6Class
#' @export
Metric <- R6::R6Class("Metric",
    public = list(
        name = NULL,
        fn = NULL,
        type = NULL,
        description = NULL,

        #' @description
        #' Create a new Metric object.
        #' @param name Character. The name of the metric.
        #' @param fn Function. The function used to calculate the metric.
        #' @param type binary/continuous/both.
        #' @param description Longer name.
        initialize = function(name, fn, type, description) {
            checkmate::assert_string(name, null.ok = FALSE)
            checkmate::assert_function(fn, null.ok = FALSE)
            checkmate::assert_choice(type, c("binary", "continuous", "both"), null.ok = FALSE)
            checkmate::assert_string(description, null.ok = FALSE)
            self$name <- name
            self$fn <- fn
            self$type <- type
            self$description <- description
        },

        #' @description
        #' Calculate the metric value.
        #' @param m glm model
        #' @param m_lm lm model (family='gaussian')
        #' @param p probabilities from glm
        #' @param y response
        #' @param liab_prep data input to r2_liability
        #' @param n number of observations
        #' @param w weights
        #' @return The calculated metric value.
        calculate = function(m = NULL,
                             m_lm = NULL,
                             p = NULL,
                             y = NULL,
                             liab_prep = NULL,
                             n = NULL,
                             w = NULL) {
            self$fn(
                m = m,
                m_lm = m_lm,
                p = p,
                y = y,
                liab_prep = liab_prep,
                n = n,
                w = w
            )
        }
    )
)


###' Metric Registry
# metric_registry <- new.env(parent = emptyenv())


#' Register a new metric
#'
#' @param metrics_list List of metrics
#' @param name Character. Name of the metric.
#' @param fn Function. The metric calculation function.
#' @param type binary/continuous/both.
#' @param description Description.
#' @export
register_metric <- function(metrics_list, name, fn, type, description) {
    metrics_list[[name]] <- Metric$new(name, fn, type, description)
    metrics_list
}


#' Filter a list of metrics by type
#'
#' @param metrics_list List of metrics
#' @param model_type binary/continuous/both.
#' @return A list of all registered metrics.
#' @export
filter_metrics <- function(metrics_list, model_type) {
    checkmate::assert_choice(model_type, c("binary", "continuous", "both"))
    l <- list()
    for (m in metrics_list) {
        if (m$type == model_type || m$type == "both") {
            l[[m$name]] <- m
        }
    }
    l
}


#' Get name description map
#'
#' @keywords internal
get_name_description_map <- function(metrics_list) {
    names <- c()
    descriptions <- c()
    for (m in metrics_list) {
        names <- c(names, m$name)
        descriptions <- c(descriptions, m$description)
    }
    data.table(name = names, description = descriptions)
}


#' Initialize default metrics
#'
#' @keywords internal
get_default_metrics <- function() {
    metrics <- list()
    metrics <- register_metric(metrics, "brier", brier_score, type = "binary", description = "Brier ~ score")
    metrics <- register_metric(metrics, "roc_auc", roc_auc, type = "binary", description = "ROC ~ AUC")
    metrics <- register_metric(metrics, "r2_nagelkerke", r2_nagelkerke, type = "binary", description = "Nagelkerke ~ italic(R)^2")
    metrics <- register_metric(metrics, "r2_liability", r2_liability, type = "binary", description = "Liability ~ scale ~ italic(R)^2")
    metrics <- register_metric(metrics, "aic", aic, type = "both", description = "AIC")
    metrics <- register_metric(metrics, "r2", r2, type = "continuous", description = "italic(R)^2")
    metrics <- register_metric(metrics, "mse", brier_score, type = "continuous", description = "MSE") # brier=mse for binary
    metrics <- register_metric(metrics, "mae", mae, type = "continuous", description = "MAE")
    metrics
}
