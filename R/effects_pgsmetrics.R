#' Create coefficient table containing PGS effect size adjusted for covariates. 
#' 
#' @param data A data.table containing necessary data: PGS, dependent variable, and covariates.
#' @param pgs Character vector. Metrics are computed for these pgs. Must exist in `data`.
#' @param dep Character. Name of the dependent variable in `data`. Binary variable must be 0/1 (otherwise continuous assumed).
#' @param covars Character vector. Names of covariate to adjust for. Must be in columns in `data`.
#' @param missing Warns if missing values in `data`. Use missing="drop" to drop rows with missing values. However, we recommend imputing missing values instead. Default: "warn".
#'
#' @return A data.frame with calculated covariate-adjusted effect size of PGS. 
#'   \item{coef_table}{A data.frame with calculated covariate-adjusted effect size of PGS.}
#'   
#' @examples
#' # Simulate data
#' data <- simulate_data(n = 1000, n_pgs = 3)
#'
#' # raw PGS coefficients 
#' coef_table <- coef.pgsmetrics(data,
#'     pgs = c("pgs1", "pgs2", "pgs3"),
#'     dep = "y_bin",
#'     covars = c("age", "sex")
#' )
#' # print(coef_table)
#'
#' # Obtain coefficients for top decile of PGS
#' data$pgs1.top10 <- 1*(data$pgs1 >= quantile(data$pgs1, 0.9, na.rm = TRUE))
#' 
#' coef <- coef.pgsmetrics(data, 
#'     pgs = c("pgs1", "pgs1.top10"), 
#'     dep = "y_bin", 
#'     covars = c("age", "sex"),
#'  )
#' # print(coef)
#' 
#' @export
effects_pgsmetrics <- function(
    data,
    pgs,
    dep,
    covars = NULL,
    missing = "warn",
    report_covars = FALSE) { 
  
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }
  
  ##  Verify columns 
  if(any(duplicated(pgs))) {
    warning("Duplicate PGS found - removing duplicates")
    pgs <- pgs[!duplicated(pgs)]
  }
  
  if(any(duplicated(covars))) {
    warning("Duplicate covariates - removing duplicates")
    covars <- covars[!duplicated(covars)]
  }
  
  if(length(dep) > 1) {
    stop("Multiple Dependent variables are not allowed")
  }

  required <- c(pgs, dep, covars)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("The following variables are not present in 'data': %s",
                 paste(missing_cols, collapse = ", ")),
      call. = FALSE)
  }

  data <- handle_missing_values(data, pgs, dep, covars, missing)
  
  quant <- determine_outcome_type(data[[dep]])
  if(quant){
    fam <- "gaussian"
  } else {
    fam <- "binomial"
  }
  
  ## Estimate effect sizes of PGS adjusted for covar
  ## - optionally include covar effect sizes
  
  npgs <- length(pgs)
  coef_list <- list()

  if (npgs == 0) {
    stop("Must have at least one PGS column.")
  }

  for(i in 1:npgs){
    
    full_pgs <- as.formula(
      paste(dep, "~", paste(c(covars, pgs[i]) , collapse = " + "))
    )
    
    fit_full <- glm(full_pgs, data=data, family=fam)
    if(report_covars) {
      index <- which(grepl(paste(c(pgs[i], covars), collapse="|"), names(fit_full$coefficients)))
    } else {
      index <- grep(pgs[i], names(fit_full$coefficients))
    }

    ##  Extract terms & add columns for model (PGS name) & all terms
    ct <- as.data.frame(
      summary(fit_full)$coefficients[index, , drop = FALSE]
    )
    
    ct$Term <- rownames(ct)
    rownames(ct) <- NULL
    ct$Model <- pgs[i]
    
    coef_list[[i]] <- ct
  }
  coef_table <- do.call('rbind', coef_list) 
  
  return(coef_table)
}




