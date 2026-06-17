
test_that("effects_pgsmetrics accepts numeric outcome and categorical PGS", {
  d <- simulate_data(n = 100)

  d$pgscat <- ifelse(
    d$pgs1 <= 0.25, "<0.25",
    ifelse(d$pgs1 <= 0.5, "(0.25, 0.5]", ">0.5")
  )

  expect_no_error(
    effects_pgsmetrics(
      d,
      pgs = "pgscat",
      dep = "y",
      covars = c("age", "sex")
    )
  )
})

test_that("effects_pgsmetrics rejects non-0/1 binary outcome", {
  d <- simulate_data(n = 100)
  d$y_bin <- d$y_bin + 1

  expect_error(
    effects_pgsmetrics(d,
                       pgs = "pgs1",
                       dep = "y_bin"
    ),
    "Binary outcome variable must be coded as 0/1"
  )
})

test_that("effects_pgsmetrics validates PGS arguments", {
    d <- simulate_data(n = 100)
    
    expect_error(
      effects_pgsmetrics(
        d,
        pgs = "pgs10",
        dep = "y_bin",
        covars = c("age", "sex")
      ),
      "not present"
    )
    
    expect_error(
      effects_pgsmetrics(
        d,
        pgs = c("pgs1", NA),
        dep = "y_bin",
        covars = c("age", "sex")
      ),
      "not present"
    )
    
    expect_warning(
      res <- effects_pgsmetrics(
        d,
        pgs = c("pgs1", "pgs1"),
        dep = "y_bin",
        covars = c("age", "sex")
      ),
      "Duplicate PGS"
    )
    
    expect_equal(nrow(res), 1)
    expect_equal(res$Model, "pgs1")
  })
  
test_that("effects_pgsmetrics validates covariate arguments", {
    d <- simulate_data(n = 100)
    
    expect_error(
      effects_pgsmetrics(
        d,
        pgs = "pgs1",
        dep = "y_bin",
        covars = c("age", NA)
      ),
      "not present"
    )
    
    expect_error(
      effects_pgsmetrics(
        d,
        pgs = "pgs1",
        dep = "y_bin",
        covars = "missing_covariate"
      ),
      "not present"
    )
    
    expect_warning(
      effects_pgsmetrics(
        d,
        pgs = "pgs1",
        dep = "y_bin",
        covars = c("age", "age")
      ),
      "Duplicate covariates"
    )
    
    expect_no_error(
      effects_pgsmetrics(
        d,
        pgs = "pgs1",
        dep = "y_bin",
        covars = NULL,
        missing = "warn",
        report_covars = TRUE
      )
    )
  })
  
test_that("effects_pgsmetrics validates dependent variable argument", {
    d <- simulate_data(n = 100)
    
    expect_error(
      effects_pgsmetrics(
        d,
        pgs = "pgs1",
        dep = c("dep1", "dep2"),
        covars = c("age", "sex")
      ),
      "Multiple Dependent variables"
    )
    
    expect_error(
      effects_pgsmetrics(
        d,
        pgs = "pgs1",
        dep = c("dep1", NA),
        covars = c("age", "sex")
      ),
      "Multiple Dependent variables"
    )
    
    expect_error(
      effects_pgsmetrics(
        d,
        pgs = "pgs1",
        dep = "missing_dep",
        covars = c("age", "sex")
      ),
      "not present"
    )
  })
  
test_that("effects_pgsmetrics returns expected columns for continuous outcome", {
  set.seed(789)
  
  n <- 500
  d <- data.table(
    covar1 = rnorm(n),
    pgs1 = rnorm(n),
    pgs2 = rnorm(n),
    dep = rnorm(n)
  )
  
  expected_cols <- c(
    "Model", "Term", "Estimate", "Std. Error", "t value", "Pr(>|t|)"
  )
  
  r_standard <- effects_pgsmetrics(
    data = d,
    pgs = c("pgs1", "pgs2"),
    dep = "dep",
    covars = "covar1",
    report_covars = FALSE
  )
  
  r_with_covars <- effects_pgsmetrics(
    data = d,
    pgs = c("pgs1", "pgs2"),
    dep = "dep",
    covars = "covar1",
    report_covars = TRUE
  )
  
  expect_equal(names(r_standard), expected_cols)
  expect_equal(names(r_with_covars), expected_cols)
  
  expect_equal(nrow(r_standard), 2)
  expect_true(all(r_standard$Term %in% c("pgs1", "pgs2")))
  
  expect_equal(nrow(r_with_covars), 4)
  expect_true(all(c("pgs1", "pgs2", "covar1") %in% r_with_covars$Term))
})

test_that("effects_pgsmetrics returns expected columns for binary outcome", {
  set.seed(789)
  
  n <- 500
  d <- data.table(
    covar1 = rnorm(n),
    pgs1 = rnorm(n),
    pgs2 = rnorm(n),
    dep = rnorm(n)
  )
  d$dep_bin <- as.numeric(d$dep > 0)
  
  expected_cols <- c(
    "Model", "Term", "Estimate", "Std. Error", "z value", "Pr(>|z|)"
  )
  
  r_standard <- effects_pgsmetrics(
    data = d,
    pgs = c("pgs1", "pgs2"),
    dep = "dep_bin",
    covars = "covar1",
    report_covars = FALSE
  )
  
  r_with_covars <- effects_pgsmetrics(
    data = d,
    pgs = c("pgs1", "pgs2"),
    dep = "dep_bin",
    covars = "covar1",
    report_covars = TRUE
  )
  
  expect_equal(names(r_standard), expected_cols)
  expect_equal(names(r_with_covars), expected_cols)
  
  expect_equal(nrow(r_standard), 2)
  expect_true(all(r_standard$Term %in% c("pgs1", "pgs2")))
  
  expect_equal(nrow(r_with_covars), 4)
  expect_true(all(c("pgs1", "pgs2", "covar1") %in% r_with_covars$Term))
})

test_that("effects_pgsmetrics handles missing = drop", {
  d <- simulate_data(n = 100)
  
  d$y[1:10] <- NA
  
  expect_message(
    res <- effects_pgsmetrics(
      d,
      pgs = c("pgs1", "pgs2"),
      dep = "y",
      covars = NULL,
      missing = "drop"
    ),
    "Removed 10 rows"
  )
  
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
})

test_that("effects_pgsmetrics rejects binary outcome not coded 0/1", {
  d <- simulate_data(n = 100)
  
  d$y_bin <- d$y_bin + 1
  
  expect_error(
    effects_pgsmetrics(
      d,
      pgs = c("pgs1", "pgs2"),
      dep = "y_bin"
    ),
    "Binary outcome variable must be coded as 0/1"
  )
})

test_that("effects_pgsmetrics rejects single-value outcome", {
  d <- simulate_data(n = 100)
  
  d$y_bin <- 1
  
  expect_error(
    effects_pgsmetrics(
      d,
      pgs = "pgs1",
      dep = "y_bin",
      covars = c("age", "sex"),
      report_covars = TRUE
    ),
    "Outcome variable must have >1 values"
  )
})

test_that("effects_pgsmetrics rejects empty PGS vector", {
  d <- simulate_data(n = 100)
  
  expect_error(
    effects_pgsmetrics(
      d,
      pgs = character(0),
      dep = "y_bin",
      covars = c("age", "sex")
    ),
    "at least one PGS"
  )
})

test_that("effects_pgsmetrics output identifies model and term correctly", {
  d <- simulate_data(n = 100)
  
  res <- effects_pgsmetrics(
    d,
    pgs = c("pgs1", "pgs2"),
    dep = "y_bin",
    covars = c("age", "sex"),
    report_covars = FALSE
  )
  
  expect_equal(res$Model, c("pgs1", "pgs2"))
  expect_equal(res$Term, c("pgs1", "pgs2"))
})

test_that("effects_pgsmetrics report_covars returns covariate rows per model", {
  d <- simulate_data(n = 100)
  
  res <- effects_pgsmetrics(
    d,
    pgs = c("pgs1", "pgs2"),
    dep = "y_bin",
    covars = c("age", "sex"),
    report_covars = TRUE
  )
  
  expect_true(all(c("pgs1", "pgs2") %in% res$Model))
  
  expect_equal(
    table(res$Model)[["pgs1"]],
    3
  )
  
  expect_equal(
    table(res$Model)[["pgs2"]],
    3
  )
  
  expect_true(all(c("age", "sex") %in% res$Term))
})
