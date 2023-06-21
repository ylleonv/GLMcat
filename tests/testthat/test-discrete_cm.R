# Test if discrete_cm works with a known example using logistic CDF
test_that("discrete_cm works with a known example using logistic CDF", {
  expected_loglik <- -185.9149
  tolerance <- 1e-5

  result <- discrete_cm(
    formula = choice ~ hinc[air] + psize[air] + gc + ttme,
    case_id = "indv",
    alternatives = "mode",
    reference = "car",
    data = TravelChoice,
    alternative_specific = c("gc", "ttme"),
    cdf = "logistic"
  )

  actual_loglik <- as.numeric(logLik(result))

  expect_equal(actual_loglik, expected_loglik, tolerance = tolerance)
})

# Test if the 'df of the model' is correctly defined
test_that("Check correct definition of the 'df of the model'", {
  expected_df <- 7
  tolerance <- 1e-5

  result <- discrete_cm(
    formula = choice ~ hinc[air] + psize[air] + gc + ttme,
    case_id = "indv",
    alternatives = "mode",
    reference = "car",
    data = TravelChoice,
    alternative_specific = c("gc", "ttme"),
    cdf = list("student", 1)
  )

  actual_df <- result$`df of the model`

  expect_equal(actual_df, expected_df, tolerance = tolerance)
})

# Test if the 'student' distribution returns the same result as 'cauchy' when df = 1
test_that("Check that the 'student' is well defined, when df = 1 it should return the same result as 'cauchy'", {
  tolerance <- 1e-5
  result_student <- discrete_cm(
    formula = choice ~ hinc[air] + psize[air] + gc + ttme,
    case_id = "indv",
    alternatives = "mode",
    reference = "car",
    data = TravelChoice,
    alternative_specific = c("gc", "ttme"),
    cdf = list("student", 1)
  )

  result_cauchy <- discrete_cm(
    formula = choice ~ hinc[air] + psize[air] + gc + ttme,
    case_id = "indv",
    alternatives = "mode",
    reference = "car",
    data = TravelChoice,
    alternative_specific = c("gc", "ttme"),
    cdf = "cauchy"
  )

  expect_equal(as.numeric(logLik(result_student)),
               as.numeric(logLik(result_cauchy)),
               tolerance = tolerance)
})

# Test if different coefficients are obtained when normalization is set to TRUE
test_that("Check different results for the coefficients when normalization is set to TRUE", {
  result_normalization <- discrete_cm(
    formula = choice ~ hinc[air] + psize[air] + gc + ttme,
    case_id = "indv",
    alternatives = "mode",
    reference = "car",
    data = TravelChoice,
    alternative_specific = c("gc", "ttme"),
    cdf = list("student", 1),
    normalization = 0.8
  )

  result_default <- discrete_cm(
    formula = choice ~ hinc[air] + psize[air] + gc + ttme,
    case_id = "indv",
    alternatives = "mode",
    reference = "car",
    data = TravelChoice,
    alternative_specific = c("gc", "ttme"),
    cdf = list("student", 1)
  )

  coef_normalization <- summary(result_normalization, normalized = TRUE)$coefficients[1, 1]
  coef_default <- summary(result_default)$coefficients[1, 1]

  expect_false(coef_normalization == coef_default)
})
