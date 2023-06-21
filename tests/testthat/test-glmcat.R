# Test if glmcat works with ADJACENT gompertz CDF
test_that("glmcat works with ADJACENT gompertz CDF", {
  expected_loglik <- -278.9917
  tolerance <- 1e-5

  result <- glmcat(
    formula = Level ~ Age,
    data = DisturbedDreams,
    categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
    ratio = "adjacent",
    cdf = "gompertz"
  )

  actual_loglik <- as.numeric(logLik(result))

  expect_equal(actual_loglik, expected_loglik, tolerance = tolerance)
})

# Test if glmcat works with CUMULATIVE cauchy CDF
test_that("glmcat works with CUMULATIVE cauchy CDF", {
  expected_loglik <- -279.8666
  tolerance <- 1e-5

  result <- glmcat(
    formula = Level ~ Age,
    data = DisturbedDreams,
    categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
    ratio = "cumulative",
    cdf = "cauchy"
  )

  actual_loglik <- as.numeric(logLik(result))

  expect_equal(actual_loglik, expected_loglik, tolerance = tolerance)
})

# Test if glmcat works with SEQUENTIAL normal CDF
test_that("glmcat works with SEQUENTIAL normal CDF", {
  expected_loglik <- -280.5465
  tolerance <- 1e-5

  result <- glmcat(
    formula = Level ~ Age,
    data = DisturbedDreams,
    categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
    ratio = "sequential",
    cdf = "normal"
  )

  actual_loglik <- as.numeric(logLik(result))

  expect_equal(actual_loglik, expected_loglik, tolerance = tolerance)
})

# Test if glmcat works with reference logistic CDF
test_that("glmcat works with cumulative logistic CDF", {
  expected_loglik <- -278.4682
  tolerance <- 1e-5

  result <- glmcat(
    formula = Level ~ Age,
    data = DisturbedDreams,
    ratio = "cumulative",
    cdf = "logistic"
  )

  actual_loglik <- as.numeric(logLik(result))

  expect_equal(actual_loglik, expected_loglik, tolerance = tolerance)
})

# Test if glmcat assumes logistic as default CDF if not specified
test_that("glmcat assumes logistic as default CDF if not specified", {
  result <- glmcat(
    formula = Level ~ Age,
    data = DisturbedDreams,
    ref_category = "Very.severe",
    ratio = "cumulative"
  )

  expect_true(result$cdf[1] == "logistic")
})

# Test if glmcat throws an error when ratio is not specified
test_that("glmcat throws an error when ratio is not specified", {
  expect_error(glmcat(
    formula = Level ~ Age,
    data = DisturbedDreams,
    ref_category = "Very.severe"
  ))
})

library(ordinal)

# Test if glmcat produces equivalent results to ordinal package
test_that("glmcat produces equivalent results to ordinal package", {
  wine <- ordinal::wine

  fm1 <- clm(rating ~ temp * contact, data = wine)
  glm_cum1 <- glmcat(rating ~ temp * contact, data = wine, ratio = "cumulative")

  expected_loglik <- as.numeric(logLik(fm1))
  actual_loglik <- as.numeric(logLik(glm_cum1))
  tolerance <- 1e-5

  expect_equal(actual_loglik, expected_loglik, tolerance = tolerance)
})
