test_that("glmcat works (known example)", {
  expect_equal(
    as.numeric(logLik(glmcat(formula = Level ~ Age, data = DisturbedDreams,
                  categories_order = c("Not.severe", "Severe.1", "Severe.2",
                                       "Very.severe"), ratio = "adjacent", cdf = "gompertz"))),
    -278.9917,
    tolerance=1e-5)
})

library(ordinal)

test_that("glmcat works (comparison with ordinal package)", {

  fm1 <- clm(rating ~ temp * contact, data = wine)
  glm_cum1 <- glmcat(rating ~ temp * contact, data = wine,
                     ratio = "cumulative")
  expect_equal(
    as.numeric(logLik(fm1)), as.numeric(logLik(glm_cum1)),
    tolerance=1e-5)
})



