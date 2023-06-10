test_that("discrete_cm works (known example)", {
  logistic_car_alt <- discrete_cm(formula = choice ~
                                    hinc[air] + psize[air] + gc + ttme, case_id = "indv",
                                  alternatives = "mode", reference = "car", data = TravelChoice,
                                  alternative_specific = c("gc", "ttme"), cdf = "logistic")
  expect_equal(as.numeric(logLik(logistic_car_alt)), -185.9149, tolerance=1e-5)
})
