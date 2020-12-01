## -----------------------------------------------------------------------------
# devtools::load_all()
library(GLMcat)

## -----------------------------------------------------------------------------
data("TravelChoice")
head(TravelChoice)
str(TravelChoice)

## -----------------------------------------------------------------------------
exp_8.4 <- Discrete_CM(
  formula = choice ~ hinc + gc + invt,
  case_id = "indv",
  alternatives = "mode",
  reference = "air",
  data = TravelChoice,
  alternative_specific = c("gc", "invt"),
  distribution = "logistic")

## -----------------------------------------------------------------------------
summary(exp_8.4)

## -----------------------------------------------------------------------------
(constant_model <- Discrete_CM(
  formula = choice ~ 1 ,
  case_id = "indv",
  alternatives = "mode",
  reference = c("air", "train", "bus", "car"),
  data = TravelChoice,
  distribution = "logistic"
))

(car_0 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = c("air", "train", "bus", "car"),
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "logistic"
))

## -----------------------------------------------------------------------------
(air_30 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "air",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 3
))

## -----------------------------------------------------------------------------
(bus_30 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "bus",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 30
))

## -----------------------------------------------------------------------------
(car_.2 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "car",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 0.2
))

## -----------------------------------------------------------------------------
(train_1.35 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "train",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 1.35
))

## -----------------------------------------------------------------------------
(table_4 <- Discrete_CM(
  formula = choice ~ ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "car",
  alternative_specific = "ttme",
  data = TravelChoice,
  distribution = "logistic"
))
# AIC(table_4)

## -----------------------------------------------------------------------------
(car_.2 <- Discrete_CM(
  formula = choice ~ hinc[air] + psize[air] + gc + ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "car",
  alternative_specific = c("gc", "ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = 0.45
))
# logLik(car_.2)

## -----------------------------------------------------------------------------
(car_8 <- Discrete_CM(
  formula = choice ~  ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "car",
  alternative_specific = c("ttme"),
  data = TravelChoice,
  distribution = "logistic"
))
# logLik(car_8)

## -----------------------------------------------------------------------------
(car_.45 <- Discrete_CM(
  formula = choice ~  ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "car",
  alternative_specific = c("ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = .45
))
# logLik(car_.45)

## -----------------------------------------------------------------------------
(car_.05 <- Discrete_CM(
  formula = choice ~  ttme,
  case_id = "indv",
  alternatives = "mode",
  reference = "car",
  alternative_specific = c("ttme"),
  data = TravelChoice,
  distribution = "student",
  freedom_degrees = .05
))
# logLik(car_.05)

## -----------------------------------------------------------------------------
dat_p <- TravelChoice[1:8,]

## -----------------------------------------------------------------------------
# (adj_1 <- Discrete_CM(
#   formula = choice ~ hinc + psize + gc + ttme,
#   case_id = "indv",
#   alternatives = "mode",
#   reference = c("train"),
#   alternative_specific = c("gc", "ttme"),
#   data = dat_p,
#   distribution = "logistic",
#   ratio = "adjacent"
# ))
# 
# (adj_2 <- Discrete_CM(
#   formula = choice ~ hinc[car] + psize[air] + gc + ttme,
#   case_id = "indv",
#   alternatives = "mode",
#   reference = c("train"),
#   alternative_specific = c("gc", "ttme"),
#   data = dat_p,
#   distribution = "logistic",
#   ratio = "adjacent"
# ))
# 
# (adj_3 <- Discrete_CM(
#   formula = choice ~ hinc[air] + psize[bus] + gc + ttme,
#   case_id = "indv",
#   alternatives = "mode",
#   reference = c("air"),
#   alternative_specific = c("gc", "ttme"),
#   data = dat_p,
#   distribution = "logistic",
#   ratio = "adjacent"
# ))

## -----------------------------------------------------------------------------
# (adj_4 <- Discrete_CM(
#   formula = choice ~ hinc[air] + psize[bus] + gc + ttme,
#   case_id = "indv",
#   alternatives = "mode",
#   reference = c("air"),
#   alternative_specific = c("gc", "ttme"),
#   data = TravelChoice,
#   distribution = "logistic",
#   ratio = "adjacent"
# ))
# 
# (adj_5 <- Discrete_CM(
#   formula = choice ~ hinc[air] + psize[bus] + gc + ttme,
#   case_id = "indv",
#   alternatives = "mode",
#   reference = c("air"),
#   alternative_specific = c("gc", "ttme"),
#   data = TravelChoice,
#   distribution = "cauchit",
#   ratio = "adjacent"
# ))
# 
# (ref_4 <- Discrete_CM(
#   formula = choice ~ hinc[air] + psize[bus] + gc + ttme,
#   case_id = "indv",
#   alternatives = "mode",
#   reference = c("air"),
#   alternative_specific = c("gc", "ttme"),
#   data = TravelChoice,
#   distribution = "logistic",
#   ratio = "reference"
# ))
# 
# (ref_5 <- Discrete_CM(
#   formula = choice ~ hinc[air] + psize[bus] + gc + ttme,
#   case_id = "indv",
#   alternatives = "mode",
#   reference = c("air"),
#   alternative_specific = c("gc", "ttme"),
#   data = TravelChoice,
#   distribution = "cauchit",
#   ratio = "reference"
# ))

