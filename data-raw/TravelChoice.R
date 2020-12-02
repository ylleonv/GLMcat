## code to prepare `TravelChoice`

TravelChoice1 <- read.csv("~/Downloads/TableF18-2.csv")
indv <- rep(1:210, 4)
indv <- sort(indv)
choice <- rep(c("air", "train", "bus", "car"), 210)
mode <- TravelChoice1$MODE == 1

TravelChoice <- as.data.frame(cbind(indv, choice, mode, TravelChoice1[, -1]))
colnames(TravelChoice) <- c("indv", "mode", "choice", "ttme", "invc", "invt", "gc", "hinc", "psize")
TravelChoice$indv <- as.factor(TravelChoice$indv)

# usethis::use_data(TravelChoice,compress = "bzip2", overwrite = T)
