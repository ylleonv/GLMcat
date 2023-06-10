## code to prepare `disturbed_dreams`

disturbed_dreams <- read.csv("~/Desktop/Test package/data/Severity of Disturbed Dreams.csv")

# Wide to long
library(tidyr)
dreams_d1 <- gather(disturbed_dreams, Level, Total, Not.severe:Very.severe)

# Grouped to ungrouped
library(vcdExtra)
DisturbedDreams <- expand.dft(dreams_d1, freq = "Total")

DisturbedDreams$Level <- as.ordered(DisturbedDreams$Level)


# usethis::use_data(DisturbedDreams, compress = "bzip2", overwrite = TRUE)
