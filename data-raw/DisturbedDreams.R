## code to prepare `disturbed_dreams`

disturbed_dreams <- read.csv("~/Desktop/Test package/data/Severity of Disturbed Dreams.csv")

# Wide to long
library(tidyr)
dreams_d1 <- gather(dreams_d, Level, Total, Not.severe:Very.severe)

# Grouped to ungrouped
library(vcdExtra)
disturbed_dreams <- expand.dft(disturbed_dreams, freq = "Total")

# usethis::use_data("disturbed_dreams")
