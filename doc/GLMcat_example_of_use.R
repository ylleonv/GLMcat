## -----------------------------------------------------------------------------
library(GLMcat)
# library(gtools); library(vcdExtra); library(grid)
library(gridExtra)
library(gtools)
library(dplyr)
# #library(stringr)
library(ggplot2)
library(tidyr)

## -----------------------------------------------------------------------------
data("DisturbedDreams")
summary(DisturbedDreams)
str(DisturbedDreams)

## -----------------------------------------------------------------------------
mod1 <- GLMcat(formula = Level ~ Age,
  categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  data = DisturbedDreams, distribution = "logistic")

## -----------------------------------------------------------------------------
summary(mod1)

## -----------------------------------------------------------------------------
nobs_glmcat(mod1)

## -----------------------------------------------------------------------------
coef(mod1)

## -----------------------------------------------------------------------------
logLik(mod1)

## -----------------------------------------------------------------------------
AIC(mod1)

## -----------------------------------------------------------------------------
BIC(mod1)

## -----------------------------------------------------------------------------
predict_glmcat(mod1, data = DisturbedDreams[1:5,], type = "prob")

## -----------------------------------------------------------------------------
predict_glmcat(mod1, data = DisturbedDreams[1:5,], type = "linear.predict_glmcator")

## -----------------------------------------------------------------------------
predict_glmcat(mod1, data = DisturbedDreams[1:5,], type = "cum.prob")

## -----------------------------------------------------------------------------
mod2 <- GLMcat(formula = Level ~ Age, proportional_effects = "Age",
  categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  data = DisturbedDreams, distribution = "logistic")

summary(mod2)
logLik(mod2)

## -----------------------------------------------------------------------------
predict_glmcat(mod2, data = DisturbedDreams[1:5,], type = "prob")

## -----------------------------------------------------------------------------
# predict_glmcat(mod2, data = DisturbedDreams[1:5,], type = "linear.predict_glmcator")

## -----------------------------------------------------------------------------
predict_glmcat(mod2, data = DisturbedDreams[1:5,], type = "cum.prob")

## -----------------------------------------------------------------------------
mod3 <- GLMcat(formula = Level ~ Age,
  categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  data = DisturbedDreams, distribution = "cauchit")
summary(mod3)
logLik(mod3)

## -----------------------------------------------------------------------------
mod1_a <- GLMcat(formula = Level ~ Age,
  categories_order = c("Not.severe", "Severe.1", "Very.severe", "Severe.2"),
  data = DisturbedDreams, distribution = "logistic")
summary(mod1_a)
logLik(mod1_a)

## -----------------------------------------------------------------------------
mod2_a <- GLMcat(formula = Level ~ Age, proportional_effects = "Age",
  categories_order = c("Not.severe", "Severe.1", "Very.severe", "Severe.2"),
  data = DisturbedDreams, distribution = "logistic")
summary(mod2_a)
logLik(mod2_a)

## -----------------------------------------------------------------------------
mod3_a <- GLMcat(formula = Level ~ Age,
  categories_order = c("Not.severe", "Severe.1", "Very.severe", "Severe.2"),
  data = DisturbedDreams, distribution = "cauchit")
summary(mod3_a)
logLik(mod3_a)

## -----------------------------------------------------------------------------
adj_logi <- GLMcat(formula = Level ~ Age, ratio = "adjacent", distribution = "logistic",
  categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  data = DisturbedDreams)
logLik(adj_logi)
summary(adj_logi)

## -----------------------------------------------------------------------------
matrix(c(1,-1,0,0,1,-1,0,0,1), nrow = 3)

## -----------------------------------------------------------------------------
adj_cauchit <- GLMcat(formula = Level ~ Age, ratio = "adjacent", distribution = "cauchit",
  categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  data = DisturbedDreams)
logLik(adj_cauchit)
summary(adj_cauchit)

## -----------------------------------------------------------------------------
adj_cauchit_rev <- GLMcat(formula = Level ~ Age, ratio = "adjacent", distribution = "cauchit",
  categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
  data = DisturbedDreams)
logLik(adj_cauchit_rev)
summary(adj_cauchit_rev)

## -----------------------------------------------------------------------------
adj_gompertz <- GLMcat(formula = Level ~ Age, ratio = "adjacent", distribution = "gompertz",
  categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  proportional_effects = c("(Intercept)", "Age"),
  data = DisturbedDreams)
logLik(adj_gompertz)
summary(adj_gompertz)

## -----------------------------------------------------------------------------
adj_gompertz_rev <- GLMcat(formula = Level ~ Age, ratio = "adjacent", distribution = "gompertz",
  categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
  proportional_effects = c("(Intercept)", "Age"),
  data = DisturbedDreams)
logLik(adj_gompertz_rev)
summary(adj_gompertz_rev)

## -----------------------------------------------------------------------------
adj_gumbel <- GLMcat(formula = Level ~ Age, ratio = "adjacent", distribution = "gumbel",
  categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
  proportional_effects = c("(Intercept)", "Age"),
  data = DisturbedDreams)
logLik(adj_gumbel)
summary(adj_gumbel)

## -----------------------------------------------------------------------------
cum_gompertz <- GLMcum(formula = Level ~ Age,
  categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  proportional_effects = "Age",
  data = DisturbedDreams,
  distribution = "gompertz")
cum_gompertz$`Log-likelihood`
(cum_gompertz$coefficients)

## -----------------------------------------------------------------------------
seq_gompertz <- GLMcat(formula = Level ~ Age, ratio = "adjacent", distribution = "gompertz",
  categories_order = c("Very.severe", "Severe.2", "Severe.1", "Not.severe"),
  proportional_effects = c("(Intercept)", "Age"),
  data = DisturbedDreams)
logLik(seq_gompertz)
summary(seq_gompertz)

knitr::opts_chunk$set(fig.width=9, fig.height=5.5,
  message=FALSE, warning = FALSE, echo = FALSE,
  tidy = TRUE, tidy.opts = list(comment = FALSE))

## -----------------------------------------------------------------------------
#

## -----------------------------------------------------------------------------
all_permutations <- permutations(
  v = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  repeats.allowed = F, n = 4, r = 4
)

plot <- list()

for (dist in c("logistic", "normal", "cauchit", "student")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    l <- GLMcat(formula = Level ~ Age, ratio = "reference",
      distribution = dist, freedom_degrees = 3,
      categories_order = all_permutations[element, ],
      # proportional_effects = c("(Intercept)", "Age"),
      data = DisturbedDreams)
    Log_lik_Vec[element] <- logLik(l)
  }
  Log_lik_Vec
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Reference, ", dist, ", complete")
  plot[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}

plot_adj <- grid.arrange(plot[["logistic"]], plot[["normal"]],
  plot[["cauchit"]], plot[["student"]],
  ncol = 2, nrow = 2)

# plot_adj

## -----------------------------------------------------------------------------
plot <- list()
for (dist in c("logistic", "normal", "cauchit", "student")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    l <- GLMcat(formula = Level ~ Age, ratio = "reference",
      distribution = dist,
      categories_order = all_permutations[element, ],
      proportional_effects = c("(Intercept)", "Age"),
      data = DisturbedDreams)
    Log_lik_Vec[element] <- logLik(l)
  }
  Log_lik_Vec
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Reference, ", dist, ", proportional")
  plot[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}

plot_adj <- grid.arrange(plot[["logistic"]], plot[["normal"]],
  plot[["cauchit"]], plot[["student"]],
  ncol = 2, nrow = 2
)

# plot_adj

## -----------------------------------------------------------------------------
all_permutations <- permutations(
  v = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  repeats.allowed = F, n = 4, r = 4
)
plot <- list()
for (dist in c("logistic", "normal", "cauchit", "gompertz")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    l <- GLMcat(formula = Level ~ Age, ratio = "adjacent",
      distribution = dist,
      categories_order = all_permutations[element, ],
      # proportional_effects = c("(Intercept)", "Age"),
      data = DisturbedDreams)
    Log_lik_Vec[element] <- logLik(l)
  }
  Log_lik_Vec
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Adjacent, ", dist, ", complete")
  plot[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}

plot_adj <- grid.arrange(plot[["logistic"]], plot[["normal"]],
  plot[["cauchit"]], plot[["gompertz"]],
  ncol = 2, nrow = 2
)

## -----------------------------------------------------------------------------
all_permutations <- permutations(
  v = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
  repeats.allowed = F, n = 4, r = 4
)
plot <- list()
for (dist in c("logistic", "normal", "cauchit", "gompertz")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    l <- GLMcat(formula = Level ~ Age, ratio = "adjacent",
      distribution = dist,
      categories_order = all_permutations[element, ],
      proportional_effects = c("(Intercept)", "Age"),
      data = DisturbedDreams)
    Log_lik_Vec[element] <- logLik(l)
  }
  Log_lik_Vec
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Adjacent, ", dist, ", proportional")
  plot[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}

plot_adj <- grid.arrange(plot[["logistic"]], plot[["normal"]],
  plot[["cauchit"]], plot[["gompertz"]],
  ncol = 2, nrow = 2
)

## -----------------------------------------------------------------------------
plot3 <- list()
for (dist in c("logistic", "normal", "cauchit", "gompertz")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    skip_to_next <- FALSE
    tryCatch(
      {
        l <- GLMcat(formula = Level ~ Age, ratio = "sequential",
          distribution = dist,
          categories_order = all_permutations[element, ],
          # proportional_effects = c("(Intercept)", "Age"),
          data = DisturbedDreams)
        Log_lik_Vec[element] <- logLik(l)
      },
      error = function(e) {
        Log_lik_Vec[element] <- NA
        skip_to_next <<- TRUE
      }
    )
    if (skip_to_next) {
      next
    } else {
      Log_lik_Vec[element] <- l$`Log-likelihood`
    }
  }
  Log_lik_Vec[element] <- l$`Log-likelihood`
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Sequential, ", dist, ", complete")
  plot3[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}
plot_cum_com <- grid.arrange(plot3[["logistic"]], plot3[["normal"]],
  plot3[["cauchit"]], plot3[["gompertz"]],
  ncol = 2, nrow = 2
)

## -----------------------------------------------------------------------------
plot3 <- list()
for (dist in c("logistic", "normal", "cauchit", "gompertz")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    skip_to_next <- FALSE
    tryCatch(
      {
        l <- GLMcat(formula = Level ~ Age, ratio = "sequential",
          distribution = dist,
          categories_order = all_permutations[element, ],
          proportional_effects = c("(Intercept)", "Age"),
          data = DisturbedDreams)

      },
      error = function(e) {
        Log_lik_Vec[element] <- logLik(l)
        skip_to_next <<- TRUE
      }
    )
    if (skip_to_next) {
      next
    } else {
      Log_lik_Vec[element] <- l$`Log-likelihood`
    }
  }
  Log_lik_Vec[element] <- l$`Log-likelihood`
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Sequential, ", dist, ", proportional")
  plot3[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}
plot_cum_com <- grid.arrange(plot3[["logistic"]], plot3[["normal"]],
  plot3[["cauchit"]], plot3[["gompertz"]],
  ncol = 2, nrow = 2
)

## -----------------------------------------------------------------------------
for (dist in c("logistic", "normal", "cauchit", "gompertz")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    l <- GLMcum(
      formula = Level ~ Age,
      categories_order = all_permutations[element, ],
      data = DisturbedDreams,
      distribution = dist
    )
    Log_lik_Vec[element] <- l$`Log-likelihood`
  }
  Log_lik_Vec
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Cumulative, ", dist, ", complete")
  plot[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}
plot_cum <- grid.arrange(plot[["logistic"]], plot[["normal"]],
  plot[["cauchit"]], plot[["gompertz"]],
  ncol = 2, nrow = 2
)

## -----------------------------------------------------------------------------
for (dist in c("logistic", "normal", "cauchit", "gompertz")) {
  Log_lik_Vec <- NA
  for (element in 1:nrow(all_permutations)) {
    l <- GLMcum(
      formula = Level ~ Age,
      categories_order = all_permutations[element, ],
      proportional_effects = c("Age"),
      data = DisturbedDreams,
      distribution = dist
    )
    Log_lik_Vec[element] <- l$`Log-likelihood`
  }
  Log_lik_Vec
  all_permutations_num <- permutations(v = c("1", "2", "3", "4"), repeats.allowed = F, n = 4, r = 4)
  names <- as.data.frame(all_permutations_num) %>% unite("z", remove = FALSE, sep = "")
  to_plot <- data.frame("LogLik" = Log_lik_Vec, "Permutation" = names[, 1])
  to_plot$LogLik <- round(to_plot$LogLik, digits = 4)
  to_plot$Permutation <- as.factor(to_plot$Permutation)
  to_plot$Distribution <- dist
  groups <- data.frame(gn = 1:length(unique(to_plot$LogLik)), LogLik = unique(to_plot$LogLik))
  to_plot <- left_join(to_plot, groups)
  to_plot <- to_plot[to_plot$LogLik >= -1000, ]
  title <- paste0("Cumulative, ", dist, ", Proportional")
  plot[[dist]] <- to_plot %>%
    arrange(-LogLik) %>%
    mutate(Permutation = factor(Permutation, levels = (Permutation))) %>%
    ggplot(aes(x = Permutation, y = LogLik)) +
    geom_point() +
    geom_line(aes(group = gn)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(title) +
    xlab("") +
    ylab("")
}
plot_cum <- grid.arrange(plot[["logistic"]], plot[["normal"]],
  plot[["cauchit"]], plot[["gompertz"]],
  ncol = 2, nrow = 2
)

# source('~/Desktop/Test package/GLMcat_examples.R', echo=TRUE)
# source('C:/Users/leonvelasco/Google Drive/PhD most updated/Test package 3008/Test package/GLMcat_examples.R', echo=TRUE)
# rmarkdown::render("C:/Users/leonvelasco/Google Drive/PhD most updated/Test package 3008/Test package/GLMcat_examples.R", "pdf_document")
# rmarkdown::render("/home/leonvelasco/Downloads/GLMcat_examples.R", "pdf_document")

