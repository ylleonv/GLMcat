
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GLMcat

<!-- badges: start -->

<!-- badges: end -->

The goal of GLMcat is to …

## Installation

You can install the released version of GLMcat from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GLMcat")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ylleonv/GLMcat")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(GLMcat)
data(DisturbedDreams)
Multinomial_mod <- GLMcat(formula = Level ~ Age, categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
data = DisturbedDreams, distribution = "logistic")
```

The model `summary` :

``` r
summary(Multinomial_mod)
#>                         Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) Not.severe -2.454438   0.845593 -2.9026    0.0037 ** 
#> (Intercept) Severe.1   -0.554640   0.891015 -0.6225    0.5336    
#> (Intercept) Severe.2   -1.124641   0.916512 -1.2271    0.2198    
#> Age Not.severe          0.309988   0.078044  3.9719 7.129e-05 ***
#> Age Severe.1            0.059972   0.085817  0.6988    0.4847    
#> Age Severe.2            0.112281   0.086843  1.2929    0.1960    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
