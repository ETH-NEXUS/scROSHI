
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scROSHI

<!-- badges: start -->

<!-- badges: end -->

Robust supervised hierarchical identification of single cells

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ETH-NEXUS/scROSHI")
```

## Example

This is a basic example for the Wilcox test function:

``` r
library(scROSHI)
f_my_wilcox_test(rnorm(10,1,2),rnorm(10,5,2))
#> [1] 1
```
