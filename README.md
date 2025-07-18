# psme: Penalized Splines Mixed-Effects Models

`psme` is an R package for fitting penalized splines mixed-effects models for large longitudinal datasets. The core computational routines are implemented in C++ and interfaced with R. The package is designed to be fast and memory-efficient, and can handle datasets with millions of observations.

## Installation

You can install the development version of `psme` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("ZheyuanLi/psme")
```

## Example

Here is a simple example of how to use `psme` to fit a penalized splines mixed-effects model:

```r
require(psme)

## A simple example
set.seed(123)
n <- 100
x <- runif(n)
y <- 2 * x + sin(2 * pi * x) + rnorm(n, 0, 0.5)
dat <- data.frame(x = x, y = y)
fit <- psme(y ~ s(x), data = dat)
summary(fit$lme4.fit)
```

## Documentation

To generate the documentation for this package, you will need to have `roxygen2` installed. You can install it from CRAN with:

```r
install.packages("roxygen2")
```

Once `roxygen2` is installed, you can generate the documentation by running the following command in your R console:

```r
roxygen2::roxygenise()
```
