\name{psme}

\alias{psme}

\title{
Fit penalized splines mixed-effects models
}

\description{
Fit penalized splines mixed-effects models by leveraging the known equivalence between penalized splines and mixed-effects models. \code{psme} reparameterizes each penalized spline term specified in a \emph{mgcv}-style model formula as linear mixed-effects model terms, and fits such equivalent model using \emph{lme4} as the computational engine.
}

\usage{
psme(mgcv.form, data, knots = NULL)
}

\arguments{
\item{mgcv.form}{a \emph{mgcv}-style model formula that specify a penalized splines mixed-effects model in the fashion of penalized splines additive models.}
\item{data}{a data frame with all variables in the formula, and any rows with NA must be mannually removed.}
\item{knots}{an optional named list providing knots.}
}

\details{
A penalized splines mixed-effects model is a special case of additive models, most suitable for modelling longitudinal data. Such model typically contains some population smooth functions and a group of subject smooth functions, both nonlinear in a variable (for example, time). The function relies on \emph{mgcv} for constructing smooth functions as penalized splines, and users should specify these terms via a \emph{mgcv}-style model formula. Examples of population smooth functions include a single smooth function like \code{s(x)} and a factor `by' smooth term like \code{s(x, by = f)}, where the smooth function in \code{x} changes with the levels of a factor variable \code{f} (for example, gender). Subject smooth functions are set up as a factor-smooth interaction term like \code{s(x, f, bs = "fs")}. The model may also involve smooth functions of other covariate variables constructed using \code{s}, \code{te} and \code{ti}, as well as any parametric terms that are legitimate in a \code{lm} formula. In addition, users can customize the class and basis dimension of each penalized spline, as they can when using \emph{mgcv}.

To fit a penalized splines mixed-effects model, the function transforms each penalized spline into a combination of fixed and random effects, producing an equivalent linear mixed-effects model. This model is then estimated using \emph{lme4}'s low-level functions. Finally, the function backtransforms the estimated fixed-effects and predicted random-effects to the basis coefficients of the original penalized splines, so that the users can extract and plot the estimated smooth functions. In particular, the transformation from a penalized splines mixed-effects model to a linear mixed-effects model is performed with great caution to preserve the sparsity of the design matrices. As a result, the function is able to handle large longitudinal datasets with millions of observations.
}

\value{
a list containing:
\item{pform}{The parametric component of the model formula.}
\item{pcoef}{The estimated coefficients of the parametric component of the model.}
\item{smooth}{A list of constructed and estimated smooth functions corresponding to the smooth terms in the model formula.}
\item{lme4.fit}{The raw output of \emph{lme4}'s low-level functions.}
}

\seealso{
\code{\link[mgcv]{gam}}, \code{\link[mgcv]{gamm}}
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
library(psme)

## simulate a toy dataset of 50 subjects, each with a random quadratic trajectory
n.subjects <- 50
y <- x <- vector("list", n.subjects)
## a subject has 5 to 10 random observations
n <- sample(5:10, size = n.subjects, replace = TRUE)
a <- rnorm(n.subjects, mean = 1, sd = 0.2)
b <- rnorm(n.subjects, mean = 0, sd = 0.3)
c <- rnorm(n.subjects, mean = 0, sd = 0.4)
for (i in 1:n.subjects) {
  x_i <- sort(runif(n[i], -1, 1))
  y_i <- a[i] * x_i ^ 2 + b[i] * x_i + c[i]
  x[[i]] <- x_i; y[[i]] <- y_i
}
x <- unlist(x)
y <- unlist(y)
## add noise to the true trajectories
y <- y + rnorm(length(y), sd = 0.1)
## compose a dataset in long format
id <- rep.int(1:n.subjects, n)
dat <- data.frame(x = x, y = y, id = as.factor(id))

## fit a penalized splines mixed-effects model
mgcv.form <- y ~ s(x, bs = 'ps', k = 8) +
                 s(x, id, bs = 'fs', xt = 'ps', k = 8, m = c(2, 1))

\donttest{
fit <- psme(mgcv.form, data = dat)

## evaluate population smooth function and subject smooth functions
xp <- seq.int(-1, 1, length.out = 51)
pop.smooth <- EvalSmooth(fit$smooth[[1]], new.x = xp)
sub.smooth <- EvalSmooth(fit$smooth[[2]], new.x = xp)
smooth <- pop.smooth + sub.smooth + fit$pcoef[["(Intercept)"]]
matplot(xp, smooth, type = "l", lty = 1, xlab = "x", ylab = "fitted trajectories")
points(dat, pch = 20, col = 8)
}
}

