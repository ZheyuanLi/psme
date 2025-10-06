\name{EvalSmooth}

\alias{EvalSmooth}

\title{
Evaluate a fitted smooth term in a fitted model
}

\description{
A helper function to evaluate the fitted population and subject smooth functions in a fitted penalized splines mixed-effects model.
}

\usage{
EvalSmooth(mgcv.smooth, new.x)
}

\arguments{
\item{mgcv.smooth}{A single smooth term taken from the \code{smooth} list of the fitted penalized splines mixed-effects model.}
\item{new.x}{A numeric vector of new covariate values at which to evaluate the smooth term.}
}

\value{
These are two common cases of population smooth functions. For a single smooth function like \code{s(x)}, the function returns a vector of predicted values. For a factor `by' smooth term like \code{s(x, by = f)}, the function returns a matrix of predicted values, where each column corresponds to a factor level.

For subject smooth functions specified as a factor-smooth interaction term \code{s(x, f, bs = "fs")}, the function returns a matrix of predicted values, where each column corresponds to an individual subject.
}

\seealso{
\code{\link{psme}}
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
## see examples for the psme() function
}
