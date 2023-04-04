\name{psme}

\alias{psme}

\title{
Penalized Splines Mixed-Effects Models
}

\description{
Fit a Gaussian additive model using lme4
}

\usage{
psme(mgcv.form, data, knots = NULL)
}

\arguments{
\item{mgcv.form}{a mgcv-style model formula.}
\item{data}{a data frame with all variables in the formula (note: please remove rows with NA).}
\item{knots}{an optional named list providing knots.}
}

\value{
to do
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
require(psme)
}
