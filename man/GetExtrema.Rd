\name{GetExtrema}

\alias{GetExtrema}

\title{
Identify all the local extrema of a curve with a grid search
}

\description{
Approximately identify all the local extrema of a 1-dimensional curve \eqn{y = f(x)}, by finding all the local minima and all the local maxima in a sequence of \eqn{(x, y)} values, discretely observed or sampled from the curve.
}

\usage{
GetExtrema(x, y)
}

\arguments{
\item{x}{A numeric vector of x-coordinates for a curve.}
\item{y}{A numeric vector of y-coordinates for a curve.}
}

\value{
A list with four elements:
\item{min.x}{The x-coordinates of the local minima.}
\item{min.y}{The y-coordinates of the local minima.}
\item{max.x}{The x-coordinates of the local maxima.}
\item{max.y}{The y-coordinates of the local maxima.}
}

\seealso{
\code{\link{EvalSmooth}}
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
x <- seq(0, 2 * pi, length.out = 100)
y <- sin(x)
ex <- GetExtrema(x, y)
plot(x, y, type = "l")
points(ex$min.x, ex$min.y, pch = 19, col = 2)
points(ex$max.x, ex$max.y, pch = 19, col = 3)
}

