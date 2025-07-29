\name{psme}
\alias{psme}
\alias{EvalSmooth}
\alias{GetExtrema}

\title{Fit Penalized splines linear mixed effects models}

\description{
Fits penalized splines linear mixed effects models by leveraging the known equivalence between penalized splines and mixed effects models. \code{psme} re-parameterizes penalized spline smooth terms from \code{mgcv}-style formulas as linear mixed effects models using \code{lme4} as the computational engine. This allows smoothing parameters to be estimated as variance components.
}

\usage{
psme(mgcv.form, data, knots = NULL)

EvalSmooth(mgcv.smooth, new.x)

GetExtrema(x, y)
}

\arguments{
  \item{mgcv.form}{A formula following \code{mgcv} syntax, specifying fixed effects and smooth terms via \code{s()}. For example, \code{y ~ x1 + s(x2)}.}
  \item{data}{A data frame containing all variables in the formula. Rows with \code{NA} values should be removed beforehand, as they are not handled by the function.}
  \item{knots}{An optional list of knot locations for the smooth terms, passed directly to \code{\link[mgcv]{smoothCon}}.}
  \item{mgcv.smooth}{A single smooth object taken from the \code{smooth} list of a fitted \code{psme} object.}
  \item{new.x}{A numeric vector of new data points at which to evaluate the smooth.}
  \item{x}{A numeric vector of x-coordinates for a curve.}
  \item{y}{A numeric vector of y-coordinates for a curve.}
}

\details{
The \code{psme} function transforms each penalized spline into a combination of fixed and random effects, which \code{lme4} can then estimate. A penalized spline with model matrix \eqn{X} and penalty matrix \eqn{S} has a penalty \eqn{\lambda \beta'S\beta}. This penalty is treated as a variance component in a mixed model.

The strength of the penalty is controlled by a smoothing parameter, which is estimated as a variance component within the mixed-effects framework. This often involves an L2 penalty on the differences of adjacent B-spline coefficients, which encourages smoother functions by penalizing "wiggliness".

The fitting process involves several key steps:
\enumerate{
  \item \strong{Formula Parsing}: The function uses \code{\link[mgcv]{interpret.gam}} to parse the formula.
  
  \item \strong{Smooth Term Processing}: Each smooth term is processed based on its type:
  \itemize{
    \item \strong{Singleton Smooth} (e.g., \code{s(x)}).
    \item \strong{Factor-Smooth Interaction} (e.g., \code{s(x, f, bs="fs")}).
    \item \strong{Factor \code{by} Variable} (e.g., \code{s(x, by=f)}).
  }
  
  \item \strong{Model Fitting}: A new formula is constructed for \code{lme4}, and the model is fitted using \code{lme4}'s low-level functions to accommodate the custom basis from \code{mgcv}.
  
  \item \strong{Coefficient Recovery}: Coefficients from the \code{lme4} fit are transformed back to the original spline basis.
}

\code{EvalSmooth} is a utility for predicting from a single fitted smooth term. For more complex predictions (e.g., with \code{by} variables), use \code{\link[mgcv]{PredictMat}} directly.

\code{GetExtrema} is a helper to find local minima and maxima of a curve.
}

\value{
\code{psme} returns a list containing:
\item{pform}{The parsed fixed effects formula used in the \code{lme4} model.}
\item{pcoef}{A named vector of fixed effects coefficients from the \code{lme4} fit.}
\item{smooth}{A list, where each element corresponds to an original smooth term from \code{mgcv.form}. Each element contains details about the smooth, including its type, and importantly, its estimated \code{coefficients} (the value of the smooth function at the basis points).}
\item{lme4.fit}{The raw fitted model object from \code{lme4}, an object of class \code{merMod}.}

\code{EvalSmooth} returns a vector of predicted values.

\code{GetExtrema} returns a list with four elements: \code{min.x}, \code{min.y}, \code{max.x}, and \code{max.y}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\seealso{
\code{\link[mgcv]{gam}}, \code{\link[lme4]{lmer}}
}

\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{splines}

\examples{
# A simple example to go here.
# require(psme)
}
