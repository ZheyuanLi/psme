#' Convert a smooth specification to a mixed-effects model representation
#'
#' This function takes a smooth specification from `mgcv` and converts it into a
#' mixed-effects model representation that can be used with `lme4`.
#'
#' @param object An object of class `mgcv.smooth`.
#' @return A list containing the random effects design matrix (`Xr`), the fixed
#'   effects design matrix (`Xf`), and the penalty matrix (`P`).
#' @keywords internal
ss2me <- function (object) {
  if (!inherits(object, "mgcv.smooth")) {
    stop("'object' is not a singleton s() term!")
  }
  X <- object$X
  S <- object$S[[1]]
  n <- nrow(X)
  p <- ncol(X)
  r <- as.integer(object$rank)
  f <- as.integer(object$null.space.dim)
  R <- chol.default(base::crossprod(X))
  Omega <- forwardsolve(R, S, upper.tri = TRUE, transpose = TRUE)
  Omega <- forwardsolve(R, t.default(Omega), upper.tri = TRUE, transpose = TRUE)
  ei <- eigen(Omega, symmetric = TRUE)
  D <- ei$values[1:r]
  U <- ei$vectors
  P <- backsolve(R, U)
  Z <- forwardsolve(R, t.default(X), upper.tri = TRUE, transpose = TRUE)
  Z <- t.default(crossprod(U, Z))
  d <- 1 / sqrt(D)
  P[, 1:r] <- P[, 1:r] * rep(d, each = p)
  Xr <- Z[, 1:r] * rep(d, each = n)
  if (f) {
    Xf <- Z[, (r + 1):p, drop = FALSE]
  } else {
    Xf <- NULL
  }
  if (f) {
    list(Xr = Xr, Xf = Xf, P = P)
  } else {
    list(Xr = Xr, P = P)
  }
}

#' Create a design matrix for a factor interaction
#'
#' This function creates a design matrix for the interaction between a continuous
#' variable and a factor.
#'
#' @param X A numeric matrix.
#' @param fctr A factor variable.
#' @param split.by.level A logical value indicating whether to split the design
#'   matrix by the levels of the factor.
#' @return A sparse design matrix or a list of sparse design matrices.
#' @keywords internal
Xfctr <- function (X, fctr, split.by.level = FALSE) {
  if (!is.factor(fctr)) stop("'fctr' is not a factor variable!")
  nlev <- nlevels(fctr)
  n <- nrow(X)
  p <- ncol(X)
  .i <- seq.int(0, n - 1L)
  .j <- as.integer(fctr) - 1L
  i <- rep.int(.i, p)
  j <- c(outer(p * .j, 0L:(p - 1L), `+`))
  .M <- methods::new("dgTMatrix", i = i, j = j, x = as.double(X), Dim = c(n, p * nlev))
  M <- methods::as(.M, "CsparseMatrix")
  if (!split.by.level) return(M)
  Mlev <- vector("list", nlev)
  ind <- 1:p
  for (k in 1:nlev) {
    Mlev[[k]] <- M[, ind, drop = FALSE]
    ind <- ind + p
  }
  Mlev
}

#' Convert a factor smooth interaction to a random-effects representation
#'
#' This function takes a factor smooth interaction from `mgcv` and converts it
#' into a random-effects representation that can be used with `lme4`.
#'
#' @param object An object of class `mgcv.smooth`.
#' @param fctr A factor variable.
#' @return A list containing the random effects design matrices (`Xr`) and the
#'   penalty matrix (`P`).
#' @keywords internal
fs2re <- function (object, fctr) {
  if (!is.null(object$Xf)) f <- ncol(object$Xf) else f <- 0L
  Xr <- vector("list", 1L + f)
  Xr[[1L]] <- Xfctr(object$Xr, fctr)
  if (f) {
    M <- Xfctr(object$Xf, fctr)
    ind <- seq.int(1L, by = f, length.out = nlevels(fctr))
    for (k in 1:f) {
      Xr[[k + 1L]] <- M[, ind, drop = FALSE]
      ind <- ind + 1L
    }
  }
  list(Xr = Xr, P = object$P)
}

#' Convert a factor-by smooth to a mixed-effects model representation
#'
#' This function takes a factor-by smooth from `mgcv` and converts it into a
#' mixed-effects model representation that can be used with `lme4`.
#'
#' @param object An object of class `mgcv.smooth`.
#' @param fctr A factor variable.
#' @return A list containing the random effects design matrices (`Xr`), the fixed
#'   effects design matrix (`Xf`), and the penalty matrix (`P`).
#' @keywords internal
ssfby2me <- function (object, fctr) {
  Xr <- Xfctr(object$Xr, fctr, split.by.level = TRUE)
  if (!is.null(object$Xf)) {
    list(Xr = Xr, Xf = Xfctr(object$Xf, fctr), P = object$P)
  } else {
    list(Xr = Xr, P = object$P)
  }
}

#' Fit a Penalized Splines Mixed-Effects Model
#'
#' This function fits a Gaussian additive model using `lme4`. It is intended for
#' large longitudinal datasets.
#'
#' @param mgcv.form A `mgcv`-style model formula.
#' @param data A data frame containing the variables in the formula. Note: please
#'   remove rows with `NA`.
#' @param knots An optional named list providing knots.
#' @return A list containing the following components:
#'   \item{pform}{The formula for the parametric part of the model.}
#'   \item{pcoef}{The coefficients for the parametric part of the model.}
#'   \item{smooth}{A list of smooth objects.}
#'   \item{lme4.fit}{The fitted `lme4` model.}
#' @author Zheyuan Li \email{zheyuan.li@@bath.edu}
#' @examples
#' require(psme)
#' ## A simple example
#' set.seed(123)
#' n <- 100
#' x <- runif(n)
#' y <- 2 * x + sin(2 * pi * x) + rnorm(n, 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- psme(y ~ s(x), data = dat)
#' summary(fit$lme4.fit)
#' @export
psme <- function (mgcv.form, data, knots = NULL) {
  n.obs <- nrow(data)
  parsed.mgcv.form <- mgcv::interpret.gam(mgcv.form)
  n.smooths <- length(parsed.mgcv.form$smooth.spec)
  if (n.smooths == 0L) {
    stop("There is no smooth term in the formula!")
  }
  smooth <- vector("list", n.smooths)
  for (i in 1:n.smooths) {
    spec <- parsed.mgcv.form$smooth.spec[[i]]
    if (inherits(spec, c("tensor.smooth.spec", "t2.smooth.spec"))) {
      stop("Sorry, but we can not process tensor product smooths at the momemnt!")
    }
    if (inherits(spec, "fs.smooth.spec")) {
      if (spec$dim == 1L) {
        stop("'fs' smooth should involve two or more variables!")
      }
      fctr.var <- which(unlist(lapply(data[spec$term], is.factor), use.names = FALSE))
      if (length(fctr.var) != 1L) {
        stop("'fs' smooth should have (only) one factor variable!")
      }
      spec$fctr <- spec$term[fctr.var]
      base.spec <- spec
      base.spec$term <- base.spec$term[-fctr.var]
      base.spec$dim <- base.spec$dim - 1
      if (is.null(base.spec$xt)) {
        class(base.spec) <- "tp.smooth.spec"
      } else {
        class(base.spec) <- paste0(base.spec$xt, ".smooth.spec")
      }
      base.spec$xt <- NULL
      base.spec$tensor.possible <- NULL
      sm <- mgcv::smooth.construct2(base.spec, data, knots)
      sm$lme4 <- fs2re(ss2me(sm), data[[spec$fctr]])
      sm$type <- "factor.smooth.interaction"
    } else if (is.factor(data[[spec$by]])) {
      by.var <- spec$by
      spec$by <- "NA"
      sm <- mgcv::smoothCon(spec, data, knots, TRUE, FALSE)[[1L]]
      sm$factor.by <- by.var
      sm$lme4 <- ssfby2me(ss2me(sm), data[[by.var]])
      sm$type <- "factor.by"
    } else {
      sm <- mgcv::smoothCon(spec, data, knots, TRUE, FALSE)[[1L]]
      sm$lme4 <- ss2me(sm)
      sm$type <- "singleton"
    } ## end if else
    smooth[[i]] <- sm
  }  ## end for
  lme4.form <- deparse1(parsed.mgcv.form$pf)
  kf <- 0L
  kr <- 0L
  Xf <- Xr <- Z <- list()
  for (i in 1:n.smooths) {
    sm <- smooth[[i]]
    if (sm$type == "singleton") {
      nm <- character(1L)
      kr <- kr + 1L
      nm[1L] <- Xr.label <- paste0("Xr.", kr)
      Xr[[kr]] <- rep_len(as.factor(1:ncol(sm$lme4$Xr)), n.obs)
      lme4.form <- sprintf("%s + (1 | %s)", lme4.form, Xr.label)
      Z[[kr]] <- methods::as(sm$lme4$Xr, "CsparseMatrix")
      new.lme4.component <- list(sm$lme4$Xr)
    } else if (sm$type %in% c("factor.smooth.interaction", "factor.by")) {
      nm <- character(length(sm$lme4$Xr))
      for (j in 1:length(sm$lme4$Xr)) {
        kr <- kr + 1L
        nm[j] <- Xr.label <- paste0("Xr.", kr)
        Xr[[kr]] <- rep_len(as.factor(1:ncol(sm$lme4$Xr[[j]])), n.obs)
        lme4.form <- sprintf("%s + (1 | %s)", lme4.form, Xr.label)
        Z[[kr]] <- sm$lme4$Xr[[j]]
      }
      new.lme4.component <- sm$lme4$Xr
    }
    if (!is.null(sm$lme4$Xf)) {
      kf <- kf + 1L
      Xf.label <- paste0("Xf.", kf)
      nm <- c(nm, Xf.label)
      Xf[[kf]] <- as.matrix(sm$lme4$Xf)
      lme4.form <- sprintf("%s + %s", lme4.form, Xf.label)
      new.lme4.component <- c(new.lme4.component, list(sm$lme4$Xf))
    }
    new.lme4.component <- c(new.lme4.component, list(sm$lme4$P))
    names(new.lme4.component) <- c(nm, "P")
    smooth[[i]]$lme4 <- new.lme4.component
  }  ## end for
  lme4.form <- stats::formula(lme4.form)
  names(Z) <- names(Xr) <- sprintf("Xr.%d", 1:kr)
  if (kf) names(Xf) <- sprintf("Xf.%d", 1:kf)
  lme4.vars <- c(data[all.vars(parsed.mgcv.form$pf)], Xf, Xr)
  ctrl <- lme4::lmerControl(check.nobs.vs.nlev = "warning", check.nobs.vs.nRE = "warning")
  parsed.lme4.form <- lme4::lFormula(lme4.form, lme4.vars, REML = TRUE, control = ctrl)
  reTrms <- parsed.lme4.form$reTrms
  Xr.names <- names(reTrms$flist)
  reTrms$Ztlist <- lapply(Z[Xr.names], Matrix::t)
  reTrms$Zt <- do.call(rbind, reTrms$Ztlist)
  reTrms$cnms <- Map(paste, reTrms$cnms, sub("Xr.", "", Xr.names, fixed = TRUE), sep = ".")
  parsed.lme4.form$reTrms <- reTrms
  devianceFunction <- do.call(lme4::mkLmerDevfun, parsed.lme4.form)
  optimizerOutput <- lme4::optimizeLmer(devianceFunction)
  lme4.fit <- lme4::mkMerMod(rho = environment(devianceFunction),
                             opt = optimizerOutput,
                             reTrms = parsed.lme4.form$reTrms,
                             fr = parsed.lme4.form$fr)
  bf <- lme4.fit@beta
  br <- lme4.fit@u
  theta <- lme4.fit@theta
  X.terms <- attr(stats::terms(parsed.lme4.form$fr), "varnames.fixed")[-1]
  X <- lme4.fit@pp$X
  X.names <- dimnames(X)[[2]]
  X.assign <- attr(X, "assign")
  for (i in 1:n.smooths) {
    sm <- smooth[[i]]
    nm <- names(sm$lme4)
    kr <- length(grep("Xr.", nm, fixed = TRUE))
    coef.r <- vector("list", kr)
    for (k in 1:kr) {
      Xr.label <- nm[k]
      j <- match(Xr.label, Xr.names)
      ind <- seq.int(reTrms$Gp[j] + 1L, reTrms$Gp[j + 1L])
      coef.r[[k]] <- br[ind] * theta[j]
    }
    if (sm$type == "singleton") {
      coef.r <- matrix(coef.r[[1L]])
    } else if (sm$type == "factor.smooth.interaction") {
      nlev <- nlevels(data[[sm$fctr]])
      coef.r <- lapply(coef.r, matrix, ncol = nlev)
      coef.r <- do.call(rbind, coef.r)
    } else {
      coef.r <- do.call(cbind, coef.r)
    }
    if (length(sm$lme4) == kr + 2L) {
      Xf.label <- names(sm$lme4)[kr + 1L]
      j <- match(Xf.label, X.terms)
      ind <- X.assign == j
      coef.f <- matrix(bf[ind], ncol = ncol(coef.r))
      omit <- !ind
      X.names <- X.names[omit]
      X.assign <- X.assign[omit]
      bf <- bf[omit]
      X.terms <- X.terms[-j]
    } else {
      coef.f <- matrix(0, nrow = 0, ncol = ncol(coef.r))
    }
    smooth[[i]]$coefficients <- sm$lme4$P %*% rbind(coef.r, coef.f)
  }
  names(bf) <- X.names
  list(pform = parsed.mgcv.form$pf, pcoef = bf, smooth = smooth, lme4.fit = lme4.fit)
}

#' Evaluate a smooth term
#'
#' This function evaluates a smooth term at new data points.
#'
#' @param mgcv.smooth A smooth object from a `psme` fit.
#' @param new.x A vector or matrix of new data points.
#' @return A vector or matrix of the evaluated smooth term.
#' @export
EvalSmooth <- function (mgcv.smooth, new.x) {
  dat <- data.frame(x = new.x)
  names(dat) <- mgcv.smooth$term
  X <- mgcv::PredictMat(mgcv.smooth, dat)
  b <- mgcv.smooth$coefficients
  y <- X %*% b
  if (ncol(y) == 1L) c(y) else y
}

#' Find the extrema of a curve
#'
#' This function finds the local minima and maxima of a curve.
#'
#' @param x A numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates.
#' @return A list with the following components:
#'   \item{min.x}{The x-coordinates of the local minima.}
#'   \item{min.y}{The y-coordinates of the local minima.}
#'   \item{max.x}{The x-coordinates of the local maxima.}
#'   \item{max.y}{The y-coordinates of the local maxima.}
#' @export
GetExtrema <- function (x, y) {
  turn <- diff.default(sign(diff.default(y)))
  maxInd <- which(turn == -2) + 1L
  minInd <- which(turn == 2) + 1L
  list(min.x = x[minInd], min.y = y[minInd],
       max.x = x[maxInd], max.y = y[maxInd])
}
