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
  Qt <- forwardsolve(R, t.default(X), upper.tri = TRUE, transpose = TRUE)
  Z <- t.default(crossprod(U, Qt))
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

ssfby2me <- function (object, fctr) {
  Xr <- Xfctr(object$Xr, fctr, split.by.level = TRUE)
  if (!is.null(object$Xf)) {
    list(Xr = Xr, Xf = Xfctr(object$Xf, fctr), P = object$P)
  } else {
    list(Xr = Xr, P = object$P)
  }
}

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

EvalSmooth <- function (mgcv.smooth, new.x) {
  dat <- data.frame(x = new.x)
  names(dat) <- mgcv.smooth$term
  X <- mgcv::PredictMat(mgcv.smooth, dat)
  b <- mgcv.smooth$coefficients
  y <- X %*% b
  if (ncol(y) == 1L) c(y) else y
}

GetExtrema <- function (x, y) {
  turn <- diff.default(sign(diff.default(y)))
  maxInd <- which(turn == -2) + 1L
  minInd <- which(turn == 2) + 1L
  list(min.x = x[minInd], min.y = y[minInd],
       max.x = x[maxInd], max.y = y[maxInd])
}

