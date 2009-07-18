# ====================================================================
#           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
#              Oak Ridge National Laboratory, Oak Ridge TN.
#        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
#                            Srikanth Yoginath.
#     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
#     Contact: Guruprasad Kora; (865) 576-6210; koragh@ornl.gov
#                 Computer Science and Mathematics Division
#             Oak Ridge National Laboratory, Oak Ridge TN 37831 
#                   (C) 2004 All Rights Reserved
#
#                              NOTICE
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby granted
# provided that the above copyright notice appear in all copies and
# that both the copyright notice and this permission notice appear in
# supporting documentation.
#
# Neither the Oak Ridge National Laboratory nor the Authors make any
# representations about the suitability of this software for any
# purpose.  This software is provided ``as is'' without express or
# implied warranty.
#
# RScaLAPACK (http://www.aspect-sdm.org/Parallel-R) was funded
# as part of the Scientific Data Management Center
# (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
# Scientific Discovery through Advanced Computing (DOE SciDAC) program
# (http://www.scidac.org ). 
# ======================================================================
# Function PRCOMP using scalapack functions to evaluate linear algebra problems
# uses sla.svd

sla.prcomp <- function(x,  retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, NPROWS=0, NPCOLS=0, MB=16)
{

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	dimX = dim(x);

	if ( dimX[1] < NPROWS || dimX[1] < NPCOLS || dimX[2] < NPROWS || dimX[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Check block size for proper distribution
	MB = sla.checkBlockSize(dim(x), MB)


    x <- as.matrix(x)
    x <- scale(x, center = center, scale = scale.)
    s <- sla.svd(x, nu = 0, NPROWS=NPROWS, NPCOLS=NPCOLS, MB=MB)
        s$vt=t(s$vt)
    if (!is.null(tol)) {
        rank <- sum(s$d > (s$d[1] * tol))
        if (rank < ncol(x))
            s$v <- s$v[, 1:rank, drop = FALSE]
    }
    s$d <- s$d/sqrt(max(1, nrow(x) - 1))
    dimnames(s$v) <- list(colnames(x), paste("PC", seq(len = ncol(s$v)),
        sep = ""))
    r <- list(sdev = s$d, rotation = s$v)
    if (retx)
        r$x <- x %*% s$v
    class(r) <- "prcomp"
    r
}

# Function PRINCOMP using scalapack functions to evaluate linear algebra problems
# uses sla.eigen
sla.princomp <- function(x, ...) UseMethod("sla.princomp")
                                                                                
## use formula to allow update() to be used.
sla.princomp.formula <- function(formula, data = NULL, subset, na.action, NPROWS=0, NPCOLS=0, MB=16, ...)
{


    mt <- terms(formula, data = data)
    if(attr(mt, "response") > 0) stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    na.act <- attr(mf, "na.action")
    x <- model.matrix(mt, mf)
    res <- sla.princomp.default(x,NPROWS=NPROWS,NPCOLS=NPCOLS,MB=MB, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("princomp")
    res$call <- cl
    if(!is.null(na.act)) {
        res$na.action <- na.act
        if(!is.null(sc <- res$scores))
            res$scores <- napredict(na.act, sc)
    }
    res
}

sla.princomp.default <-
    function(x, cor = FALSE, scores = TRUE, covmat = NULL,
             subset = rep(TRUE, nrow(as.matrix(x))),
               NPROWS=0, NPCOLS=0, MB=16, ...)
{

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	dimX = dim(x);

	if ( dimX[1] < NPROWS || dimX[1] < NPCOLS || dimX[2] < NPROWS || dimX[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Check block size for proper distribution
	MB = sla.checkBlockSize(dim(x), MB)

    cl <- match.call()
    cl[[1]] <- as.name("princomp")
    z <- if(!missing(x)) as.matrix(x)[subset, , drop = FALSE]
    if (is.list(covmat)) {
        if(any(is.na(match(c("cov", "n.obs"), names(covmat)))))
            stop("covmat is not a valid covariance list")
        cv <- covmat$cov
        n.obs <- covmat$n.obs
        cen <- covmat$center
    } else if(is.matrix(covmat)) {
        cv <- covmat
        n.obs <- NA
        cen <- NULL
    } else if(is.null(covmat)){
        dn <- dim(z)
        if(dn[1] < dn[2])
            stop("princomp can only be used with more units than variables")
        covmat <- cov.wt(z)             # returns list, cov() does not
        n.obs <- covmat$n.obs
        cv <- covmat$cov * (1 - 1/n.obs)# for S-PLUS compatibility
        cen <- covmat$center
    } else stop("covmat is of unknown type")
    if (cor) {
        sds <- sqrt(diag(cv))
        cv <- cv/(sds %o% sds)
    }
     edc <- sla.eigen(cv, NPROWS=NPROWS, NPCOLS=NPCOLS,MB=MB)

#    edc <- eigen(cv, symmetric = TRUE)
#    ev <- edc$values
    ev <- rev(edc$values)
    if (any(neg <- ev < 0)) { # S-PLUS sets all := 0
        ## 9 * : on Solaris found case where 5.59 was needed (MM)
        if (any(ev[neg] < - 9 * .Machine$double.eps * ev[1]))
            stop("covariance matrix is not non-negative definite")
        else
            ev[neg] <- 0
    }
    cn <- paste("Comp.", 1:ncol(cv), sep = "")
    names(ev) <- cn
    dimnames(edc$vectors) <- if(missing(x))
        list(dimnames(cv)[[2]], cn) else list(dimnames(x)[[2]], cn)
    sdev <- sqrt(ev)
    sc <- if (cor) sds else rep(1, ncol(cv))
    names(sc) <- colnames(cv)
    scr <- if (scores && !missing(x))
        scale(z, center = TRUE, scale = sc) %*% edc$vectors
    if (is.null(cen)) cen <- rep(NA, nrow(cv))
    edc <- list(sdev = sdev,

loadings = structure(edc$vectors, class="loadings"),
                center = cen, scale = sc, n.obs = n.obs,
                scores = scr, call = cl)
    ## The Splus function also return list elements factor.sdev,
    ## correlations and coef, but these are not documented in the help.
    ## coef seems to equal load.  The Splus function also returns list
    ## element terms which is not supported here.
    class(edc) <- "princomp"
    edc
}

# Function FACTANAL using scalapack functions to evaluate linear algebra problems
# uses sla.svd and sla.eigen

sla.factanal<-
    function (x, factors, data = NULL, covmat = NULL, n.obs = NA,
             subset, na.action, start = NULL,
              scores = c("none", "regression", "Bartlett"),
              rotation = "sla.varimax",
              control = NULL, NPROWS=0, NPCOLS=0, MB=16, ...)
{

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	if ( is.matrix(x) )
	{

		dimX = dim(x);

		if ( dimX[1] < NPROWS || dimX[1] < NPCOLS || dimX[2] < NPROWS || dimX[2] < NPCOLS )
		{
			stop("Input matrix is smaller than process grid");
		}
		# Check block size for proper distribution
		MB = sla.checkBlockSize(dim(x), MB)
	}



    sortLoadings <- function(Lambda)
    {
        cn <- colnames(Lambda)
        Phi <- attr(Lambda, "covariance")
        ssq <- apply(Lambda, 2, function(x) -sum(x^2))
        Lambda <- Lambda[, order(ssq), drop = FALSE]
        colnames(Lambda) <- cn
        neg <- colSums(Lambda) < 0
        Lambda[, neg] <- -Lambda[, neg]
        if(!is.null(Phi)) {
            unit <- ifelse(neg, -1, 1)
            attr(Lambda, "covariance") <-
                unit %*% Phi[order(ssq), order(ssq)] %*% unit
        }
        Lambda
    }
    cl <- match.call()
    na.act <- NULL
    if (is.list(covmat)) {
        if (any(is.na(match(c("cov", "n.obs"), names(covmat)))))
            stop("covmat is not a valid covariance list")
        cv <- covmat$cov
        n.obs <- covmat$n.obs
        have.x <- FALSE
    }
    else if (is.matrix(covmat)) {
        cv <- covmat
        have.x <- FALSE
    }
    else if (is.null(covmat)) {
        if(missing(x)) stop("neither x nor covmat supplied")
        have.x <- TRUE
        if(inherits(x, "formula")) {
            mt <- terms(x, data = data)
            if(attr(mt, "response") > 0)
                stop("response not allowed in formula")
            attr(mt, "intercept") <- 0
            mf <- match.call(expand.dots = FALSE)
            names(mf)[names(mf) == "x"] <- "formula"
            mf$factors <- mf$covmat <- mf$scores <- mf$start <-
                mf$rotation <- mf$control <- mf$... <- NULL
            mf[[1]] <- as.name("model.frame")
            mf <- eval(mf, parent.frame())
            na.act <- attr(mf, "na.action")
            z <- model.matrix(mt, mf)
        } else {
            z <- as.matrix(x)
            if(!missing(subset)) z <- z[subset, , drop = FALSE]
        }
        covmat <- cov.wt(z)
        cv <- covmat$cov
        n.obs <- covmat$n.obs
    }
    else stop("covmat is of unknown type")
    scores <- match.arg(scores)
    if(scores != "none" && !have.x)
        stop("requested scores without an x matrix")
    sds <- sqrt(diag(cv))
    cv <- cv/(sds %o% sds)
    p <- ncol(cv)
    dof <- 0.5 * ((p - factors)^2 - p - factors)
    if(dof < 0)
        stop(paste(factors, "factors is too many for", p, "variables"))

    cn <- list(nstart = 1, trace = FALSE, lower = 0.005)
    cn[names(control)] <- control
    more <- list(...)[c("nstart", "trace", "lower", "opt", "rotate")]
    if(length(more)) cn[names(more)] <- more

    if(is.null(start)) {
        #start <- (1 - 0.5*factors/p)/diag(solve(cv))
        start <- (1 - 0.5*factors/p)/diag(sla.solve(cv,NPROWS=NPROWS,NPCOLS=NPCOLS,MB=MB))
        if((ns <- cn$nstart) > 1)
            start <- cbind(start, matrix(runif(ns-1), p, ns-1, byrow=TRUE))
    }
    start <- as.matrix(start)
    if(nrow(start) != p) stop(paste("start must have", p, "rows"))
    nc <- ncol(start)
    if(nc < 1) stop("no starting values supplied")

    best <- Inf
    for (i in 1:nc) {
        nfit <- sla.factanal.fit.mle(cv, factors, start[, i],
                                 max(cn$lower, 0), cn$opt, NPROWS=NPROWS, NPCOLS=NPCOLS, MB=MB)
        if(cn$trace)
            cat("start", i, "value:", format(nfit$criteria[1]),
                "uniqs:", format(as.vector(round(nfit$uniquenesses, 4))), "\n")
        if(nfit$converged && nfit$criteria[1] < best) {
            fit <- nfit
            best <- fit$criteria[1]
        }
    }
    if(best == Inf) stop("Unable to optimize from these starting value(s)")
    load <- fit$loadings
    if(rotation != "none") {
        rot <- do.call(rotation, c(list(load), cn$rotate))
        load <- if(is.list(rot)) rot$loadings else rot
    }
    fit$loadings <- sortLoadings(load)
    class(fit$loadings) <- "loadings"
    fit$na.action <- na.act
    if(have.x && scores != "none") {
        Lambda <- fit$loadings
        zz <- scale(z, TRUE, TRUE)
        switch(scores,
               regression = {
                   #sc <- zz %*% solve(cv, Lambda)
                   sc <- zz %*% sla.solve(cv, Lambda,NPROWS=NPROWS,NPCOLS=NPCOLS,MB=MB)
                   if(!is.null(Phi <- attr(Lambda, "covariance")))
                       sc <- sc %*% Phi
               },
               Bartlett = {
                   d <- 1/fit$uniquenesses
                   tmp <- t(Lambda * d)
                   #sc <- t(solve(tmp %*% Lambda, tmp %*% t(zz)))
                   sc <- t(sla.solve(tmp %*% Lambda, tmp %*% t(zz),NPROWS=NPROWS,NPCOLS=NPCOLS,MB=MB))
               })
        rownames(sc) <- rownames(z)
        colnames(sc) <- colnames(Lambda)
        if(!is.null(na.act)) sc <- napredict(na.act, sc)
        fit$scores <- sc
    }
    if(!is.na(n.obs) && dof > 0) {
        fit$STATISTIC <- (n.obs - 1 - (2 * p + 5)/6 -
                     (2 * factors)/3) * fit$criteria["objective"]
        fit$PVAL <- pchisq(fit$STATISTIC, dof, lower.tail = FALSE)
    }
    fit$n.obs <- n.obs
    fit$call <- cl
    fit
}


sla.factanal.fit.mle <-
function(cmat, factors, start=NULL, lower = 0.005, control = NULL,NPROWS=0,NPCOLS=0,MB=16, ...)
{

FAout <- function(Psi, S, q)
    {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE)
        #E <- sla.eigen(Sstar)
        L <- E$vectors[, 1:q, drop = FALSE]
        load <- L %*% diag(sqrt(pmax(E$values[1:q] - 1, 0)), q)
        diag(sqrt(Psi)) %*% load
    }
    FAfn <- function(Psi, S, q)
    {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE, only.values = TRUE)
        #E <- sla.eigen(Sstar)
        e <- E$values[-(1:q)]
        e <- sum(log(e) - e) - q + nrow(S)
##        print(round(c(Psi, -e), 5))  # for tracing
        -e
    }
    FAgr <- function(Psi, S, q)
    {
        sc <- diag(1/sqrt(Psi))
        Sstar <- sc %*% S %*% sc
        E <- eigen(Sstar, symmetric = TRUE)
        #E <- sla.eigen(Sstar)
        L <- E$vectors[, 1:q, drop = FALSE]
        load <- L %*% diag(sqrt(pmax(E$values[1:q] - 1, 0)), q)
        load <- diag(sqrt(Psi)) %*% load
        g <- load %*% t(load) + diag(Psi) - S
        diag(g)/Psi^2
    }
    p <- ncol(cmat)
    if(is.null(start))
        #start <- (1 - 0.5*factors/p)/diag(solve(cmat))
        start <- (1 - 0.5*factors/p)/diag(sla.solve(cmat,NPROWS=NPROWS,NPCOLS=NPCOLS,MB=MB))
    res <- optim(start, FAfn, FAgr, method = "L-BFGS-B",
                 lower = lower, upper = 1,
                 control = c(list(fnscale=1,
                 parscale = rep(0.01, length(start))), control),
                 q = factors, S = cmat)
    Lambda <- FAout(res$par, cmat, factors)
    dimnames(Lambda) <- list(dimnames(cmat)[[1]],
                             paste("Factor", 1:factors, sep = ""))
    p <- ncol(cmat)
    dof <- 0.5 * ((p - factors)^2 - p - factors)
    un <- res$par
    names(un) <- colnames(cmat)
    class(Lambda) <- "loadings"
    ans <- list(converged = res$convergence == 0,
                loadings = Lambda, uniquenesses = un,
                correlation = cmat,
                criteria = c(objective = res$value, counts = res$counts),
                factors = factors, dof = dof, method = "mle")
    class(ans) <- "factanal"
    ans
}



sla.varimax <- function(x, normalize = TRUE, eps = 1e-5,NPROWS=0,NPCOLS=0,MB=16)
{


    nc <- ncol(x)
    if(nc < 2) return(x)
    if(normalize) {
        sc <- sqrt(drop(apply(x, 1, function(x) sum(x^2))))
        x <- x/sc
    }
    p <- nrow(x)
    TT <- diag(nc)
    d <- 0
    for(i in 1:1000) {
        z <- x %*% TT
        B  <- t(x) %*% (z^3 - z %*% diag(drop(rep(1, p) %*% z^2))/p)
#        sB <- La.svd(B)
        sB <- sla.svd(B,NPROWS=NPROWS,NPCOLS=NPCOLS,MB=MB)
        TT <- sB$u %*% sB$vt
        dpast <- d
        d <- sum(sB$d)
        if(d < dpast * (1 + eps)) break
    }
    z <- x %*% TT
    if(normalize) z <- z * sc
    dimnames(z) <- dimnames(x)
    list(loadings = z, rotmat = TT)
}

sla.promax <- function(x, m = 4, NPROWS=0, NPCOLS=0, MB=16)
{

# Check for illigal values

	if ( NPROWS < 0 )
	{
		stop("NPROWS cannot be a negative value.");
	}

	if ( NPCOLS < 0 )
	{
		stop("NPCOLS cannot be a negative value.");
	}

	if ( MB < 0 )
	{
		stop("MB cannot be a negative value.");
	}

	dimX = dim(x);

	if ( dimX[1] < NPROWS || dimX[1] < NPCOLS || dimX[2] < NPROWS || dimX[2] < NPCOLS )
	{
		stop("Input matrix is smaller than process grid");
	}

# Check block size for proper distribution
	MB = sla.checkBlockSize(dim(x), MB)

    if(ncol(x) < 2) return(x)
    dn <- dimnames(x)
    xx <- sla.varimax(x)
    x <- xx$loadings
    Q <- x * abs(x)^(m-1)
    U <- lm.fit(x, Q)$coefficients
    d <- diag(sla.solve(t(U) %*% U,NPROWS=NPROWS,NPCOLS=NPCOLS,MB=MB))
#    d <- diag(solve(t(U) %*% U))
    U <- U %*% diag(sqrt(d))
    dimnames(U) <- NULL
    z <- x %*% U
    U <- xx$rotmat %*% U
    dimnames(z) <- dn
    list(loadings = z, rotmat = U)
}

