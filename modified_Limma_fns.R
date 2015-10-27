# modifications:
#   added argument test_coefs: - list of vectors of coefficients, each specifying a group of parameters to test jointly
#   This is used only when a matrix of weights for each observation is provided.
#   In this case:
#       a new value, ts_sub is calculated, which are components of t-statistics except not divided by the residual SD
#       a new value, fs_sub is calculated for each gene x element of test_coefs which are the F statistics corresponding to this test, except not divided by the residual squeezeVar
#       ts_sub and fs_sub are used to calculate ts and Fs in ebayes_mod
lm.series_mod = function (M, design = NULL, ndups = 1, spacing = 1, weights = NULL,test_coefs= NULL) 
{
    M <- as.matrix(M)
    narrays <- ncol(M)
    if (is.null(design)) 
        design <- matrix(1, narrays, 1)
    else design <- as.matrix(design)
    nbeta <- ncol(design)
    coef.names <- colnames(design)
    if (is.null(coef.names)) 
        coef.names <- paste("x", 1:nbeta, sep = "")
    if (!is.null(weights)) {
        weights <- asMatrixWeights(weights, dim(M))
        weights[weights <= 0] <- NA
        M[!is.finite(weights)] <- NA
    }
    if (ndups > 1) {
        M <- unwrapdups(M, ndups = ndups, spacing = spacing)
        design <- design %x% rep(1, ndups)
        if (!is.null(weights)) 
            weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
    }
    ngenes <- nrow(M)
    stdev.unscaled <- beta <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M), 
        coef.names))
    NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights, 
        "arrayweights")))
    if (NoProbeWts) {
        if (is.null(weights)) 
            fit <- lm.fit(design, t(M))
        else {
            print('adsf')
            fit <- lm.wfit(design, t(M), weights[1, ])
            fit$weights <- NULL
        }
        if (fit$df.residual > 0) {
            if (is.matrix(fit$effects)) 
                fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 
                  1):narrays, , drop = FALSE]^2))
            else fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 
                1):narrays]^2))
        }
        else fit$sigma <- rep(NA, ngenes)
        fit$fitted.values <- fit$residuals <- fit$effects <- NULL
        fit$coefficients <- t(fit$coefficients)
        fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
        est <- fit$qr$pivot[1:fit$qr$rank]
        dimnames(fit$cov.coefficients) <- list(coef.names[est], 
            coef.names[est])
        stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)), 
            ngenes, fit$qr$rank, byrow = TRUE)
        fit$stdev.unscaled <- stdev.unscaled
        fit$df.residual <- rep.int(fit$df.residual, ngenes)
        dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
        fit$pivot <- fit$qr$pivot
        return(fit)
    }
    ts_sub <- beta <- stdev.unscaled
    fs_sub <- array(NA,c(ngenes,length(test_coefs)))
    sigma <- rep(NA, ngenes)
    df.residual <- rep(0, ngenes)

    for (i in 1:ngenes) {
        y <- as.vector(M[i, ])
        obs <- is.finite(y)
        if (sum(obs) > 0) {
            X <- design[obs, , drop = FALSE]
            y <- y[obs]
            if (is.null(weights)) {
                out <- lm.fit(X, y)
            }else {
                w <- as.vector(weights[i, obs])
                out <- lm.wfit(X, y, w)
            }
            est <- !is.na(out$coef)
            beta[i, ] <- out$coef
            cov.coefficients <- chol2inv(out$qr$qr, size = out$rank)
            stdev.unscaled[i, est] <- sqrt(diag(cov.coefficients))
            df.residual[i] <- out$df.residual
            ts_sub[i, est] <- beta[i, est] / stdev.unscaled[i, est]
            if(!is.null(test_coefs)){
                for(j in 1:length(test_coefs)){
                    coefs = test_coefs[[j]]
                    res = classifyTestsF(t(matrix(ts_sub[i,coefs])),cov2cor(cov.coefficients)[coefs,coefs],fstat.only = TRUE)
                    fs_sub[i,j] <- res       
                }
            }
            if (df.residual[i] > 0) 
                sigma[i] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
        }
    }
    QR <- qr(design)
    cov.coef <- chol2inv(QR$qr, size = QR$rank)
    est <- QR$pivot[1:QR$rank]
    dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
    list(coefficients = beta, stdev.unscaled = stdev.unscaled, 
        sigma = sigma, df.residual = df.residual, cov.coefficients = cov.coef, 
        pivot = QR$pivot, rank = QR$rank,ts_sub = ts_sub,fs_sub = fs_sub)
}


# modifications: 
#   passes this test_coefficients to the modified fitting function lm.series_mod
#   adds field test_coef_ranks to the returned fit object to be used for df1 downstream
lmFit_mod = function (object, design = NULL, ndups = 1, spacing = 1, block = NULL, 
    correlation, weights = NULL, method = "ls", test_coefs = NULL,...) 
{
    y <- getEAWP(object)
    if (is.null(design)) 
        design <- y$design
    if (is.null(design)) 
        design <- matrix(1, ncol(y$exprs), 1)
    else {
        design <- as.matrix(design)
        if (mode(design) != "numeric") 
            stop("design must be a numeric matrix")
        if (nrow(design) != ncol(y$exprs)) 
            stop("row dimension of design doesn't match column dimension of data object")
    }
    ne <- nonEstimable(design)
    if (!is.null(ne)) 
        cat("Coefficients not estimable:", paste(ne, collapse = " "), 
            "\n")
    if (missing(ndups) && !is.null(y$printer$ndups)) 
        ndups <- y$printer$ndups
    if (missing(spacing) && !is.null(y$printer$spacing)) 
        spacing <- y$printer$spacing
    if (missing(weights) && !is.null(y$weights)) 
        weights <- y$weights
    method <- match.arg(method, c("ls", "robust"))
    if (ndups > 1) {
        if (!is.null(y$probes)) 
            y$probes <- uniquegenelist(y$probes, ndups = ndups, 
                spacing = spacing)
        if (!is.null(y$Amean)) 
            y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean), 
                ndups = ndups, spacing = spacing), na.rm = TRUE)
    }
    if (method == "robust") 
        fit <- mrlm(y$exprs, design = design, ndups = ndups, 
            spacing = spacing, weights = weights, ...)
    else if (ndups < 2 && is.null(block)){
        fit <- lm.series_mod(y$exprs, design = design, ndups = ndups, 
            spacing = spacing, weights = weights,test_coefs=test_coefs)
    } else {
        if (missing(correlation)) 
            stop("the correlation must be set, see duplicateCorrelation")
        fit <- gls.series(y$exprs, design = design, ndups = ndups, 
            spacing = spacing, block = block, correlation = correlation, 
            weights = weights, ...)
    }
    if (NCOL(fit$coef) > 1) {
        i <- is.na(fit$coef)
        i <- apply(i[, 1] == i[, -1, drop = FALSE], 1, all)
        n <- sum(!i)
        if (n > 0) 
            warning("Partial NA coefficients for ", n, " probe(s)", 
                call. = FALSE)
    }
    fit$genes <- y$probes
    fit$Amean <- y$Amean
    fit$method <- method
    fit$design <- design
    fit$test_coef_ranks = sapply(test_coefs,length)
    new("MArrayLM", fit)
}

# modifications:
#   after calculating s2.post, uses this to scale the newly calculated ts_sub and fs_sub to form 
#       t statistics ts and F statistics Fs for each pre-specified set of coefficients 
ebayes_mod = function (fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
    trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) 
{
    coefficients <- fit$coefficients
    stdev.unscaled <- fit$stdev.unscaled
    sigma <- fit$sigma
    df.residual <- fit$df.residual
    if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || 
        is.null(df.residual)) 
        stop("No data, or argument is not a valid lmFit object")
    if (all(df.residual == 0)) 
        stop("No residual degrees of freedom in linear model fits")
    if (all(!is.finite(sigma))) 
        stop("No finite residual standard deviations")
    if (trend) {
        covariate <- fit$Amean
        if (is.null(covariate)) 
            stop("Need Amean component in fit to estimate trend")
    }
    else {
        covariate <- NULL
    }
    out <- squeezeVar(sigma^2, df.residual, covariate = covariate, 
        robust = robust, winsor.tail.p = winsor.tail.p)
    out$s2.prior <- out$var.prior
    out$s2.post <- out$var.post
    out$var.prior <- out$var.post <- NULL
    out$t <- coefficients/stdev.unscaled/sqrt(out$s2.post)
    df.total <- df.residual + out$df.prior
    df.pooled <- sum(df.residual, na.rm = TRUE)
    df.total <- pmin(df.total, df.pooled)
    out$df.total <- df.total
    out$p.value <- 2 * pt(-abs(out$t), df = df.total)
    var.prior.lim <- stdev.coef.lim^2/median(out$s2.prior)
    out$var.prior <- tmixture.matrix(out$t, stdev.unscaled, df.total, 
        proportion, var.prior.lim)
    if (any(is.na(out$var.prior))) {
        out$var.prior[is.na(out$var.prior)] <- 1/out$s2.prior
        warning("Estimation of var.prior failed - set to default value")
    }
    r <- rep(1, NROW(out$t)) %o% out$var.prior
    r <- (stdev.unscaled^2 + r)/stdev.unscaled^2
    t2 <- out$t^2
    Infdf <- out$df.prior > 10^6
    if (any(Infdf)) {
        kernel <- t2 * (1 - 1/r)/2
        if (any(!Infdf)) {
            t2.f <- t2[!Infdf]
            r.f <- r[!Infdf]
            df.total.f <- df.total[!Infdf]
            kernel[!Infdf] <- (1 + df.total.f)/2 * log((t2.f + 
                df.total.f)/(t2.f/r.f + df.total.f))
        }
    }
    else kernel <- (1 + df.total)/2 * log((t2 + df.total)/(t2/r + 
        df.total))
    out$lods <- log(proportion/(1 - proportion)) - log(r)/2 + 
        kernel


    out$ts <- fit$ts_sub / sqrt(out$s2.post)
    out$Fs <- fit$fs_sub / out$s2.post
    out
}

# modifies eBayes to calculate F statistics for the pre-specified lists of coefficients.
#   calls modified ebayes_mod which calculates the new values ts and Fs.
#   also calculates p-values for the new Fs
#   the Fs are the goal. The ts should be identical to the t values calculate normally.
eBayes_mod = function (fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
    trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) 
{
    if (trend) 
        if (is.null(fit$Amean)) 
            stop("Need Amean component in fit to estimate trend")
    eb <- ebayes_mod(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, 
        trend = trend, robust = robust, winsor.tail.p = winsor.tail.p)
    fit$df.prior <- eb$df.prior
    fit$s2.prior <- eb$s2.prior
    fit$var.prior <- eb$var.prior
    fit$proportion <- proportion
    fit$s2.post <- eb$s2.post
    fit$t <- eb$t
    fit$df.total <- eb$df.total
    fit$p.value <- eb$p.value
    fit$lods <- eb$lods
    if (!is.null(fit$design) && is.fullrank(fit$design)) {
        F.stat <- classifyTestsF(fit, fstat.only = TRUE)
        fit$F <- as.vector(F.stat)
        df1 <- attr(F.stat, "df1")
        df2 <- attr(F.stat, "df2")
        if (df2[1] > 1e+06) 
            fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)
        else fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
    }
    if(length(fit$test_coef_ranks) > 0) {
        fit$ts = eb$ts
        fit$Fs = eb$Fs
        fit$Fs.p.values = array(NA,dim(fit$Fs))
        for(j in 1:ncol(fit$Fs)){
            fit$Fs.p.values[,j] <- pf(fit$Fs[,j],fit$test_coef_ranks[j],fit$df.prior+fit$df.residual,lower.tail=F)
        }
        colnames(fit$Fs) = colnames(fit$Fs.p.values) = names(fit$test_coef_ranks)
    }
    fit
}


TopTableF2 = function (fit, number = 10, genelist = fit$genes, adjust.method = "BH", 
    sort.by = "F", p.value = 1, lfc = 0) 
{
    if (is.null(fit$coefficients)) 
        stop("Coefficients not found in fit")
    M <- as.matrix(fit$coefficients)
    rn <- rownames(M)
    if (is.null(colnames(M))) 
        colnames(M) <- paste("Coef", 1:ncol(M), sep = "")
    Amean <- fit$Amean
    Fstat <- fit$F
    Fp <- fit$F.p.value
    if (is.null(Fstat)) 
        stop("F-statistics not found in fit")
    if (!is.null(genelist) && is.null(dim(genelist))) 
        genelist <- data.frame(ProbeID = genelist, stringsAsFactors = FALSE)
    if (is.null(rn)) 
        rn <- 1:nrow(M)
    else if (anyDuplicated(rn)) {
        if (is.null(genelist)) 
            genelist <- data.frame(ID = rn, stringsAsFactors = FALSE)
        else if ("ID" %in% names(genelist)) 
            genelist$ID0 <- rn
        else genelist$ID <- rn
        rn <- 1:nrow(M)
    }
    sort.by <- match.arg(sort.by, c("F", "none"))
    adj.P.Value <- p.adjust(Fp, method = adjust.method)
    if (lfc > 0 || p.value < 1) {
        if (lfc > 0) 
            big <- rowSums(abs(M) > lfc, na.rm = TRUE) > 0
        else big <- TRUE
        if (p.value < 1) {
            sig <- adj.P.Value <= p.value
            sig[is.na(sig)] <- FALSE
        }
        else sig <- TRUE
        keep <- big & sig
        if (!all(keep)) {
            M <- M[keep, , drop = FALSE]
            rn <- rn[keep]
            Amean <- Amean[keep]
            Fstat <- Fstat[keep]
            Fp <- Fp[keep]
            genelist <- genelist[keep, , drop = FALSE]
            adj.P.Value <- adj.P.Value[keep]
        }
    }
    if (nrow(M) < number) 
        number <- nrow(M)
    if (number < 1) 
        return(data.frame())
    if (sort.by == "F") 
        o <- order(Fp, decreasing = FALSE)[1:number]
    else o <- 1:number
    if (is.null(genelist)) 
        tab <- data.frame(M[o, , drop = FALSE])
    else tab <- data.frame(genelist[o, , drop = FALSE], M[o, 
        , drop = FALSE])
    tab$AveExpr <- Amean[o]
    tab <- data.frame(tab, F = Fstat[o], P.Value = Fp[o], adj.P.Val = adj.P.Value[o])
    rownames(tab) <- rn[o]
    tab
}