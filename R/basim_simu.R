#' Simulate batch data from Beta-Binomial
#'
#' @param bsize batch size
#' @param p mean of the Beta distribution a/(a+b)
#' @param rho inter-class correlation 1/(a+b+1)
#'
#' @return Simulated binary outcomes by batches
#' @export
#'
baSimuBetaBin <- function(bsize, p, rho = 0, seed = NULL, nreps = 1000000) {
    if (!is.null(seed))
        set.seed(seed);

    if (0 == rho) {
        bps <- rep(p, nreps);
    } else {
        alpha  <- p * (1/rho - 1);
        beta   <- (1-p) * (1/rho - 1);
        bps    <- rbeta(nreps, alpha, beta);
    }

    rst      <- rbinom(bsize * nreps, size = 1, prob = bps);
    dim(rst) <- c(nreps, bsize);

    ## return
    list(y = rst, ps = bps);
}


#' Simulate random error
#'
#'
#' @details Skewed distribution: Prameterization: mu, phi; RNBINOM uses (n,p)
#'     with: phi = n, mu = n(1-p)/p; Mean: mu = n(1-p)/p #; Variances:
#'     mu(1+mu/phi) = n (1-p)/p^2
#'
#' @export
#'
baSimuError <- function(n,
                        error.type = c("normal", "skewed"),
                        sig = 1, mu = 0,
                        skew.mu    = NULL,
                        skew.phi   = NULL,
                        skew.noise = 0.001) {

    type <- match.arg(error.type);
    rst <- switch(type,
                  normal = {rnorm(n, mu, sig)},
                  skewed = {
                      ## smu <- skew.n * (1-skew.p) / skew.p;
                      skew.p    <- skew.phi / (skew.mu + skew.phi)
                      skew.va   <- skew.phi * (1-skew.p) / skew.p^2;
                      rst       <- rnbinom(n = n, size = skew.phi, prob = skew.p);
                      rst       <- rst - skew.mu + rnorm(n, mu, skew.noise);
                      rst       <- rst/sqrt(skew.va + skew.noise^2)*sig;
                  });

    rst
}


#' Simulate errors for all batches
#'
#'
#'
#' @export
#'
baSimuBatchError <- function(bsize,
                             par.err = list(gamma   = list(error.type = "normal", sig = 1),
                                            delta   = list(error.type = "normal", sig = 1),
                                            epsilon = list(error.type = "normal", sig = 1)),
                             nreps = 1000000) {

    n.tot <- bsize * nreps;

    ## batch effects
    beffs <- NULL;
    for (ef in c("gamma", "delta")) {
        cur.par <- par.err[[ef]];
        cur.eff <- do.call(baSimuError, c(n = nreps, cur.par));
        beffs   <- cbind(beffs, cur.eff);
    }

    epsilon <- do.call(baSimuError, c(n = n.tot,
                                      par.err[["epsilon"]]));

    rst <- list(gamma   = beffs[,1],
                delta   = beffs[,2],
                epsilon = epsilon);

    class(rst) <- "ClsBaErr";
    rst;
}



#' Simulate t-cell counts using Negative binomial
#'
#' @export
#'
baSimuTcell <- function(bsize, par.err, par.other, nreps = 1000000, ...) {
    ## par.others
    u0   <- par.other["u0"];
    u1   <- par.other["u1"];
    v    <- par.other["v"];
    beta <- par.other["beta"];

    ## baseline error and oucome
    err.base <- baSimuBatchError(bsize = bsize, par.err = par.err, nreps = nreps);
    u0.eps   <- u0     + err.base$epsilon;
    log.mu0  <- u0.eps + rep(err.base$gamma, each = bsize);
    log.phi0 <- v      + rep(err.base$delta, each = bsize);
    y0       <- baRNB(log.mu0, log.phi0);

    ## post treatment error and oucome
    err.post <- baSimuBatchError(bsize = bsize, par.err = par.err, nreps = nreps);
    log.mu1  <- u1 + beta*u0.eps + rep(err.post$gamma, each = bsize) + err.post$epsilon;
    log.phi1 <- v  + rep(err.post$delta, each = bsize);
    y1       <- baRNB(log.mu1, log.phi1);
    ry       <- baGetOutcome(cbind(y0, y1), ...);

    ## return
    rst <- list(y0 = matrix(y0, nrow = nreps, ncol = bsize, byrow = TRUE),
                y1 = matrix(y1, nrow = nreps, ncol = bsize, byrow = TRUE),
                y  = matrix(ry, nrow = nreps, ncol = bsize, byrow = TRUE));
}

#' Get cut off of the outcome to get given response rates
#'
#'
#' @export
#'
baGetCuts <- function(par.err, par.other, rates, f.simu = baSimuTcell, ..., seed = NULL) {

    if (!is.null(seed))
        set.seed(seed);

    true.pts <- f.simu(bsize     = 1,
                       par.err   = par.err,
                       par.other = par.other,
                       ...)$y;

    cut.y <- unname(quantile(true.pts, probs = 1 - rates));
    rst   <- cbind(rate = cut.ps, cuts = cut.y);
    return(rst);
}


#' Get CV, batch effect variance ratio and ICC
#'
#'
#'
#' @export
#'
#'
baGetBvrCvIcc <- function(par.err, par.other, rates, f.simu = baSimuTcell, ...) {
    par.e2         <- par.err;
    par.e2$gamma   <- list(error.type = "normal", sig = 0);
    par.e2$delta   <- list(error.type = "normal", sig = 0);

    cur.y.wb <- f.simu(bsize = 1, par.err, par.other, ...);
    cur.y.wo <- f.simu(bsize = 1, par.e2,  par.other, ...);

    v.wb     <- var(cur.y.wb$y0, na.rm = TRUE);
    v.wo     <- var(cur.y.wo$y0, na.rm = TRUE);
    bvr      <- v.wb / v.wo;
    cv       <- 100 * sqrt(v.wb - v.wo)/mean(cur.y.wo$y0);

    ## icc depends on the binary outcomes
    cuty     <- quantile(cur.y.wb$y, probs = 1 - rates);
    y.wb.ba  <- f.simu(bsize = 10, par.err, par.other, ...)$y;

    icc <- NULL;
    for (i in 1:length(cuty)) {
        icc <- c(icc, bacICC(y.wb.ba > cuty[i]));
    }

    list(par.err   = par.err,
         par.other = par.other,
         bvr       = bvr,
         cv        = unname(cv),
         icc       = cbind(rates = rates,
                           cuts  = unname(cuty),
                           icc   = icc));
}
