#' Simulate batch data from Beta-Binomial
#'
#' @param batch.sizes number of subjects in each batch
#' @param p mean of the Beta distribution a/(a+b)
#' @param rho inter-class correlation 1/(a+b+1)
#'
#' @return Simulated binary outcomes by batches
#' @export
#'
baSimuBetaBin <- function(batch.sizes, p, rho = 0, seed = NULL) {

    if (!is.null(seed))
        set.seed(seed);

    nbatch <- length(batch.sizes);

    if (0 == rho) {
        bps <- rep(p, nbatch);
    } else {
        alpha  <- p * (1/rho - 1);
        beta   <- (1-p) * (1/rho - 1);
        bps    <- rbeta(nbatch, alpha, beta);
    }

    rst <- apply(cbind(batch.sizes, bps), 1,
                 function(x) rbinom(x[1], size = 1, prob = x[2]));

    ## return
    list(y   = as.numeric(rst),
         ind = baBatInd(batch.sizes),
         ps  = bps);
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
baSimuBatchError <- function(batch.sizes,
                             par.err = list(gamma   = list(error.type = "normal", sig = 1),
                                            delta   = list(error.type = "normal", sig = 1),
                                            epsilon = list(error.type = "normal", sig = 1))) {

    n.batch <- length(batch.sizes);
    n.tot   <- sum(batch.sizes);

    ## batch effects
    beffs <- NULL;
    for (ef in c("gamma", "delta")) {
        cur.par <- par.err[[ef]];
        cur.eff <- do.call(baSimuError, c(n = n.batch, cur.par));
        beffs   <- cbind(beffs, cur.eff);
    }

    epsilon <- do.call(baSimuError, c(n = n.tot,
                                      par.err[["epsilon"]]));

    rst <- list(gamma   = beffs[,1],
                delta   = beffs[,2],
                epsilon = epsilon,
                bis     = baBatInd(batch.sizes));

    class(rst) <- "ClsBaErr";
    rst;
}



#' Simulate t-cell counts using Negative binomial
#'
#' @export
#'
baSimuTcell <- function(batch.sizes, par.err, par.other, ...) {
    ntot <- sum(batch.sizes);

    ## par.others
    u0   <- par.other["u0"];
    u1   <- par.other["u1"];
    v    <- par.other["v"];
    beta <- par.other["beta"];

    ## baseline error and oucome
    err.base <- baSimuBatchError(batch.sizes, par.err = par.err);
    u0.eps   <- u0 + err.base$epsilon;
    log.mu0  <- u0.eps + err.base$gamma[err.base$bis];
    log.phi0 <- v + err.base$delta[err.base$bis];
    y0       <- baRNB(log.mu0, log.phi0);

    ## post treatment error and oucome
    err.post <- baSimuBatchError(batch.sizes, par.err = par.err);
    log.mu1  <- u1 + beta*u0.eps + err.post$gamma[err.post$bis] + err.post$epsilon;
    log.phi1 <- v + err.post$delta[err.post$bis];
    y1       <- baRNB(log.mu1, log.phi1);
    ry       <- baGetOutcome(cbind(y0, y1), ...);

    ## return
    rst <- list(y0 = y0, y1 = y1, y = ry);
}

#' Get cut off of the outcome to get given response rates
#'
#'
#' @export
#'
baGetCuts <- function(par.err, par.other, rates, n.reps = 100000,
                      f.simu = baSimuTcell, ..., seed = NULL) {

    if (!is.null(seed))
        set.seed(seed);

    true.pts <- f.simu(rep(1, n.reps),
                       par.err = par.err,
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
baGetBvrCvIcc <- function(par.err, par.other, rates, n.reps = 100000, f.simu = baSimuTcell, ...) {
    par.e2         <- par.err;
    par.e2$gamma   <- list(error.type = "normal", sig = 0);
    par.e2$delta   <- list(error.type = "normal", sig = 0);

    cur.y.wb <- f.simu(rep(1, n.reps), par.err, par.other, ...);
    cur.y.wo <- f.simu(rep(1, n.reps), par.e2,  par.other, ...);

    v.wb     <- var(cur.y.wb$y0, na.rm = TRUE);
    v.wo     <- var(cur.y.wo$y0, na.rm = TRUE);
    bvr      <- v.wb / v.wo;
    cv       <- 100 * sqrt(v.wb - v.wo)/mean(cur.y.wo$y0);

    ## icc depends on the binary outcomes
    cuty     <- quantile(cur.y.wb$y, probs = 1 - rates);
    y.wb.ba  <- f.simu(rep(10, n.reps), par.err, par.other, ...)$y;

    icc <- NULL;
    for (i in 1:length(cuty)) {
        icc <- c(icc, bacICC(y.wb.ba > cuty[i], rep(10, n.reps)));
    }

    list(par.err   = par.err,
         par.other = par.other,
         bvr       = bvr,
         cv        = unname(cv),
         icc       = cbind(rates = rates,
                           cuts  = unname(cuty),
                           icc   = icc));
}
