#' Get designs for BetaBinomial designs
#'
#'
#'
#' @export
#'
baDesignBetaBin <- function(par.design, rho = 0, bsizes = 1, nmin = 10, nmax = 100, ...) {
    rho <- rep(rho, 2);
    yp0 <- baSimuBetaBin(max(bsizes), p = par.design$P0, rho = rho[1], ...)$y;
    yp1 <- baSimuBetaBin(max(bsizes), p = par.design$P1, rho = rho[2], ...)$y;

    rst <- NULL;
    for (i in 1:length(bsizes)) {
        cur.rst        <- bacSimonDesign(yp0, yp1, nmax, nmin, bsizes[i],
                                        par.design$ALPHA, 1 - par.design$POWER);
        cur.rst        <- data.frame(cur.rst);
        cur.rst$design <- c("Optimal", "Minimax");
        cur.rst$bsize  <- bsizes[i];
        rst            <- rbind(rst, cur.rst);
    }

    rst$rho0  <- rho[1];
    rst$rho1  <- rho[2];
    rst$p0    <- par.design$P0;
    rst$p1    <- par.design$P1;
    rst$alpha <- par.design$ALPHA;
    rst$beta  <- par.design$BETA;
    rst
}

#' Get type I error and power given y0, y1, n1, r1, n, r
#'
#'
#'
#' @export
#'
baGetRejRate <- function(yp0, yp1, n1, r1, nt, r, bsize = 1) {

    stopifnot(ncol(yp0) >= bsize &
              ncol(yp1) >= bsize);

    cumu0 <- bacCumProb(yp0, nt, n1, bsize);
    cumu1 <- bacCumProb(yp1, nt, n1, bsize);
    cp0   <- bacSimonSingle(cumu0, n1, r1, nt, r);
    cp1   <- bacSimonSingle(cumu1, n1, r1, nt, r);

    c(type1 = cp0[3], power = cp1[3],
      en0   = cp0[1], en1   = cp1[1],
      pet0  = cp0[2], pet1  = cp1[2]);
}


## #' Get actuarial type I and II error in Sargent's method
## baSaAlphaBeta <- function(bsizes, r, p0, p1, rho = 0) {
##     v0    <- baSaGetVar(bsizes, p0, rho);
##     v1    <- baSaGetVar(bsizes, p1, rho);
##     n     <- sum(bsizes);
##     m0    <- n * p0;
##     m1    <- n * p1;

##     alpha <- 1 - pnorm( (r + 1.5 - m0)/sqrt(v0) );
##     beta  <- pnorm( (r + 1.5 - m1)/sqrt(v1) );

##     c(alpha, beta);
## }


