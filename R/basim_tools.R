
#' Get batch indices
#'
#' @export
#'
baBatInd <- function(batch.sizes) {
    nbatch <- length(batch.sizes);
    bis    <- apply(cbind(1:nbatch, batch.sizes),
                    1,
                    function(x) rep(x[1], x[2]));
    bis    <- as.numeric(unlist(bis));
}


#' Define outcome based on longitudinal outcomes
#'
#'
#' @export
#'
baGetOutcome <- function(mat.y, type = c("ratio", "change", "base")) {
    type <- match.arg(type);
    switch(type,
           "ratio"  = mat.y[,ncol(mat.y)] / mat.y[,1],
           "change" = mat.y[,ncol(mat.y)] / mat.y[,1],
           "base"   = mat.y[,1]);
}


#' Simulate NB using mu and phi instead of n and phi
#'
#'
#' @export
#'
baRNB <- function(mu, phi, n = NULL, is.log = TRUE) {

    if (is.null(n))
        n <- length(mu);

    if (is.log) {
        mu  <- exp(mu);
        phi <- exp(phi);
    }

    rnbinom(n, size = phi, prob = phi/(mu + phi));
}



