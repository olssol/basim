#' Get designs
#'
#'
#'
#' @export
#'
baGetDesign  <- function(par.design, lst.setting = NULL, design = c("simon", "truth", "betabinom"),
                         f.simu = baSimuTcell, rho = NULL, adj.rho = 1, bsize = 1,
                         nmin = 10, nmax = 100, nreps = 10000, ...) {

    f.dta <- function(rst, ...) {
        data.frame(r1     = rst[1],
                   n1     = rst[2],
                   r      = rst[3],
                   n      = rst[4],
                   en     = rst[5],
                   pet    = rst[6],
                   ...);
    }

    design <- match.arg(design);
    if ("simon" == design) {
        rst <- ph2simon(pu  = par.design$P0,
                        pa  = par.design$P1,
                        ep1 = par.design$ALPHA,
                        ep2 = 1 - par.design$POWER);
        rst <- rst$out[which.min(rst$out[,5]),];
        rst <- f.dta(rst,
                     design = design,
                     rhoadj = NA,
                     rho0   = NA,
                     rho1   = NA,
                     bsize  = NA);

    } else {
        cur.icc <- lst.setting$icc;
        bsizes  <- rep(bsize, (ceiling(nmax/bsize)+2)*bsize);

        if ("betabinom" == design) {
            if (is.null(rho)) {
                rho0 <- adj.rho * cur.icc[which(par.design$P0 == cur.icc[,"rates"]), "icc"];
                rho1 <- adj.rho * cur.icc[which(par.design$P1 == cur.icc[,"rates"]), "icc"];
            } else {
                rho0 <- rho[1];
                rho1 <- rho[min(2, length(rho))];
            }

            yp0 <- sapply(1:nreps, function(x) baSimuBetaBin(bsizes, p = par.design$P0, rho = rho0)$y);
            yp1 <- sapply(1:nreps, function(x) baSimuBetaBin(bsizes, p = par.design$P1, rho = rho1)$y);
            yp0 <- t(yp0);
            yp1 <- t(yp1);

            rst <- bacSimonDesign(yp0, yp1, nmax, nmin, bsize, par.design$ALPHA, 1-par.design$POWER);

            rst <- f.dta(rst[1,],
                         design = design,
                         rhoadj = adj.rho,
                         rho0   = rho0,
                         rho1   = rho1,
                         bsize  = bsize);

        } else if ("truth" == design) {
            cut0 <- cur.icc[which(par.design$P0 == cur.icc[,"rates"]), "cuts"];
            cut1 <- cur.icc[which(par.design$P1 == cur.icc[,"rates"]), "cuts"];
            ys   <- sapply(1:nreps,
                           function(x) f.simu(bsizes,
                                              par.err = lst.setting$par.err,
                                              par.other = lst.setting$par.other, ...)$y);

            yp0 <- t(ys > cut0);
            yp1 <- t(ys > cut1);
            rst <- bacSimonDesign(yp0, yp1, nmax, nmin, bsize, par.design$ALPHA, 1-par.design$POWER);

            rst <- f.dta(rst[1,],
                         design = design,
                         rhoadj = NA,
                         rho0   = NA,
                         rho1   = NA,
                         bsize  = bsize);
        }
    }

    rst <- cbind(rst, par.design);
    row.names(rst) <- NULL;
    rst
}
