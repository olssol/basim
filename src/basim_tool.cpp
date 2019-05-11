#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]


//' Get InterClass Correlation
//' @export
// [[Rcpp::export]]
double bacICC(NumericMatrix ys) {
    int    i, j, k, nc, nr;
    double cor, vx;
    double mu = 0, num = 0, denom = 0;
    double rst;

    nc = ys.ncol();
    nr = ys.nrow();

    // get mu
    mu = 0;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
          mu += ys(i,j);
        }
    }
    mu /= (nr * nc);

    //get icc
    for (i = 0; i < nr; i++) {
        vx  = 0;
        for (j = 0; j < nc; j++) {
          vx += pow(ys(i,j)-mu, 2);
        }
        denom += vx / nc;

        cor = 0;
        for (j = 0; j < nc-1; j++) {
          for (k = j+1; k < nc; k++) {
            cor += (ys(i,j) - mu) * (ys(i,k) - mu);
          }
        }

        num += cor/nc/(nc-1)*2;
    }

    rst = num/denom;
    return(rst);
}


//' Get batch sizes to get the total n
//' @export
// [[Rcpp::export]]
NumericVector baBatches(int n, int bsize) {
  int nb = floor(n/bsize), nl = n % bsize;
  int nt = nb;
  int i;

  if (nl > 0) {
    nt++;
  }

  NumericVector rst(nt);
  for (i = 0; i < nb; i++) {
    rst(i) = bsize;
  }

  if (nl > 0)
    rst(nt - 1) = nl;

  return(rst);
}

//' Get variance given batch sizes, response rate and ICC
//'  @export
// [[Rcpp::export]]
double baSaGetVar(NumericVector bsizes, double p, double rho) {
  int  nb  = bsizes.size();
  double rst = 0, pp = p * (1-p);
  int  i;

  for (i = 0; i < nb; i++) {
    rst += bsizes(i) * pp * (1 + (bsizes(i) - 1) * rho);
  }

  return(rst);
}


//' Get actuarial type I and II error in Sargent method
//' @export
// [[Rcpp::export]]
NumericVector baSaAlphaBeta(NumericVector bsizes, int r, double p0, double p1, double rho0, double rho1) {
  double v0, v1, m0, m1, t0, t1;
  int  nb, n = 0, i;

  nb = bsizes.size();
  for (i = 0; i < nb; i++) {
    n += bsizes(i);
  }

  v0 = baSaGetVar(bsizes, p0, rho0);
  v1 = baSaGetVar(bsizes, p1, rho1);
  m0 = n * p0;
  m1 = n * p1;
  t0 = (r + 1.5 - m0)/sqrt(v0);
  t1 = (r + 1.5 - m1)/sqrt(v1);

  NumericVector rst = NumericVector::create(_["alpha"] = 1 - R::pnorm(t0, 0, 1, 1, 0),
                                            _["beta"]  = R::pnorm(t1, 0, 1, 1, 0));

  return(rst);
}

//' Get Design type I and II error using Sargent method
//' @export
// [[Rcpp::export]]
NumericVector baSaDesign(int nmin, int bsize, double alpha, double beta, double p0, double p1,
                         double rho0, double rho1) {
  int flag = 0, n = nmin - 1, r;
  NumericVector ab(2);
  NumericVector bsizes;

  while (0 == flag) {
    n++;
    bsizes = baBatches(n, bsize);
    for (r = floor(n*p0); r <= n; r++) {
      ab = baSaAlphaBeta(bsizes, r, p0, p1, rho0, rho1);

      if (ab(0) < alpha & ab(1) < beta) {
        flag = 1;
        break;
      }
    }
  }

  ab = baSaAlphaBeta(bsizes, r, p0, p1, 0, 0);
  return(ab);
}
