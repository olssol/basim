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
