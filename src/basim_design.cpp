#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

NumericVector baBatches(int n, int bsize);


//' PDF of single batch
//' @export
// [[Rcpp::export]]
NumericMatrix bacBatchFreq(NumericMatrix y) {
  int           nrep  = y.nrow();
  int           bsize = y.ncol();
  NumericMatrix rst(bsize, bsize+1);

  //initial
  std::fill(rst.begin(), rst.end(), 0);

  int i, j, r;
  NumericMatrix::iterator ir;

  // counts
  for (i = 0; i < nrep; i++) {
    r = 0;
    for (j = 0; j < bsize; j++) {
      r += y(i,j);
      rst(j, r)++;
    }
  }

  for (ir = rst.begin(); ir != rst.end(); ++ir) {
    *ir /= nrep;
  }

  return(rst);
}


//' Probability of (n = n, r = 0:n)
//' @export
// [[Rcpp::export]]
NumericMatrix bacProb(int n, int bsize, NumericMatrix pmat) {
  NumericVector bs = baBatches(n, bsize);
  NumericVector tmp(n+1);
  NumericMatrix rst(n+1,2);

  int i, j, k, r;

  //initial
  std::fill(rst.begin(), rst.end(), 0);

  //PDF
  for (i = 0; i <= bs(0); i++) {
    rst(i,0) = pmat(bs(0) - 1, i);
  }
  r = bs(0);

  //continue
  for (i = 1; i < bs.size(); i++) {
    std::fill(tmp.begin(), tmp.end(), 0);
    for (j = 0; j <= r; j++) {
      for (k = 0; k <= bs(i); k++) {
        tmp(j+k) += pmat(bs(i) - 1, k) * rst(j,0);
      }
    }
    r += bs(i);
    for (j = 0; j <= r; j++) {
      rst(j,0) = tmp(j);
    }
  }

  //CDF
  rst(0,1) =  rst(0,0);
  for (i = 1; i < n+1; i++) {
    rst(i,1) = rst(i,0) + rst(i-1,1);
  }

  return(rst);
}

void tkProb(NumericMatrix y, int nmax, int nmin,  int bsize,
            NumericMatrix prob, NumericMatrix cumu) {
  NumericMatrix pmat = bacBatchFreq(y);
  NumericMatrix tmp;
  int i,j;

  //init
  std::fill(prob.begin(), prob.end(), 0);
  std::fill(cumu.begin(), cumu.end(), 0);

  for (i = 1; i <= nmax; i++) {
    tmp = bacProb(i, bsize, pmat);

    for (j = 0; j <= i; j++) {
      prob(i - 1, j) = tmp(j,0);
      cumu(i - 1, j) = tmp(j,1);
    }
  } 
}

//' Prepare probabilities
//' @return A matrix with 2*nmax rows. Row 1-nmax is marginal P(x = r|n)
//'         Row nmax-2*nmax P(x<=r|n);
//' @export
// [[Rcpp::export]]
NumericMatrix bacCumProb(NumericMatrix y, int nmax, int nmin,  int bsize) {
  NumericMatrix cumu(2*nmax, nmax+1);
  NumericMatrix pmat = bacBatchFreq(y);
  NumericMatrix tmp;
  int i,j;

  //init
  std::fill(cumu.begin(), cumu.end(), 0);

  for (i = 1; i <= nmax; i++) {
    tmp = bacProb(i, bsize, pmat);
    for (j = 0; j <= i; j++) {
      cumu(i - 1, j)        = tmp(j,0);
      cumu(nmax+(i - 1), j) = tmp(j,1);
    }
  }

  return(cumu);
}

//' Single Simon 2-stage
//' @export
// [[Rcpp::export]]
NumericVector bacSimonSingle(NumericMatrix cumu, int n1, int r1, int n, int r) {
  NumericVector rst(3);
  int    nmax = cumu.nrow() / 2;
  int    n2   = n - n1;
  double stop1, stop, en, x;

  stop1  = cumu(nmax + n1 - 1, r1);
  en     = stop1 * n1 + (1-stop1) * n;  
  stop   = stop1;
  for (x = r1+1; x <= fmin(r, n1); x++) {
    stop += cumu(n1 - 1, x) * cumu(nmax + n2 - 1, r - x);
  }

  rst(0) = en;
  rst(1) = stop1;
  rst(2) = 1 - stop;

  return(rst);
}

//' Simon's two-stage design 
//' @export
// [[Rcpp::export]]
NumericMatrix bacSimonDesign(NumericMatrix y0, NumericMatrix y1, int nmax, int nmin,
                             int bsize, double alpha, double beta) {

  NumericMatrix rst(2, 8);
  NumericVector cp0, cp1;
  int           n1, n, r1, r, step = 5;
  NumericMatrix cumu0, cumu1;

  //init
  std::fill(rst.begin(), rst.end(), 99999);

  cumu0 = bacCumProb(y0, nmax, nmin, bsize);
  cumu1 = bacCumProb(y1, nmax, nmin, bsize);

  for (n = nmin; n <= nmax; n=n+step) {
    for (n1 = 5; n1 < n-1; n1++) {
      for (r = n; r >= 0; r--) {
        for (r1 = fmin(r,n1); r1 >= 0; r1--) {
          cp0 = bacSimonSingle(cumu0, n1, r1, n, r);
          if (cp0(2) > alpha |
              cp0(0) > rst(0,4))
            break;

          cp1 = bacSimonSingle(cumu1, n1, r1, n, r);
          if (ceil(cp1(2)*1000)/1000 < 1 - beta)
            continue;

          //accelerated n
          if (step > 1) {
            n    = n - step;
            step = 1;
          }

          rst(0, 0) = r1;
          rst(0, 1) = n1;
          rst(0, 2) = r;
          rst(0, 3) = n;
          rst(0, 4) = cp0(0);
          rst(0, 5) = cp0(1);
          rst(0, 6) = cp0(2);
          rst(0, 7) = cp1(2);

          if (n <= rst(1, 3)) {
            rst(1,_) = rst(0,_);
          }
        }
      }
    }
  }

  colnames(rst) = CharacterVector::create("r1", "n1", "r", "n",
                                          "en0", "pet0",
                                          "type1", "power");
  return(rst);
}
