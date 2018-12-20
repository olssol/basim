#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

//' Simon's two-stage rules 
//' @export
// [[Rcpp::export]]
void bacSimonSingle0(NumericMatrix y, int n1, int r1, int n, int r, int bsize, NumericVector rst) {

  // NumericVector rst = NumericVector::create(_["nres"]  = 0.0,
  //                                           _["nenr"]  = 0.0,
  //                                           _["mean"]  = 0.0,
  //                                           _["earl"]  = 0.0,
  //                                           _["rej"]   = 0.0);

  int rep, i, sy, inx;

  //start of 2nd stage
  inx  = bsize * floor(n1/bsize);
  inx += bsize * (n1 % bsize > 0);

  for (rep = 0; rep < y.nrow(); rep++) {
    sy = 0;
    //stage 1
    for (i = 0; i < n1; i++) {
      sy += y(rep,i);
    }

    //stage 2
    if (sy > r1) {
      for (i = inx; i < n; i++) {
        sy += y(rep,i);
      }
    }

    rst(0) += sy;
    rst(1) += i;
    rst(2) += 1.0 * sy / i;
    rst(3) += (i == n1);
    rst(4) += (sy > r);
  }

  //summary
  for (i = 0; i < rst.size(); i++) {
    rst[i] /= y.nrow();
  } 

  //return(rst);
}

void tlSetVisOpt0(int i, int j, NumericMatrix visited, NumericVector optimal,
                  NumericMatrix y0, NumericMatrix y1, int n1, int n, int bsize,
                  double alpha, double beta, NumericVector cp0, NumericVector cp1) {

  if (1 == visited(i,j))
    return;

  // visit current
  visited(i,j) = 1;

  double        rej0, rej1;
  int           l, m;

  bacSimonSingle0(y0, n1, i, n, i+j, bsize, cp0);
  bacSimonSingle0(y1, n1, i, n, i+j, bsize, cp1);
  rej0 = round(cp0["rej"]*1000)/1000;
  rej1 = round(cp1["rej"]*1000)/1000;

  //Rcout << i << ":" << j << "--" << rej0 << "-" << rej1 << "\n";

  if (rej0 > alpha) {
    for (l = 0; l <= i; l++) {
      for (m = 0; m <= j; m++) {
        visited(l,m) = 1;
      }
    }
  }

  if (rej1 < 1-beta) {
    for (l = i; l < n1+1; l++) {
      for (m = j; m < n-n1+1; m++) {
        visited(l,m) = 1;
      }
    }
  }

  if (rej0 <= alpha  &
      rej1 >= 1-beta &
      cp0["nenr"] < optimal(0)) {
    optimal(0) = cp0["nenr"];
    optimal(1) = cp0["earl"];
    optimal(2) = cp0["rej"];
    optimal(3) = cp1["rej"];
    optimal(4) = i;
    optimal(5) = i+j;
  }
}

//' Search r1 and r given n1 and n
//' @export
// [[Rcpp::export]]
void bacSimonSearchR0(NumericMatrix y0, NumericMatrix y1, int n1, int n,
                     int bsize, double alpha, double beta, NumericVector optimal) {

  NumericVector cp0 = NumericVector::create(_["nres"]  = 0.0,
                                            _["nenr"]  = 0.0,
                                            _["mean"]  = 0.0,
                                            _["earl"]  = 0.0,
                                            _["rej"]   = 0.0);
  NumericVector cp1 = clone(cp0); 
  int           i, j, k, n2 = n-n1;

  NumericMatrix visited(n1+1, n2+1);

  //initial visited
  std::fill(visited.begin(), visited.end(), 0);

  for (k = 0; k < fmin(n1+1, n2+1); k = k+2) {
    tlSetVisOpt0(k,    k,    visited, optimal, y0, y1, n1, n, bsize, alpha, beta, cp0, cp1);
    tlSetVisOpt0(n1-k, n2-k, visited, optimal, y0, y1, n1, n, bsize, alpha, beta, cp0, cp1);
    tlSetVisOpt0(n1-k, k,    visited, optimal, y0, y1, n1, n, bsize, alpha, beta, cp0, cp1);
    tlSetVisOpt0(k,    n2-k, visited, optimal, y0, y1, n1, n, bsize, alpha, beta, cp0, cp1);
  }


  //Rf_PrintValue(visited);

  // all combinations
  for (i = n1; i >= 0; i--) {
    for (j = n2; j >= 0; j--) {
      tlSetVisOpt0(i, j, visited, optimal, y0, y1, n1, n, bsize, alpha, beta, cp0, cp1);
    }
  }

  //Rf_PrintValue(optimal);
}


//' Simon's two-stage design 
//' @export
// [[Rcpp::export]]
NumericMatrix bacSimonDesign0(NumericMatrix y0, NumericMatrix y1, int nmax, int nmin,
                             int bsize, double alpha, double beta) {

  NumericMatrix rst(2, 8);
  NumericVector optimal(6, 999999.0);

  //initial 
  std::fill(rst.begin(), rst.end(), 999999.0);

  int n1, n, step = 5;
  for (n = nmin; n <= nmax; n = n+step) {
    Rcout << "n = " << n << "\n";
    for (n1 = 5; n1 < n-1; n1++) {
      bacSimonSearchR0(y0, y1, n1, n, bsize, alpha, beta, optimal);

      //accelerated n
      if (step > 1 & optimal[0] < 999999) {
        n    = n - step;
        step = 1;
      }

      if (optimal(0) < rst(0,4)) {
        rst(0,0) = optimal(4);
        rst(0,1) = n1;
        rst(0,2) = optimal(5);
        rst(0,3) = n;
        rst(0,4) = optimal(0);
        rst(0,5) = optimal(1);
        rst(0,6) = optimal(2);
        rst(0,7) = optimal(3);

        //minimax
        if (n < rst(1, 1)) {
          rst(1,_) = rst(0,_);
        }
      }
    }
  }

  colnames(rst) = CharacterVector::create("r1", "n1", "r", "n",
                                          "en0", "pet0",
                                          "type1", "power");
  return(rst);
}

