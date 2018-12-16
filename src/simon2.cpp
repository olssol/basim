#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

//' Simon's two-stage rules 
//' @export
// [[Rcpp::export]]
NumericVector bacSimonSingle(NumericMatrix y, int n1, int r1, int n, int r, int bsize) {

  NumericVector rst = NumericVector::create(_["nres"]  = 0.0,
                                            _["nenr"]  = 0.0,
                                            _["mean"]  = 0.0,
                                            _["earl"]  = 0.0,
                                            _["rej"]   = 0.0);

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

    rst[0] += sy;
    rst[1] += i;
    rst[2] += 1.0 * sy / i;
    rst[3] += (i == n1);
    rst[4] += (sy > r);
  }

  //summary
  for (i = 0; i < rst.size(); i++) {
    rst[i] /= y.nrow();
  } 

  return(rst);
}

NumericMatrix tlSetVis(NumericMatrix visited, int i, int j, int nr, int nc,
                       double rej0, double rej1, double alpha, double beta) {

  NumericMatrix rst(visited);
  int l, m;

  //current visited
  rst(i,j) = 1;

  if (round(rej0*1000)/1000 > alpha) {
    for (l = 0; l <= i; l++) {
      for (m = 0; m <= j; m++) {
        rst(l,m) = 1;
      }
    }
  }

  if (round(rej1*1000)/1000 < 1-beta) {
    for (l = i; l < nr; l++) {
      for (m = j; m < nc; m++) {
        rst(l,m) = 1;
      }
    }
  }

  return(rst);
}

NumericVector tlSetOpt(NumericVector optimal, NumericVector cp0, NumericVector cp1,
                       int i, int j, double alpha, double beta) {

  NumericVector rst(optimal);

  if (round(cp0["rej"]*1000)/1000 <= alpha  &
      round(cp1["rej"]*1000)/1000 >= 1-beta &
      cp0["nenr"] < optimal(0)) {

    rst(0) = cp0["nenr"];
    rst(1) = cp0["earl"];
    rst(2) = cp0["rej"];
    rst(3) = cp1["rej"];
    rst(4) = i;
    rst(5) = i+j;
  }

  return(rst);
}

void tlSetVisOpt(int i, int j, NumericMatrix visited, NumericVector optimal,
                 NumericMatrix y0, NumericMatrix y1, int n1, int n, int bsize,
                 double alpha, double beta) {

  if (1 == visited(i,j))
    return;

  // visit current
  visited(i,j) = 1;

  NumericVector cp0, cp1;
  double        rej0, rej1;
  int           l, m;

  cp0  = bacSimonSingle(y0, n1, i, n, i+j, bsize);
  cp1  = bacSimonSingle(y1, n1, i, n, i+j, bsize);
  rej0 = round(cp0["rej"]*1000)/1000;
  rej1 = round(cp1["rej"]*1000)/1000;

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
NumericVector bacSimonSearchR(NumericMatrix y0, NumericMatrix y1, int n1, int n,
                              int bsize, double alpha, double beta) {

  NumericVector optimal(6, 999999.0);
  NumericVector cp0(5), cp1(5);
  int           i, j, k, n2 = n-n1;

  NumericMatrix visited(n1+1, n2+1);

  //initial visited
  std::fill(visited.begin(), visited.end(), 0);

  // diagonal only 
  // Rf_PrintValue(diagij);

  for (k = 0; k < fmin(n1+1, n2+1); k++) {
    tlSetVisOpt(k,    k,    visited, optimal, y0, y1, n1, n, bsize, alpha, beta);
    tlSetVisOpt(n1-k, n2-k, visited, optimal, y0, y1, n1, n, bsize, alpha, beta);
    tlSetVisOpt(n1-k, k,    visited, optimal, y0, y1, n1, n, bsize, alpha, beta);
    tlSetVisOpt(k,    n2-k, visited, optimal, y0, y1, n1, n, bsize, alpha, beta);
  }

  // all combinations
  for (i = 0; i < n1+1; i++) {
    for (j = 0; j < n2+1; j++) {
      tlSetVisOpt(i, j, visited, optimal, y0, y1, n1, n, bsize, alpha, beta);
    }
  }

  return(optimal);
}


//' Simon's two-stage design 
//' @export
// [[Rcpp::export]]
NumericMatrix bacSimonDesign(NumericMatrix y0, NumericMatrix y1, int nmax, int nmin,
                             int bsize, double alpha, double beta) {

  NumericMatrix rst(2, 8);
  NumericVector optimal;

  //initial 
  std::fill(rst.begin(), rst.end(), 999999.0);

  int n1, n;
  for (n = nmin; n <= nmax; n++) {
    Rcout << "n = " << n << "\n";
    for (n1 = 5; n1 < n-1; n1++) {
      optimal = bacSimonSearchR(y0, y1, n1, n, bsize, alpha, beta);

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



//' Get InterClass Correlation
//' @export
// [[Rcpp::export]]
double bacICC(NumericVector ys, NumericVector bsizes) {
    int    i, j, k;
    double sumx, cor, vx;
    double mu = 0, num = 0, denom = 0;
    int    nbatch = bsizes.size(), cur;
    double rst;

    // get mu
    cur = 0;
    for (i = 0; i < nbatch; i++) {
        sumx = 0;
        for (j = 0; j < bsizes(i); j++) {
            sumx += ys[cur + j];
        }

        mu  += sumx / bsizes(i);
        cur += bsizes(i);
    }
    mu /= nbatch;

    //get icc
    cur = 0;
    for (i = 0; i < nbatch; i++) {
        vx  = 0;
        for (j = 0; j < bsizes(i); j++) {
            vx += pow(ys[cur+j]-mu, 2);
        }
        denom += vx / bsizes(i);

        cor = 0;
        for (j = 0; j < bsizes(i)-1; j++) {
            for (k = j+1; k < bsizes(i); k++) {
                cor += (ys[j+cur] - mu) * (ys[k+cur] - mu);
            }
        }

        num += cor/bsizes(i)/(bsizes(i)-1)*2;
        cur += bsizes(i);
    }

    rst = num/denom;
    return(rst);
}
