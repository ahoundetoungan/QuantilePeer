/* GENERAL NOTATIONS
 * y       : is the vector of the outcome values
 * X       : is the matrix of the explanatory variables. Add the intercept 
 *           if it is included in the model. The intercept will not be added automatically. 
 * G       : is the network matrix List. That is G[r] is the subnetwork of the group r. 
 *           Gs(i,j) = measures the intensity of the outgoing link from i to j. 
 * igroup  : is the matrix of groups indexes. The data should be organized by group. 
 *           igroup[s,] is a 2-dimension vector ans gives the first and the last rows
 *           for the group s.
 * tau     : is a quantile level
 * stau    : set of values for tau
 * sm      : set of values for m; for type 7 quantile m = 1 - p
 * ntau    : the number of tau
 * ngroup  : is the number of groups.
 * theta   : is the vector of parameters ordered as follow: lambda_tau and explanatory variables (intercept is included)
 * n       : The sample size.
 * nvec    : Vector of agents in each subnet
 * tol     : A tolerance value for the iterative method solving the game. 
 * maxit   : The maximum number of iterations of the iterative method. If this
 *           number is reached, the algorithm stops and the last y is used as the solution. maxit
 *           is important for numerical reasons if tol is too small. 
 */

// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]

#include <RcppArmadillo.h>
//#define NDEBUG
#include <RcppNumerical.h>
#include <RcppEigen.h>

typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
typedef Eigen::Map<Eigen::VectorXd> MapVect;

using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// gP vector of gij for gij>0, this is sorted from the smallest based of yj of peers
// IyP includes the indexes of sorterd yj of peers in the original y
// cumsumgP cumulative sum of gP
// This function computes liste of gP, IyP, cumsumgP
//[[Rcpp::export]]
Rcpp::List fgPIyP(const arma::vec& y,
                  Rcpp::List& G,
                  const arma::vec d,
                  const arma::mat& igroup,
                  const arma::vec& nvec,
                  const int& ngroup,
                  const int& n){
  Rcpp::List lgP(n), lcumsumgP(n), lIyP(n);
  int k(-1), n1, n2;
  arma::vec yr, csgPi, Gri;
  arma::uvec ipGri, tp;
  for(int r(0); r < ngroup; ++ r){
    n1           = igroup(r, 0);
    n2           = igroup(r, 1);
    arma::mat Gr = G[r];
    yr           = y.subvec(n1, n2);
    
    for(int i(0); i < nvec(r); ++ i){
      ++ k;
      if(d(k) == 0) continue;
      Gri = Gr.row(i).t();
      ipGri = arma::find(Gri > 0);
      Gri   = Gri.elem(ipGri);
      Gri   = Gri/d(k);
      tp    = arma::sort_index(yr.elem(ipGri));
      Gri   = Gri.elem(tp);
      tp    = arma::join_cols(ipGri.elem(tp) + n1, arma::zeros<arma::uvec>(1));
      csgPi = arma::cumsum(Gri);
      Gri          = arma::join_cols(Gri, arma::ones<arma::vec>(1));
      lgP[k]       = Gri;
      lcumsumgP[k] = csgPi;
      lIyP[k]      = tp;
    }
  }
  return Rcpp::List::create(_["gP"] = lgP, _["cumsumgP"] = lcumsumgP, _["IyP"] = lIyP);
}

// This function computes the weights and the indexes to dertermine quantiles
// weight for y_(pi + 1) in w2 (matrix n*ntau) 
// index for y_(pi) in pi1 (matrix n*ntau)
// index for y_(pi + 1) in pi2 (matrix n*ntau) 
//[[Rcpp::export]]
void fQWeightIndex(arma::mat& w2, 
                   arma::umat& pi1,
                   arma::umat& pi2,
                   Rcpp::List& lgP, 
                   Rcpp::List& lcumsumgP, 
                   Rcpp::List& lIyP,
                   const arma::vec& d,
                   const arma::vec& stau,
                   const arma::vec& sm,
                   const int& n,
                   const int& ntau){
  for(int i(0); i < n; ++ i){
    if(d(i) == 0) continue;
    arma::vec gP = lgP[i];
    arma::vec cumsumgP = lcumsumgP[i];
    arma::uvec IyP = lIyP[i];
    arma::vec w2i(ntau);
    arma::uvec pii1(ntau);
    
    for(int k(0); k < ntau; ++ k){
      int l(sum(cumsumgP <= stau(k)));
      double tp1(0); if(l > 0) tp1 = cumsumgP(l - 1);
      double tp2(l - 1 + (stau(k) - tp1)/gP(l) + sm(k));
      pii1(k) = floor(tp2);
      w2i(k)   = tp2 - pii1(k);
    }
    w2.col(i)  = w2i;
    pi1.col(i) = IyP.elem(pii1);
    pi2.col(i) = IyP.elem(pii1 + 1);
  }
}

// This function computes Qtau(V) and Qtau y
//[[Rcpp::export]]
Rcpp::List fQtauyVars(const arma::vec& y,
                      const arma::mat Vars,
                      Rcpp::List& G,
                      const arma::vec d,
                      const arma::mat& igroup,
                      const arma::vec& nvec,
                      const arma::vec& stau,
                      const arma::vec& sm,
                      const int& ngroup,
                      const int& colVars,
                      const int& n,
                      const int& ntau){
  // compute gP, cumsumgP, and IyP
  Rcpp::List tp = fgPIyP(y, G, d, igroup, nvec, ngroup, n);
  Rcpp::List lgP       = tp["gP"];
  Rcpp::List lcumsumgP = tp["cumsumgP"];
  Rcpp::List lIyP      = tp["IyP"];
  
  // compute w2, pi1, and pi2
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w2, pi1, pi2, lgP, lcumsumgP, lIyP, d, stau, sm, n, ntau);
  
  Rcpp::List QtV(colVars + 1);
  
  // Qtau(y)
  arma::mat QtVk(n, ntau);
  for(int i(0); i < n; ++ i){
    QtVk.row(i) = arma::trans(y.elem(pi1.col(i))%(1 - w2.col(i)) +  y.elem(pi2.col(i))%w2.col(i));
  }
  QtV[0]        = QtVk;
  
  // Qtau(Vars)
  for(int k(0); k < colVars; ++ k){
    arma::vec Vk(Vars.col(k));
    for(int i(0); i < n; ++ i){
      QtVk.row(i) = arma::trans(Vk.elem(pi1.col(i))%(1 - w2.col(i)) +  Vk.elem(pi2.col(i))%w2.col(i));
    }
    QtV[k + 1]    = QtVk;
  }
  return QtV;
}

// The same function for a single variable
//[[Rcpp::export]]
arma::mat fQtauy(const arma::vec& y,
                 Rcpp::List& G,
                 const arma::vec d,
                 const arma::mat& igroup,
                 const arma::vec& nvec,
                 const arma::vec& stau,
                 const arma::vec& sm,
                 const int& ngroup,
                 const int& n,
                 const int& ntau){
  // compute gP, cumsumgP, and IyP
  Rcpp::List tp = fgPIyP(y, G, d, igroup, nvec, ngroup, n);
  Rcpp::List lgP       = tp["gP"];
  Rcpp::List lcumsumgP = tp["cumsumgP"];
  Rcpp::List lIyP      = tp["IyP"];
  
  // compute w2, pi1, and pi2
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w2, pi1, pi2, lgP, lcumsumgP, lIyP, d, stau, sm, n, ntau);
  
  arma::mat Qty(ntau, n);
  for(int i(0); i < n; ++ i){
    Qty.col(i) = y.elem(pi1.col(i))%(1 - w2.col(i)) +  y.elem(pi2.col(i))%w2.col(i);
  }
  return Qty.t();
}

// Best response function
// talpha is alpha tilde
//[[Rcpp::export]]
arma::vec BR(const arma::vec& talpha,
             const arma::mat& Qtauy,
             const arma::vec& lambdatau){
  return talpha + Qtauy*lambdatau;
}

// Nash Equilibrium
// y is initial solution
//[[Rcpp::export]]
int fNashE(arma::vec& y,
           Rcpp::List& G,
           const arma::vec d,
           const arma::vec& talpha,
           const arma::vec& lambdatau,
           const arma::mat& igroup,
           const arma::vec& nvec,
           const arma::vec& stau,
           const arma::vec& sm,
           const int& ngroup,
           const int& n,
           const int& ntau,
           const double& tol = 1e-10,
           const int& maxit  = 500){
  int t(0);
  computeBR: ++t;
  
  // Compute Qtauy
  arma::mat Qtauy = fQtauy(y, G, d, igroup, nvec, stau, sm, ngroup, n, ntau);
  
  // New y
  arma::vec yst   = BR(talpha, Qtauy, lambdatau);
  
  // check convergence
  double dist     = max(arma::abs((yst - y)/(y + 1e-50)));
  y               = yst;
  if (dist > tol && t < maxit) goto computeBR;
  return t; 
}