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

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
//#define NDEBUG
// #include <RcppNumerical.h>
// #include <RcppEigen.h>

// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;

// using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// floorP and ceilP are floor and ceil functions with the precision that
// decimal digit longer than k are ignored.
int floorP(const long double& x, const unsigned int& k = 10){
  long tp1(pow(10, k));
  long tp2(round(x*tp1));
  return tp2/tp1;
}

int ceilP(const long double& x, const unsigned int& k = 10){
  long tp1(pow(10, k));
  long tp2(round(x*tp1));
  long tp3(tp2/tp1);
  if (tp3*tp1 < tp2) return tp3 + 1;
  return tp3;
}


// gP vector of gij for gij>0, this is sorted from the smallest based of yj of peers
// IyP includes the indexes of sorterd yj of peers in the original y
// cumsumgP cumulative sum of gP
// This function computes list of gP, IyP, cumsumgP
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
      
      // Index
      tp    = ipGri.elem(tp) + n1;
      tp    = arma::join_cols(arma::join_cols(tp.head(1), tp), tp.tail(1));
      
      // cumsum
      csgPi = arma::join_cols(arma::zeros<arma::vec>(1), arma::cumsum(Gri));
      csgPi.elem(arma::find(csgPi >= 1)).ones(); 
      
      // Weights
      Gri          = arma::join_cols(Gri, arma::ones<arma::vec>(1));
      lgP[k]       = Gri;
      lcumsumgP[k] = csgPi;
      lIyP[k]      = tp;
    }
  }
  return Rcpp::List::create(_["gP"] = lgP, _["cumsumgP"] = lcumsumgP, _["IyP"] = lIyP);
}

// This function computes the weights and the indexes to dertermine quantiles
// weight for y_(pi + 1) in w1 and w2 (matrix n*ntau) 
// index for y_(pi) in pi1 (matrix n*ntau)
// index for y_(pi + 1) in pi2 (matrix n*ntau) 
void fQWeightIndex(arma::mat& w1,
                   arma::mat& w2, 
                   arma::umat& pi1,
                   arma::umat& pi2,
                   Rcpp::List& lgP, 
                   Rcpp::List& lcumsumgP, 
                   Rcpp::List& lIyP,
                   const arma::vec& d,
                   const arma::vec& stau,
                   const int& n,
                   const int& ntau,
                   const int type){
  if(type == 1){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + (stau(k) - tp1)/gP(l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = ceilP(tp2 - pii1(k));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 2){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + (stau(k) - tp1)/gP(l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = 0.5*(1 + ceilP(tp2 - pii1(k)));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 3){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l - 0.5 + (stau(k) - tp1)/gP(l)); 
        pii1(k)  = floorP(tp2); if (pii1(k) < 0) pii1(k) = 0;//this is j in HYNDMAN and FAN (1996)
        if((pii1(k)%2 == 0) & (tp2 >= 1)){
          w2i(k) = ceilP(tp2 - pii1(k));
        } else{
          w2i(k) = 1;
        }
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 4){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + (stau(k) - tp1)/gP(l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 5){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + 0.5 + (stau(k) - tp1)/gP(l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 6){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + (stau(k) - tp1)/gP(l) + stau(k)); 
        if(stau(k) >= 1) tp2 = l;
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
        // if(i == 1) cout<<d(i)<<endl;
        // if(i == 1) cout<<arma::trans(cumsumgP - stau(k))<<endl;
        // if(i == 1) cout<<tp1<<endl;
        // if(i == 1) cout<<tp2<<endl;
        // if(i == 1) cout<<pii1(k)<<endl;
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 7){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + (stau(k) - tp1)/gP(l) + 1 - stau(k)); //here tp2 is necessarily >= 0
        pii1(k) = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)  = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 8){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + (stau(k) - tp1)/gP(l) + (stau(k) + 1.0)/3.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
  
  if(type == 9){
    for(int i(0); i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec gP(lgP[i]);
      arma::vec cumsumgP(lcumsumgP[i]);
      arma::uvec IyP(lIyP[i]);
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(cumsumgP <= stau(k)) - 1);
        long double tp1(cumsumgP(l));
        long double tp2(l + (stau(k) - tp1)/gP(l) + stau(k)/4.0 + 3.0/8.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = IyP.elem(pii1);
      pi2.col(i) = IyP.elem(pii1 + 1);
    }
  }
}

// This function computes Qtau(y) 
//[[Rcpp::export]]
arma::mat fQtauy(const arma::vec& y,
                 Rcpp::List& G,
                 const arma::vec d,
                 const arma::mat& igroup,
                 const arma::vec& nvec,
                 const arma::vec& stau,
                 const int& ngroup,
                 const int& n,
                 const int& ntau,
                 const int& type){
  // compute gP, cumsumgP, and IyP
  Rcpp::List tp = fgPIyP(y, G, d, igroup, nvec, ngroup, n);
  Rcpp::List lgP       = tp["gP"];
  Rcpp::List lcumsumgP = tp["cumsumgP"];
  Rcpp::List lIyP      = tp["IyP"];
  
  // compute w1, w2, pi1, and pi2
  arma::mat w1(ntau, n, arma::fill::zeros); 
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w1, w2, pi1, pi2, lgP, lcumsumgP, lIyP, d, stau, n, ntau, type);
  arma::mat Qty(ntau, n, arma::fill::zeros);
  for(int i(0); i < n; ++ i){
    Qty.col(i) = y.elem(pi1.col(i))%w1.col(i) +  y.elem(pi2.col(i))%w2.col(i);
  }
  return Qty.t();
}

// The same function with indices
//[[Rcpp::export]]
Rcpp::List fQtauyWithIndex(const arma::vec& y,
                           Rcpp::List& G,
                           const arma::vec d,
                           const arma::mat& igroup,
                           const arma::vec& nvec,
                           const arma::vec& stau,
                           const int& ngroup,
                           const int& n,
                           const int& ntau,
                           const int& type){
  // compute gP, cumsumgP, and IyP
  Rcpp::List tp = fgPIyP(y, G, d, igroup, nvec, ngroup, n);
  Rcpp::List lgP       = tp["gP"];
  Rcpp::List lcumsumgP = tp["cumsumgP"];
  Rcpp::List lIyP      = tp["IyP"];
  
  // compute w1, w2, pi1, and pi2
  arma::mat w1(ntau, n, arma::fill::zeros); 
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w1, w2, pi1, pi2, lgP, lcumsumgP, lIyP, d, stau, n, ntau, type);
  
  arma::mat Qty(ntau, n);
  for(int i(0); i < n; ++ i){
    Qty.col(i) = y.elem(pi1.col(i))%w1.col(i) +  y.elem(pi2.col(i))%w2.col(i);
  }
  
  return Rcpp::List::create(_["qy"] = Qty.t(), _["pi1"] = pi1.t(), _["pi2"] = pi2.t(), 
                            _["w1"] = w1.t(), _["w2"] = w2.t());
}

// The same function only indices
//[[Rcpp::export]]
Rcpp::List fQtauyIndex(const arma::vec& y,
                       Rcpp::List& G,
                       const arma::vec d,
                       const arma::mat& igroup,
                       const arma::vec& nvec,
                       const arma::vec& stau,
                       const int& ngroup,
                       const int& n,
                       const int& ntau,
                       const int& type){
  // compute gP, cumsumgP, and IyP
  Rcpp::List tp = fgPIyP(y, G, d, igroup, nvec, ngroup, n);
  Rcpp::List lgP       = tp["gP"];
  Rcpp::List lcumsumgP = tp["cumsumgP"];
  Rcpp::List lIyP      = tp["IyP"];
  
  // compute w1, w2, pi1, and pi2
  arma::mat w1(ntau, n, arma::fill::zeros); 
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w1, w2, pi1, pi2, lgP, lcumsumgP, lIyP, d, stau, n, ntau, type);
  return Rcpp::List::create(_["pi1"] = pi1.t(), _["pi2"] = pi2.t(), 
                            _["w1"] = w1.t(), _["w2"] = w2.t());
}

// Convert index into sparse matrix
//[[Rcpp::export]]
arma::sp_mat fIndexMat(const arma::uvec& pi1,
                       const arma::uvec& pi2,
                       const arma::vec& w1,
                       const arma::vec& w2,
                       const int& n){
  arma::sp_mat out(n, n);
  for (int i(0); i < n; ++ i){
    out(i, pi1(i))  = w1(i);
    out(i, pi2(i)) += w2(i);
  }
  return out;
}

// Product Weights Variables
//[[Rcpp::export]]
arma::mat fProdWVI(const arma::sp_mat& W,
                   const arma::mat& V,
                   const int& distance = 1){
  arma::mat tp(W*V);
  if (distance == 1) return tp;
  arma::mat out(tp);
  for (int k(1); k < distance; ++ k){
    tp  = W*tp;
    out = arma::join_rows(out, tp);
  }
  return out;
}

// // Compute expectation of y
// // talpha is alpha tilde
// // sWl is the sum of W*lambda
// //[[Rcpp::export]]
// arma::vec fEy(const arma::sp_mat& sWl, 
//               const arma::vec& talpha,
//               const int& n){
//   return arma::spsolve(arma::speye<arma::sp_mat>(n, n) - sWl, talpha);
// }

// Best response function
// talpha is alpha tilde
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
           const int& ngroup,
           const int& n,
           const int& ntau,
           const int& type = 7,
           const double& tol = 1e-10,
           const int& maxit  = 500){
  int t(0);
  computeBR: ++t;
  
  // Compute Qtauy
  arma::mat Qtauy = fQtauy(y, G, d, igroup, nvec, stau, ngroup, n, ntau, type);
  
  // New y
  arma::vec yst   = BR(talpha, Qtauy, lambdatau);
  
  // check convergence
  double dist     = max(arma::abs((yst - y)/(y + 1e-50)));
  y               = yst;
  if (dist > tol && t < maxit) goto computeBR;
  return t; 
}

//  This function compute optimal instruments (reduced form model)
//[[Rcpp::export]]
arma::mat optins_red(const arma::vec& beta,
                     const arma::vec& y,
                     Rcpp::List& G,
                     const arma::mat& X,
                     const arma::vec d,
                     const arma::mat& igroup,
                     const arma::vec& nvec,
                     const arma::vec& stau,
                     const int& ngroup,
                     const int& n,
                     const int& ntau,
                     const int& type,
                     const int& Kx,
                     const double& tol,
                     const int& maxit) {
  arma::vec talpha(X*beta.tail(Kx));
  arma::vec Ey = y;
  fNashE(Ey, G, d, talpha, beta.head(ntau), igroup, nvec, stau, ngroup, n, 
                 ntau, type, tol, maxit);
  return fQtauy(Ey, G, d, igroup, nvec, stau, ngroup, n, ntau, type);
}

//  This function compute optimal instruments (reduced form model)
//[[Rcpp::export]]
arma::mat optins_struc(const arma::vec& beta,
                     const arma::vec& y,
                     Rcpp::List& G,
                     const arma::mat& X,
                     const arma::vec d,
                     const arma::mat& igroup,
                     const arma::vec& nvec,
                     const arma::vec& stau,
                     const arma::uvec& nIs,
                     const int& ngroup,
                     const int& n,
                     const int& ntau,
                     const int& type,
                     const int& Kx,
                     const double& tol,
                     const int& maxit) {
  arma::vec talpha(X*beta.tail(Kx)); talpha.elem(nIs) *= (1 - beta(0));
  arma::vec Ey = y;
  fNashE(Ey, G, d, talpha, beta.subvec(1, ntau), igroup, nvec, stau, ngroup, n, 
         ntau, type, tol, maxit);
  return fQtauy(Ey, G, d, igroup, nvec, stau, ngroup, n, ntau, type);
}
