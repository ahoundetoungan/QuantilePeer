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

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

//===============================================================================================//
////////////////////////////////////// Structures ////////////////////////////////////////////////
//===============================================================================================//
//===============================================================================================//
struct gPIyP {
  std::vector<arma::vec> gP;
  std::vector<arma::vec> cumsumgP;
  std::vector<arma::uvec> IyP;
  
  gPIyP(int n) {
    gP.resize(n);
    cumsumgP.resize(n);
    IyP.resize(n);
  }
};

struct gPIyP_EIGEN {
  std::vector<Eigen::ArrayXd> gP;
  std::vector<Eigen::ArrayXd> cumsumgP;
  std::vector<Eigen::ArrayXi> IyP;
  
  gPIyP_EIGEN(int n) {
    gP.resize(n);
    cumsumgP.resize(n);
    IyP.resize(n);
  }
};

//===============================================================================================//
//===============================================================================================//

// floorP and ceilP are floor and ceil functions with the precision that
// decimal digit longer than k are ignored.
int floorP(const long double& x, const unsigned int& k = 10){
  long long int tp1(std::pow(10, k));
  long long int tp2(std::round(x * tp1));
  return tp2/tp1;
}

int ceilP(const long double& x, const unsigned int& k = 10){
  long long int tp1(std::pow(10, k));
  long long int tp2(std::round(x * tp1));
  int tp3(tp2/tp1);
  if (tp3*tp1 < tp2) return tp3 + 1;
  return tp3;
}

// This function set nthreads
//[[Rcpp::export]]
int fnthreads(const int& nthreads) {
#if defined(_OPENMP)
  return nthreads;
#else
  if (nthreads > 1) {
    Rf_warning("OpenMP is not available. Sequential processing is used.");
  }
  return 1;
#endif
}


// gP vector of gij for gij>0, this is sorted from the smallest based of yj of peers
// IyP includes the indexes of sorterd yj of peers in the original y
// cumsumgP cumulative sum of gP
// This function computes list of gP, IyP, cumsumgP
gPIyP fgPIyP(const arma::vec& y,
             const std::vector<arma::mat>& G,
             const arma::vec& d,
             const arma::uvec& igroup,
             const arma::uvec& group, // group index
             const arma::uvec& groupidx, // position in the group
             const arma::uvec& nvec,
             const int& ngroup,
             const int& n,
             const int& nthreads){
  gPIyP out(n);
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel
{
  arma::vec Gri;
  arma::uvec ipGri, tp;
  int i, r, n1, nfr;
#pragma omp for
  for(int k = 0; k < n; ++ k){
    if(d(k) < 1e-50) continue;
    i  = groupidx(k); 
    r  = group(k);
    n1 = igroup(r);

      Gri   = G[r].row(i).t();
      ipGri = arma::find(Gri > 1e-50);
      Gri   = Gri.elem(ipGri) / d(k);
      tp    = arma::sort_index(y.elem(ipGri + n1));
      Gri   = Gri.elem(tp);
      
      // Index
      tp    = ipGri.elem(tp) + n1;
      nfr   = tp.n_elem;
      arma::uvec idx(nfr + 2);
      idx(0)             = tp(0);
      idx.subvec(1, nfr) = tp;
      idx(nfr + 1)       = tp(nfr - 1);
      
      // cumsum
      arma::vec csgPi(nfr + 1);
      csgPi(0)        = 0;
      csgPi.tail(nfr) = arma::cumsum(Gri);
      
      // Weights
      arma::vec gP(nfr + 1);
      gP.head(nfr) = Gri;
      gP(nfr)      = 1;
      
      out.gP[k] = gP;
      out.cumsumgP[k] = csgPi;
      out.IyP[k]      = idx;
    }
}
#else
  arma::vec Gri;
  arma::uvec ipGri, tp;
  int i, r, n1, nfr;
  for(int k = 0; k < n; ++ k){
    if(d(k) < 1e-50) continue;
    i  = groupidx(k); 
    r  = group(k);
    n1 = igroup(r);
    
    Gri   = G[r].row(i).t();
    ipGri = arma::find(Gri > 1e-50);
    Gri   = Gri.elem(ipGri) / d(k);
    tp    = arma::sort_index(y.elem(ipGri + n1));
    Gri   = Gri.elem(tp);
    
    // Index
    tp    = ipGri.elem(tp) + n1;
    nfr   = tp.n_elem;
    arma::uvec idx(nfr + 2);
    idx(0)             = tp(0);
    idx.subvec(1, nfr) = tp;
    idx(nfr + 1)       = tp(nfr - 1);
    
    // cumsum
    arma::vec csgPi(nfr + 1);
    csgPi(0)        = 0;
    csgPi.tail(nfr) = arma::cumsum(Gri);
    
    // Weights
    arma::vec gP(nfr + 1);
    gP.head(nfr) = Gri;
    gP(nfr)      = 1;
    
    out.gP[k] = gP;
    out.cumsumgP[k] = csgPi;
    out.IyP[k]      = idx;
  }
#endif
  return out;
}


gPIyP_EIGEN fgPIyP_EIGEN(const Eigen::ArrayXd& y,
                         const std::vector<Eigen::ArrayXXd>& G,
                         const Eigen::ArrayXd& d,
                         const Eigen::ArrayXi& igroup,
                         const Eigen::ArrayXi& group, // group index
                         const Eigen::ArrayXi& groupidx, // position in the group
                         const Eigen::ArrayXi& nvec,
                         const int& ngroup,
                         const int& n,
                         const int& nthreads){
  gPIyP_EIGEN out(n);
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for(int k = 0; k < n; ++ k){
    if(d(k) < 1e-50) continue;
    int i(groupidx(k)), r(group(k)), n1(igroup(r));
    
    std::vector<int> ipGriVec;
    for (int j(0); j < nvec(r); ++ j){
      if (G[r](i, j) > 1e-50) {
        ipGriVec.push_back(j);
      }
    }
    int nfr(ipGriVec.size());
    Eigen::ArrayXi ipGri = Eigen::Map<Eigen::ArrayXi>(ipGriVec.data(), nfr);
    Eigen::ArrayXd Gri(G[r](i, ipGri).transpose() / d(k));
    
    Eigen::ArrayXi tp(Eigen::ArrayXi::LinSpaced(nfr, 0, nfr - 1));
    std::sort(tp.data(), tp.data() + nfr,
              [&y, &ipGri, &n1](int a, int b) {
                return y(ipGri(a) + n1) < y(ipGri(b) + n1);
              });
    Gri   = Gri(tp);
    
    // Index
    tp    = ipGri(tp) + n1;
    Eigen::ArrayXi idx(nfr + 2);
    idx << tp(0), tp, tp(nfr - 1);
    
    // cumsum
    Eigen::ArrayXd csgPi(Eigen::ArrayXd::Zero(nfr + 1));
    for (int j(0); j < nfr; ++ j) {
      csgPi(j + 1) = std::min(1.0, csgPi(j) + Gri(j));
    }
    
    // Weights
    Eigen::ArrayXd gP(nfr + 1);
    gP << Gri, 1;
    out.gP[k] = gP;
    out.cumsumgP[k] = csgPi;
    out.IyP[k]      = idx;
  }
#else
  for(int k = 0; k < n; ++ k){
    if(d(k) < 1e-50) continue;
    int i(groupidx(k)), r(group(k)), n1(igroup(r));
    
    std::vector<int> ipGriVec;
    for (int j(0); j < nvec(r); ++ j){
      if (G[r](i, j) > 1e-50) {
        ipGriVec.push_back(j);
      }
    }
    int nfr(ipGriVec.size());
    Eigen::ArrayXi ipGri = Eigen::Map<Eigen::ArrayXi>(ipGriVec.data(), nfr);
    Eigen::ArrayXd Gri(G[r](i, ipGri).transpose() / d(k));
    
    Eigen::ArrayXi tp(Eigen::ArrayXi::LinSpaced(nfr, 0, nfr - 1));
    std::sort(tp.data(), tp.data() + nfr,
              [&y, &ipGri, &n1](int a, int b) {
                return y(ipGri(a) + n1) < y(ipGri(b) + n1);
              });
    Gri   = Gri(tp);
    
    // Index
    tp    = ipGri(tp) + n1;
    Eigen::ArrayXi idx(nfr + 2);
    idx << tp(0), tp, tp(nfr - 1);
    
    // cumsum
    Eigen::ArrayXd csgPi(Eigen::ArrayXd::Zero(nfr + 1));
    for (int j(0); j < nfr; ++ j) {
      csgPi(j + 1) = std::min(1.0, csgPi(j) + Gri(j));
    }
    
    // Weights
    Eigen::ArrayXd gP(nfr + 1);
    gP << Gri, 1;
    out.gP[k] = gP;
    out.cumsumgP[k] = csgPi;
    out.IyP[k]      = idx;
  }
#endif
  return out;
}

// This function computes the weights and the indexes to dertermine quantiles
// weight for y_(pi + 1) in w1 and w2 (matrix n*ntau) 
// index for y_(pi) in pi1 (matrix n*ntau)
// index for y_(pi + 1) in pi2 (matrix n*ntau) 
void fQWeightIndex(arma::mat& w1,
                   arma::mat& w2, 
                   arma::umat& pi1,
                   arma::umat& pi2,
                   std::vector<arma::vec>& lgP, 
                   std::vector<arma::vec>& lcumsumgP, 
                   std::vector<arma::uvec>& lIyP,
                   const arma::vec& d,
                   const arma::vec& stau,
                   const int& n,
                   const int& ntau,
                   const int& type,
                   const int& nthreads){
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
  if(type == 1){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = ceilP(tp2 - pii1(k));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 2){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = 0.5*(1 + ceilP(tp2 - pii1(k)));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 3){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l - 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); if (pii1(k) < 0) pii1(k) = 0;//this is j in HYNDMAN and FAN (1996)
        if((pii1(k)%2 == 0) & (tp2 >= 1)){
          w2i(k) = ceilP(tp2 - pii1(k));
        } else{
          w2i(k) = 1;
        }
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 4){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 5){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 6){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)); 
        if(stau(k) >= 1) tp2 = l;
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 7){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + 1 - stau(k)); //here tp2 is necessarily >= 0
        pii1(k) = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)  = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 8){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + (stau(k) + 1.0)/3.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 9){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)/4.0 + 3.0/8.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
#else
  if(type == 1){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = ceilP(tp2 - pii1(k));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 2){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = 0.5*(1 + ceilP(tp2 - pii1(k)));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 3){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l - 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); if (pii1(k) < 0) pii1(k) = 0;//this is j in HYNDMAN and FAN (1996)
        if((pii1(k)%2 == 0) & (tp2 >= 1)){
          w2i(k) = ceilP(tp2 - pii1(k));
        } else{
          w2i(k) = 1;
        }
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 4){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 5){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 6){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)); 
        if(stau(k) >= 1) tp2 = l;
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 7){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + 1 - stau(k)); //here tp2 is necessarily >= 0
        pii1(k) = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)  = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 8){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + (stau(k) + 1.0)/3.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
  
  if(type == 9){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      arma::vec w2i(ntau);
      arma::uvec pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l(sum(lcumsumgP[i] <= stau(k)) - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)/4.0 + 3.0/8.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i].elem(pii1);
      pi2.col(i) = lIyP[i].elem(pii1 + 1);
    }
  }
#endif
}

void fQWeightIndex_EIGEN(Eigen::ArrayXXd& w1,
                         Eigen::ArrayXXd& w2, 
                         Eigen::ArrayXXi& pi1,
                         Eigen::ArrayXXi& pi2,
                         std::vector<Eigen::ArrayXd>& lgP, 
                         std::vector<Eigen::ArrayXd>& lcumsumgP, 
                         std::vector<Eigen::ArrayXi>& lIyP,
                         const Eigen::ArrayXd& d,
                         const Eigen::ArrayXd& stau,
                         const int& n,
                         const int& ntau,
                         const int& type,
                         const int& nthreads){
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
  if(type == 1){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = ceilP(tp2 - pii1(k));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 2){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = 0.5*(1 + ceilP(tp2 - pii1(k)));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 3){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l - 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); if (pii1(k) < 0) pii1(k) = 0;//this is j in HYNDMAN and FAN (1996)
        if((pii1(k)%2 == 0) & (tp2 >= 1)){
          w2i(k) = ceilP(tp2 - pii1(k));
        } else{
          w2i(k) = 1;
        }
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 4){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 5){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 6){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)); 
        if(stau(k) >= 1) tp2 = l;
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 7){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + 1 - stau(k)); //here tp2 is necessarily >= 0
        pii1(k) = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)  = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 8){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + (stau(k) + 1.0)/3.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 9){
#pragma omp parallel for
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)/4.0 + 3.0/8.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
#else
  if(type == 1){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = ceilP(tp2 - pii1(k));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 2){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = 0.5*(1 + ceilP(tp2 - pii1(k)));
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 3){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l - 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); if (pii1(k) < 0) pii1(k) = 0;//this is j in HYNDMAN and FAN (1996)
        if((pii1(k)%2 == 0) & (tp2 >= 1)){
          w2i(k) = ceilP(tp2 - pii1(k));
        } else{
          w2i(k) = 1;
        }
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 4){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 5){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + 0.5 + (stau(k) - tp1)/lgP[i](l)); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 6){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)); 
        if(stau(k) >= 1) tp2 = l;
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 7){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + 1 - stau(k)); //here tp2 is necessarily >= 0
        pii1(k) = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)  = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 8){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + (stau(k) + 1.0)/3.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
  
  if(type == 9){
    for(int i = 0; i < n; ++ i){
      if(d(i) == 0) continue;
      Eigen::ArrayXd w2i(ntau);
      Eigen::ArrayXi pii1(ntau);
      
      for(int k(0); k < ntau; ++ k){
        int l((lcumsumgP[i] <= stau(k)).sum() - 1);
        long double tp1(lcumsumgP[i](l));
        long double tp2(l + (stau(k) - tp1)/lgP[i](l) + stau(k)/4.0 + 3.0/8.0); 
        pii1(k)  = floorP(tp2); //this is j in HYNDMAN and FAN (1996)
        w2i(k)   = tp2 - pii1(k);
      }
      w1.col(i)  = 1 - w2i;
      w2.col(i)  = w2i;
      pi1.col(i) = lIyP[i](pii1);
      pi2.col(i) = lIyP[i](pii1 + 1);
    }
  }
#endif
}

// This function computes Qtau(y) 
//[[Rcpp::export]]
arma::mat fQtauy(const arma::vec& y,
                 const std::vector<arma::mat>& G,
                 const arma::vec& d,
                 const arma::uvec& igroup,
                 const arma::uvec& group, // group index
                 const arma::uvec& groupidx, // position in the group
                 const arma::uvec& nvec,
                 const arma::vec& stau,
                 const int& ngroup,
                 const int& n,
                 const int& ntau,
                 const int& type,
                 const int& nthreads){
  // compute gP, cumsumgP, and IyP
  gPIyP tp = fgPIyP(y, G, d, igroup, group, groupidx, nvec, ngroup, n, nthreads);
  
  // compute w1, w2, pi1, and pi2
  arma::mat w1(ntau, n, arma::fill::zeros); 
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w1, w2, pi1, pi2, tp.gP, tp.cumsumgP, tp.IyP, d, stau, n, ntau, type, nthreads);
  arma::mat Qty(ntau, n);
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y.elem(pi1.col(i))%w1.col(i) +  y.elem(pi2.col(i))%w2.col(i);
  }
#else
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y.elem(pi1.col(i))%w1.col(i) +  y.elem(pi2.col(i))%w2.col(i);
  }
#endif
  return Qty.t();
}


//[[Rcpp::export]]
Eigen::ArrayXXd fQtauy_EIGEN(const Eigen::ArrayXd& y,
                             const std::vector<Eigen::ArrayXXd>& G,
                             const Eigen::ArrayXd& d,
                             const Eigen::ArrayXi& igroup,
                             const Eigen::ArrayXi& group, // group index
                             const Eigen::ArrayXi& groupidx, // position in the group
                             const Eigen::ArrayXi& nvec,
                             const Eigen::ArrayXd& stau,
                             const int& ngroup,
                             const int& n,
                             const int& ntau,
                             const int& type,
                             const int& nthreads){
  // compute gP, cumsumgP, and IyP
  gPIyP_EIGEN tp = fgPIyP_EIGEN(y, G, d, igroup, group, groupidx, nvec, ngroup, n, nthreads);
  
  // compute w1, w2, pi1, and pi2
  Eigen::ArrayXXd w1(Eigen::ArrayXXd::Zero(ntau, n)), w2(Eigen::ArrayXXd::Zero(ntau, n)); 
  Eigen::ArrayXXi pi1(Eigen::ArrayXXi::Zero(ntau, n)), pi2(Eigen::ArrayXXi::Zero(ntau, n));
  fQWeightIndex_EIGEN(w1, w2, pi1, pi2, tp.gP, tp.cumsumgP, tp.IyP, d, stau, n, ntau, type, nthreads);
  Eigen::ArrayXXd Qty(ntau, n);
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y(pi1.col(i)) * w1.col(i) +  y(pi2.col(i)) * w2.col(i);
  }
#else
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y(pi1.col(i)) * w1.col(i) +  y(pi2.col(i)) * w2.col(i);
  }
#endif
  return Qty.transpose();
}

// The same function with indices
//[[Rcpp::export]]
Rcpp::List fQtauyWithIndex(const arma::vec& y,
                           const std::vector<arma::mat>& G,
                           const arma::vec& d,
                           const arma::uvec& igroup,
                           const arma::uvec& group, // group index
                           const arma::uvec& groupidx, // position in the group
                           const arma::uvec& nvec,
                           const arma::vec& stau,
                           const int& ngroup,
                           const int& n,
                           const int& ntau,
                           const int& type,
                           const int& nthreads){
  // compute gP, cumsumgP, and IyP
  gPIyP tp = fgPIyP(y, G, d, igroup, group, groupidx, nvec, ngroup, n, nthreads);
  
  // compute w1, w2, pi1, and pi2
  arma::mat w1(ntau, n, arma::fill::zeros); 
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w1, w2, pi1, pi2, tp.gP, tp.cumsumgP, tp.IyP, d, stau, n, ntau, type, nthreads);
  arma::mat Qty(ntau, n);
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y.elem(pi1.col(i))%w1.col(i) +  y.elem(pi2.col(i))%w2.col(i);
  }
#else
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y.elem(pi1.col(i))%w1.col(i) +  y.elem(pi2.col(i))%w2.col(i);
  }
#endif

  return Rcpp::List::create(Rcpp::_["qy"]  = Qty.t(), 
                            Rcpp::_["pi1"] = pi1.t(), 
                            Rcpp::_["pi2"] = pi2.t(), 
                            Rcpp::_["w1"]  = w1.t(), 
                            Rcpp::_["w2"]  = w2.t());
}

//[[Rcpp::export]]
Rcpp::List fQtauyWithIndex_EIGEN(const Eigen::ArrayXd& y,
                                 const std::vector<Eigen::ArrayXXd>& G,
                                 const Eigen::ArrayXd& d,
                                 const Eigen::ArrayXi& igroup,
                                 const Eigen::ArrayXi& group, // group index
                                 const Eigen::ArrayXi& groupidx, // position in the group
                                 const Eigen::ArrayXi& nvec,
                                 const Eigen::ArrayXd& stau,
                                 const int& ngroup,
                                 const int& n,
                                 const int& ntau,
                                 const int& type,
                                 const int& nthreads){
  // compute gP, cumsumgP, and IyP
  gPIyP_EIGEN tp = fgPIyP_EIGEN(y, G, d, igroup, group, groupidx, nvec, ngroup, n, nthreads);
  
  // compute w1, w2, pi1, and pi2
  Eigen::ArrayXXd w1(Eigen::ArrayXXd::Zero(ntau, n)), w2(Eigen::ArrayXXd::Zero(ntau, n)); 
  Eigen::ArrayXXi pi1(Eigen::ArrayXXi::Zero(ntau, n)), pi2(Eigen::ArrayXXi::Zero(ntau, n));
  fQWeightIndex_EIGEN(w1, w2, pi1, pi2, tp.gP, tp.cumsumgP, tp.IyP, d, stau, n, ntau, type, nthreads);
  Eigen::ArrayXXd Qty(ntau, n);
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y(pi1.col(i)) * w1.col(i) +  y(pi2.col(i)) * w2.col(i);
  }
#else
  for(int i = 0; i < n; ++ i){
    Qty.col(i) = y(pi1.col(i)) * w1.col(i) +  y(pi2.col(i)) * w2.col(i);
  }
#endif
  
  return Rcpp::List::create(Rcpp::_["qy"]  = Qty.transpose(), 
                            Rcpp::_["pi1"] = pi1.transpose(), 
                            Rcpp::_["pi2"] = pi2.transpose(), 
                            Rcpp::_["w1"]  = w1.transpose(), 
                            Rcpp::_["w2"]  = w2.transpose());
}

// The same function only indices
//[[Rcpp::export]]
Rcpp::List fQtauyIndex(const arma::vec& y,
                       const std::vector<arma::mat>& G,
                       const arma::vec& d,
                       const arma::uvec& igroup,
                       const arma::uvec& group, // group index
                       const arma::uvec& groupidx, // position in the group
                       const arma::uvec& nvec,
                       const arma::vec& stau,
                       const int& ngroup,
                       const int& n,
                       const int& ntau,
                       const int& type,
                       const int& nthreads){
  // compute gP, cumsumgP, and IyP
  gPIyP tp = fgPIyP(y, G, d, igroup, group, groupidx, nvec, ngroup, n, nthreads);
  
  // compute w1, w2, pi1, and pi2
  arma::mat w1(ntau, n, arma::fill::zeros); 
  arma::mat w2(ntau, n, arma::fill::zeros); 
  arma::umat pi1(ntau, n, arma::fill::zeros), pi2(ntau, n, arma::fill::zeros);
  fQWeightIndex(w1, w2, pi1, pi2, tp.gP, tp.cumsumgP, tp.IyP, d, stau, n, ntau, type, nthreads);
  return Rcpp::List::create(Rcpp::_["pi1"] = pi1.t(), 
                            Rcpp::_["pi2"] = pi2.t(), 
                            Rcpp::_["w1"]  = w1.t(), 
                            Rcpp::_["w2"]  = w2.t());
}

//[[Rcpp::export]]
Rcpp::List fQtauyIndex_EIGEN(const Eigen::ArrayXd& y,
                             const std::vector<Eigen::ArrayXXd>& G,
                             const Eigen::ArrayXd& d,
                             const Eigen::ArrayXi& igroup,
                             const Eigen::ArrayXi& group, // group index
                             const Eigen::ArrayXi& groupidx, // position in the group
                             const Eigen::ArrayXi& nvec,
                             const Eigen::ArrayXd& stau,
                             const int& ngroup,
                             const int& n,
                             const int& ntau,
                             const int& type,
                             const int& nthreads){
  // compute gP, cumsumgP, and IyP
  gPIyP_EIGEN tp = fgPIyP_EIGEN(y, G, d, igroup, group, groupidx, nvec, ngroup, n, nthreads);
  
  // compute w1, w2, pi1, and pi2
  Eigen::ArrayXXd w1(Eigen::ArrayXXd::Zero(ntau, n)), w2(Eigen::ArrayXXd::Zero(ntau, n)); 
  Eigen::ArrayXXi pi1(Eigen::ArrayXXi::Zero(ntau, n)), pi2(Eigen::ArrayXXi::Zero(ntau, n));
  fQWeightIndex_EIGEN(w1, w2, pi1, pi2, tp.gP, tp.cumsumgP, tp.IyP, d, stau, n, ntau, type, nthreads);
  
  return Rcpp::List::create(Rcpp::_["pi1"] = pi1.transpose(), 
                            Rcpp::_["pi2"] = pi2.transpose(), 
                            Rcpp::_["w1"]  = w1.transpose(), 
                            Rcpp::_["w2"]  = w2.transpose());
}

// Convert index into sparse matrix
// [[Rcpp::export]]
arma::sp_mat fIndexMat_ARMA(const arma::uvec& pi1,
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

// Convert index into sparse matrix
// [[Rcpp::export]]
Eigen::SparseMatrix<double> fIndexMat(const Eigen::ArrayXi& pi1,
                                            const Eigen::ArrayXi& pi2,
                                            const Eigen::ArrayXd& w1,
                                            const Eigen::ArrayXd& w2,
                                            const int& n) {
  std::vector<Eigen::Triplet<double>> Tpl;
  
  for (int i = 0; i < n; ++i) {
    Tpl.emplace_back(i, pi1(i), w1(i));
    Tpl.emplace_back(i, pi2(i), w2(i));
  }
  
  // Build sparse matrix
  Eigen::SparseMatrix<double> out(n, n);
  out.setFromTriplets(Tpl.begin(), Tpl.end());
  
  return out;
}

// Product Weights Variables
//[[Rcpp::export]]
arma::mat fProdWVI(const arma::sp_mat& W,
                   const arma::mat& V,
                   const int& distance = 1){
  int n(W.n_rows), kV(V.n_cols);
  arma::mat out(n, kV * distance);
  if ((kV * distance) == 0) {
    return out;
  }
  arma::mat tp = W * V;
  out.cols(0, kV - 1) = tp;
  for (int k(1); k < distance; ++ k){
    tp  = W*tp;
    out.cols(k * kV, (k + 1) * kV - 1) = tp;
  }
  return out;
}

// [[Rcpp::export]]
Eigen::MatrixXd fProdWVI_EIGEN(const Eigen::SparseMatrix<double>& W,
                               const Eigen::MatrixXd& V,
                               const int distance = 1) {
  int n(V.rows()), kV(V.cols());
  Eigen::MatrixXd out(n, kV * distance);
  if ((kV * distance) == 0) {
    return out;
  }
  Eigen::MatrixXd tp = W * V;
  out.block(0, 0, n, kV) = tp;
  for (int k = 1; k < distance; ++k) {
    tp = W * tp;
    out.block(0, k * kV, n, kV) = tp;
  }
  return out;
}


// Nash Equilibrium
// y is initial solution
//[[Rcpp::export]]
int fNashE(arma::vec& y,
           const std::vector<arma::mat>& G,
           const arma::vec& d,
           const arma::vec& talpha,
           const arma::vec& lambdatau,
           const arma::uvec& igroup,
           const arma::uvec& group, // group index
           const arma::uvec& groupidx, // position in the group
           const arma::uvec& nvec,
           const arma::vec& stau,
           const int& ngroup,
           const int& n,
           const int& ntau,
           const int& type = 7,
           const double& tol = 1e-10,
           const int& maxit  = 500,
           const int& nthreads = 1){
  int t(0);
  computeBR: ++t;
  
  // Compute Qtauy
  arma::mat Qtauy = fQtauy(y, G, d, igroup, group, groupidx, nvec, stau, ngroup,
                           n, ntau, type, nthreads);
  
  // New y
  arma::vec yst   = talpha + Qtauy*lambdatau;
  
  // check convergence
  double dist     = max(arma::abs((yst - y)/(y + 1e-50)));
  y               = yst;
  if (dist > tol && t < maxit) goto computeBR;
  return t; 
}

//[[Rcpp::export]]
Rcpp::List fNashE_EIGEN(const Eigen::ArrayXd& y0,
                        const std::vector<Eigen::ArrayXXd>& G,
                        const Eigen::ArrayXd& d,
                        const Eigen::VectorXd& talpha,
                        const Eigen::VectorXd& lambdatau,
                        const Eigen::ArrayXi& igroup,
                        const Eigen::ArrayXi& group, // group index
                        const Eigen::ArrayXi& groupidx, // position in the group
                        const Eigen::ArrayXi& nvec,
                        const Eigen::ArrayXd& stau,
                        const int& ngroup,
                        const int& n,
                        const int& ntau,
                        const int& type = 7,
                        const double& tol = 1e-10,
                        const int& maxit  = 500,
                        const int& nthreads = 1){
  Eigen::ArrayXd y(y0);
  int t(0);
  computeBR: ++t;
  
  // Compute Qtauy
  Eigen::MatrixXd Qtauy = fQtauy_EIGEN(y, G, d, igroup, group, groupidx, nvec, stau, 
                                       ngroup, n, ntau, type, nthreads);
  
  // New y
  Eigen::ArrayXd yst    = talpha + Qtauy * lambdatau;
  
  // check convergence
  double dist     = ((yst - y).abs()/(y.abs() + 1e-50)).maxCoeff();
  y               = yst;
  if (dist > tol && t < maxit) goto computeBR;
  return Rcpp::List::create(Rcpp::_["y"] = y, Rcpp::_["t"] = t); 
}


// Optimal (good) instrument 
//[[Rcpp::export]]
arma::mat simInstrqpeer(const arma::vec& y,
                        const arma::mat& qy,
                        const arma::mat& X,
                        const std::vector<arma::mat>& G,
                        const arma::vec& d,
                        const arma::uvec& igroup,
                        const arma::uvec& group, // group index
                        const arma::uvec& groupidx, // position in the group
                        const arma::vec& estimate,
                        const arma::uvec& nIs,
                        const arma::uvec& nvec,
                        const arma::vec& stau,
                        const int& boot,
                        const bool& fixedeffects,
                        const bool& structural,
                        const unsigned int& nthreads,
                        const unsigned int& seed,
                        const int& type = 7,
                        const double& tol = 1e-10,
                        const int& maxit  = 500) {
  int n(y.n_elem), ntau(stau.n_elem), Kx(X.n_cols), ngroup(nvec.n_elem);
  
  // residuals
  arma::vec eps(y - qy * estimate.subvec(structural, ntau - 1 + structural));
  
  // remove X*beta
  arma::vec Xbeta(X * estimate.tail(Kx));
  if (structural) {
    Xbeta(nIs) *= (1 - estimate(0));
  }
  eps          -= Xbeta;
  if (structural) {
    eps(nIs)   /= (1 - estimate(0));
  }
  
  // strata setup
  std::vector<int> strata;
  strata.push_back(0);
  if (fixedeffects) {
    for (int m = 0; m < ngroup; ++m) {
      strata.push_back(strata.back() + nvec(m));
    }
  } else {
    strata.push_back(n);
  }
  int nstrata = strata.size() - 1;
  
  // output
  arma::mat out;
  
  // If it can be run in parallel
#if defined(_OPENMP)
  // list of instruments
  std::vector<arma::mat> listQtauy(nthreads);
  
  // bootstrap
  omp_set_num_threads(nthreads);
#pragma omp parallel
{
  unsigned int tid = omp_get_thread_num();
  // Initialize RNG with seed_seq for thread safety
  std::vector<unsigned int> seq_data = {seed, tid};
  std::seed_seq seq(seq_data.begin(), seq_data.end());
  std::mt19937 rng(seq);
  
  listQtauy[tid].resize(n, ntau);
  listQtauy[tid].zeros();
  
#pragma omp for
  for (int b = 0; b < boot; ++b) {
    arma::vec epsb = eps;
    for (int s = 0; s < nstrata; ++s) {
      std::shuffle(epsb.memptr() + strata[s], epsb.memptr() + strata[s + 1], rng);
    }
    if (structural) {
      epsb(nIs) *= (1 - estimate(0));
    }
    arma::vec talpha(Xbeta + epsb);
    
    arma::vec yy(y);
    int t(0);
    computeBR: ++t;
    
    // Compute Qtauy
    arma::mat Qtauy = fQtauy(yy, G, d, igroup, group, groupidx, nvec, stau, 
                             ngroup, n, ntau, type, 1);
    
    // New y
    arma::vec yyst = talpha + Qtauy * estimate.subvec(structural, ntau - 1 + structural);
    
    // check convergence
    double dist     = max(abs(yyst - yy)/(abs(yy) + 1e-50));
    yy              = yyst;
    if (dist > tol && t < maxit) goto computeBR;
    listQtauy[tid] += Qtauy;
  }
}

  out = listQtauy[0];
  for (unsigned int k = 1; k < nthreads; ++k){
    out += listQtauy[k];
  }

#else
  // Initialize RNG with seed_seq for thread safety
  std::mt19937 rng(seed);
  
  out.resize(n, ntau);
  out.zeros();
  
  for (int b = 0; b < boot; ++b) {
    arma::vec epsb = eps;
    for (int s = 0; s < nstrata; ++s) {
      std::shuffle(epsb.memptr() + strata[s], epsb.memptr() + strata[s + 1], rng);
    }
    if (structural) {
      epsb(nIs) *= (1 - estimate(0));
    }
    arma::vec talpha(Xbeta + epsb);
    
    arma::vec yy(y);
    int t(0);
    computeBR: ++t;
    
    // Compute Qtauy
    arma::mat Qtauy = fQtauy(yy, G, d, igroup, group, groupidx, nvec, stau, 
                             ngroup, n, ntau, type, 1);
    
    // New y
    arma::vec yyst = talpha + Qtauy * estimate.subvec(structural, ntau - 1 + structural);
    
    // check convergence
    double dist = max(abs(yyst - yy)/(abs(yy) + 1e-50));
    yy          = yyst;
    if (dist > tol && t < maxit) goto computeBR;
    out        += Qtauy;
  }
#endif
return out / boot;
}

//[[Rcpp::export]]
Eigen::ArrayXXd simInstrqpeer_EIGEN(const Eigen::VectorXd& y,
                                    const Eigen::MatrixXd& qy,
                                    const Eigen::MatrixXd& X,
                                    const std::vector<Eigen::ArrayXXd>& G,
                                    const Eigen::ArrayXd& d,
                                    const Eigen::ArrayXXi& igroup,
                                    const Eigen::ArrayXi& group, // group index
                                    const Eigen::ArrayXi& groupidx, // position in the group
                                    const Eigen::VectorXd& estimate,
                                    const Eigen::ArrayXi& nIs,
                                    const Eigen::ArrayXi& nvec,
                                    const Eigen::ArrayXd& stau,
                                    const int& boot,
                                    const bool& structural,
                                    const bool& fixedeffects,
                                    const unsigned int& nthreads,
                                    const unsigned int& seed,
                                    const int& type = 7,
                                    const double& tol = 1e-10,
                                    const int& maxit  = 500) {
  int n(y.size()), ntau(stau.size()), Kx(X.cols()), ngroup(nvec.size());
  
  // residuals
  Eigen::ArrayXd eps(y - qy * estimate.segment(structural, ntau));
  
  // remove X*beta
  Eigen::ArrayXd Xbeta(X * estimate.segment(structural + ntau, Kx));
  if (structural){
    Xbeta(nIs) *= (1 - estimate(0)); 
  }
  eps          -= Xbeta;
  if (structural){
    eps(nIs)   /= (1 - estimate(0));
  }
  
  // strata setup
  std::vector<int> strata;
  strata.push_back(0);
  if (fixedeffects) {
    for (int m = 0; m < ngroup; ++m) {
      strata.push_back(strata.back() + nvec(m));
    }
  } else {
    strata.push_back(n);
  }
  int nstrata = strata.size() - 1;
  
  // output 
  Eigen::MatrixXd out;
  
  // bootstrap
#if defined(_OPENMP)
  
  // list of instruments
  std::vector<Eigen::MatrixXd> listQtauy(nthreads);
  
  omp_set_num_threads(nthreads);
#pragma omp parallel
{
  unsigned int tid = omp_get_thread_num();
  // Initialize RNG with seed_seq for thread safety
  std::vector<unsigned int> seq_data = {seed, tid};
  std::seed_seq seq(seq_data.begin(), seq_data.end());
  std::mt19937 rng(seq);
  
  listQtauy[tid].resize(n, ntau);           
  listQtauy[tid].setZero();
  
  Eigen::ArrayXd epsb, talpha, yy, yyst;
  Eigen::MatrixXd Qtauy;
  int t;
  double dist;
  
#pragma omp for
  for (int b = 0; b < boot; ++b) {
    epsb = eps;
    for (int s = 0; s < nstrata; ++s) {
      std::shuffle(epsb.begin() + strata[s], epsb.begin() + strata[s + 1], rng);
    }
    if (structural){
      epsb(nIs) *= (1 - estimate(0));
    }
    talpha     = Xbeta + epsb;
    yy         = y;
    t          = 0;
    computeBR: ++t;
    
    // Compute Qtauy
    Qtauy = fQtauy_EIGEN(yy, G, d, igroup, group, groupidx, nvec, stau, 
                         ngroup, n, ntau, type, 1);
    
    // New y
    yyst   = talpha.matrix() + Qtauy * estimate.segment(structural, ntau);
    
    // check convergence
    dist   = ((yyst - yy).abs()/(yy.abs() + 1e-50)).maxCoeff();
    yy              = yyst;
    if (dist > tol && t < maxit) goto computeBR;
    listQtauy[tid] += Qtauy;
  }
}

  out = listQtauy[0];
  for (unsigned int k = 1; k < nthreads; ++k){
    out += listQtauy[k];
  }

#else
  // Initialize RNG with seed_seq for thread safety
  std::mt19937 rng(seed);
  
  out.resize(n, ntau);           
  out.setZero();
  
  Eigen::ArrayXd epsb, talpha, yy, yyst;
  Eigen::MatrixXd Qtauy;
  int t;
  double dist;
  
  for (int b = 0; b < boot; ++b) {
    epsb = eps;
    for (int s = 0; s < nstrata; ++s) {
      std::shuffle(epsb.begin() + strata[s], epsb.begin() + strata[s + 1], rng);
    }
    if (structural){
      epsb(nIs) *= (1 - estimate(0));
    }
    talpha     = Xbeta + epsb;
    yy         = y;
    t          = 0;
    computeBR: ++t;
    
    // Compute Qtauy
    Qtauy = fQtauy_EIGEN(yy, G, d, igroup, group, groupidx, nvec, stau, 
                         ngroup, n, ntau, type, 1);
    
    // New y
    yyst   = talpha.matrix() + Qtauy * estimate.segment(structural, ntau);
    
    // check convergence
    dist   = ((yyst - yy).abs()/(yy.abs() + 1e-50)).maxCoeff();
    yy              = yyst;
    if (dist > tol && t < maxit) goto computeBR;
    out += Qtauy;
  }
#endif
return out / boot;
}


// Equilibrium of the standard LIM model
//[[Rcpp::export]]
void fylim_ARMA(arma::vec& y,
           arma::vec& Gy,
           const std::vector<arma::mat>& G,
           const arma::vec& talpha,
           const arma::uvec& igroup,
           const int& ngroup,
           const double& lambda,
           const int& nthreads) {
  //loop over group
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel
{
  arma::mat Am;
  arma::vec ym;
  
#pragma omp for
  for (int m = 0; m < ngroup; ++m) {
    Am = -lambda * G[m];
    Am.diag() += 1.0;
    ym = arma::solve(Am, talpha.subvec(igroup(m), igroup(m + 1) - 1));
    y.rows(igroup(m), igroup(m + 1) - 1)  = ym;
    Gy.rows(igroup(m), igroup(m + 1) - 1) = G[m] * ym;
  }
}
#else
  arma::mat Am;
  arma::vec ym;

  for (int m = 0; m < ngroup; ++m) {
    Am = -lambda * G[m];
    Am.diag() += 1.0;
    ym = arma::solve(Am, talpha.subvec(igroup(m), igroup(m + 1) - 1));
    y.rows(igroup(m), igroup(m + 1) - 1)  = ym;
    Gy.rows(igroup(m), igroup(m + 1) - 1) = G[m] * ym;
  }
#endif
}

//[[Rcpp::export]]
Rcpp::List fylim(const std::vector<Eigen::MatrixXd>& G,
                       const Eigen::VectorXd& talpha,
                       const Eigen::ArrayXi& igroup,
                       const Eigen::ArrayXi& nvec,
                       const int& ngroup,
                       const double& lambda,
                       const int& n,
                       const int& nthreads) {
  Eigen::VectorXd y(n), Gy(n);
  //loop over group
#if defined(_OPENMP)
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for (int m = 0; m < ngroup; ++ m) {
    int nm(nvec(m));
    Eigen::MatrixXd Am(-lambda * G[m]);
    Am.diagonal().array() += 1;
    y.segment(igroup(m), nm)  = Am.colPivHouseholderQr().solve(talpha.segment(igroup(m), nm));
    Gy.segment(igroup(m), nm) = G[m] * y.segment(igroup(m), nm);
  }
#else
  for (int m = 0; m < ngroup; ++ m) {
    int nm(nvec(m));
    Eigen::MatrixXd Am(-lambda * G[m]);
    Am.diagonal().array() += 1;
    y.segment(igroup(m), nm)  = Am.colPivHouseholderQr().solve(talpha.segment(igroup(m), nm));
    Gy.segment(igroup(m), nm) = G[m] * y.segment(igroup(m), nm);
  }
#endif
  return Rcpp::List::create(Rcpp::_["y"] = y, Rcpp::_["Gy"] = Gy);
}

// //  This function compute optimal instruments (reduced form model)
// //[[Rcpp::export]]
// arma::mat optins_red(const arma::vec& beta,
//                      const arma::vec& y,
//                      const std::vector<arma::mat>& G,
//                      const arma::mat& X,
//                      const arma::vec& d,
//                      const arma::mat& igroup,
//                      const arma::vec& nvec,
//                      const arma::vec& stau,
//                      const int& ngroup,
//                      const int& n,
//                      const int& ntau,
//                      const int& type,
//                      const int& Kx,
//                      const double& tol,
//                      const int& maxit) {
//   arma::vec talpha(X*beta.tail(Kx));
//   arma::vec Ey = y;
//   fNashE(Ey, G, d, talpha, beta.head(ntau), igroup, nvec, stau, ngroup, n, 
//          ntau, type, tol, maxit);
//   return fQtauy(Ey, G, d, igroup, nvec, stau, ngroup, n, ntau, type);
// }
// 
// //  This function compute optimal instruments (reduced form model)
// //[[Rcpp::export]]
// arma::mat optins_struc(const arma::vec& beta,
//                        const arma::vec& y,
//                        const std::vector<arma::mat>& G,
//                        const arma::mat& X,
//                        const arma::vec& d,
//                        const arma::mat& igroup,
//                        const arma::vec& nvec,
//                        const arma::vec& stau,
//                        const arma::uvec& nIs,
//                        const int& ngroup,
//                        const int& n,
//                        const int& ntau,
//                        const int& type,
//                        const int& Kx,
//                        const double& tol,
//                        const int& maxit) {
//   arma::vec talpha(X*beta.tail(Kx)); talpha.elem(nIs) *= (1 - beta(0));
//   arma::vec Ey = y;
//   fNashE(Ey, G, d, talpha, beta.subvec(1, ntau), igroup, nvec, stau, ngroup, n, 
//          ntau, type, tol, maxit);
//   return fQtauy(Ey, G, d, igroup, nvec, stau, ngroup, n, ntau, type);
// }