// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;
using namespace std;

// Assigning folds to groups
//[[Rcpp::export]]
std::vector<Eigen::ArrayXi> fassignfold(const Eigen::ArrayXi& nvec,
                                        const int& nfold) {
  int ngroup(nvec.size());
  Eigen::ArrayXi ncs(ngroup + 1); ncs(0) = 0;
  
  // Fold for each group
  int fold;
  std::vector<Eigen::ArrayXi> lfold(nfold); 
  Eigen::ArrayXi foldsize(Eigen::ArrayXi::Zero(nfold));
  for (int s(0); s < ngroup; ++s) {
    ncs(s + 1) = ncs(s) + nvec(s);
    int minsize(foldsize.minCoeff());
    for (int k(0); k < nfold; ++k) {
      if (foldsize(k) == minsize) {
        fold         = k;
        foldsize(k) += nvec(s);
        break;
      }
    }
    
    lfold[fold].conservativeResize(lfold[fold].size() + nvec(s));
    lfold[fold].tail(nvec(s)) = Eigen::ArrayXi::LinSpaced(nvec(s), ncs(s), ncs(s + 1) - 1); ;
  }
  return lfold;
}



// The first stage to compute instruments
//[[Rcpp::export]]
Rcpp::List  fpredinst(const Eigen::MatrixXd& endo,
                      const Eigen::MatrixXd& ins,
                      const std::vector<Eigen::ArrayXi>& fold,
                      const int& nthreads) {
  int nfold(fold.size()), n(endo.rows()), nendo(endo.cols());
  
  // complement of the fold
  std::vector<Eigen::ArrayXi> cfold(nfold);
  for (int k = 0; k < nfold; ++k) {
    for (int l = 0; l < nfold; ++l) {
      if (k != l) {
        int nl(fold[l].size());
        cfold[k].conservativeResize(cfold[k].size() + nl);
        cfold[k].tail(nl) = fold[l];
      }
    }
  }
  
  // Predict instruments and residuals
  Eigen::MatrixXd pins(n, nendo);
  
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for (int k = 0; k < nfold; ++k) {
    Eigen::MatrixXd ZZ = ins(cfold[k], Eigen::all).transpose() * ins(cfold[k], Eigen::all);
    Eigen::MatrixXd ZV = ins(cfold[k], Eigen::all).transpose() * endo(cfold[k], Eigen::all);
    Eigen::MatrixXd pi = ZZ.householderQr().solve(ZV);
    pins(fold[k], Eigen::all) = ins(fold[k], Eigen::all) * pi;
  }
#else
  for (int k = 0; k < nfold; ++k) {
    Eigen::MatrixXd ZZ = ins(cfold[k], Eigen::all).transpose() * ins(cfold[k], Eigen::all);
    Eigen::MatrixXd ZV = ins(cfold[k], Eigen::all).transpose() * endo(cfold[k], Eigen::all);
    Eigen::MatrixXd pi = ZZ.householderQr().solve(ZV);;
    pins(fold[k], Eigen::all) = ins(fold[k], Eigen::all) * pi;
  }
#endif
  Eigen::MatrixXd u = endo - pins;
  
  return Rcpp::List::create(Rcpp::_["ins"] = pins, Rcpp::_["resid"] = u);
}
