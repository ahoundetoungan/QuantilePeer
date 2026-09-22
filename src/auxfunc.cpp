// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
//#define NDEBUG
// #include <RcppNumerical.h>
#include <RcppEigen.h>
// typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;
// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;

// using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;


// This function removes columns to obtain full rank matrices
//[[Rcpp::export]]
Eigen::Array<bool, Eigen::Dynamic, 1> fcheckrankEigen(const Eigen::MatrixXd& X, const double& tol = 1e-10) {
  int n(X.rows());
  Eigen::RowVectorXd m(X.colwise().mean());
  Eigen::RowVectorXd s(((X.rowwise() - m).array().square().colwise().sum() / n).sqrt());
  m = (s.array() < tol).select(0, m);
  s = (s.array() < tol).select(1, s);
  Eigen::MatrixXd U((X.rowwise() - m).array().rowwise() / s.array());
  U = U.transpose()*U/n;
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(U);
  Eigen::MatrixXd R(qr.matrixQR().topRows(U.cols()));
  // std::cout<<R.diagonal().transpose()<<std::endl;
  return R.diagonal().array().abs() > tol;
}

// Demean
//[[Rcpp::export]]
arma::mat Demean(arma::mat X,
                 const arma::uvec& igroup,
                 const int & ngroup) {
  for (int r(0); r < ngroup; ++ r) {
    int n1(igroup(r)), n2(igroup(r + 1) - 1);
    X.rows(n1, n2).each_row() -= arma::mean(X.rows(n1, n2), 0);
  }
  return X;
}

// Demean for the structural model
//[[Rcpp::export]]
arma::mat Demean_separate(arma::mat X,
                          const arma::uvec& igroup,
                          const Rcpp::List& LIs,
                          const Rcpp::List& LnIs,
                          const int & ngroup,
                          const int& n) {
  for (int r(0); r < ngroup; ++ r) {
    arma::uvec Isr = LIs[r], nIsr = LnIs[r];
    // For isolated
    if (Isr.n_elem > 0) {
      arma::mat Xr(X.rows(Isr));
      X.rows(Isr) = Xr.each_row() - arma::mean(Xr, 0);
    }
    // For non-isolated
    if (nIsr.n_elem > 0) {
      arma::mat Xr(X.rows(nIsr));
      X.rows(nIsr) = Xr.each_row() - arma::mean(Xr, 0);
    }
  }
  return X;
}


// data for the diagnostic function
//[[Rcpp::export]]
Rcpp::List fdatadiagnostic(arma::vec& y,
                           arma::mat& endo,
                           arma::mat& X,
                           arma::mat& ins,
                           const arma::vec& theta,
                           const arma::uvec& idX1,
                           const arma::uvec& idX2,
                           const arma::uvec& igroup,
                           const arma::uvec& nIs,
                           const Rcpp::List& LIs,
                           const Rcpp::List& LnIs,
                           const int& n,
                           const int & ngroup,
                           const int& ntau,
                           const bool& struc,
                           const std::string& FE) {
  if (struc) {
    arma::vec xb(X.cols(idX1)*theta.elem(idX1 + ntau + 1));
    X   = arma::join_rows(X.cols(idX2), xb);
  }
  
  if (FE == "join") {
    y    = Demean(y, igroup, ngroup);
    endo = Demean(endo, igroup, ngroup);
    X    = Demean(X, igroup, ngroup);
    ins  = Demean(ins, igroup, ngroup);
  } else if(FE == "separate") {
    y    = Demean_separate(y, igroup, LIs, LnIs, ngroup, n);
    endo = Demean_separate(endo, igroup, LIs, LnIs, ngroup, n);
    X    = Demean_separate(X, igroup, LIs, LnIs, ngroup, n);
    ins  = Demean_separate(ins, igroup, LIs, LnIs, ngroup, n);
  }
  
  arma::uvec nvecnIs(ngroup);
  if (struc) {
    y      = y.elem(nIs);
    endo   = endo.rows(nIs);
    X      = X.rows(nIs);
    ins    = ins.rows(nIs);
  }
  
  return Rcpp::List::create(_["y"] = y, _["endo"] = endo, _["X"] = X, _["ins"] = ins);
}


// Assigning folds to groups
// This function assigning folds to subnetworks
//[[Rcpp::export]]
Eigen::ArrayXi fassignfold(const Eigen::ArrayXi& subnetwork,
                           const int& nfold,
                           const unsigned long long& seed) {
  
  // sizes
  int n = subnetwork.size();
  std::unordered_set<int> usubnetwork(subnetwork.data(), subnetwork.data() + n);
  int R(usubnetwork.size());
  if (R < 2) {
    Rcpp::stop("There is only one subnetwork.");
  }
  
  // group should take values from 0 to ngroup - 1 with possible duplication
  if ((subnetwork.minCoeff() != 0) || (subnetwork.maxCoeff() != (R - 1))) {
    Rcpp::stop("Subnetwork should take values from 0, 1, 2, ... without jumps.");
  }
  
  // Number of data per subnetwork
  Eigen::ArrayXi nvec(Eigen::ArrayXi::Zero(R));
  for (int i(0); i < n; ++i) {
    nvec(subnetwork(i)) += 1;
  }
  
  // Shuffle order
  std::mt19937 rng(seed);
  std::vector<int> netorder(R);
  std::iota(netorder.begin(), netorder.end(), 0);
  std::shuffle(netorder.begin(), netorder.end(), rng);
  
  // Fold for each subnetwork
  Eigen::ArrayXi fold(R);
  Eigen::ArrayXi foldsize(Eigen::ArrayXi::Zero(nfold));
  for (int r : netorder) {
    int minsize(foldsize.minCoeff());
    for (int k(0); k < nfold; ++k) {
      if (foldsize(k) == minsize) {
        fold(r)      = k;
        foldsize(k) += nvec(r);
        break;
      }
    }
  }
  return fold(subnetwork);
}

