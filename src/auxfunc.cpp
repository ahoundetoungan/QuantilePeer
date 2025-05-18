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
                 const arma::mat& igroup,
                 const int & ngroup) {
  for (int r(0); r < ngroup; ++ r) {
    int n1(igroup(r, 0)), n2(igroup(r, 1));
    X.rows(n1, n2).each_row() -= arma::mean(X.rows(n1, n2), 0);
  }
  return X;
}

// Demean for the structural model
//[[Rcpp::export]]
arma::mat Demean_separate(arma::mat X,
                          const arma::mat& igroup,
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
                           const arma::mat& igroup,
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
    ins = arma::join_rows(ins, xb);
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

// //[[Rcpp::export]]
// arma::uvec fcheckrankcpp(const arma::mat& X, const double& tol = 1e-10) {
//   arma::rowvec s(arma::stddev(X, 1, 0)), m(arma::mean(X, 0));
//   m.elem(arma::find(s < tol)).zeros();
//   s.elem(arma::find(s < tol)).ones();
//   arma::mat U((X.each_row() - m).each_row()/s);
//   U = U.t()*U/U.n_rows; 
//   arma::mat Q, R;
//   arma::qr(Q, R, U);
//   return arma::find(arma::abs(R.diag()) > tol);
// }
// 
// 
// arma::uvec fcheckrankcpp_this_is_slower(const arma::mat& X) {
//   arma::rowvec m(arma::mean(X, 0)), s(arma::stddev(X, 1, 0));
//   s.elem(arma::find(s == 0)).ones();
//   arma::mat U((X.each_row() - m).each_row()/s);
//   U = U.t()*U/X.n_rows;
// 
//   int Kx(X.n_cols), Rout(arma::rank(U));
//   unsigned int tp(1);
//   arma::uvec out(arma::ones<arma::uvec>(Kx));
//   if (Rout < Kx) {
//     for (int k(0); k < Kx; ++ k) {
//       arma::mat Up(U.cols(arma::find(out.head(k + 1) == 1)));
//       if (arma::rank(Up) != tp) {
//         out(k) = 0;
//       } else {
//         ++ tp;
//       }
//     }
//   }
//   return arma::find(out == 1);
// }
// 
// 
// // This function computes the matrix that will be used tfor the qr decomp
// //[[Rcpp::export]]
// arma::mat fmatforrankarma(const arma::mat& X, const double& tol = 1e-10) {
//   arma::rowvec s(arma::stddev(X, 1, 0)), m(arma::mean(X, 0));
//   m.elem(arma::find(s < tol)).zeros();
//   s.elem(arma::find(s < tol)).ones();
//   arma::mat U((X.each_row() - m).each_row()/s);
//   return U.t()*U/U.n_rows;
// }
// 
// //[[Rcpp::export]]
// Eigen::MatrixXd fmatforrank(const Eigen::MatrixXd& X, const double& tol = 1e-10) {
//   int n(X.rows());
//   Eigen::RowVectorXd m(X.colwise().mean());
//   Eigen::RowVectorXd s(((X.rowwise() - m).array().square().colwise().sum() / n).sqrt());
//   m = (s.array() < tol).select(0, m);
//   s = (s.array() < tol).select(1, s);
//   Eigen::MatrixXd U((X.rowwise() - m).array().rowwise() / s.array());
//   return U.transpose()*U/n;
// }
// 
// 
// //[[Rcpp::export]]
// Eigen::MatrixXd Demean(Eigen::MatrixXd& X,
//                        const Eigen::MatrixXi& igroup,
//                        const int & ngroup) {
//   for (int r(0); r < ngroup; ++ r) {
//     int n1(igroup(r, 0)), n2(igroup(r, 1));
//     X(Eigen::seq(n1, n2), Eigen::all).rowwise() -= X(Eigen::seq(n1, n2), Eigen::all).colwise().mean();
//   }
//   return X;
// }

// // Demean for the structural model
// //[[Rcpp::export]]
// Eigen::MatrixXd Demean_separate(Eigen::MatrixXd& X,
//                                 const Eigen::MatrixXi& igroup,
//                                 const LVectorXi& LIs,
//                                 const LVectorXi& LnIs,
//                                 const int & ngroup,
//                                 const int& n) {
//   for (int r(0); r < ngroup; ++ r) {
//     Eigen::VectorXi Isr(LIs[r]), nIsr(LnIs[r]);
//     // For isolated
//     if (Isr.size() > 0) {
//       X(Isr, Eigen::all).rowwise() -= X(Isr, Eigen::all).colwise().mean();
//     }
//     // For non-isolated
//     if (nIsr.size() > 0) {
//       X(nIsr, Eigen::all).rowwise() -= X(nIsr, Eigen::all).colwise().mean();
//     }
//   }
//   return X;
// }

// // from unconstrained parameters to constrained parameters
// // The first lambda is the sum of lambda2
// //[[Rcpp::export]]
// arma::vec flambda(const arma::vec& lambdatilde,
//                   const arma::vec& linf,
//                   const arma::vec& lsup,
//                   const int& ntau) {
//   arma::vec out(ntau + 1, arma::fill::zeros);
//   out(0) = (lsup(0) - linf(0))/(1 + exp(-lambdatilde(0))) + linf(0);
//   double sl(0);
//   for (int k(0); k < ntau; ++ k) {
//     double ak(linf(k + 1)*(1 - sl)), bk(lsup(k + 1)*(1 - sl));
//     out(k + 1) = (bk - ak)/(1 + exp(-lambdatilde(k + 1))) + ak;
//     sl        += abs(out(k + 1));
//   }
//   return out;
// }
// 
// // from constrained parameters to unconstrained parameters
// //[[Rcpp::export]]
// arma::vec flambdatilde(const arma::vec& lambda,
//                        const arma::vec& linf,
//                        const arma::vec& lsup,
//                        const int& ntau) {
//   arma::vec out(ntau + 1, arma::fill::zeros);
//   out(0) = log(lambda(0) - linf(0)) - log(lsup(0) - lambda(0));
//   double sl(1);
//   for (int k(0); k < ntau; ++ k) {
//     double ak(linf(k + 1)*sl), bk(lsup(k + 1)*sl);
//     out(k + 1) = log(lambda(k + 1) - ak) - log(bk - lambda(k + 1));
//     sl        -= abs(lambda(k + 1));
//   }
//   return out;
// }
// 
// // derivative flambda
// arma::mat fdlambda(const arma::vec& lambdatilde,
//                    const arma::vec& lambda,
//                    const arma::rowvec& linf,
//                    const arma::rowvec& lsup,
//                    const int& ntau) {
//   arma::mat out(ntau + 1, ntau + 1, arma::fill::zeros);
//   out(0, 0) = (lsup(0) - linf(0))*exp(lambdatilde(0))/pow(1 + exp(lambdatilde(0)), 2);
//   arma::rowvec da(arma::zeros<arma::rowvec>(ntau));
//   double sl(0);
//   for (int k(0); k < ntau; ++ k) {
//     double tp(exp(-lambdatilde(k + 1)));
//     double ak(linf(k + 1)*(1 - sl)), bk(lsup(k + 1)*(1 - sl));
//     out(k + 1, k + 1)          = (bk - ak)*tp/pow(1 + tp, 2);
//     
//     arma::rowvec dak(linf(k + 1)*da), dbk(lsup(k + 1)*da);
//     out.row(k + 1).tail(ntau) += (dak - dbk)/(1 + tp) - dak;
//     da                        += (out.row(k + 1).tail(ntau))*(linf(k + 1) + lsup(k + 1));
//     sl                        += abs(lambda(k + 1));
//   }
//   return out;
// }
// 
// // moment function g
// //[[Rcpp::export]]
// arma::mat g(const arma::vec& theta, 
//             Rcpp::List& x){
//   arma::vec y = x["y"];
//   arma::mat qy = x["qy"];
//   arma::mat X = x["X"];
//   arma::mat ins = x["ins"];
//   arma::uvec Is = x["Is"];
//   arma::uvec nIs = x["nIs"];
//   arma::vec linf = x["linf"];
//   arma::vec lsup = x["lsup"];
//   int n = x["n"];
//   int Kx = x["Kx"];
//   int ntau = x["ntau"];
//   
//   arma::vec lambda(flambda(theta.head(ntau + 1), linf, lsup, ntau));
//   arma::vec xb(X*theta.tail(Kx));
//   arma::vec e(n);
//   // Isolated
//   e.elem(Is)  = y.elem(Is) - xb.elem(Is);
//   // Nonisolated
//   e.elem(nIs) = y.elem(nIs) - qy.rows(nIs)*lambda.tail(ntau) - xb.elem(nIs)*(1 - lambda(0));
//   return ins.each_col()%e;
// }
// 
// // Jocobian of the moment function g
// //[[Rcpp::export]]
// arma::mat dg(const arma::vec& theta, 
//              Rcpp::List& x){
//   arma::vec y = x["y"];
//   arma::mat qy = x["qy"];
//   arma::mat X = x["X"];
//   arma::mat ins = x["ins"];
//   arma::uvec Is = x["Is"];
//   arma::uvec nIs = x["nIs"];
//   arma::vec linf = x["linf"];
//   arma::vec lsup = x["lsup"];
//   int n = x["n"];
//   int Kx = x["Kx"];
//   int ntau = x["ntau"];
//   int Kins = x["Kins"];
//   arma::vec lambda(flambda(theta.head(ntau + 1), linf, lsup, ntau));
//   arma::vec xb(X*theta.tail(Kx));
//   arma::mat out(Kins, 1 + ntau + Kx);
//   
//   // Derivative with respect to theta0
//   out.col(0) = -arma::trans(ins.rows(nIs))*xb.elem(nIs);
//   // Derivative with respect to theta1 to theta(ntau)
//   out.cols(1, ntau) = -arma::trans(arma::sum(qy.rows(nIs), 0));
//   // Derivative with respect to beta
//   out.tail_cols(Kx) = -arma::trans(arma::sum(X.rows(Is), 0) + arma::sum(X.rows(nIs), 0)*(1 - lambda(0)));
//   return out.t()/n;
// }

// //[[Rcpp::export]]
// Rcpp::List fdatadiagnostic(arma::vec& y,
//                            arma::mat& endo,
//                            arma::mat& X,
//                            arma::mat& ins,
//                            const arma::vec& theta,
//                            const arma::uvec& idX1,
//                            const arma::uvec& idX2,
//                            const arma::mat& igroup,
//                            const arma::uvec& Is,
//                            const arma::uvec& nIs,
//                            const int& n,
//                            const int & ngroup,
//                            const int& ntau,
//                            const bool& struc,
//                            const std::string& FE) {
//   if (struc) {
//     arma::vec xb(X.cols(idX1)*theta.elem(idX1 + ntau + 1));
//     if (idX2.n_elem == 0) {
//       X = xb;
//     } else {
//       X = arma::join_rows(xb, X.cols(idX2));
//     }
//     ins = arma::join_rows(xb, ins);
//   }
//   
//   if (FE == "join") {
//     y    = Demean(y, igroup, ngroup);
//     endo = Demean(endo, igroup, ngroup);
//     X    = Demean(X, igroup, ngroup);
//     ins  = Demean(ins, igroup, ngroup);
//   } else if(FE == "separate") {
//     y    = Demean_separate(y, igroup, Is, ngroup, n);
//     endo = Demean_separate(endo, igroup, Is, ngroup, n);
//     X    = Demean_separate(X, igroup, Is, ngroup, n);
//     ins  = Demean_separate(ins, igroup, Is, ngroup, n);
//   }
//   
//   arma::uvec nvecnIs(ngroup);
//   if (struc) {
//     y      = y.elem(nIs);
//     endo   = endo.rows(nIs);
//     X      = X.rows(nIs);
//     ins    = ins.rows(nIs);
//     arma::uvec isnIs(arma::ones<arma::uvec>(n));
//     isnIs.elem(Is).zeros();
//     for (int r(0); r < ngroup; ++ r) {
//       int n1(igroup(r, 0)), n2(igroup(r, 1));
//       nvecnIs(r) = sum(isnIs.subvec(n1, n2));
//     }
//   }
//   
//   return Rcpp::List::create(_["nvecnIs"] = nvecnIs, _["y"] = y, _["endo"] = endo, _["X"] = X, _["ins"] = ins);
// }
