// [[Rcpp::depends(RcppArmadillo)]]

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


// This function removes columns to obtain full rank matrices
//[[Rcpp::export]]
arma::uvec fcheckrank(const arma::mat& X) {
  int Kx(X.n_cols), Rout(arma::rank(X));
  arma::uvec out(arma::ones<arma::uvec>(Kx));
  if (Rout < Kx) {
    for (int k(0); k < Kx; ++ k) {
      arma::uvec outst(out); 
      outst(k) = 0;
      arma::mat X1(X.cols(arma::find(out == 1)));
      arma::mat X2(X.cols(arma::find(outst == 1)));
      int Routst(arma::rank(X2));
      
      int drank1 = X1.n_cols - Rout;
      int drank2 = X2.n_cols - Routst;
      if (drank2 < drank1){
        out(k) = 0;
        Rout   = Routst;
      } 
    }
  }
  return out;
}

// demean
//[[Rcpp::export]]
arma::mat demean(arma::mat X,
                 const arma::mat& igroup,
                 const int & ngroup) {
  for (int r(0); r < ngroup; ++ r) {
    int n1(igroup(r, 0)), n2(igroup(r, 1));
    X.rows(n1, n2).each_row() -= arma::mean(X.rows(n1, n2), 0);
  }
  return X;
}

// demean for the structural model
//[[Rcpp::export]]
arma::mat demean_separate(arma::mat X,
                          const arma::mat& igroup,
                          const arma::uvec& Is,
                          const int & ngroup,
                          const int& n) {
  arma::uvec isnIs = arma::ones<arma::uvec>(n);
  isnIs.elem(Is).zeros();
  for (int r(0); r < ngroup; ++ r) {
    int n1(igroup(r, 0)), n2(igroup(r, 1));
    arma::mat Xm(X.rows(n1, n2));
    arma::uvec tp(isnIs.subvec(n1, n2));
    // For isolated
    arma::uvec itp(arma::find(tp == 0));
    if (itp.n_elem > 0) {
      arma::mat Xm1 = Xm.rows(itp);
      Xm1.each_row() -= arma::mean(Xm1, 0);
      Xm.rows(itp) = Xm1;
    }
    // For non-isolated
    itp = arma::find(tp != 0);
    if (itp.n_elem > 0) {
      arma::mat Xm1 = Xm.rows(itp);
      Xm1.each_row() -= arma::mean(Xm1, 0);
      Xm.rows(itp) = Xm1;
    }
    X.rows(n1, n2) = Xm;
  }
  return X;
}



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
