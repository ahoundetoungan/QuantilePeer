#include <RcppArmadillo.h>
//#define NDEBUG
// #include <RcppNumerical.h>
// #include <RcppEigen.h>
// 
// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;

// using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// This function implements the GMM estimator (reduced-form)
//[[Rcpp::export]]
Rcpp::List fgmm_red(const arma::vec& y,
                    const arma::mat& V,
                    const arma::mat& ins,
                    arma::mat& W,
                    const arma::mat& igroup,
                    const int& ngroup,
                    const int& Kx,
                    const int& Kins, 
                    const int& ntau,
                    const int& n,
                    const int& HAC = 0,
                    const bool& iv = true){
  arma::mat ZV(ins.t()*V), Zy(ins.t()*y);
  if (iv) {
    W = arma::inv(ins.t()*ins/n);
  }
  arma::mat VZW(ZV.t()*W), VZWZV(VZW*ZV);
  
  // estimate
  arma::vec parms(arma::solve(VZWZV, VZW*Zy));
  
  // variance
  arma::vec e(y - V*parms);
  arma::mat VZe(Kins, Kins, arma::fill::zeros);
  if (HAC == 0) {
    VZe    = (sum(e%e)/(n - Kx - ntau))*(ins.t()*ins);
  }
  if (HAC == 1) {
    arma::mat Ze(ins.each_col()%e);
    VZe    = Ze.t()*Ze;
  }
  if (HAC == 2) {
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), n2(igroup(r, 1));
      arma::vec tp(ins.rows(n1, n2).t()*e.subvec(n1, n2));
      VZe += tp*tp.t();
    }
  }
  arma::mat Vpa(arma::solve(VZWZV, VZW*VZe)); // inv(V'ZWZ'V)*V'ZW*Var(Ze)
  Vpa = arma::solve(VZWZV, VZW*Vpa.t());//inv(V'ZWZ'V)*V'ZW*Var(Ze)*WZ'V*inv(V'ZWZ'V)
  
  // overidentification
  arma::vec Ze(ins.t()*e);
  double stat = sum(Ze.t()*arma::solve(VZe, Ze));
  
  // // OLS
  // arma::mat VV(V.t()*V);
  // arma::vec OLS(arma::solve(VV, V.t()*y)), eOLS = y - V*OLS;
  // arma::mat VVe(Kx, Kx, arma::fill::zeros);
  // if (cluster) {
  //   for (int r(0); r < ngroup; ++ r) {
  //     int n1(igroup(r, 0)), n2(igroup(r, 1));
  //     arma::vec tp(V.rows(n1, n2).t()*eOLS.subvec(n1, n2));
  //     VVe += tp*tp.t();
  //   }
  // } else {
  //   arma::mat Ve(V.each_col()%eOLS);
  //   VVe    = Ve.t()*Ve;
  // } 
  // arma::mat VOLS(arma::solve(VV, VVe)); //inv(V'V)*Var(V'e)
  // VOLS = arma::solve(VV, VOLS.t());//inv(V'V)*Var(V'e)*inv(V'V)
  // 
  // // Hausman Wu
  // arma::vec dparms = parms - OLS;
  // double haus = sum(dparms.t()*arma::solve(Vpa - VOLS, dparms));
  return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, _["VZe"] = VZe, 
                            _["Overident"] = stat, _["df"] = Kins - Kx - ntau);
}

  
// This function implements the GMM estimator (Structural model)
//[[Rcpp::export]]
Rcpp::List fgmm_struc(const arma::vec& y,
                       const arma::mat& X,
                       const arma::mat& qy,
                       const arma::mat& ins,
                       arma::mat& W1,
                       arma::mat& W2,
                       const arma::mat& igroup,
                       const arma::uvec& nIs,
                       const arma::uvec& Is,
                       const int& ngroup,
                       const int& Kins,
                       const int& Kx,
                       const int& ntau,
                       const int& n,
                       const int& HAC = 0,
                       const bool& iv = true){
  int n_iso(Is.n_elem), n_niso(n - n_iso);
  // First stage
  arma::vec y1(y.elem(Is));
  arma::mat X1(X.rows(Is)), XX1(X1.t()*X1);
  if (iv) {
    W1 = arma::inv(XX1/n_iso);
  }
  arma::mat XXW1(XX1*W1), XXWXX1(XXW1*XX1);
  arma::vec b(arma::solve(XXWXX1, XXW1*X1.t()*y1));
  
  // Second stage
  arma::vec Xb(X*b), Xb1(Xb.elem(Is)), Xb2(Xb.elem(nIs)), y2(y.elem(nIs));
  arma::mat X2(X.rows(nIs)), V2(arma::join_rows(qy.rows(nIs), Xb2)), 
  Z2(arma::join_rows(Xb2, ins.rows(nIs))), ZV2(Z2.t()*V2), ZZ2(Z2.t()*Z2); 
  if (iv) {
    W2 = arma::inv(ZZ2/n_niso);
  }
  arma::mat VZW2(ZV2.t()*W2), VZWZV2(VZW2*ZV2);
  arma::vec Zy2(Z2.t()*y2);
  arma::vec lambda(arma::solve(VZWZV2, VZW2*Zy2));
  
  // Variance
  arma::vec  e1(y1 - Xb1), e2(y2 - V2*lambda);
  arma::mat H(Kx + ntau + 1, Kx + Kins + 1, arma::fill::zeros);
  H.submat(0, 0, Kx - 1, Kx - 1) = XXW1;
  H.submat(Kx, Kx, Kx + ntau, Kx + Kins) = VZW2;
  
  arma::mat dF(Kx + Kins + 1, Kx + ntau + 1, arma::fill::zeros);
  dF.submat(0, 0, Kx - 1, Kx - 1) = -XX1;
  dF.submat(Kx, 0, Kx, Kx - 1) = (e2.t()*X2);
  dF.submat(Kx, 0, Kx + Kins, Kx - 1) -= (Z2.t()*X2*lambda(ntau));
  dF.submat(Kx, Kx, Kx + Kins, Kx + ntau) = -ZV2;
  
  arma::mat VF(Kx + Kins + 1, Kx + Kins + 1, arma::fill::zeros);
  if (HAC == 0) {
    VF.submat(0, 0, Kx - 1, Kx - 1) = (sum(e1%e1)/(n_iso - Kx))*XX1;
    VF.submat(Kx, Kx, Kx + Kins, Kx + Kins) = (sum(e2%e2)/(n_niso - ntau - 1))*ZZ2;
  }
  if (HAC == 1) {
    arma::mat Xe1(X1.each_col()%e1);
    arma::mat Ze2(Z2.each_col()%e2);
    VF.submat(0, 0, Kx - 1, Kx - 1) = Xe1.t()*Xe1;
    VF.submat(Kx, Kx, Kx + Kins, Kx + Kins) = Ze2.t()*Ze2;
  }
  if (HAC == 2) {
    X1   = X; X1.rows(nIs).zeros();
    e1   = y - Xb; e1.elem(nIs).zeros();
    Z2   = arma::join_rows(Xb, ins); Z2.rows(Is).zeros();
    arma::vec Vl(arma::join_rows(qy, Xb)*lambda);
    e2   = y - Vl; e2.elem(Is).zeros();
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), n2(igroup(r, 1));
      arma::vec tp(arma::join_cols(X1.rows(n1, n2).t()*e1.subvec(n1, n2),
                                   Z2.rows(n1, n2).t()*e2.subvec(n1, n2)));
      VF += tp*tp.t();
    }
  }
  arma::mat HdF(H*dF);//H * dF
  arma::mat Vpa(arma::solve(HdF, H*VF*H.t())); // inv(H * dF) * H * Var(F) * H'
  Vpa = arma::solve(HdF, Vpa.t());//inv(H * dF) * H * Var(F) * H' * inv(H * dF)'
  // overidentification
  arma::vec F2(Z2.t()*e2);
  arma::mat VF1(VF.submat(0, 0, Kx - 1, Kx - 1)), VF2(VF.submat(Kx, Kx, Kx + Kins, Kx + Kins));
  double stat = sum(F2.t()*arma::solve(VF2, F2));
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa, 
                            _["VF1"] = VF1, _["VF2"] = VF2, _["Overident"] = stat, 
                              _["df"] = Kins - ntau);
}




        