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
                    const int& Kest, 
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
  arma::vec yhat(V*parms), e(y - yhat);
  arma::mat VZe(Kins, Kins, arma::fill::zeros);
  double s2(R_NaN);
  if (HAC == 0) {
    s2     = sum(e%e)/(n - Kest);
    VZe    = s2*(ins.t()*ins);
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
  return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, _["VZe"] = VZe, _["Overident"] = stat, 
                            _["df"] = Kins - Kx - ntau, _["yhat"] = yhat, _["sigma2"] = s2);
}



// This function implements the GMM estimator (Structural model)
//[[Rcpp::export]]
Rcpp::List fgmm_struc(const arma::vec& y,
                      const arma::mat& X,
                      const arma::mat& qy,
                      const arma::mat& ins,
                      arma::mat& W1,
                      arma::mat& W2,
                      const arma::uvec& idX1,
                      const arma::uvec& idX2,
                      const int& Kx1,
                      const int& Kx2,
                      const arma::mat& igroup,
                      const arma::uvec& nIs,
                      const arma::uvec& Is,
                      const int& ngroup,
                      const int& Kins,
                      const int& Kx,
                      const int& ntau,
                      const int& n,
                      const int& Kest,
                      const int& HAC = 0,
                      const bool& iv = true){
  int n_iso(Is.n_elem), n_niso(n - n_iso);
  // First stage
  arma::vec y1(y.elem(Is));
  arma::mat X1(X.rows(Is)); X1 = X1.cols(idX1);
  arma::mat XX1(X1.t()*X1);
  if (iv) {
    W1 = arma::inv(XX1/n_iso);
  }
  arma::mat XXW1(XX1*W1), XXWXX1(XXW1*XX1);
  arma::vec b(arma::solve(XXWXX1, XXW1*X1.t()*y1));
  
  // Second stage
  arma::vec Xb(X.cols(idX1)*b), Xb1(Xb.elem(Is)), Xb2(Xb.elem(nIs)), y2(y.elem(nIs));
  arma::mat X2(X.rows(nIs));
  arma::mat X21(X2.cols(idX1));
  X2 = X2.cols(idX2);
  arma::mat V2(arma::join_rows(arma::join_rows(Xb2, qy.rows(nIs)), X2)), 
  Z2(arma::join_rows(Xb2, ins.rows(nIs))), ZV2(Z2.t()*V2), ZZ2(Z2.t()*Z2); 
  if (iv) {
    W2 = arma::inv(ZZ2/n_niso);
  }
  arma::mat VZW2(ZV2.t()*W2), VZWZV2(VZW2*ZV2);
  arma::vec Zy2(Z2.t()*y2);
  arma::vec lambda(arma::solve(VZWZV2, VZW2*Zy2));
  
  // Variance
  arma::vec y1hat(Xb1), y2hat(V2*lambda), e1(y1 - y1hat), e2(y2 - y2hat), yhat(n);
  yhat.elem(Is) = y1hat; yhat.elem(nIs) = y2hat;
  arma::vec e = y - yhat;
  arma::mat H(Kx + ntau + 1, Kx1 + Kins + 1, arma::fill::zeros);
  H.submat(0, 0, Kx1 - 1, Kx1 - 1) = XXW1;
  H.submat(Kx1, Kx1, Kx + ntau, Kx1 + Kins) = VZW2;
  
  arma::mat dF(Kx1 + Kins + 1, Kx + ntau + 1, arma::fill::zeros);
  dF.submat(0, 0, Kx1 - 1, Kx1 - 1) = -XX1;
  dF.submat(Kx1, 0, Kx1, Kx1 - 1) = (e2.t()*X21);
  dF.submat(Kx1, 0, Kx1 + Kins, Kx1 - 1) -= (Z2.t()*X21*lambda(0));
  dF.submat(Kx1, Kx1, Kx1 + Kins, Kx + ntau) = -ZV2;
  
  arma::mat VF(Kx1 + Kins + 1, Kx1 + Kins + 1, arma::fill::zeros);
  double s2(R_NaN);
  if (HAC == 0) {
    s2     = sum(e%e)/(n - Kest);
    VF.submat(0, 0, Kx1 - 1, Kx1 - 1) = s2*XX1;
    VF.submat(Kx1, Kx1, Kx1 + Kins, Kx1 + Kins) = s2*ZZ2;
  }
  if (HAC == 1) {
    arma::mat Xe1(X1.each_col()%e1);
    arma::mat Ze2(Z2.each_col()%e2);
    VF.submat(0, 0, Kx1 - 1, Kx1 - 1) = Xe1.t()*Xe1;
    VF.submat(Kx1, Kx1, Kx1 + Kins, Kx1 + Kins) = Ze2.t()*Ze2;
  }
  if (HAC == 2) {
    X1   = X.cols(idX1); X1.rows(nIs).zeros();
    Z2   = arma::join_rows(Xb, ins); Z2.rows(Is).zeros();
    e2   = y - yhat;
    // e1.elem(nIs).zeros();
    // e2.elem(Is).zeros();
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), n2(igroup(r, 1));
      arma::vec tp(arma::join_cols(X1.rows(n1, n2).t()*e2.subvec(n1, n2),
                                   Z2.rows(n1, n2).t()*e2.subvec(n1, n2)));
      VF += tp*tp.t();
    }
  }
  arma::mat HdF(H*dF);//H * dF
  arma::mat Vpa(arma::solve(HdF, H*VF*H.t())); // inv(H * dF) * H * Var(F) * H'
  Vpa = arma::solve(HdF, Vpa.t());//inv(H * dF) * H * Var(F) * H' * inv(H * dF)'
  // overidentification
  arma::vec F2(Z2.t()*e2);
  arma::mat VF1(VF.submat(0, 0, Kx1 - 1, Kx1 - 1)), VF2(VF.submat(Kx1, Kx1, Kx1 + Kins, Kx1 + Kins));
  double stat = sum(F2.t()*arma::solve(VF2, F2));
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa, _["VF1"] = VF1, 
                            _["VF2"] = VF2, _["Overident"] = stat, _["df"] = Kins - ntau - Kx2, 
                              _["yhat"] = yhat, _["sigma2"] = s2);
}

// This function return the structural parameters using the GMM estimates
//[[Rcpp::export]]
Rcpp::List fStructParam(const arma::vec& param,
                        const arma::mat& covp,
                        const arma::uvec& idX1,
                        const arma::uvec& idX2,
                        const int& ntau,
                        const int& Kx,
                        const int& Kx1,
                        const int& Kx2) {
  arma::uvec idx(1 + ntau + Kx);
  idx.head(ntau + 1) = arma::linspace<arma::uvec>(Kx1, Kx1 + ntau, ntau + 1);
  idx.elem(ntau + 1 + idX1) = arma::linspace<arma::uvec>(0, Kx1 - 1, Kx1);
  if (Kx2 > 0) {
    idx.elem(ntau + 1 + idX2) = arma::linspace<arma::uvec>(Kx1 + ntau + 1, Kx, Kx2);
  }
  arma::vec theta(param.elem(idx));
  arma::mat covt(covp.cols(idx));
  covt = covt.rows(idx);
  arma::mat R(arma::eye<arma::mat>(1 + ntau + Kx, 1 + ntau + Kx));
  
  // diagonal elements
  arma::vec Rd(arma::ones<arma::vec>(1 + ntau + Kx));
  Rd(0)    = -1;
  if (Kx2 > 0) {
    Rd.elem(ntau + 1 + idX2) /= param(Kx1);
  }
  R.diag() = Rd;
  
  // first column
  if (Kx2 > 0) {
    arma::vec R0(arma::zeros<arma::vec>(1 + ntau + Kx));
    R0(0)    = -1;
    R0.elem(ntau + 1 + idX2) = -param.tail(Kx2)/pow(param(Kx1), 2);
    R.col(0) = R0;
  }
  
  // theta and covariance
  theta.elem(ntau + 1 + idX2) /= theta(0);
  theta(0)                     = 1 - theta(0);
  covt       = R * covt * R.t();
  return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covt);
}

// This function estimates F stats and predict endogenous variables
// This assumes homoskedasticity
//[[Rcpp::export]]
Rcpp::List fFstathomo(const arma::mat& y,
                  const arma::mat& Xc,
                  const arma::mat& Xu) {
  int n(y.n_rows), ku(Xu.n_cols), kc(Xc.n_cols), df1(ku - kc), df2(n - ku);
  
  // constrained models
  arma::mat bc(arma::solve(Xc.t()*Xc, Xc.t()*y));
  arma::mat rc(y - Xc*bc);
  arma::rowvec ssrc(arma::sum(arma::square(rc), 0));
  
  // unconstrained models
  arma::mat bu(arma::solve(Xu.t()*Xu, Xu.t()*y));
  arma::mat ru(y - Xu*bu);
  arma::rowvec ssru(arma::sum(arma::square(ru), 0));
  
  arma::rowvec F((ssrc - ssru)*df2/(ssru*df1));
  return Rcpp::List::create(_["F"] = F, _["df1"] = df1, _["df2"] = df2, _["ru"] = ru);
}

// Same function without assuming homoskedasticity
//[[Rcpp::export]]
Rcpp::List fFstat(const arma::mat& y,
                  const arma::mat& X,
                  const arma::uvec& index,
                  const arma::mat& igroup,
                  const int& ngroup,
                  const int& HAC = 0) {
  int n(y.n_rows), K(X.n_cols), df1(index.n_elem), df2(n - K), S(y.n_cols);
  
  arma::mat XX(X.t()*X), iXX(arma::inv(XX));
  arma::mat b(iXX*X.t()*y), e(y - X*b);
  arma::vec F(S);
  for (int s(0); s < S; ++ s) {
    arma::mat V(K, K, arma::fill::zeros);
    if (HAC == 0) {
      V = sum(e.col(s)%e.col(s))*XX/(n - K);
    }
    if (HAC == 1) {
      arma::mat Xe(X.each_col()%e.col(s));
      V = Xe.t()*Xe;
    }
    if (HAC == 2) {
      for (int r(0); r < ngroup; ++ r) {
        int n1(igroup(r, 0)), n2(igroup(r, 1));
        arma::vec tp(X.rows(n1, n2).t()*e.submat(n1, s, n2, s));
        V += tp*tp.t();
      }
    }
    V   = iXX*V*iXX; V = V.rows(index); V = V.cols(index);
    arma::vec bs(b.col(s)); bs = bs.elem(index);
    F(s) = arma::sum(bs % arma::solve(V, bs))/df1;
  }
  return Rcpp::List::create(_["F"] = F, _["df1"] = df1, _["df2"] = df2, _["ru"] = e);
}

