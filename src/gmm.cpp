// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
//#define NDEBUG
// #include <RcppNumerical.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/KroneckerProduct>
// 
// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;

// using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// This function implements the GMM estimator (reduced-form)
//[[Rcpp::export]]
Rcpp::List fgmm_red(const Eigen::VectorXd& y,
                    const Eigen::MatrixXd& V,
                    const Eigen::MatrixXd& ins,
                    Eigen::MatrixXd& W,
                    const Eigen::MatrixXd& igroup,
                    const int& ngroup,
                    const int& Kx,
                    const int& Kins, 
                    const int& ntau,
                    const int& n,
                    const int& Kest, 
                    const int& HAC = 0,
                    const bool& iv = true){
  Eigen::MatrixXd ZV(ins.transpose()*V), Zy(ins.transpose()*y);
  if (iv) {
    W = (ins.transpose()*ins/n).inverse();
  }
  Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
  
  // estimate
  Eigen::VectorXd parms(VZWZV.colPivHouseholderQr().solve(VZW*Zy));
  
  // variance
  Eigen::VectorXd yhat(V*parms);
  Eigen::ArrayXd e(y - yhat);
  Eigen::MatrixXd VZe(Eigen::MatrixXd::Zero(Kins, Kins));
  double s2(R_NaN);
  if (HAC == 0) {
    s2     = e.square().sum()/(n - Kest);
    VZe    = s2*(ins.transpose()*ins);
  }
  if (HAC == 1) {
    Eigen::MatrixXd Ze(ins.array().colwise()*e);
    VZe    = Ze.transpose()*Ze;
  }
  if (HAC == 2) {
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), n2(igroup(r, 1));
      Eigen::VectorXd tp(ins(Eigen::seq(n1, n2), Eigen::all).transpose()*e(Eigen::seq(n1, n2)).matrix());
      VZe += tp*tp.transpose();
    }
  }
  // Eigen::MatrixXd tp(VZWZV.colPivHouseholderQr().solve(VZW*VZe)); // inv(V'ZWZ'V)*V'ZW*Var(Ze)
  // Eigen::MatrixXd Vpa(VZWZV.colPivHouseholderQr().solve(VZW*tp.transpose()));//inv(V'ZWZ'V)*V'ZW*Var(Ze)*WZ'V*inv(V'ZWZ'V)
  Eigen::MatrixXd iHdF((VZWZV).inverse()); 
  Eigen::MatrixXd HVFH(VZW*VZe*VZW.transpose());
  Eigen::MatrixXd Vpa(iHdF * HVFH * iHdF.transpose());
  
  // overidentification
  Eigen::VectorXd Ze(ins.transpose()*e.matrix());
  double stat = (Ze.array()*(VZe.colPivHouseholderQr().solve(Ze)).array()).sum();
  
  // criterion
  double cri, BIC, AIC, HQIC;
  if (HAC == 2) {
    cri  = Ze.dot(W*Ze)/ngroup;
    BIC  = cri - (Kins - Kx - ntau)*log(ngroup);
    AIC  = cri - 2*(Kins - Kx - ntau);
    HQIC = cri - 2.01*(Kins - Kx - ntau)*log(log(ngroup));
  } else{
    cri  = Ze.dot(W*Ze)/n;
    BIC  = cri - (Kins - Kx - ntau)*log(n);
    AIC  = cri - 2*(Kins - Kx - ntau);
    HQIC = cri - 2.01*(Kins - Kx - ntau)*log(log(n));
  }
  Rcpp::List critLis = Rcpp::List::create(_["criterion"] = cri,
                                          _["BIC"]       = BIC,
                                          _["AIC"]       = AIC,
                                          _["HQIC"]      = HQIC);
    
    return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, _["VZe"] = VZe, _["Overident"] = stat, 
                              _["df"] = Kins - Kx - ntau, _["yhat"] = yhat, _["sigma2"] = s2, _["W"] = W,
                                _["criterion"] = critLis);
}

Rcpp::List fgmm_redARMA(const arma::vec& y,
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
  return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, _["VZe"] = VZe, _["Overident"] = stat, 
                            _["df"] = Kins - Kx - ntau, _["yhat"] = yhat, _["sigma2"] = s2);
}

// This function implements the GMM estimator (Structural model)
//[[Rcpp::export]]
Rcpp::List fgmm_struc(const Eigen::VectorXd& y,
                      const Eigen::MatrixXd& X,
                      const Eigen::MatrixXd& qy,
                      const Eigen::MatrixXd& ins,
                      Eigen::MatrixXd& W1,
                      Eigen::MatrixXd& W2,
                      const Eigen::VectorXi& idX1,
                      const Eigen::VectorXi& idX2,
                      const int& Kx1,
                      const int& Kx2,
                      const Eigen::MatrixXi& igroup,
                      const Eigen::VectorXi& nIs,
                      const Eigen::VectorXi& Is,
                      const int& ngroup,
                      const int& ngroup2,
                      const int& Kins,
                      const int& Kx,
                      const int& ntau,
                      const int& n,
                      const int& Kest1,
                      const int& Kest2,
                      const int& HAC = 0,
                      const bool& iv = true){
  int n_iso(Is.size()), n_niso(n - n_iso);
  // First stage
  Eigen::VectorXd y1(y(Is));
  Eigen::MatrixXd X1(X(Is, idX1));
  Eigen::MatrixXd XX1(X1.transpose()*X1);
  if (iv) {
    W1 = (XX1/n_iso).inverse();
  }
  Eigen::MatrixXd XXW1(XX1*W1), XXWXX1(XXW1*XX1);
  Eigen::VectorXd b(XXWXX1.colPivHouseholderQr().solve(XXW1*X1.transpose()*y1));
  
  // Second stage
  Eigen::VectorXd Xb(X(Eigen::all, idX1)*b), Xb1(Xb(Is)), Xb2(Xb(nIs)), y2(y(nIs));
  Eigen::MatrixXd X2(X(nIs, idX2));
  Eigen::MatrixXd X21(X(nIs, idX1));
  Eigen::MatrixXd V2(n_niso, 1 + ntau + Kx2);
  V2 << Xb2, qy(nIs, Eigen::all), X2;
  Eigen::MatrixXd Z2(n_niso, 1 + Kins);
  Z2 << Xb2, ins(nIs, Eigen::all);
  Eigen::MatrixXd ZV2(Z2.transpose()*V2), ZZ2(Z2.transpose()*Z2); 
  if (iv) {
    W2 = (ZZ2/n_niso).inverse();
  }
  Eigen::MatrixXd VZW2(ZV2.transpose()*W2), VZWZV2(VZW2*ZV2);
  Eigen::VectorXd Zy2(Z2.transpose()*y2);
  Eigen::VectorXd lambda(VZWZV2.colPivHouseholderQr().solve(VZW2*Zy2));
  // cout<<Xb.array().sum()<<endl;
  // cout<<V2.array().sum()<<endl;
  // cout<<VZWZV2.array().sum()<<endl;
  // Variance
  Eigen::VectorXd y1hat(Xb1), y2hat(V2*lambda), yhat(n);
  yhat(Is) = y1hat; yhat(nIs) = y2hat;
  Eigen::ArrayXd e1(y1 - y1hat), e2(y2 - y2hat);
  // Eigen::ArrayXd e = y - yhat;
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kx + ntau + 1, Kx1 + Kins + 1));
  H.block(0, 0, Kx1, Kx1)  = XXW1;
  H.block(Kx1, Kx1, Kx2 + ntau + 1, Kins + 1) = VZW2;
  // cout<<H<<endl;
  
  Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(Kx1 + Kins + 1, Kx + ntau + 1));
  dF.block(0, 0, Kx1, Kx1) = XX1;
  dF(Eigen::seqN(Kx1, Kins + 1), Eigen::all) << (Z2.transpose()*X21*lambda(0)), ZV2;
  
  Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(Kx1 + Kins + 1, Kx1 + Kins + 1));
  double s21(R_NaN), s22(R_NaN);
  if (HAC == 0) {
    s21   = e1.square().sum()/(n_iso - Kest1);
    s22   = (e2/lambda(0)).square().sum()/(n_niso - Kest2);
    VF.block(0, 0, Kx1, Kx1) = s21*XX1;
    VF.block(Kx1, Kx1, Kins + 1, Kins + 1) = s22*ZZ2*pow(lambda(0), 2);
    // cout<<VF<<endl;
  }
  if (HAC == 1) {
    Eigen::MatrixXd Xe1(X1.array().colwise()*e1);
    Eigen::MatrixXd Ze2(Z2.array().colwise()*e2);
    VF.block(0, 0, Kx1, Kx1) = Xe1.transpose()*Xe1;
    VF.block(Kx1, Kx1, Kins + 1, Kins + 1) = Ze2.transpose()*Ze2;
  }
  if (HAC == 2) {
    X1   = Eigen::MatrixXd::Zero(n, Kx1);
    X1(Is, Eigen::all) = X(Is, idX1);
    Z2   = Eigen::MatrixXd::Zero(n, 1 + Kins);
    Z2(nIs, 0) = Xb2;
    Z2(nIs, Eigen::seqN(1, Kins)) = ins(nIs, Eigen::all);
    e2   = y - yhat;
    // e1.elem(nIs).zeros();
    // e2.elem(Is).zeros();
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), n2(igroup(r, 1));
      Eigen::MatrixXd tp1(n2 - n1 + 1, Kx1 + Kins + 1);
      tp1 << X1(Eigen::seq(n1, n2), Eigen::all), Z2(Eigen::seq(n1, n2), Eigen::all);
      Eigen::VectorXd tp2(tp1.transpose()*e2.matrix().segment(n1, n2));
      VF += tp2*tp2.transpose();
    }
  }
  // Eigen::MatrixXd HdF(H*dF);//H * dF
  // Eigen::MatrixXd tp(HdF.colPivHouseholderQr().solve(H*VF*H.transpose())); // inv(H * dF) * H * Var(F) * H'
  // Eigen::MatrixXd Vpa(HdF.colPivHouseholderQr().solve(tp.transpose()));//inv(H * dF) * H * Var(F) * H' * inv(H * dF)'
  Eigen::MatrixXd iHdF((H*dF).inverse()); 
  Eigen::MatrixXd HVFH(H*VF*H.transpose());
  Eigen::MatrixXd Vpa(iHdF * HVFH * iHdF.transpose());
  // cout<<Vpa<<endl;
  
  // overidentification
  Eigen::VectorXd F2(Z2.transpose()*e2.matrix());
  Eigen::MatrixXd VF1(VF.block(0, 0, Kx1, Kx1)), VF2(VF.block(Kx1, Kx1, Kins + 1, Kins + 1));
  double stat = (F2.array()*(VF2.colPivHouseholderQr().solve(F2)).array()).sum();
  
  // criterion
  double cri, BIC, AIC, HQIC;
  if (HAC == 2) {
    cri  = F2.dot(W2*F2)/ngroup2;
    BIC  = cri - (Kins - Kx2 - ntau)*log(ngroup2);
    AIC  = cri - 2*(Kins - Kx2 - ntau);
    HQIC = cri - 2.01*(Kins - Kx2 - ntau)*log(log(ngroup2));
  } else{
    cri  = F2.dot(W2*F2)/n_niso;
    BIC  = cri - (Kins - Kx2 - ntau)*log(n_niso);
    AIC  = cri - 2*(Kins - Kx2 - ntau);
    HQIC = cri - 2.01*(Kins - Kx2 - ntau)*log(log(n_niso));
  }
  Rcpp::List critLis = Rcpp::List::create(_["criterion"] = cri,
                                          _["BIC"]       = BIC,
                                          _["AIC"]       = AIC,
                                          _["HQIC"]      = HQIC);
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa, _["VF1"] = VF1, 
                            _["VF2"] = VF2, _["Overident"] = stat, _["df"] = Kins - ntau - Kx2, 
                              _["yhat"] = yhat, _["sigma21"] = s21, _["sigma22"] = s22, _["W1"] = W1, 
                                _["W2"] = W2, _["criterion"] = critLis);
}


// Not exported, No longer used.
Rcpp::List fgmm_strucARMA(const arma::vec& y,
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
  dF.submat(0, 0, Kx1 - 1, Kx1 - 1) = XX1;
  dF.submat(Kx1, 0, Kx1 + Kins, Kx1 - 1) = (Z2.t()*X21*lambda(0));
  dF.submat(Kx1, Kx1, Kx1 + Kins, Kx + ntau) = ZV2;
  
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
// Transformations were applied to estimate a function of parameter
// Here we return the true parameters
//[[Rcpp::export]]
Rcpp::List fStructParam(const arma::vec& param,
                        const arma::mat& covp,
                        const arma::uvec& idX1,
                        const arma::uvec& idX2,
                        const int& ntau,
                        const int& Kx,
                        const int& Kx1,
                        const int& Kx2, 
                        const bool& COV = true) {
  arma::uvec idx(1 + ntau + Kx);
  idx.head(ntau + 1)          = arma::linspace<arma::uvec>(Kx1, Kx1 + ntau, ntau + 1);
  idx.elem(ntau + 1 + idX1)   = arma::linspace<arma::uvec>(0, Kx1 - 1, Kx1);
  if (Kx2 > 0) {
    idx.elem(ntau + 1 + idX2) = arma::linspace<arma::uvec>(Kx1 + ntau + 1, Kx + ntau, Kx2);
  }
  
  // Theta
  arma::vec theta(param.elem(idx));
  theta.elem(ntau + 1 + idX2) /= theta(0);
  theta(0)                     = 1 - theta(0);
  
  // Covariance
  if (COV) {
    arma::mat covt(covp.cols(idx));
    covt      = covt.rows(idx);
    arma::mat R(arma::eye<arma::mat>(1 + ntau + Kx, 1 + ntau + Kx));
    
    // diagonal elements
    arma::vec Rd(arma::ones<arma::vec>(1 + ntau + Kx));
    Rd(0)     = -1;
    if (Kx2 > 0) {
      Rd.elem(ntau + 1 + idX2) /= param(Kx1);
    }
    R.diag()  = Rd;
    
    // first column
    if (Kx2 > 0) {
      arma::vec R0(arma::zeros<arma::vec>(1 + ntau + Kx));
      R0(0)    = -1;
      R0.elem(ntau + 1 + idX2) = -param.tail(Kx2)/pow(param(Kx1), 2);
      R.col(0) = R0;
    }
    
    // theta and covariance
    covt       = R * covt * R.t();
    return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covt);
  } else {
    return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covp);
  }
  
}

// This function estimates F stats and predict endogenous variables
// This assumes homoskedasticity
//[[Rcpp::export]]
Rcpp::List fFstathomo(const Eigen::MatrixXd& y,
                      const Eigen::MatrixXd& Xc,
                      const Eigen::MatrixXd& Xu) {
  int n(y.rows()), ku(Xu.cols()), kc(Xc.cols()), df1(ku - kc), df2(n - ku);
  
  // constrained models
  Eigen::MatrixXd bc((Xc.transpose()*Xc).colPivHouseholderQr().solve(Xc.transpose()*y));
  Eigen::MatrixXd rc(y - Xc*bc);
  Eigen::ArrayXXd ssrc(rc.array().square().colwise().sum());
  
  // unconstrained models
  Eigen::MatrixXd bu((Xu.transpose()*Xu).colPivHouseholderQr().solve(Xu.transpose()*y));
  Eigen::MatrixXd ru(y - Xu*bu);
  Eigen::ArrayXXd ssru(ru.array().square().colwise().sum());
  
  Eigen::ArrayXXd F((ssrc - ssru)*df2/(ssru*df1));
  return Rcpp::List::create(_["F"] = F, _["df1"] = df1, _["df2"] = df2, _["ru"] = ru);
}

//[[Rcpp::export]]
Rcpp::List fFstathomoARMA(const arma::mat& y,
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
Rcpp::List fFstat(const Eigen::MatrixXd& y,
                  const Eigen::MatrixXd& X,
                  const Eigen::VectorXi& index,
                  const Eigen::MatrixXd& igroup,
                  const int& ngroup,
                  const int& HAC = 0) {
  int n(y.rows()), K(X.cols()), df1(index.size()), df2(n - K), S(y.cols());
  
  Eigen::MatrixXd XX(X.transpose()*X), iXX(XX.inverse());
  Eigen::MatrixXd b(iXX*X.transpose()*y);
  Eigen::ArrayXXd e(y - X*b);
  Eigen::VectorXd F(S);
  for (int s(0); s < S; ++ s) {
    Eigen::MatrixXd V(Eigen::MatrixXd::Zero(K, K));
    if (HAC == 0) {
      V = (e.col(s).square().sum())*XX/(n - K);
    }
    if (HAC == 1) {
      Eigen::MatrixXd Xe((X.array().colwise()*e.col(s)));
      V = Xe.transpose()*Xe;
    }
    if (HAC == 2) {
      for (int r(0); r < ngroup; ++ r) {
        int n1(igroup(r, 0)), n2(igroup(r, 1));
        Eigen::VectorXd tp(X(Eigen::seq(n1, n2), Eigen::all).transpose()*e(Eigen::seq(n1, n2), s).matrix());
        V += tp*tp.transpose();
      }
    }
    Eigen::MatrixXd tp(iXX*V*iXX);
    V    = tp(index, index);
    Eigen::VectorXd bs(b(index, s));
    F(s) = (bs.array() * (V.colPivHouseholderQr().solve(bs)).array()).sum()/df1;
  }
  return Rcpp::List::create(_["F"] = F, _["df1"] = df1, _["df2"] = df2, _["ru"] = e);
}

//[[Rcpp::export]]
Rcpp::List fFstatARMA(const arma::mat& y,
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

// Computes sqrt of matrices
Eigen::MatrixXd matrixSqrt(const Eigen::MatrixXd& A) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  Eigen::VectorXd sqrt_evals = es.eigenvalues().array().sqrt();
  return es.eigenvectors() * sqrt_evals.asDiagonal() * es.eigenvectors().transpose();
}

// This function computes KP stat
//[[Rcpp::export]]
Rcpp::List fKPstat(const Eigen::MatrixXd& qy_,
                   const Eigen::MatrixXd& X,
                   const Eigen::MatrixXd& Z_,
                   const arma::uvec& index,
                   const Eigen::MatrixXd& igroup,
                   const int& HAC = 0) {
  // cout<<qy_.cols()<<" "<<qy_.rows()<<endl;
  // cout<<X.cols()<<" "<<X.rows()<<endl;
  // cout<<Z_.cols()<<" "<<Z_.rows()<<endl;
  Eigen::MatrixXd iXX((X.transpose() * X).inverse());
  Eigen::MatrixXd qy(qy_ - X * iXX * X.transpose() * qy_);
  Eigen::MatrixXd Z(Z_(Eigen::all, index) - X * iXX * X.transpose() * Z_(Eigen::all, index));
  // Eigen::MatrixXd qy(qy_);
  // Eigen::MatrixXd Z(Z_);
  int n(qy.rows()), ntau(qy.cols()), l(Z.cols()), ngroup(igroup.rows());
  Eigen::MatrixXd ZZ(Z.transpose() * Z);
  Eigen::MatrixXd iZZ(ZZ.inverse());
  Eigen::MatrixXd Zqy(Z.transpose() * qy);
  
  // estimator
  Eigen::MatrixXd Pi(Zqy.transpose() * iZZ);
  Eigen::VectorXd pi(Pi.reshaped(l * ntau, 1)); // Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(l, l), Zqy.transpose()) * ZZ.inverse().reshaped(l*l, 1)
  
  // vec(Ze)
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(l * ntau, l * ntau));
  for (int s1(0); s1 < l; ++ s1) {
    for (int s2(0); s2 < ntau; ++ s2) {
      R(s1 * ntau + s2, s2 * l + s1) = 1;
    }
  }
  
  Eigen::MatrixXd eps(qy - Z * Pi.transpose());
  Eigen::MatrixXd vecZe(n, l*ntau);
  for (int s(0); s < ntau; ++ s) {
    vecZe.block(0, s*l, n, l) = (Z.array().colwise()*eps.col(s).array()).matrix();
  }
  
  // Variance of vec(Ze), covqy and covz
  Eigen::MatrixXd VvecZe(Eigen::MatrixXd::Zero(l*ntau, l*ntau)),
  Eee(Eigen::MatrixXd::Zero(ntau, ntau)),
  Ezz(Eigen::MatrixXd::Zero(l, l));
  if (HAC <= 1) {
    VvecZe = vecZe.transpose() * vecZe;
    Eee    = eps.transpose() * eps;
    Ezz    = Z.transpose() * Z;
  } else if (HAC == 2) {
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), n2(igroup(r, 1));
      Eigen::VectorXd tp(vecZe(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix());
      VvecZe += tp * tp.transpose();
      tp      = eps(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix();
      Eee    += tp * tp.transpose();
      tp      = Z(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix();
      Ezz    += tp * tp.transpose();
    }
  }
  
  // Variance of pi
  Eigen::MatrixXd H(R * Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(ntau, ntau), iZZ));
  Eigen::MatrixXd varpi(H * VvecZe * H.transpose()); // O(1/n)
  
  // normalisation
  Eigen::LLT<Eigen::MatrixXd> tpF(ZZ * Ezz.colPivHouseholderQr().solve(ZZ)), tpG(Eee.inverse());
  Eigen::MatrixXd F(tpF.matrixL()); // O(sqrt(n))
  Eigen::MatrixXd G(tpG.matrixL().transpose()); // O(1/sqrt(n))
  F *= sqrt(ngroup); // O(n)
  
  // Theta and its variance
  Eigen::MatrixXd Theta(G * Pi * F.transpose());
  Eigen::VectorXd theta(Theta.reshaped(l*ntau, 1));
  Eigen::MatrixXd FG(Eigen::KroneckerProduct(F, G));
  Eigen::MatrixXd vartheta(FG * varpi * FG.transpose());
  // cout << vartheta << endl;
  
  // SDV decomposition of Theta
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Theta, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd U = svd.matrixU(); //ntau * ntau
  Eigen::VectorXd d = svd.singularValues();
  Eigen::MatrixXd ddiag = d.asDiagonal();
  Eigen::MatrixXd D(ntau, l);
  D << ddiag, Eigen::MatrixXd::Zero(ntau, l - ntau); //l*ntau
  Eigen::MatrixXd V = svd.matrixV(); //ntau * ntau
  
  //U12, U22, V12, V22
  int q(ntau - 1);
  Eigen::MatrixXd U12(U.block(0, q, q, ntau - q));
  Eigen::MatrixXd U22(U.block(q, q, ntau - q, ntau - q));
  Eigen::MatrixXd V12(V.block(0, q, q, l - q));
  Eigen::MatrixXd V22(V.block(q, q, l - q, l - q));
  
  // Aqper and Bqper
  Eigen::MatrixXd U12U22(ntau, ntau - q), V12V22(l, l - q);
  U12U22 << U12, U22;
  V12V22 << V12, V22;
  
  Eigen::MatrixXd Aper(U12U22 * U22.colPivHouseholderQr().solve(matrixSqrt(U22 * U22.transpose())));
  Eigen::MatrixXd Bper((V12V22 * V22.colPivHouseholderQr().solve(matrixSqrt(V22 * V22.transpose()))).transpose());
  
  // lambda and its varianve
  Eigen::MatrixXd BAper(Eigen::KroneckerProduct(Bper, Aper.transpose()));
  Eigen::VectorXd lambda (BAper * theta);
  Eigen::MatrixXd varlambda (BAper * vartheta * BAper.transpose());
  
  // statistic
  double stat((lambda.transpose() * varlambda.colPivHouseholderQr().solve(lambda))(0, 0));
  return Rcpp::List::create(_["stat"] = stat, _["df"] = (ntau - q)*(l - q));
}

// This function adds spillover effects to fStructParam which only have total effects and conformity
//[[Rcpp::export]]
Rcpp::List fStructParamFull(const arma::vec& param,
                            const arma::mat& covp,
                            const int& ntau,
                            const int& Kx1,
                            const int& Kx2,
                            const int& quantile,
                            const int& ces = false) {
  int Kx(Kx1 + Kx2);
  arma::vec theta;
  arma::mat covt;
  
  if (ces) {
    arma::mat R(arma::zeros<arma::mat>(4 + Kx, 3 + Kx));
    R.submat(4, 3, 3 + Kx, 2 + Kx).eye();
    // rho
    R(0, 0) =  1;
    // Spillover : index 2 - index 1
    R(1, 1) = -1;
    R(1, 2) =  1;
    // Conformity index 1
    R(2, 1) =  1; 
    // total index 2
    R(3, 2) =  1; 
    
    // theta
    theta = R*param;
    covt  = R*covp*R.t();
  } else{
    arma::mat R(arma::zeros<arma::mat>(2 + quantile + ntau + Kx, 1 + ntau + Kx));
    R.submat(2 + quantile, 1, 1 + quantile + ntau + Kx, ntau + Kx).eye();
    // Spillover sum(index 1:ntau) - index 0
    R(0, 0)                 = -1;
    R.submat(0, 1, 0, ntau).ones();
    // Conformity index 0
    R(1, 0)                 =  1; 
    // total if quantile sum(index 1:ntau)
    if (quantile == 1) {
      R.submat(2, 1, 2, ntau).ones();
    }
    
    // theta
    theta = R*param;
    covt  = R*covp*R.t();
  }
  
  return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covt);
}

// This function add total effects 
//[[Rcpp::export]]
Rcpp::List fParamFull(const arma::vec& param,
                      const arma::mat& covp,
                      const int& ntau,
                      const int& Kx1,
                      const int& Kx2) {
  int Kx(Kx1 + Kx2);
  arma::mat R(arma::zeros<arma::mat>(1 + ntau + Kx, ntau + Kx));
  
  // total if quantile sum(index 0:(ntau - 1))
  R.submat(0, 0, 0, ntau - 1).ones();
  R.submat(1, 0, ntau + Kx, ntau + Kx - 1).eye();
  
  // theta
  arma::vec theta(R*param);
  arma::mat covt(R*covp*R.t());
  return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covt);
}

// DC test adapted to many variables, But does not make sense,
// Rcpp::List fDCstat(const Eigen::MatrixXd& qy,
//                    const Eigen::MatrixXd& Z,
//                    const Eigen::MatrixXd& igroup,
//                    const int& ngroup,
//                    const int& HAC = 0) {
//   int n(qy.rows()), ntau(qy.cols()), l(Z.cols());
//   Eigen::MatrixXd ZZ(Z.transpose() * Z);
//   Eigen::MatrixXd Zqy(Z.transpose() * qy);
//   Eigen::MatrixXd eps(qy - Z * ZZ.colPivHouseholderQr().solve(Zqy));
//   
//   // vec(Ze)
//   Eigen::MatrixXd vecZe(n, l*ntau);
//   for (int s(0); s < ntau; ++ s) {
//     vecZe.block(0, s*l, n, l) = (Z.array().colwise()*eps.col(s).array()).matrix();
//   }
//   
//   // Variance of vec(Ze)
//   Eigen::MatrixXd VvecZe(Eigen::MatrixXd::Zero(l*ntau, l*ntau));
//   if (HAC <= 1) {
//     VvecZe =vecZe.transpose()*vecZe;
//   } else if (HAC == 2) {
//     for (int r(0); r < ngroup; ++ r) {
//       int n1(igroup(r, 0)), n2(igroup(r, 1));
//       Eigen::VectorXd tp(vecZe(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix());
//       VvecZe += tp*tp.transpose();
//     }
//   }
//   
//   // stat
//   Eigen::VectorXd vecZqy(Zqy.reshaped(l*ntau, 1));
//   double stat = (vecZqy.transpose() * VvecZe.selfadjointView<Eigen::Lower>().ldlt().solve(vecZqy))(0, 0) / ntau;
//   
//   // critical values
//   Eigen::VectorXd CV30(30), CV20(30), CV10(30), CV05(30);
//   CV30 << 12.05, 9.57, 8.53, 7.92, 7.51, 7.21, 6.98, 6.80, 6.65, 6.52,
//           6.41, 6.32, 6.24, 6.16, 6.10, 6.04, 5.99, 5.94, 5.89, 5.85,
//           5.81, 5.78, 5.74, 5.71, 5.68, 5.66, 5.63, 5.61, 5.58, 5.56;
//   
//   CV20 << 15.06, 12.17, 10.95, 10.23, 9.75, 9.40, 9.14, 8.92, 8.74, 8.59,
//           8.47, 8.36, 8.26, 8.17, 8.10, 8.03, 7.96, 7.91, 7.85, 7.80,
//           7.76, 7.72, 7.68, 7.64, 7.61, 7.57, 7.54, 7.51, 7.49, 7.46;
//   
//   CV10 << 23.11, 19.29, 17.67, 16.72, 16.08, 15.62, 15.26, 14.97, 14.73, 14.53,
//           14.36, 14.21, 14.08, 13.96, 13.86, 13.77, 13.68, 13.60, 13.53, 13.46,
//           13.40, 13.35, 13.29, 13.24, 13.20, 13.15, 13.11, 13.07, 13.04, 13.00;
//   
//   CV05 << 37.42, 32.32, 30.13, 28.85, 27.98, 27.35, 26.86, 26.47, 26.15, 25.87,
//           25.64, 25.44, 25.26, 25.10, 24.96, 24.83, 24.71, 24.60, 24.50, 24.41,
//           24.33, 24.25, 24.18, 24.11, 24.05, 23.98, 23.93, 23.87, 23.82, 23.77;
//   
//   if (ntau <= 30) {
//     Eigen::VectorXd CV(4);
//     CV << CV30[ntau - 1], CV20[ntau - 1], CV10[ntau - 1], CV05[ntau - 1];
//     return Rcpp::List::create(
//       _["stat"] = stat,
//       _["cv"] = CV,
//       _["cv_names"] = Rcpp::CharacterVector::create("Bias tolerance 30%", "  20%", "  10%", "  5%")
//     );
//   } else {
//     Rcpp::CharacterVector CV = Rcpp::CharacterVector::create(
//       "< 5.56", "< 7.46", "< 13.00", "< 23.77"
//     );
//     return Rcpp::List::create(
//       _["stat"] = stat,
//       _["cv"] = CV,
//       _["cv_names"] = Rcpp::CharacterVector::create("Bias tolerance 30%", "  20%", "  10%", "  5%")
//     );
//   }
// }