// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppEigen)]]
#include <RcppArmadillo.h>
//#define NDEBUG
// #include <RcppNumerical.h>
#include <random>
#include <RcppEigen.h>
#include <unsupported/Eigen/KroneckerProduct>
#include <progress.hpp>
#include <progress_bar.hpp>
// 
// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;

#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

// using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// generalized inverse
Eigen::MatrixXd ginv_gmm(const Eigen::MatrixXd& A) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  
  const Eigen::VectorXd& singularValues = svd.singularValues();
  Eigen::VectorXd singularValuesInv = singularValues;
  
  // Compute robust tolerance
  double tol = std::numeric_limits<double>::epsilon()
    * std::max(A.rows(), A.cols())
    * singularValues(0);  // σ_max
    
    for (int i = 0; i < singularValues.size(); ++i) {
      if (singularValues(i) > tol)
        singularValuesInv(i) = 1.0 / singularValues(i);
      else
        singularValuesInv(i) = 0.0;
    }
    
    Eigen::MatrixXd Dplus = singularValuesInv.asDiagonal();
    return svd.matrixV() * Dplus * svd.matrixU().transpose();
}


// Computes sqrt of matrices
Eigen::MatrixXd matrixSqrt(const Eigen::MatrixXd& A) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  Eigen::VectorXd sqrt_evals = es.eigenvalues().array().sqrt();
  return es.eigenvectors() * sqrt_evals.asDiagonal() * es.eigenvectors().transpose();
}

// This function implements the GMM estimator (reduced-form)
//[[Rcpp::export]]
Rcpp::List fgmm_red(const Eigen::VectorXd& y,
                    const Eigen::MatrixXd& V,
                    const Eigen::MatrixXd& ins,
                    Eigen::MatrixXd& W,
                    const Eigen::ArrayXi& igroup,
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
      int n1(igroup(r)), n2(igroup(r + 1) - 1);
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
  double stat = Ze.dot(ginv_gmm(VZe) * Ze);
  
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
                      const Eigen::ArrayXi& igroup,
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
  Eigen::MatrixXd Z2(ins(nIs, Eigen::all));
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
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kx + ntau + 1, Kx1 + Kins));
  H.block(0, 0, Kx1, Kx1)  = XXW1;
  H.block(Kx1, Kx1, Kx2 + ntau + 1, Kins) = VZW2;
  // cout<<H<<endl;
  
  Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(Kx1 + Kins, Kx + ntau + 1));
  dF.block(0, 0, Kx1, Kx1) = XX1;
  dF(Eigen::seqN(Kx1, Kins), Eigen::all) << (Z2.transpose()*X21*lambda(0)), ZV2;
  
  Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(Kx1 + Kins, Kx1 + Kins));
  double s21(R_NaN), s22(R_NaN);
  if (HAC == 0) {
    s21   = e1.square().sum()/(n_iso - Kest1);
    s22   = (e2/lambda(0)).square().sum()/(n_niso - Kest2);
    VF.block(0, 0, Kx1, Kx1) = s21*XX1;
    VF.block(Kx1, Kx1, Kins, Kins) = s22*ZZ2*pow(lambda(0), 2);
    // cout<<VF<<endl;
  }
  if (HAC == 1) {
    Eigen::MatrixXd Xe1(X1.array().colwise()*e1);
    Eigen::MatrixXd Ze2(Z2.array().colwise()*e2);
    VF.block(0, 0, Kx1, Kx1) = Xe1.transpose()*Xe1;
    VF.block(Kx1, Kx1, Kins, Kins) = Ze2.transpose()*Ze2;
  }
  if (HAC == 2) {
    X1   = Eigen::MatrixXd::Zero(n, Kx1);
    X1(Is, Eigen::all) = X(Is, idX1);
    Z2   = Eigen::MatrixXd::Zero(n, Kins);
    Z2(nIs, Eigen::all) = ins(nIs, Eigen::all);
    e2   = y - yhat;
    // e1.elem(nIs).zeros();
    // e2.elem(Is).zeros();
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r)), nr(igroup(r + 1) - n1);
      Eigen::MatrixXd tp1(nr, Kx1 + Kins);
      tp1 << X1(Eigen::seqN(n1, nr), Eigen::all), Z2(Eigen::seqN(n1, nr), Eigen::all);
      Eigen::VectorXd tp2(tp1.transpose()*e2.matrix().segment(n1, nr));
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
  Eigen::MatrixXd VF1(VF.block(0, 0, Kx1, Kx1)), VF2(VF.block(Kx1, Kx1, Kins, Kins));
  double stat = F2.dot(ginv_gmm(VF2) * F2);
  
  // criterion
  double cri, BIC, AIC, HQIC;
  if (HAC == 2) {
    cri  = F2.dot(W2*F2)/ngroup2;
    BIC  = cri - (Kins - Kx2 - 1 - ntau)*log(ngroup2);
    AIC  = cri - 2*(Kins - Kx2 - 1 - ntau);
    HQIC = cri - 2.01*(Kins - Kx2 - 1 - ntau)*log(log(ngroup2));
  } else{
    cri  = F2.dot(W2*F2)/n_niso;
    BIC  = cri - (Kins - Kx2 - 1 - ntau)*log(n_niso);
    AIC  = cri - 2*(Kins - Kx2 - 1 - ntau);
    HQIC = cri - 2.01*(Kins - Kx2 - 1 - ntau)*log(log(n_niso));
  }
  Rcpp::List critLis = Rcpp::List::create(_["criterion"] = cri,
                                          _["BIC"]       = BIC,
                                          _["AIC"]       = AIC,
                                          _["HQIC"]      = HQIC);
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa, _["VF1"] = VF1, 
                            _["VF2"] = VF2, _["Overident"] = stat, _["df"] = Kins - Kx2 - 1 - ntau, 
                              _["yhat"] = yhat, _["sigma21"] = s21, _["sigma22"] = s22, _["W1"] = W1, 
                                _["W2"] = W2, _["criterion"] = critLis);
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
//[[Rcpp::export]]
Rcpp::List fFstat(const Eigen::MatrixXd& y,
                  const Eigen::MatrixXd& X,
                  const Eigen::VectorXi& index,
                  const Eigen::ArrayXd& igroup,
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
        int n1(igroup(r)), n2(igroup(r + 1) - 1);
        Eigen::VectorXd tp(X(Eigen::seq(n1, n2), Eigen::all).transpose()*e(Eigen::seq(n1, n2), s).matrix());
        V += tp*tp.transpose();
      }
    }
    Eigen::MatrixXd tp(iXX*V*iXX);
    V    = tp(index, index);
    Eigen::VectorXd bs(b(index, s));
    F(s) = bs.dot(ginv_gmm(V) * bs) / df1;
  }
  return Rcpp::List::create(_["F"] = F, _["df1"] = df1, _["df2"] = df2, _["ru"] = e);
}

// This function computes KP stat
//[[Rcpp::export]]
Rcpp::List fKPstat(const Eigen::MatrixXd& qy,
                   const Eigen::MatrixXd& Z,
                   const arma::uvec& index,
                   const Eigen::ArrayXd& igroup,
                   const int& HAC = 0) {
  int n(qy.rows()), ntau(qy.cols()), l(index.size()), Kins(Z.cols()), 
  Kx(Kins - l), ngroup(igroup.size() - 1);
  Eigen::MatrixXd ZZ(Z.transpose() * Z);
  Eigen::MatrixXd iZZ(ZZ.inverse());
  Eigen::MatrixXd Zqy(Z.transpose() * qy);
  
  // estimator and residuals
  Eigen::MatrixXd Pi(Zqy.transpose() * iZZ);
  Eigen::MatrixXd eps(qy - Z * Pi.transpose());
  Eigen::MatrixXd selPi = Pi(Eigen::all, index);
  Eigen::MatrixXd selZZ = ZZ(index, Eigen::all);
  Eigen::VectorXd pi(selPi.reshaped(l * ntau, 1)); // Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(l, l), Zqy.transpose()) * ZZ.inverse().reshaped(l*l, 1)
  
  // vec(Ze)
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(l * ntau, Kins * ntau));
  for (int s1(0); s1 < l; ++ s1) {
    for (int s2(0); s2 < ntau; ++ s2) {
      R(s1 * ntau + s2, s2 * Kins + s1 + Kx) = 1;
    }
  }
  
  
  Eigen::MatrixXd vecZe(n, Kins*ntau);
  for (int s(0); s < ntau; ++ s) {
    vecZe.block(0, s*Kins, n, Kins) = Z.array().colwise()*eps.col(s).array();
  }
  
  // Variance of vec(Ze), covqy and covz
  Eigen::MatrixXd VvecZe(Eigen::MatrixXd::Zero(Kins*ntau, Kins*ntau));
  if (HAC <= 1) {
    VvecZe = vecZe.transpose() * vecZe;
  } else if (HAC == 2) {
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r)), n2(igroup(r + 1) - 1);
      Eigen::VectorXd tp(vecZe(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix());
      VvecZe += tp * tp.transpose();
    }
  }
  
  // Variance of pi
  Eigen::MatrixXd H(R * Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(ntau, ntau), iZZ));
  Eigen::MatrixXd varpi(H * VvecZe * H.transpose()); // O(1/n)
  
  // covz and cove
  Eigen::MatrixXd dZ    = Z.array().rowwise() - Z.array().colwise().mean();
  Eigen::MatrixXd dqy   = qy.array().rowwise() - qy.array().colwise().mean();
  Eigen::MatrixXd covz  = dZ.transpose() * dZ / ngroup;
  Eigen::MatrixXd covqy = dqy.transpose() * dqy / ngroup;
  
  // normalisation
  Eigen::LLT<Eigen::MatrixXd> tpF(selZZ * covz.colPivHouseholderQr().solve(selZZ.transpose())), tpG(covqy.inverse());
  Eigen::MatrixXd F(tpF.matrixL()); // O(sqrt(n))
  Eigen::MatrixXd G(tpG.matrixL().transpose()); // O(1/sqrt(n))
  
  // Theta and its variance
  Eigen::MatrixXd Theta(G * selPi * F.transpose());
  Eigen::VectorXd theta(Theta.reshaped(l*ntau, 1));
  Eigen::MatrixXd FG(Eigen::kroneckerProduct(F, G));
  Eigen::MatrixXd vartheta(FG * varpi * FG.transpose());
  // cout << F << endl;
  // cout << G << endl;
  // cout << Pi << endl;
  // cout << vartheta << endl;
  // Until this, replicate using bootstrap
  
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
  Eigen::MatrixXd BAper(Eigen::kroneckerProduct(Bper, Aper.transpose()));
  Eigen::VectorXd lambda (BAper * theta);
  Eigen::MatrixXd varlambda (BAper * vartheta * BAper.transpose());
  // cout << theta << endl;
  // cout << varlambda << endl;
  // statistic
  double stat = lambda.dot(ginv_gmm(varlambda) * lambda);
  
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


/////////////////////////////////////// Bootstrap
// This function implements the GMM estimator (reduced-form), for one iteration of bootstrap
void fgmm_red_bootcoef0(Eigen::VectorXd& theta,
                        Eigen::VectorXd& Ze,
                        double& Jstat,
                        const std::vector<Eigen::MatrixXd>& LZV,
                        const std::vector<Eigen::VectorXd>& LZy,
                        const std::vector<Eigen::MatrixXd>& LZZ,
                        Eigen::MatrixXd& W,
                        const int& ngroup,
                        const int& Kins,
                        const int& Kv,
                        const bool& iv,
                        const bool& overident){
  Eigen::MatrixXd ZV = Eigen::MatrixXd::Zero(Kins, Kv);
  Eigen::MatrixXd ZZ = Eigen::MatrixXd::Zero(Kins, Kins);
  Eigen::VectorXd Zy = Eigen::VectorXd::Zero(Kins);
  
  // Find iso and nisolated
  for (int s = 0; s < ngroup; ++s) {
    ZV += LZV[s];
    Zy += LZy[s];
    if (iv) ZZ += LZZ[s];
  }
  if (iv) {
    W = (ZZ / ngroup).ldlt().solve(Eigen::MatrixXd::Identity(Kins, Kins));
  }
  Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
  
  // estimate
  theta = VZWZV.colPivHouseholderQr().solve(VZW*Zy);
  Ze    =  (Zy - ZV * theta);
  if (overident) {
    Eigen::MatrixXd VZe = Eigen::MatrixXd::Zero(Kins, Kins);
    for (int s = 0; s < ngroup; ++s) {
      Eigen::VectorXd Zes(LZy[s] -  LZV[s] * theta);
      VZe += Zes * Zes.transpose();
    }
    Jstat = Ze.dot(ginv_gmm(VZe) * Ze);
  }
}


void fgmm_red_bootcoef(Eigen::VectorXd& theta,
                       double& Jstat,
                       const Eigen::VectorXd& Ze0,
                       const std::vector<Eigen::MatrixXd>& LZV,
                       const std::vector<Eigen::VectorXd>& LZy,
                       const std::vector<Eigen::MatrixXd>& LZZ,
                       Eigen::MatrixXd& W,
                       const Eigen::ArrayXi& sgroup,
                       const int& ngroup,
                       const int& Kins,
                       const int& Kv,
                       const bool& iv,
                       const bool& overident){
  Eigen::MatrixXd ZV = Eigen::MatrixXd::Zero(Kins, Kv);
  Eigen::MatrixXd ZZ = Eigen::MatrixXd::Zero(Kins, Kins);
  Eigen::VectorXd Zy = Eigen::VectorXd::Zero(Kins);
  
  // Find iso and nisolated
  for (int s = 0; s < ngroup; ++s) {
    ZV += LZV[sgroup(s)];
    Zy += LZy[sgroup(s)];
    if (iv) ZZ += LZZ[sgroup(s)];
  }
  if (iv) {
    W = (ZZ / ngroup).ldlt().solve(Eigen::MatrixXd::Identity(Kins, Kins));
  }
  Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
  
  // estimate
  Eigen::ColPivHouseholderQR <Eigen::MatrixXd> qr(VZWZV);
  theta = qr.solve(VZW*Zy);
  
  // J distribution
  if (overident) {
    Eigen::VectorXd Zyc     = Zy - Ze0;
    Eigen::VectorXd theta_c = qr.solve(VZW * Zyc);
    Eigen::VectorXd Ze      =  Zyc - ZV * theta_c;
    Eigen::MatrixXd VZe = Eigen::MatrixXd::Zero(Kins, Kins);
    for (int s = 0; s < ngroup; ++s) {
      Eigen::VectorXd Zes(LZy[sgroup(s)] -  LZV[sgroup(s)] * theta_c - Ze0 / ngroup);
      VZe += Zes * Zes.transpose();
    }
    Jstat = Ze.dot(ginv_gmm(VZe) * Ze);
  }
}


//[[Rcpp::export]]
Rcpp::List fgmm_red_boot(const Eigen::VectorXd& y,
                         const Eigen::MatrixXd& V,
                         const Eigen::MatrixXd& ins,
                         const Eigen::MatrixXd& W,
                         const Eigen::ArrayXi& igroup,
                         const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                         const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                         const int& ngroup,
                         const int& Kx,
                         const int& Kins, 
                         const int& ntau,
                         const int& Kest, 
                         const bool& iv,
                         const int& boot,
                         const int& nthreads,
                         const unsigned long long seed,
                         const bool& print,
                         const bool& overident){
  int Kv(ntau + Kx);
  
  // Compute blocks
  std::vector<Eigen::MatrixXd> LZV(ngroup);
  std::vector<Eigen::VectorXd> LZy(ngroup);
  std::vector<Eigen::MatrixXd> LZZ(ngroup);
  
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#pragma omp parallel for
  for (int s = 0; s < ngroup; ++s) {
    auto idx = Eigen::seq(igroup(s), igroup(s + 1) - 1);
    LZV[s] = ins(idx, Eigen::all).transpose() * V(idx, Eigen::all);
    LZy[s] = ins(idx, Eigen::all).transpose() * y(idx);
    if (iv) LZZ[s] = ins(idx, Eigen::all).transpose() * ins(idx, Eigen::all);
  }
#else
  for (int s = 0; s < ngroup; ++s) {
    auto idx = Eigen::seq(igroup(s), igroup(s + 1) - 1);
    LZV[s] = ins(idx, Eigen::all).transpose() * V(idx, Eigen::all);
    LZy[s] = ins(idx, Eigen::all).transpose() * y(idx);
    if (iv) LZZ[s] = ins(idx, Eigen::all).transpose() * ins(idx, Eigen::all);
  }
#endif
  
  Progress Prog(boot + 1, print);
  
  // First iteration
  Eigen::MatrixXd W0(W);
  Eigen::VectorXd theta, Ze;
  double stat = R_NaN;
  fgmm_red_bootcoef0(theta, Ze, stat, LZV, LZy, LZZ, W0, ngroup, Kins, Kv, 
                     iv, overident);
  Prog.increment();
  
  // Where to save ltheta
  Eigen::MatrixXd ltheta(ntau + Kx, boot);
  Eigen::ArrayXd lstat(boot);
  
  //setup parallel settings
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  std::mt19937 rng(seed + tid * 7919); 
  
  Eigen::VectorXd theta_loc;
  Eigen::MatrixXd W_loc;
  double stat_loc;
  
#pragma omp for
  for (int k = 0; k < boot; ++ k) {
    W_loc = W;
    // Select subnets
    Eigen::ArrayXi sgroup(ngroup);
    std::uniform_int_distribution<int> unidist(0, ngroup - 1);
    for (int s = 0; s < ngroup; ++ s) {
      sgroup(s)   =  unidist(rng);
    }
    
    fgmm_red_bootcoef(theta_loc, stat_loc, Ze, LZV, LZy, LZZ, W_loc, sgroup,
                      ngroup, Kins, Kv, iv, overident);
    ltheta.col(k) = theta_loc;
    if (overident) {
      lstat(k)    = stat_loc;
    }
#pragma omp critical
    Prog.increment();
  }
}
#else
std::mt19937 rng(seed); 

Eigen::VectorXd theta_loc;
Eigen::MatrixXd W_loc;
double stat_loc;
for (int k = 0; k < boot; ++ k) {
  W_loc = W;
  // Select subnets
  Eigen::ArrayXi sgroup(ngroup);
  std::uniform_int_distribution<int> unidist(0, ngroup - 1);
  for (int s = 0; s < ngroup; ++ s) {
    sgroup(s)   =  unidist(rng);
  }
  
  fgmm_red_bootcoef(theta_loc, stat_loc, Ze, LZV, LZy, LZZ, W_loc, sgroup,
                    ngroup, Kins, Kv, iv, overident);
  ltheta.col(k) = theta_loc;
  if (overident) {
    lstat(k)    = stat_loc;
  }
  Prog.increment();
}
#endif

Eigen::VectorXd yhat(V * theta);

// Variance of theta
Eigen::ArrayXd mtheta(ltheta.array().rowwise().mean());
Eigen::MatrixXd dtheta(ltheta.array().colwise() - mtheta);
Eigen::MatrixXd Vpa(dtheta * dtheta.transpose() / (boot - 1));

// overidentification
double pval = R_NaN;
if (overident) {
  pval = (lstat > stat).count() / static_cast<double>(boot);
}

// criterion
double cri, BIC, AIC, HQIC;
cri  = Ze.dot(W0 * Ze);
BIC  = cri - (Kins - Kx - ntau)*log(ngroup);
AIC  = cri - 2*(Kins - Kx - ntau);
HQIC = cri - 2.01*(Kins - Kx - ntau)*log(log(ngroup));
Rcpp::List critLis = Rcpp::List::create(_["criterion"] = cri,
                                        _["BIC"]       = BIC,
                                        _["AIC"]       = AIC,
                                        _["HQIC"]      = HQIC);

return Rcpp::List::create(_["parms"] = theta, _["Vpa"] = Vpa, _["Overident.stat"] = stat, _["Overident.pvalue"] = pval,
                          _["df"] = Kins - Kx - ntau, _["yhat"] = yhat, _["sigma2"] = R_NaN, _["W"] = W0,
                            _["criterion"] = critLis);
}

// Weak instrument Kp test with bootstrap
// This estimates the pi
void fKPstat_bootCoef0(Eigen::MatrixXd& selPi,
                       Eigen::MatrixXd& selZZ,
                       Eigen::MatrixXd& covz,
                       Eigen::MatrixXd& covqy,
                       const std::vector<Eigen::MatrixXd>& LZZ,
                       const std::vector<Eigen::MatrixXd>& LZqy,
                       const Eigen::MatrixXd& qy,
                       const Eigen::MatrixXd& Z,
                       const Eigen::ArrayXi& index,
                       const int& ngroup,
                       const int& Kins,
                       const int& ntau,
                       const int& boot) {
  Eigen::MatrixXd ZZ  = Eigen::MatrixXd::Zero(Kins, Kins);
  Eigen::MatrixXd Zqy = Eigen::MatrixXd::Zero(Kins, ntau);
  for (int s = 0; s < ngroup; ++s) {
    ZZ  += LZZ[s];
    Zqy += LZqy[s];
  }
  
  Eigen::MatrixXd Pi = ZZ.colPivHouseholderQr().solve(Zqy);
  selPi = Pi(index, Eigen::all).transpose();
  selZZ = ZZ(index, Eigen::all);
  
  Eigen::MatrixXd dZ  = Z.array().rowwise() - Z.array().colwise().mean();
  Eigen::MatrixXd dqy = qy.array().rowwise() - qy.array().colwise().mean();
  covz  = dZ.transpose() * dZ  / ngroup;
  covqy = dqy.transpose() * dqy / ngroup;
}



Eigen::VectorXd fKPstat_bootCoef(const std::vector<Eigen::MatrixXd>& LZZ,
                                 const std::vector<Eigen::MatrixXd>& LZqy,
                                 const Eigen::ArrayXi& index,
                                 const Eigen::ArrayXi& sgroup,
                                 const int& ngroup,
                                 const int& Kins,
                                 const int& l,
                                 const int& ntau) {
  Eigen::MatrixXd ZZ  = Eigen::MatrixXd::Zero(Kins, Kins);
  Eigen::MatrixXd Zqy = Eigen::MatrixXd::Zero(Kins, ntau);
  for (int s = 0; s < ngroup; ++s) {
    ZZ  += LZZ[sgroup(s)];
    Zqy += LZqy[sgroup(s)];
  }
  
  // Pi
  Eigen::MatrixXd Pi  = ZZ.colPivHouseholderQr().solve(Zqy);
  return Pi(index, Eigen::all).transpose().reshaped(l * ntau, 1); // pi as a vector
}


//[[Rcpp::export]]
Rcpp::List fKPstat_boot(const Eigen::MatrixXd& qy,
                        const Eigen::MatrixXd& Z,
                        const Eigen::ArrayXi& index,
                        const Eigen::ArrayXi& igroup,
                        const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                        const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                        const int& ngroup,
                        const int& boot,
                        const int& nthreads,
                        const unsigned long long seed,
                        const bool& print) {
  int ntau(qy.cols()), Kins(Z.cols()), l(index.size());
  Eigen::MatrixXd selPi, selZZ, covz, covqy;
  Eigen::MatrixXd lpi(l * ntau, boot);
  
  {
    // Compute blocks
    std::vector<Eigen::MatrixXd> LZZ(ngroup);
    std::vector<Eigen::MatrixXd> LZqy(ngroup);
    
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#pragma omp parallel for
    for (int s = 0; s < ngroup; ++s) {
      auto idx = Eigen::seq(igroup(s), igroup(s + 1) - 1);
      LZZ[s]  = Z(idx, Eigen::all).transpose() * Z(idx, Eigen::all);
      LZqy[s] = Z(idx, Eigen::all).transpose() * qy(idx, Eigen::all);
    }
#else
    for (int s = 0; s < ngroup; ++s) {
      auto idx = Eigen::seq(igroup(s), igroup(s + 1) - 1);
      LZZ[s]  = Z(idx, Eigen::all).transpose() * Z(idx, Eigen::all);
      LZqy[s] = Z(idx, Eigen::all).transpose() * qy(idx, Eigen::all);
    }
#endif
    
    // For the main sample
    Progress Prog(boot + 1, print);
    fKPstat_bootCoef0(selPi, selZZ, covz, covqy, LZZ, LZqy, qy, Z, index,
                      ngroup, Kins, ntau, boot);
    Prog.increment();
    
    //setup parallel settings
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  std::mt19937 rng(seed + tid * 7919); 
  
#pragma omp for
  for (int k = 0; k < boot; ++ k) {
    // Select subnets
    Eigen::ArrayXi sgroup(ngroup);
    std::uniform_int_distribution<int> unidist(0, ngroup - 1);
    for (int s = 0; s < ngroup; ++ s) {
      sgroup(s)   =  unidist(rng);
    }
    lpi.col(k) = fKPstat_bootCoef(LZZ, LZqy, index, sgroup, ngroup, Kins, l, ntau);
    
#pragma omp critical
    Prog.increment();
  }
}
#else
std::mt19937 rng(seed); 
for (int k = 0; k < boot; ++ k) {
  // Select subnets
  Eigen::ArrayXi sgroup(ngroup);
  std::uniform_int_distribution<int> unidist(0, ngroup - 1);
  for (int s = 0; s < ngroup; ++ s) {
    sgroup(s)   =  unidist(rng);
  }
  lpi.col(k) = fKPstat_bootCoef(LZZ, LZqy, index, sgroup, ngroup, Kins, l, ntau);
  Prog.increment();
}
#endif
  }
  
  // Variance of pi
  Eigen::ArrayXd mpi(lpi.array().rowwise().mean());
  Eigen::MatrixXd dpi(lpi.array().colwise() - mpi);
  Eigen::MatrixXd varpi(dpi * dpi.transpose() / (boot - 1));
  
  // normalisation
  Eigen::LLT<Eigen::MatrixXd> tpF(selZZ * ginv_gmm(covz) * selZZ.transpose()), tpG(covqy.inverse());
  Eigen::MatrixXd F(tpF.matrixL()); // O(1/sqrt(n))
  Eigen::MatrixXd G(tpG.matrixL().transpose()); // O(1)
  
  
  // Theta and its variance
  Eigen::MatrixXd Theta(G * selPi * F.transpose());
  Eigen::VectorXd theta(Theta.reshaped(l*ntau, 1));
  Eigen::MatrixXd FG(Eigen::kroneckerProduct(F, G));
  Eigen::MatrixXd vartheta(FG * varpi * FG.transpose());
  // cout << vartheta << endl;
  // Until this, replicate using bootstrap
  
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
  Eigen::MatrixXd BAper(Eigen::kroneckerProduct(Bper, Aper.transpose()));
  Eigen::VectorXd lambda (BAper * theta);
  Eigen::MatrixXd varlambda (BAper * vartheta * BAper.transpose());
  
  // statistic
  double stat = lambda.dot(ginv_gmm(varlambda) * lambda);
  
  return Rcpp::List::create(_["stat"] = stat, _["df"] = (ntau - q)*(l - q));
}
