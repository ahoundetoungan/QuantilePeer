// [[Rcpp::depends(RcppArmadillo, RcppNumerical, RcppEigen)]]
#include <RcppArmadillo.h>
//#define NDEBUG
#include <RcppNumerical.h>
#include <RcppEigen.h>

using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// Covariance matrix of two theta using two different sets of instruments, reduced form
//[[Rcpp::export]]
Eigen::MatrixXd Cov2ThetaRed(const Eigen::MatrixXd& V1,
                             const Eigen::MatrixXd& Z1,
                             const Eigen::MatrixXd& W1,
                             const Eigen::ArrayXd& e1,
                             const int& Kest1,
                             const Eigen::MatrixXd& V2,
                             const Eigen::MatrixXd& Z2,
                             const Eigen::MatrixXd& W2,
                             const Eigen::ArrayXd& e2,
                             const int& Kest2,
                             const int& ngroup, //common to both models
                             const Eigen::ArrayXi& cumsn, //common to both models
                             const int& HAC = 0) { //common to both models
  int Kins1(Z1.cols()), Kins2(Z2.cols()), Kins(Kins1 + Kins2), n(V1.rows()), K1(V1.cols()), K2(V2.cols()), K(K1+K2);
  Eigen::MatrixXd VZ1(V1.transpose()*Z1), VZW1(VZ1*W1), VZWZV1(VZW1*VZ1.transpose());
  Eigen::MatrixXd VZ2(V2.transpose()*Z2), VZW2(VZ2*W2), VZWZV2(VZW2*VZ2.transpose());
  Eigen::MatrixXd VZWZV(Eigen::MatrixXd::Zero(K, K)), VZW(Eigen::MatrixXd::Zero(K, Kins));
  VZWZV.block(0, 0, K1, K1)       = VZWZV1;
  VZWZV.block(K1, K1, K2, K2)     = VZWZV2;
  VZW.block(0, 0, K1, Kins1)      = VZW1;
  VZW.block(K1, Kins1, K2, Kins2) = VZW2;
  
  // variance of Ze
  Eigen::MatrixXd VZe(Eigen::MatrixXd::Zero(Kins, Kins));
  if (HAC == 0) {
    double s2((e1.square().sum() + e2.square().sum())/(n - Kest1 - Kest2));
    Eigen::MatrixXd ZZ1(Z1.transpose()*Z1), ZZ2(Z2.transpose()*Z2), Z1Z2(Z1.transpose()*Z2);
    VZe.block(0, 0, Kins1, Kins1)         = s2*ZZ1;
    VZe.block(0, Kins1, Kins1, Kins2)     = s2*Z1Z2;
    VZe.block(Kins1, 0, Kins2, Kins1)     = s2*Z1Z2.transpose();
    VZe.block(Kins1, Kins1, Kins2, Kins2) = s2*ZZ2;
    // cout<<s21<<endl;
    // cout<<s22<<endl;
    // cout<<s1s2<<endl;
  }
  if (HAC == 1) {
    Eigen::MatrixXd Ze(n, Kins1 + Kins2);
    Ze << (Z1.array().colwise()*e1).matrix(), (Z2.array().colwise()*e2).matrix();
    VZe    = Ze.transpose()*Ze;
  }
  if (HAC == 2) {
    Eigen::ArrayXXd Ze(n, Kins1 + Kins2);
    Ze << Z1.array().colwise()*e1, Z2.array().colwise()*e2;
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::VectorXd tp(Ze(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix());
      VZe += tp*tp.transpose();
    }
  }
  Eigen::MatrixXd iHdF((VZWZV).inverse()); 
  Eigen::MatrixXd HVFH(VZW*VZe*VZW.transpose());
  return iHdF * HVFH * iHdF.transpose();
}


// Covariance matrix of two theta using two different sets of instruments, Structural form
//[[Rcpp::export]]
Eigen::MatrixXd Cov2ThetaStruc(const Eigen::MatrixXd& X1,
                               const Eigen::MatrixXd& qy1,
                               const Eigen::MatrixXd& Z1,
                               const Eigen::MatrixXd& W11,
                               const Eigen::MatrixXd& W21,
                               const Eigen::ArrayXd& e1,
                               const Eigen::VectorXd theta1,
                               const int& Kest11,
                               const int& Kest21,
                               const Eigen::MatrixXd& X2,
                               const Eigen::MatrixXd& qy2,
                               const Eigen::MatrixXd& Z2,
                               const Eigen::MatrixXd& W12,
                               const Eigen::MatrixXd& W22,
                               const Eigen::ArrayXd& e2,
                               const Eigen::VectorXd theta2,
                               const int& Kest12,
                               const int& Kest22,
                               const Eigen::ArrayXi& idX1, //common to both models
                               const Eigen::ArrayXi& idX2, //common to both models
                               const Eigen::ArrayXi& nIs, //common to both models
                               const Eigen::ArrayXi& Is, //common to both models
                               const int& ngroup, //common to both models
                               const Eigen::ArrayXi& cumsn, //common to both models
                               const int& HAC = 0) {//common to both models
  int Kins1(Z1.cols()), Kins2(Z2.cols()), Kins(Kins1 + Kins2), n(X1.rows()), K1(idX1.size()), K2(idX2.size()), K(K1 + K2),
  ntau(qy1.cols()), n_iso(Is.size()), n_niso(n - n_iso);
  
  // First stage
  Eigen::MatrixXd X11(X1(Is, idX1)), XX11(X11.transpose()*X11), XXW11(XX11*W11), XXWXX11(XXW11*XX11);
  Eigen::MatrixXd X12(X2(Is, idX1)), XX12(X12.transpose()*X12), XXW12(XX12*W12), XXWXX12(XXW12*XX12);
  Eigen::VectorXd b1(theta1(1 + ntau + idX1)), b2(theta2(1 + ntau + idX1));
  
  // Second stage
  Eigen::VectorXd Xb1(X1(Eigen::all, idX1)*b1), Xb11(Xb1(Is)), Xb21(Xb1(nIs));
  Eigen::VectorXd Xb2(X2(Eigen::all, idX1)*b2), Xb12(Xb2(Is)), Xb22(Xb2(nIs));
  
  Eigen::MatrixXd X21(X1(nIs, idX2)),  X211(X1(nIs, idX1)), V21(n_niso, 1 + ntau + K2);
  V21 << Xb21, qy1(nIs, Eigen::all), X21;
  Eigen::MatrixXd X22(X2(nIs, idX2)),  X212(X2(nIs, idX1)), V22(n_niso, 1 + ntau + K2);
  V22 << Xb22, qy2(nIs, Eigen::all), X22;
  
  Eigen::MatrixXd Z21(n_niso, 1 + Kins1);
  Z21 << Xb21, Z1(nIs, Eigen::all);
  Eigen::MatrixXd Z22(n_niso, 1 + Kins2);
  Z22 << Xb22, Z2(nIs, Eigen::all);
  
  Eigen::MatrixXd ZV21(Z21.transpose()*V21), ZZ21(Z21.transpose()*Z21);
  Eigen::MatrixXd ZV22(Z22.transpose()*V22), ZZ22(Z22.transpose()*Z22);
  
  Eigen::MatrixXd VZW21(ZV21.transpose()*W21), VZWZV21(VZW21*ZV21);
  Eigen::MatrixXd VZW22(ZV22.transpose()*W22), VZWZV22(VZW22*ZV22);
  
  Eigen::ArrayXd e11(e1(Is)), e21(e1(nIs)*(1 - theta1(0)));
  Eigen::ArrayXd e12(e2(Is)), e22(e2(nIs)*(1 - theta2(0)));
  
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(2*K + 2*ntau + 2, 2*K1 + Kins + 2));
  H.block(0, 0, K1, K1)  = XXW11;
  H.block(K1, K1, K2 + ntau + 1, Kins1 + 1) = VZW21;
  H.block(K + ntau + 1, K1 + Kins1 + 1, K1, K1)  = XXW12;
  H.block(K + K1 + ntau + 1, 2*K1 + Kins1 + 1, K2 + ntau + 1, Kins2 + 1) = VZW22;
  
  Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(2*K1 + Kins + 2, 2*K + 2*ntau + 2));
  dF.block(0, 0, K1, K1) = XX11;
  dF.block(K1, 0, Kins1 + 1, K + ntau + 1) << (Z21.transpose()*X211*(1 - theta1(0))), ZV21;
  dF.block(K1 + Kins1 + 1, K + ntau + 1, K1, K1) = XX12;
  dF.block(2*K1 + Kins1 + 1, K + ntau + 1, Kins2 + 1, K + ntau + 1) << (Z22.transpose()*X212*(1 - theta2(0))), ZV22;
  
  Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(2*K1 + Kins + 2, 2*K1 + Kins + 2));
  // cout<<K1<<endl;
  // cout<<Kins<<endl;
  if (HAC == 0) {
    double s21((e11.square().sum() + e12.square().sum())/(n_iso - Kest11 - Kest12));
    double s22((e21.square().sum() + e22.square().sum())/(n_niso - Kest21 - Kest22));
    Eigen::MatrixXd X11X12(X11.transpose()*X12), Z21Z22(Z21.transpose()*Z22);
    VF.block(0, 0, K1, K1) = s21*XX11;
    VF.block(0, K1 + Kins1 + 1, K1, K1) = s21*X11X12;
    VF.block(K1, K1, Kins1 + 1, Kins1 + 1) = s22*ZZ21;
    VF.block(K1, 2*K1 + Kins1 + 1, Kins1 + 1, Kins2 + 1) = s22*Z21Z22;
    VF.block(K1 + Kins1 + 1, 0, K1, K1) = s21*X11X12.transpose();
    VF.block(K1 + Kins1 + 1, K1 + Kins1 + 1, K1, K1) = s21*XX12;
    VF.block(2*K1 + Kins1 + 1, K1, Kins2 + 1, Kins1 + 1) = s22*Z21Z22.transpose();
    VF.block(2*K1 + Kins1 + 1, 2*K1 + Kins1 + 1, Kins2 + 1, Kins2 + 1) = s22*ZZ22;
    // cout<<s211<<endl;
    // cout<<s221<<endl;
    // cout<<s212<<endl;
    // cout<<s222<<endl;
    // cout<<s1s2iso<<endl;
    // cout<<s1s2niso<<endl;
    // cout<<VF<<endl;
    // cout << VF.eigenvalues().array() <<endl;
  }
  if (HAC == 1) {
    Eigen::MatrixXd Xe11(X11.array().colwise()*e11), Ze21(Z21.array().colwise()*e21);
    Eigen::MatrixXd Xe12(X12.array().colwise()*e12), Ze22(Z22.array().colwise()*e22);
    Eigen::MatrixXd Xe11Xe12(Xe11.transpose()*Xe12), Ze21Ze22(Ze21.transpose()*Ze22);
    VF.block(0, 0, K1, K1) = Xe11.transpose()*Xe11;
    VF.block(0, K1 + Kins1 + 1, K1, K1) = Xe11Xe12;
    VF.block(K1, K1, Kins1 + 1, Kins1 + 1) = Ze21.transpose()*Ze21;
    VF.block(K1, 2*K1 + Kins1 + 1, Kins1 + 1, Kins2 + 1) = Ze21Ze22;
    VF.block(K1 + Kins1 + 1, 0, K1, K1) = Xe11Xe12.transpose();
    VF.block(K1 + Kins1 + 1, K1 + Kins1 + 1, K1, K1) = Xe12.transpose()*Xe12;
    VF.block(2*K1 + Kins1 + 1, K1, Kins2 + 1, Kins1 + 1) = Ze21Ze22.transpose();
    VF.block(2*K1 + Kins1 + 1, 2*K1 + Kins1 + 1, Kins2 + 1, Kins2 + 1) = Ze22.transpose()*Ze22;
  }
  if (HAC == 2) {
    X11   = Eigen::MatrixXd::Zero(n, K1);
    X11(Is, Eigen::all) = X1(Is, idX1);
    X12   = Eigen::MatrixXd::Zero(n, K1);
    X12(Is, Eigen::all) = X2(Is, idX1);
    Z21   = Eigen::MatrixXd::Zero(n, 1 + Kins1);
    Z21(nIs, Eigen::all) << Xb21, Z1(nIs, Eigen::all);
    Z22   = Eigen::MatrixXd::Zero(n, 1 + Kins2);
    Z22(nIs, Eigen::all) << Xb22, Z2(nIs, Eigen::all);
    Eigen::VectorXd est1(Eigen::VectorXd::Zero(n)), est2(Eigen::VectorXd::Zero(n));
    est1(Is)  = e11; est1(nIs) = e21;
    est2(Is)  = e12; est2(nIs) = e22;
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::MatrixXd tp11(n2 - n1 + 1, K1 + Kins1 + 1), tp12(n2 - n1 + 1, K1 + Kins2 + 1);
      tp11 << X11(Eigen::seq(n1, n2), Eigen::all), Z21(Eigen::seq(n1, n2), Eigen::all);
      tp12 << X12(Eigen::seq(n1, n2), Eigen::all), Z22(Eigen::seq(n1, n2), Eigen::all);
      Eigen::VectorXd tp2(2*K1 + Kins + 2);
      tp2 << tp11.transpose()*est1.segment(n1, n2), tp12.transpose()*est2.segment(n1, n2);
      VF += tp2*tp2.transpose();
    }
  }
  Eigen::MatrixXd iHdF((H*dF).inverse()); 
  Eigen::MatrixXd HVFH(H*VF*H.transpose());
  Eigen::MatrixXd Vpa(iHdF * HVFH * iHdF.transpose());
  
  // Find correct order
  Eigen::MatrixXd R1(Eigen::MatrixXd::Zero(1 + ntau + K, 1 + ntau + K));
  Eigen::MatrixXd R2(Eigen::MatrixXd::Zero(1 + ntau + K, 1 + ntau + K));
  R1(0, K1) = -1; R2(0, K1) = -1;
  R1.block(1, K1 + 1, ntau, ntau) = Eigen::MatrixXd::Identity(ntau, ntau);
  R2.block(1, K1 + 1, ntau, ntau) = Eigen::MatrixXd::Identity(ntau, ntau);
  R1(1 + ntau + idX1, Eigen::seqN(0, K1)) = Eigen::MatrixXd::Identity(K1, K1);
  R2(1 + ntau + idX1, Eigen::seqN(0, K1)) = Eigen::MatrixXd::Identity(K1, K1);
  if (K2 > 0) {
    R1.block(0, K1 + ntau + 1, 1, K2) = -theta1(1 + ntau + idX2)/(1 - theta1(0));
    R2.block(0, K1 + ntau + 1, 1, K2) = -theta2(1 + ntau + idX2)/(1 - theta2(0));
    R1(1 + ntau + idX2, Eigen::seqN(1 + ntau + K1, K2)) = Eigen::MatrixXd::Identity(K2, K2)/(1 - theta1(0));
    R2(1 + ntau + idX2, Eigen::seqN(1 + ntau + K1, K2)) = Eigen::MatrixXd::Identity(K2, K2)/(1 - theta2(0));
  }
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(2 + 2*ntau + 2*K, 2 + 2*ntau + 2*K));
  R.block(0, 0, 1 + ntau + K, 1 + ntau + K) = R1;
  R.block(1 + ntau + K, 1 + ntau + K, 1 + ntau + K, 1 + ntau + K) = R2;
  return R * Vpa * R.transpose();
}

// Optimization
class OptimTest: public MFuncGrad
{
private:
  const Eigen::VectorXd& thetahat;
  const Eigen::MatrixXd& Sigma;
  const Eigen::VectorXd& a;
  const int K;
public:
  OptimTest(const Eigen::VectorXd& thetahat_, 
            const Eigen::MatrixXd& Sigma_,
            const Eigen::VectorXd& a_,
            const int K_) : 
  thetahat(thetahat_), 
  Sigma(Sigma_),
  a(a_),
  K(K_) {}
  double f_grad(Constvec& beta, Refvec grad)
  {
    Eigen::VectorXd theta(beta);
    Eigen::MatrixXd tp1(Eigen::MatrixXd::Identity(K, K));
    for (int i(0); i < (K - 1); ++i) {
      theta(i + 1)  = theta(i) + a(i)*exp(beta(i + 1));
      tp1(i + 1, Eigen::seqN(0, i + 2)) << 1, (a(Eigen::seqN(0, i + 1)).array()*beta.segment(1, i + 1).array().exp()).matrix().transpose();
    }
    theta = theta.array().min(1e3).max(-1e3).matrix();
    // cout<<theta.transpose()<<endl;
    Eigen::VectorXd tp2(Sigma.colPivHouseholderQr().solve(theta - thetahat));
    double f(tp2.dot(theta - thetahat));
    grad = (2 * tp1.transpose() * tp2);
    // cout<<f<<endl;
    return f;
  }
};

Rcpp::List fOptimTest(const Eigen::VectorXd& thetahat, 
                      const Eigen::MatrixXd& Sigma,
                      const Eigen::VectorXd& a,
                      const int& K,
                      const int &maxit,
                      const double& eps_f, 
                      const double& eps_g){
  OptimTest f(thetahat, Sigma, a, K);
  Eigen::VectorXd beta(Eigen::VectorXd::Zero(K));
  double fopt;
  int status = optim_lbfgs(f, beta, fopt, maxit, eps_f, eps_g);
  Eigen::VectorXd theta(beta);
  for (int i(0); i < (K - 1); ++i) {
    theta(i + 1)  = theta(i) + a(i)*exp(beta(i + 1));
  }
  return Rcpp::List::create(_["lambda"] = theta, _["minimum"] = fopt, _["status"] = status);
}


// [[Rcpp::export]]
Rcpp::List fTestMonotone(const Eigen::VectorXd& thetahat,
                         const Eigen::MatrixXd& Sigma,
                         const Eigen::VectorXd& a,
                         const Eigen::MatrixXd& thetasimu,
                         const int& Boot = 1e4,
                         const int& maxit = 1e6,
                         const double& eps_f = 1e-9,
                         const double& eps_g = 1e-9){
  int K(thetahat.size());
  Rcpp::List tpopt = fOptimTest(thetahat, Sigma, a, K, maxit, eps_f, eps_g);
  
  // Bootstrap
  // arma::mat Listtheta(arma::randn(K, Boot)), Sigmamat(Sigma.data(), K, K, /*copy_aux_mem=*/false, /*strict=*/true);
  // arma::vec thetahatvec(thetahat.data(), K, /*copy_aux_mem=*/false, /*strict=*/true);
  // Listtheta = arma::chol(Sigmamat).t()*Listtheta*0;
  // Listtheta.each_col() += thetahatvec;
  // Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> thetasimu(Listtheta.memptr(), K, Boot);
  Eigen::ArrayXi countPos(Boot);
  for (int b(0); b < Boot; ++ b) {
    // cout<<b<<endl;
    Rcpp::List tp = fOptimTest(thetasimu.col(b), Sigma, a, K, maxit, eps_f, eps_g);
    Eigen::ArrayXd thetab = tp["lambda"];
    Eigen::ArrayXi tp1((((thetab.segment(1, K - 1) - thetab.segment(0, K - 1)) * a.array())
                         > (1e5 * (eps_f + eps_g))).cast<int>());
    // cout<<thetab.transpose()<<endl;
    // cout<<tp1.transpose()<<endl;
    // cout<<tp1.sum()<<endl;
    countPos(b) = tp1.sum();
  }
  return Rcpp::List::create(_["optim"] = tpopt, _["count"] = countPos);
}