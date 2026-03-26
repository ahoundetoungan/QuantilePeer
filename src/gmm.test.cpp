// [[Rcpp::depends(RcppArmadillo, RcppProgress, RcppNumerical, RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
//#define NDEBUG
#include <random>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <RcppNumerical.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/KroneckerProduct>

#if defined(_OPENMP)
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

// generalized inverse
Eigen::MatrixXd ginv(const Eigen::MatrixXd& A) {
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
Eigen::MatrixXd matrixSqrt2(const Eigen::MatrixXd& A) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  Eigen::VectorXd sqrt_evals = es.eigenvalues().array().sqrt();
  return es.eigenvectors() * sqrt_evals.asDiagonal() * es.eigenvectors().transpose();
}

// Covariance matrix of two theta using two different sets of instruments, reduced form
//[[Rcpp::export]]
Rcpp::List Cov2ThetaRed(const Eigen::MatrixXd& Z1,
                        const Eigen::MatrixXd& W1,
                        const Eigen::ArrayXd& e1,
                        const Eigen::VectorXd& theta1,
                        const Eigen::MatrixXd& Z2,
                        const Eigen::MatrixXd& W2,
                        const Eigen::ArrayXd& e2,
                        const Eigen::VectorXd theta2,
                        const Eigen::MatrixXd& X, //common to both models
                        const Eigen::MatrixXd& qy,//common to both models
                        const int& Kest, //common to both models
                        const int& ngroup, //common to both models
                        const Eigen::ArrayXi& cumsn, //common to both models
                        const int& HAC = 0, //common to both models
                        const bool& full = false) { 
  int ntau(qy.cols()), Kins1(Z1.cols()), Kins2(Z2.cols()), 
  Kins(Kins1 + Kins2), n(X.rows()), K(X.cols()), Kv(K + ntau);
  
  Eigen::MatrixXd V(n, Kv);
  V << qy, X;
  
  Eigen::MatrixXd VZ1(V.transpose()*Z1), VZW1(VZ1*W1), VZWZV1(VZW1*VZ1.transpose());
  Eigen::MatrixXd VZ2(V.transpose()*Z2), VZW2(VZ2*W2), VZWZV2(VZW2*VZ2.transpose());
  Eigen::MatrixXd VZWZV(Eigen::MatrixXd::Zero(2*Kv, 2*Kv)), VZW(Eigen::MatrixXd::Zero(2*Kv, Kins));
  VZWZV.block(0, 0, Kv, Kv)        = VZWZV1;
  VZWZV.block(Kv, Kv, Kv, Kv)      = VZWZV2;
  VZW.block(0, 0, Kv, Kins1)       = VZW1;
  VZW.block(Kv, Kins1, Kv, Kins2)  = VZW2;
  
  // variance of Ze
  Eigen::MatrixXd VZe(Eigen::MatrixXd::Zero(Kins, Kins));
  if (HAC == 0) {
    double s2((e1.square().sum() + e2.square().sum())/(2*n - 2*Kest));
    Eigen::MatrixXd ZZ1(Z1.transpose()*Z1), ZZ2(Z2.transpose()*Z2), Z1Z2(Z1.transpose()*Z2);
    VZe.block(0, 0, Kins1, Kins1)         = s2*ZZ1;
    VZe.block(0, Kins1, Kins1, Kins2)     = s2*Z1Z2;
    VZe.block(Kins1, 0, Kins2, Kins1)     = s2*Z1Z2.transpose();
    VZe.block(Kins1, Kins1, Kins2, Kins2) = s2*ZZ2;
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
  
  // theta 
  Eigen::VectorXd theta(2*Kv);
  theta << theta1, theta2;
  
  // test all parameters of just quantiles
  int df(ntau);
  if (full) { 
    df = Kv; 
  }
  
  // R matrix
  Eigen::ArrayXi itheta(df);
  itheta << Eigen::ArrayXi::LinSpaced(df, 1, df);
  
  // R for selected parameters
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(df, 2*Kv));
  R(Eigen::all, Eigen::seqN(0, df))  = Eigen::MatrixXd::Identity(df, df);
  R(Eigen::all, Eigen::seqN(Kv, df)) = -Eigen::MatrixXd::Identity(df, df);
  
  // R for full 
  Eigen::MatrixXd Rfull(Eigen::MatrixXd::Zero(Kv, 2*Kv));
  Rfull(Eigen::all, Eigen::seqN(0, Kv))  = Eigen::MatrixXd::Identity(Kv, Kv);
  Rfull(Eigen::all, Eigen::seqN(Kv, Kv)) = -Eigen::MatrixXd::Identity(Kv, Kv);
  
  // test statistic
  Eigen::VectorXd Rtheta(R*theta);
  Eigen::MatrixXd RiHdF(R*iHdF);
  Eigen::MatrixXd VarRtheta(RiHdF * HVFH * RiHdF.transpose());
  Eigen::VectorXd stat(Rtheta.transpose() * ginv(VarRtheta) * Rtheta);
  
  Eigen::VectorXd Rfulltheta(Rfull*theta);
  Eigen::MatrixXd RfulliHdF(Rfull*iHdF);
  Eigen::MatrixXd VarRfulltheta(RfulliHdF * HVFH * RfulliHdF.transpose());
  
  return Rcpp::List::create(_["stat"]    = stat, 
                            _["df"]      = df, 
                            _["dtheta"]  = Rfulltheta,
                            _["Vdtheta"] = VarRfulltheta,
                            _["itheta"]  = itheta);
}


// Covariance matrix of two theta using two different sets of instruments, Structural form
//[[Rcpp::export]]
Rcpp::List Cov2ThetaStruc(const Eigen::MatrixXd& Z1,
                          const Eigen::MatrixXd& W21,
                          const Eigen::ArrayXd& e1,
                          const Eigen::VectorXd& theta1,
                          const Eigen::MatrixXd& Z2,
                          const Eigen::MatrixXd& W22,
                          const Eigen::ArrayXd& e2,
                          const Eigen::VectorXd& theta2,
                          const Eigen::MatrixXd& X, //common to both models
                          const Eigen::MatrixXd& qy, //common to both models
                          const Eigen::MatrixXd& W1, //common to both models
                          const int& Kest1, //common to both models
                          const int& Kest2, //common to both models
                          const Eigen::ArrayXi& idX1, //common to both models
                          const Eigen::ArrayXi& idX2, //common to both models
                          const Eigen::ArrayXi& nIs, //common to both models
                          const Eigen::ArrayXi& Is, //common to both models
                          const int& ngroup, //common to both models
                          const Eigen::ArrayXi& cumsn, //common to both models
                          const int& HAC = 0, //common to both models
                          const bool& full = false) {
  int Kins1(Z1.cols()), Kins2(Z2.cols()), Kins(Kins1 + Kins2), n(X.rows()), 
  K1(idX1.size()), K2(idX2.size()), K(K1 + K2),
  ntau(qy.cols()), n_iso(Is.size()), n_niso(n - n_iso);
  
  // First stage
  Eigen::MatrixXd X1(X(Is, idX1)), XX1(X1.transpose()*X1), XXW1(XX1*W1), XXWXX1(XXW1*XX1);
  Eigen::VectorXd b(theta1(1 + ntau + idX1));// b2(theta2(1 + ntau + idX1));
  
  // Second stage
  Eigen::MatrixXd X21(X(nIs, idX1)), X22(X(nIs, idX2)), V2(n_niso, 1 + ntau + K2);
  V2 << X(nIs, idX1)*b, qy(nIs, Eigen::all), X22;
  
  Eigen::MatrixXd Z21(Z1(nIs, Eigen::all));
  Eigen::MatrixXd Z22(Z2(nIs, Eigen::all));
  
  Eigen::MatrixXd ZV21(Z21.transpose()*V2), ZZ21(Z21.transpose()*Z21);
  Eigen::MatrixXd ZV22(Z22.transpose()*V2), ZZ22(Z22.transpose()*Z22);
  
  Eigen::MatrixXd VZW21(ZV21.transpose()*W21), VZWZV21(VZW21*ZV21);
  Eigen::MatrixXd VZW22(ZV22.transpose()*W22), VZWZV22(VZW22*ZV22);
  
  Eigen::ArrayXd eiso(e1(Is)), e21(e1(nIs)*(1 - theta1(0))), e22(e2(nIs)*(1 - theta2(0)));
  
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(K + K2 + 2*ntau + 2, K1 + Kins));
  H.block(0, 0, K1, K1)  = XXW1;
  H.block(K1, K1, K2 + ntau + 1, Kins1) = VZW21;
  H.block(K + ntau + 1, K1 + Kins1, K2 + ntau + 1, Kins2) = VZW22;
  
  
  Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(K1 + Kins, K + K2 + 2*ntau + 2));
  dF.block(0, 0, K1, K1) = XX1;
  dF.block(K1, 0, Kins1, K + ntau + 1) << Z21.transpose()*X21*(1 - theta1(0)), ZV21;
  dF.block(K1 + Kins1, 0, Kins2, K1) << (Z22.transpose()*X21*(1 - theta2(0)));
  dF.block(K1 + Kins1, K + ntau + 1, Kins2, K2 + ntau + 1) << ZV22;
  
  Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(K1 + Kins, K1 + Kins));
  if (HAC == 0) {
    double s21(eiso.square().sum()/(n_iso - Kest1));
    double s22((e21.square().sum() + e22.square().sum())/(2*n_niso - 2*Kest2));
    Eigen::MatrixXd Z21Z22(Z21.transpose()*Z22);
    VF.block(0, 0, K1, K1) = s21*XX1;
    VF.block(K1, K1, Kins1, Kins1) = s22*ZZ21;
    VF.block(K1, K1 + Kins1, Kins1, Kins2) = s22*Z21Z22;
    VF.block(K1 + Kins1, K1, Kins2, Kins1) = s22*Z21Z22.transpose();
    VF.block(K1 + Kins1, K1 + Kins1, Kins2, Kins2) = s22*ZZ22;
  }
  if (HAC == 1) {
    Eigen::MatrixXd Xe1(X1.array().colwise()*eiso), Ze21(Z21.array().colwise()*e21);
    Eigen::MatrixXd Ze22(Z22.array().colwise()*e22), Ze21Ze22(Ze21.transpose()*Ze22);
    VF.block(0, 0, K1, K1) = Xe1.transpose()*Xe1;
    VF.block(K1, K1, Kins1, Kins1) = Ze21.transpose()*Ze21;
    VF.block(K1, K1 + Kins1, Kins1, Kins2) = Ze21Ze22;
    VF.block(K1 + Kins1, K1, Kins2, Kins1) = Ze21Ze22.transpose();
    VF.block(K1 + Kins1, K1 + Kins1, Kins2, Kins2) = Ze22.transpose()*Ze22;
  }
  if (HAC == 2) {
    X1   = Eigen::MatrixXd::Zero(n, K1);
    X1(Is, Eigen::all) = X(Is, idX1);
    Z21   = Eigen::MatrixXd::Zero(n, Kins1);
    Z21(nIs, Eigen::all) = Z1(nIs, Eigen::all);
    Z22   = Eigen::MatrixXd::Zero(n, Kins2);
    Z22(nIs, Eigen::all) = Z2(nIs, Eigen::all);
    Eigen::VectorXd eps1(Eigen::VectorXd::Zero(n)), eps2(Eigen::VectorXd::Zero(n));
    eps1(Is)  = eiso; eps1(nIs) = e21;
    eps2(Is)  = eiso; eps2(nIs) = e22;
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::VectorXd tp(K1 + Kins);
      tp << X1(Eigen::seq(n1, n2), Eigen::all).transpose() * eps1.segment(n1, n2), 
            Z21(Eigen::seq(n1, n2), Eigen::all).transpose() * eps1.segment(n1, n2),
            Z22(Eigen::seq(n1, n2), Eigen::all).transpose() * eps2.segment(n1, n2);
      VF += tp * tp.transpose();
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
    R1(1 + ntau + idX2, K1) = -theta1(1 + ntau + idX2)/(1 - theta1(0));
    R2(1 + ntau + idX2, K1) = -theta2(1 + ntau + idX2)/(1 - theta2(0));
    R1(1 + ntau + idX2, Eigen::seqN(1 + ntau + K1, K2)) = Eigen::MatrixXd::Identity(K2, K2)/(1 - theta1(0));
    R2(1 + ntau + idX2, Eigen::seqN(1 + ntau + K1, K2)) = Eigen::MatrixXd::Identity(K2, K2)/(1 - theta2(0));
  }
  
  Eigen::ArrayXi seltp(1 + ntau + K2);
  seltp << Eigen::ArrayXi::LinSpaced(1 + ntau, 0, ntau), 1 + ntau + idX2;
  
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(2 + 2*ntau + 2*K2, 2 + 2*ntau + K + K2));
  R.block(0, 0, 1 + ntau + K2, 1 + ntau + K) = R1(seltp, Eigen::all);
  R.block(1 + ntau + K2, 1 + ntau + K, 1 + ntau + K2, 1 + ntau + K2) = 
    R2(seltp, Eigen::seqN(K1, 1 + ntau + K2));
  
  // theta 
  Eigen::VectorXd theta(2*K2 + 2*ntau + 2);
  theta << theta1(Eigen::seqN(0, 1 + ntau)), theta1(1 + ntau + idX2), 
           theta2(Eigen::seqN(0, 1 + ntau)), theta2(1 + ntau + idX2);
  
  // test all parameters of just quantiles
  int df(ntau);
  if (full) { 
    df = 1 + ntau + K2; 
  }
  
  // R matrix for the test
  Eigen::ArrayXi itheta(df);
  Eigen::MatrixXd Rt(Eigen::MatrixXd::Zero(df, 2 + 2*ntau + 2*K2));
  if (full) { 
    Rt <<  Eigen::MatrixXd::Identity(df, df), -Eigen::MatrixXd::Identity(df, df);
    itheta << Eigen::ArrayXi::LinSpaced(1 + ntau, 1, 1 + ntau),  2 + ntau + idX2;
  } else {
    Rt(Eigen::all, Eigen::seqN(1, ntau)) = Eigen::MatrixXd::Identity(df, df);
    Rt(Eigen::all, Eigen::seqN(2 + ntau + K2, ntau)) = -Eigen::MatrixXd::Identity(df, df);
    itheta << Eigen::ArrayXi::LinSpaced(ntau, 2, 1 + ntau);
  }
  
  // R for full 
  Eigen::MatrixXd Rfull(Eigen::MatrixXd::Zero(1 + ntau + K2, 2 + 2*ntau + 2*K2));
  Rfull <<  Eigen::MatrixXd::Identity(1 + ntau + K2, 1 + ntau + K2), 
            -Eigen::MatrixXd::Identity(1 + ntau + K2, 1 + ntau + K2);
  
  // test statistic
  Eigen::VectorXd Rtheta(Rt * theta);
  Eigen::MatrixXd RtR(Rt * R);
  Eigen::MatrixXd VarRtheta(RtR * Vpa * RtR.transpose());
  Eigen::VectorXd stat(Rtheta.transpose() * ginv(VarRtheta) * Rtheta);
  
  
  Eigen::VectorXd Rfulltheta(Rfull * theta);
  Eigen::MatrixXd RfullR(Rfull * R);
  Eigen::MatrixXd VarRfulltheta(RfullR * Vpa * RfullR.transpose());
  
  return Rcpp::List::create(_["stat"]    = stat, 
                            _["df"]      = df, 
                            _["dtheta"]  = Rfulltheta,
                            _["Vdtheta"] = VarRfulltheta,
                            _["itheta"]  = itheta);
}

// This function implement test for validity of Z2 using the fact that Z1 is valid
//[[Rcpp::export]]
Rcpp::List validZ2SarganRed(const Eigen::MatrixXd& Z1,
                            const Eigen::MatrixXd& W1,
                            const Eigen::ArrayXd& e1,
                            const Eigen::VectorXd& theta1,
                            const Eigen::MatrixXd& Z2,
                            const Eigen::MatrixXd& W2,
                            const Eigen::ArrayXd& e2,
                            const Eigen::VectorXd& theta2,
                            const Eigen::MatrixXd& X, //common to both models
                            const Eigen::MatrixXd& qy,//common to both models
                            const int& Kest, //common to both models
                            const int& ngroup, //common to both models
                            const Eigen::ArrayXi& cumsn, //common to both models
                            const int& HAC = 0,
                            const bool& full = false) {
  int ntau(qy.cols()), Kins2(Z2.cols()), n(X.rows());
  Eigen::MatrixXd H(Z2 - Z1*(Z1.transpose() * Z1).colPivHouseholderQr().solve(Z1.transpose() * Z2));
  Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
  int df(lu.rank());
  
  // variance of He
  Eigen::MatrixXd VHe(Eigen::MatrixXd::Zero(Kins2, Kins2));
  if (HAC == 0) {
    double s2(e1.square().sum()/(n - Kest));
    VHe = s2 * (H.transpose() * H);
  }
  if (HAC == 1) {
    Eigen::MatrixXd He(n, Kins2);
    He << (H.array().colwise()*e1).matrix();
    VHe    = He.transpose()*He;
  }
  if (HAC == 2) {
    Eigen::ArrayXXd He(n, Kins2);
    He << H.array().colwise()*e1;
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::VectorXd tp(He(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix());
      VHe += tp*tp.transpose();
    }
  }
  
  Eigen::VectorXd He(H.transpose() * e1.matrix());
  double stat((He.transpose() * ginv(VHe) * He)(0, 0));
  return Rcpp::List::create(_["stat"]   = stat, 
                            _["df"]     = df, 
                            _["theta1"] = theta1.head(ntau),
                            _["theta2"] = theta2.head(ntau));
}


//[[Rcpp::export]]
Rcpp::List validZ2SarganStruc(const Eigen::MatrixXd& Z1,
                              const Eigen::MatrixXd& W21,
                              const Eigen::ArrayXd& e1,
                              const Eigen::VectorXd theta1,
                              const Eigen::MatrixXd& Z2,
                              const Eigen::MatrixXd& W22,
                              const Eigen::ArrayXd& e2,
                              const Eigen::VectorXd& theta2,
                              const Eigen::MatrixXd& X, //common to both models
                              const Eigen::MatrixXd& qy, //common to both models
                              const Eigen::MatrixXd& W1, //common to both models
                              const int& Kest1, //common to both models
                              const int& Kest2, //common to both models
                              const Eigen::ArrayXi& idX1, //common to both models
                              const Eigen::ArrayXi& idX2, //common to both models
                              const Eigen::ArrayXi& nIs, //common to both models
                              const Eigen::ArrayXi& Is, //common to both models
                              const int& ngroup, //common to both models
                              const Eigen::ArrayXi& cumsn, //common to both models
                              const int& HAC = 0,
                              const bool& full = false) {
  int Kins2(Z2.cols()), n(X.rows()), ntau(qy.cols()), n_iso(Is.size()), n_niso(n - n_iso);
  
  // First stage
  Eigen::VectorXd b(theta1(1 + ntau + idX1));
  
  // Instruments
  Eigen::MatrixXd Z21(Z1(nIs, Eigen::all));
  Eigen::MatrixXd Z22(Z2(nIs, Eigen::all));
  
  Eigen::MatrixXd H(Z22 - Z21*(Z21.transpose() * Z21).colPivHouseholderQr().solve(Z21.transpose() * Z22));
  Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
  int df(lu.rank());
  
  // variance of He
  Eigen::MatrixXd VHe(Eigen::MatrixXd::Zero(Kins2, Kins2));
  if (HAC == 0) {
    double s2(e1(nIs).square().sum()/(n_niso - Kest2));
    VHe = s2 * (H.transpose() * H);
  }
  if (HAC == 1) {
    Eigen::MatrixXd He((H.array().colwise() * e1(nIs)).matrix());
    VHe    = He.transpose()*He;
  }
  if (HAC == 2) {
    Eigen::ArrayXXd He(Eigen::MatrixXd::Zero(n, Kins2));
    He(nIs, Eigen::all) << H.array().colwise() * e1(nIs);
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::VectorXd tp(He(Eigen::seq(n1, n2), Eigen::all).colwise().sum().matrix());
      VHe += tp*tp.transpose();
    }
  }
  
  Eigen::VectorXd He(H.transpose() * e1(nIs).matrix());
  double stat((He.transpose() * ginv(VHe) * He)(0, 0));
  return Rcpp::List::create(_["stat"]   = stat, 
                            _["df"]     = df, 
                            _["theta1"] = theta1(Eigen::seqN(1, ntau)),
                            _["theta2"] = theta2(Eigen::seqN(1, ntau)));
}


// Optimization for monotonicity test
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
                         const int& Boot,
                         const int& maxit,
                         const double& eps_f,
                         const double& eps_g){
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
    countPos(b) = tp1.sum();
  }
  return Rcpp::List::create(_["optim"] = tpopt, _["count"] = countPos);
}


// Encompassing test
// Estimation one delta given selected groups
Eigen::VectorXd fEncompassingStrucDeltaCoef(const Eigen::VectorXd& y,
                                            const Eigen::MatrixXd& qy1,
                                            const Eigen::MatrixXd& Z1,
                                            const Eigen::MatrixXd& qy2,
                                            const Eigen::MatrixXd& Z2,
                                            const Eigen::MatrixXd& X, //common to both models
                                            const Eigen::ArrayXi& idX1, //common to both models
                                            const Eigen::ArrayXi& idX2, //common to both models
                                            const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                                            const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                                            const int& ngroup,
                                            const int& ntau1,
                                            const int& ntau2,
                                            const int& K2,
                                            const int& Kins1,
                                            const int& Kins2,
                                            const int& Kv1,
                                            const int& Kv2,
                                            const bool& iv1,
                                            const bool& iv2){
  // Find iso and nisolated
  int n_iso(0), n_niso(0);
  for (int s = 0; s < ngroup; ++s) {
    n_iso  += LIs[s].size();
    n_niso += LnIs[s].size();
  }
  
  Eigen::ArrayXi Is(n_iso), nIs(n_niso);
  int pos_iso(0), pos_niso(0);
  for (int s = 0; s < ngroup; ++s) {
    int n1(LIs[s].size());
    if (n1 > 0) {
      Is.segment(pos_iso, n1) = LIs[s];
      pos_iso += n1;
    }
    
    int n2(LnIs[s].size());
    if (n2 > 0) {
      nIs.segment(pos_niso, n2) = LnIs[s];
      pos_niso += n2;
    }
  }
  
  // First stage
  Eigen::VectorXd b(X(Is, idX1).colPivHouseholderQr().solve(y(Is)));
  
  // Variables and Instruments
  Eigen::VectorXd Xb2(X(nIs, idX1) * b);
  Eigen::MatrixXd V21(n_niso, Kv1), V22(n_niso, Kv2);
  V21 << Xb2, qy1(nIs, Eigen::all), X(nIs, idX2);
  V22 << Xb2, qy2(nIs, Eigen::all), X(nIs, idX2);
  
  Eigen::MatrixXd Z21(Z1(nIs, Eigen::all));
  Eigen::MatrixXd Z22(Z2(nIs, Eigen::all));
  
  Eigen::MatrixXd ZV21(Z21.transpose()*V21);
  Eigen::MatrixXd ZV22(Z22.transpose()*V22);
  
  // Weights
  Eigen::MatrixXd W21 = Eigen::MatrixXd::Identity(Kins1, Kins1);
  Eigen::MatrixXd W22 = Eigen::MatrixXd::Identity(Kins2, Kins2);
  if (iv1) {
    W21 = (Z21.transpose()*Z21).inverse();
  }
  if (iv2) {
    W22 = (Z22.transpose()*Z22).inverse();
  }
  
  Eigen::VectorXd y2(y(nIs));
  Eigen::MatrixXd VZW21(ZV21.transpose()*W21), VZWZV21(VZW21*ZV21);
  Eigen::MatrixXd VZW22(ZV22.transpose()*W22), VZWZV22(VZW22*ZV22);
  Eigen::VectorXd lambda1(VZWZV21.colPivHouseholderQr().solve(VZW21 * Z21.transpose() * y2));
  Eigen::ArrayXd e21(y2 - V21 * lambda1);
  
  return (VZWZV22).colPivHouseholderQr().solve(VZW22 * Z22.transpose() * e21.matrix());
}


// Estimating delta and its variance using bootstreap
Rcpp::List fEncompassingStrucDelta(const Eigen::VectorXd& y,
                                   const Eigen::MatrixXd& qy1,
                                   const Eigen::MatrixXd& Z1,
                                   const Eigen::MatrixXd& qy2,
                                   const Eigen::MatrixXd& Z2,
                                   const Eigen::MatrixXd& X, //common to both models
                                   const Eigen::ArrayXi& idX1, //common to both models
                                   const Eigen::ArrayXi& idX2, //common to both models
                                   const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                                   const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                                   const int& ngroup,
                                   const bool& iv1,
                                   const bool& iv2, //common to both models
                                   const int& boot,
                                   const int& nthreads,
                                   const unsigned long long seed,
                                   const bool& print) { 
  // Fixed variables
  int ntau1(qy1.cols()), ntau2(qy2.cols()), K2(idX2.size()), 
  Kins1(Z1.cols()), Kins2(Z2.cols()),
  Kv1(1 + ntau1 + K2), Kv2(1 + ntau2 + K2);
  
  Progress Prog(boot + 1, print);
  
  // Where to save ldelta
  Eigen::MatrixXd ldelta(Kv2, boot + 1);
  ldelta.col(0) = fEncompassingStrucDeltaCoef(y, qy1, Z1, qy2, Z2, X, idX1,
             idX2, LnIs, LIs, ngroup, ntau1, ntau2, K2, Kins1, Kins2,
             Kv1, Kv2, iv1, iv2);
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
    std::vector<Eigen::ArrayXi> LnIs_boot(ngroup);
    std::vector<Eigen::ArrayXi> LIs_boot(ngroup);
    std::uniform_int_distribution<int> unidist(0, ngroup - 1);
    for (int s = 0; s < ngroup; ++ s) {
      int sboot    =  unidist(rng);
      LIs_boot[s]  = LIs[sboot];
      LnIs_boot[s] = LnIs[sboot];
    }
    
    ldelta.col(k + 1) = fEncompassingStrucDeltaCoef(y, qy1, Z1, qy2, Z2, X, idX1,
               idX2, LnIs_boot, LIs_boot, ngroup, ntau1, ntau2, K2, Kins1, Kins2,
               Kv1, Kv2, iv1, iv2);
#pragma omp critical
    Prog.increment();
  }
}
#else
std::mt19937 rng(seed); 
for (int k = 0; k < boot; ++ k) {
  // Select subnets
  std::vector<Eigen::ArrayXi> LnIs_boot(ngroup);
  std::vector<Eigen::ArrayXi> LIs_boot(ngroup);
  std::uniform_int_distribution<int> unidist(0, ngroup - 1);
  for (int s = 0; s < ngroup; ++ s) {
    int sboot    =  unidist(rng);
    LIs_boot[s]  = LIs[sboot];
    LnIs_boot[s] = LnIs[sboot];
  }
  
  ldelta.col(k + 1) = fEncompassingStrucDeltaCoef(y, qy1, Z1, qy2, Z2, X, idX1,
             idX2, LnIs_boot, LIs_boot, ngroup, ntau1, ntau2, K2, Kins1, Kins2,
             Kv1, Kv2, iv1, iv2);
  Prog.increment();
}
#endif

Eigen::ArrayXd mdelta(ldelta(Eigen::all, Eigen::seqN(1, boot)).array().rowwise().mean());
Eigen::MatrixXd ddelta(ldelta(Eigen::all, Eigen::seqN(1, boot)).array().colwise() - mdelta);
Eigen::MatrixXd Vdelta(ddelta * ddelta.transpose() / (boot - 1));
return Rcpp::List::create(_["delta"]  = ldelta.col(0), 
                          _["mdelta"] = mdelta,
                          _["Vdelta"] = Vdelta);
}

// Encompassing test KP method
//[[Rcpp::export]]
Rcpp::List fEncompassingStruc(const Eigen::VectorXd& y,
                              const Eigen::MatrixXd& qy1,
                              const Eigen::MatrixXd& Z1,
                              const int& Kest11,
                              const int& Kest21,
                              const Eigen::MatrixXd& qy2,
                              const Eigen::MatrixXd& Z2,
                              const int& Kest12,
                              const int& Kest22,
                              const Eigen::MatrixXd& X, //common to both models
                              const Eigen::ArrayXi& idX1, //common to both models
                              const Eigen::ArrayXi& idX2, //common to both models
                              const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                              const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                              const int& ngroup,
                              const bool& iv1,
                              const bool& iv2, //common to both models
                              const int& boot,
                              const int& nthreads,
                              const unsigned long long seed,
                              const bool& print,
                              const bool& full) {
  int ntau2(qy2.cols()), Kv2(1 + ntau2 + idX2.size());
  // Delta and its variable
  Eigen::VectorXd delta;
  Eigen::VectorXd mdelta;
  Eigen::MatrixXd Vardelta;
  {
    Rcpp::List tp = fEncompassingStrucDelta(y, qy1, Z1, qy2, Z2, X, idX1, idX2, 
                                            LnIs, LIs, ngroup, iv1, iv2, boot, 
                                            nthreads, seed, print);
    delta     = tp["delta"];
    mdelta    = tp["mdelta"];
    Vardelta  = tp["Vdelta"];
  }
  
  // R matrix
  int Kdeltasel(ntau2);
  if (full) {
    Kdeltasel = Kv2;
  }
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(Kdeltasel, Kv2));
  if (full) {
    R(Eigen::all, Eigen::seqN(0, Kdeltasel)) = Eigen::MatrixXd::Identity(Kdeltasel, Kdeltasel);
  } else {
    R(Eigen::all, Eigen::seqN(1, Kdeltasel)) = Eigen::MatrixXd::Identity(Kdeltasel, Kdeltasel);
  }
  
  // delta select
  Eigen::VectorXd deltasel(R * delta);
  Eigen::MatrixXd Vardeltasel(R * Vardelta * R.transpose());
  
  int df1(Kdeltasel), df2(y.size() - Kest12 - Kest22);
  Eigen::ArrayXi itheta(df1);
  if (full) {
    R(Eigen::all, Eigen::seqN(0, df1)) = Eigen::MatrixXd::Identity(df1, df1);
    itheta << Eigen::ArrayXi::LinSpaced(1 + ntau2, 1, 1 + ntau2),  2 + ntau2 + idX2;
  } else {
    R(Eigen::all, Eigen::seqN(1, df1)) = Eigen::MatrixXd::Identity(df1, df1);
    itheta << Eigen::ArrayXi::LinSpaced(ntau2, 2, 1 + ntau2);
  }
  
  // statistic
  // KP
  double statKP(deltasel.dot(ginv(Vardeltasel)*deltasel));
  // K
  double statF(statKP / Kdeltasel);
  
  return Rcpp::List::create(_["KP.stat"] = statKP,
                            _["KP.df"]   = df1,
                            _["F.stat"]  = statF,
                            _["F.df1"]   = df1, 
                            _["F.df2"]   = df2,  
                            _["delta"]   = delta,
                            _["Vdelta"]  = Vardelta,
                            _["mdelta"]  = mdelta,
                            _["itheta"]  = itheta);
}



// Same function for the reduxed form model
Eigen::VectorXd fEncompassingRedDeltaCoef(const Eigen::VectorXd& y,
                                          const Eigen::MatrixXd& qy1,
                                          const Eigen::MatrixXd& Z1,
                                          const Eigen::MatrixXd& qy2,
                                          const Eigen::MatrixXd& Z2,
                                          const Eigen::MatrixXd& X, //common to both models
                                          const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                                          const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                                          const int& ngroup,
                                          const int& ntau1,
                                          const int& ntau2,
                                          const int& K,
                                          const int& Kv1,
                                          const int& Kv2,
                                          const int& Kins1,
                                          const int& Kins2,
                                          const bool& iv1,
                                          const bool& iv2){
  // Find iso and nisolated
  int n = 0;
  for (int s = 0; s < ngroup; ++s){
    n += LIs[s].size() + LnIs[s].size();
  }
  
  Eigen::ArrayXi sel(n);
  int pos = 0;
  for (int s = 0; s < ngroup; ++s) {
    int n1 = LIs[s].size();
    if (n1 > 0) {
      sel.segment(pos, n1) = LIs[s];
      pos += n1;
    }
    
    int n2 = LnIs[s].size();
    if (n2 > 0) {
      sel.segment(pos, n2) = LnIs[s];
      pos += n2;
    }
  }
  
  // Variables and Instruments
  Eigen::MatrixXd V1(n, Kv1), V2(n, Kv2);
  V1 << qy1(sel, Eigen::all), X(sel, Eigen::all);
  V2 << qy2(sel, Eigen::all), X(sel, Eigen::all);
  Eigen::MatrixXd Z1sel = Z1(sel, Eigen::all);
  Eigen::MatrixXd Z2sel = Z2(sel, Eigen::all);
  Eigen::VectorXd ysel = y(sel);
  
  // Weights
  Eigen::MatrixXd W1 = Eigen::MatrixXd::Identity(Kins1, Kins1);
  Eigen::MatrixXd W2 = Eigen::MatrixXd::Identity(Kins2, Kins2);
  if (iv1) {
    W1 = (Z1sel.transpose()*Z1sel).inverse();
  }
  if (iv2) {
    W2 = (Z2sel.transpose()*Z2sel).inverse();
  }
  
  Eigen::MatrixXd ZV1(Z1sel.transpose() * V1), 
  VZW1(ZV1.transpose() * W1), VZWZV1(VZW1 * ZV1);
  Eigen::MatrixXd ZV2(Z2sel.transpose() * V2), 
  VZW2(ZV2.transpose() * W2), VZWZV2(VZW2 * ZV2);
  Eigen::VectorXd Z1y(Z1sel.transpose() * ysel);
  
  Eigen::VectorXd lambda1(VZWZV1.colPivHouseholderQr().solve(VZW1 * Z1y));
  // cout << lambda1.transpose() << endl;
  Eigen::VectorXd e1(ysel - V1 * lambda1);
  
  Eigen::VectorXd Z2e1(Z2sel.transpose() * e1);
  
  return VZWZV2.colPivHouseholderQr().solve(VZW2 * Z2e1);
}

Rcpp::List fEncompassingRedDelta(const Eigen::VectorXd& y,
                                 const Eigen::MatrixXd& qy1,
                                 const Eigen::MatrixXd& Z1,
                                 const Eigen::MatrixXd& qy2,
                                 const Eigen::MatrixXd& Z2,
                                 const Eigen::MatrixXd& X, //common to both models
                                 const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                                 const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                                 const int& ngroup,
                                 const bool& iv1,
                                 const bool& iv2, //common to both models
                                 const int& boot,
                                 const int& nthreads,
                                 const unsigned long long seed,
                                 const bool& print) { 
  // Fixed variables
  int ntau1(qy1.cols()), ntau2(qy2.cols()), K(X.cols()), 
  Kv1(ntau1 + K), Kv2(ntau2 + K),
  Kins1(Z1.cols()), Kins2(Z2.cols());
  
  Progress Prog(boot + 1, print);
  
  // Where to save ldelta
  Eigen::MatrixXd ldelta(Kv2, boot + 1);
  ldelta.col(0) = fEncompassingRedDeltaCoef(y, qy1, Z1, qy2, Z2, X, LnIs, LIs, ngroup,
             ntau1, ntau2, K, Kv1, Kv2, Kins1, Kins2, iv1, iv2);
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
    std::vector<Eigen::ArrayXi> LnIs_boot(ngroup);
    std::vector<Eigen::ArrayXi> LIs_boot(ngroup);
    std::uniform_int_distribution<int> unidist(0, ngroup - 1);
    for (int s = 0; s < ngroup; ++ s) {
      int sboot    =  unidist(rng);
      LIs_boot[s]  = LIs[sboot];
      LnIs_boot[s] = LnIs[sboot];
    }
    
    ldelta.col(k + 1) = fEncompassingRedDeltaCoef(y, qy1, Z1, qy2, Z2, X, LnIs_boot, LIs_boot, ngroup,
               ntau1, ntau2, K, Kv1, Kv2, Kins1, Kins2, iv1, iv2);
#pragma omp critical
    Prog.increment();
  }
}
#else
std::mt19937 rng(seed); 
for (int k = 0; k < boot; ++ k) {
  // Select subnets
  std::vector<Eigen::ArrayXi> LnIs_boot(ngroup);
  std::vector<Eigen::ArrayXi> LIs_boot(ngroup);
  std::uniform_int_distribution<int> unidist(0, ngroup - 1);
  for (int s = 0; s < ngroup; ++ s) {
    int sboot    =  unidist(rng);
    LIs_boot[s]  = LIs[sboot];
    LnIs_boot[s] = LnIs[sboot];
  }
  
  ldelta.col(k + 1) = fEncompassingRedDeltaCoef(y, qy1, Z1, qy2, Z2, X, LnIs_boot, LIs_boot, ngroup,
             ntau1, ntau2, K, Kv1, Kv2, Kins1, Kins2, iv1, iv2);
#pragma omp critical
}
#endif

Eigen::ArrayXd mdelta(ldelta(Eigen::all, Eigen::seqN(1, boot)).array().rowwise().mean());
Eigen::MatrixXd ddelta(ldelta(Eigen::all, Eigen::seqN(1, boot)).array().colwise() - mdelta);
Eigen::MatrixXd Vdelta(ddelta * ddelta.transpose() / (boot - 1));
return Rcpp::List::create(_["delta"]  = ldelta.col(0), 
                          _["mdelta"] = mdelta,
                          _["Vdelta"] = Vdelta);
}

//[[Rcpp::export]]
Rcpp::List fEncompassingRed(const Eigen::VectorXd& y,
                            const Eigen::MatrixXd& qy1,
                            const Eigen::MatrixXd& Z1,
                            const int& Kest1,
                            const Eigen::MatrixXd& qy2,
                            const Eigen::MatrixXd& Z2,
                            const int& Kest2,
                            const Eigen::MatrixXd& X, //common to both models
                            const std::vector<Eigen::ArrayXi>& LnIs, //common to both models
                            const std::vector<Eigen::ArrayXi>& LIs,  //common to both models
                            const int& ngroup,
                            const bool& iv1,
                            const bool& iv2, //common to both models
                            const int& boot,
                            const int& nthreads,
                            const unsigned long long seed,
                            const bool& print,
                            const bool& full) {
  int ntau2(qy2.cols()), Kv2(ntau2 + X.cols());
  // Delta and its variable
  Eigen::VectorXd delta;
  Eigen::VectorXd mdelta;
  Eigen::MatrixXd Vardelta;
  {
    Rcpp::List tp = fEncompassingRedDelta(y, qy1, Z1, qy2, Z2, X, 
                                          LnIs, LIs, ngroup, iv1, iv2, 
                                          boot, nthreads, seed, print);
    delta     = tp["delta"];
    mdelta    = tp["mdelta"];
    Vardelta  = tp["Vdelta"];
  }
  
  // R matrix
  int Kdeltasel(ntau2);
  if (full) {
    Kdeltasel = Kv2;
  }
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(Kdeltasel, Kv2));
  R(Eigen::all, Eigen::seqN(0, Kdeltasel)) = Eigen::MatrixXd::Identity(Kdeltasel, Kdeltasel);
  
  // delta select
  Eigen::VectorXd deltasel(R * delta);
  // variance
  Eigen::MatrixXd Vardeltasel(R * Vardelta * R.transpose());
  
  
  int df1(Kdeltasel), df2(y.rows() - Kest2);
  Eigen::ArrayXi itheta(df1);
  itheta << Eigen::ArrayXi::LinSpaced(df1, 1, df1);
  // statistic
  // KP
  double statKP(deltasel.dot(ginv(Vardeltasel)*deltasel));
  // F stat
  double statF(statKP / Kdeltasel);
  return Rcpp::List::create(_["KP.stat"] = statKP,
                            _["KP.df"]   = df1,
                            _["F.stat"]  = statF,
                            _["F.df1"]   = df1, 
                            _["F.df2"]   = df2,  
                            _["delta"]   = delta,
                            _["Vdelta"]  = Vardelta,
                            _["mdelta"]  = mdelta,
                            _["itheta"]  = itheta);
}




// // Encompassing test with projection onto Z
// //[[Rcpp::export]]
// Rcpp::List fEncompassingStruc2(const Eigen::MatrixXd& X1,
//                                const Eigen::MatrixXd& qy1,
//                                const Eigen::MatrixXd& Z1,
//                                const Eigen::MatrixXd& W11,
//                                const Eigen::MatrixXd& W21,
//                                const Eigen::VectorXd& e1,
//                                const Eigen::VectorXd theta1,
//                                const Eigen::ArrayXi& idX11,
//                                const Eigen::ArrayXi& idX21,
//                                const int& Kest11,
//                                const int& Kest21,
//                                const Eigen::MatrixXd& X2,
//                                const Eigen::MatrixXd& qy2,
//                                const Eigen::MatrixXd& Z2,
//                                const Eigen::MatrixXd& W12,
//                                const Eigen::MatrixXd& W22,
//                                const Eigen::VectorXd& e2,
//                                const Eigen::VectorXd theta2,
//                                const Eigen::ArrayXi& idX12,
//                                const Eigen::ArrayXi& idX22,
//                                const int& Kest12,
//                                const int& Kest22,
//                                const Eigen::ArrayXi& nIs, //common to both models
//                                const Eigen::ArrayXi& Is, //common to both models
//                                const int& ngroup, //common to both models
//                                const Eigen::ArrayXi& cumsn, //common to both models
//                                const int& HAC = 0, //common to both models
//                                const bool& full = false) {
//   int Kins1(Z1.cols()), Kins2(Z2.cols()), n(X1.rows()), K21(idX21.size()), 
//   K22(idX22.size()), ntau1(qy1.cols()), ntau2(qy2.cols()), n_iso(Is.size()), n_niso(n - n_iso);
//   
//   // First stage
//   Eigen::MatrixXd X11(X1(Is, idX11)), XX11(X11.transpose()*X11), XXW11(XX11*W11), XXWXX11(XXW11*XX11);
//   Eigen::MatrixXd X12(X2(Is, idX12)), XX12(X12.transpose()*X12), XXW12(XX12*W12), XXWXX12(XXW12*XX12);
//   Eigen::VectorXd b1(theta1(1 + ntau1 + idX11)), b2(theta2(1 + ntau2 + idX12));
//   
//   // Second stage
//   Eigen::VectorXd Xb1(X1(Eigen::all, idX11)*b1), Xb11(Xb1(Is)), Xb21(Xb1(nIs));
//   Eigen::VectorXd Xb2(X2(Eigen::all, idX12)*b2), Xb12(Xb2(Is)), Xb22(Xb2(nIs));
//   
//   Eigen::MatrixXd X21(X1(nIs, idX21)),  X211(X1(nIs, idX11)), V21(n_niso, 1 + ntau1 + K21);
//   V21 << Xb21, qy1(nIs, Eigen::all), X21;
//   Eigen::MatrixXd X22(X2(nIs, idX22)),  X212(X2(nIs, idX12)), V22(n_niso, 1 + ntau2 + K22);
//   V22 << Xb22, qy2(nIs, Eigen::all), X22;
//   
//   Eigen::MatrixXd Z21(n_niso, 1 + Kins1);
//   Z21 << Xb21, Z1(nIs, Eigen::all);
//   Eigen::MatrixXd Z22(n_niso, 1 + Kins2);
//   Z22 << Xb22, Z2(nIs, Eigen::all);
//   
//   //residuals
//   Eigen::VectorXd e21(e1(nIs)*(1 - theta1(0)));
//   // Eigen::VectorXd e22(e2(nIs)*(1 - theta2(0)));
//   
//   //H
//   Eigen::MatrixXd Z2V2(Z22.transpose() * V22);
//   Eigen::MatrixXd V2Z2W2(Z2V2.transpose()*W22);
//   Eigen::MatrixXd Z2Z1iZZ1(((Z21.transpose() * Z21).colPivHouseholderQr().solve(Z21.transpose() * Z22)).transpose());
//   Eigen::MatrixXd H((V2Z2W2*Z2V2).colPivHouseholderQr().solve(V2Z2W2 * Z2Z1iZZ1));
//   // cout<<((V2Z2W2*Z2V2).colPivHouseholderQr().solve(V2Z2W2*(Z22.transpose()*e21))).transpose()<<endl;
//   
//   // Variance Ze1
//   Eigen::MatrixXd Z1e1((Z21.array().colwise()*e21.array()).matrix());
//   Eigen::MatrixXd VZ1e1(Eigen::MatrixXd::Zero(Kins2 + 1, Kins2 + 1));
//   if (HAC <= 1) {
//     VZ1e1 = Z1e1.transpose()*Z1e1;
//   }
//   if (HAC == 2) {
//     Eigen::MatrixXd Z(Eigen::MatrixXd::Zero(n, 1 + Kins1));
//     Z(nIs, Eigen::all) << Xb21, Z1(nIs, Eigen::all);
//     Eigen::VectorXd e(Eigen::VectorXd::Zero(n));
//     e(nIs) = e21;
//     for (int r(0); r < ngroup; ++ r) {
//       int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
//       Eigen::VectorXd tp(Z(Eigen::seq(n1, n2), Eigen::all).transpose()*e.segment(n1, n2));
//       VZ1e1 += tp*tp.transpose();
//     }
//   }
//   
//   // R matrix
//   int df(ntau2);
//   if (full) { 
//     df = 1 + ntau2 + K22; 
//   }
//   Eigen::MatrixXd R(Eigen::MatrixXd::Zero(df, 1 + ntau2 + K22));
//   if (full) {
//     R(Eigen::all, Eigen::seqN(0, df)) = Eigen::MatrixXd::Identity(df, df);
//   } else {
//     R(Eigen::all, Eigen::seqN(1, df)) = Eigen::MatrixXd::Identity(df, df);
//   }
//   
//   
//   // statistic and its variance
//   Eigen::VectorXd HZ1e1(H*Z21.transpose()*e21);
//   Eigen::VectorXd RHZ1e1(R*HZ1e1);
//   Eigen::MatrixXd RH(R*H);
//   Eigen::MatrixXd Vstat(RH*VZ1e1*RH.transpose());
//   Eigen::MatrixXd stat(RHZ1e1.transpose()*ginv(Vstat)*RHZ1e1);
//   
//   return Rcpp::List::create(_["stat"] = stat, _["df"] = df, _["diff"] = HZ1e1);
// }
// 
// 
// //[[Rcpp::export]]
// Rcpp::List fEncompassingRed2(const Eigen::MatrixXd& X1,
//                              const Eigen::MatrixXd& qy1,
//                              const Eigen::MatrixXd& Z1,
//                              const Eigen::MatrixXd& W1,
//                              const Eigen::VectorXd& e1,
//                              const Eigen::VectorXd theta1,
//                              const int& Kest1,
//                              const Eigen::MatrixXd& X2,
//                              const Eigen::MatrixXd& qy2,
//                              const Eigen::MatrixXd& Z2,
//                              const Eigen::MatrixXd& W2,
//                              const Eigen::VectorXd& e2,
//                              const Eigen::VectorXd theta2,
//                              const int& Kest2,
//                              const int& ngroup, //common to both models
//                              const Eigen::ArrayXi& cumsn, //common to both models
//                              const int& HAC = 0, //common to both models
//                              const bool& full = false) {
//   int Kins2(Z2.cols()), n(X1.rows()), K1(X1.cols()), K2(X2.cols()), 
//   ntau1(qy1.cols()), ntau2(qy2.cols());
//   
//   Eigen::MatrixXd V1(n, ntau1 + K1), V2(n, ntau2 + K2);
//   V1 << qy1, X1;
//   V2 << qy2, X2;
//   
//   //residuals
//   // Eigen::VectorXd e21(e1(nIs)*(1 - theta1(0)));
//   // Eigen::VectorXd e22(e2(nIs)*(1 - theta2(0)));
//   
//   //H
//   Eigen::MatrixXd Z2V2(Z2.transpose()*V2);
//   Eigen::MatrixXd V2Z2W2(Z2V2.transpose()*W2);
//   Eigen::MatrixXd Z2Z1iZZ1(((Z1.transpose() * Z1).colPivHouseholderQr().solve(Z1.transpose() * Z2)).transpose());
//   Eigen::MatrixXd H((V2Z2W2*Z2V2).colPivHouseholderQr().solve(V2Z2W2 * Z2Z1iZZ1));
//   
//   // Variance Ze1
//   Eigen::MatrixXd Z1e1((Z1.array().colwise()*e1.array()).matrix());
//   Eigen::MatrixXd VZ1e1(Eigen::MatrixXd::Zero(Kins2, Kins2));
//   if (HAC <= 1) {
//     VZ1e1 = Z1e1.transpose()*Z1e1;
//   }
//   if (HAC == 2) {
//     for (int r(0); r < ngroup; ++ r) {
//       int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
//       Eigen::VectorXd tp(Z1(Eigen::seq(n1, n2), Eigen::all).transpose()*e1.segment(n1, n2));
//       VZ1e1 += tp*tp.transpose();
//     }
//   }
//   
//   // R matrix
//   int df(ntau2);
//   if (full) { 
//     df = ntau2 + K2; 
//   }
//   Eigen::MatrixXd R(Eigen::MatrixXd::Zero(df, ntau2 + K2));
//   R(Eigen::all, Eigen::seqN(0, df)) = Eigen::MatrixXd::Identity(df, df);
//   
//   // statistic and its variance
//   Eigen::VectorXd HZ1e1(H*Z1.transpose()*e1);
//   Eigen::VectorXd RHZ1e1(R*HZ1e1);
//   Eigen::MatrixXd RH(R*H);
//   Eigen::MatrixXd Vstat(RH*VZ1e1*RH.transpose());
//   Eigen::MatrixXd stat(RHZ1e1.transpose()*ginv(Vstat)*RHZ1e1);
//   
//   return Rcpp::List::create(_["stat"] = stat, _["df"] = df, _["diff"] = HZ1e1);
// }
// 
// 
