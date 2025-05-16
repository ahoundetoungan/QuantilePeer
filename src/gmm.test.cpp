// [[Rcpp::depends(RcppArmadillo, RcppNumerical, RcppEigen)]]
#include <RcppArmadillo.h>
//#define NDEBUG
#include <RcppNumerical.h>
#include <RcppEigen.h>

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
Rcpp::List Cov2ThetaRed(const Eigen::MatrixXd& X1,
                        const Eigen::MatrixXd& qy1,
                        const Eigen::MatrixXd& Z1,
                        const Eigen::MatrixXd& W1,
                        const Eigen::ArrayXd& e1,
                        const Eigen::VectorXd theta1,
                        const int& Kest1,
                        const Eigen::MatrixXd& X2,
                        const Eigen::MatrixXd& qy2,
                        const Eigen::MatrixXd& Z2,
                        const Eigen::MatrixXd& W2,
                        const Eigen::ArrayXd& e2,
                        const Eigen::VectorXd theta2,
                        const int& Kest2,
                        const int& ngroup, //common to both models
                        const Eigen::ArrayXi& cumsn, //common to both models
                        const int& HAC = 0, //common to both models
                        const bool& full = false) { 
  int ntau1(qy1.cols()), ntau2(qy2.cols()), Kins1(Z1.cols()), Kins2(Z2.cols()), 
  Kins(Kins1 + Kins2), n(X1.rows()), K1(X1.cols()), K2(X2.cols()), 
  Kv1(K1 + ntau1), Kv2(K2 + ntau2), Kv(Kv1 + Kv2);
  
  Eigen::MatrixXd V1(n, Kv1), V2(n, Kv2);
  V1 << qy1, X1;
  V2 << qy2, X2;
  
  Eigen::MatrixXd VZ1(V1.transpose()*Z1), VZW1(VZ1*W1), VZWZV1(VZW1*VZ1.transpose());
  Eigen::MatrixXd VZ2(V2.transpose()*Z2), VZW2(VZ2*W2), VZWZV2(VZW2*VZ2.transpose());
  Eigen::MatrixXd VZWZV(Eigen::MatrixXd::Zero(Kv, Kv)), VZW(Eigen::MatrixXd::Zero(Kv, Kins));
  VZWZV.block(0, 0, Kv1, Kv1)       = VZWZV1;
  VZWZV.block(Kv1, Kv1, Kv2, Kv2)   = VZWZV2;
  VZW.block(0, 0, Kv1, Kins1)       = VZW1;
  VZW.block(Kv1, Kins1, Kv2, Kins2) = VZW2;
  
  // variance of Ze
  Eigen::MatrixXd VZe(Eigen::MatrixXd::Zero(Kins, Kins));
  if (HAC == 0) {
    double s2((e1.square().sum() + e2.square().sum())/(2*n - Kest1 - Kest2));
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
  Eigen::VectorXd theta(Kv);
  theta << theta1, theta2;
  
  // test all parameters of just quantiles
  int df(ntau1);
  if (full) { 
    df = Kv1; 
  }
  
  // R matrix
  Eigen::ArrayXi itheta(df);
  itheta << Eigen::ArrayXi::LinSpaced(df, 1, df);
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(df, Kv));
  R(Eigen::all, Eigen::seqN(0, df))   = Eigen::MatrixXd::Identity(df, df);
  R(Eigen::all, Eigen::seqN(Kv1, df)) = -Eigen::MatrixXd::Identity(df, df);
  
  // test statistic
  Eigen::VectorXd Rtheta(R*theta);
  Eigen::MatrixXd RiHdF(R*iHdF);
  Eigen::MatrixXd VarRtheta(RiHdF * HVFH * RiHdF.transpose());
  Eigen::VectorXd stat(Rtheta.transpose() * ginv(VarRtheta) * Rtheta);
  
  return Rcpp::List::create(_["stat"]    = stat, 
                            _["df"]      = df, 
                            _["dtheta"]  = Rtheta,
                            _["Vdtheta"] = VarRtheta,
                            _["itheta"]  = itheta);
}


// Covariance matrix of two theta using two different sets of instruments, Structural form
//[[Rcpp::export]]
Rcpp::List Cov2ThetaStruc(const Eigen::MatrixXd& X1,
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
                          const int& HAC = 0, //common to both models
                          const bool& full = false) {
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
    double s21((e11.square().sum() + e12.square().sum())/(2*n_iso - Kest11 - Kest12));
    double s22((e21.square().sum() + e22.square().sum())/(2*n_niso - Kest21 - Kest22));
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
  
  // theta 
  Eigen::VectorXd theta(2*K + 2*ntau + 2);
  theta << theta1, theta2;
  
  // test all parameters of just quantiles
  int df(ntau);
  if (full) { 
    df = 1 + ntau + K2; 
  }
  
  // R matrix for the test
  Eigen::ArrayXi itheta(df);
  Eigen::MatrixXd Rt(Eigen::MatrixXd::Zero(df, 2*K + 2*ntau + 2));
  if (full) { 
    Rt(Eigen::seqN(0, 1 + ntau), Eigen::seqN(0, 1 + ntau)) = Eigen::MatrixXd::Identity(1 + ntau, 1 + ntau);
    Rt(Eigen::seqN(0, 1 + ntau), Eigen::seqN(1 + ntau + K, 1 + ntau)) = -Eigen::MatrixXd::Identity(1 + ntau, 1 + ntau);
    if (K2 > 0) {
      Rt(Eigen::seqN(1 + ntau, K2), 1 + ntau + idX2) = Eigen::MatrixXd::Identity(K2, K2);
      Rt(Eigen::seqN(1 + ntau, K2), 2 + 2*ntau + K + idX2) = -Eigen::MatrixXd::Identity(K2, K2); 
    }
    itheta << Eigen::ArrayXi::LinSpaced(1 + ntau, 1, 1 + ntau),  1 + ntau + idX2;
  } else {
    Rt(Eigen::all, Eigen::seqN(1, ntau)) = Eigen::MatrixXd::Identity(ntau, ntau);
    Rt(Eigen::all, Eigen::seqN(2 + ntau + K, ntau)) = -Eigen::MatrixXd::Identity(ntau, ntau);
    itheta << Eigen::ArrayXi::LinSpaced(ntau, 2, 1 + ntau);
  }
  
  // test statistic
  Eigen::VectorXd Rtheta(Rt * theta);
  Eigen::MatrixXd RtR(Rt * R);
  Eigen::MatrixXd VarRtheta(RtR * Vpa * RtR.transpose());
  Eigen::VectorXd stat(Rtheta.transpose() * ginv(VarRtheta) * Rtheta);
  
  return Rcpp::List::create(_["stat"]    = stat, 
                            _["df"]      = df, 
                            _["dtheta"]  = Rtheta,
                            _["Vdtheta"] = VarRtheta,
                            _["itheta"]  = itheta);
}

// This function implement test for validity of Z2 using the fact that Z1 is valid
//[[Rcpp::export]]
Rcpp::List validZ2SarganRed(const Eigen::MatrixXd& X1,
                            const Eigen::MatrixXd& qy1,
                            const Eigen::MatrixXd& Z1,
                            const Eigen::MatrixXd& W1,
                            const Eigen::ArrayXd& e1,
                            const Eigen::VectorXd theta1,
                            const int& Kest1,
                            const Eigen::MatrixXd& X2,
                            const Eigen::MatrixXd& qy2,
                            const Eigen::MatrixXd& Z2,
                            const Eigen::MatrixXd& W2,
                            const Eigen::ArrayXd& e2,
                            const Eigen::VectorXd theta2,
                            const int& Kest2,
                            const int& ngroup, //common to both models
                            const Eigen::ArrayXi& cumsn, //common to both models
                            const int& HAC = 0, //common to both models
                            const bool& full = false) {
  int ntau1(qy1.cols()), ntau2(qy2.cols()), Kins2(Z2.cols()), n(X1.rows());
  Eigen::MatrixXd H(Z2 - Z1*(Z1.transpose() * Z1).colPivHouseholderQr().solve(Z1.transpose() * Z2));
  Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
  int df(lu.rank());
  
  // variance of He
  Eigen::MatrixXd VHe(Eigen::MatrixXd::Zero(Kins2, Kins2));
  if (HAC == 0) {
    double s2(e1.square().sum()/(n - Kest1));
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
                            _["theta1"] = theta1.head(ntau1),
                            _["theta2"] = theta2.head(ntau2));
}


//[[Rcpp::export]]
Rcpp::List validZ2SarganStruc(const Eigen::MatrixXd& X1,
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
                              const int& HAC = 0, //common to both models
                              const bool& full = false) {
  int Kins1(Z1.cols()), Kins2(Z2.cols()), n(X1.rows()), ntau(qy1.cols()), n_iso(Is.size()), n_niso(n - n_iso);
  
  // First stage
  Eigen::VectorXd b1(theta1(1 + ntau + idX1)), b2(theta2(1 + ntau + idX1));
  
  // Instruments
  Eigen::VectorXd Xb1(X1(Eigen::all, idX1)*b1), Xb21(Xb1(nIs));
  Eigen::VectorXd Xb2(X2(Eigen::all, idX1)*b2), Xb22(Xb2(nIs));
  Eigen::MatrixXd Z21(n_niso, 1 + Kins1);
  Z21 << Xb21, Z1(nIs, Eigen::all);
  Eigen::MatrixXd Z22(n_niso, 1 + Kins2);
  Z22 << Xb22, Z2(nIs, Eigen::all);
  
  Eigen::MatrixXd H(Z22 - Z21*(Z21.transpose() * Z21).colPivHouseholderQr().solve(Z21.transpose() * Z22));
  Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
  int df(lu.rank());
  
  // variance of He
  Eigen::MatrixXd VHe(Eigen::MatrixXd::Zero(Kins2 + 1, Kins2 + 1));
  if (HAC == 0) {
    double s2(e1(nIs).square().sum()/(n_niso - Kest21));
    VHe = s2 * (H.transpose() * H);
  }
  if (HAC == 1) {
    Eigen::MatrixXd He((H.array().colwise() * e1(nIs)).matrix());
    VHe    = He.transpose()*He;
  }
  if (HAC == 2) {
    Eigen::ArrayXXd He(Eigen::MatrixXd::Zero(n, Kins2 + 1));
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


// Encompassing test
// Estimating delta and its variance
Rcpp::List fEncompassingStrucDelta(const Eigen::MatrixXd& X1,
                                   const Eigen::MatrixXd& qy1,
                                   const Eigen::MatrixXd& Z1,
                                   const Eigen::MatrixXd& W11,
                                   const Eigen::MatrixXd& W21,
                                   const Eigen::VectorXd& e1,
                                   const Eigen::VectorXd theta1,
                                   const Eigen::ArrayXi& idX11,
                                   const Eigen::ArrayXi& idX21,
                                   const int& Kest11,
                                   const int& Kest21,
                                   const Eigen::MatrixXd& X2,
                                   const Eigen::MatrixXd& qy2,
                                   const Eigen::MatrixXd& Z2,
                                   const Eigen::MatrixXd& W12,
                                   const Eigen::MatrixXd& W22,
                                   const Eigen::VectorXd& e2,
                                   const Eigen::VectorXd theta2,
                                   const Eigen::ArrayXi& idX12,
                                   const Eigen::ArrayXi& idX22,
                                   const int& Kest12,
                                   const int& Kest22,
                                   const Eigen::ArrayXi& nIs, //common to both models
                                   const Eigen::ArrayXi& Is, //common to both models
                                   const int& ngroup, //common to both models
                                   const Eigen::ArrayXi& cumsn, //common to both models
                                   const int& HAC = 0) { //common to both models
  int n(X1.rows()), n_iso(Is.size()), n_niso(n - n_iso), 
  ntau1(qy1.cols()), ntau2(qy2.cols()),
  K11(idX11.size()), K12(idX12.size()), K21(idX21.size()), K22(idX22.size()), 
  K1(K11 + K12), 
  Kins1(Z1.cols()), Kins2(Z2.cols()), Kins(Kins1 + Kins2),
  Kv1(1 + ntau1 + K21), Kv2(1 + ntau2 + K22), Kv(Kv1 + Kv2);
  
  // First stage
  Eigen::MatrixXd X11(X1(Is, idX11)), XX11(X11.transpose()*X11), XXW11(XX11*W11), XXWXX11(XXW11*XX11);
  Eigen::MatrixXd X12(X2(Is, idX12)), XX12(X12.transpose()*X12), XXW12(XX12*W12), XXWXX12(XXW12*XX12);
  Eigen::VectorXd b1(theta1(1 + ntau1 + idX11)), b2(theta2(1 + ntau2 + idX12));
  
  // Variables and Instruments
  Eigen::VectorXd Xb1(X1(Eigen::all, idX11)*b1), Xb11(Xb1(Is)), Xb21(Xb1(nIs));
  Eigen::VectorXd Xb2(X2(Eigen::all, idX12)*b2), Xb12(Xb2(Is)), Xb22(Xb2(nIs));
  
  Eigen::MatrixXd X21(X1(nIs, idX21)), X221(X1(nIs, idX21)),  X211(X1(nIs, idX11)), V21(n_niso, Kv1);
  V21 << Xb21, qy1(nIs, Eigen::all), X21;
  Eigen::MatrixXd X22(X2(nIs, idX22)), X222(X2(nIs, idX22)),  X212(X2(nIs, idX12)), V22(n_niso, Kv2);
  V22 << Xb22, qy2(nIs, Eigen::all), X22;
  
  Eigen::MatrixXd Z21(n_niso, 1 + Kins1);
  Z21 << Xb21, Z1(nIs, Eigen::all);
  Eigen::MatrixXd Z22(n_niso, 1 + Kins2);
  Z22 << Xb22, Z2(nIs, Eigen::all);
  
  Eigen::MatrixXd ZV21(Z21.transpose()*V21);
  Eigen::MatrixXd ZV22(Z22.transpose()*V22);
  
  Eigen::MatrixXd VZW21(ZV21.transpose()*W21), VZWZV21(VZW21*ZV21);
  Eigen::MatrixXd VZW22(ZV22.transpose()*W22), VZWZV22(VZW22*ZV22);
  
  Eigen::ArrayXd e11(e1(Is)), e21(e1(nIs)*(1 - theta1(0)));
  Eigen::ArrayXd e12(e2(Is)), e22(e2(nIs)*(1 - theta2(0)));
  Eigen::VectorXd delta((VZWZV22).colPivHouseholderQr().solve(VZW22 * Z22.transpose() * e21.matrix()));
  
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(K1 + Kv, K1 + Kins + 2));
  H.block(0, 0, K11, K11)  = XXW11;
  H.block(K11, K11, Kv1, Kins1 + 1) = VZW21;
  H.block(K11 + Kv1, K11 + Kins1 + 1, K12, K12)  = XXW12;
  H.block(K1 + Kv1, K1 + Kins1 + 1, Kv2, Kins2 + 1) = VZW22;
  
  Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(K1 + Kins + 2, K1 + Kv));
  dF.block(0, 0, K11, K11) = XX11;
  dF.block(K11, 0, Kins1 + 1, K11 + Kv1) << (Z21.transpose()*X211*(1 - theta1(0))), ZV21;
  dF.block(K11 + Kins1 + 1, K11 + Kv1, K12, K12) = XX12;
  dF(Eigen::seqN(K1 + Kins1 + 1, Kins2 + 1) , Eigen::all) << Z22.transpose()*X211*(1 - theta1(0)),
                                                             Z22.transpose()*V21, 
                                                             Z22.transpose()*X212*(1 - theta2(0)), 
                                                             ZV22;
  
  Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(K1 + Kins + 2, K1 + Kins + 2));
  if (HAC <= 1) {
    Eigen::MatrixXd Xe11(X11.array().colwise()*e11), Ze21(Z21.array().colwise()*e21);
    Eigen::MatrixXd Xe12(X12.array().colwise()*e12), Z22edelta(Z22.array().colwise()*(e21 - (V22 * delta).array()));
    Eigen::MatrixXd Xe11Xe12(Xe11.transpose()*Xe12), Ze21Z22delta(Ze21.transpose()*Z22edelta);
    VF.block(0, 0, K11, K11) = Xe11.transpose()*Xe11;
    VF.block(0, K11 + Kins1 + 1, K11, K12) = Xe11Xe12;
    VF.block(K11, K11, Kins1 + 1, Kins1 + 1) = Ze21.transpose()*Ze21;
    VF.block(K11, K1 + Kins1 + 1, Kins1 + 1, Kins2 + 1) = Ze21Z22delta;
    VF.block(K11 + Kins1 + 1, 0, K12, K11) = Xe11Xe12.transpose();
    VF.block(K11 + Kins1 + 1, K11 + Kins1 + 1, K12, K12) = Xe12.transpose()*Xe12;
    VF.block(K1 + Kins1 + 1, K11, Kins2 + 1, Kins1 + 1) = Ze21Z22delta.transpose();
    VF.block(K1 + Kins1 + 1, K1 + Kins1 + 1, Kins2 + 1, Kins2 + 1) = Z22edelta.transpose()*Z22edelta;
  }
  if (HAC == 2) {
    X11   = Eigen::MatrixXd::Zero(n, K11);
    X11(Is, Eigen::all) = X1(Is, idX11);
    X12   = Eigen::MatrixXd::Zero(n, K12);
    X12(Is, Eigen::all) = X2(Is, idX12);
    Z21   = Eigen::MatrixXd::Zero(n, 1 + Kins1);
    Z21(nIs, Eigen::all) << Xb21, Z1(nIs, Eigen::all);
    Z22   = Eigen::MatrixXd::Zero(n, 1 + Kins2);
    Z22(nIs, Eigen::all) << Xb22, Z2(nIs, Eigen::all);
    Eigen::VectorXd est1(Eigen::VectorXd::Zero(n)), est2(Eigen::VectorXd::Zero(n));
    est1(Is)  = e11; est1(nIs) = e21;
    est2(Is)  = e12; est2(nIs) = e21.matrix() - V22 * delta;
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::MatrixXd tp11(n2 - n1 + 1, K11 + Kins1 + 1), tp12(n2 - n1 + 1, K12 + Kins2 + 1);
      tp11 << X11(Eigen::seq(n1, n2), Eigen::all), Z21(Eigen::seq(n1, n2), Eigen::all);
      tp12 << X12(Eigen::seq(n1, n2), Eigen::all), Z22(Eigen::seq(n1, n2), Eigen::all);
      Eigen::VectorXd tp2(K1 + Kins + 2);
      tp2 << tp11.transpose()*est1.segment(n1, n2), tp12.transpose()*est2.segment(n1, n2);
      VF += tp2*tp2.transpose();
    }
  }
  
  Eigen::MatrixXd iHdF((H*dF).inverse()); 
  Eigen::MatrixXd HVFH(H*VF*H.transpose());
  Eigen::MatrixXd Vdelta(iHdF * HVFH * iHdF.transpose());
  // cout<<Vdelta<<endl;
  
  return Rcpp::List::create(_["delta"]  = delta, 
                            _["Vdelta"] = Vdelta(Eigen::seqN(K1 + Kv1, Kv2), Eigen::seqN(K1 + Kv1, Kv2)));
}

Rcpp::List fEncompassingRedDelta(const Eigen::MatrixXd& X1,
                                 const Eigen::MatrixXd& qy1,
                                 const Eigen::MatrixXd& Z1,
                                 const Eigen::MatrixXd& W1,
                                 const Eigen::VectorXd& e1,
                                 const Eigen::VectorXd theta1,
                                 const int& Kest1,
                                 const Eigen::MatrixXd& X2,
                                 const Eigen::MatrixXd& qy2,
                                 const Eigen::MatrixXd& Z2,
                                 const Eigen::MatrixXd& W2,
                                 const Eigen::VectorXd& e2,
                                 const Eigen::VectorXd theta2,
                                 const int& Kest2,
                                 const int& ngroup, //common to both models
                                 const Eigen::ArrayXi& cumsn, //common to both models
                                 const int& HAC = 0) { //common to both models
  int n(X1.rows()), ntau1(qy1.cols()), ntau2(qy2.cols()),
  K1(X1.cols()), K2(X2.cols()),  
  Kins1(Z1.cols()), Kins2(Z2.cols()), Kins(Kins1 + Kins2),
  Kv1(K1 + ntau1), Kv2(K2 + ntau2), Kv(Kv1 + Kv2);
  
  // Variables and Instruments
  Eigen::MatrixXd V1(n, Kv1), V2(n, Kv2);
  V1 << qy1, X1;
  V2 << qy2, X2;
  
  Eigen::MatrixXd ZV1(Z1.transpose() * V1), VZW1(ZV1.transpose() * W1), VZWZV1(VZW1 * ZV1);
  Eigen::MatrixXd ZV2(Z2.transpose() * V2), VZW2(ZV2.transpose() * W2), VZWZV2(VZW2 * ZV2);
  Eigen::VectorXd delta(VZWZV2.colPivHouseholderQr().solve(VZW2 * Z2.transpose() * e1));
  
  Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kv, Kins));
  H.block(0, 0, Kv1, Kins1) = VZW1;
  H.block(Kv1, Kins1, Kv2, Kins2) = VZW2;
  
  Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(Kins, Kv));
  dF.block(0, 0, Kins1, Kv1) = ZV1;
  dF(Eigen::seqN(Kins1, Kins2) , Eigen::all) << Z2.transpose() * V1, ZV2;
  
  Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(Kins, Kins));
  if (HAC <= 1) {
    Eigen::MatrixXd Ze(n, Kins);
    Ze << (Z1.array().colwise() * e1.array()).matrix(), (Z2.array().colwise() * (e1 - V2 * delta).array()).matrix();
    VF    = Ze.transpose()*Ze;
    // cout<<e1.array().square().sum()/(n - Kest1)<<endl;
  }
  if (HAC == 2) {
    Eigen::ArrayXXd Ze(n, Kins);
    Ze << (Z1.array().colwise() * e1.array()).matrix(), (Z2.array().colwise() * (e1 - V2 * delta).array()).matrix();
    for (int r(0); r < ngroup; ++ r) {
      int n1(cumsn(r)), n2(cumsn(r + 1) - 1);
      Eigen::VectorXd tp(Ze(Eigen::seq(n1, n2), Eigen::all).array().colwise().sum().matrix());
      VF += tp*tp.transpose();
    }
  }
  
  Eigen::MatrixXd iHdF((H*dF).inverse());
  Eigen::MatrixXd HVFH(H*VF*H.transpose());
  Eigen::MatrixXd Vdelta(iHdF * HVFH * iHdF.transpose());
  
  return Rcpp::List::create(_["delta"]  = delta,
                            _["Vdelta"] = Vdelta(Eigen::seqN(Kv1, Kv2), Eigen::seqN(Kv1, Kv2)));
}

// Encompassing test KP method
//[[Rcpp::export]]
Rcpp::List fEncompassingStrucKP(const Eigen::MatrixXd& X1,
                                const Eigen::MatrixXd& qy1,
                                const Eigen::MatrixXd& Z1,
                                const Eigen::MatrixXd& W11,
                                const Eigen::MatrixXd& W21,
                                const Eigen::VectorXd& e1,
                                const Eigen::VectorXd theta1,
                                const Eigen::ArrayXi& idX11,
                                const Eigen::ArrayXi& idX21,
                                const int& Kest11,
                                const int& Kest21,
                                const Eigen::MatrixXd& X2,
                                const Eigen::MatrixXd& qy2,
                                const Eigen::MatrixXd& Z2,
                                const Eigen::MatrixXd& W12,
                                const Eigen::MatrixXd& W22,
                                const Eigen::VectorXd& e2,
                                const Eigen::VectorXd theta2,
                                const Eigen::ArrayXi& idX12,
                                const Eigen::ArrayXi& idX22,
                                const int& Kest12,
                                const int& Kest22,
                                const Eigen::ArrayXi& nIs, //common to both models
                                const Eigen::ArrayXi& Is, //common to both models
                                const int& ngroup, //common to both models
                                const Eigen::ArrayXi& cumsn, //common to both models
                                const int& HAC = 0, //common to both models
                                const bool& full = false) {
  int ntau2(qy2.cols()), Kv2(1 + ntau2 + idX22.size());
  // Delta and its variable
  Eigen::VectorXd deltaFull;
  Eigen::MatrixXd VardeltaFull;
  {
    Rcpp::List tp = fEncompassingStrucDelta(X1, qy1, Z1, W11, W21, e1, theta1, idX11, idX21, Kest11,
                                            Kest21, X2, qy2, Z2, W12, W22, e2, theta2, idX12, idX22,
                                            Kest12, Kest22, nIs, Is, ngroup, cumsn, HAC);
    deltaFull     = tp["delta"];
    VardeltaFull  = tp["Vdelta"];
    // cout<<VardeltaFull<<endl;
  }
  
  
  // R matrix
  int Kdelta(ntau2);
  if (full) {
    Kdelta = Kv2;
  }
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(Kdelta, Kv2));
  if (full) {
    R(Eigen::all, Eigen::seqN(0, Kdelta)) = Eigen::MatrixXd::Identity(Kdelta, Kdelta);
  } else {
    R(Eigen::all, Eigen::seqN(1, Kdelta)) = Eigen::MatrixXd::Identity(Kdelta, Kdelta);
  }
  
  // delta
  Eigen::VectorXd delta(R * deltaFull);
  
  // variance
  Eigen::MatrixXd Vardelta(R * VardeltaFull * R.transpose());
  // cout<<"KP - delta"<<endl;
  // cout<<delta.transpose()<<endl;
  // cout<<"KP - Vardelta"<<endl;
  // cout<<Vardelta<<endl;
  
  
  // Normalization
  Eigen::LLT<Eigen::MatrixXd> tpOmega(Vardelta);
  Eigen::MatrixXd Omega(tpOmega.matrixL());
  Eigen::VectorXd deltatilde(Omega.colPivHouseholderQr().solve(delta));
  // cout<<"KP - deltatilde"<<endl;
  // cout<<deltatilde.transpose()<<endl;
  // cout<<"KP - Vardeltatilde"<<endl;
  // cout<<Omega.inverse() * Vardelta * Omega.inverse().transpose()<<endl;
  // Test if variance is I
  // Eigen::MatrixXd Vtp(Omega.colPivHouseholderQr().solve(Vardelta) * (Omega.transpose()).inverse());
  // cout << Vtp << endl;
  // cout << delta.transpose() << endl;
  // cout << Vardelta.diagonal().array().sqrt().transpose() << endl;
  // cout << deltatilde.transpose() << endl;
  
  // SDV decomposition of Theta
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(deltatilde.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd U = svd.matrixU(); //ntau * ntau
  Eigen::VectorXd d = svd.singularValues();
  Eigen::MatrixXd ddiag = d.asDiagonal();
  Eigen::MatrixXd D(1, Kdelta);
  D << ddiag, Eigen::MatrixXd::Zero(1, Kdelta - 1); 
  Eigen::MatrixXd V = svd.matrixV(); 
  // cout << U << endl;
  // cout << "**" << endl;
  // cout << D << endl;
  // cout << "**" << endl;
  // cout << V << endl;
  
  // Aqper and Bqper
  Eigen::MatrixXd Aper(matrixSqrt2(U * U.transpose()));
  Eigen::MatrixXd Bper((matrixSqrt2(V * V.transpose())).transpose());
  
  
  // lambda and its varianve
  Eigen::MatrixXd BAper(Eigen::KroneckerProduct(Bper, Aper.transpose()));
  Eigen::VectorXd lambda (BAper * deltatilde);
  Eigen::MatrixXd varlambda (BAper * BAper.transpose());
  
  // statistic
  double stat((lambda.transpose() * varlambda.colPivHouseholderQr().solve(lambda))(0, 0));
  return Rcpp::List::create(_["stat"]  = stat, 
                            _["df"]    = Kdelta);
}


//[[Rcpp::export]]
Rcpp::List fEncompassingRedKP(const Eigen::MatrixXd& X1,
                              const Eigen::MatrixXd& qy1,
                              const Eigen::MatrixXd& Z1,
                              const Eigen::MatrixXd& W1,
                              const Eigen::VectorXd& e1,
                              const Eigen::VectorXd theta1,
                              const int& Kest1,
                              const Eigen::MatrixXd& X2,
                              const Eigen::MatrixXd& qy2,
                              const Eigen::MatrixXd& Z2,
                              const Eigen::MatrixXd& W2,
                              const Eigen::VectorXd& e2,
                              const Eigen::VectorXd theta2,
                              const int& Kest2,
                              const int& ngroup, //common to both models
                              const Eigen::ArrayXi& cumsn, //common to both models
                              const int& HAC = 0, //common to both models
                              const bool& full = false) {
  int ntau2(qy2.cols()), Kv2(ntau2 + X2.cols());
  // Delta and its variable
  Eigen::VectorXd deltaFull;
  Eigen::MatrixXd VardeltaFull;
  {
    Rcpp::List tp = fEncompassingRedDelta(X1, qy1, Z1, W1, e1, theta1, Kest1,
                                          X2, qy2, Z2, W2, e2, theta2, Kest2,
                                          ngroup, cumsn, HAC);
    deltaFull     = tp["delta"];
    VardeltaFull  = tp["Vdelta"];
  }
  
  // R matrix
  int Kdelta(ntau2);
  if (full) {
    Kdelta = Kv2;
  }
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(Kdelta, Kv2));
  R(Eigen::all, Eigen::seqN(0, Kdelta)) = Eigen::MatrixXd::Identity(Kdelta, Kdelta);
  
  // delta
  Eigen::VectorXd delta(R * deltaFull);
  
  // variance
  Eigen::MatrixXd Vardelta(R * VardeltaFull * R.transpose());
  
  // Normalization
  Eigen::LLT<Eigen::MatrixXd> tpOmega(Vardelta);
  Eigen::MatrixXd Omega(tpOmega.matrixL());
  Eigen::VectorXd deltatilde(Omega.colPivHouseholderQr().solve(delta));
  
  // Test if variance is I
  // Eigen::MatrixXd Vtp(Omega.colPivHouseholderQr().solve(Vardelta) * (Omega.transpose()).inverse());
  // cout << Vtp << endl;
  // cout<< deltatilde.transpose() << endl;
  
  // SDV decomposition of Theta
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(deltatilde.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd U = svd.matrixU(); //ntau * ntau
  Eigen::VectorXd d = svd.singularValues();
  Eigen::MatrixXd ddiag = d.asDiagonal();
  Eigen::MatrixXd D(1, Kdelta);
  D << ddiag, Eigen::MatrixXd::Zero(1, Kdelta - 1); 
  Eigen::MatrixXd V = svd.matrixV(); 
  // cout << U << endl;
  // cout << "**" << endl;
  // cout << D << endl;
  // cout << "**" << endl;
  // cout << V << endl;
  
  // Aqper and Bqper
  Eigen::MatrixXd Aper(matrixSqrt2(U * U.transpose()));
  Eigen::MatrixXd Bper((matrixSqrt2(V * V.transpose())).transpose());
  
  
  // lambda and its varianve
  Eigen::MatrixXd BAper(Eigen::KroneckerProduct(Bper, Aper.transpose()));
  Eigen::VectorXd lambda (BAper * deltatilde);
  Eigen::MatrixXd varlambda (BAper * BAper.transpose());
  
  // statistic
  double stat((lambda.transpose() * varlambda.colPivHouseholderQr().solve(lambda))(0, 0));
  return Rcpp::List::create(_["stat"]  = stat,
                            _["df"]    = Kdelta);
}

// Encompassing test, F method
//[[Rcpp::export]]
Rcpp::List fEncompassingStruc(const Eigen::MatrixXd& X1,
                              const Eigen::MatrixXd& qy1,
                              const Eigen::MatrixXd& Z1,
                              const Eigen::MatrixXd& W11,
                              const Eigen::MatrixXd& W21,
                              const Eigen::VectorXd& e1,
                              const Eigen::VectorXd theta1,
                              const Eigen::ArrayXi& idX11,
                              const Eigen::ArrayXi& idX21,
                              const int& Kest11,
                              const int& Kest21,
                              const Eigen::MatrixXd& X2,
                              const Eigen::MatrixXd& qy2,
                              const Eigen::MatrixXd& Z2,
                              const Eigen::MatrixXd& W12,
                              const Eigen::MatrixXd& W22,
                              const Eigen::VectorXd& e2,
                              const Eigen::VectorXd theta2,
                              const Eigen::ArrayXi& idX12,
                              const Eigen::ArrayXi& idX22,
                              const int& Kest12,
                              const int& Kest22,
                              const Eigen::ArrayXi& nIs, //common to both models
                              const Eigen::ArrayXi& Is, //common to both models
                              const int& ngroup, //common to both models
                              const Eigen::ArrayXi& cumsn, //common to both models
                              const int& HAC = 0, //common to both models
                              const bool& full = false) {
  int ntau2(qy2.cols()), Kv2(1 + ntau2 + idX22.size());
  // Delta and its variable
  Eigen::VectorXd deltaFull;
  Eigen::MatrixXd VardeltaFull;
  {
    Rcpp::List tp = fEncompassingStrucDelta(X1, qy1, Z1, W11, W21, e1, theta1, idX11, idX21, Kest11,
                                            Kest21, X2, qy2, Z2, W12, W22, e2, theta2, idX12, idX22,
                                            Kest12, Kest22, nIs, Is, ngroup, cumsn, HAC);
    deltaFull     = tp["delta"];
    VardeltaFull  = tp["Vdelta"];
  }
  
  
  // R matrix
  int df1(ntau2), df2(nIs.size() - Kest22);
  if (full) {
    df1 = Kv2;
  }
  Eigen::ArrayXi itheta(df1);
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(df1, Kv2));
  if (full) {
    R(Eigen::all, Eigen::seqN(0, df1)) = Eigen::MatrixXd::Identity(df1, df1);
    itheta << Eigen::ArrayXi::LinSpaced(1 + ntau2, 1, 1 + ntau2),  1 + ntau2 + idX22;
  } else {
    R(Eigen::all, Eigen::seqN(1, df1)) = Eigen::MatrixXd::Identity(df1, df1);
    itheta << Eigen::ArrayXi::LinSpaced(ntau2, 2, 1 + ntau2);
  }
  
  // delta
  Eigen::VectorXd delta(R * deltaFull);
  
  // variance
  Eigen::MatrixXd Vardelta(R * VardeltaFull * R.transpose());
  // cout<<"F - delta"<<endl;
  // cout<<delta.transpose()<<endl;
  // cout<<"F - Vardelta"<<endl;
  // cout<<Vardelta<<endl;
  
  // statistic
  double stat((delta.transpose()*ginv(Vardelta)*delta/R.rows())(0, 0));
  
  return Rcpp::List::create(_["stat"]   = stat,
                            _["df1"]    = df1, 
                            _["df2"]    = df2, 
                            _["delta"]  = delta,
                            _["Vdelta"] = Vardelta,
                            _["itheta"] = itheta);
}


//[[Rcpp::export]]
Rcpp::List fEncompassingRed(const Eigen::MatrixXd& X1,
                            const Eigen::MatrixXd& qy1,
                            const Eigen::MatrixXd& Z1,
                            const Eigen::MatrixXd& W1,
                            const Eigen::VectorXd& e1,
                            const Eigen::VectorXd theta1,
                            const int& Kest1,
                            const Eigen::MatrixXd& X2,
                            const Eigen::MatrixXd& qy2,
                            const Eigen::MatrixXd& Z2,
                            const Eigen::MatrixXd& W2,
                            const Eigen::VectorXd& e2,
                            const Eigen::VectorXd theta2,
                            const int& Kest2,
                            const int& ngroup, //common to both models
                            const Eigen::ArrayXi& cumsn, //common to both models
                            const int& HAC = 0, //common to both models
                            const bool& full = false) {
  int ntau2(qy2.cols()), Kv2(ntau2 + X2.cols());
  // Delta and its variable
  Eigen::VectorXd deltaFull;
  Eigen::MatrixXd VardeltaFull;
  {
    Rcpp::List tp = fEncompassingRedDelta(X1, qy1, Z1, W1, e1, theta1, Kest1,
                                          X2, qy2, Z2, W2, e2, theta2, Kest2,
                                          ngroup, cumsn, HAC);
    deltaFull     = tp["delta"];
    VardeltaFull  = tp["Vdelta"];
  }
  
  // R matrix
  int df1(ntau2), df2(qy2.rows() - Kest2);
  if (full) {
    df1 = Kv2;
  }
  Eigen::ArrayXi itheta(df1);
  itheta << Eigen::ArrayXi::LinSpaced(df1, 1, df1);
  Eigen::MatrixXd R(Eigen::MatrixXd::Zero(df1, Kv2));
  R(Eigen::all, Eigen::seqN(0, df1)) = Eigen::MatrixXd::Identity(df1, df1);
  
  // delta
  Eigen::VectorXd delta(R * deltaFull);
  
  // variance
  Eigen::MatrixXd Vardelta(R * VardeltaFull * R.transpose());
  
  // statistic
  double stat((delta.transpose()*ginv(Vardelta)*delta/R.rows())(0, 0));
  
  return Rcpp::List::create(_["stat"]   = stat,
                            _["df1"]    = df1, 
                            _["df2"]    = df2,  
                            _["delta"]  = delta,
                            _["Vdelta"] = Vardelta,
                            _["itheta"] = itheta);
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
