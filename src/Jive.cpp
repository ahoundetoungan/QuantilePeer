// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
//#define NDEBUG
// #include <RcppNumerical.h>
#include <RcppEigen.h>
// 
// typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
// typedef Eigen::Map<Eigen::VectorXd> MapVect;

// using namespace Numer;
using namespace Rcpp;
using namespace arma;
using namespace std;

////////////// Reduced-form models with independent errors
// This function implements the JIVE estimator (reduced-form) when errors are independent
//[[Rcpp::export]]
Rcpp::List fJIVE_redInd(const Eigen::VectorXd& y,
                           const Eigen::MatrixXd& V,
                           const Eigen::MatrixXd& ins,
                           const int& Kx,
                           const int& Kins,
                           const int& ntau,
                           const int& n,
                           const int& Kest, 
                           const int& HAC  = 0,
                           const bool& COV = true){
  int Kv(Kx + ntau);
  Eigen::MatrixXd ZZ(ins.transpose()*ins), iZZ(ZZ.inverse());
  //Projection matrix (pseudo) Pii
  Eigen::MatrixXd tZ(ins*iZZ);
  Eigen::VectorXd Pii((tZ.array()*ins.array()).rowwise().sum());
  Eigen::MatrixXd hV(V.array().colwise()*Pii.array()); // PiiXi
  Eigen::MatrixXd tV(V.array().colwise()/(1 - Pii.array())); // inv(1 - Pii)Xi
  Eigen::VectorXd ty(y.array()/(1 - Pii.array())); // inv(1 - Pii)yi
  Eigen::MatrixXd VtZ(V.transpose()*tZ);
  
  // estimate
  Eigen::MatrixXd tpH(VtZ*ins.transpose() - hV.transpose());
  Eigen::MatrixXd H(tpH*tV);
  Eigen::VectorXd parms(H.colPivHouseholderQr().solve(tpH*ty));
  
  // Covariance computation
  Eigen::VectorXd yhat(V*parms);
  Eigen::ArrayXd eh(y - yhat);
  double s2(R_NaN);
  if (HAC == 0) {
    s2 = eh.square().sum()/(n - Kest);
  }
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::ArrayXd teh(eh.array()/(1 - Pii.array()));
    Eigen::MatrixXd Zteh(ins.array().colwise()*teh), hVteh(hV.array().colwise()*teh), Vteh(V.array().colwise()*teh);
    Eigen::MatrixXd Zteh2(Zteh.transpose()*Zteh), hVtehZteh(hVteh.transpose()*Zteh);
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv, Kv));
    Eigen::VectorXd tp1, tp2;
    for (int k(0); k < Kins; ++ k) {
      tp1    = Vteh.transpose()*(tZ.col(k).array().square()).matrix();
      tp2    = Vteh.transpose()*(ins.col(k).array().square()).matrix();
      tp0   += (tp1*tp2.transpose()); 
      for (int l(k + 1); l < Kins; ++ l) {
        tp1  = Vteh.transpose()*(tZ.col(k).array()*tZ.col(l).array()).matrix();
        tp2  = Vteh.transpose()*(ins.col(k).array()*ins.col(l).array()).matrix();
        tp0 += (2*tp1*tp2.transpose()); 
      }
    }
    Eigen::MatrixXd tp3(hVtehZteh*VtZ.transpose());
    Eigen::MatrixXd Sigma(VtZ*Zteh2*VtZ.transpose() - tp3 - tp3.transpose() + tp0);
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, // _["Overident"] = stat,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}


// This function implements the JIVE2 estimator (reduced-form) when errors are independent
//[[Rcpp::export]]
Rcpp::List fJIVE2_redInd(const Eigen::VectorXd& y,
                            const Eigen::MatrixXd& V,
                            const Eigen::MatrixXd& ins,
                            const int& Kx,
                            const int& Kins,
                            const int& ntau,
                            const int& n,
                            const int& Kest, 
                            const int& HAC  = 0,
                            const bool& COV = true){
  int Kv(Kx + ntau);
  Eigen::MatrixXd ZZ(ins.transpose()*ins), iZZ(ZZ.inverse());
  //Projection matrix (pseudo) Pii
  Eigen::MatrixXd tZ(ins*iZZ);
  Eigen::VectorXd Pii((tZ.array()*ins.array()).rowwise().sum());
  Eigen::MatrixXd hV(V.array().colwise()*Pii.array()); 
  Eigen::MatrixXd VtZ(V.transpose()*tZ);
  
  // estimate
  Eigen::MatrixXd tpH(VtZ*ins.transpose() - hV.transpose());
  Eigen::MatrixXd H(tpH*V);
  Eigen::VectorXd parms(H.colPivHouseholderQr().solve(tpH*y));
  
  // Covariance computation
  Eigen::VectorXd yhat(V*parms);
  Eigen::ArrayXd eh(y - yhat);
  double s2(R_NaN);
  if (HAC == 0) {
    s2 = eh.square().sum()/(n - Kest);
  }
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::MatrixXd Zeh(ins.array().colwise()*eh), hVeh(hV.array().colwise()*eh), Veh(V.array().colwise()*eh);
    Eigen::MatrixXd Zeh2(Zeh.transpose()*Zeh), hVehZeh(hVeh.transpose()*Zeh);
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv, Kv));
    Eigen::VectorXd tp1, tp2;
    for (int k(0); k < Kins; ++ k) {
      tp1    = Veh.transpose()*(tZ.col(k).array().square()).matrix();
      tp2    = Veh.transpose()*(ins.col(k).array().square()).matrix();
      tp0   += (tp1*tp2.transpose()); 
      for (int l(k + 1); l < Kins; ++ l) {
        tp1  = Veh.transpose()*(tZ.col(k).array()*tZ.col(l).array()).matrix();
        tp2  = Veh.transpose()*(ins.col(k).array()*ins.col(l).array()).matrix();
        tp0 += (2*tp1*tp2.transpose()); 
      }
    }
    Eigen::MatrixXd tp3(hVehZeh*VtZ.transpose());
    Eigen::MatrixXd Sigma(VtZ*Zeh2*VtZ.transpose() - tp3 - tp3.transpose() + tp0);
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, // _["Overident"] = stat,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}

////////////// Reduced-form models with dependent errors
// This function implements the JIVE estimator (reduced-form) when errors are clustered
//[[Rcpp::export]]
Rcpp::List fJIVE_redClu(const Eigen::VectorXd& y,
                           const Eigen::MatrixXd& V,
                           const Eigen::MatrixXd& ins,
                           const Eigen::MatrixXi& igroup,
                           const Eigen::VectorXi& nvec,
                           const int& ngroup,
                           const int& Kx,
                           const int& Kins,
                           const int& ntau,
                           const int& n,
                           const int& Kest,
                           const bool& COV){
  int Kv(Kx + ntau);
  Eigen::MatrixXd ZZ(ins.transpose()*ins), iZZ(ZZ.inverse());
  
  Eigen::MatrixXd tZ(ins*iZZ), hV(n, Kv), tV(n, Kv);
  Eigen::VectorXd ty(n);
  Rcpp::List LiImPkk(ngroup);
  for (int r(0); r < ngroup; ++ r) {
    int n1(igroup(r, 0)), nr(nvec(r));
    Eigen::MatrixXd Vr(V(Eigen::seqN(n1, nr), Eigen::all));
    Eigen::MatrixXd Pkk(tZ(Eigen::seqN(n1, nr), Eigen::all)*ins(Eigen::seqN(n1, nr), Eigen::all).transpose());
    Eigen::MatrixXd iImPkk((Eigen::MatrixXd::Identity(nr, nr) - Pkk).inverse());
    
    LiImPkk[r] = iImPkk;
    hV(Eigen::seqN(n1, nr), Eigen::all) = Pkk*Vr;
    tV(Eigen::seqN(n1, nr), Eigen::all) = iImPkk*Vr;
    ty.segment(n1, nr)                  = iImPkk*y.segment(n1, nr);
  }
  Eigen::MatrixXd VtZ(V.transpose()*tZ);
  
  // estimate
  Eigen::MatrixXd tpH(VtZ*ins.transpose() - hV.transpose());
  Eigen::MatrixXd H(tpH*tV);
  Eigen::VectorXd parms(H.colPivHouseholderQr().solve(tpH*ty));
  
  // Covariance computation
  Eigen::VectorXd yhat(V*parms);
  Eigen::VectorXd eh(y - yhat);
  double s2(R_NaN);
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::MatrixXd Zteh2(Eigen::MatrixXd::Zero(Kins, Kins)); 
    Eigen::MatrixXd hVtehZteh(Eigen::MatrixXd::Zero(Kv, Kins));
    Eigen::MatrixXd Zteh(Kins, ngroup), tehtZ(ngroup, Kins);
    Rcpp::List LVtZ(ngroup), LZV(ngroup);
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), nr(nvec(r));
      Eigen::MatrixXd iImPkk = LiImPkk[r];
      Eigen::MatrixXd tZr(tZ(Eigen::seqN(n1, nr), Eigen::all));
      Eigen::MatrixXd Zr(ins(Eigen::seqN(n1, nr), Eigen::all));
      Eigen::MatrixXd Vr(V(Eigen::seqN(n1, nr), Eigen::all));
      Eigen::VectorXd tehr(iImPkk*eh.segment(n1, nr));
      Eigen::VectorXd Ztehr(Zr.transpose()*tehr);
      
      Zteh2       += (Ztehr*Ztehr.transpose());
      hVtehZteh   += ((hV(Eigen::seqN(n1, nr), Eigen::all).transpose()*tehr)*Ztehr.transpose());
      Eigen::MatrixXd VtZr(Vr.transpose()*tZr);
      Eigen::MatrixXd ZVr(Zr.transpose()*Vr);
      Zteh.col(r)  = Ztehr;
      tehtZ.row(r) = (tehr.transpose()*tZr);
      LVtZ[r]      = VtZr;
      LZV[r]       = ZVr;
    }
    // Compute the double sum
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv, Kv));
    for (int r(0); r < ngroup; ++ r) {
      Eigen::MatrixXd VtZr   = LVtZ[r];
      Eigen::MatrixXd ZVr    = LZV[r];
      // tp0   += (VtZr*Zteh.col(r)*tehtZ.row(r)*ZVr);
      for (int s(0); s < ngroup; ++ s) {
        Eigen::MatrixXd ZVs  = LZV[s];
        tp0 += (VtZr*Zteh.col(s)*tehtZ.row(r)*ZVs);
      }
    }
    Eigen::MatrixXd tp3(hVtehZteh*VtZ.transpose());
    Eigen::MatrixXd Sigma(VtZ*Zteh2*VtZ.transpose() - tp3 - tp3.transpose() + tp0);
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, // _["Overident"] = stat,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}

// This function implements the JIVE2 estimator (reduced-form) when errors are clustered
//[[Rcpp::export]]
Rcpp::List fJIVE2_redClu(const Eigen::VectorXd& y,
                            const Eigen::MatrixXd& V,
                            const Eigen::MatrixXd& ins,
                            const Eigen::MatrixXd& igroup,
                            const Eigen::VectorXi& nvec,
                            const int& ngroup,
                            const int& Kx,
                            const int& Kins,
                            const int& ntau,
                            const int& n,
                            const int& Kest,
                            const bool& COV){
  int Kv(Kx + ntau);
  Eigen::MatrixXd ZZ(ins.transpose()*ins), iZZ(ZZ.inverse());
  
  Eigen::MatrixXd tZ(ins*iZZ), hV(n, Kv);
  for (int r(0); r < ngroup; ++ r) {
    int n1(igroup(r, 0)), nr(nvec(r));
    Eigen::MatrixXd Pkk((tZ(Eigen::seqN(n1, nr), Eigen::all))*(ins(Eigen::seqN(n1, nr), Eigen::all)).transpose());
    hV(Eigen::seqN(n1, nr), Eigen::all) = Pkk*V(Eigen::seqN(n1, nr), Eigen::all);
  }
  Eigen::MatrixXd VtZ(V.transpose()*tZ);
  
  // estimate
  Eigen::MatrixXd tpH(VtZ*ins.transpose() - hV.transpose());
  Eigen::MatrixXd H(tpH*V);
  Eigen::VectorXd parms(H.colPivHouseholderQr().solve(tpH*y));
  
  // Covariance computation
  Eigen::VectorXd yhat(V*parms);
  Eigen::VectorXd eh(y - yhat);
  double s2(R_NaN);
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::MatrixXd Zeh2(Eigen::MatrixXd::Zero(Kins, Kins)); 
    Eigen::MatrixXd hVehZeh(Eigen::MatrixXd::Zero(Kv, Kins));
    Eigen::MatrixXd Zeh(Kins, ngroup), ehtZ(ngroup, Kins);
    Rcpp::List LVtZ(ngroup), LZV(ngroup);
    for (int r(0); r < ngroup; ++ r) {
      int n1(igroup(r, 0)), nr(nvec(r));
      Eigen::MatrixXd tZr(tZ(Eigen::seqN(n1, nr), Eigen::all));
      Eigen::MatrixXd Zr(ins(Eigen::seqN(n1, nr), Eigen::all));
      Eigen::MatrixXd Vr(V(Eigen::seqN(n1, nr), Eigen::all));
      Eigen::MatrixXd hVr(hV(Eigen::seqN(n1, nr), Eigen::all));
      Eigen::VectorXd ehr(eh.segment(n1, nr));
      
      Eigen::VectorXd Zehr(Zr.transpose()*ehr);
      Eigen::VectorXd hVehr(hVr.transpose()*ehr);
      Eigen::RowVectorXd ehtZr(ehr.transpose()*tZr);
      
      Zeh2       += (Zehr*Zehr.transpose());
      hVehZeh    += (hVehr*Zehr.transpose());
      Eigen::MatrixXd VtZr(Vr.transpose()*tZr);
      Eigen::MatrixXd ZVr(Zr.transpose()*Vr);
      Zeh.col(r)  = Zehr;
      ehtZ.row(r) = ehtZr;
      LVtZ[r]     = VtZr;
      LZV[r]      = ZVr;
    }
    // Compute the double sum
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv, Kv));
    for (int r(0); r < ngroup; ++ r) {
      Eigen::MatrixXd VtZr   = LVtZ[r];
      Eigen::MatrixXd ZVr    = LZV[r];
      // tp0   += (VtZr*Zeh.col(r)*ehtZ.row(r)*ZVr);
      for (int s(0); s < ngroup; ++ s) {
        Eigen::MatrixXd ZVs  = LZV[s];
        tp0 += (VtZr*Zeh.col(s)*ehtZ.row(r)*ZVs);
      }
    }
    Eigen::MatrixXd tp3(hVehZeh*VtZ.transpose());
    Eigen::MatrixXd Sigma(VtZ*Zeh2*VtZ.transpose() - tp3 - tp3.transpose() + tp0);
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  return Rcpp::List::create(_["parms"] = parms, _["Vpa"] = Vpa, // _["Overident"] = stat,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}

////////////// Structural models with independent errors
// This function implements the JIVE estimator for the structural model when errors are independent
//[[Rcpp::export]]
Rcpp::List fJIVE_strucInd(const Eigen::VectorXd& y,
                             const Eigen::MatrixXd& X,
                             const Eigen::MatrixXd& qy,
                             const Eigen::MatrixXd& ins,
                             const Eigen::VectorXi& idX1,
                             const Eigen::VectorXi& idX2,
                             const Eigen::VectorXi& nIs,
                             const Eigen::VectorXi& Is,
                             const int& Kx1,
                             const int& Kx2,
                             const int& Kins,
                             const int& ntau,
                             const int& n,
                             const int& Kest, 
                             const int& HAC  = 0,
                             const bool& COV = true){
  int n_iso(Is.size()), n_niso(n - n_iso), Kv2(1 + ntau + Kx2), Kv(Kx1 + Kv2);
  // First stage
  Eigen::VectorXd y1(y(Is));
  Eigen::MatrixXd X1(X(Is, idX1));
  Eigen::MatrixXd XX1(X1.transpose()*X1);
  Eigen::VectorXd b(XX1.colPivHouseholderQr().solve(X1.transpose()*y1));
  
  // Second stage
  Eigen::VectorXd Xb(X(Eigen::all, idX1)*b), Xb1(Xb(Is)), Xb2(Xb(nIs)), y2(y(nIs));
  Eigen::MatrixXd X2(X(nIs, idX2));
  Eigen::MatrixXd X21(X(nIs, idX1));
  Eigen::MatrixXd V2(n_niso, 1 + ntau + Kx2);
  V2 << Xb2, qy(nIs, Eigen::all), X2;
  Eigen::MatrixXd Z2(n_niso, 1 + Kins);
  Z2 << Xb2, ins(nIs, Eigen::all);
  Eigen::MatrixXd ZZ2(Z2.transpose()*Z2), iZZ2(ZZ2.inverse()); 
  //Projection matrix (pseudo) Pii
  Eigen::MatrixXd tZ2(Z2*iZZ2);
  Eigen::VectorXd Pii2((tZ2.array()*Z2.array()).rowwise().sum());
  Eigen::MatrixXd hV2(V2.array().colwise()*Pii2.array()); // PiiXi
  Eigen::MatrixXd tV2(V2.array().colwise()/(1 - Pii2.array())); // inv(1 - Pii)Xi
  Eigen::VectorXd ty2(y2.array()/(1 - Pii2.array())); // inv(1 - Pii)yi
  Eigen::MatrixXd VtZ2(V2.transpose()*tZ2);
  Eigen::MatrixXd tX21(X21.array().colwise()/(1 - Pii2.array()));
  
  Eigen::MatrixXd tpH(VtZ2*Z2.transpose() - hV2.transpose());
  Eigen::MatrixXd H2(tpH*tV2);
  Eigen::VectorXd lambda(H2.colPivHouseholderQr().solve(tpH*ty2));
  
  // Covariance computation
  Eigen::VectorXd y1hat(Xb1), y2hat(V2*lambda), yhat(n);
  yhat(Is) = y1hat; yhat(nIs) = y2hat;
  Eigen::ArrayXd e1(y1 - y1hat), e2(y2 - y2hat);
  Eigen::ArrayXd e = y - yhat;
  double s2(R_NaN);
  if (HAC == 0) {
    s2 = e.square().sum()/(n - Kest);
  }
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kv, Kv));
    H.block(0, 0, Kx1, Kx1)   = XX1;
    H(Eigen::seqN(Kx1, Kv2), Eigen::all) << lambda(0)*tpH*tX21, H2;
    Eigen::MatrixXd Sigma(Eigen::MatrixXd::Zero(Kv, Kv));
    if (HAC == 0) {
      Sigma.block(0, 0, Kx1, Kx1) = s2*XX1;
    } else {
      Eigen::MatrixXd Xe1(X1.array().colwise()*e1);
      Sigma.block(0, 0, Kx1, Kx1) = Xe1.transpose()*Xe1;
    }
    Eigen::ArrayXd teh2(e2.array()/(1 - Pii2.array()));
    Eigen::MatrixXd Zteh2(Z2.array().colwise()*teh2), hVteh2(hV2.array().colwise()*teh2), Vteh2(V2.array().colwise()*teh2);
    Eigen::MatrixXd Zteh22(Zteh2.transpose()*Zteh2), hVtehZteh2(hVteh2.transpose()*Zteh2);
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv2, Kv2));
    Eigen::VectorXd tp1, tp2;
    for (int k(0); k < (1 + Kins); ++ k) {
      tp1    = Vteh2.transpose()*(tZ2.col(k).array().square()).matrix();
      tp2    = Vteh2.transpose()*(Z2.col(k).array().square()).matrix();
      tp0   += (tp1*tp2.transpose()); 
      for (int l(k + 1); l < (1 + Kins); ++ l) {
        tp1  = Vteh2.transpose()*(tZ2.col(k).array()*tZ2.col(l).array()).matrix();
        tp2  = Vteh2.transpose()*(Z2.col(k).array()*Z2.col(l).array()).matrix();
        tp0 += (2*tp1*tp2.transpose()); 
      }
    }
    Eigen::MatrixXd tp3(hVtehZteh2*VtZ2.transpose());
    Sigma.block(Kx1, Kx1, Kv2, Kv2) = (VtZ2*Zteh22*VtZ2.transpose() - tp3 - tp3.transpose() + tp0);
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}


// This function implements the JIVE estimator for the structural model when errors are independent
//[[Rcpp::export]]
Rcpp::List fJIVE2_strucInd(const Eigen::VectorXd& y,
                              const Eigen::MatrixXd& X,
                              const Eigen::MatrixXd& qy,
                              const Eigen::MatrixXd& ins,
                              const Eigen::VectorXi& idX1,
                              const Eigen::VectorXi& idX2,
                              const Eigen::VectorXi& nIs,
                              const Eigen::VectorXi& Is,
                              const int& Kx1,
                              const int& Kx2,
                              const int& Kins,
                              const int& ntau,
                              const int& n,
                              const int& Kest, 
                              const int& HAC  = 0,
                              const bool& COV = true){
  int n_iso(Is.size()), n_niso(n - n_iso), Kv2(1 + ntau + Kx2), Kv(Kx1 + Kv2);
  // First stage
  Eigen::VectorXd y1(y(Is));
  Eigen::MatrixXd X1(X(Is, idX1));
  Eigen::MatrixXd XX1(X1.transpose()*X1);
  Eigen::VectorXd b(XX1.colPivHouseholderQr().solve(X1.transpose()*y1));
  
  // Second stage
  Eigen::VectorXd Xb(X(Eigen::all, idX1)*b), Xb1(Xb(Is)), Xb2(Xb(nIs)), y2(y(nIs));
  Eigen::MatrixXd X2(X(nIs, idX2));
  Eigen::MatrixXd X21(X(nIs, idX1));
  Eigen::MatrixXd V2(n_niso, 1 + ntau + Kx2);
  V2 << Xb2, qy(nIs, Eigen::all), X2;
  Eigen::MatrixXd Z2(n_niso, 1 + Kins);
  Z2 << Xb2, ins(nIs, Eigen::all);
  Eigen::MatrixXd ZZ2(Z2.transpose()*Z2), iZZ2(ZZ2.inverse()); 
  //Projection matrix (pseudo) Pii
  Eigen::MatrixXd tZ2(Z2*iZZ2);
  Eigen::VectorXd Pii2((tZ2.array()*Z2.array()).rowwise().sum());
  Eigen::MatrixXd hV2(V2.array().colwise()*Pii2.array()); // PiiXi
  Eigen::MatrixXd VtZ2(V2.transpose()*tZ2);
  
  Eigen::MatrixXd tpH(VtZ2*Z2.transpose() - hV2.transpose());
  Eigen::MatrixXd H2(tpH*V2);
  Eigen::VectorXd lambda(H2.colPivHouseholderQr().solve(tpH*y2));
  
  // Covariance computation
  Eigen::VectorXd y1hat(Xb1), y2hat(V2*lambda), yhat(n);
  yhat(Is) = y1hat; yhat(nIs) = y2hat;
  Eigen::ArrayXd e1(y1 - y1hat), e2(y2 - y2hat);
  Eigen::ArrayXd e = y - yhat;
  double s2(R_NaN);
  if (HAC == 0) {
    s2 = e.square().sum()/(n - Kest);
  }
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kv, Kv));
    Eigen::MatrixXd Sigma(Eigen::MatrixXd::Zero(Kv, Kv));
    H.block(0, 0, Kx1, Kx1)   = XX1;
    H(Eigen::seqN(Kx1, Kv2), Eigen::all) << lambda(0)*tpH*X21, H2;
    if (HAC == 0) {
      Sigma.block(0, 0, Kx1, Kx1) = s2*XX1;
    } else {
      Eigen::MatrixXd Xe1(X1.array().colwise()*e1);
      Sigma.block(0, 0, Kx1, Kx1) = Xe1.transpose()*Xe1;
    }
    Eigen::MatrixXd Ze2(Z2.array().colwise()*e2), hVe2(hV2.array().colwise()*e2), Ve2(V2.array().colwise()*e2);
    Eigen::MatrixXd Ze22(Ze2.transpose()*Ze2), hVeZe2(hVe2.transpose()*Ze2);
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv2, Kv2));
    Eigen::VectorXd tp1, tp2;
    for (int k(0); k < (1 + Kins); ++ k) {
      tp1    = Ve2.transpose()*(tZ2.col(k).array().square()).matrix();
      tp2    = Ve2.transpose()*(Z2.col(k).array().square()).matrix();
      tp0   += (tp1*tp2.transpose());
      for (int l(k + 1); l < (1 + Kins); ++ l) {
        tp1  = Ve2.transpose()*(tZ2.col(k).array()*tZ2.col(l).array()).matrix();
        tp2  = Ve2.transpose()*(Z2.col(k).array()*Z2.col(l).array()).matrix();
        tp0 += (2*tp1*tp2.transpose()); 
      }
    }
    Eigen::MatrixXd tp3(hVeZe2*VtZ2.transpose());
    Sigma.block(Kx1, Kx1, Kv2, Kv2) = (VtZ2*Ze22*VtZ2.transpose() - tp3 - tp3.transpose() + tp0);
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}

////////////// Structural models with clustered errors
// This function implements the JIVE estimator for the structural model when errors are clustered
//[[Rcpp::export]]
Rcpp::List fJIVE_strucClu(const Eigen::VectorXd& y,
                             const Eigen::MatrixXd& X,
                             const Eigen::MatrixXd& qy,
                             const Eigen::MatrixXd& ins,
                             const Eigen::VectorXi& idX1,
                             const Eigen::VectorXi& idX2,
                             const Eigen::VectorXi& nIs,
                             const Eigen::VectorXi& Is,
                             Rcpp::List LnIs,
                             Rcpp::List LIs,
                             const Eigen::MatrixXi& igroup,
                             const Eigen::VectorXi& nvecnIs,
                             const Eigen::VectorXi& hasnIs,
                             const Eigen::VectorXi& hasIs,
                             const int& ngroup,
                             const int& Kx1,
                             const int& Kx2,
                             const int& Kins,
                             const int& ntau,
                             const int& n,
                             const int& Kest, 
                             const bool& COV = true){
  int Kv2(1 + ntau + Kx2), Kv(Kx1 + Kv2);
  // First stage
  Eigen::VectorXd y1(y(Is));
  Eigen::MatrixXd X1(X(Is, idX1));
  Eigen::MatrixXd XX1(X1.transpose()*X1);
  Eigen::VectorXd b(XX1.colPivHouseholderQr().solve(X1.transpose()*y1));
  
  // Second stage
  Eigen::VectorXd Xb(X(Eigen::all, idX1)*b);
  Eigen::MatrixXd V(n, Kv2);
  V << Xb, qy, X(Eigen::all, idX2);
  Eigen::MatrixXd Z(n, 1 + Kins);
  Z << Xb, ins;
  Eigen::MatrixXd ZZ2(Z(nIs, Eigen::all).transpose()*Z(nIs, Eigen::all)), iZZ2(ZZ2.inverse()); 
  //Projection matrix (pseudo) Pii
  Eigen::MatrixXd tZ(Z*iZZ2), hV(n, Kv2), tV(n, Kv2), tX1(n, Kx1);
  Eigen::VectorXd ty(n);
  Rcpp::List LiImPkk(ngroup);
  for (int r(0); r < ngroup; ++ r) {
    if (hasnIs(r) == 1) {
      int nrnIso(nvecnIs(r));
      Eigen::VectorXi nIsr = LnIs[r];
      Eigen::MatrixXd Vr(V(nIsr, Eigen::all));
      Eigen::MatrixXd Pkk(tZ(nIsr, Eigen::all)*Z(nIsr, Eigen::all).transpose());
      Eigen::MatrixXd iImPkk((Eigen::MatrixXd::Identity(nrnIso, nrnIso) - Pkk).inverse());
      
      LiImPkk[r]            = iImPkk;
      hV(nIsr, Eigen::all)  = Pkk*Vr;
      tV(nIsr, Eigen::all)  = iImPkk*Vr;
      ty(nIsr)              = iImPkk*y(nIsr);
      tX1(nIsr, Eigen::all) = iImPkk*X(nIsr, idX1);
    }
  }
  Eigen::MatrixXd VtZ2(V(nIs, Eigen::all).transpose()*tZ(nIs, Eigen::all));
  
  // estimate
  Eigen::MatrixXd tpH(VtZ2*Z(nIs, Eigen::all).transpose() - hV(nIs, Eigen::all).transpose());
  Eigen::MatrixXd H2(tpH*tV(nIs, Eigen::all));
  Eigen::VectorXd lambda(H2.colPivHouseholderQr().solve(tpH*ty(nIs)));
  
  // Covariance computation
  Eigen::VectorXd y1hat(Xb(Is)), y2hat(V(nIs, Eigen::all)*lambda), yhat(n);
  yhat(Is) = y1hat; yhat(nIs) = y2hat;
  Eigen::VectorXd eh = y - yhat;
  double s2(R_NaN);
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kv, Kv));
    H.block(0, 0, Kx1, Kx1)   = XX1;
    H(Eigen::seqN(Kx1, Kv2), Eigen::all) << lambda(0)*tpH*tX1(nIs, Eigen::all), H2;
    
    Eigen::MatrixXd Sigma1(Eigen::MatrixXd::Zero(Kx1, Kx1));
    Eigen::MatrixXd Sigma12(Eigen::MatrixXd::Zero(Kx1, Kv2));
    Eigen::MatrixXd Zteh22(Eigen::MatrixXd::Zero(Kins + 1, Kins + 1)); 
    Eigen::MatrixXd hVtehZteh2(Eigen::MatrixXd::Zero(Kv2, Kins + 1));
    Eigen::MatrixXd Zteh(Kins + 1, ngroup), tehtZ(ngroup, Kins + 1), Xeh(Kx1, ngroup);
    Rcpp::List LVtZ(ngroup), LZV(ngroup);
    for (int r(0); r < ngroup; ++ r) {
      if (hasIs(r) == 1) {
        Eigen::VectorXi Isr = LIs[r];
        Eigen::VectorXd Xehr(X(Isr, idX1).transpose()*eh(Isr));
        Xeh.col(r) = Xehr;
        Sigma1    += (Xehr*Xehr.transpose());
      }
      if (hasnIs(r) == 1) {
        // int n1(igroup(r, 0)), nr(nvec(r));
        Eigen::VectorXi nIsr = LnIs[r];
        Eigen::MatrixXd iImPkk = LiImPkk[r];
        Eigen::MatrixXd tZr(tZ(nIsr, Eigen::all));
        Eigen::MatrixXd Zr(Z(nIsr, Eigen::all));
        Eigen::MatrixXd Vr(V(nIsr, Eigen::all));
        Eigen::MatrixXd hVr(hV(nIsr, Eigen::all));
        Eigen::VectorXd tehr(iImPkk*eh(nIsr));
        Eigen::VectorXd Ztehr(Zr.transpose()*tehr);

        Zteh22       += (Ztehr*Ztehr.transpose());
        hVtehZteh2   += ((hVr.transpose()*tehr)*Ztehr.transpose());
        Eigen::MatrixXd VtZr(Vr.transpose()*tZr);
        Eigen::MatrixXd ZVr(Zr.transpose()*Vr);
        Zteh.col(r)  = Ztehr;
        tehtZ.row(r) = (tehr.transpose()*tZr);
        LVtZ[r]      = VtZr;
        LZV[r]       = ZVr;
        if (hasIs(r) == 1) {
          Eigen::VectorXi Isr = LIs[r];
          Sigma12   += (Xeh.col(r)*(tehr.transpose()*(Zr*VtZr.transpose() - hVr)));
        }
      }
    }

    // Compute the double sum in Sigma2
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv2, Kv2));
    for (int r(0); r < ngroup; ++ r) {
      if (hasnIs(r) == 1) {
        Eigen::MatrixXd VtZr   = LVtZ[r];
        Eigen::MatrixXd ZVr    = LZV[r];
        // tp0   += (VtZr*Zteh.col(r)*tehtZ.row(r)*ZVr); not symmetric by group
        for (int s(0); s < ngroup; ++ s) {
          if (hasnIs(s) == 1) {
            Eigen::MatrixXd ZVs  = LZV[s];
            tp0 += (VtZr*Zteh.col(s)*tehtZ.row(r)*ZVs);
          }
        }
      }
    }
    Eigen::MatrixXd tp3(hVtehZteh2*VtZ2.transpose());
    Eigen::MatrixXd Sigma2(VtZ2*Zteh22*VtZ2.transpose() - tp3 - tp3.transpose() + tp0);
    Eigen::MatrixXd Sigma(Kv, Kv);
    Sigma << Sigma1, Sigma12,
             Sigma12.transpose(), Sigma2;
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}


// This function implements the JIVE2 estimator for the structural model when errors are clustering
//[[Rcpp::export]]
Rcpp::List fJIVE2_strucClu(const Eigen::VectorXd& y,
                             const Eigen::MatrixXd& X,
                             const Eigen::MatrixXd& qy,
                             const Eigen::MatrixXd& ins,
                             const Eigen::VectorXi& idX1,
                             const Eigen::VectorXi& idX2,
                             const Eigen::VectorXi& nIs,
                             const Eigen::VectorXi& Is,
                             Rcpp::List LnIs,
                             Rcpp::List LIs,
                             const Eigen::MatrixXi& igroup,
                             const Eigen::VectorXi& hasnIs,
                             const Eigen::VectorXi& hasIs,
                             const int& ngroup,
                             const int& Kx1,
                             const int& Kx2,
                             const int& Kins,
                             const int& ntau,
                             const int& n,
                             const int& Kest, 
                             const bool& COV = true){
  int Kv2(1 + ntau + Kx2), Kv(Kx1 + Kv2);
  // First stage
  Eigen::VectorXd y1(y(Is));
  Eigen::MatrixXd X1(X(Is, idX1));
  Eigen::MatrixXd XX1(X1.transpose()*X1);
  Eigen::VectorXd b(XX1.colPivHouseholderQr().solve(X1.transpose()*y1));
  
  // Second stage
  Eigen::VectorXd Xb(X(Eigen::all, idX1)*b);
  Eigen::MatrixXd V(n, Kv2);
  V << Xb, qy, X(Eigen::all, idX2);
  Eigen::MatrixXd Z(n, 1 + Kins);
  Z << Xb, ins;
  Eigen::MatrixXd ZZ2(Z(nIs, Eigen::all).transpose()*Z(nIs, Eigen::all)), iZZ2(ZZ2.inverse()); 
  //Projection matrix (pseudo) Pii
  Eigen::MatrixXd tZ(Z*iZZ2), hV(n, Kv2);
  for (int r(0); r < ngroup; ++ r) {
    if (hasnIs(r) == 1) {
      Eigen::VectorXi nIsr = LnIs[r];
      hV(nIsr, Eigen::all)  = (tZ(nIsr, Eigen::all)*Z(nIsr, Eigen::all).transpose())*(V(nIsr, Eigen::all));
    }
  }
  Eigen::MatrixXd VtZ2(V(nIs, Eigen::all).transpose()*tZ(nIs, Eigen::all));
  
  // estimate
  Eigen::MatrixXd tpH(VtZ2*Z(nIs, Eigen::all).transpose() - hV(nIs, Eigen::all).transpose());
  Eigen::MatrixXd H2(tpH*V(nIs, Eigen::all));
  Eigen::VectorXd lambda(H2.colPivHouseholderQr().solve(tpH*y(nIs)));
  
  // Covariance computation
  Eigen::VectorXd y1hat(Xb(Is)), y2hat(V(nIs, Eigen::all)*lambda), yhat(n);
  yhat(Is) = y1hat; yhat(nIs) = y2hat;
  Eigen::VectorXd eh = y - yhat;
  double s2(R_NaN);
  Eigen::MatrixXd Vpa(Eigen::MatrixXd::Zero(Kv, Kv));
  if (COV) {
    Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kv, Kv));
    Eigen::MatrixXd Sigma(Kv, Kv);
    H.block(0, 0, Kx1, Kx1)   = XX1;
    H(Eigen::seqN(Kx1, Kv2), Eigen::all) << lambda(0)*tpH*X(nIs, Eigen::all), H2;
    
    Eigen::MatrixXd Sigma1(Eigen::MatrixXd::Zero(Kx1, Kx1));
    Eigen::MatrixXd Sigma12(Eigen::MatrixXd::Zero(Kx1, Kv2));
    Eigen::MatrixXd Zeh22(Eigen::MatrixXd::Zero(Kins + 1, Kins + 1)); 
    Eigen::MatrixXd hVehZeh2(Eigen::MatrixXd::Zero(Kv2, Kins + 1));
    Eigen::MatrixXd Zeh(Kins + 1, ngroup), ehtZ(ngroup, Kins + 1), Xeh(Kx1, ngroup);
    Rcpp::List LVtZ(ngroup), LZV(ngroup);
    for (int r(0); r < ngroup; ++ r) {
      if (hasIs(r) == 1) {
        Eigen::VectorXi Isr = LIs[r];
        Eigen::VectorXd Xehr(X(Isr, idX1).transpose()*eh(Isr));
        Xeh.col(r) = Xehr;
        Sigma1    += (Xehr*Xehr.transpose());
      }
      if (hasnIs(r) == 1) {
        // int n1(igroup(r, 0)), nr(nvec(r));
        Eigen::VectorXi nIsr = LnIs[r];
        Eigen::MatrixXd tZr(tZ(nIsr, Eigen::all));
        Eigen::MatrixXd Zr(Z(nIsr, Eigen::all));
        Eigen::MatrixXd Vr(V(nIsr, Eigen::all));
        Eigen::MatrixXd hVr(hV(nIsr, Eigen::all));
        Eigen::VectorXd ehr(eh(nIsr));
        Eigen::VectorXd Zehr(Zr.transpose()*ehr);
        
        Zeh22       += (Zehr*Zehr.transpose());
        hVehZeh2    += ((hVr.transpose()*ehr)*Zehr.transpose());
        Eigen::MatrixXd VtZr(Vr.transpose()*tZr);
        Eigen::MatrixXd ZVr(Zr.transpose()*Vr);
        Zeh.col(r)   = Zehr;
        ehtZ.row(r)  = (ehr.transpose()*tZr);
        LVtZ[r]      = VtZr;
        LZV[r]       = ZVr;
        if (hasIs(r) == 1) {
          Eigen::VectorXi Isr = LIs[r];
          Sigma12   += (Xeh.col(r)*(ehr.transpose()*(Zr*VtZr.transpose() - hVr)));
        }
      }
    }
    
    // Compute the double sum in Sigma2
    Eigen::MatrixXd tp0(Eigen::MatrixXd::Zero(Kv2, Kv2));
    for (int r(0); r < ngroup; ++ r) {
      if (hasnIs(r) == 1) {
        Eigen::MatrixXd VtZr   = LVtZ[r];
        Eigen::MatrixXd ZVr    = LZV[r];
        // tp0   += (VtZr*Zeh.col(r)*ehtZ.row(r)*ZVr);
        for (int s(0); s < ngroup; ++ s) {
          if (hasnIs(s) == 1) {
            Eigen::MatrixXd ZVs  = LZV[s];
            tp0 += (VtZr*Zeh.col(s)*ehtZ.row(r)*ZVs);
          }
        }
      }
    }
    Eigen::MatrixXd tp3(hVehZeh2*VtZ2.transpose());
    Eigen::MatrixXd Sigma2(VtZ2*Zeh22*VtZ2.transpose() - tp3 - tp3.transpose() + tp0);
    Sigma << Sigma1, Sigma12,
             Sigma12.transpose(), Sigma2;
    Eigen::MatrixXd tp4(H.colPivHouseholderQr().solve(Sigma)); // inv(H)*Sigma;
    Vpa      = H.colPivHouseholderQr().solve(tp4.transpose());// inv(H)*Sigma*inv(H');
  }
  
  return Rcpp::List::create(_["beta"] = b, _["lambda"] = lambda, _["Vpa"] = Vpa,
                            _["yhat"] = yhat, _["sigma2"] = s2);
}
