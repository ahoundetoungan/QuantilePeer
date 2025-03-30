// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace arma;
using namespace Eigen;


///////////////////////////////////////////////
///////////////////////////////////////////////
////////////////// Functions //////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////

///////////////////////////////////////////////
///////////////// Prep of Data ////////////////
///////////////////////////////////////////////

// [[Rcpp::export]]
arma::mat fCESdata(const arma::mat& X, 
                   const arma::vec& y, 
                   const arma::vec& z, 
                   const Rcpp::List& G,
                   const Rcpp::List& friendindex,
                   const arma::umat& igroup,
                   const arma::uvec& frzeroy, 
                   const arma::uvec& frzeroz, 
                   const Rcpp::List& lIs,//in selection 
                   const Rcpp::List& lnIs,//in selection 
                   const IntegerVector& nvec, 
                   const Eigen::MatrixXd& yFMiMa,
                   const Eigen::MatrixXd& zFMiMa,
                   const int& n,
                   const int& Kx, 
                   const int& ngroup,
                   const double& rho, 
                   const int& FEnum,
                   const bool& deriv) {
  vec Gy(n, fill::zeros), dGy(n, fill::zeros), Gz(n, fill::zeros), dGz(n, fill::zeros), ddGz(n, fill::zeros);
  
  for (int r(0); r < ngroup; r++) {
    // Extract data for group r
    int nm(nvec(r)), n1(igroup(r, 0)), n2(igroup(r, 1));
    arma::mat Gr = G[r];
    arma::mat Xm(X.rows(n1, n2));
    arma::vec ym(y.subvec(n1, n2)), zm(z.subvec(n1, n2));
    Rcpp::List frindexm(friendindex[r]);
    
    for (int i(0); i < nm; i++) {
      uvec fri = frindexm[i];
      if (!fri.is_empty()) {
        if (rho == R_PosInf) {
          Gy[n1 + i] = max(ym.elem(fri));
          Gz[n1 + i] = max(zm.elem(fri)); 
        } else if (rho == R_NegInf) {
          Gy[n1 + i] = min(ym.elem(fri));
          Gz[n1 + i] = min(zm.elem(fri)); 
        } else if (rho < 0) {
          arma::vec Gri(Gr.row(i).t());
          
          if (frzeroy[n1 + i] == 0) {
            arma::vec ymf(ym.elem(fri)/yFMiMa(n1 + i, 0)); // y of friends but scaled
            arma::vec tpy(Gri.elem(fri) % pow(ymf, rho));
            double stpy(sum(tpy));
            Gy(n1 + i)  = yFMiMa(n1 + i, 0)*pow(stpy, 1.0 / rho);
            if (deriv) {
              dGy(n1 + i) = Gy(n1 + i)*(sum(tpy % log(ymf))/(rho * stpy) - log(stpy)/pow(rho, 2));// Derivative of Gy for gradient computation
            }
          }
          
          if (frzeroz[n1 + i] == 0) {
            arma::vec zmf(zm.elem(fri)/zFMiMa(n1 + i, 0)); // z of friends but scaled
            arma::vec tpz(Gri.elem(fri) % pow(zmf, rho));
            double stpz(sum(tpz));
            Gz(n1 + i)  = zFMiMa(n1 + i, 0)*pow(stpz, 1.0 / rho);
            double tp(sum(tpz % log(zmf)));
            dGz(n1 + i)  = Gz(n1 + i)*(tp/(rho * stpz) - log(stpz)/pow(rho, 2));
            // if ((n1 + i) == 4) {
            //   cout<< frzeroz[n1 + i] <<endl;
            //   cout<< stpz <<endl;
            //   cout<< pow(stpz, 1.0 / rho) <<endl;
            //   cout<< Gz(n1 + i) <<endl;
            //   cout<< zm.elem(fri) <<endl;
            // }
            if (deriv) {
              ddGz(n1 + i) = Gz(n1 + i)*(2*log(stpz)/pow(rho, 3) - 2*tp/(stpz * pow(rho, 2)) + 
                (sum(tpz % pow(log(zmf), 2))*stpz - pow(tp, 2))/(rho * pow(stpz, 2)) + 
                pow(tp/(rho * stpz) - log(stpz)/pow(rho, 2), 2));
            }
          }
          
        } else {
          arma::vec Gri(Gr.row(i).t());
          arma::vec ymf(ym.elem(fri)/yFMiMa(n1 + i, 1)); // y of friends but scaled
          arma::vec tpy(Gri.elem(fri) % pow(ymf, rho));
          double stpy(sum(tpy));
          Gy(n1 + i)  = yFMiMa(n1 + i, 1)*pow(stpy, 1.0 / rho);
          if (deriv) {
            dGy(n1 + i) = Gy(n1 + i)*(sum(tpy % log(ymf))/(rho * stpy) - log(stpy)/pow(rho, 2));// Derivative of Gy for gradient computation
          }
          
          arma::vec zmf(zm.elem(fri)/zFMiMa(n1 + i, 1)); // z of friends but scaled
          arma::vec tpz(Gri.elem(fri) % pow(zmf, rho));
          double stpz(sum(tpz));
          Gz(n1 + i)  = zFMiMa(n1 + i, 1)*pow(stpz, 1.0 / rho);
          if (frzeroz[n1 + i] == 0) {
            double tp(sum(tpz % log(zmf)));
            dGz(n1 + i)  = Gz(n1 + i)*(tp/(rho * stpz) - log(stpz)/pow(rho, 2));
            // if ((n1 + i) == 7081) {
            //   cout<< frzeroz[n1 + i] <<endl;
            //   cout<< Gz(n1 + i) <<endl;
            //   cout<< zm.elem(fri) <<endl;
            //   cout<< stpz <<endl;
            // }
            if (deriv) {
              ddGz(n1 + i) = Gz(n1 + i)*(2*log(stpz)/pow(rho, 3) - 2*tp/(stpz * pow(rho, 2)) + 
                (sum(tpz % pow(log(zmf), 2))*stpz - pow(tp, 2))/(rho * pow(stpz, 2)) + 
                pow(tp/(rho * stpz) - log(stpz)/pow(rho, 2), 2));
            }
          }
        }
      }
    }
  }
  
  arma::mat data(n, Kx + 6, fill::zeros);
  data.cols(0, Kx - 1) = X;
  data.col(Kx)         = y;
  data.col(Kx + 1)     = Gy;
  data.col(Kx + 2)     = Gz;
  data.col(Kx + 3)     = dGz;
  data.col(Kx + 4)     = dGy;
  data.col(Kx + 5)     = ddGz;
  if (FEnum == 1) {
    for (int r(0); r < ngroup; ++ r) {
      arma::uvec Isr = lIs[r], nIsr = lnIs[r];
      if (!Isr.is_empty() || !nIsr.is_empty()) {
        arma::uvec IsnIsr = arma::join_cols(Isr, nIsr);
        arma::mat tp(data.rows(IsnIsr));
        data.rows(IsnIsr) = tp.each_row() - arma::mean(tp, 0);
      }
    }
  } else if (FEnum == 2) {
    for (int r(0); r < ngroup; ++ r) {
      arma::uvec Isr = lIs[r], nIsr = lnIs[r];
      // For isolated
      if (!Isr.is_empty()) {
        arma::mat tp(data.rows(Isr));
        data.rows(Isr) = tp.each_row() - arma::mean(tp, 0);
      }
      // For non-isolated
      if (!nIsr.is_empty()) {
        arma::mat tp(data.rows(nIsr));
        data.rows(nIsr) = tp.each_row() - arma::mean(tp, 0);
      }
    }
  }
  return data;
}


///////////////////////////////////////////////////////////
////////////////////////// fOLS ///////////////////////////
///////////////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::VectorXd fOLS(const Eigen::MatrixXd& data,
                     const Eigen::ArrayXi& idX1,
                     const arma::uvec& Is, 
                     const int& Kx) {
  Eigen::MatrixXd X(data(Is, idX1));
  return (X.transpose()*X).colPivHouseholderQr().solve(X.transpose()*data(Is, Kx));
}

///////////////////////////////////////////////////////////
//////////////////// Compute Moments //////////////////////
///////////////////////////////////////////////////////////
Eigen::ArrayXXd fCESmoments(const Eigen::VectorXd& theta, 
                            const Eigen::ArrayXXd& data,
                            const Eigen::ArrayXi& sel, //row to be selected
                            const arma::uvec& nIs,
                            const Eigen::ArrayXi& idX1,
                            const Eigen::ArrayXi& idX2,
                            const bool& structural,
                            const int& Kx,
                            const int& Kx2,
                            const int& nniso,
                            const int& nst) {
  //
  // Computes the moments of isolated and non-isolated individuals
  //
  Eigen::VectorXd beta = theta.tail(Kx);
  double lambda, lambda2;
  if (structural) {
    lambda2 = theta(1); 
    lambda  = theta(2);
  } else {
    lambda2 = 0; 
    lambda  = theta(1);
  }
  
  // this is how data are organized in dataiso and dataniso
  // X: 0 to Kx - 1
  // y: Kx
  // Gy: Kx + 1
  // Gz: Kx + 2
  // dGz: Kx + 3
  
  // epsilon
  Eigen::ArrayXXd X(data(Eigen::all, Eigen::seq(0, Kx - 1)));
  Eigen::ArrayXd eps;
  Eigen::ArrayXXd Z;
  if (structural) {
    Eigen::ArrayXd Xb1((X(Eigen::all, idX1).matrix()*beta(idX1)).array());
    Z.resize(nniso, Kx2 + 3);
    if (Kx2 > 0) {
      eps = data(nIs, Kx) - Xb1(nIs)*(1 - lambda2) - data(nIs, Kx + 1)*lambda - (X(nIs, idX2).matrix()*beta(idX2)).array();
      Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3), X(nIs, idX2);
    } else {
      eps = data(nIs, Kx) - Xb1(nIs)*(1 - lambda2) - data(nIs, Kx + 1)*lambda;
      Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3);
    }
  } else {
    Z.resize(nst, Kx + 2);
    Z << data(sel, Kx + 2), data(sel, Kx + 3), X(sel, Eigen::all);
    eps = data(sel, Kx) - data(sel, Kx + 1)*lambda - (X(sel, Eigen::all).matrix()*beta).array();
  }
  
  // Moments
  Eigen::ArrayXXd mom(Z.colwise()*eps);
  return mom;
}

////////////////////////////////////////////////////////////////////////////////
//////////////////// GMM Objective Function Given Moments //////////////////////
////////////////////////////////////////////////////////////////////////////////
double fCESgmm(const Eigen::ArrayXXd& mom,
               const Eigen::MatrixXd& W) {
  //
  // GMM objective function following Arelano and Meghir (1992)
  //
  
  // Compute means for non-isolated individuals
  Eigen::VectorXd mmom(mom.colwise().mean().matrix());
  
  // Compute objective function
  return (mmom.transpose() * W * mmom)(0, 0);
}

///////////////////////////////////////////////////////////////////////////////////////
//////////////////// Full Parameters Given non concentrated ones //////////////////////
//////// This main GMM function in Boucher et al. 2024 ECMA, when e == 0 //////////////
///////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fCESgmmparms(const double& rho, 
                        const Eigen::VectorXd& beta1,
                        const Eigen::ArrayXi& idX1,
                        const Eigen::ArrayXi& idX2,
                        const arma::mat& X,
                        const arma::vec& y,
                        const arma::vec& z,
                        const Rcpp::List& G,
                        const Rcpp::List& friendindex,
                        const arma::umat& igroup,
                        const arma::uvec& frzeroy,
                        const arma::uvec& frzeroz,
                        const Rcpp::List& lIs,//in selection
                        const Rcpp::List& lnIs,//in selection
                        const arma::uvec& nIs,
                        const arma::uvec& Is,
                        const Eigen::ArrayXi sel, //row to be selected
                        const Eigen::MatrixXd& W,
                        const IntegerVector& nvec,
                        const Eigen::MatrixXd& yFMiMa,
                        const Eigen::MatrixXd& zFMiMa,
                        const int& n,
                        const int& nniso,
                        const int& niso,
                        const int& nst,
                        const int& Kx,
                        const int& Kx2,
                        const int& ngroup,
                        const int& FEnum,
                        const int& Kest1,
                        const int& Kest2,
                        const int& Kest,
                        const bool& structural,
                        const int& HAC = 0, 
                        const bool& COV = true) {
  int Kx1(Kx - Kx2);
  arma::mat tp(fCESdata(X, y, z, G, friendindex, igroup, frzeroy, frzeroz, lIs, lnIs, nvec, yFMiMa, zFMiMa,
                        n, Kx, ngroup, rho, FEnum, COV));
  Eigen::Map<Eigen::MatrixXd> data(tp.memptr(), n, tp.n_cols);
  // this is how data are organized in data
  // X: 0 to Kx - 1
  // y: Kx
  // Gy: Kx + 1
  // Gz: Kx + 2
  // dGz: Kx + 3
  // dGy: Kx + 4
  // ddGz: Kx + 5
  
  // Extract variables
  Eigen::MatrixXd Xmat(data(Eigen::all, Eigen::seq(0, Kx - 1)));
  Eigen::VectorXd theta;
  Eigen::MatrixXd Vpa;
  double s2(R_NaN), s21(R_NaN), s22(R_NaN);
  
  if (structural) {
    theta.resize(Kx + 3);
    theta(3 + idX1) = beta1;
    theta(0)         = rho;
    
    Eigen::MatrixXd X1(Xmat(Eigen::all, idX1).matrix());
    Eigen::MatrixXd XX1(X1(Is, Eigen::all).transpose()*X1(Is, Eigen::all));
    Eigen::VectorXd Xb1(X1*beta1);
    Eigen::MatrixXd Z(nniso, Kx2 + 3), V(nniso, Kx2 + 2);
    if (Kx2 > 0) {
      Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3), Xmat(nIs, idX2);
      V << Xb1(nIs), data(nIs, Kx + 1), Xmat(nIs, idX2);
    } else {
      Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3);
      V << Xb1(nIs), data(nIs, Kx + 1);
    }
    Eigen::VectorXd y(data(nIs, Kx));
    
    Eigen::MatrixXd ZZ(Z.transpose()*Z);
    Eigen::MatrixXd ZV(Z.transpose()*V);
    Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
    Eigen::VectorXd Zy(Z.transpose()*y);
    Eigen::VectorXd tp(VZWZV.colPivHouseholderQr().solve(VZW*Zy));
    
    // fill theta other component of theta
    theta(1) = 1 - tp(0); // because the variable coefficient is 1 - lambda2
    theta(2) = tp(1);
    if (Kx2 > 0) {
      theta(3 + idX2) = tp.tail(Kx2);
    }
    
    // Variance
    if (COV) {
      Eigen::VectorXd y1hat(Xb1(Is)), y2hat(V*tp);
      Eigen::ArrayXd e1(data(Is, Kx) - y1hat), e2(y - y2hat);
      
      Eigen::MatrixXd dV(nniso, Kx2 + 3), dtp(Eigen::MatrixXd::Zero(Kx2 + 3, Kx2 + 3));
      dV << data(nIs, Kx + 4)*tp(1), V;
      dtp(0, 1) = (e2.matrix().transpose()*data(nIs, Kx + 3))(0, 0);
      dtp(0, 2) = (e2.matrix().transpose()*data(nIs, Kx + 5))(0, 0);
      Eigen::MatrixXd deZ(dV.transpose()*Z - dtp);
      
      Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kx + 3, Kx + 3));
      H.block(0, 0, Kx1, Kx1)  = Eigen::MatrixXd::Identity(Kx1, Kx1);
      H.block(Kx1, Kx1, Kx2 + 3, Kx2 + 3) = deZ*W;
      
      Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(Kx + 3, Kx + 3));
      dF.block(0, 0, Kx1, Kx1) = XX1;
      dF(Eigen::seq(Kx1, Kx + 2), Eigen::all) << Z.transpose()*X1(nIs, Eigen::all)*tp(0),
                                                 Z.transpose()*data(nIs, Kx + 4)*tp(1), ZV;
      
      Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(Kx + 3, Kx + 3));
      if (HAC == 0) {
        s21   = e1.square().sum()/(niso - Kest1);
        s22   = (e2/tp(0)).square().sum()/(nniso - Kest2);
        VF.block(0, 0, Kx1, Kx1) = s21*XX1;
        VF.block(Kx1, Kx1, Kx2 + 3, Kx2 + 3) = s22*ZZ*pow(tp(0), 2);
      }
      if (HAC == 1) {
        Eigen::MatrixXd Xe1(X1(nIs, Eigen::all).array().colwise()*e1);
        Eigen::MatrixXd Ze2(Z.array().colwise()*e2);
        VF.block(0, 0, Kx1, Kx1) = Xe1.transpose()*Xe1;
        VF.block(Kx1, Kx1, Kx2 + 3, Kx2 + 3) = Ze2.transpose()*Ze2;
      }
      if (HAC == 2) {
        Eigen::MatrixXd Xall1(Eigen::MatrixXd::Zero(n, Kx1)), Zall2(Eigen::MatrixXd::Zero(n, Kx2 + 3));
        Eigen::VectorXd e(Eigen::VectorXd::Zero(n));
        
        Xall1(Is, Eigen::all) = X1(Is, Eigen::all);
        Zall2(nIs, Eigen::all) = Z;
        e(Is) = e1;
        e(nIs) = e2;
        
        for (int r(0); r < ngroup; ++ r) {
          int n1(igroup(r, 0)), n2(igroup(r, 1));
          Eigen::MatrixXd tp1(n2 - n1 + 1, Kx + 3);
          tp1 << X1(Eigen::seq(n1, n2), Eigen::all), Z(Eigen::seq(n1, n2), Eigen::all);
          Eigen::VectorXd tp2(tp1.transpose()*e.matrix().segment(n1, n2));
          VF += tp2*tp2.transpose();
        }
      }
      Eigen::MatrixXd iHdF((H*dF).inverse()); 
      Eigen::MatrixXd HVFH(H*VF*H.transpose());
      Vpa = iHdF * HVFH * iHdF.transpose();
    }
  } else {
    theta.resize(Kx + 2);
    theta(0)   = rho;
    
    Eigen::MatrixXd Z(nst, Kx + 2), V(nst, Kx + 1);
    Z << data(sel, Kx + 2), data(sel, Kx + 3), Xmat(sel, Eigen::all);
    V << data(sel, Kx + 1), Xmat(sel, Eigen::all);
    Eigen::VectorXd y(data(sel, Kx));
    
    Eigen::MatrixXd ZV(Z.transpose()*V);
    Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
    Eigen::VectorXd Zy(Z.transpose()*y);
    Eigen::VectorXd tp(VZWZV.colPivHouseholderQr().solve(VZW*Zy));
    
    // fill theta other component of theta
    theta.tail(Kx + 1) = tp;
    
    // Variance
    if (COV) {
      Eigen::VectorXd yhat(V*tp);
      Eigen::ArrayXd e(y - yhat);
      
      Eigen::MatrixXd dV(nst, Kx + 2), dtp(Eigen::MatrixXd::Zero(Kx + 2, Kx + 2));
      dV << data(sel, Kx + 4)*tp(0), V;
      dtp(0, 0) = (e.matrix().transpose()*data(sel, Kx + 3))(0, 0);
      dtp(0, 1) = (e.matrix().transpose()*data(sel, Kx + 5))(0, 0);
      Eigen::MatrixXd deZ(-dV.transpose()*Z + dtp);
      
      Eigen::MatrixXd H(deZ*W);
      Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(Kx + 2, Kx + 2));
      dF << Z.transpose()*data(sel, Kx + 4)*tp(0), ZV;
      
      Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(Kx + 2, Kx + 2));
      if (HAC == 0) {
        s2     = e.square().sum()/(nst - Kest);
        VF     = s2*(Z.transpose()*Z);
      }
      if (HAC == 1) {
        Eigen::MatrixXd Ze(Z.array().colwise()*e);
        VF     = Ze.transpose()*Ze;
      }
      if (HAC == 2) {
        Eigen::MatrixXd Zall(Eigen::MatrixXd::Zero(n, Kx + 2));
        Zall(sel, Eigen::all) = Z;
        Eigen::VectorXd eall(Eigen::VectorXd::Zero(n));
        eall(sel) = e;
        for (int r(0); r < ngroup; ++ r) {
          int n1(igroup(r, 0)), n2(igroup(r, 1));
          Eigen::VectorXd tp1(Zall(Eigen::seq(n1, n2), Eigen::all).transpose()*eall(Eigen::seq(n1, n2)));
          VF += tp1*tp1.transpose();
        }
      }
      Eigen::MatrixXd iHdF((H*dF).inverse()); 
      Eigen::MatrixXd HVFH(H*VF*H.transpose());
      // cout<<H<<endl;
      // cout<<"*****"<<endl;
      // cout<<dF<<endl;
      // cout<<"*****"<<endl;
      // cout<<VF<<endl;
      Vpa = iHdF * HVFH * iHdF.transpose();
    }
  }
  // cout<<theta.transpose()<<endl;
  return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = Vpa, 
                            _["sigma21"] = s21, _["sigma22"] = s22, _["sigma2"] = s2);
}

///////////////////////////////////////////////////////////////////////////////////
//////////////////// GMM Objective Function Given Parameters //////////////////////
//////// This main GMM function in Boucher et al. 2024 ECMA, when e == 1 //////////
///////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double fCESgmmobj(const double& rho, 
                  const Eigen::VectorXd& beta1,
                  const Eigen::ArrayXi& idX1,
                  const Eigen::ArrayXi& idX2,
                  const arma::mat& X,
                  const arma::vec& y,
                  const arma::vec& z,
                  const Rcpp::List& G,
                  const Rcpp::List& friendindex,
                  const arma::umat& igroup,
                  const arma::uvec& frzeroy,
                  const arma::uvec& frzeroz,
                  const Rcpp::List& lIs,//in selection
                  const Rcpp::List& lnIs,//in selection
                  const arma::uvec& nIs,
                  const Eigen::ArrayXi sel, //row to be selected
                  const Eigen::MatrixXd& W,
                  const IntegerVector& nvec,
                  const Eigen::MatrixXd& yFMiMa,
                  const Eigen::MatrixXd& zFMiMa,
                  const int& n,
                  const int& nniso,
                  const int& nst,
                  const int& Kx,
                  const int& Kx2,
                  const int& ngroup,
                  const int& FEnum,
                  const bool& structural) {
  arma::mat tp(fCESdata(X, y, z, G, friendindex, igroup, frzeroy, frzeroz, lIs, lnIs, nvec, yFMiMa, zFMiMa,
                        n, Kx, ngroup, rho, FEnum, false));
  Eigen::Map<Eigen::MatrixXd> data(tp.memptr(), n, tp.n_cols);
  // cout<<rho<<endl;
  // cout<<tp.col(3).t()<<endl;
  // cout<<"-----"<<endl;
  // this is how data are organized in data
  // X: 0 to Kx - 1
  // y: Kx
  // Gy: Kx + 1
  // Gz: Kx + 2
  // dGz: Kx + 3
  
  // Extract variables
  Eigen::MatrixXd Xmat(data(Eigen::all, Eigen::seq(0, Kx - 1)));
  Eigen::VectorXd theta;
  if (structural) {
    theta.resize(Kx + 3);
    theta(3 + idX1) = beta1;
    theta(0)         = rho;
    
    Eigen::VectorXd Xb1(Xmat(Eigen::all, idX1).matrix()*beta1);
    Eigen::MatrixXd Z(nniso, Kx2 + 3), V(nniso, Kx2 + 2);
    if (Kx2 > 0) {
      Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3), Xmat(nIs, idX2);
      V << Xb1(nIs), data(nIs, Kx + 1), Xmat(nIs, idX2);
    } else {
      Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3);
      V << Xb1(nIs), data(nIs, Kx + 1);
    }
    Eigen::VectorXd y(data(nIs, Kx));
    
    Eigen::MatrixXd ZV(Z.transpose()*V);
    Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
    Eigen::VectorXd Zy(Z.transpose()*y);
    Eigen::VectorXd tp(VZWZV.colPivHouseholderQr().solve(VZW*Zy));
    
    // fill theta other component of theta
    theta(1) = 1 - tp(0); // because the variable coefficient is 1 - lambda2
    theta(2) = tp(1);
    if (Kx2 > 0) {
      theta(3 + idX2) = tp.tail(Kx2);
    }
  } else {
    theta.resize(Kx + 2);
    theta(0) = rho;
    
    Eigen::MatrixXd Z(nst, Kx + 2), V(nst, Kx + 1);
    Z << data(sel, Kx + 2), data(sel, Kx + 3), Xmat(sel, Eigen::all);
    V << data(sel, Kx + 1), Xmat(sel, Eigen::all);
    Eigen::VectorXd y(data(sel, Kx));
    Eigen::MatrixXd ZV(Z.transpose()*V);
    Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
    Eigen::VectorXd Zy(Z.transpose()*y);
    Eigen::VectorXd tp(VZWZV.colPivHouseholderQr().solve(VZW*Zy));
    // fill theta other component of theta
    theta.tail(Kx + 1) = tp;
  }
  // objective
  return fCESgmm(fCESmoments(theta, data, sel, nIs, idX1, idX2, structural, Kx, Kx2, nniso, nst), W);
}



///////////////////////////////////////////////////////
//////////////////// GMM Weights //////////////////////
///////////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::MatrixXd fCESWeight(const Eigen::VectorXd& theta, 
                           const Eigen::ArrayXXd& data,
                           const Eigen::ArrayXi& sel, //row to be selected
                           const arma::uvec& nIs,
                           const Eigen::ArrayXi& idX1,
                           const Eigen::ArrayXi& idX2,
                           const bool& structural,
                           const int& Kx,
                           const int& Kx2,
                           const int& nniso,
                           const int& nst) {
  // moments
  Eigen::MatrixXd mom(fCESmoments(theta, data, sel, nIs, idX1, idX2, structural, Kx, Kx2, nniso, nst));
  return ((mom.transpose() * mom) / nst).inverse();
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// GMM when rho is fixed /////////////////////////////
///////////////// In this case dataiso and dataniso are given //////////////////
///////////// This function will be used to find best initial rho //////////////
///// This corresponds to cobj in Boucher et al. 2024 ECMA, when eloc == 0 /////
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Rcpp::List fCESgmmrhoparms(const double& rho,
                           const Eigen::VectorXd& beta1,
                           const Eigen::MatrixXd& data,
                           const Eigen::ArrayXi sel, //row to be selected
                           const arma::uvec& nIs,
                           const arma::uvec& Is,
                           const Eigen::ArrayXi& idX1,
                           const Eigen::ArrayXi& idX2,
                           const Eigen::MatrixXi& igroup,
                           const int& ngroup,
                           const int& Kx,
                           const int& Kx2,
                           const int& nniso,
                           const int& niso,
                           const int& n,
                           const int& nst,
                           const int& Kest1,
                           const int& Kest2,
                           const int& Kest,
                           const bool& structural,
                           const int& rhoinf = 0,
                           const int& HAC = 0, 
                           const bool& COV = true) {
  int Kx1(Kx - Kx2);
  Eigen::MatrixXd Xmat(data(Eigen::all, Eigen::seq(0, Kx - 1)));
  Eigen::VectorXd theta;
  Eigen::MatrixXd Vpa;
  double s2(R_NaN), s21(R_NaN), s22(R_NaN);
  
  if (structural) {
    theta.resize(Kx + 3);
    theta(3 + idX1) = beta1;
    theta(0)        = rho;
    
    Eigen::MatrixXd X1(Xmat(Eigen::all, idX1).matrix());
    Eigen::MatrixXd XX1(X1(Is, Eigen::all).transpose()*X1(Is, Eigen::all));
    Eigen::VectorXd Xb1(X1*beta1);
    Eigen::MatrixXd Z(nniso, Kx2 + 3 - rhoinf), V(nniso, Kx2 + 2);
    if (Kx2 > 0) {
      if (rhoinf == 0) {
        Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3), Xmat(nIs, idX2);
      } else {
        Z << Xb1(nIs), data(nIs, Kx + 2), Xmat(nIs, idX2);
      }
      V << Xb1(nIs), data(nIs, Kx + 1), Xmat(nIs, idX2);
    } else {
      if (rhoinf == 0) {
        Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3);
      } else {
        Z << Xb1(nIs), data(nIs, Kx + 2);
      }
      V << Xb1(nIs), data(nIs, Kx + 1);
    }
    Eigen::VectorXd y(data(nIs, Kx));
    Eigen::MatrixXd ZZ(Z.transpose()*Z);
    Eigen::MatrixXd W((ZZ/nniso).inverse());
    // W.setIdentity();
    
    Eigen::MatrixXd ZV(Z.transpose()*V);
    Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
    Eigen::VectorXd Zy(Z.transpose()*y);
    Eigen::VectorXd tp(VZWZV.colPivHouseholderQr().solve(VZW*Zy));
    
    // fill theta other component of theta
    theta(1) = 1 - tp(0); // because the variable coefficient is 1 - lambda2
    theta(2) = tp(1);
    if (Kx2 > 0) {
      theta(3 + idX2) = tp.tail(Kx2); //should be divided by 1 - theta(1) later
    }
    
    // Variance
    if (COV) {
      Eigen::VectorXd y1hat(Xb1(Is)), y2hat(V*tp);
      Eigen::ArrayXd e1((data(Is, Kx) - y1hat).array()), e2((y - y2hat).array());
      Eigen::MatrixXd H(Eigen::MatrixXd::Zero(Kx + 2, Kx + 3 - rhoinf));
      H.block(0, 0, Kx1, Kx1)  = Eigen::MatrixXd::Identity(Kx1, Kx1);
      H.block(Kx1, Kx1, Kx2 + 2, Kx2 + 3 - rhoinf) = VZW;
      
      Eigen::MatrixXd dF(Eigen::MatrixXd::Zero(Kx + 3 - rhoinf, Kx + 2));
      dF.block(0, 0, Kx1, Kx1) = XX1;
      dF(Eigen::seq(Kx1, Kx + 2), Eigen::all) << (Z.transpose()*X1(nIs, Eigen::all)*tp(0)), ZV;
      
      Eigen::MatrixXd VF(Eigen::MatrixXd::Zero(Kx + 3 - rhoinf, Kx + 3 - rhoinf));
      if (HAC == 0) {
        s21   = e1.square().sum()/(niso - Kest1);
        s22   = (e2/tp(0)).square().sum()/(nniso - Kest2);
        VF.block(0, 0, Kx1, Kx1) = s21*XX1;
        VF.block(Kx1, Kx1, Kx2 + 3 - rhoinf, Kx2 + 3 - rhoinf) = s22*ZZ*pow(tp(0), 2);
      }
      if (HAC == 1) {
        Eigen::MatrixXd Xe1(X1(nIs, Eigen::all).array().colwise()*e1);
        Eigen::MatrixXd Ze2(Z.array().colwise()*e2);
        VF.block(0, 0, Kx1, Kx1) = Xe1.transpose()*Xe1;
        VF.block(Kx1, Kx1, Kx2 + 3 - rhoinf, Kx2 + 3 - rhoinf) = Ze2.transpose()*Ze2;
      }
      if (HAC == 2) {
        Eigen::MatrixXd Xall1(Eigen::MatrixXd::Zero(n, Kx1)), Zall2(Eigen::MatrixXd::Zero(n, Kx2 + 3));
        Eigen::VectorXd e(Eigen::VectorXd::Zero(n));
        
        Xall1(Is, Eigen::all) = X1(Is, Eigen::all);
        Zall2(nIs, Eigen::all) = Z;
        e(Is) = e1;
        e(nIs) = e2;
        
        for (int r(0); r < ngroup; ++ r) {
          int n1(igroup(r, 0)), n2(igroup(r, 1));
          Eigen::MatrixXd tp1(n2 - n1 + 1, Kx + 3 - rhoinf);
          tp1 << X1(Eigen::seq(n1, n2), Eigen::all), Z(Eigen::seq(n1, n2), Eigen::all);
          Eigen::VectorXd tp2(tp1.transpose()*e.matrix().segment(n1, n2));
          VF += tp2*tp2.transpose();
        }
      }
      Eigen::MatrixXd iHdF((H*dF).inverse()); 
      Eigen::MatrixXd HVFH(H*VF*H.transpose());
      Vpa = iHdF * HVFH * iHdF.transpose();
    }
  } else {
    theta.resize(Kx + 2);
    theta(0) = rho;
    
    Eigen::MatrixXd Z(nst, Kx + 2 - rhoinf), V(nst, Kx + 1);
    if (rhoinf == 0) {
      Z << data(sel, Kx + 2), data(sel, Kx + 3), Xmat(sel, Eigen::all);
    } else {
      Z << data(sel, Kx + 2), Xmat(sel, Eigen::all);
    }
    V << data(sel, Kx + 1), Xmat(sel, Eigen::all);
    Eigen::MatrixXd ZZ(Z.transpose()*Z);
    Eigen::VectorXd y(data(sel, Kx));
    Eigen::MatrixXd W((ZZ/nst).inverse());
    // W.setIdentity();
    
    Eigen::MatrixXd ZV(Z.transpose()*V);
    Eigen::MatrixXd VZW(ZV.transpose()*W), VZWZV(VZW*ZV);
    Eigen::VectorXd Zy(Z.transpose()*y);
    Eigen::VectorXd tp(VZWZV.colPivHouseholderQr().solve(VZW*Zy));
    // cout<<Z.sum()<<endl;
    // cout<<V.sum()<<endl;
    // cout<<y.sum()<<endl;
    // cout<<ZZ.sum()<<endl;
    // cout<<W.sum()<<endl;
    // fill theta other component of theta
    theta.tail(Kx + 1) = tp;
    
    // Variance
    if (COV) {
      Eigen::VectorXd yhat(V*tp);
      Eigen::ArrayXd e((y - yhat).array());
      Eigen::MatrixXd VZe(Eigen::MatrixXd::Zero(Kx + 2 - rhoinf, Kx + 2 - rhoinf));
      if (HAC == 0) {
        s2     = e.square().sum()/(nst - Kest);
        VZe    = s2*(ZZ);
      }
      if (HAC == 1) {
        Eigen::MatrixXd Ze(Z.array().colwise()*e);
        VZe    = Ze.transpose()*Ze;
      }
      if (HAC == 2) {
        Eigen::MatrixXd Zall(Eigen::MatrixXd::Zero(n, Kx + 2 - rhoinf));
        Zall(sel, Eigen::all) = Z;
        Eigen::VectorXd eall(Eigen::VectorXd::Zero(n));
        eall(sel) = y - yhat;
        for (int r(0); r < ngroup; ++ r) {
          int n1(igroup(r, 0)), n2(igroup(r, 1));
          Eigen::VectorXd tp1(Zall(Eigen::seq(n1, n2), Eigen::all).transpose()*eall(Eigen::seq(n1, n2)));
          VZe += tp1*tp1.transpose();
        }
      }
      Eigen::MatrixXd iHdF((VZWZV).inverse());
      Eigen::MatrixXd HVFH(VZW*VZe*VZW.transpose());
      Vpa = iHdF * HVFH * iHdF.transpose();
    }
  }
  
  return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = Vpa, 
                            _["sigma21"] = s21, _["sigma22"] = s22, _["sigma2"] = s2);
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////// GMM when rho is fixed /////////////////////////////
///////////////// In this case dataiso and dataniso are given //////////////////
///////////// This function will be used to find best initial rho //////////////
///// This corresponds to cobj in Boucher et al. 2024 ECMA, when eloc == 1 /////
////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double fCESgmmrhoobj(const Eigen::VectorXd& theta,
                     const Eigen::MatrixXd& data,
                     const Eigen::ArrayXi sel, //row to be selected
                     const arma::uvec& nIs,
                     const Eigen::ArrayXi& idX1,
                     const Eigen::ArrayXi& idX2,
                     const int& Kx,
                     const int& Kx2,
                     const int& nniso,
                     const int& nst,
                     const bool& structural,
                     const int& rhoinf = 0) {
  Eigen::MatrixXd Xmat(data(Eigen::all, Eigen::seq(0, Kx - 1))), W;
  if (structural) {
    Eigen::VectorXd Xb1(Xmat(Eigen::all, idX1).matrix()*theta(idX1));
    Eigen::MatrixXd Z(nniso, Kx2 + 3 - rhoinf);
    if (Kx2 > 0) {
      if (rhoinf == 0) {
        Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3), Xmat(nIs, idX2);
      } else {
        Z << Xb1(nIs), data(nIs, Kx + 2), Xmat(nIs, idX2);
      }
    } else {
      if (rhoinf == 0) {
        Z << Xb1(nIs), data(nIs, Kx + 2), data(nIs, Kx + 3);
      } else {
        Z << Xb1(nIs), data(nIs, Kx + 2);
      }
    }
    W = (Z.transpose()*Z/nniso).inverse();
    // W.setIdentity();
  } else {
    Eigen::MatrixXd Z(nst, Kx + 2 - rhoinf);
    if (rhoinf == 0) {
      Z << data(sel, Kx + 2), data(sel, Kx + 3), Xmat(sel, Eigen::all);
    } else {
      Z << data(sel, Kx + 2), Xmat(sel, Eigen::all);
    }
    W = (Z.transpose()*Z/nst).inverse();
    // W.setIdentity();
  }
  return fCESgmm(fCESmoments(theta, data, sel, nIs, idX1, idX2, structural, Kx, Kx2, nniso, nst), W);
}

// This function return the structural parameters using the GMM estimates
// Transformations were applied to estimate a function of parameter
// Here we return the true parameters
//[[Rcpp::export]]
Rcpp::List fCESParam(const arma::vec& param,
                     const arma::mat& covp,
                     const arma::uvec& idX1,
                     const arma::uvec& idX2,
                     const int& Kx,
                     const int& Kx2, 
                     const int& estimrho,
                     const bool& structural,
                     const bool& COV = true) {
  if (structural) {
    int Kx1(Kx - Kx2);
    // Theta
    arma::vec theta(param);
    if (Kx2 > 0) {
      theta.elem(3 + idX2) /= (1 - theta(1));
    }
    
    // Covariance
    if (COV) {
      int Ke(estimrho + 2);
      arma::uvec idx(Ke + Kx);
      idx.head(Ke)          = arma::linspace<arma::uvec>(Kx1, Kx1 + Ke - 1, Ke);
      idx.elem(Ke + idX1)   = arma::linspace<arma::uvec>(0, Kx1 - 1, Kx1);
      if (Kx2 > 0) {
        idx.elem(Ke + idX2) = arma::linspace<arma::uvec>(Kx1 + Ke, Kx + Ke - 1, Kx2);
      }
      
      arma::mat tp(covp.cols(idx));
      tp        = tp.rows(idx);
      arma::mat covt(arma::zeros<arma::mat>(3 + Kx, 3 + Kx));
      covt.submat(1 - estimrho, 1 - estimrho, 2 + Kx, 2 + Kx) = tp;
      arma::mat R(arma::eye<arma::mat>(3 + Kx, 3 + Kx));
      
      // diagonal elements
      arma::vec Rd(arma::ones<arma::vec>(3 + Kx));
      Rd(1)     = -1;
      if (Kx2 > 0) {
        Rd.elem(3 + idX2) /= (1 - theta(1));
      }
      R.diag()  = Rd;
      
      // second column
      if (Kx2 > 0) {
        arma::vec R0(arma::zeros<arma::vec>(3 + Kx));
        R0(1)    = -1;
        R0.elem(3 + idX2) = -theta.elem(3 + idX2)/pow(1 - theta(1), 2);
        R.col(1) = R0;
      }
      
      // theta and covariance
      covt       = R * covt * R.t();
      return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covt);
    } else {
      return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covp);
    }
  } else{
    // theta
    arma::vec theta(param);
    // Covariance
    if (COV) {
      arma::mat covt(arma::zeros<arma::mat>(2 + Kx, 2 + Kx));
      covt.submat(1 - estimrho, 1 - estimrho, 1 + Kx, 1 + Kx) = covp;
      return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covt);
    } else {
      return Rcpp::List::create(_["theta"] = theta, _["Vpa"] = covp);
    }
  }
}

// Best response function for the CES
// talpha is alpha tilde
arma::vec BRCES(const arma::vec& y, 
                const Rcpp::List& G,
                const Rcpp::List& friendindex,
                const arma::umat& igroup,
                const arma::uvec& frzeroy, 
                const IntegerVector& nvec, 
                const Eigen::MatrixXd& yFMiMa,
                const int& n,
                const int& ngroup,
                const arma::vec& talpha,
                const double& lambda,
                const double& rho) {
  vec Gy(n, fill::zeros);
  for (int r(0); r < ngroup; r++) {
    // Extract data for group r
    int nm(nvec(r)), n1(igroup(r, 0)), n2(igroup(r, 1));
    arma::mat Gr = G[r];
    arma::vec ym(y.subvec(n1, n2));
    Rcpp::List frindexm(friendindex[r]);
    
    for (int i(0); i < nm; i++) {
      uvec fri = frindexm[i];
      if (!fri.is_empty()) {
        if (rho == R_PosInf) {
          Gy[n1 + i] = max(ym.elem(fri));
        } else if (rho == R_NegInf) {
          Gy[n1 + i] = min(ym.elem(fri));
        } else if (rho < 0) {
          arma::vec Gri(Gr.row(i).t());
          if (frzeroy[n1 + i] == 0) {
            arma::vec ymf(ym.elem(fri)/yFMiMa(n1 + i, 0)); // y of friends but scaled
            arma::vec tpy(Gri.elem(fri) % pow(ymf, rho));
            double stpy(sum(tpy));
            Gy(n1 + i)  = yFMiMa(n1 + i, 0)*pow(stpy, 1.0 / rho);
          }
        } else {
          arma::vec Gri(Gr.row(i).t());
          arma::vec ymf(ym.elem(fri)/yFMiMa(n1 + i, 1)); // y of friends but scaled
          arma::vec tpy(Gri.elem(fri) % pow(ymf, rho));
          double stpy(sum(tpy));
          Gy(n1 + i)  = yFMiMa(n1 + i, 1)*pow(stpy, 1.0 / rho);
        }
      }
    }
  }
  return talpha + lambda*Gy;
}


// Nash Equilibrium for the CES case
// y is initial solution
//[[Rcpp::export]]
int fNashECES(arma::vec& y,
              Rcpp::List& G,
              const arma::vec& talpha,
              const double& lambda,
              const double& rho,
              const Rcpp::List& friendindex,
              const arma::umat& igroup,
              const arma::uvec& frzeroy, 
              const IntegerVector& nvec,
              const Eigen::MatrixXd& yFMiMa,
              const int& ngroup,
              const int& n,
              const double& tol = 1e-10,
              const int& maxit  = 500){
  int t(0);
  computeBR: ++t;
  // New y
  arma::vec yst(BRCES(y, G, friendindex, igroup, frzeroy, nvec, yFMiMa, n, ngroup, talpha, lambda, rho));
  
  // check convergence
  double dist     = max(arma::abs((yst - y)/(y + 1e-50)));
  y               = yst;
  if (dist > tol && t < maxit) goto computeBR;
  return t; 
}