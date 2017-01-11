#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include "linear_algebra_addon.hpp"
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicg(const MatrixXcd& A, const MatrixXcd& B, const double& tol, const int& itermax)
{
    // Rashedi et al 2016, On short recurrence Krylov type methods for linear systems with many right-hand sides
  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B-A*X;
  MatrixXcd R_hat= R; // or Rhat= R.conjugate();
  MatrixXcd P= R;
  MatrixXcd P_hat= R_hat;


  for(int k= 0; k < itermax; ++k){
      MatrixXcd Q= A*P;
      MatrixXcd Q_hat=A.adjoint()*P_hat;

      MatrixXcd alpha= (P_hat.adjoint()*Q).partialPivLu().solve(R_hat.adjoint()*R);
      MatrixXcd alpha_hat= (P.adjoint()*Q_hat).partialPivLu().solve(R.adjoint()*R_hat);

      X= X+P*alpha;
      MatrixXcd Rnew= R-Q*alpha;
      MatrixXcd Rnew_hat= R_hat-Q_hat*alpha_hat;

      MatrixXcd beta= (R_hat.adjoint()*R).partialPivLu().solve(Rnew_hat.adjoint()*Rnew);
      MatrixXcd beta_hat= (R.adjoint()*R_hat).partialPivLu().solve(Rnew.adjoint()*Rnew_hat);

      double err= R.norm()/Bnorm;
      cout << "bl_bicg: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;

      R= Rnew;
      R_hat= Rnew_hat;

      P= R+P*beta;
      P_hat= R_hat+P_hat*beta_hat;
  }

  if((A*X-B).norm()/Bnorm > 10*tol){
      cerr << "bl_bicg did not converge to solution within error tolerance !" << endl;
     // exit(EXIT_FAILURE);
  }

  return X;
}
