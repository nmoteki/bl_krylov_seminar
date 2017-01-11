#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <chrono>
#include <fstream>
#include <tuple>
#include <Eigen/Dense>
#include "bl_krylov_solvers.hpp"
#include "linear_algebra_addon.hpp"
using namespace std;
using namespace Eigen;
using namespace std::chrono;

int main(int argc, char *argv[]){

    double tol= 1e-5;
    int itermax= 100;

    std::srand((unsigned int) time(0));


    int m= atoi(argv[1]); // dimension of matrix A
    int s= atoi(argv[2]); // number of RHS vector
    double cond;
    if(argc == 4) cond= atof(argv[3]); // condition number of A (only for demo3)

    MatrixXcd B= MatrixXcd::Random(m,s);


    // demo 1

/*    {
        auto A= MatrixXcd::Random(m,m);

        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicg ..." << endl;
            MatrixXcd X= bl_bicg(A, B, tol, itermax); // applicable to general matrix A but needs A.adjoint()
        cout << " bl_bicg relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicg computation time= " << duration_cast<milliseconds>((end - start)).count() << " milliseconds" << endl ;
        cout << "bl_bicg computation time per RHS= " << duration_cast<milliseconds>((end - start)).count()/s << " milliseconds" << endl << endl<< endl;
    }
*/


    // demo 2
/*    {
        auto A= MatrixXcd::Random(m,m);

        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicg_rq ..." << endl;
            MatrixXcd X= bl_bicg_rq(A, B, tol, itermax); // applicable to general matrix A but needs A.adjoint()
        cout << " bl_bicg_rq relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicg_rq computation time= " << duration_cast<milliseconds>((end - start)).count() << " milliseconds" << endl ;
        cout << "bl_bicg_rq computation time per RHS= " << duration_cast<milliseconds>((end - start)).count()/s << " milliseconds" << endl << endl<< endl;
    }
*/


    // demo 3
    {
        auto A1= MatrixXcd::Random(m,m);
        auto A2= MatrixXcd::Random(m,m);
        MatrixXcd Q1,Q2,R;
        tie(Q1,R)=qr_reduced(A1); // random unitary matrix Q1
        tie(Q2,R)=qr_reduced(A2); // random unitary matrix Q2
        VectorXd vec(m);
        double d= pow(cond,1.0/(m-1));
        for(int i=0; i<m; ++i) vec(i)=pow(d,-m/2.0+i); // set singular value

        auto D= vec.asDiagonal();
        auto A=Q1*D*Q2; // construct A from SVD

        cout << "condition number of A = " << vec.maxCoeff()/vec.minCoeff() << endl;

        time_point<steady_clock> start= steady_clock::now();
        cout << "solving AX=B using bl_bicg_rq ..." << endl;
            MatrixXcd X= bl_bicg_rq(A, B, tol, itermax); // applicable to general matrix A but needs A.adjoint()
        cout << " bl_bicg_rq relative error: " << (A*X-B).norm()/B.norm() << endl << endl;
        auto end= steady_clock::now();
        cout << "bl_bicg_rq computation time= " << duration_cast<milliseconds>((end - start)).count() << " milliseconds" << endl ;
        cout << "bl_bicg_rq computation time per RHS= " << duration_cast<milliseconds>((end - start)).count()/s << " milliseconds" << endl << endl<< endl;

    }





}
