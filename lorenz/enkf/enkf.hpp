#ifndef READ_FILE_HPP
#define READ_FILE_HPP

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>

using namespace Eigen;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;

class EnsembleKalmanFilter
{
    private:
    int M_dim_x,M_dim_z,M_N;
    double M_dt;
    MyMatrix (* M_hx)(MyMatrix * x);
    MyMatrix (* M_fx)(double dt,MyMatrix * x);
    MyMatrix* M_x,M_K,M_z,M_S,M_SI;
    MyMatrix* M_Q,M_R,M_P,M_mean,M_mean_z;
    MyMatrix* M_x_prior,M_P_prior,M_x_post,M_P_post;
    x, P, dim_z, dt, N, hx, fx
    EnsembleKalmanFilter(MyMatrix *x,MyMatrix *P, int dim_z,double dt, int N, MyMatrix (* hx)(MyMatrix * x),MyMatrix (* fx)(double dt,MyMatrix * x))
    {
        M_dt=dt;
        M_x=x;
        M_P=P;
        M_dim_z=dim_z;
        M_N=N;
        M_hx=hx;
        M_fx=fx;
        

    }




}