
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
    MyMatrix* M_x,*M_K,*M_z,*M_S,*M_SI;
    MyMatrix* M_Q,*M_R,*M_P,*M_mean,*M_mean_z;
    MyMatrix* M_x_prior,*M_P_prior,*M_x_post,*M_P_post;
    public:
        EnsembleKalmanFilter(double dim_x,double dim_z,MyMatrix *x,MyMatrix *P,MyMatrix *Q,MyMatrix *R,double dt, int N, MyMatrix (* hx)(MyMatrix * x),MyMatrix (* fx)(double dt,MyMatrix * x))
        {
            M_dim_x=dim_x;
            M_dim_z=dim_z;
            M_dt=dt;
            M_N=N;
            M_x=new MyMatrix(M_dim_x,1);
            M_z=new MyMatrix(M_dim_z,1);
            M_P=new MyMatrix(M_dim_x,M_dim_x);
            M_K=new MyMatrix(M_dim_x,M_dim_z);
            M_S=new MyMatrix(M_dim_z,M_dim_z);
            M_SI=new MyMatrix(M_dim_z,M_dim_z);
            M_Q=new MyMatrix(M_dim_x,M_dim_x);
            M_R=new MyMatrix(M_dim_z,M_dim_z);
            M_x=x;
            M_P=P;
            M_Q=Q;
            M_R=R;
            M_hx=hx;
            M_fx=fx;
        }
    int get_dim_x() const
    {
        return M_dim_x;
    }
    int get_dim_z() const
    {
        return M_dim_z;
    }
    int get_N() const
    {
        return M_N;
    }
    double get_dt() const
    {
        return M_dt;
    }
    MyMatrix* get_x() const
    {
        return M_x;
    }
    MyMatrix* get_z() const
    {
        return M_z;
    }
    MyMatrix* get_P() const
    {
        return M_P;
    }
    MyMatrix* get_Q() const
    {
        return M_Q;
    }
    MyMatrix* get_R() const
    {
        return M_R;
    }
        
};

