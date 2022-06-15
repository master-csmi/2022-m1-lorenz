
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
#include <EigenRand/EigenRand>
using namespace Eigen;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;

class EnsembleKalmanFilter
{
    private:
        int M_dim_x,M_dim_z,M_N;
        double M_dt;
        MyMatrix (* M_hx)(MyMatrix  x);
        MyMatrix (* M_fx)(double dt,MyMatrix x);
        MyMatrix M_x,M_K,M_z;
        MyMatrix M_Q,M_R,M_P,M_mean,M_mean_z;
        MyMatrix M_x_prior,M_P_prior,M_x_post,M_P_post;
        MyMatrix M_sigmas;
    public:
        EnsembleKalmanFilter(double dim_x,double dim_z,MyMatrix x,MyMatrix P,double dt, int N, MyMatrix ( *hx)(MyMatrix  x),MyMatrix ( *fx)(double dt,MyMatrix  x))
        {
            M_dim_x=dim_x;
            M_dim_z=dim_z;
            M_dt=dt;
            M_N=N;

            M_x=MyMatrix(M_dim_x,1);
            M_z=MyMatrix(M_dim_z,1);

            M_x_prior=MyMatrix(M_dim_x,1);
            M_x_post=MyMatrix(M_dim_x,1);

            M_P=MyMatrix(M_dim_x,M_dim_x);
            M_P_prior=MyMatrix(M_dim_x,M_dim_x);
            M_P_post=MyMatrix(M_dim_x,M_dim_x);


            M_Q=MyMatrix(M_dim_x,M_dim_x);
            M_R=MyMatrix(M_dim_z,M_dim_z);
            M_K=MyMatrix(M_dim_x,M_dim_z);
            
            M_mean=MyMatrix(M_dim_x,1);
            M_mean_z=MyMatrix(M_dim_z,1);

            //M_sigmas=MyMatrix(N,M_dim_x);


            M_x=x;
            M_x_post=x;
            M_x_prior=x;
            M_P=P;
            M_P_post=P;
            M_P_post=P;
            M_hx=hx;
            M_fx=fx;

            //Initialise
            std::mt19937_64 urng(1);
            auto gen1 = Rand::makeMvNormalGen(x, P);
            M_sigmas = gen1.generate(urng, N);

        }
    // accesseurs
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
    MyMatrix get_x() const
    {
        return M_x;
    }
    MyMatrix get_z() const
    {
        return M_z;
    }
    MyMatrix get_P() const
    {
        return M_P;
    }
    MyMatrix get_Q() const
    {
        return M_Q;
    }
    MyMatrix get_R() const
    {
        return M_R;
    }
    // mutateur
    void set_dim_x(int dim_x) 
    {
        M_dim_x=dim_x;
    }
    void set_dim_z(int dim_z)
    {
        M_dim_z=dim_z;
    }
    void set_N(int N) 
    {
        M_N=N;
    }
    void set_dt(double dt)
    {
        M_dt=dt;
    }
    void set_x(MyMatrix x)
    {
        M_x=x;
    }
    void set_z(MyMatrix z) 
    {
        M_z=z;
    }
    void set_P(MyMatrix P) 
    {
        M_P=P;
    }
    void set_Q(MyMatrix Q) 
    {
        M_Q=Q;
    }
    void set_R(MyMatrix R) 
    {
        M_R=R;
    }
    void update(MyMatrix z)
    {
        /*sigmas_h=MyMatrix(M_N,M_dim_z);
        for(int i=0;i<M_n;i++)
        {
            sigmas_h[i]=M_hx()
        }
        */




    }
    void predict()
    {





    }
};

