#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <sstream>
#include<time.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <EigenRand/EigenRand>

using namespace Eigen;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;

MyMatrix mean(MyMatrix M,int a);
MyMatrix f_lorenz(double t,MyMatrix X,MyMatrix p);
MyMatrix hx_model(MyMatrix x);

MyMatrix read_sensor_model(int index,MyMatrix lorenz_1);
MyMatrix fx(MyMatrix X,double t,double dt,MyMatrix p);
MyMatrix fx_2(double dt,MyMatrix X);

MyMatrix RK4(int dim_x,MyMatrix p,MyMatrix X0,int N,double T);