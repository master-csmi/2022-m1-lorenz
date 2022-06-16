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

MyMatrix mean(MyMatrix M,int a)
{
    MyMatrix Mean;
    int dim_col=M.cols();
    int dim_row=M.rows();
    
    if (a==0)
    {
        Mean=MyMatrix::Zero(dim_col,1);
        for(int i=0;i<dim_row;i++)
        {
            for(int j=0;j<dim_col;j++)
            {
                Mean(j,0)+=M(i,j);
            }
        }
        Mean=Mean/dim_row;
    }
    else
    {
        Mean=MyMatrix(dim_col,1);
        for(int i=0;i<dim_col;i++)
        {
            for(int j=0;j<dim_row;j++)
            {
                Mean(i,0)+=M(i,j);
            }
        }
    }
    return Mean;
}