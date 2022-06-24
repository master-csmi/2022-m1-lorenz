#include"enkf_fct.hpp"
void essaie()
{
    MyMatrix A=MyMatrix::Zero(1,1);
    std::cout << "A \n  "<<A<<std::endl;
}
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
MyMatrix f_lorenz(double t,MyMatrix X,MyMatrix p)
{
    MyMatrix f;
    f=MyMatrix(3,1);
    f(0)=p(0)*(X(1)-X(0));
    f(1)=X(0)*(p(2)-X(2))-X(1);
    f(2)=X(0)*X(1)-p(1)*X(2);
    return f; 
}
MyMatrix hx_model(MyMatrix x)
{
    return x;
}

MyMatrix read_sensor_model(int index,MyMatrix lorenz_1)
{
    MyMatrix z;
    int dim_z=lorenz_1.cols();
    z=MyMatrix::Zero(dim_z,1);
    z=lorenz_1.row(index).transpose();
    return z;
} 
MyMatrix fx(MyMatrix X,double t,double dt,MyMatrix p)
{
    
    int dim_x=X.cols();
    MyMatrix X_p,K1,K2,K3,K4;
    X_p=MyMatrix(dim_x,1);
    K1=MyMatrix(dim_x,1);
    K2=MyMatrix(dim_x,1);
    K3=MyMatrix(dim_x,1);
    K4=MyMatrix(dim_x,1);
    K1=f_lorenz(t,X,p);
    K2=f_lorenz(t+dt/2, X + 1./2. * K1 * dt,p);
    K3=f_lorenz(t+dt/2,X + 1./2. * K2 * dt,p);
    K4=f_lorenz(t+dt, X+ K3 * dt,p);
    X_p=X+(dt/6.* (K1+(2.*K2)+(2.*K3)+K4));
    return X_p;
}
MyMatrix fx_2(double dt,MyMatrix X)
{
    double t=dt;
    MyMatrix p{{12.},{6.},{12.}};
    int dim_x=X.cols();
    MyMatrix X_p;
    X_p=MyMatrix::Zero(dim_x,1);
    MyMatrix K1=f_lorenz(t,X,p);
    MyMatrix K2=f_lorenz(t+dt/2, X + 1./2. * K1 * dt,p);
    MyMatrix K3=f_lorenz(t+dt/2,X + 1./2. * K2 * dt,p);
    MyMatrix K4=f_lorenz(t+dt, X+ K3 * dt,p);
    X_p=X+(dt/6.* (K1+(2.*K2)+(2.*K3)+K4));
    return X_p;
}

MyMatrix RK4(int dim_x,MyMatrix p,MyMatrix X0,int N,double T)
{
    double dt=T/N;
    MyMatrix X,K1,K2,K3,K4,tab_t,X_n_1,X_n;
    X=MyMatrix(N+1,dim_x);
    K1=MyMatrix(dim_x,1);
    X_n=MyMatrix(dim_x,1);
    K2=MyMatrix(dim_x,1);
    K3=MyMatrix(dim_x,1);
    K4=MyMatrix(dim_x,1);
    tab_t=MyMatrix(N+1,1);
    X.row(0)=X0.transpose();
    tab_t(0)=0;

    double t=0;
    for(int n=1;n<N+1;n++)
    {
        K1=f_lorenz(t,X.row(n-1).transpose() ,p);
        K2=f_lorenz(t+dt/2, X.row(n-1).transpose() + 1./2. * K1 * dt,p);
        K3=f_lorenz(t+dt/2,X.row(n-1).transpose()  + 1./2. * K2 * dt,p);
        K4=f_lorenz(t+dt, X.row(n-1).transpose() + K3 * dt,p);
        tab_t(n)=t+dt;
        X.row(n)=X.row(n-1).transpose() +(dt/6.* (K1+(2.*K2)+(2.*K3)+K4));
        
    }
    return X;
    
}