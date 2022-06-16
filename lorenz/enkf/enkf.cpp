#include "enkf.hpp"


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
    for(int i=0;i<dim_z;i++)
    {
        
        z(i,0)=lorenz_1(index,i);
    }
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
    X(0,0)=X0(0);
    X(0,1)=X0(1);
    X(0,2)=X0(2);
    tab_t(0)=0;

    double t=0;
    for(int n=1;n<N+1;n++)
    {
        X_n_1=MyMatrix(dim_x,1);
        X_n_1(0,0)=X(n-1,0);
        X_n_1(1,0)=X(n-1,1);
        X_n_1(2,0)=X(n-1,2);
        K1=f_lorenz(t,X_n_1,p);
        K2=f_lorenz(t+dt/2, X_n_1 + 1./2. * K1 * dt,p);
        K3=f_lorenz(t+dt/2,X_n_1 + 1./2. * K2 * dt,p);
        K4=f_lorenz(t+dt, X_n_1+ K3 * dt,p);
        tab_t(n)=t+dt;
        X_n=X_n_1+(dt/6.* (K1+(2.*K2)+(2.*K3)+K4));
        X(n,0)=X_n(0,0);
        X(n,1)=X_n(1,0);
        X(n,2)=X_n(2,0);
    }
    return X;
    
}
int main()
{
    
    MyMatrix lorenz_1,lorenz_2;
    double T,dt;
    int N,N2;
    MyMatrix p{{12.},{6.},{12.}};
    MyMatrix X_0{{-10.},{-10.},{25.}}; 
    MyMatrix P{{0.1,0,0},{0,0.1,0},{0,0,0.1}};
    MyMatrix Q{{0.1,0,0},{0,0.1,0},{0,0,0.1}};
    MyMatrix R{{0.001,0,0},{0,0.001,0},{0,0,0.001}};
    T=1;
    double t=0;
    N=int(T/0.01);
    
    std::cout << "p:  "<<p<<std::endl;
    std::cout << "X_0:  "<<X_0<<std::endl;
    lorenz_1=RK4(3,p,X_0,N,T);
    std::cout << "lorenz1 \n  "<<lorenz_1<<std::endl;
    N2=int(T/0.1);
    dt=0.1;
    lorenz_2=RK4(3,p,X_0,N2,T);
    std::cout << "lorenz2 \n  "<<lorenz_2<<std::endl;
    EnsembleKalmanFilter E1(3,3,X_0,P,dt, 20, &hx_model,&fx_2);
    E1.set_Q(Q);
    E1.set_R(R);
    int M=10;
    int index=M;
    double time=0;
    while(time<1-dt)
    {
        MyMatrix z=read_sensor_model(index,lorenz_1);
        E1.predict();
        E1.update(z);
        index+=10;
        time+=dt;
    }

    

    

    



   /*
   //Rand::P8_mt19937_64 urng{ 1 };
   //Rand::P8_mt19937_64_32  urng{ 1 };

   //Rand::Vmt19937_64 urng{ 1 };
   //Rand::Pmt19937_64 urng{ 1 };
    //std::mt19937_64 urng(1);
    
    Vector4f mean{ 0, 1, 2, 3 };
    Matrix4f cov;
    cov << 1, 1, 0, 0,
           1, 2, 0, 0,
           0, 0, 3, 1,
           0, 0, 1, 2;
    // constructs MvNormalGen with Scalar=float, Dim=4
   
    auto gen1 = Rand::makeMvNormalGen(mean, cov);
    std::mt19937 urng(12345);
    Matrix<float, 4, -1> samples = gen1.generate(urng, 10);
    std::cout << "Sample: \n"<<samples<< std::endl;
    std::mt19937 urng2(12345);
    std::cout << urng2() << std::endl;
    */
    return 0;
}