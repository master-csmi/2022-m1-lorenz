#include "enkf.hpp"

EnsembleKalmanFilter::EnsembleKalmanFilter(double dim_x,double dim_z,MyMatrix x,MyMatrix P,double dt, int N, MyMatrix ( *hx)(MyMatrix  x),MyMatrix ( *fx)(double dt,MyMatrix  x))
{
    M_dim_x=dim_x;
    M_dim_z=dim_z;
    M_dt=dt;
    M_N=N;

    M_x=MyMatrix::Zero(M_dim_x,1);
    M_z=MyMatrix::Zero(M_dim_z,1);

    M_x_prior=MyMatrix::Zero(M_dim_x,1);
    M_x_post=MyMatrix::Zero(M_dim_x,1);

    M_P=MyMatrix::Zero(M_dim_x,M_dim_x);
    M_P_prior=MyMatrix::Zero(M_dim_x,M_dim_x);
    M_P_post=MyMatrix::Zero(M_dim_x,M_dim_x);


    M_Q=MyMatrix::Zero(M_dim_x,M_dim_x);
    M_R=MyMatrix::Zero(M_dim_z,M_dim_z);
    M_K=MyMatrix::Zero(M_dim_x,M_dim_z);
            
    M_mean=MyMatrix::Zero(M_dim_x,1);
    M_mean_z=MyMatrix::Zero(M_dim_z,1);

    M_sigmas=MyMatrix::Zero(N,M_dim_x);

    M_x=x;
    M_x_post=x;
    M_x_prior=x;
    M_P=P;
    M_P_post=P;
    M_P_post=P;
    M_hx=hx;
    M_fx=fx;
            
    //Initialise
    std::mt19937_64 urng( time(NULL) );
    //Rand::P8_mt19937_64 urng{ time(NULL)  };
    auto gen1 = Rand::makeMvNormalGen(x, P);
    M_sigmas = (gen1.generate(urng, N)).transpose();
}

// accesseurs
int EnsembleKalmanFilter::get_dim_x() const
{
    return M_dim_x;
}
int EnsembleKalmanFilter::get_dim_z() const
{
    return M_dim_z;
}
int EnsembleKalmanFilter::get_N() const
{
    return M_N;
}
double EnsembleKalmanFilter::get_dt() const
{
    return M_dt;
}
MyMatrix EnsembleKalmanFilter::get_x() const
{
    return M_x;
}
MyMatrix EnsembleKalmanFilter::get_z() const
{
    return M_z;
}
MyMatrix EnsembleKalmanFilter::get_P() const
{
    return M_P;
}
MyMatrix EnsembleKalmanFilter::get_Q() const
{
    return M_Q;
}
MyMatrix EnsembleKalmanFilter::get_R() const
{
    return M_R;
}
MyMatrix EnsembleKalmanFilter::get_sigmas()const
{
    return M_sigmas;
}
// mutateur
void EnsembleKalmanFilter::set_dim_x(int dim_x) 
{
    M_dim_x=dim_x;
}
void EnsembleKalmanFilter::set_dim_z(int dim_z)
{
    M_dim_z=dim_z;
}
void EnsembleKalmanFilter::set_N(int N) 
{
    M_N=N;
}
void EnsembleKalmanFilter::set_dt(double dt)
{
    M_dt=dt;
}
void EnsembleKalmanFilter::set_x(MyMatrix x)
{
    M_x=x;
}
void EnsembleKalmanFilter::set_z(MyMatrix z) 
{
    M_z=z;
}
void EnsembleKalmanFilter::set_P(MyMatrix P) 
{
    M_P=P;
}
void EnsembleKalmanFilter::set_Q(MyMatrix Q) 
{
    M_Q=Q;
}
void EnsembleKalmanFilter::set_R(MyMatrix R) 
{
    M_R=R;
}
void EnsembleKalmanFilter::set_sigmas(MyMatrix sigmas)
{
    M_sigmas=sigmas;
}
void EnsembleKalmanFilter::update(MyMatrix z)
{
        
    MyMatrix sigmas_h=MyMatrix::Zero(M_N,M_dim_z);
    for(int i=0;i<M_N;i++)
    {
        sigmas_h.row(i)=M_hx(M_sigmas.row(i).transpose()).transpose();

    }
    MyMatrix z_mean=MyMatrix::Zero(M_dim_z,1);
    z_mean=mean(sigmas_h,0);



    MyMatrix P_zz=MyMatrix::Zero(M_dim_z,M_dim_z);
    for(int i=0;i<M_N;i++)
    {
        P_zz+=(sigmas_h.row(i).transpose()-z_mean)*(sigmas_h.row(i).transpose()-z_mean).transpose();
    }
    P_zz=(P_zz/(M_N-1))+M_R;




    MyMatrix P_xz=MyMatrix::Zero(M_dim_x,M_dim_z);
    for(int i=0;i<M_N;i++)
    {
        P_xz+=(M_sigmas.row(i).transpose()-M_x)*((sigmas_h.row(i).transpose()-z_mean).transpose());
    }
    P_xz=(P_xz/(M_N-1));
    P_zz.inverse();
    MyMatrix M_K=P_xz*(P_zz.inverse());



    std::mt19937_64 urng3( time(NULL) );
    auto gen3 = Rand::makeMvNormalGen(M_mean_z, M_R);
    MyMatrix e_r = (gen3.generate(urng3, M_N)).transpose();
    for(int i=0;i<M_N;i++)
    {
            
        M_sigmas.row(i)+=(M_K*(z+e_r.row(i).transpose()-sigmas_h.row(i).transpose())).transpose();
    }
    M_x=mean(M_sigmas,0);
    M_P=M_P-((M_K*P_zz)*M_K.transpose());
    M_z=z;
    M_x_post=M_x;
    M_P_post=M_P;
}
void EnsembleKalmanFilter::predict()
{
   
    for (int i=0;i<M_N;i++)
    {
        M_sigmas.row(i)=M_fx(M_dt,M_sigmas.row(i).transpose()).transpose();
    }
    std::mt19937_64 urng2( time(NULL) );
    auto gen2 = Rand::makeMvNormalGen(M_mean, M_Q);
    MyMatrix e = (gen2.generate(urng2, M_N)).transpose();
    M_sigmas+=e;
    MyMatrix P_=MyMatrix::Zero(M_dim_x,M_dim_x);
    for(int i=0;i<M_N;i++)
    {
        P_+=(M_sigmas.row(i).transpose()-M_x)*(M_sigmas.row(i).transpose()-M_x).transpose();
    }
    M_P=(P_/(M_N-1));
    M_x_prior=M_x;
    M_P_prior=M_P;
}
void EnsembleKalmanFilter::creation_csv(MyMatrix lorenz1,double dt_1,MyMatrix lorenz2,double dt_2,MyMatrix etat,double dt_3)
{
    int dim_col=lorenz1.cols();
    int dim_row=lorenz1.rows();
    std :: ofstream ofile("../../../examples/cpp/enkf/plotlorenz1.csv",std :: ios :: out|std :: ios :: app );
    if ( ofile ) 
    {
        double time=0 ;
        for (int i=0;i<dim_row;i++)
        {
            ofile << time <<" ";
                
            for (int j=0;j<dim_col-1;j++)
            {
                ofile <<lorenz1(i,j)<<" ";
             }
            ofile <<lorenz1(i,dim_col-1)<<std :: endl;
           
            time=time+dt_1;   
        }
    }
    ofile.close ();
    int dim_col2=lorenz2.cols();
    int dim_row2=lorenz2.rows();
    std :: ofstream oofile("../../../examples/cpp/enkf/plotlorenz2.csv",std :: ios :: out|std :: ios :: app );
    if ( oofile )
    {
        double time2=0 ;
        for (int i=0;i<dim_row2;i++)
        {
            oofile << time2 <<" ";
                    
            for (int j=0;j<dim_col2-1;j++)
            {
                oofile <<lorenz2(i,j)<<" ";
            }
            oofile <<lorenz2(i,dim_col2-1)<<std :: endl;
           
            time2=time2+dt_2;   
        }
    }
    oofile.close ();
    int dim_col3=etat.cols();
    int dim_row3=etat.rows();
    std :: ofstream ooofile("../../../examples/cpp/enkf/plotetat.csv",std :: ios :: out|std :: ios :: app );
    if ( ooofile ) 
    {
        double time3=0 ;
        for (int i=0;i<dim_row3;i++)
        {
            ooofile << time3 <<" ";
                    
            for (int j=0;j<dim_col3-1;j++)
            {
                ooofile <<etat(i,j)<<" ";
            }
            ooofile <<etat(i,dim_col3-1)<<std :: endl;
           
            time3=time3+dt_3;   
        }
    }
    ooofile.close ();
}

