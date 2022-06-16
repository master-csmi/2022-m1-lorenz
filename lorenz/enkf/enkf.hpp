#include "enkf_fct.cpp"
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
            auto gen1 = Rand::makeMvNormalGen(x, P);
            M_sigmas = (gen1.generate(urng, N)).transpose();
            
            

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
    MyMatrix get_sigmas()const
    {
        return M_sigmas;
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
    void set_sigmas(MyMatrix sigmas)
    {
        M_sigmas=sigmas;
    }
    void update(MyMatrix z)
    {
        MyMatrix sigmas_h,sigmas_x,sigmas_z,z_mean,P_zz,P_xz,sigma,sigma_2;
        MyMatrix K;
        sigmas_h=MyMatrix::Zero(M_N,M_dim_z);
        for(int i=0;i<M_N;i++)
        {
            sigmas_x=MyMatrix::Zero(M_dim_x,1);
            for(int j=0;j<M_dim_x;j++)
            {
                sigmas_x(j,0)=M_sigmas(i,j);
            }
            sigmas_z=M_hx(sigmas_x);
            for(int j=0;j<M_dim_z;j++)
            {
                sigmas_h(i,j)=sigmas_z(j,0);
            }
        }
        z_mean=MyMatrix::Zero(M_dim_z,1);
        //std::cout << "sigmas_h\n "<<sigmas_h<<std::endl;
        z_mean=mean(sigmas_h,0);
        //std::cout << "z_meam\n "<<z_mean<<std::endl;
        P_zz=MyMatrix::Zero(M_dim_z,M_dim_z);
        sigma=MyMatrix::Zero(M_dim_z,1);
        for(int i=0;i<M_N;i++)
        {
            for(int j=0;j<M_dim_z;j++)
            {
                sigma(j,0)=sigmas_h(i,j);
            }
            sigma=sigma-z_mean;
            P_zz+=sigma*sigma.transpose();
        }
        P_zz=(P_zz/(M_N-1))+M_R;
        //std::cout << "P_zz\n "<<P_zz<<std::endl;
        P_xz=MyMatrix::Zero(M_dim_x,M_dim_z);
        sigma=MyMatrix::Zero(M_dim_z,1);
        sigma_2=MyMatrix::Zero(M_dim_x,1);
        for(int i=0;i<M_N;i++)
        {
            for(int j=0;j<M_dim_z;j++)
            {
                sigma(j,0)=sigmas_h(i,j);
            }
            for(int j=0;j<M_dim_x;j++)
            {
                sigma_2(j,0)=M_sigmas(i,j);
            }
            P_xz+=(sigma_2-M_x)*((sigma-z_mean).transpose());
        }
        P_xz=(P_xz/(M_N-1));
        //std::cout << "P_xz\n "<<P_xz<<std::endl;
        P_zz.inverse();
        M_K=P_xz*(P_zz.inverse());
        //std::cout << "K\n "<<M_K<<std::endl;
        std::mt19937_64 urng3( time(NULL) );
        auto gen3 = Rand::makeMvNormalGen(M_mean_z, M_R);
        MyMatrix e_r = (gen3.generate(urng3, M_N)).transpose();
        //std::cout << "e_r\n "<<e_r<<std::endl;
        MyMatrix e_r_i=MyMatrix::Zero(M_dim_z,1);
        MyMatrix sigma_h_i=MyMatrix::Zero(M_dim_z,1);
        for(int i=0;i<M_N;i++)
        {
            for(int j=0;j<M_dim_z;j++)
            {
                sigma_h_i(j,0)=sigmas_h(i,j);
                e_r_i(j,0)=e_r(i,j);
            }
            MyMatrix sigma_i=M_K*(z+e_r_i-sigma_h_i);
            for(int j=0;j<M_dim_x;j++)
            {
                M_sigmas(i,j)+=sigma_i(j,0);
                
            }
        }
        //std::cout << "M_sigmas\n "<<M_sigmas<<std::endl;
        M_x=mean(M_sigmas,0);
        //std::cout << "M_x\n "<<M_x<<std::endl;
        M_P=M_P-((M_K*P_zz)*M_K.transpose());
        //std::cout << "M_P\n "<<M_P<<std::endl;
        M_z=z;
        M_x_post=M_x;
        M_P_post=M_P;
        std::cout << "M_x_post\n "<<M_x_post<<std::endl;
        //std::cout << "M_P_post\n "<<M_P_post<<std::endl;
    }
    void predict()
    {
       
        MyMatrix sigmas_x;
        for (int i=0;i<M_N;i++)
        {
            sigmas_x=MyMatrix::Zero(M_dim_x,1);
            for(int j=0;j<M_dim_x;j++)
            {
                sigmas_x(j,0)=M_sigmas(i,j);
            }
            
            sigmas_x=M_fx(M_dt,sigmas_x);
            
            for(int j=0;j<M_dim_x;j++)
            {
                M_sigmas(i,j)=sigmas_x(j,0);
            }
        }
        
        std::mt19937_64 urng2( time(NULL) );
        auto gen2 = Rand::makeMvNormalGen(M_mean, M_Q);
        MyMatrix e = (gen2.generate(urng2, M_N)).transpose();
        M_sigmas+=e;
        MyMatrix P_=MyMatrix::Zero(M_dim_x,M_dim_x);
        MyMatrix sigma=MyMatrix::Zero(M_dim_x,1);
        for(int i=0;i<M_N;i++)
        {
            for(int j=0;j<M_dim_x;j++)
            {
               
                sigma(j,0)=M_sigmas(i,j);
            }
            sigma=sigma-M_x;
            P_+=sigma*sigma.transpose();
        }
        M_P=(P_/(M_N-1));
        M_x_prior=M_x;
        M_P_prior=M_P;
    }
};

