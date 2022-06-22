#include "enkf_fct.hpp"
/**
 * @file enkf.hpp
 * @author Melissa AYDOGDU
 * @brief The enkf class.
 * @version 1.0
 * @date 2022-06-22
 */
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
        EnsembleKalmanFilter(double dim_x,double dim_z,MyMatrix x,MyMatrix P,double dt, int N, MyMatrix ( *hx)(MyMatrix  x),MyMatrix ( *fx)(double dt,MyMatrix  x));
    
    
        // accesseurs
        int get_dim_x()const;
        int get_dim_z()const;
        int get_N()const;
        double get_dt()const;
        MyMatrix get_x()const;
        MyMatrix get_z()const;
        MyMatrix get_P()const;
        MyMatrix get_Q()const;
        MyMatrix get_R()const;
        MyMatrix get_sigmas()const;
        
        // mutateur
        void set_dim_x(int dim_x);
        void set_dim_z(int dim_z);
        void set_N(int N) ;
        void set_dt(double dt);
        void set_x(MyMatrix x);
        void set_z(MyMatrix z) ;
        void set_P(MyMatrix P);
        void set_Q(MyMatrix Q);
        void set_R(MyMatrix R);
        void set_sigmas(MyMatrix sigmas);



        void update(MyMatrix z);
        void predict();
        void creation_csv(MyMatrix lorenz1,double dt_1,MyMatrix lorenz2,double dt_2,MyMatrix etat,double dt_3);
 };

