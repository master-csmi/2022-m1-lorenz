#include "enkf_fct.hpp"
/**
 * @file enkf.hpp
 * @author Melissa AYDOGDU
 * @brief The enkf class.
 * @version 1.0
 * @date 2022-06-22
 */

/**
 * @brief An enkf class to implement the ensemble kalman filter method.
 *
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
        /**
         * @brief Get the M_dim_x attibute.
         * 
         * @return type:int The dimensions of the state matrix.
         */
        int get_dim_x()const;

        /**
         * @brief Get the M_dim_z attibute.
         * 
         * @return type:int The dimensions of the observations matrix.
         */
        int get_dim_z()const;

        /**
         * @brief Get the M_N attibute.
         * 
         * @return type:int Number of sigma points (ensembles).
         */
        int get_N()const;

        /**
         * @brief Get the M_dt attibute.
         * 
         * @return type:double Time step in seconds.
         */
        double get_dt()const;

        /**
         * @brief Get the M_x attibute.
         * 
         * @return type:MyMatrix State mean matrix.
         */
        MyMatrix get_x()const;

        /**
         * @brief Get the M_z attibutes.
         * 
         * @return type:MyMatrix Last measurement used in update().
         */
        MyMatrix get_z()const;

        /**
         * @brief Get the M_P attibutes.
         * 
         * @return type:MyMatrix State covariance matrix.
         */
        MyMatrix get_P()const;

        /**
         * @brief Get the M_Q attibutes.
         * 
         * @return type:MyMatrix Process noise matrix.
         */
        MyMatrix get_Q()const;

        /**
         * @brief Get the M_R attibutes.
         * 
         * @return type:MyMatrix  Measurement noise matrix
         */
        MyMatrix get_R()const;

        /**
         * @brief Get the M_sigmas attibutes.
         * 
         * @return MyMatrix 
         */
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


        /**
         * @brief Add a new measurement (z) to the kalman filter
         * 
         * @param z An eigen matrix with the new observations.
         */
        void update(MyMatrix z);

        /**
         * @brief Predict next position.
         * 
         */
        void predict();

        /**
         * @brief Create 3 csv files.
         * 
         * @param lorenz1 An eigen matrix with the model.
         * @param dt_1 An eigen matrix with the discretisation for the model.
         * @param lorenz2 An eigen matrix with the observations.
         * @param dt_2 An eigen matrix with the discretisation for the observations.
         * @param etat An eigen matrix with the analyse state.
         * @param dt_3 An eigen matrix with the discretisation for analyse state.
         */
        void creation_csv(MyMatrix lorenz1,double dt_1,MyMatrix lorenz2,double dt_2,MyMatrix etat,double dt_3);
 };

