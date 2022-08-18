
#include <boost/test/unit_test.hpp>
//#include <../../src/cpp/enkf/enkf.hpp>
// #include <../../src/cpp/enkf/enkf_fct.hpp>
// #include <../../src/cpp/enkf/enkf.cpp>
// #include <../../src/cpp/enkf/enkf_fct.cpp>
#include <enkf/enkf_fct.hpp>
#include <enkf/enkf.cpp>
#include <enkf/enkf_fct.cpp>

MyMatrix fx_harmonique(double dt,MyMatrix X,int nbr_exhantillon,double dt_2)
{
    double t=dt;
    MyMatrix p{{2}};
    int dim_x=X.cols();
    MyMatrix X_p;
    X_p=MyMatrix::Zero(dim_x,1);
    MyMatrix K1=f_ocillateur(t,X,p);
    MyMatrix K2=f_ocillateur(t+dt/2, X + 1./2. * K1 * dt,p);
    MyMatrix K3=f_ocillateur(t+dt/2,X + 1./2. * K2 * dt,p);
    MyMatrix K4=f_ocillateur(t+dt, X+ K3 * dt,p);
    X_p=X+(dt/6.* (K1+(2.*K2)+(2.*K3)+K4));
    return X_p;
}
MyMatrix hx_harmonique(MyMatrix x)
{
    MyMatrix X{{x(0,0)}};
    return X;
}

MyMatrix read_sensor_harmonique(double t)
{
    MyMatrix X;
    X(0)=2*cos(2*t);
    return X;
} 


BOOST_AUTO_TEST_SUITE(enkf)

BOOST_AUTO_TEST_CASE(test_enkf_1)
{
    using namespace Eigen;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;
    MyMatrix lorenz_1,lorenz_2;
    double T,dt;
    int N,N2;
    MyMatrix p{{12.},{6.},{12.}};
    MyMatrix X_0{{-10.},{-10.},{25.}}; 
    MyMatrix P{{0,0,0},{0,0,0},{0,0,0}};
    MyMatrix Q{{0.001,0,0},{0,0.001,0},{0,0,0.001}};
    MyMatrix R{{0,0,0},{0,0,0},{0,0,0}};
    T=1;
    N=int(T/0.001);
    
    std::cout << "p:  "<<p<<std::endl;
    std::cout << "X_0:  "<<X_0<<std::endl;
    lorenz_1=RK4(3,p,X_0,N,T,&f_lorenz);
    std::cout << "lorenz1 \n  "<<lorenz_1<<std::endl;
    N2=int(T/0.01);
    dt=0.01;
    lorenz_2=RK4(3,p,X_0,N2,T,&f_lorenz);
    std::cout << "lorenz2 \n  "<<lorenz_2<<std::endl;

    EnsembleKalmanFilter E1(3,3,X_0,P,dt, 20, &hx_model,&fx_2);
    E1.set_Q(Q);
    E1.set_R(R);
    int M=10;
    int index=M;
    double time=0;

    MyMatrix etat,tab_temps;
    etat=MyMatrix::Zero(lorenz_2.rows(),lorenz_2.cols());
    tab_temps=MyMatrix::Zero(lorenz_2.rows(),1);


    
    etat.row(0)=((E1.get_x()).col(0)).transpose();
    tab_temps(0,0)=time;
    int compteur_=1;
    while(time<T)
    {
        MyMatrix z=read_sensor_model(index,lorenz_1);
        E1.predict();
        E1.update(z);
        index+=M;
        time+=dt;
        etat.row(compteur_)=((E1.get_x()).col(0)).transpose();
        tab_temps(compteur_,0)=time;
        compteur_+=1;

    }

    int taille_tab_etat=etat.rows();
    int compteur=0;
    for (int i=0;i<taille_tab_etat;i++)
    {
        if (((lorenz_2.row(i)(0)-etat.row(i)(0))<1e-2) and ((lorenz_2.row(i)(1)-etat.row(i)(1))<1e-2) and ((lorenz_2.row(i)(2)-etat.row(i)(2))<1e2) )
        {
            compteur+=1;
        }
    }
    BOOST_CHECK(compteur == taille_tab_etat);
    
}
BOOST_AUTO_TEST_CASE(test_enkf_2)
{
    using namespace Eigen;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;
    MyMatrix lorenz_1,lorenz_2;
    double T,dt;
    int N,N2;
    MyMatrix p{{12.},{6.},{12.}};
    MyMatrix X_0{{-10.},{-10.},{25.}}; 
    MyMatrix P{{0,0,0},{0,0,0},{0,0,0}};
    MyMatrix Q{{0,0,0},{0,0,0},{0,0,0}};
    MyMatrix R{{0.001,0,0},{0,0.001,0},{0,0,0.001}};
    T=1;
    N=int(T/0.01);
    
    std::cout << "p:  "<<p<<std::endl;
    std::cout << "X_0:  "<<X_0<<std::endl;
    lorenz_1=RK4(3,p,X_0,N,T,&f_lorenz);
    std::cout << "lorenz1 \n  "<<lorenz_1<<std::endl;
    N2=int(T/0.1);
    dt=0.1;
    lorenz_2=RK4(3,p,X_0,N2,T,&f_lorenz);
    std::cout << "lorenz2 \n  "<<lorenz_2<<std::endl;

    EnsembleKalmanFilter E1(3,3,X_0,P,dt, 20, &hx_model,&fx_2);
    E1.set_Q(Q);
    E1.set_R(R);
    int M=10;
    int index=M;
    double time=0;

    MyMatrix etat,tab_temps;
    etat=MyMatrix::Zero(lorenz_2.rows(),lorenz_2.cols());
    tab_temps=MyMatrix::Zero(lorenz_2.rows(),1);


    
    etat.row(0)=((E1.get_x()).col(0)).transpose();
    tab_temps(0,0)=time;
    int compteur_=1;
    while(time<1-dt)
    {
        MyMatrix z=read_sensor_model(index,lorenz_1);
        E1.predict();
        E1.update(z);
        index+=M;
        time+=dt;
        etat.row(compteur_)=((E1.get_x()).col(0)).transpose();
        tab_temps(compteur_,0)=time;
        compteur_+=1;

    }
    
    
    int taille_tab_etat=etat.rows();
    int compteur=0;
    for (int i=0;i<taille_tab_etat;i++)
    {
        if (((lorenz_2.row(i)(0)-etat.row(i)(0))<1e-10) and ((lorenz_2.row(i)(1)-etat.row(i)(1))<1e-10) and ((lorenz_2.row(i)(2)-etat.row(i)(2))<1e-10) )
        {
            compteur+=1;
        }
    }
    
    BOOST_CHECK(compteur == taille_tab_etat);
    
}
BOOST_AUTO_TEST_CASE(test_enkf_3)
{
    using namespace Eigen;
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;
    MyMatrix harmonique_1,harmonique_2;

    int N2;
    MyMatrix p{{2}};
    MyMatrix X_0{{2},{0}}; 
    MyMatrix P{{0,0},{0,0}};
    MyMatrix Q{{0,0},{0,0}};
    MyMatrix R{{0.001}};
    double Pe = 2*M_PI/2;
    double dt= Pe/20;
    double T = 3*Pe;
    int N = int(round(T/dt));
    
    std::cout << "p:  "<<p<<std::endl;
    std::cout << "X_0:  "<<X_0<<std::endl;
    harmonique_1=RK4(2,p,X_0,N,T,&f_ocillateur);
    std::cout << "harmonique \n  "<<harmonique_1<<std::endl;
    
    EnsembleKalmanFilter E2(2,1,X_0,P,dt, 20, &hx_harmonique,&fx_harmonique);
    E2.set_Q(Q);
    E2.set_R(R);
    double time=0;

    MyMatrix etat,tab_temps;
    etat=MyMatrix::Zero(harmonique_1.rows(),1);
    tab_temps=MyMatrix::Zero(harmonique_1.rows(),1);


    
    etat.row(0)=((E2.get_x()).col(0)).transpose();
    tab_temps(0,0)=time;
    int compteur_=1;
    while(time<T-dt)
    {
        MyMatrix z{{2*cos(2*time)}};
        E2.predict();
        E2.update(z);
        time+=dt;
        etat.row(compteur_)=((E2.get_x()).col(0)).transpose();
        tab_temps(compteur_,0)=time;
        compteur_+=1;

    }
    
    
    int taille_tab_etat=etat.rows();
    int compteur=0;
    for (int i=0;i<taille_tab_etat;i++)
    {
        if (((etat.row(i)(0)-harmonique_1.row(i)(0))<1e-10))
        {
            compteur+=1;
        }
    }
    BOOST_CHECK(compteur == taille_tab_etat);
}

BOOST_AUTO_TEST_SUITE_END()