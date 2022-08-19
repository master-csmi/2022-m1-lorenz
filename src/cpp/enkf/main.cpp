#include "enkf.hpp"



int main()
{
    
    MyMatrix lorenz_1,lorenz_2;
    double T,dt;
    int N,N2;
    MyMatrix p1{{12.},{6.},{12.}};
    MyMatrix p2{{10.},{6.},{10.}};
    MyMatrix X_0{{-10.},{10.},{25.}}; 
    MyMatrix P{{0.1,0,0},{0,0.1,0},{0,0,0.1}};
    MyMatrix Q{{0.1,0,0},{0,0.1,0},{0,0,0.1}};
    MyMatrix R{{0.001,0,0},{0,0.001,0},{0,0,0.001}};
    T=4;
    double t=0;
    N=int(T/0.1);
    
    std::cout << "p_obs:  "<<p1<<std::endl;
    std::cout << "p_mod:  "<<p2<<std::endl;
    std::cout << "X_0:  "<<X_0<<std::endl;
    lorenz_1=RK4(3,p1,X_0,N,T,&f_lorenz);
    std::cout << "lorenz1 \n  "<<lorenz_1<<std::endl;
    N2=int(T/0.1);
    dt=0.1;
    lorenz_2=RK4(3,p2,X_0,N2,T,&f_lorenz);
    std::cout << "lorenz2 \n  "<<lorenz_2<<std::endl;

    EnsembleKalmanFilter E1(3,3,X_0,P,dt, 20, &hx_model,&fx_2);
    E1.set_Q(Q);
    E1.set_R(R);
    int M=1;
    int index=M;
    double time=0;

    MyMatrix etat,tab_temps;
    etat=MyMatrix::Zero(lorenz_2.rows(),lorenz_2.cols());
    tab_temps=MyMatrix::Zero(lorenz_2.rows(),1);

    int dim_col=lorenz_2.cols();
    
    etat.row(0)=((E1.get_x()).col(0)).transpose();
    tab_temps(0,0)=time;
    int compteur=1;
    while(time<T)
    {
        MyMatrix z=read_sensor_model(index,lorenz_1);
        E1.predict();
        E1.update(z);
        index+=M;
        time+=dt;
        etat.row(compteur)=((E1.get_x()).col(0)).transpose();
        tab_temps(compteur,0)=time;
        compteur+=1;

    }
    std::cout << "etat\n  "<<etat<<std::endl;

    E1.creation_csv(lorenz_1,0.1,lorenz_2,0.1,etat,0.1); 
    return 0;
}