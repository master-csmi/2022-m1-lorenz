#include <iostream>
#include <cmath>
#include <Eigen/Dense>

template<class T> using Vector = Eigen::Matrix<T,1,Eigen::Dynamic>;

int main(){
    int n_proc=3;
    double t0=0.;
    double T=2.;
    double dt_G=0.1;
    double dt_F=0.01;
    // double dt_G=0.15;
    // double dt_F=0.015;
    // double dt_G=0.15;
    // double dt_F=0.015*5;
    // double dt_G=0.15;
    // double dt_F=0.015*2;

    int nb_t_G = (T-t0)/dt_G;
    std::cout << "nb_t_G : " << nb_t_G << std::endl;
    
    Vector<int> tab_nb_t_G_p(n_proc);
    for(int i=0; i<n_proc; i++){
        tab_nb_t_G_p(i) = nb_t_G/n_proc;
    }
    int reste = nb_t_G-tab_nb_t_G_p.sum();
    int index=0;
    for(int i=0; i<reste; i++){
        if(index==n_proc-1)
            index=0;
        tab_nb_t_G_p(index)++;
        index++;
    }
    std::cout << "tab_nb_t_G_p : " << tab_nb_t_G_p << std::endl;
    
    Vector<double> temps_G(n_proc+1);
    temps_G(0) = t0;
    for(int i=1; i<n_proc+1; i++){
        temps_G(i) = temps_G(i-1) + tab_nb_t_G_p(i-1) * dt_G;
    }
    std::cout << "temps_G : " << temps_G << std::endl;

    Vector<int> tab_nb_t_F_p(n_proc);
    int nb_t_F=(T-t0)/dt_F;
    std::cout << "nb_t_F : " << nb_t_F << std::endl;

    int rapport = dt_G/dt_F;
    std::cout << "rapport : " << rapport << std::endl;
    for(int i=0; i<n_proc-1; i++){
        tab_nb_t_F_p(i) = tab_nb_t_G_p(i) * rapport;
    }
    tab_nb_t_F_p(n_proc-1) = nb_t_F-tab_nb_t_F_p.head(n_proc-1).sum(); 
    std::cout << "tab_nb_t_F_p : " << tab_nb_t_F_p << std::endl;

    Vector<double> temps_F(n_proc+1);
    temps_F(0) = t0;
    for(int i=1; i<n_proc+1; i++){
        temps_F(i) = temps_F(i-1) + tab_nb_t_F_p(i-1) * dt_F;
    }
    std::cout << "temps_F : " << temps_F << std::endl;

    assert(nb_t_G==tab_nb_t_G_p.sum());
    assert(nb_t_F==tab_nb_t_F_p.sum());

    return 0;
}