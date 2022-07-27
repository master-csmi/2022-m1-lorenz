#include "data_assim.cpp"
#include "../enkf/enkf.cpp"
#include "../enkf/enkf_fct.cpp"

int main()
{
    std::string date_heure="2022-01-08 20:00:00";
    int nbr_d_obs=12;
    int nbr_model=73;
    MyMatrix obs=read_obs("meraki_results.csv",date_heure,nbr_d_obs);
    MyMatrix model=read_model("heat_model.csv",nbr_model);
    std::cout << "model\n  "<<model<<std::endl;
    std::cout << "obs\n  "<<obs<<std::endl;

    MyMatrix X_0{{20.},{20.},{20},{20.},{20.},{20},{20.},{20.},{20},{20.}}; 
    MyMatrix P{{0.1,0,0,0,0,0,0,0,0,0},{0,0.1,0,0,0,0,0,0,0,0},{0,0,0.1,0,0,0,0,0,0,0},{0,0,0,0.1,0,0,0,0,0,0},{0,0,0,0,0.1,0,0,0,0,0},
                {0,0,0,0,0,0.1,0,0,0,0},{0,0,0,0,0,0,0.1,0,0,0},{0,0,0,0,0,0,0,0.1,0,0},{0,0,0,0,0,0,0,0,0.1,0},{0,0,0,0,0,0,0,0,0,0.1}};
    MyMatrix Q{{0.1,0,0,0,0,0,0,0,0,0},{0,0.1,0,0,0,0,0,0,0,0},{0,0,0.1,0,0,0,0,0,0,0},{0,0,0,0.1,0,0,0,0,0,0},{0,0,0,0,0.1,0,0,0,0,0},
                {0,0,0,0,0,0.1,0,0,0,0},{0,0,0,0,0,0,0.1,0,0,0},{0,0,0,0,0,0,0,0.1,0,0},{0,0,0,0,0,0,0,0,0.1,0},{0,0,0,0,0,0,0,0,0,0.1}};
    MyMatrix R{{0.01,0,0,0,0,0,0,0,0,0},{0,0.01,0,0,0,0,0,0,0,0},{0,0,0.01,0,0,0,0,0,0,0},{0,0,0,0.01,0,0,0,0,0,0},{0,0,0,0,0.01,0,0,0,0,0},
                {0,0,0,0,0,0.01,0,0,0,0},{0,0,0,0,0,0,0.01,0,0,0},{0,0,0,0,0,0,0,0.01,0,0},{0,0,0,0,0,0,0,0,0.01,0},{0,0,0,0,0,0,0,0,0,0.01}};
    int T=12;
    double t=0;
    double N=int(T/0.1); //model
    double dt=0.1; //model
    
    int N2=int(T/1); //obs
    int dt2=1; //obs
    

    EnsembleKalmanFilter E1(10,10,X_0,P,dt, 20, &hx_heat,&fx_heat);
    E1.set_Q(Q);
    E1.set_R(R);
    int M=6;
    int index=M;
    double time=0;

    MyMatrix etat,tab_temps;
    etat=MyMatrix::Zero(model.rows(),model.cols());
    tab_temps=MyMatrix::Zero(model.rows(),1);

    int dim_col=model.cols();
    
    etat.row(0)=((E1.get_x()).col(0)).transpose();
    tab_temps(0,0)=time;
    int compteur=1;





}