//#include "data_assim.cpp"
#include "data_assim.hpp"
#include <enkf/enkf_fct.hpp>
#include <enkf/enkf.hpp>
#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/mean.hpp>



//#include <enkf/enkf_fct.cpp>
template <typename ToolboxType>
MyMatrix
runToolboxSimulation( std::shared_ptr<ToolboxType> toolbox,double t,MyMatrix X,int nbr_echantillon,double dt)
{
    using namespace Feel;
    double dt_initial=700;
    double dt_final=3300;

    toolbox->init();
    //toolbox->timeStepBase()->setTimeInitial(dt_initial);
    //toolbox->timeStepBase()->setTimeFinal(dt_final);
    //toolbox->timeStepBase()->setTimeInitial(t);
    //toolbox->timeStepBase()->setTimeFinal(t+2*dt);
    toolbox->printAndSaveInfo();
    auto u=toolbox->spaceTemperature()->element();
    std::cout << "t\n  "<<t<<std::endl;
    std::cout << "t+2*dt\n  "<<t+dt<<std::endl;


    MyMatrix X_;
    X_=MyMatrix::Ones(10,1);
    //X_*=10;

    
    //MyMatrix coord{{2.8,1.5,0.7},{4.8,0.8,0.7 },{1.3,3.5,0.7 },{4.8,1.5,0.5},{4.8,3.0,0.7},{3.7,3.5,0.7},{1.3,0.0,0.7},{0.0,0.5,0.7},{3.7,0.0,0.7},{0,2.0,0.7}};
    //std::cout << "coord2\n  "<<coord<<std::endl;
    //double rayon=0.20;
    auto [coord,rayon] = read_coord("/home/u4/csmi/2020/aydogdu/2022-m1-lorenz/src/cpp/feelpp_heat/csv/coord.csv");
    std::string rep="/data/home/aydogdu/feelppdb/data_assim_heat.e/";
    rep+=std::to_string(nbr_echantillon);
    rep+="/np_1/heat.ts/temperature/temperature-6.h5";

    u.load(_path=rep,_type="hdf5");
    for (int i=0;i<10;i++)
    {
        auto gaussien=cst(1./(rayon(i)*sqrt(2*pi)))*exp(-((Px()-cst(coord(i,0)))*(Px()-cst(coord(i,0)))+((Py()-cst(coord(i,1)))*(Py()-cst(coord(i,1))))+(Pz()-cst(coord(i,2)))*(Pz()-cst(coord(i,2))))/(2*(cst(rayon(i))*cst(rayon(i)))));
        u.on(_range=elements(toolbox->mesh()),_expr=cst(X(i))*gaussien);



    }
    std::cout << "print_1\n  "<<std::endl;
    toolbox->timeStepBdfTemperature()->unknown(0)=u;
    std::cout << "bool\n  "<<toolbox->timeStepBase()->isFinished()<<std::endl;
    for ( toolbox->startTimeStep() ; !toolbox->timeStepBase()->isFinished(); toolbox->updateTimeStep() )
    {
        std::cout << "print_2\n  "<<std::endl;
        if (toolbox->worldComm().isMasterRank())
        {
            std::cout << "============================================================\n";
            std::cout << "time simulation: " << toolbox->time() << "s \n";
            std::cout << "============================================================\n";
        }
        
        
        
        toolbox->solve();
        toolbox->exportResults();
    }
    MyMatrix X_dt;
    X_dt=MyMatrix::Zero(10,1);
    
    for (int i=0;i<10;i++)
    {
        auto mean_exp=((Px()-cst(coord(i,0)))*Px()-cst(coord(i,0))+(Py()-cst(coord(i,1)))*(Py()-cst(coord(i,1)))+(Pz()-cst(coord(i,2)))*(Pz()-cst(coord(i,2)))-
        cst(rayon(i))*cst(rayon(i)))<=0;
        auto A=integrate(_range=elements(toolbox->mesh()),_expr=mean_exp).evaluate()(0,0);
        X_dt(i)=integrate(_range=elements(toolbox->mesh()),_expr=idv(toolbox->fieldTemperature())*mean_exp/cst(A)).evaluate()(0,0);

    }
    std::cout<<"X_dt"<<X_dt<<std::endl;
    return X_dt;
}




template <typename ToolboxType>
MyMatrix
runHeatSimulation(double t,MyMatrix X,int nbr_echantillon,double dt)
{
    using namespace Feel;
    std::string rep="/data/home/aydogdu/feelppdb/data_assim_heat.e/";
    rep+=std::to_string(nbr_echantillon);
    auto heat = ToolboxType::New(_prefix="heat",_repository=rep);
    std::cout << "Toolbox initial"<<std::endl;
    return runToolboxSimulation( heat,t,X,nbr_echantillon,dt);
    //std::cout << "Toolbox initial"<<std::endl;
    //return ;

}
MyMatrix fx_heat(double t,MyMatrix X,int nbr_echantillon,double dt)
{
    using namespace Feel;
    typedef Feel::FeelModels::Heat< Simplex<3,1>,Lagrange<1, Scalar,Continuous,PointSetFekete> > model_type;
    return runHeatSimulation< model_type>(t,X,nbr_echantillon,dt);
}


int main(int argc,char **argv)
{
    using namespace Feel;
	po::options_description heatoptions( "heat options" );
    heatoptions.add( toolboxes_options("heat") );

    Environment env( _argc=argc, _argv=argv,
                     _desc=heatoptions,
                     _about=about(_name="toolboxes_heat",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    typedef Feel::FeelModels::Heat< Simplex<3,1>,Lagrange<1, Scalar,Continuous,PointSetFekete> > model_type;
    


    




    std::string date_heure="2022-01-08 20:00:00";
    int nbr_d_obs=12;
    int nbr_model=73;
    std::string path=std::filesystem::current_path();
    std::cout << "path\n  "<<path<<std::endl;
    //MyMatrix obs=read_obs(path+"/csv/meraki_results.csv",date_heure,nbr_d_obs);
    //MyMatrix model=read_model(path+"/csv/heat_model.csv",nbr_model);
    MyMatrix obs=read_obs("/home/u4/csmi/2020/aydogdu/2022-m1-lorenz/src/cpp/feelpp_heat/csv/meraki_results.csv",date_heure,nbr_d_obs);
    std::cout << "obs\n  "<<obs<<std::endl;
    MyMatrix model=read_model("/home/u4/csmi/2020/aydogdu/2022-m1-lorenz/src/cpp/feelpp_heat/csv/heat_model.csv",nbr_model);
    std::cout << "model\n  "<<model<<std::endl;
    

    //MyMatrix coord;
    //MyMatrix rayon;
    //auto [coord,rayon] = read_coord(path+"../..//csv/coord.csv");
    
    
    
    




    //MyMatrix X_0{{20.},{20.},{20},{20.},{20.},{20},{20.},{20.},{20},{20.}}; 
    MyMatrix X_0=MyMatrix::Ones(10,1);
    X_0*=293.15;
    MyMatrix P{{0.1,0,0,0,0,0,0,0,0,0},{0,0.1,0,0,0,0,0,0,0,0},{0,0,0.1,0,0,0,0,0,0,0},{0,0,0,0.1,0,0,0,0,0,0},{0,0,0,0,0.1,0,0,0,0,0},
                {0,0,0,0,0,0.1,0,0,0,0},{0,0,0,0,0,0,0.1,0,0,0},{0,0,0,0,0,0,0,0.1,0,0},{0,0,0,0,0,0,0,0,0.1,0},{0,0,0,0,0,0,0,0,0,0.1}};
    //MyMatrix Q{{0.1,0,0,0,0,0,0,0,0,0},{0,0.1,0,0,0,0,0,0,0,0},{0,0,0.1,0,0,0,0,0,0,0},{0,0,0,0.1,0,0,0,0,0,0},{0,0,0,0,0.1,0,0,0,0,0},
    //            {0,0,0,0,0,0.1,0,0,0,0},{0,0,0,0,0,0,0.1,0,0,0},{0,0,0,0,0,0,0,0.1,0,0},{0,0,0,0,0,0,0,0,0.1,0},{0,0,0,0,0,0,0,0,0,0.1}};
    MyMatrix R{{0.01,0,0,0,0,0,0,0,0,0},{0,0.01,0,0,0,0,0,0,0,0},{0,0,0.01,0,0,0,0,0,0,0},{0,0,0,0.01,0,0,0,0,0,0},{0,0,0,0,0.01,0,0,0,0,0},
                {0,0,0,0,0,0.01,0,0,0,0},{0,0,0,0,0,0,0.01,0,0,0},{0,0,0,0,0,0,0,0.01,0,0},{0,0,0,0,0,0,0,0,0.01,0},{0,0,0,0,0,0,0,0,0,0.01}};
    MyMatrix Q=MyMatrix::Ones(10,10);
    Q*=0.5;
    MyMatrix identity;
    identity=MyMatrix::Identity(10,10);
    Q+=identity;
    int T=12;
    double t=0;
    double N=int(T/0.1); //model
    double dt=0.1; //model
    
    int N2=int(T/1); //obs
    int dt2=1; //obs
    double dim_x=10;
    double dim_z=10;

    EnsembleKalmanFilter E1(dim_x,dim_z,X_0,P,dt2, 2, &hx_heat,&fx_heat);
    E1.set_Q(Q);
    E1.set_R(R);
    double time=0;

    MyMatrix etat,tab_temps;
    etat=MyMatrix::Zero(obs.rows(),obs.cols());
    tab_temps=MyMatrix::Zero(model.rows(),1);

    int dim_col=model.cols();
    
    etat.row(0)=((E1.get_x()).col(0)).transpose();
    tab_temps(0,0)=time;
    int compteur=1;
    //while(time<T-dt2)
    while(time<1)
    {
        MyMatrix z=read_sensor_model(compteur,obs);
        E1.predict_2();
        E1.update(z);
        //std::cout << "z:\n  "<<z<<std::endl;
        time+=dt2;
        etat.row(compteur)=((E1.get_x()).col(0)).transpose();
        tab_temps(compteur,0)=time;
        compteur+=1;
        std::cout << "etat_chaque_iter\n  "<<etat<<std::endl;

    }
    std::cout << "etat\n  "<<etat<<std::endl;
    //E1.creation_csv(lorenz_1,0.01,lorenz_2,0.1,etat,0.1); 
    MyMatrix model_after_data_ass=read_model("heat_model.csv",nbr_model);
    for (int i=1;i<nbr_d_obs;i++)
    {
        model_after_data_ass.row(i*6)=etat.row(i);
    }
    //std::cout << "model after correction\n  "<<model_after_data_ass<<std::endl;
    return 0;
}