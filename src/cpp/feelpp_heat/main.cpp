//#include "data_assim.cpp"
#include "data_assim.hpp"
#include <enkf/enkf_fct.hpp>
#include <enkf/enkf.hpp>
#include <feel/feelmodels/heat/heat.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/mean.hpp>



//#include <enkf/enkf_fct.cpp>
template <typename ToolboxType>
int
runToolboxSimulation( std::shared_ptr<ToolboxType> toolbox )
{
    using namespace Feel;
    double dt_initial=700;
    double dt_final=3300;

    toolbox->init();
    toolbox->timeStepBase()->setTimeInitial(dt_initial);
    toolbox->timeStepBase()->setTimeFinal(dt_final);
    toolbox->printAndSaveInfo();
    auto u=toolbox->spaceTemperature()->element();

    MyMatrix X;
    X=MyMatrix::Ones(10,1);
    X*=10;


    MyMatrix coord{{2.8,1.5,0.7},{4.8,0.8,0.7 },{1.3,3.5,0.7 },{4.8,1.5,0.5},{4.8,3.0,0.7},{3.7,3.5,0.7},{1.3,0.0,0.7},{0.0,0.5,0.7},{3.7,0.0,0.7},{0,2.0,0.7}};
    double rayon=0.20;

    u.load(_path="/data/home/aydogdu/feelppdb/data_assim_heat.e/1/np_1/heat.ts/temperature/temperature-6.h5",_type="hdf5");
    for (int i=0;i<10;i++)
    {
        auto gaussien=exp(-(pow(Px()-cst(coord(i,0)),2)+pow(Py()-cst(coord(i,1)),2)+pow(Pz()-cst(coord(i,2)),2))/(2*pow(cst(rayon),2)));
        u.on(_range=elements(toolbox->mesh()),_expr=cst(X(i))*gaussien);



    }
    toolbox->timeStepBdfTemperature()->unknown(0)=u;
    for ( toolbox->startTimeStep() ; !toolbox->timeStepBase()->isFinished(); toolbox->updateTimeStep() )
    {
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
        auto mean_exp=(pow(Px()-cst(coord(i,0)),2)+pow(Py()-cst(coord(i,1)),2)+pow(Pz()-cst(coord(i,2)),2)-pow(cst(rayon),2))<=0;
        auto A=integrate(_range=elements(toolbox->mesh()),_expr=mean_exp).evaluate()(0,0);
        std::cout <<"A"<< A<<std::endl;
        X_dt(i)=integrate(_range=elements(toolbox->mesh()),_expr=idv(toolbox->fieldTemperature())*mean_exp/cst(A)).evaluate()(0,0);
        std::cout << X_dt(i)<<std::endl;

    }
        

    return !toolbox->checkResults();
}




template <typename ToolboxType>
int
runHeatSimulation()
{
    using namespace Feel;
    std::string rep="/data/home/aydogdu/feelppdb/data_assim_heat.e";
    rep+="/1";
    auto heat = ToolboxType::New(_prefix="heat",_repository=rep);
    std::cout << "Toolbox initial"<<std::endl;
    return runToolboxSimulation( heat );
    //std::cout << "Toolbox initial"<<std::endl;
    //return ;

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
    int status=0;
    status = runHeatSimulation<model_type>();


    std::string date_heure="2022-01-08 20:00:00";
    int nbr_d_obs=12;
    int nbr_model=73;
    std::string path=std::filesystem::current_path();
    std::cout << "path\n  "<<path<<std::endl;
    MyMatrix obs=read_obs(path+"/csv/meraki_results.csv",date_heure,nbr_d_obs);
    MyMatrix model=read_model(path+"/csv/heat_model.csv",nbr_model);
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
    double dim_x=10;
    double dim_z=10;

    EnsembleKalmanFilter E1(dim_x,dim_z,X_0,P,dt2, 20, &hx_heat,&fx_heat);
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
    while(time<T-dt2)
    {
        MyMatrix z=read_sensor_model(compteur,obs);
        E1.predict();
        E1.update(z);
        //std::cout << "z:\n  "<<z<<std::endl;
        time+=dt2;
        etat.row(compteur)=((E1.get_x()).col(0)).transpose();
        tab_temps(compteur,0)=time;
        compteur+=1;

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