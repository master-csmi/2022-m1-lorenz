#include"data_assim.hpp"
#include <feel/feelmodels/heat/heat.hpp>
//#include <enkf/enkf_fct.cpp>


void split(const std::string &chaine, char delimiteur, std::vector<std::string> &elements)
{
    std::stringstream ss(chaine);
    std::string sousChaine;
    while (getline(ss, sousChaine, delimiteur))
    {
        elements.push_back(sousChaine);
    }
}
std::vector<std::string> split(const std::string &chaine, char delimiteur) 
{
    std::vector<std::string> elements;
    split(chaine, delimiteur, elements);
    return elements;
}
MyMatrix read_obs(std::string donnée,std::string date_heure,int nbr_d_obs)
{
    MyMatrix MT=MyMatrix::Zero(nbr_d_obs,10);
    int compteur=0;
    std :: ifstream ifile(donnée,std :: ios ::in);
    if (ifile.good())
    {
        std::string str;
        getline(ifile, str);
        while(getline(ifile, str)) 
        {
            std::istringstream ss(str);
            int num;
            std::string a;
            while(ss >> num)
            {
                std::vector<std::string> x=split(str, ',');
                if(x[0]==date_heure)
                {
                    while(compteur != nbr_d_obs)
                    {
                        MT(compteur,0)=stod(x[3]);
                        MT(compteur,1)=stod(x[1]);
                        MT(compteur,2)=stod(x[4]);
                        MT(compteur,3)=stod(x[2]);
                        MT(compteur,4)=stod(x[5]);
                        MT(compteur,5)=stod(x[6]);
                        MT(compteur,6)=stod(x[7]);
                        MT(compteur,7)=stod(x[8]);
                        MT(compteur,8)=stod(x[9]);
                        MT(compteur,9)=stod(x[10]);
                        getline(ifile, str);
                        x=split(str, ',');
                        compteur+=1;
                    }
                    return MT;
                }
            }
        }
    }
    return MT;
}
MyMatrix read_model(std::string donnée,int nbr_model)
{
    MyMatrix MT=MyMatrix::Zero(nbr_model,10);
    int compteur=0;
    std :: ifstream ifile(donnée  ,std :: ios ::in);

    if (ifile.good())
    {
        std::string str;
        getline(ifile, str);
        std::cout << "ss\n  "<<str<<std::endl;
        while(getline(ifile, str)) 
        {
            std::istringstream ss(str);
            int num;
            std::string a;
            while(ss >> num)
            {
                std::vector<std::string> x=split(str, ',');
                std::cout << "x\n  "<<x[0]<<std::endl;
                MT(compteur,0)=stod(x[0])-273.15;
                MT(compteur,1)=stod(x[1])-273.15;
                MT(compteur,2)=stod(x[2])-273.15;
                MT(compteur,3)=stod(x[3])-273.15;
                MT(compteur,4)=stod(x[4])-273.15;
                MT(compteur,5)=stod(x[5])-273.15;
                MT(compteur,6)=stod(x[6])-273.15;
                MT(compteur,7)=stod(x[7])-273.15;
                MT(compteur,8)=stod(x[8])-273.15;
                MT(compteur,9)=stod(x[9])-273.15;
                
            }
            compteur+=1;
        }
    }
    return MT;
}
MyMatrix hx_heat(MyMatrix x)
{
    return x;
}

MyMatrix read_sensor_heat(int index,MyMatrix obs)
{
    MyMatrix z;
    int dim_z=obs.cols();
    z=MyMatrix::Zero(dim_z,1);
    z=obs.row(index).transpose();
    return z;
} 
MyMatrix fx_heat(double t,MyMatrix X,int nbr_echantillon)
{
    
    
    
   return X;
    /*
    int nbr_model=73;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    using namespace Feel;
    using Feel::cout;
    po::options_description laplacianoptions( "Laplacian options" );
    laplacianoptions.add_options()
        ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
        ;

    Environment env( _argc=argc, _argv=argv,
                   _desc=laplacianoptions,
                   _about=about(_name="qs_laplacian",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));

    tic();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    toc("loadMesh");

    tic();
    auto Vh = Pch<2>( mesh ); 
    auto u = Vh->element("u"); 
    auto mu = expr(soption(_name="functions.mu")); // diffusion term 
    auto f = expr( soption(_name="functions.f"), "f" ); 
    auto r_1 = expr( soption(_name="functions.a"), "a" ); // Robin left hand side expression 
    auto r_2 = expr( soption(_name="functions.b"), "b" ); // Robin right hand side expression 
    auto n = expr( soption(_name="functions.c"), "c" ); // Neumann expression 
    auto solution = expr( checker().solution(), "solution" ); 
    auto g = checker().check()?solution:expr( soption(_name="functions.g"), "g" ); 
    auto v = Vh->element( g, "g" ); 
    toc("Vh");

    tic();
    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=f*id(v));
    l+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_2*id(v));
    l+=integrate(_range=markedfaces(mesh,"Neumann"), _expr=n*id(v));
    toc("l");

    tic();
    auto a = form2( _trial=Vh, _test=Vh);
    tic();
    a = integrate(_range=elements(mesh),
                  _expr=mu*inner(gradt(u),grad(v)) );
    toc("a");
    a+=integrate(_range=markedfaces(mesh,"Robin"), _expr=r_1*idt(u)*id(v));
    a+=on(_range=markedfaces(mesh,"Dirichlet"), _rhs=l, _element=u, _expr=g );
    //! if no markers Robin Neumann or Dirichlet are present in the mesh then
    //! impose Dirichlet boundary conditions over the entire boundary
    if ( !mesh->hasAnyMarker({"Robin", "Neumann","Dirichlet"}) )
        a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=g );
    toc("a");

    tic();
    //! solve the linear system, find u s.t. a(u,v)=l(v) for all v
    if ( !boption( "no-solve" ) )
        a.solve(_rhs=l,_solution=u);
    toc("a.solve");

    tic();
    auto e = exporter( _mesh=mesh );
    e->addRegions();
    e->add( "uh", u );
    if ( checker().check() )
    {
        v.on(_range=elements(mesh), _expr=solution );
        e->add( "solution", v );
    }
    e->save();
    toc("Exporter");

    
    std::string path=std::filesystem::current_path();
    MyMatrix model=read_model(path+"/csv/heat_model.csv",nbr_model);
    int dim=model.cols();
    MyMatrix m=MyMatrix::Zero(dim,1);
    m=model.row(t*12).transpose();
    return m;
    */
}

