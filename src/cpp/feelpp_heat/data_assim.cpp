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
                        MT(compteur,0)=std::stod(x[3]);
                        MT(compteur,1)=std::stod(x[1]);
                        MT(compteur,2)=std::stod(x[4]);
                        MT(compteur,3)=std::stod(x[2]);
                        MT(compteur,4)=std::stod(x[5]);
                        MT(compteur,5)=std::stod(x[6]);
                        MT(compteur,6)=std::stod(x[7]);
                        MT(compteur,7)=std::stod(x[8]);
                        MT(compteur,8)=std::stod(x[9]);
                        MT(compteur,9)=std::stod(x[10]);
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
        while(getline(ifile, str)) 
        {
            std::istringstream ss(str);
            int num;
            std::string a;
            while(ss >> num)
            {
                std::vector<std::string> x=split(str, ',');
                MT(compteur,0)=std::stod(x[0])-273.15;
                MT(compteur,1)=std::stod(x[1])-273.15;
                MT(compteur,2)=std::stod(x[2])-273.15;
                MT(compteur,3)=std::stod(x[3])-273.15;
                MT(compteur,4)=std::stod(x[4])-273.15;
                MT(compteur,5)=std::stod(x[5])-273.15;
                MT(compteur,6)=std::stod(x[6])-273.15;
                MT(compteur,7)=std::stod(x[7])-273.15;
                MT(compteur,8)=std::stod(x[8])-273.15;
                MT(compteur,9)=std::stod(x[9])-273.15;
                
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
    MyMatrix Kelv;
    int dim_z=obs.cols();
    z=MyMatrix::Zero(dim_z,1);
    Kelv=MyMatrix::Ones(dim_z,1);
    Kelv*=273.15;
    z=obs.row(index).transpose();
    return z+Kelv;
} 


std::tuple<MyMatrix, MyMatrix>  read_coord(std::string donnée)
{
    MyMatrix coord=MyMatrix::Zero(10,3);
    MyMatrix rayon=MyMatrix::Zero(10,1);
    int compteur=0;
    std :: ifstream ifile(donnée  ,std :: ios ::in);

    if (ifile.good())
    {
        std::string str;
        while(getline(ifile, str)) 
        {
            std::istringstream ss(str);
            int num;
            std::string a;
            while(ss >> num)
            {
                std::vector<std::string> x=split(str, ',');
                coord(compteur,0)=std::stod(x[0]);
                coord(compteur,1)=std::stod(x[1]);
                coord(compteur,2)=std::stod(x[2]);
                rayon(compteur,0)=std::stod(x[3]);
            }
            compteur+=1;
        }
    }
    return {coord, rayon};
}