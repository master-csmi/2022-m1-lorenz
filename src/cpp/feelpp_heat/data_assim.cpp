#include"data_assim.hpp"
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
    std::cout << "prob\n  "<<donnée<<std::endl;
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
MyMatrix fx_heat(double t,MyMatrix X)
{
    int nbr_model=73;
    MyMatrix model=read_model("heat_model.csv",nbr_model);
    int dim=model.cols();
    MyMatrix m=MyMatrix::Zero(dim,1);
    m=model.row(t*12).transpose();
    return m;
}