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
MyMatrix read_sensor_heat(std::string donnée,std::string date_heure,int nbr_d_obs)
{
    MyMatrix MT=MyMatrix::Zero(nbr_d_obs,10);
    int compteur=0;
    std :: ifstream ifile(donnée  ,std :: ios ::in);
    if (ifile.good())
    {
        std::string str;
        ifile;
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
                        
                        //std::cout << "str1\n  "<<str<<std::endl;
                        //std::cout << "x1\n  "<<x[3]<<std::endl;

                        getline(ifile, str);
                        //std::cout << "str2\n  "<<str<<std::endl;
                        //std::istringstream ss(str);
                        x=split(str, ',');
                        //std::cout << "x2\n  "<<x[3]<<std::endl;
                        compteur+=1;
                    }
                    return MT;
                }
            

            }

        }
        
    
    }
    return MT;

}