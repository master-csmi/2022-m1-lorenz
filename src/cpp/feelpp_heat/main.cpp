#include "data_assim.cpp"

int main()
{
    std::string date_heure="2022-01-08 20:00:00";
    int nbr_d_obs=10;
    int nbr_model=5;
    MyMatrix obs=read_obs("meraki_results.csv",date_heure,nbr_d_obs);
    MyMatrix model=read_model("heat_model.csv",nbr_model);
    std::cout << "model\n  "<<model<<std::endl;
    std::cout << "obs\n  "<<obs<<std::endl;
}