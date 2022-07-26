#include "data_assim.cpp"

int main()
{
    std::string a="2022-01-08 20:00:00";
    int nbr_d_obs=10;
    read_sensor_heat("meraki_results.csv",a,nbr_d_obs);
}