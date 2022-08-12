#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include<time.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <filesystem>

using namespace Eigen;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;

MyMatrix read_obs(std::string donnée,std::string date_heure,int nbr_d_obs);
MyMatrix read_model(std::string donnée,int nbr_model);
MyMatrix hx_heat(MyMatrix x);
MyMatrix read_sensor_heat(int index,MyMatrix obs);
MyMatrix fx_heat(double t,MyMatrix X,int nbr_echantillon);

//MyMatrix hx_heat();
