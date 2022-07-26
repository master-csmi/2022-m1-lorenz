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
#include <EigenRand/EigenRand>

using namespace Eigen;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;

MyMatrix read_obs(std::string donnée,std::string date,std::string heure);
MyMatrix read_model(std::string donnée,int nbr_model);

//MyMatrix hx_heat();
