// #include <parareal/write_csv.hpp>
#include "write_csv.hpp"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>

#include <stdio.h>
#include <dirent.h>

void delete_old_files(){
    DIR *theFolder = opendir("examples/cpp/parareal/data");
    struct dirent *next_file;
    char filepath[256];

    while ( (next_file = readdir(theFolder)) != NULL ){
        sprintf(filepath, "%s/%s", "examples/cpp/parareal/data", next_file->d_name);
        remove(filepath);
    }
    closedir(theFolder);
}

void write_sol_k(int k, Vector<double> t, Matrix sol_k){
    std::ostringstream filename; 
    filename << "examples/cpp/parareal/data/solution_" << k << ".csv";

    std::ofstream ofile(filename.str());
    for(int i=0; i < sol_k.rows(); i++){
        ofile << t(i) << ", ";
        for(int j=0; j < sol_k.cols()-1; j++){
            ofile << sol_k(i,j) << ", ";
        }
        ofile << sol_k(i,sol_k.cols()-1) << std::endl;
    }
    ofile.close();
}

void write_X0_k(int k, Vector<double> times, Matrix X0_k){
    std::ostringstream filename; 
    filename << "examples/cpp/parareal/data/init_pts_" << k << ".csv";

    std::ofstream ofile(filename.str());
    for(int i=0; i < X0_k.rows(); i++){
        ofile << times(i) << ", ";
        for(int j=0; j < X0_k.cols()-1; j++){
            ofile << X0_k(i,j) << ", ";
        }
        ofile << X0_k(i,X0_k.cols()-1) << std::endl;
    }
    ofile.close();
}

void write_sol_ex(Vector<double> t, Matrix sol_ex){
    std::ostringstream filename; 
    filename << "examples/cpp/parareal/data/solution_ex.csv";
    std::ofstream ofile(filename.str());

    for(int i=0; i < sol_ex.rows(); i++){
        ofile << t(i) << ", ";
        for(int j=0; j < sol_ex.cols()-1; j++){
            ofile << sol_ex(i,j) << ", ";
        }
        ofile << sol_ex(i,sol_ex.cols()-1) << std::endl;
    }
    ofile.close();
}