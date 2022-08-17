#include "parareal.hpp"
#include "utils.hpp"
#include "write_csv.hpp"
#include <mpi.h>
#include <ctime>
#include <cmath>
#include <iostream>

// RK4

// int main (){
//     double t0 = 0.;
//     double T = 500.;
//     double dt_F = 0.01;
//     double gamma[3] = {10.,8./3,28.};

//     // Lorenz

//     Vector<double> X0(3); 
//     for(int i=0; i<3; i++){
//         X0[i] = 5.;
//     }

//     Matrix sol_k = RK4(X0,dt_F,t0,static_cast<int>((T-t0)/dt_F),lorenz,gamma);
    
//     std::cout << "avant : " << sol_k.row(sol_k.rows()-2).leftCols(3) << std::endl;
//     std::cout << sol_k.row(sol_k.rows()-1).leftCols(3) << std::endl;

//     return 0;
// }

// Lorenz

int main (){
    MPI_Init ( nullptr , nullptr );
    int world_rank, n_proc ;
    MPI_Comm_rank ( MPI_COMM_WORLD , & world_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & n_proc );

    // if(world_rank==0)
    //     delete_old_files();

    bool write = false;

    double t0 = 0.;
    double T = 500.;
    double dt_G = 0.1;
    double dt_F = 0.01;

    double gamma[3] = {10.,8./3,28.};

    // Lorenz

    // Vector<double> X0(3); X0 << 5., 5., 5.;
    // Matrix sol_k = parareal(X0, t0, T, lorenz, dt_G, dt_F, gamma, world_rank, n_proc, write);

    // Big Lorenz with 100*3 equations

    Vector<double> X0(3*100); 
    for(int i=0; i<3*100; i++){
        X0[i] = 5.;
    }

    time_t start_time, final_time, time_passed;
    if(!write){
        time(&start_time);
    }

    Matrix sol_k = parareal(X0, t0, T, big_lorenz, dt_G, dt_F, gamma, world_rank, n_proc, write);
    if(world_rank==0){
        std::cout << sol_k.rows() << std::endl;
    }

    if(!write){
        time(&final_time);

        if(world_rank==0){
            time_passed = final_time - start_time;
            std::cout << time_passed << std::endl;
            std::cout << "Time : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
        }
    }

    if(world_rank==0){
        // std::cout << "avant : " << sol_k.row(sol_k.rows()-2).leftCols(3) << std::endl;
        // std::cout << sol_k.row(sol_k.rows()-1).leftCols(3) << std::endl;
        std::cout << sol_k.leftCols(3) << std::endl;
    }

    MPI_Finalize ();

    return 0;
}


//oscillator

// int main (){
//     MPI_Init ( nullptr , nullptr );
//     int world_rank, n_proc ;
//     MPI_Comm_rank ( MPI_COMM_WORLD , & world_rank );
//     MPI_Comm_size ( MPI_COMM_WORLD , & n_proc );

//     double gamma[3] = {5.,-1./5.,M_PI/2.};
//     Vector<double> X0(2); X0 << 0.,1.;

//     double t0 = 0.;
//     double T = 20.;
//     double dt_G = 0.01;
//     double dt_F = 0.001;

//     bool write = false;

//     // if(world_rank==0)
//     //     delete_old_files();

//     time_t start_time, final_time, time_passed;
//     if(!write){
//         time(&start_time);
//         // std::cout << "Start : " << start_time << std::endl;
//     }

//     Matrix sol_k = parareal(X0, t0, T, oscillator, dt_G, dt_F, gamma, world_rank, n_proc, write);

//     if(!write){
//         time(&final_time);
//         // std::cout << "Final : " <<  final_time << std::endl;

//         if(world_rank==0){
//             time_passed = final_time - start_time;
//             std::cout << time_passed << std::endl;
//             std::cout << "Time : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
//         }
//     }

//     MPI_Finalize ();

//     return 0;
// }