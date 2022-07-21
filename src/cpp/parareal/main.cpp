#include "parareal.hpp"
#include "utils.hpp"
#include "write_csv.hpp"
#include <mpi.h>
#include <ctime>
#include <cmath>
// #include <feel/feel.hpp> //tic/toc ?

int main (){
    MPI_Init ( nullptr , nullptr );
    int world_rank, n_proc ;
    MPI_Comm_rank ( MPI_COMM_WORLD , & world_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & n_proc );

    double gamma[3] = {5.,-1./5.,M_PI/2.};
    Vector<double> X0(2); X0 << 0.,1.;

    double t0 = 0.;
    double T = 20.;
    double dt_G = 0.01;
    double dt_F = 0.001;

    bool write = true;

    // if(world_rank==0)
    //     delete_old_files();

    time_t start_time, final_time, time_passed;
    if(!write){
        time(&start_time);
        std::cout << "Start : " << start_time << std::endl;
    }

    Matrix sol_k = parareal(X0, t0, T, oscillator, dt_G, dt_F, gamma, world_rank, n_proc, write);

    if(!write){
        time(&final_time);
        std::cout << "Final : " <<  final_time << std::endl;

        if(world_rank==0){
            time_passed = final_time - start_time;
            std::cout << time_passed << std::endl;
            std::cout << "Time : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
        }
    }

    MPI_Finalize ();

    return 0;
}

// int main (){
//     MPI_Init ( nullptr , nullptr );
//     int world_rank, n_proc ;
//     MPI_Comm_rank ( MPI_COMM_WORLD , & world_rank );
//     MPI_Comm_size ( MPI_COMM_WORLD , & n_proc );

//     double gamma[3] = {10.,8./3,28.};
//     Vector<double> X0(3); X0 << 5., 5., 5.;

//     double t0 = 0.;
//     double T = 20.;
//     double dt_G = 0.1;
//     double dt_F = 0.01;

//     bool write = true;

//     // if(world_rank==0)
//     //     delete_old_files();

//     time_t start_time, final_time, time_passed;
//     if(!write){
//         time(&start_time);
//         std::cout << "Start : " << start_time << std::endl;
//     }

//     Matrix sol_k = parareal(X0, t0, T, lorenz, dt_G, dt_F, gamma, world_rank, n_proc, write);

//     if(!write){
//         time(&final_time);
//         std::cout << "Final : " <<  final_time << std::endl;

//         if(world_rank==0){
//             time_passed = final_time - start_time;
//             std::cout << time_passed << std::endl;
//             std::cout << "Time : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
//         }
//     }

//     MPI_Finalize ();

//     return 0;
// }