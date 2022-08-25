#include "parareal.hpp"
#include "utils.hpp"
#include "write_csv.hpp"
#include <mpi.h>
#include <ctime>
#include <cmath>
#include <iostream>
#define RK4_bool 0

// Lorenz

#if 1
int main (){
    MPI_Init ( nullptr , nullptr );
    int world_rank, n_proc ;
    MPI_Comm_rank ( MPI_COMM_WORLD , & world_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & n_proc );

    // if(world_rank==0)
    //     delete_old_files();

    bool write = false;

    double t0 = 0.;
    double T = 10000.;
    double dt_G = 0.01;
    double dt_F = 0.001;

    double gamma[3] = {10.,8./3,28.};

    // Big Lorenz (which is Lorenz with duplicate*3 equations)
    // if duplicate=1 => Lorenz

    int duplicate = 100;
    Vector<double> X0(3*duplicate); 
    for(int i=0; i<3*duplicate; i++){
        X0[i] = 5.;
    }

    time_t start_time, final_time, time_passed;
    if(!write)
        time(&start_time);

    Vector<double> sol = parareal(X0, t0, T, big_lorenz, dt_G, dt_F, gamma, world_rank, n_proc, write);
    if(world_rank==0){
        std::cout << "Solution :" << std::endl;
        std::cout << sol.leftCols(3) << std::endl;
    }

    if(!write){
        time(&final_time);

        if(world_rank==0){
            time_passed = final_time - start_time;
            std::cout << time_passed << std::endl;
            std::cout << "Time : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
        }
    }

    // RK4
    #if RK4_bool
    if(n_proc==1){
            if(!write)
        time(&start_time);
        Matrix sol_rk4 = RK4(X0,dt_F,t0,static_cast<int>((T-t0)/dt_F),big_lorenz,gamma);
        if(!write){
            time(&final_time);
            time_passed = final_time - start_time;
            std::cout << "Time - RK4 : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
        }
        std::cout << "RK4 : " << sol_rk4.row(sol_rk4.rows()-1).leftCols(3) << std::endl;
    }
    #endif

    MPI_Finalize ();

    return 0;
}
#endif

//oscillator

#if 0
int main (){
    MPI_Init ( nullptr , nullptr );
    int world_rank, n_proc ;
    MPI_Comm_rank ( MPI_COMM_WORLD , & world_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & n_proc );

    double gamma[3] = {5.,-1./5.,M_PI/2.};

    double t0 = 0.;
    double T = 10000.;
    double dt_G = 0.01;
    double dt_F = 0.001;

    bool write = true;

    int duplicate = 200;
    Vector<double> X0(2*duplicate); 
    for(int i=0; i<2*duplicate; i++){
        if(i%2==0)
            X0[i] = 0.;
        else if(i%2==1)
            X0[i] = 1.;
    }

    // if(world_rank==0)
    //     delete_old_files();

    // Speed-up ?

    time_t start_time, final_time, time_passed;
    if(!write){
        time(&start_time);
        // std::cout << "Start : " << start_time << std::endl;
    }

    Vector<double> sol = parareal(X0, t0, T, big_oscillator, dt_G, dt_F, gamma, world_rank, n_proc, write);
    if(world_rank==0){
        std::cout << "Solution :" << std::endl;
        std::cout << sol.leftCols(2) << std::endl;
    }

    if(!write){
        time(&final_time);
        // std::cout << "Final : " <<  final_time << std::endl;

        if(world_rank==0){
            time_passed = final_time - start_time;
            std::cout << time_passed << std::endl;
            std::cout << "Time : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
        }
    }

    // RK4
    #if RK4_bool
    if(n_proc==1){
            if(!write)
        time(&start_time);
        Matrix sol_rk4 = RK4(X0,dt_F,t0,static_cast<int>((T-t0)/dt_F),big_oscillator,gamma);
        if(!write){
            time(&final_time);
            time_passed = final_time - start_time;
            std::cout << "Time - RK4 : " << std::floor(time_passed/60) << "m" << time_passed%60 << "s" << std::endl;
        }
        std::cout << "RK4 : " << sol_rk4.row(sol_rk4.rows()-1).leftCols(2) << std::endl;
    }
    #endif

    MPI_Finalize ();

    return 0;
}
#endif