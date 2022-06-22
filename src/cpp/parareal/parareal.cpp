#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <mpi.h>
#include <parareal/parareal.hpp>
#include <parareal/utils.hpp>

Matrix parareal(Vector X0_t0, double t0, double T, Vector prob(double, 
        Vector, int, double*), double dt_G, double dt_F, double* gamma, 
        int world_rank, const int n_proc, bool write_csv){

    int dim = static_cast<int>(X0_t0.cols());

    int size_verif;
    if(world_rank==0){
        size_verif = (T-t0)/dt_F;
        // std::cout << "(T-t0)/dt_F : " << size_verif << std::endl;
    }

    /*
        we initialize the time intervals for each processing units
    */

    // [t_0,...,t_P]
    Vector times;
    // [t_0,...,t_{P-1}]
    Vector times_head;
    // [t_1,...,t_P]
    Vector times_tail;
    if(world_rank==0){
        times = compute_times(t0,T,dt_G,n_proc);
        // std::cout << times << std::endl;
        times_head = times.head(n_proc);
        times_tail = times.tail(n_proc);
    }
    
    // starting time and final time on each interval
    double t_j, t_jp;

    MPI_Scatter(times_head.data(), 1, MPI_DOUBLE, &t_j, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(times_tail.data(), 1, MPI_DOUBLE, &t_jp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // number of time step on each interval
    int nb_t_p = static_cast<int>((t_jp-t_j)/dt_F);
    if((t_jp-t_j)/dt_F - nb_t_p >1e-6)
        nb_t_p++;
    // std::cout << world_rank << " : " << (t_jp-t_j)/dt_F << std::endl;
    // std::cout << world_rank << " : " << nb_t_p << std::endl;

    /*
        Usefull for communication
    */
    int nEltByProc[n_proc];
    MPI_Gather(&nb_t_p, 1, MPI_INTEGER, nEltByProc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    int displacement[n_proc];
    if(world_rank==0){
        for(int i=0; i<n_proc; i++){
            nEltByProc[i] = nEltByProc[i] * dim;
        }

        displacement[0] = 0;
        for(int i=1; i<n_proc; i++){
            displacement[i] = displacement[i-1]+nEltByProc[i-1];
        }
    }

    /*
        Iteration k = 0
    */

    // initial point of the process j
    Vector X0_k_j(dim);
    // initial points for each process ( tab of X0_k_j )
    Matrix X0_k(n_proc,dim); 
    // final point of the coarse integrator on each interval
    // starting with the initial point X0_k_j
    Vector coarse_k_j(dim);
    // final point of the coarse integrator on each interval
    // starting with the initial points in X0_k
    Matrix coarse_k(n_proc-1,dim); 

    /*
        1st step : compute the X0s for each tj (in sequential)
        (using the coarse integrator)
    */
    if(world_rank==0){
        X0_k.row(0) = X0_t0;
        X0_k_j = X0_t0;

        for (int j=1; j<n_proc; j++){
            coarse_k_j = RK4(X0_k.row(j-1),dt_G,times[j-1],nEltByProc[j-1]/dim,prob,gamma).bottomRows<1>();
            MPI_Send(coarse_k_j.data(), dim, MPI_DOUBLE, j, j, MPI_COMM_WORLD);

            X0_k.row(j) = coarse_k_j;
            coarse_k.row(j-1) = coarse_k_j;
        }
    }
    // final coarse point of proc j is the initial point of proc j+1
    if(world_rank > 0)
        MPI_Recv(X0_k_j.data(), dim, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /*
        2nd step : compute the solution on each interval
        (using the fine integrator)
    */

    // fine solution on each interval
    Matrix fine = RK4(X0_k_j,dt_F,t_j,nb_t_p,prob,gamma);
    // std::cout << "fine : " << fine.size()/dim << std::endl;
    
    // final point of the fine integrator on each interval
    Vector fine_k_j = fine.bottomRows<1>();
    
    // final point of the fine integrator for each interval
    Matrix fine_k(n_proc, dim);
    MPI_Gather(fine_k_j.data(), dim, MPI_DOUBLE, fine_k.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // to save the solution
    // std::cout << "rank=" << world_rank << " ; nb_t_p=" << nb_t_p << std::endl;
    // std::cout << "rank=" << world_rank << " ; fine.rows()=" << fine.rows() << std::endl;

    // assert(nb_t_p==static_cast<int>(fine.rows()));

    int nb_t;
    MPI_Reduce(&nb_t_p,&nb_t,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    // if(world_rank==0)
    //     assert(size_verif==nb_t);

    // fine solution between t0 and T ( only proc 0 )
    Matrix sol_k(nb_t,dim);
    if(world_rank==0)
        assert(size_verif==nb_t);
    MPI_Gatherv(fine.data(), static_cast<int>(fine.size()), MPI_DOUBLE, sol_k.data(), nEltByProc, displacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Matrix solution(nb_t,dim);
    Matrix init_pts(n_proc,dim);
    if(write_csv){
        if(world_rank==0){
            // std::cout << nb_t << std::endl;
            solution = sol_k;
            solution = solution.reshaped<Eigen::RowMajor>(1,solution.cols()*solution.rows());
            // solution.resize(1,solution.cols()*solution.rows());

            init_pts = X0_k;
            init_pts = init_pts.reshaped<Eigen::RowMajor>(1,init_pts.cols()*init_pts.rows());
            // init_pts.resize(1,init_pts.cols()*init_pts.rows());
        }
    }

    // Following iterations (until the solution converge)

    int k = 1;
    bool converge = false;

    Matrix X0_kp(n_proc,dim);
    X0_kp = Matrix::Zero(n_proc,dim);

    Vector to_send(dim);

    while(not converge){             
        if(world_rank==0){
            X0_kp.row(0) = X0_k.row(0);
            X0_k_j = X0_k.row(0);
            for(int j=1; j<n_proc; j++){
                coarse_k_j = RK4(X0_kp.row(j-1),dt_G,times[j-1],nEltByProc[j-1]/dim,prob,gamma).bottomRows<1>();
                X0_kp.row(j) = coarse_k_j + fine_k.row(j-1) - coarse_k.row(j-1);
                to_send = X0_k.row(j);
                MPI_Send(to_send.data(), dim, MPI_DOUBLE, j, j, MPI_COMM_WORLD);
                coarse_k.row(j-1) = coarse_k_j;
            }
        }

        if(world_rank>0)
            MPI_Recv(X0_k_j.data(), dim, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        fine = RK4(X0_k_j,dt_F,t_j,nb_t_p,prob,gamma);

        fine_k_j = fine.bottomRows<1>();

        MPI_Gather(fine_k_j.data(), dim, MPI_DOUBLE, fine_k.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(world_rank==0)
            converge = sol_converge(X0_k,X0_kp);
        
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&converge, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

        if(world_rank==0)
            X0_k = X0_kp; 

        // to save the solution
        MPI_Gatherv(fine.data(), static_cast<int>(fine.size()), MPI_DOUBLE, sol_k.data(), nEltByProc, displacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(write_csv){
            if(world_rank==0){
                solution.conservativeResize(solution.rows()+1,solution.cols());
                init_pts.conservativeResize(init_pts.rows()+1,init_pts.cols());

                solution.row(k) = sol_k.reshaped<Eigen::RowMajor>(1,sol_k.rows()*sol_k.cols());
                init_pts.row(k) = X0_k.reshaped<Eigen::RowMajor>(1,X0_k.rows()*X0_k.cols());
            }
        }

        k++;
    }

    if(world_rank==0){
        // std::cout << "sol : " <<  solution << std::endl << std::endl;
        // std::cout << "init_pts : " <<  init_pts << std::endl << std::endl;
        std::cout << k << " itérations" << std::endl;
    }

    #if 0
    if(write_csv){
        if(world_rank==0)
            csv_files_lorenz(t0, T, dt_F, times, nEltByProc, k, solution, init_pts);
    }
    #endif
    
    return sol_k;
}