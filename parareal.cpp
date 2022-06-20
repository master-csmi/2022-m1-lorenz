#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <sstream>

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> Vector;

Vector lorenz(double t, Vector X, int dim, double* gamma){
    Vector sol(X.cols());
    sol << gamma[0] * (X[1]-X[0]), X[0] * (gamma[2]-X[2])-X[1],
         X[0]*X[1]-gamma[1]*X[2];

    return sol;
}

Matrix RK4(Vector X0, double dt, double t0, double T, Vector prob(double, Vector, 
        int, double*), double* gamma){

    int dim = X0.cols();
    Matrix X(1,dim); X << X0;

    double t = t0;
    Vector K1(dim), K2(dim), K3(dim), K4(dim);
    Vector X_prec(dim);
    
    while ( (t+dt)<=T or std::abs(t+dt-T)<1e-6){ 
        X_prec = X.bottomRows<1>();
        K1=lorenz(t, X_prec, dim, gamma);
        K2=lorenz(t+dt/2., X_prec + 1./2. * K1 * dt, dim, gamma);
        K3=lorenz(t+dt/2., X_prec + 1./2. * K2 * dt, dim, gamma);
        K4=lorenz(t+dt, X_prec+ K3 * dt, dim, gamma);

        X.conservativeResize(X.rows()+1, X.cols());
        X.row(X.rows()-1) = X_prec + dt/6.* (K1+2.*K2+2.*K3+K4);

        t+=dt;
    }

    return X;
}

Vector compute_times(double t0, double T, double dt_G, int P){
    // time between t_j and t_{j+1}
    double dt_P = (T-t0)/P;
    // nb points
    int nb_pts = dt_P/dt_G;

    // t_j exact
    Vector times_exact(P+1);
    times_exact[0] = t0;
    for(int j=1; j<P+1; j++){
        times_exact[j] = times_exact[j-1] + dt_P;
    }
    times_exact[P]=T;

    // t_j approach (to be a multiple of dt_G)
    double exact, approach;
    Vector times(P+1);
    times[0] = t0;
    for(int j=1; j<P+1; j++){
        exact=times_exact[j];
        approach=times[j-1]+dt_G*nb_pts;
        if(exact-approach<=abs(exact-(approach+dt_G))) { // lower rounding
            times[j] = times[j-1] + dt_G*nb_pts;
        }
        else { // upper rounding
            times[j] = times[j-1] + dt_G*(nb_pts+1);
        }
    }
    times[P]=T;

    return times;
}

bool sol_converge(Matrix X0_k, Matrix X0_knext, double eps=1e-9){
    Matrix diff = X0_knext-X0_k;
    return diff.cwiseAbs().maxCoeff() < eps;
}

void csv_files_lorenz(double t0, double T, double dt_F, Vector times,
        int* nb_tj, int nb_iter, Matrix* solution, Matrix* init_pts){

    // time between t0 and T for the fine integrator
    // t = np.arange(t0,T+1e-6,dt_F);
    // int nb = (T-t0)/dt_F;
    // std::cout << nb << std::endl;
    // Vector t(nb); t << t0;
    // // while ( t(t.rows()-1)<=T or std::abs(t(t.rows()-1)-T)<1e-6){
    // //     t.conservativeResize(t.cols()+1);
    // //     t(t.rows()-1) = t(t.rows()-2)+dt_F;
    // // }
    // // for(int i=1; i<nb; i++)
    // //     t(i)=t(i-1)+dt_F;

    // std::cout << t << std::endl;

    // std::ostringstream filename; 
    // filename << "../resultats/fichiers_resultats/solution1D_instat_activ_desactiv.csv";
    // std::ofstream ofile(filename.str());
}

Matrix parareal(Vector X0_t0, double t0, double T, Vector prob(double, 
        Vector, int, double*), double dt_G, double dt_F, double* gamma, 
        int world_rank, int n_proc, bool write_csv = true){

    int dim = X0_t0.cols();

    Vector times;
    Matrix X0_k(n_proc,dim); 
    Matrix coarse_k(n_proc-1,dim); // pas de valeur pour j=0

    Vector X0_k_j(dim);
    double t_j, t_jp;

    // Iteration k = 0

    Vector X0_j(dim);
    if(world_rank==0){
        // we initialize the time intervals for each processing units
        times = compute_times(t0,T,dt_G,n_proc);

        // 1st step : compute the X0s for each tj (in sequential)
        // We use the coarse integrator

        X0_k.row(0) = X0_t0;
        X0_k_j = X0_t0;
        for (int j=1; j<n_proc; j++){
            X0_j = RK4(X0_k.row(j-1),dt_G,times[j-1],times[j],prob,gamma).bottomRows<1>();
            X0_k.row(j) = X0_j;
            MPI_Send(X0_j.data(), dim, MPI_DOUBLE, j, j, MPI_COMM_WORLD);
            coarse_k.row(j-1) = X0_j;
        }
    }

    if(world_rank > 0)
        MPI_Recv(X0_k_j.data(), dim, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // MPI_Scatter(times.head(n_proc).data(), 1, MPI_DOUBLE, &t_j, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Scatter(times.tail(n_proc).data(), 1, MPI_DOUBLE, &t_jp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(world_rank==0){
        t_j = times(0);
        t_jp = times(1);
        for(int j=1; j<n_proc; j++){
            MPI_Send(&times(j), 1, MPI_DOUBLE, j, j, MPI_COMM_WORLD);
            MPI_Send(&times(j+1), 1, MPI_DOUBLE, j, n_proc+j, MPI_COMM_WORLD);
        }
    }
    if(world_rank > 0){
        MPI_Recv(&t_j, 1, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&t_jp, 1, MPI_DOUBLE, 0, n_proc+world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    // 2nd step : compute the solution on each interval
    // We use the fine integrator
    
    Matrix fine = RK4(X0_k_j,dt_F,t_j,t_jp,prob,gamma);
    
    Vector fine_k_j = fine.bottomRows<1>();
    
    Matrix fine_k(n_proc, dim);
    MPI_Gather(fine_k_j.data(), dim, MPI_DOUBLE, fine_k.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // to save the solution

    int nEltByProc[n_proc];
    int displacement[n_proc];
    int nb_t_p = fine.rows();
    int nb_t = 0;
    MPI_Gather(&nb_t_p, 1, MPI_INTEGER, nEltByProc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    if(world_rank==0){
        for(int i=0; i<n_proc; i++){
            nb_t+=nEltByProc[i];
            nEltByProc[i] = nEltByProc[i] * dim;
        }
        displacement[0] = 0;
        for(int i=1; i<n_proc; i++){
            displacement[i] = displacement[i-1]+nEltByProc[i-1];
        }
    }
    
    Matrix sol_k(nb_t,dim);
    MPI_Gatherv(fine.data(), fine.size(), MPI_DOUBLE, sol_k.data(), nEltByProc, displacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Matrix solution(nb_t,dim);
    Matrix init_pts(n_proc,dim);
    if(write_csv){
        if(world_rank==0){
            // solution = new Matrix(nb_t,dim);
            solution = sol_k;
            // solution.reshaped<Eigen::RowMajor>(1,solution.cols()*solution.rows());
            solution.resize(1,solution.cols()*solution.rows());

            // init_pts = new Matrix(n_proc,dim);
            init_pts = X0_k;
            // init_pts.reshaped<Eigen::RowMajor>(1,init_pts.cols()*init_pts.rows());
            init_pts.resize(1,init_pts.cols()*init_pts.rows());
        }
    }

    // Following iterations (until the solution converge)

    int k = 1;
    bool converge = false;

    Matrix X0_kp(n_proc,dim);
    X0_kp = Matrix::Zero(n_proc,dim);

    Vector coarse_k_j(dim);
    Vector to_send(dim);

    Matrix* sol_k_resize;

    while(not converge){             
        if(world_rank==0){
            X0_kp.row(0) = X0_k.row(0);
            X0_k_j = X0_k.row(0);
            for(int j=1; j<n_proc; j++){
                coarse_k_j = RK4(X0_kp.row(j-1),dt_G,times[j-1],times[j],prob,gamma).bottomRows<1>();
                X0_kp.row(j) = coarse_k_j + fine_k.row(j-1) - coarse_k.row(j-1);
                to_send = X0_k.row(j);
                MPI_Send(to_send.data(), dim, MPI_DOUBLE, j, j, MPI_COMM_WORLD);
                coarse_k.row(j-1) = coarse_k_j;
            }
        }

        if(world_rank>0)
            MPI_Recv(X0_k_j.data(), dim, MPI_DOUBLE, 0, world_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        fine = RK4(X0_k_j,dt_F,t_j,t_jp,prob,gamma);

        fine_k_j = fine.bottomRows<1>();

        MPI_Gather(fine_k_j.data(), dim, MPI_DOUBLE, fine_k.data(), dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(world_rank==0)
            converge = sol_converge(X0_k,X0_kp);
        
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&converge, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

        if(world_rank==0)
            X0_k = X0_kp; 

        // to save the solution

        MPI_Gatherv(fine.data(), fine.size(), MPI_DOUBLE, sol_k.data(), nEltByProc, displacement, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(write_csv){
            if(world_rank==0){
                #if 1
                solution.conservativeResize(solution.rows()+1,solution.cols());
                sol_k_resize = new Matrix(sol_k.rows(),sol_k.cols());
                *sol_k_resize = sol_k;
                sol_k_resize->resize(1,sol_k_resize->rows()*sol_k_resize->cols());
                solution.row(k) = *sol_k_resize;
                delete sol_k_resize;
                #endif
                // solution.row(k) = sol_k.reshaped<Eigen::RowMajor>(1,sol_k.rows()*sol_k.cols());

                #if 1
                init_pts.conservativeResize(init_pts.rows()+1,init_pts.cols());
                init_pts.row(k) = X0_k.reshaped<Eigen::RowMajor>(1,X0_k.rows()*X0_k.cols());
                #endif
            }
        }

        k++;
    }

    if(world_rank==0){
        std::cout << "sol : " <<  solution << std::endl << std::endl;
        std::cout << "init_pts : " <<  init_pts << std::endl << std::endl;
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


int main (){
    MPI_Init ( nullptr , nullptr );
    int world_rank, n_proc ;
    MPI_Comm_rank ( MPI_COMM_WORLD , & world_rank );
    MPI_Comm_size ( MPI_COMM_WORLD , & n_proc );

    double gamma[3] = {10.,8./3,28.};
    Vector X0(3); X0 << 5., 5., 5.;
    double t0 = 0.;
    double T = 2.;
    double dt_G = 0.1;
    double dt_F = 0.01;

    parareal(X0, t0, T, lorenz, dt_G, dt_F, gamma, world_rank, n_proc);

    MPI_Finalize ();


    return 0;
}