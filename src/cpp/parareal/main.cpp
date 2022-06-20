#include <parareal/parareal.hpp>

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