#include <Eigen/Dense>

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> Vector;

Vector lorenz(double t, Vector X, int dim, double* gamma);

Matrix* RK4(Vector X0, double dt, double t0, double T, Vector prob(double, Vector, 
        int, double*), double* gamma);

Vector* compute_times(double t0, double T, double dt_G, int P);

bool sol_converge(Matrix* X0_k, Matrix* X0_knext, double eps=1e-9);

Matrix* parareal(Vector X0_t0, double t0, double T, Vector prob(double, 
        Vector, int, double*), double dt_G, double dt_F, double* gamma, 
        int world_rank, int n_proc, bool write_csv = true);