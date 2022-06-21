#include <Eigen/Dense>

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> Vector;

Matrix parareal(Vector X0_t0, double t0, double T, Vector prob(double, 
        Vector, int, double*), double dt_G, double dt_F, double* gamma, 
        int world_rank, int n_proc, bool write_csv=true);