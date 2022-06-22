#include <Eigen/Dense>

/**
 * @file parareal.hpp
 * @author Frédérique Lecourtier
 * @brief 
 * @version 1.0
 * @date 2022-06-22
 */

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using Vector = Eigen::Matrix<double,1,Eigen::Dynamic>;

Matrix parareal(Vector X0_t0, double t0, double T, Vector prob(double, 
        Vector, int, double*), double dt_G, double dt_F, double* gamma, 
        int world_rank, int n_proc, bool write_csv=true);