#include <Eigen/Dense>

/**
 * @file utils.hpp
 * @author Frédérique Lecourtier
 * @brief Utilities for the parareal algorithm.
 * @version 1.0
 * @date 2022-06-22
 */

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using Vector = Eigen::Matrix<double,1,Eigen::Dynamic>;

Vector lorenz(double t, Vector X, int dim, double* gamma);

Matrix RK4(Vector X0, double dt, double t0, int nb_t, Vector prob(double, Vector, 
        int, double*), double* gamma);

Vector compute_times(double t0, double T, double dt_G, int P);

bool sol_converge(Matrix X0_k, Matrix X0_knext, double eps=1e-9);