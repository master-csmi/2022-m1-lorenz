#include <Eigen/Dense>

/**
 * @file utils.hpp
 * @author Frédérique Lecourtier
 * @brief Utilities for the parareal algorithm.
 * @version 1.0
 * @date 2022-06-22
 */

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
template<class T> using Vector = Eigen::Matrix<T,1,Eigen::Dynamic>;

Vector<double> lorenz(double t, Vector<double> X, int dim, double* gamma);

Matrix RK4(Vector<double> X0, double dt, double t0, int nb_t, Vector<double> prob(double, Vector<double>, 
        int, double*), double* gamma);

Vector<double> compute_times(double t0, double T, double dt_G, double dt_F, int n_proc,
        Vector<int>* tab_nb_t_G_p, Vector<int>* tab_nb_t_F_p);

bool sol_converge(Matrix X0_k, Matrix X0_knext, double eps=1e-9);