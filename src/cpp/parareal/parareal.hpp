#include <Eigen/Dense>

/**
 * @file parareal.hpp
 * @author Frédérique Lecourtier
 * @brief 
 * @version 1.0
 * @date 2022-06-22
 */

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
template<class T> using Vector = Eigen::Matrix<T,1,Eigen::Dynamic>;

Matrix parareal(Vector<double> X0_t0, double t0, double T, Vector<double> prob(double, 
        Vector<double>, int, double*), double dt_G, double dt_F, double* gamma, 
        int world_rank, int n_proc, bool write_csv=true);