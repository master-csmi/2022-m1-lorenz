#include <Eigen/Dense>

/**
 * @file parareal.hpp
 * @author Frédérique Lecourtier
 * @brief The parareal method.
 */

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
template<class T> using Vector = Eigen::Matrix<T,1,Eigen::Dynamic>;

/**
 * @brief Parareal method.
 * 
 * @param X0_t0 Initial point.
 * @param t0 Starting time.
 * @param T Final time.
 * @param prob Function which represent the ODE.
 * @param dt_G Coarse time step.
 * @param dt_F Fine time step.
 * @param gamma System parameters.
 * @param world_rank Process rank.
 * @param n_proc Number of processes.
 * @param write_csv True to write the solutions in a csv file. Defaults to True.
 * @return Matrix The solution by using the parareal method.
 */
Matrix parareal(Vector<double> const& X0_t0, double t0, double T, 
                Vector<double> prob(double, Vector<double> const&, double*), 
                double dt_G, double dt_F, double* gamma, 
                int world_rank, int n_proc, bool write_csv=true);