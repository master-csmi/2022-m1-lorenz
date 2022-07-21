#include <Eigen/Dense>

/**
 * @file utils.hpp
 * @author Frédérique Lecourtier
 * @brief Utilities for the parareal algorithm.
 */

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
template<class T> using Vector = Eigen::Matrix<T,1,Eigen::Dynamic>;

/**
 * @brief Compute the value of the Harmonic oscillator.
 * 
 * @param t Current time.
 * @param X Current point.
 * @param gamma Parameters.
 * @return Vector<double> The value of the Harmonic oscillator.
 */
Vector<double> oscillator(double t, Vector<double> const& X, double* gamma);

/**
 * @brief Compute the exact solution at the current time of the Harmonic oscillator.
 * 
 * @param t Current time.
 * @param gamma Parameters.
 * @return double The exact value the Harmonic oscillator.
 */
double oscillator_ex(double t, double* gamma);

/**
 * @brief Compute the value of the Lorenz system.
 * 
 * @param t Current time.
 * @param X Current point.
 * @param gamma System parameters.
 * @return Vector<double> The value of the Lorenz system.
 */
Vector<double> lorenz(double t, Vector<double> const& X, double* gamma);

/**
 * @brief Runge Kutta order 4 method.
 * 
 * @param X0 Initial point of the system.
 * @param dt Time step.
 * @param t0 Starting time.
 * @param nb_t Final time.
 * @param prob Function which represent the ODE.
 * @param gamma System parameters.
 * @return Matrix Solution calculated by the method.
 */
Matrix RK4(Vector<double> const& X0, double dt, double t0, int nb_t, 
        Vector<double> prob(double, Vector<double> const&, double*), double* gamma);

/**
 * @brief Compute time used to the system resolution.
 * 
 * @param t0 Starting time.
 * @param T Final time.
 * @param dt_G Coarse time step.
 * @param dt_F Fine time step.
 * @param n_proc Number of processes.
 * @param tab_nb_t_G_p Number of coarse time step for each process.
 * @param tab_nb_t_F_p Number of coarse time step for each process.
 * @return Vector<double> Times used to the system resolution.
 */
Vector<double> compute_times(double t0, double T, double dt_G, double dt_F, int n_proc,
        Vector<int>* tab_nb_t_G_p, Vector<int>* tab_nb_t_F_p);

/**
 * @brief Check if initial points converge.
 * 
 * @param X0_k Initial points at the previous time.
 * @param X0_knext Initial points at the current time.
 * @param eps Error tolerance. Defaults to 1e-9.
 * @return true 
 * @return false 
 */
bool sol_converge(Matrix const& X0_k, Matrix const& X0_knext, double eps=1e-9);