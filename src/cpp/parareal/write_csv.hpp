#include <Eigen/Dense>

/**
 * @file write_csv.hpp
 * @author Frédérique Lecourtier
 * @brief To create csv files xith the solution for each iteration.
 */

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
template<class T> using Vector = Eigen::Matrix<T,1,Eigen::Dynamic>;

/**
 * @brief Delete old data files.
 * 
 */
void delete_old_files();

/**
 * @brief Write the solution at iteration k.
 * 
 * @param k Iteration.
 * @param t All fine time step between t0 and T.
 * @param sol_k Soltuion at the iteration k.
 */
void write_sol_k(int k, Vector<double> t, Matrix sol_k);

/**
 * @brief Write the initial points at iteration k.
 * 
 * @param k Iteration.
 * @param times Time used to the system resolution.
 * @param X0_k Initial points at iteration k.
 */
void write_X0_k(int k, Vector<double> times, Matrix X0_k);

/**
 * @brief Write the exact solution.
 * 
 * @param t Time used to the system resolution.
 * @param sol_ex Exact solution.
 */
void write_sol_ex(Vector<double> t, Matrix sol_ex);