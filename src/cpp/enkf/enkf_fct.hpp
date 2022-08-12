#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <sstream>
#include<time.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <EigenRand/EigenRand>

/**
 * @file enkf_fct.hpp
 * @author Melissa AYDOGDU
 * @brief Some fonction to apply the enkf method .
 * @version 1.0
 * @date 2022-06-22
 */

using namespace Eigen;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;

int add(int i,int j);
/**
 * @brief A function that returns an average based on the columns (a=0) or rows (a=1)
 * 
 * @param M An eigen matrix containing numbers whose mean is desired.
 * @param a Axis or axes along which the means are computed.
 * @return Returns a new eigen matrix containing the mean values.
 */
MyMatrix mean(MyMatrix M,int a);

/**
 * @brief Computes the Lorenz system for a given time, point and parameters.
 * 
 * @param t An double which represents time.
 * @param X An eigen matrix that represents the positions.
 * @param p An eigen matrix that represents the points.
 * @return An eigen matrix with the next positions.
 */
MyMatrix f_lorenz(double t,MyMatrix X,MyMatrix p);


MyMatrix f_ocillateur(double t,MyMatrix X,MyMatrix p);


/**
 * @brief Measurement function. Convert state x into a measurement.
 * 
 * @param x An eigen matrix with the state x.
 * @return Return an eigen matrix in order to switch from the model space to the observation space.
 */
MyMatrix hx_model(MyMatrix x);

/**
 * @brief A function that will read the observation every time.
 * 
 * @param index A double that represent the index that you want take in the matrix lorenz_1.
 * @param lorenz_1 The observation matrix.
 * @return Return an eigen matrix with the position associated with the observation.
 */
MyMatrix read_sensor_model(int index,MyMatrix lorenz_1);

/**
 * @brief Will compute the next position using RK4.
 * 
 * @param X An eigen matrix with the position.
 * @param t A double representing the time.
 * @param dt A double that represent the time discretisiation.
 * @param p The parameters for the system.
 * @return Return an eigen matrix with the next position.
 */
MyMatrix fx(MyMatrix X,double t,double dt,MyMatrix p);

/**
 * @brief Will compute the next position using RK4 for a given time and parameter.
 * 
 * @param dt A double that represent the time discritisiation.
 * @param X An eigen matrix with the position.
 * @return Return an eigen matrix with the next position.
 */
MyMatrix fx_2(double dt,MyMatrix X,int nbr_echantillon);

/**
 * @brief This fonction return resolution with RK4 for the lorenz system.
 * 
 * @param dim_x A int that represent the dimension of the position.
 * @param p An eigen matrix with the parameters.
 * @param X0 An eigen matrix with the initial position.
 * @param N  An in that represents the number of discritisation.
 * @param T A double that represents the time interval
 * @return Return the matrix with the positions for each discritisition.
 */
MyMatrix RK4(int dim_x,MyMatrix p,MyMatrix X0,int N,double T,MyMatrix (*f)(double t,MyMatrix X0,MyMatrix P));