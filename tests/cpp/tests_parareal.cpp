#include <Eigen/Core>

#include <boost/test/unit_test.hpp>
#include "parareal/parareal.hpp"
#include "parareal/utils.hpp"

Vector<double> oscillator(double /*t*/, Vector<double> const& X, double* gamma){
    Vector<double> sol(X.cols());
    double w0 = gamma[0];
    sol << X[1], - pow(w0,2) * X[0];

    return sol;
}

double sol_ex(double t, double* gamma){ // gamma = (w0,x0,phi0)
    double sol;
    sol = gamma[1] * cos(gamma[0] * t + gamma[2]);
    return sol;
}

// def sol_ex(t,gamma):
//     (w0,x0,phi0)=gamma
//     sol_x=x0*np.cos(w0*t+phi0)
//     return sol_x

using Matrix = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
template<class T> using Vector = Eigen::Matrix<T,1,Eigen::Dynamic>;

BOOST_AUTO_TEST_SUITE(tests_parareal)

BOOST_AUTO_TEST_CASE(test_parareal_0){
    double gamma[3] = {10.,8./3,28.};
    Vector<double> X0(3); X0 << 5., 5., 5.;

    double t0 = 0.;
    double T = 20.;
    double dt_G = 0.01;
    double dt_F = 0.001;

    Matrix sol_k = parareal(X0, t0, T, oscillator, dt_G, dt_F, gamma, 0, 1, false);



    BOOST_CHECK(6 == 6);
}

BOOST_AUTO_TEST_SUITE_END()