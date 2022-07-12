#include <Eigen/Core>

#include <boost/test/unit_test.hpp>
#include "parareal/parareal.hpp"
#include "parareal/utils.hpp"

BOOST_AUTO_TEST_SUITE(tests_parareal)

BOOST_AUTO_TEST_CASE(test_parareal_0)
{
    typedef Eigen::Matrix<float, 3, 3> MyMatrix33f;
    typedef Eigen::Matrix<float, 3, 1> MyVector3f;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MyMatrix;

    MyMatrix33f a;
    MyVector3f v;
    MyMatrix m(10, 15);

    a = MyMatrix33f::Zero();     // fill matrix elements with zeros
    a = MyMatrix33f::Identity(); // fill matrix as Identity matrix
    v = MyVector3f::Random();    // fill matrix elements with random values

    a << 1, 2, 3,
        4, 5, 6,
        7, 8, 9;
    
    double gamma[3] = {10.,8./3,28.};
    Vector<double> X0(3); X0 << 5., 5., 5.;

    double t0 = 0.;
    double T = 200.;
    double dt_G = 0.01;
    double dt_F = 0.001;

    std::cout << lorenz(t0,X0,gamma) << std::endl;
    BOOST_CHECK(a(1,2) == 6);
}

BOOST_AUTO_TEST_SUITE_END()