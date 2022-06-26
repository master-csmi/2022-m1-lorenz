#include <Eigen/Core>

#include <Eigen/Core>
#define BOOST_TEST_MODULE test_eigen
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(test_eigen_0)
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
    BOOST_CHECK(a(1,2) == 6);
}
