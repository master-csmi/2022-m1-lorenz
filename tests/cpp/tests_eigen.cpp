#include <Eigen/Core>

#include <catch2/catch.hpp>



TEST_CASE("eigen", "[basic]")
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
    REQUIRE(a(1,2) == 6);
}
