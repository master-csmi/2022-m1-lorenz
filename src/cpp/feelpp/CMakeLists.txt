add_executable(laplacian.e laplacian.cpp)

# feelpp_add_application(laplacian.e SRCS laplacian.cpp INCLUDE_IN_ALL DEFS FEELPP_DIM=2)
# feelpp_add_testcase(cases PREFIX qs)

target_link_libraries(
  laplacian.e
  PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen
          Feelpp::feelpp)

target_link_libraries(
  laplacian.e
  PUBLIC Boost::mpi MPI::MPI_CXX)


add_executable(heat.e heat.cpp)

target_link_libraries(
  heat.e
  PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen
          Feelpp::feelpp)

target_link_libraries(
  heat.e
  PUBLIC Boost::mpi MPI::MPI_CXX)