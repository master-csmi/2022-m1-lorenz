add_executable(parareal.e main.cpp parareal.cpp utils.cpp write_csv.cpp)
target_link_libraries(
  parareal.e
  PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen
          )
target_link_libraries(
  parareal.e
  PUBLIC Boost::mpi MPI::MPI_CXX)