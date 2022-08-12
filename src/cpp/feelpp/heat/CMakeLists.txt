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

target_compile_definitions(heat.e PUBLIC FEELPP_DIM=2)
