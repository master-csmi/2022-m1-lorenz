add_library(
    data_assim_heat SHARED 
    data_assim.cpp 
    )
target_include_directories(data_assim_heat PUBLIC ${PROJECT_SOURCE_DIR}/src/cpp/EigenRand)
include_directories(${CMAKE_SOURCE_DIR}/src/cpp)
target_link_libraries(data_assim_heat PUBLIC CONAN_PKG::eigen)

add_executable(data_assim_heat.e main.cpp)
target_link_libraries(
    data_assim_heat.e
    PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen
          Feelpp::feelpp)

target_link_libraries(
    data_assim_heat.e
    PUBLIC Boost::mpi MPI::MPI_CXX)
target_include_directories(data_assim_heat.e PUBLIC ${PROJECT_SOURCE_DIR}/src/cpp/EigenRand)

target_compile_definitions(data_assim_heat.e PUBLIC FEELPP_DIM=2)