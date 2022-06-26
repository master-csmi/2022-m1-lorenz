add_executable(enkf.e main.cpp enkf.cpp enkf_fct.cpp)
#https://github.com/bab2min/EigenRand/releases 
# version 0.4.0 alpha

#include_directories("/Users/macosx/desktop/pybind11")



#include(FetchContent)
#FetchContent_Declare(
#    pybind11
#    GIT_REPOSITORY https://github.com/pybind/pybind11
#    GIT_TAG        v2.2.3
#)

#FetchContent_GetProperties(pybind11)
#if(NOT pybind11_POPULATED)
#    FetchContent_Populate(pybind11)
#    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
#endif()
# pybind11_add_module(_enkf enkf_wrapper.cpp)
target_link_libraries(
  enkf.e
  PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen)
target_include_directories(enkf.e PUBLIC ${PROJECT_SOURCE_DIR}/src/cpp/EigenRand)