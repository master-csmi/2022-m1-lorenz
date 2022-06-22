add_executable(enkf.e main.cpp enkf.cpp enkf_fct.cpp)
#https://github.com/bab2min/EigenRand/releases 
# version 0.4.0 alpha
set(CMAKE_CXX_FLAGS " -std=c++17 -I /Users/macosx/desktop/eigen.0")
target_link_libraries(
  enkf.e
  PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen)