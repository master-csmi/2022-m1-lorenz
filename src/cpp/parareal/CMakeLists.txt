add_executable(parareal main.cpp parareal.cpp utils.cpp)
target_link_libraries(
  parareal
  PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen)