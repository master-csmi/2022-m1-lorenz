add_executable(myproject main.cpp parareal.cpp utils.cpp)
target_link_libraries(
  myproject
  PRIVATE project_options
          project_warnings
          CONAN_PKG::docopt.cpp
          CONAN_PKG::fmt
          CONAN_PKG::spdlog
          CONAN_PKG::eigen)