add_subdirectory(laplacian)
add_subdirectory(heat)

add_executable(heat.j heat.cpp)

target_link_libraries(
    heat.j
    PRIVATE project_options
    project_warnings
    CONAN_PKG::docopt.cpp
    CONAN_PKG::fmt
    CONAN_PKG::spdlog
    CONAN_PKG::eigen
    Feelpp::feelpp)