

add_executable(tests_eigen tests_eigen.cpp)
target_link_libraries(tests_eigen PRIVATE project_warnings project_options CONAN_PKG::eigen Boost::unit_test_framework)
add_test(tests_eigen ${PROJECT_BINARY_DIR}/bin/tests_eigen)