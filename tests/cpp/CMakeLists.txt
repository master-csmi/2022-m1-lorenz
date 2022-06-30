add_executable(tests_parareal tests_parareal.cpp)
target_link_libraries(tests_parareal PRIVATE project_warnings project_options CONAN_PKG::eigen Boost::unit_test_framework)
add_test(tests_parareal ${PROJECT_BINARY_DIR}/bin/tests_parareal)