include_directories(${CMAKE_SOURCE_DIR}/src/cpp)
# link_directories(${CMAKE_SOURCE_DIR}/src/cpp)

add_executable(tests_enkf tests_main.cpp tests_enkf.cpp)
target_link_libraries(tests_enkf PRIVATE project_warnings project_options CONAN_PKG::eigen Boost::unit_test_framework)
add_test(tests_enkf ${PROJECT_BINARY_DIR}/bin/tests_enkf)

add_executable(tests_parareal tests_main.cpp tests_parareal.cpp)
target_link_libraries(tests_parareal PUBLIC parareal PRIVATE project_warnings project_options CONAN_PKG::eigen Boost::unit_test_framework)
add_test(tests_parareal ${PROJECT_BINARY_DIR}/bin/tests_parareal)
