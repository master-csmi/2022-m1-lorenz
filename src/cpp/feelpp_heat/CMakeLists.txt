include_directories(${CMAKE_SOURCE_DIR}/src/cpp $ENV{FEELPP_DIR}/cmake/feelpp-toolboxes.package.cmake)

add_executable(data_assim_heat.e main.cpp data_assim.cpp)
target_link_libraries(data_assim_heat.e 
    PUBLIC enkf 
    PRIVATE project_warnings project_options CONAN_PKG::eigen Feelpp::feelpp Feelpp::feelpp_toolbox_heat_lib Feelpp::feelpp_modelcore) 

#target_include_directories(data_assim_heat.e PUBLIC ${PROJECT_SOURCE_DIR}/src/cpp/EigenRand)
#feelpp_add_testcase(csv PREFIX )
