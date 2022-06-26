cmake_minimum_required(VERSION 3.16)

project(lorenz CXX)
include(cmake/StandardProjectSettings.cmake)
include(cmake/PreventInSourceBuilds.cmake)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Link this 'library' to set the c++ standard / compile-time options requested
add_library(project_options INTERFACE)
# target_compile_features(project_options INTERFACE cxx_std_17)
target_compile_features(project_options INTERFACE cxx_std_20)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)
  if(ENABLE_BUILD_WITH_TIME_TRACE)
    target_compile_options(project_options INTERFACE -ftime-trace)
  endif()
endif()

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)

# enable cache system
include(cmake/Cache.cmake)

# standard compiler warnings
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# sanitizer options if supported by compiler
include(cmake/Sanitizers.cmake)
enable_sanitizers(project_options)

# enable doxygen
include(cmake/Doxygen.cmake)
enable_doxygen()


# allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)
option(ENABLE_TESTING "Enable Test Builds" ON)



# Set up some extra Conan dependencies based on our needs before loading Conan
set(CONAN_EXTRA_REQUIRES "")
set(CONAN_EXTRA_OPTIONS "")

# set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES} eigen/3.3.9)
set(CONAN_EXTRA_REQUIRES ${CONAN_EXTRA_REQUIRES} eigen/3.4.0 )

include(cmake/Conan.cmake)
run_conan()

find_package(MPI REQUIRED)

if(NOT Boost_ARCHITECTURE)
    set(Boost_ARCHITECTURE "-x64")
endif()

set(Boost_ADDITIONAL_VERSIONS "1.61" "1.62" "1.63" "1.64" "1.65" "1.66" "1.67" "1.68" "1.69" "1.70" "1.71")
set(BOOST_COMPONENTS_REQUIRED date_time filesystem system program_options unit_test_framework mpi regex serialization )
FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${BOOST_COMPONENTS_REQUIRED})

if(ENABLE_TESTING)
  enable_testing()
  message("Building Tests. Be sure to check out test/constexpr_tests for constexpr testing")
  add_subdirectory(tests/cpp)
endif()
target_include_directories(project_options  INTERFACE ${MPI_CXX_HEADER_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)

add_subdirectory(src/cpp)
add_subdirectory(docs)