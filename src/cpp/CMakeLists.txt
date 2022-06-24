
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11)
add_subdirectory(parareal)
add_subdirectory(enkf)