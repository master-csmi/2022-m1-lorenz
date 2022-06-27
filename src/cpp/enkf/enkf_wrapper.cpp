#include "enkf.hpp"
#include <pybind11/pybind11.h>
namespace py = pybind11;



PYBIND11_MODULE(_enkf, m) {
    m.doc() = "welcome to enkf module"; 
    m.def("add", &add, "A function that adds two numbers");
    //m.def("mean", &mean, "A function that returns an average based on the columns (a=0) or rows (a=1)");
    //m.def("f_lorenz", &f_lorenz, "Computes the Lorenz system for a given time, point and parameters.");
    //m.def("hx_model", &hx_model, "Measurement function. Convert state x into a measurement.");
    //m.def("read_sensor_model", &read_sensor_model, "A function that will read the observation every time.");
    //m.def("fx", &fx, "Will compute the next position using RK4.");
    //m.def("fx_2", &fx_2, "Will compute the next position using RK4 for a given time and parameter.");
    //m.def("RK4", &RK4, "This fonction return resolution with RK4 for the lorenz system.");
   



    //py::class_<EnsembleKalmanFilter>(m, "EnsembleKalmanFilter")
    //    .def(py::init<double ,double ,MyMatrix ,MyMatrix ,double ,int ,MyMatrix (*)(MyMatrix),MyMatrix (*)(double, MyMatrix)
    //>());
        //.def("predict", &EnsembleKalmanFilter::predict)
        //.def("getName", &EnsembleKalmanFilter::update);
}       