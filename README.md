# 2022-m1-lorenz
## Description of the projet: Data assimilation for the Lorenz system

 The main goals of this project was to implement a parallel time resolution method for the Lorenz system, and to realize the data assimilation using the EnKF method. For this we also had to implement several methods to solve numerically the Lorenz system

## Authors
This project is realized during the master csmi by two Master 1 students, Lecourtier Frederique and AYDOGDU Melissa and it is managed by Cemosis which is the "Centre de Modélisation et de Simulation de Strasbourg" (Strasbourg Modeling and Simulation Center). Cemosis is hosted by the Institute of Advanced Mathematical Research (IRMA) and was created in January 2013. Cemosis relies currently on the team Modeling and Control of the IRMA. Their work is focused on the numarical simulation and mathematical modelling of different phenomena.
Logo: ![Alt](docs\presentation\images\logo-cemosis.pdf)

## Requirements
For the execution of the project it is necessary to install the following modules:
* NumPy [NumPy Documentation](https://numpy.org/doc/ "Title")
```shell
  pip install numpy
```
* FilterPy [FilterPy Documentation](https://filterpy.readthedocs.io/en/latest/kalman/EnsembleKalmanFilter.html "Title")
```shell
  pip install FilterPy==1.4.5
```
* MatPlotLib [MatPlotLib Documentation](https://matplotlib.org/stable/index.html "Title")
```shell
  pip install matplotlib
```
* Scipy [MatPlotLib Documentation](https://docs.scipy.org/doc/scipy/ "Title")
```shell
  pip install matplotlib
```

$$
\begin{aligned}
x'(t) &= \sigma(y-x)\\
y'(t) &= \\
z'(t) &= \\
\end{aligned}
$$

## Authors


## Organisation 

## Documentation

## Examples
=======
[![Python package](https://github.com/master-csmi/2022-m1-lorenz/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/master-csmi/2022-m1-lorenz/actions/workflows/python-package.yml)

