# matlab-nlopt

A Matlab wrapper package of [**NLopt Nonlinear Optimization Library**](https://github.com/stevengj/nlopt)

This project aims to create a set of NLopt-based MATLAB functions which are argument-compatible with the counterparts in Mathwork's Optimization Toolbox, namely: `nlopt.fminunc`, `nlopt.fmincon`, and `nlopt.fminbnd`.

Also this project serves as an example of how to use C++ wrapper function `mexObjectHandler` in [my matlab-mexutils project](https://github.com/hokiedsp/matlab-mexutils).

## Content
* [Build & Installation](#build-installation)
* [`nlopt.fminunc`: Unconstrained nonlinear minimization](#nloptfminunc-unconstrained-nonlinear-minimization)
* [`nlopt.fmincon`: Constrained nonlinear minimization with inequality, equality, and bound constraints](#nloptfmincon-constrained_-nonlinear-minimization-with-inequality-equality-and-bound-constraints)
* [`nlopt.fminbnd`: Bounded global nonlinear minimization](#nloptfminbnd-bounded-global-nonlinear-minimization)
* [`nlopt.options`: Object to manage NLopt options](#nloptoptions-object-to-manage-nlopt-options)
* [Useful Links](#useful-links)

## Build & Installation

This project should build and install following the standard CMake build/install procedures right out of box. However, this project has only been validated to build and install in Windows (MSVC+Ninja). *If you encounter any build/installation issue, please post it on the GitHub Issues board and I'll try to addres it. Meanwhile, if you resolved it yourself, especially on non-Windows platforms, please consider initiating a pull request. Thanks!*

Name | Default | Description
---|---|---
CMAKE_INSTALL_PREFIX | ${MATLAB_USER_DIR} | Base installation directory
MATLAB_NLOPT_TOOLBOX_DIR | "nlopt" | Installation subdirectory for Matlab-NLopt packaged functions (relative to CMAKE_INSTALL_PREFIX)
MATLAB_NLOPT_EXAMPLE_DIR | "nlopt" | Installation subdirectory for Matlab-NLopt examples (relative to CMAKE_INSTALL_PREFIX)
BUILD_NLOPT_LIBS | ON | Also build and install NLopt (as a shared library). For a Windows install, the `nlopt.dll` file is placed within MATLAB_NLOPT_TOOLBOX_DIR. If OFF, CMake's find_package() must be able to located the NLopt installation and its path is on the system PATH for MATLAB to find it. For non-Windows platform, if find_Package() locates the NLopt installation, this option is ignored.
NLOPT_LIB_DIR | "" | If NLopt already not built together and installed, specify the location of NLopt library directory
NLOPT_INCLUDE_DIR | "" | If NLopt not built together, specify the location of NLopt include directory
---

MATLAB_USER_DIR

| Platform | Path |
| Windows | $ENV{USERPROFILE}/Documents/MATLAB |
| Others | $ENV{home}/Documents/MATLAB |




The default intallation path for this project is set to `<MATLAB Users Directory>/nlopt`

## `nlopt.fminunc`: Unconstrained nonlinear minimization

## `nlopt.fmincon`: Constrained nonlinear minimization with inequality, equality, and bound constraints

## `nlopt.fminbnd`: Bounded global nonlinear minimization

## `nlopt.options`: Object to manage NLopt options

### `nlopt.options.getAlgorithms()`: Get the complete list of NLopt algorithms

## Useful Links
