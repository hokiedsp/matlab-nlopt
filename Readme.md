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

**This project has only been validated to build and install in Windows. Please let me know of your success or issues, especially on other platforms.**


## `nlopt.fminunc`: Unconstrained nonlinear minimization

## `nlopt.fmincon`: Constrained nonlinear minimization with inequality, equality, and bound constraints

## `nlopt.fminbnd`: Bounded global nonlinear minimization

## `nlopt.options`: Object to manage NLopt options

### `nlopt.options.getAlgorithms()`: Get the complete list of NLopt algorithms

## Useful Links
