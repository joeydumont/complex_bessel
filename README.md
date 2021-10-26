complex_bessel
==============

[![Join the chat at https://gitter.im/valandil/complex_bessel](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/valandil/complex_bessel?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/5354/valandil/complex_bessel.svg)](https://zenodo.org/badge/latestdoi/5354/valandil/complex_bessel)

A C++ library to evaluate Bessel functions of all kinds. More information can  be found on the [website](http://joeydumont.github.io/complex_bessel).

## Introduction

C++ library that acts as a wrapper for the Fortran subroutines developed by D.E. Amos. The library provides functionality to compute the Bessel, Hankel and Airy functions of complex argument and real order. Negative orders are implemented via the standard formulae.

We provide a shared object library and header files to be included.

## Compilation instructions

### Linux

The library uses CMake for compilation. The user should thus install CMake on their machine. On Ubuntu and other Debian-based OSes, this can be done by running
  ```bash
  sudo apt-get install cmake
  ```
On Arch Linux
  ```bash
  sudo pacman -S cmake
  ```

We suggest doing an out-of-tree build by first creating a `build/` folder, then running cmake, i.e.
  ```bash
  cmake -DCMAKE_INSTALL_PREFIX=/path/of/install/dir
 ```

The target of the complex_bessel library is exported
as `complex_bessel::complex_bessel` to a package configuration file for this
library.

After installation you can find this library in your
project `CMakeLists.txt` with:
  ```cmake
  find_package(complex_bessel)
```
or you can just put a copy of `complex_bessel` source code
into your source tree and just add it from the upper level CMakeLists.txt
  ```cmake
 add_subdirectory(complex_bessel)
```
To compile and link this library you should have C++14 and Fortran compilers
installed, and you should enable `CXX Fortran` languages in CMakeLists.txt.
After that to link with your `<target>` you can just 
  ```cmake
  target_link_library(<target> complex_bessel::complex_bessel)
```
Note that CMake will add all additional needed include files to you project compilation automatically.

To run tests you will need to use HDF5, Google Test, Boost 1.6+, and a C compiler (to link with HDF5, so `C` is present in the list of languages in tests/CMakeLists.txt)

###  Windows (experimental support)

There is currently experimental support for compiling with Visual Studio Code. The only tested configuration is with the Intel OneAPI Fortran compiler, installed with Visual Studio integration. The project compiles, but nothing else is tested. Feel free to open a PR for better Windows support.


## Other similar libraries
 
The FORTRAN library that is used as the main driver for the computation of Bessel functions is also used in
  * [`scipy.special`](https://docs.scipy.org/doc/scipy/reference/special.html)
  * MATLAB
   
[Boost](https://www.boost.org/doc/libs/1_72_0/libs/math/doc/html/math_toolkit/bessel/bessel_first.html) has its own implementation of the Bessel functions, but only supports real values for the argument.

If arbitrary precision is needed, the [`arb`](http://arblib.org/acb_hypgeom.html#bessel-functions) library supports the computation of many special functions, including Bessel functions.

## Example
![Contours of Hankel function](/tests/contours.png)
