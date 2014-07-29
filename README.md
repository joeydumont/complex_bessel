complex_bessel
==============
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.11077.png)](http://dx.doi.org/10.5281/zenodo.11077)

A C++ library to evaluate Bessel functions of all kinds. More information can be found on the [website](http://valandil.github.io/complex_bessel).

## Introduction

C++ library that acts as a wrapper for the Fortran subroutines developed by D.E. Amos. The library provides functionality to compute the Bessel, Hankel and Airy functions of complex argument and nonnegative order. Negative orders are implemented via the standard formulae.

We provide a shared object library and header files to be included.

## Compilation instructions

Since the library uses automake to generate a custom Makefile, the user must have, along with basic build tools, it installed on their system. *NIX systems usually come pre-installed with those packages. Otherwise, they should be readily available using your preferred packaging tool.

The user should then run
  ```bash
  bash autogen.sh
  ```

which uses the configure.ac file to generate a configure script. Then, running 
```bash
./configure && make && sudo make install
```

should do the trick. The library will be installed to `/usr/local/lib` by default, so make sure that your compiler looks for libraries in that folder.

## Example
![Contours of Hankel function](/tests/contours.png)
