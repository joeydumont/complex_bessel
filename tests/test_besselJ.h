#ifndef TESTS_BESSELJ_H
#define TESTS_BESSELJ_H

/*! \file tests_besselJ.h
 *
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *
 *  \date 2014-01-14
 *
 *  \brief Defines functions that test the accuracy of our Bessel functions of the first kind.
 *
 * We define a number of functions designed to tests the numerical implementation
 * of the evaluation of Bessel functions. We will try to refer to the source of the
 * tests as much as possible. We will use the following abbreviations:
 *  - AS: Abromowitz & Stegun;
 *
 * \copyright LGPL
 */

#include <armadillo>
#include <complex>
#include <complex_bessel.h>

using namespace arma;
using namespace sp_bessel;

/*! AS 9.1.76. Addition theorem for Bessel functions. */
std::complex<double>
besselJAddition1(std::complex<double> z, int kMax);

/*! AS 9.1.77. Addition theorem for Bessel functions. */
std::complex<double>
besselJAddition2(std::complex<double> z, int N, int kMax);

/*! AS 9.1.15. Wronskian for two Bessel J functions. */
std::complex<double>
besselJWronskian1(int order, std::complex<double> z);

#endif // TESTS_BESSELJ_H
