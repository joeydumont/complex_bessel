/*! \file utilities.h
 *
 * \author Joey Dumont <joey.dumont@gmail.com>
 * \author Denis Gagnon <gagnon88@gmail.com>
 *
 * \brief Definition of some useful constants and functions.
 *
 * We encapsulate some global constants in a namespace 
 * specific to the library and define some functions that
 * are of use in the evaluation of Bessel functions, but are
 * not direclty related to it.
 *
 * \copyright LGPL
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <complex>
#include <cmath>

namespace sp_bessel {

namespace constants {
  const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679; ///< Good old pi.
  const std::complex<double> i = std::complex<double>(0.0,1.0);                                                             ///< Imaginary number.
}

/// We compute the value of cos(pi*nu), paying particular
/// attention to the case where nu is an integer.
inline double cos_pi(double nu)
{
  // Detect if nu is an integer. If |nu|>1e14, the significand is saturated
  // with numbers before the decimal point, and we cannot make the difference between an
  // integer and a real number.
  double nup5 = nu + 0.5;
  if (std::floor(nup5) == nup5 && std::abs(nu) < 1e14)
    return 0.0;

  return std::cos(constants::pi*nu);
}

/// We compute the value of sin(pi*nu), paying particular 
/// attention to the case where nu is an integer.
inline double sin_pi(double nu)
{
  // Detect if nu is an integer. Same comment as above if |nu|>1e14.
  if (std::floor(nu) == nu && std::abs(nu) < 1e14)
    return 0.0;

  return std::sin(constants::pi*nu);
}

} // namespace sp_bessel

#endif // UTILITIES_H
