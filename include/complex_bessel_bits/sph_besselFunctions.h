/*! \file sph_besselFunctions.h
 *
 * \author Joey Dumont <joey.dumont@gmail.com>
 * \author Denis Gagnon <gagnon88@gmail.com>
 *
 * \brief Functions computing the spherical Bessel functions.
 *
 * We compute the spherical Bessel functions by evaluating the
 * Bessel functions with half-integer order with D.E. Amos' FORTRAN
 * implementation. 
 *
 * \copyright LGPL
 */

#ifndef SPH_BESSELFUNCTIONS_H
#define SPH_BESSELFUNCTIONS_H

#include <iostream>
#include <complex>

#include "besselFunctions.h"
#include "utilities.h"

namespace sp_bessel {

/// Near \f$z\rightarrow0\f$, the floating point division necessary would
/// destroy precision. We thus use an ascending series to compute sph_besselJ.
/// See \cite ABR65 Sec. 10.1.2.
inline std::complex<double> sph_besselJ_asc_series(double order, std::complex<double> z)
{
  return 0;
}

/// We compute the spherical Bessel function of the first kind.
inline std::complex<double> sph_besselJ(double order, std::complex<double> z)
{
  return std::sqrt(constants::pi/(2.0*z))*besselJ(order+0.5, z);
}

/// We compute the spherical Bessel function of the second kind.
inline std::complex<double> sph_besselY(double order, std::complex<double> z)
{
  return std::sqrt(constants::pi/(2.0*z))*besselY(order+0.5,z);
}

/// We compute the spherical Hankel function of the first kind.
inline std::complex<double> sph_hankelH1(double order, std::complex<double> z)
{
  return std::sqrt(constants::pi/(2.0*z))*hankelH1(order+0.5,z);
}

/// We compute the spherical Hankel function of the second kind.
inline std::complex<double> sph_hankelH2(double order, std::complex<double> z)
{
  return std::sqrt(constants::pi/(2.0*z))*hankelH2(order+0.5,z);
}

} // namespace sp_bessel

#endif // SPH_BESSELFUNCTIONS_H