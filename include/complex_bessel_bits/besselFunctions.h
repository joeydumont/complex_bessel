/*! \file besselFunctions.h
 *
 * \author Joey Dumont <joey.dumont@gmail.com>
 * \author Denis Gagnon <gagnon88@gmail.com>
 *
 * \brief Definition of the functions computing the Bessel functions.
 *
 * Using the Fortran subroutines, linked to the C language in the
 * file fortranLinkage.h, we define the functions that will be used
 * to evaluate the Bessel functions and their derivatives.
 *
 * \copyright LGPL
 */

#ifndef BESSELFUNCTIONS_H
#define BESSELFUNCTIONS_H

#include <iostream>
#include <vector>

#include "errors.h"
#include "fortranLinkage.h"
#include "utilities.h"

using namespace std::complex_literals;

namespace sp_bessel {

/*! @name Evaluation of Bessel functions.
 * We implement Amos' Fortran subroutines in C++.
 */

///@{

/*! Using a function as a template parameter, we define a function that
 * computes the derivative of the Bessel functions \f$J,\,Y,\,H^{(1,2)},\,I,\,K\f$
 * using the recurrence relations \cite ABR65 (Sects. 9.1.27/9.6.26).
 * We have to compute with scale=0, otherwise the derivative is not valid.
 */
template<std::complex<double> (*T)(double, std::complex<double>, bool, BesselErrors*)>
inline std::complex<double>
diffBessel(
  double                     order,
  std::complex<double>       z,
  int                        n,
  double                     phase,
  std::vector<BesselErrors>* errors = nullptr)
{
  // For J, Y, H1 and H2, phase = -1.
  // For I, e^(order*pi*i)K, phase = 1.
  // First term of the series.
  double p = 1.0;

  // Error handling.
  if (errors) {
    (*errors).push_back(BesselErrors());
  }
  std::complex<double> s = T(order - n, z, false, errors ? &(*errors).back() : nullptr);

  // Rest of the series
  for (int i = 1; i <= n; i++) {
    p = phase * (p * (n - i + 1)) / i; // = choose(n,k).

    if (errors) {
      (*errors).push_back(BesselErrors());
    }

    s += p * T(order - n + 2 * i, z, false, errors ? &(*errors).back() : nullptr);
  }

  return s / std::pow(2.0, n);
}

/*! Computes the Bessel functions of the first kind with the reflection formula
 * \f$J_{-\nu}(z) = (-1)^\nu J_\nu(z)\f$. \cite ABR65 Sec. 9.1.5. */
inline std::complex<double>
besselJ(double order, std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  // Input values for Fortran subroutines.
  double zr   = std::real(z);
  double zi   = std::imag(z);
  double nu   = std::abs(order);
  int    kode = scale + 1;
  int    N    = 1;

  // Output values.
  double cyr, cyi;
  int    nz, ierr;

  // External function call
  zbesj_wrap(zr, zi, nu, kode, N, &cyr, &cyi, &nz, &ierr); // Call Fortran subroutine.

  // If the input is real, then the output should be real as well.
  if (zi == 0.0 && zr >= 0.0)
    cyi = 0.0;
  std::complex<double> answer(cyr, cyi); // Placeholder for output.

  // If order is negative, then we must apply the reflection formula.
  if (order < 0.0) {
    // We prepare the rotation coefficients.
    double c = cos_pi(nu);
    double s = sin_pi(nu);

    // We compute the Bessel Y function.
    double cyrY, cyiY, cwrkr, cwrki;
    int    nzY, ierrY, kodeY(1), NY(1);

    // External function call
    zbesy_wrap(zr, zi, nu, kodeY, NY, &cyrY, &cyiY, &nzY, &cwrkr, &cwrki, &ierrY);
    std::complex<double> answerY(cyrY, cyiY);

    answer = c * answer - s * answerY;
  }

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("besselJ").at(ierr);
  }

  return answer;
}

/*! Computes the nth derivative of besselJ. */
inline std::complex<double>
besselJp(double order, std::complex<double> z, int n = 1, std::vector<BesselErrors>* errors = nullptr)
{
  return diffBessel<besselJ>(order, z, n, -1, errors);
}

/*! Computes the Bessel function of the second kind with the reflection formula
 * \f$Y_{-\nu}(z) = (-1)^\nu Y_\nu(z)\f$ \cite ABR65 Sec 9.1.5. */
inline std::complex<double>
besselY(double order, std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  // Input values for Fortran subroutines
  double zr   = std::real(z);
  double zi   = std::imag(z);
  double nu   = std::abs(order);
  int    kode = scale + 1;
  int    N    = 1;

  // Output and temporary varibles
  double cyr, cyi, cwrkr, cwrki;
  int    nz, ierr;

  // External function call
  zbesy_wrap(zr, zi, nu, kode, N, &cyr, &cyi, &nz, &cwrkr, &cwrki, &ierr); // Call Fortran subroutine.

  // In passing from C++ to FORTRAN, the exact zero becomes the numerical zero (10^(-14)).
  // The limiting form of Y_nu(z) for high order, -Gamma(nu)/pi*Re(z)^(-nu)*(1-i*nu*Im(z)/Re(z)),
  // leads to product of the form zero*infinity, which destroys numerical precision. We hence
  // manually set the imaginary part of the answer to zero is the imaginary part of the input
  // is zero.
  if (zi == 0.0 && zr >= 0.0)
    cyi = 0.0;
  std::complex<double> answer(cyr, cyi);

  // If order is negative, we must apply the reflection formula.
  if (order < 0.0) {
    // We prepare the rotation coefficients.
    double c = cos_pi(nu);
    double s = sin_pi(nu);

    // We compute the Bessel J function.
    double cyrJ, cyiJ;
    int    nzJ, ierrJ, kodeJ(1), NJ(1);

    zbesj_wrap(zr, zi, nu, kodeJ, NJ, &cyrJ, &cyiJ, &nzJ, &ierrJ);
    std::complex<double> answerJ(cyrJ, cyiJ);
    answer = s * answerJ + c * answer;
  }

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("besselY").at(ierr);
  }

  return answer;
}

/*! Computes the nth derivative of besselY. */
inline std::complex<double>
besselYp(double order, std::complex<double> z, int n = 1, std::vector<BesselErrors>* errors = nullptr)
{
  return diffBessel<besselY>(order, z, n, -1, errors);
}

/*! Computes the modified Bessel function of the first kind. Negative
 *  orders are equal to the positive ones if \nu is an integer.: \f$I_{-n}(z)=I_{n}(z)\f$.
 *  Otherwise, we have to exploit the relationship:
 *  \f$K_\nu(z) = \frac{\pi}{2}\frac{I_{-\nu}(z)-I_\nu(z)}{\sin\nu\pi}\f$
 *  \cite ABR65 \S9.6.2.
 */
inline std::complex<double>
besselI(double order, std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  // Input values for Fortran subroutines.
  double zr   = std::real(z);
  double zi   = std::imag(z);
  double nu   = std::abs(order);
  int    kode = scale + 1;
  int    N    = 1;

  // Output and temporary variables.
  double cyr, cyi;
  int    nz, ierr;

  // External function call.
  zbesi_wrap(zr, zi, nu, kode, N, &cyr, &cyi, &nz, &ierr); // Call Fortran subroutine.

  // Enforcing some conditions on the output as a function of the output.,
  if (zi == 0.0 && zr >= 0.0)
    cyi = 0.0;
  std::complex<double> answer(cyr, cyi);

  // We apply the reflection formula is order is negative.
  if (order < 0.0 && scale) {
    throw("The scale cannot be applied to a negative order, as the reflection formula involves "
          "besselK, which does not share the same scaling behaviour.");
  }

  if (order < 0.0) {
    // We prepare the reflection coefficients.
    double s = sin_pi(nu);

    // We evaluate the besselK function.
    int    nzK, ierrK;
    double cyrK, cyiK;

    // External function call.
    zbesk_wrap(zr, zi, nu, kode, N, &cyrK, &cyiK, &nzK, &ierrK);

    if (zi == 0.0 && zr >= 0.0)
      cyiK = 0.0;
    std::complex<double> answerK(cyrK, cyiK);
    answer += 2.0 / constants::pi * s * answerK;
  }

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("besselI").at(ierr);
  }

  return answer;
}

/*! Computes the nth derivative of besselI. */
inline std::complex<double>
besselIp(double order, std::complex<double> z, int n = 1, std::vector<BesselErrors>* errors = nullptr)
{
  return diffBessel<besselI>(order, z, n, 1, errors);
}

/*! Computes the modified Bessel function of the second kind. Negative
 *  orders are equal to the positive ones: \f$K_{-nu}(z)=K_{nu}(z)\f$. */
inline std::complex<double>
besselK(double order, std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  // Input values for Fortran subroutines.
  double zr   = std::real(z);
  double zi   = std::imag(z);
  double nu   = std::abs(order);
  int    kode = scale + 1;
  int    N    = 1;

  // Output and temporary variables.
  double cyr, cyi;
  int    nz, ierr;

  // External function call.
  zbesk_wrap(zr, zi, nu, kode, N, &cyr, &cyi, &nz, &ierr); // Call Fortran subroutine.

  // In passing from C++ to FORTRAN, the exact zero becomes the numerical zero (10^(-14)).
  // The limiting form of K_nu(z) for high order, Gamma(nu)/2*(z/2)^(-nu),
  // leads to product of the form zero*infinity for the imaginary part, which destroys numerical
  // precision. We hence manually set the imaginary part of the answer to zero is the imaginary part of
  // the input is zero.
  if (zi == 0.0 && zr >= 0.0)
    cyi = 0.0;
  std::complex<double> answer(cyr, cyi);

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("besselK").at(ierr);
  }

  return answer;
}

/*! This function is used only when computing the derivative of the BesselK function.
 *  It cannot be used with scale = true, but we need to provide it to obey the signature
 *  of diffBessel. We thus set it directly in the function call without using the argument.
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
inline static std::complex<double>
expBesselK(double order, std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  return std::exp(order * constants::pi * constants::i) * besselK(order, z, false, error);
}
#pragma GCC diagnostic pop

/*! Computes the nth derivative of besselK. */
inline std::complex<double>
besselKp(double order, std::complex<double> z, int n = 1, std::vector<BesselErrors>* errors = nullptr)
{
  return std::exp(-order * constants::pi * constants::i) *
         diffBessel<expBesselK>(order, z, n, 1, errors);
}

/*! Computes the Hankel function of the first kind. We also implement
 * the reflection formula \f$H^{(1)}_{-\nu}(z) = H^{(1)}_\nu(z)\exp\left(\pi\nu\imath\right)
 * \f$. */
inline std::complex<double>
hankelH1(double order, std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  // Input values.
  double zr   = std::real(z);
  double zi   = std::imag(z);
  double nu   = std::abs(order);
  int    kode = scale + 1;
  int    N    = 1;
  int    kind = 1;

  // Output values
  double cyr, cyi;
  int    nz, ierr;

  // External function call.
  zbesh_wrap(zr, zi, nu, kode, kind, N, &cyr, &cyi, &nz, &ierr);
  std::complex<double> answer(cyr, cyi);

  // Reflection formula if order is negative.
  if (order < 0.0)
    answer *= std::exp(constants::pi * nu * constants::i);

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("hankelH").at(ierr);
  }

  return answer;
}

/*! Computes the nth derivative of hankelH1. */
inline std::complex<double>
hankelH1p(double order, std::complex<double> z, int n = 1, std::vector<BesselErrors>* errors = nullptr)
{
  return diffBessel<hankelH1>(order, z, n, -1, errors);
}

/*! Computes the Hankel function of the second kind. We also implement the reflection
 * formula \f$H^{(1)}_{-\nu}(z) = H^{(1)}_\nu(z)\exp\left(-\pi\nu\imath\right)\f$. */
inline std::complex<double>
hankelH2(double order, std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  //    // Input values.
  double zr   = std::real(z);
  double zi   = std::imag(z);
  double nu   = std::abs(order);
  int    kode = scale + 1;
  int    N    = 1;
  int    kind = 2;

  // Output values
  double cyr, cyi;
  int    nz, ierr;

  // External function call.
  zbesh_wrap(zr, zi, nu, kode, kind, N, &cyr, &cyi, &nz, &ierr);
  std::complex<double> answer(cyr, cyi);

  // Reflection formula if order is negative.
  if (order < 0.0)
    answer *= std::exp(-constants::pi * nu * constants::i);

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("hankelH").at(ierr);
  }

  return answer;
}

/*! Computes the nth derivative of hankelH2.*/
inline std::complex<double>
hankelH2p(double order, std::complex<double> z, int n = 1, std::vector<BesselErrors>* errors = nullptr)
{
  return diffBessel<hankelH2>(order, z, n, -1, errors);
}

/*! Computes the complex Airy Ai(z) function. */
inline std::complex<double>
airy(std::complex<double> z, int id = 0, bool scale = false, BesselErrors* error = nullptr)
{
  // Input values.
  double zr   = std::real(z);
  double zi   = std::imag(z);
  int    kode = scale + 1;

  // Output values.
  double air, aii;
  int    nz, ierr;

  // External function call.
  zairy_wrap(zr, zi, id, kode, &air, &aii, &nz, &ierr);
  std::complex<double> answer(air, aii);

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("airy").at(ierr);
  }

  return answer;
}

/*! Computes the first derivative of airy. */
inline std::complex<double>
airyp(std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  return airy(z, 1, scale, error);
}

// Computes the complex Airy funciton Bi(z). */
inline std::complex<double>
biry(std::complex<double> z, int id = 0, bool scale = false, BesselErrors* error = nullptr)
{
  // Input values.
  double zr   = std::real(z);
  double zi   = std::imag(z);
  int    kode = scale + 1;

  // Output values.
  double bir, bii;
  int    ierr;

  // External function call.
  zbiry_wrap(zr, zi, id, kode, &bir, &bii, &ierr);
  std::complex<double> answer(bir, bii);

  std::cout << ierr;

  // If an error struct is provided, set the error number and message.
  if (error) {
    *error = errorMessages.at("biry").at(ierr);
  }

  return answer;
}

/*! Computes the first derivative of biry. */
inline std::complex<double>
biryp(std::complex<double> z, bool scale = false, BesselErrors* error = nullptr)
{
  return biry(z, 1, scale, error);
}

///@}

} // namespace sp_bessel

#endif // BESSELFUNCTIONS_H
