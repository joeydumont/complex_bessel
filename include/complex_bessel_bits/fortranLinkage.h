/*! \file fortranLinkage.h
 *
 * \author Joey Dumont <joey.dumont@gmail.com>
 * \author Denis Gagnon <gagnon88@gmail.com>
 *
 * \brief Specifies the linkage to the Fortran subroutines
 * for the Bessel functions of complex argument.
 *
 * We use the Fortran code developed by D.E. Amos \cite AMO86
 * to compute the Bessel, Hankel and Airy functions of complex arguments
 * and integer order. Only positive arguments are supported in Amos' code,
 * but we implement them using the reflection formulae.
 *
 * For real arguments of the Bessel functions, we recommend
 * using the Boost libraries, which are faster for double
 * return value and argument.
 *
 * \copyright LGPL
 */

#ifndef FORTRANLINKAGE_H
#define FORTRANLINKAGE_H

namespace sp_bessel {
/*! @name Fortran Linkage
 * We link the Fortran subroutines to our C++ code.
 * Since Fortran treats all variables by reference, we must
 * pass our C++ arguments as pointers.
 *
 * Most subroutines have flags to signal errors such as underflow
 * and overflow. Exponential scaling is also available.
 * The KODE flag is a parameter that indicates the scaling option.
 * KODE = 1 means no scaling while KODE activates the exponential scaling.
 * Integer N controls the size of the array that is returned in the Fortran code.
 * We use N=1 throughout as to return scalars.  The NZ integer counts the number of
 * elements set to zero due to underflow. We ignore it. The IERR integer
 * is for error signaling.
 */
///@{
extern "C"
{
  /*! Bessel function of the first kind. */
  extern void zbesj_wrap(double, double, double, int, int, double*, double*, int*, int*);

  /*! Bessel function of the second kind. */
  extern void
  zbesy_wrap(double, double, double, int, int, double*, double*, int*, double*, double*, int*);

  /*! Modified Bessel function of the first kind. */
  extern void zbesi_wrap(double, double, double, int, int, double*, double*, int*, int*);

  /*! Modified Bessel function of the second kind. */
  extern void zbesk_wrap(double, double, double, int, int, double*, double*, int*, int*);

  /*! Hankel function of both kinds. Kind determined by integer argument. */
  extern void zbesh_wrap(double, double, double, int, int, int, double*, double*, int*, int*);

  /*! Airy function of the first kind. */
  extern void zairy_wrap(double, double, int, int, double*, double*, int*, int*);

  /*! Airy function of the second kind. */
  extern void zbiry_wrap(double, double, int, int, double*, double*, int*);
}
///@}

} // namespace sp_bessel

#endif // FORTRANLINKAGE_H
