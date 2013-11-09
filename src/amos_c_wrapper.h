/*! \file amos_c_wrapper.h
 *
 * \author Joey Dumont <joey.dumont@gmail.com>
 * \author Denis Gagnon <gagnon88@gmail.com>
 *
 * \brief Specifies the linkage to the Fortran subroutines
 * for the Bessel functions of complex argument.
 *
 * We use the Fortran code developed by D.E. Amos \cite AMO86
 * to compute the Bessel, Hankel and Airy functions of complex arguments
 * and real order. Only positive arguments are supported in Amos' code,
 * but we implement them using the reflection formulae.
 *
 * \copyright LGPL
 */

#ifndef AMOS_C_WRAPPER_H
#define AMOS_C_WRAPPER_H

extern "C"
{
/*! Bessel function of the first kind. */
extern void zbesj_wrap(double,double,double,int,int,double*,double*,int*,int*);

/*! Bessel function of the second kind. */
extern void zbesy_wrap(double,double,double,int,int,double*,double*,int*, double*,double*,int*);

/*! Modified Bessel function of the first kind. */
extern void zbesi_wrap(double,double,double,int,int,double*,double*,int*,int*);

/*! Modified Bessel function of the second kind. */
extern void zbesk_wrap(double,double,double,int,int,double*,double*,int*,int*);

/*! Hankel function of both kinds. Kind determined by integer argument. */
extern void zbesh_wrap(double,double,double,int,int,int,double*,double*,int*,int*);

/*! Airy function of the first kind. */
extern void zairy_wrap(double,double,int,int,double*,double*,int*,int*);

/*! Airy function of the second kind. */
extern void zbiry_wrap(double,double,int,int,double*,double*,int*);
}

#endif // AMOS_C_WRAPPER_H
