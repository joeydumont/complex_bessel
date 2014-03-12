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

#include <complex>
#include <iostream>

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
extern void zbiry_wrap(double,double,int,int,double*,double*,int*,int*);
}
///@}

/*! @name Evaluation of Bessel functions.
 * We implement Amos' Fortran subroutines in C++.
 * \todo Provide error detection and signaling.  */
///@{
/*! Using a function as a template parameter, we define a function that
 * computes the derivative of the Bessel funcitons \f$J,\,Y,\,H^{(1,2)},\,I,\,K\f$
 * using the recurrence relations \cite ABR65 (Sects. 9.1.27/9.6.26). */
template <std::complex<double> (*T)(int, std::complex<double>)>
inline std::complex<double> diffBessel(int order, std::complex<double> z, int signs = 0)
{
	double sgn1(1.),sgn2(-1.);
	switch (signs)
	{
		case 1:
		// Bessel function I. 
		{
			sgn2 = 1.;
			break;
		}
		case 2:
		// Bessel function K.
		{
			sgn1 = -1.;
			sgn2 = -1.;
			break;
		}
		default:
		{
			break;
		}
	}

    return 0.5*(sgn1*T(order-1.0,z)+sgn2*T(order+1.0,z));
}

/*! Computes the Bessel functions of the first kind with the reflection formula
 * \f$J_{-\nu}(z) = (-1)^\nu J_\nu(z)\f$. \cite ABR65 Sec. 9.1.5. */
inline std::complex<double> besselJ(int order, std::complex<double> z)
{
    // Input values for Fortran subroutines.
    double zr = std::real(z);
    double zi = std::imag(z);
    double nu = std::fabs((double) order);
    int kode = 1;
    int N = 1;

    // Output values.
    double cyr,cyi;
    int nz,ierr;

    // External function call
    zbesj_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&ierr);    // Call Fortran subroutine.
    std::complex<double> answer(cyr,cyi);               // Placeholder for output.

    // If order is negative, then we must apply the reflection formula.
    if (order < 0)
    {
        answer *= (order & 1 ? -1.0 : 1.0);
    }

    // If the return code is not normal, we print the error code.
    if (ierr!=0) std::cout << "besselJ: Error code " << ierr << "." << std::endl;

    return answer;
}

/*! Computes the first derivative of besselJ. */
inline std::complex<double> besselJp(int order, std::complex<double> z)
{
	return diffBessel<besselJ>(order, z);
}

/*! Computes the Bessel function of the second kind with the reflection formula
 * \f$Y_{-\nu}(z) = (-1)^\nu Y_\nu(z)\f$ \cite ABR65 Sec 9.1.5. */
inline std::complex<double> besselY(int order, std::complex<double> z)
{
    // Input values for Fortran subroutines
    double zr = std::real(z);
    double zi = std::imag(z);
    double nu = std::fabs((double) order);
    int kode = 1;
    int N = 1;

    // Output and temporary varibles
    double cyr, cyi,cwrkr,cwrki;
    int nz, ierr;

    // External function call
    zbesy_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&cwrkr,&cwrki,&ierr); // Call Fortran subroutine.
    
    // In passing from C++ to FORTRAN, the exact zero becomes the numerical zero (10^(-14)). 
    // The limiting form of Y_nu(z) for high order, -Gamma(nu)/pi*Re(z)^(-nu)*(1-i*nu*Im(z)/Re(z)),
    // leads to product of the form zero*infinity, which destroys numerical precision. We hence
    // manually set the imaginary part of the answer to zero is the imaginary part of the input
    // is zero. 
    if (zi == 0.0) cyi=0.0;
    std::complex<double> answer(cyr,cyi);                           // Placeholder for output

    // If order is negative, we must apply the reflection formula.
    if (order < 0)
    {
        answer *= (order & 1 ? -1.0 : 1.0);
    }

    // If the return code is not normal, we print the error code.
    if (ierr!=0) std::cout << "besselY: Error code " << ierr << "." << std::endl;

    return answer;
}

/*! Computes the first derivative of besselY. */
inline std::complex<double> besselYp(int order, std::complex<double> z)
{
	return diffBessel<besselY>(order, z);
}

/*! Computes the modified Bessel function of the first kind. Negative
 *  orders are equal to the positive ones: \f$I_{-nu}(z)=I_{nu}(z)\f$. */
inline std::complex<double> besselI(int order, std::complex<double> z)
{
	// Input values for Fortran subroutines. 
	double zr = std::real(z);
	double zi = std::imag(z);
	double nu = std::fabs((double)order);
	int kode = 1;
	int N = 1;

	// Output and temporary variables.
	double cyr,cyi;
	int nz, ierr;

	// External function call. 
	zbesi_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&ierr); // Call Fortran subroutine.
	std::complex<double> answer(cyr,cyi);

	// In case of error, we print the error code. 
	if (ierr!=0) std::cout << "besselI: Error code " << ierr << "." << std::endl;

	return answer;
}

/*! Computes the first derivative of besselI. */
inline std::complex<double> besselIp(int order, std::complex<double> z)
{
	return diffBessel<besselI>(order, z, 1);
}

/*! Computes the modified Bessel function of the second kind. Negative
 *  orders are equal to the positive ones: \f$K_{-nu}(z)=K_{nu}(z)\f$. */
inline std::complex<double> besselK(int order, std::complex<double> z)
{
	// Input values for Fortran subroutines.
	double zr = std::real(z);
	double zi = std::imag(z);
	double nu = std::fabs((double) order);
	int kode = 1;
	int N = 1;

	// Output and temporary variables.
	double cyr, cyi;
	int nz, ierr;

	// External function call. 
	zbesk_wrap(zr,zi,nu,kode,N,&cyr,&cyi,&nz,&ierr); // Call Fortran subroutine.

	// In passing from C++ to FORTRAN, the exact zero becomes the numerical zero (10^(-14)). 
    // The limiting form of K_nu(z) for high order, Gamma(nu)/2*(z/2)^(-nu),
    // leads to product of the form zero*infinity for the imaginary part, which destroys numerical precision. We hence
    // manually set the imaginary part of the answer to zero is the imaginary part of the input
    // is zero.
    if (zi == 0.0) cyi = 0.0;
    std::complex<double> answer(cyr,cyi); 

    // In case of error, we print the error code.
    if (ierr!=0) std::cout << "besselK: Error code " << ierr << "." << std::endl;

    return answer;
}

/*! Computes the first derivative of besselK. */
inline std::complex<double> besselKp(int order, std::complex<double> z)
{
	return diffBessel<besselK>(order, z, 2);
}

/*! Computes the Hankel function of the first kind. We also implement
 * the reflection formula \f$H^{(1)}_{-\nu}(z) = H^{(1)}_\nu(z)\exp\left(\pi\nu\imath\right)
 * \f$. */
inline std::complex<double> besselH1(int order, std::complex<double> z)
{
    // Input values.
//    double zr = std::real(z);
//    double zi = std::imag(z);
//    double nu = std::fabs((double) order);
//    int kode = 1;
//    int N = 1;
//    int kind = 1;

//    // Output values
//    double cyr,cyi;
//    int nz,ierr;

//    // External function call.
//    zbesh_wrap(zr,zi,nu,kode,kind,N,&cyr,&cyi,&nz,&ierr);
//    std::complex<double> answer(cyr,cyi);

//    // Reflection formula if order is negative.
//    if (order < 0.0)
//    {
//        // Compute complex exponential.
//        std::complex<double> i(0.0,1.0);
//        i = std::exp(arma::datum::pi*nu*i);
//        answer *= i;
//    }

    std::complex<double> i(0.0,1.0);
    return besselJ(order,z)+i*besselY(order,z);
}

/*! Computes the first derivative of besselH1. */
inline std::complex<double> besselH1p(int order, std::complex<double> z)
{
	return diffBessel<besselH1>(order, z);
}

/*! Computes the Hankel function of the second kind. We also implement the reflection
 * formula \f$H^{(1)}_{-\nu}(z) = H^{(1)}_\nu(z)\exp\left(-\pi\nu\imath\right)\f$. */
inline std::complex<double> besselH2(int order, std::complex<double> z)
{
//    // Input values.
//    double zr = std::real(z);
//    double zi = std::imag(z);
//    double nu = std::fabs((double) order);
//    int kode = 1;
//    int N = 1;
//    int kind = 2;

//    // Output values
//    double cyr,cyi;
//    int nz,ierr;

//    // External function call.
//    zbesh_wrap(zr,zi,nu,kode,kind,N,&cyr,&cyi,&nz,&ierr);
//    std::complex<double> answer(cyr,cyi);

//    // Reflection formula if order is negative.
//    if (order < 0.0 )
//    {
//        std::complex<double> i(0.0,1.0);
//        i = std::exp(-arma::datum::pi*nu*i);
//        answer *= i;
//    }

    std::complex<double> i(0.0,1.0);
    return besselJ(order,z)-i*besselY(order,z);
}

/*! Computes the first derivative of besselH2.*/
inline std::complex<double> besselH2p(int order, std::complex<double> z)
{
	return diffBessel<besselH2>(order, z);
}


/*! Computes the complex Airy Ai(z) function. */
inline std::complex<double> airy(std::complex<double> z, int id = 0)
{
	// Input values.
	double zr = std::real(z);
	double zi = std::imag(z);
	int kode = 1;

	// Output values.
	double air, aii;
	int nz, ierr;

	// External function call. 
	zairy_wrap(zr,zi,id,kode,&air,&aii,&nz,&ierr);
	std::complex<double> answer(air,aii);

	return answer;
}

/*! Computes the first derivative of airy. */
inline std::complex<double> airyp(std::complex<double> z)
{
	return airy(z,1);
}

// Computes the complex Airy funciton Bi(z). */
inline std::complex<double> biry(std::complex<double> z, int id = 0)
{
	// Input values.
	double zr = std::real(z);
	double zi = std::imag(z);
	int kode = 1;

	// Output values.
	double bir, bii;
	int nz,ierr;

	// External function call. 
	zbiry_wrap(zr,zi,id,kode,&bir,&bii,&nz,&ierr);
	std::complex<double> answer(bir,bii);

	return answer;
}

/*! Computes the first derivative of biry. */
inline std::complex<double> biryp(std::complex<double> z)
{
	return biry(z,1);
}
///@}
}

#endif // FORTRANLINKAGE_H