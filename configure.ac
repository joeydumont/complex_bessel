# Process this file with autoconf to produce a configure script.
AC_INIT(complex_bessel, 0.1, gagnon88@gmail.com)
AM_INIT_AUTOMAKE([1.10 -Wall no-define subdir-objects])

# Required autoconf version
AC_PREREQ([2.59])

# Configuration header to generate
AC_CONFIG_HEADERS([config.h])

# Required programs
AC_PROG_CXX
AC_PROG_CC

# Check installation of Fortran 77 (gfortran is preferable)
AC_PROG_F77([gfortran f95 fl32 f77 fort77 xlf g77 f90 xlf90])
AC_PROG_FC

# Enable shared library
AC_PROG_INSTALL
AC_ENABLE_SHARED
AC_DISABLE_STATIC

# Check for archiver
AM_PROG_AR

# Check for libtool
LT_INIT

# Configuration files
AC_CONFIG_FILES([Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_OUTPUT
