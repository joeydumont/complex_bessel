# Autoconf flags
ACLOCAL_AMFLAGS = -I m4

# Name of library 
lib_LTLIBRARIES = libcomplex_bessel.la

# List of needed source files
libcomplex_bessel_la_SOURCES = src/machine.for src/zbesh.for src/amos_iso_c_fortran_wrapper.f90
libcomplex_bessel_la_LDFLAGS = -release 0.4.0

# List of needed headers
library_includedir=$(includedir)/complex_bessel
library_include_HEADERS = src/fortranLinkage.h  src/utilities.h src/besselFunctions.h

library2_includedir=$(includedir)
library2_include_HEADERS = complex_bessel.h
