# Autoconf flags
ACLOCAL_AMFLAGS = -I m4

# Name of library
lib_LTLIBRARIES = libcomplex_bessel.la

# List of needed source files
libcomplex_bessel_la_SOURCES = 	\
	src/zbesh.for 				\
	src/machine.for 			\
	src/amos_c_wrapper.f90 		

libcomplex_bessel_la_LDFLAGS = -version-info 0:0:0

# List of needed headers
library_includedir=$(includedir)/complex_bessel
library_include_HEADERS = \ 
	src/amos_c_wrapper.h

