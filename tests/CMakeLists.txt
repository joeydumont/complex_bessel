# --------------------------------------------------------------------------- #
# Author:       Joey Dumont                  <joey.dumont@gmail.com>          #
# Date created: 2020-02-24                                                    #
# Description:  CMake compilation instructions for complex_bessel tests.      #
#               It is possible to compile tests as a separate project after   #
#               complex_bessel installation into the system.
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
# CMake version
cmake_minimum_required(VERSION 3.1) #update was tested at 3.16
project(complex_bessel_tests C CXX Fortran)
# -- Dependency (Google Test)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GTEST_LIBRARIES})

# -- Required dependency: HDF5.
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set (LIBS ${LIBS} ${HDF5_LIBRARIES})

find_package(Boost 1.60.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
SET (LIBS ${LIBS} ${Boost_LIBRARIES} )

if (NOT TARGET complex_bessel)
    find_package(complex_bessel REQUIRED)
endif()
SET (LIBS ${LIBS} complex_bessel::complex_bessel )

# -- Output tests in directory
add_executable("0011" "0011.cpp")
target_link_libraries("0011" ${LIBS})
add_test(NAME "0011" COMMAND "0011")

# -- Output tests in directory
add_executable("0012" "0012.cpp")
target_link_libraries("0012" ${LIBS})
add_test(NAME "0012" COMMAND "0012")

# -- Output tests in directory
add_executable("0014" "0014.cpp")
target_link_libraries("0014" ${LIBS})
add_test(NAME "0014" COMMAND "0014")
