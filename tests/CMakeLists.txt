# --------------------------------------------------------------------------- #
# Author:       Joey Dumont                  <joey.dumont@gmail.com>          #
# Date created: 2020-02-24                                                    #
# Description:  CMake compilation instructions for complex_bessel tests       #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
# CMake version
cmake_minimum_required(VERSION 3.1)

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

# -- Output tests in directory
add_executable("0011" "0011.cpp")
include_directories(${CMAKE_SOURCE_DIR}/include)
target_link_libraries("0011" ${PROJECT_NAME})
target_link_libraries("0011" ${LIBS})
add_test(NAME "0011" COMMAND "0011")
