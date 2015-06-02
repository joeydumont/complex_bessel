# --------------------------------------------------------------------------- #
# Author:       Joey Dumont                   <joey.dumont@gmail.com>         #
# Author:       Denis Gagnon                  <gagnon88@gmail.com>            #
# Date:         2015-06-01                                                    #
# Description:  CMake compilation instructions.                               #
# ----------------------------------------------------------------------------#

# ----------------------------------------------------------------- #
# --                 Name and version of library                 -- #
# ----------------------------------------------------------------- #
project(complex_bessel)
set (complex_bessel_VERSION_MAJOR 0)
set (complex_bessel_VERSION_MINOR 6)
set (complex_bessel_VERSION_RELEASE 0)

# ----------------------------------------------------------------- #
# --               Configuration and Dependencies                -- #
# ----------------------------------------------------------------- #
# -- CMake version and installation directory.
cmake_minimum_required(VERSION 2.8)
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX /usr)
endif()

# -- Compiler config (flags and such)
enable_language (Fortran)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -march=native")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -march=native -std=c++11")

# ----------------------------------------------------------------- #
# --                  Compilation Instructions                   -- #
# ----------------------------------------------------------------- #
# -- Included files
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/include)

# -- Source files
aux_source_directory (${CMAKE_CURRENT_SOURCE_DIR}/src SRC_LIST)

# -- Build a shared library
add_library(${PROJECT_NAME} SHARED 
  ${SRC_LIST}
  "${CMAKE_CURRENT_SOURCE_DIR}/src/amos_iso_c_fortran_wrapper.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/zbesh.for"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/machine.for")

# Define the properties of the library to
# be built (version and other stuff).
set_target_properties (
  ${PROJECT_NAME}
  PROPERTIES
    VERSION   ${complex_bessel_VERSION_MAJOR}.${complex_bessel_VERSION_MINOR}.${complex_bessel_VERSION_RELEASE}
    SOVERSION ${complex_bessel_VERSION_MAJOR}.${complex_bessel_VERSION_MINOR}.${complex_bessel_VERSION_RELEASE}
)

# Link dependencies
target_link_libraries (${PROJECT_NAME} ${LIBS})

# Install directories
install (TARGETS ${PROJECT_NAME} DESTINATION lib)
install (DIRECTORY include/ DESTINATION include)

# ----------------------------------------------------------------- #
# --                   Generate Documentation                    -- #
# ----------------------------------------------------------------- #
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(
  doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)