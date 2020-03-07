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

# ----------------------------------------------------------------- #
# --                    Compiler Configuration                   -- #
# ----------------------------------------------------------------- #

# -- Compiler config (flags and such)
enable_language (Fortran)

# -- Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# C++14 Standard required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configuration for the GCC compiler.
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3 -Wall -Wextra -pedantic -march=native -mtune=native")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3 -Wall -Wextra -pedantic -march=native -mtune=native")

  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -pg -g -Wall -DNDEBUG")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -pg -g -Wall -DNDEBUG")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   -Wall -pg -O3 -DNDEBUG")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wall -pg -O3 -DNDEBUG")

# Configuration for the Intel compiler.
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -O3")
  set (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}    -O3")

  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -debug all")
  set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0 -g -debug all")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   -O3 -g -debug all")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -g -debug all")
endif()

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

add_library(${PROJECT_NAME}_static STATIC
  ${SRC_LIST}
  "${CMAKE_CURRENT_SOURCE_DIR}/src/amos_iso_c_fortran_wrapper.f90"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/zbesh.for"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/machine.for")

# Link dependencies
target_link_libraries (${PROJECT_NAME} ${LIBS})
target_link_libraries (${PROJECT_NAME}_static ${LIBS})

# Install directories
install (TARGETS ${PROJECT_NAME} DESTINATION lib)
install (DIRECTORY include/ DESTINATION include)

include(CTest)
add_subdirectory(tests)
