# --------------------------------------------------------------------------- #
# Author:       Joey Dumont                   <joey.dumont@gmail.com>         #
# Author:       Denis Gagnon                  <gagnon88@gmail.com>            #
# Date:         2015-06-01                                                    #
# Description:  CMake compilation instructions.                               #
# Author:       Konstantin Ladutenko          <kostyfisik@gmail.com>          #
# Date:         2021-03-24                                                    #
# Description:  CMake3 update, using cmake config from                        #
#               https://github.com/szaghi/PENF as a reference.                #
# ----------------------------------------------------------------------------#

# The target of the complex_bessel library is exported
# as complex_bessel::complex_bessel to a package configuration file for this
# library
#
# usage:
#     find_package(complex_bessel) or add_subdirectory(complex_bessel)
# and then for your <target>:
#     target_link_library(<target> complex_bessel::complex_bessel)
#
# the config file is generated in the build and install directories
#
# to build the tests enable the option BUILD_TESTING

cmake_minimum_required(VERSION 3.16)
project(complex_bessel VERSION 0.7.1 LANGUAGES Fortran CXX)

# set export variables
set(TARGETS_EXPORT_NAME "${PROJECT_NAME}-targets")
set(NAMESPACE "${PROJECT_NAME}::")

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

include(GNUInstallDirs)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")

###############################################################################
set(LIB complex_bessel)
add_library(${LIB} SHARED
        ${CMAKE_CURRENT_SOURCE_DIR}/src/amos_iso_c_fortran_wrapper.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/src/machine.for
        ${CMAKE_CURRENT_SOURCE_DIR}/src/zbesh.for
        )
add_library(${NAMESPACE}${LIB} ALIAS ${LIB})

# C++ includes will propagate with complex_bessel::complex_bessel target
# when linking it with some external files using cmake3.
target_include_directories(${LIB} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
        $<INSTALL_INTERFACE:include>
        )


# testing
if(${PROJECT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  set(main_project TRUE)
  option(BUILD_TESTING "Build the testing tree." OFF)
else()
  set(main_project FALSE)

  # if this is not the main project but BUILD_TESTING is set to TRUE
  # the tests for this project can be enabled by also setting
  # BUILD_TESTING_${PROJECT_NAME} to TRUE
  include(CMakeDependentOption)
  cmake_dependent_option(BUILD_TESTING_${PROJECT_NAME}
          "Build the testing tree for project ${PROJECT_NAME}." OFF
          "BUILD_TESTING;NOT main_project" OFF
          )
endif()

# build tests only by explicit request
if((main_project OR BUILD_TESTING_${PROJECT_NAME}) AND BUILD_TESTING)
  include(CTest)
  add_subdirectory(tests)
endif()

# generate package config files
set(project_config "${PROJECT_NAME}-config.cmake")
set(project_config_version "${PROJECT_NAME}-config-version.cmake")
set(cmake_files_dir "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles")
set(default_config_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
set(config_build_dir "${CMAKE_CURRENT_BINARY_DIR}/${default_config_install_dir}")
set(${PROJECT_NAME}_INSTALL_CMAKEDIR ${default_config_install_dir} CACHE PATH "Path into which the cmake files for project ${PROJECT_NAME} will be installed.")

# installation and export of targets
install(TARGETS ${LIB} EXPORT ${TARGETS_EXPORT_NAME}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT RuntimeLibraries
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT RuntimeLibraries
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )

install(EXPORT ${TARGETS_EXPORT_NAME}
        NAMESPACE
        ${NAMESPACE}
        DESTINATION
        ${${PROJECT_NAME}_INSTALL_CMAKEDIR}
        COMPONENT
        ${PROJECT_NAME}_Development
        )

install(
        FILES
        "${cmake_files_dir}/${project_config}"
        "${CMAKE_CURRENT_BINARY_DIR}/${project_config_version}"
        DESTINATION
        ${${PROJECT_NAME}_INSTALL_CMAKEDIR}
        COMPONENT
        ${PROJECT_NAME}_Development
)

install (DIRECTORY include/ DESTINATION include)

# export targets into build
export(EXPORT ${TARGETS_EXPORT_NAME}
        NAMESPACE
        ${NAMESPACE}
        FILE
        "${config_build_dir}/${TARGETS_EXPORT_NAME}.cmake"
        )

# create package config
# Variables needed by PackageConfig.cmake.in: PROJECT_NAME, TARGETS_EXPORT_NAME

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/PackageConfig.cmake.in "${cmake_files_dir}/${project_config}"
        INSTALL_DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKEDIR}
        )

configure_package_config_file(cmake/PackageConfig.cmake.in "${config_build_dir}/${project_config}"
        INSTALL_DESTINATION ${config_build_dir}
        INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        )

write_basic_package_version_file(
        ${project_config_version}
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)

#
