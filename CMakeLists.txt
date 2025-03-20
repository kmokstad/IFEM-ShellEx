cmake_minimum_required(VERSION 3.5)

project(ShellEx)
message(STATUS "Configuring project ${PROJECT_NAME}")

if(NOT IFEM_CONFIGURED)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../../cmake/Modules)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  string(APPEND CMAKE_CXX_FLAGS " ${IFEM_CXX_FLAGS}")
endif()
include_directories(${IFEM_INCLUDES})

get_filename_component(PARENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
if(IS_DIRECTORY ${PARENT_SOURCE_DIR}/Linear)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${PARENT_SOURCE_DIR}
                      ${PARENT_SOURCE_DIR}/Linear ${PARENT_SOURCE_DIR}/Beam)
else()
  message(FATAL_ERROR "Need IFEM-Elasticity as parent directory.")
endif()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT TARGET Elasticity)
  add_subdirectory(.. Elasticity)
endif()

# ANDES shell elements (3rd-party Fortran-90 code)
OPTION(IFEM_USE_ANDES "Compile with ANDES shell element?" OFF)
IF(IFEM_USE_ANDES)
  MESSAGE(STATUS "Compiling with ANDES shell element support")
  ENABLE_LANGUAGE(Fortran)
  FILE(GLOB ANDES_SOURCES ANDES/*.f90 ANDES/*.f FFlLib/FFaLib/FFaAlgebra/*.f90)
  STRING(APPEND CMAKE_Fortran_FLAGS " -cpp -frecursive -fwhole-file")
  STRING(APPEND CMAKE_CXX_FLAGS " -DHAS_ANDES")
ENDIF(IFEM_USE_ANDES)

OPTION(IFEM_USE_FFLLIB "Compile with FFlLib?" OFF)
IF(IFEM_USE_FFLLIB)
  MESSAGE(STATUS "Compiling with use of FFlLib")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/FFlLib)
  STRING(APPEND CMAKE_CXX_FLAGS " -DHAS_FFLLIB -DFFL_REFCOUNT")
  IF(IFEM_USE_ANDES)
    STRING(APPEND CMAKE_Fortran_FLAGS " -DHAS_FFLLIB")
  ENDIF(IFEM_USE_ANDES)
  ADD_SUBDIRECTORY(FFlLib)
ENDIF(IFEM_USE_FFLLIB)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
FILE(GLOB ShellEx_SOURCES *.C ../Beam/ElasticBeam.C ../Beam/BeamProperty.C ../Linear/*Driver.C)
FILE(GLOB ShellEx_HEADERS *.h ../Beam/ElasticBeam.h ../Beam/BeamProperty.h)
ADD_LIBRARY(ShellEx ${ShellEx_SOURCES} ${ANDES_SOURCES})
TARGET_LINK_LIBRARIES(ShellEx Elasticity ${IFEM_LIBRARIES})
IF(TARGET FFlLib)
  TARGET_LINK_LIBRARIES(ShellEx FFlLib)
ENDIF(TARGET FFlLib)

# Building executable
add_executable(ShellSim main.C)
target_link_libraries(ShellSim ShellEx)

# Installation
install(TARGETS ShellSim DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} \\
                      ${PROJECT_BINARY_DIR} \\
                      ${PARENT_SOURCE_DIR} \\
                      ${PARENT_SOURCE_DIR}/Beam \\
                      ${PARENT_SOURCE_DIR}/Linear")
add_doc_target(ShellSim ShellSim)

# Regression tests
enable_testing()
include(IFEMTesting)

if(IFEM_USE_ANDES AND TARGET FFlLib)
  file(GLOB SHELL_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                     ${PROJECT_SOURCE_DIR}/Test/*.reg)
endif()

if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  foreach(TESTFILE ${SHELL_TESTFILES})
    ifem_add_test(${TESTFILE} ShellSim)
  endforeach()
endif()

list(APPEND TEST_APPS ShellSim)
if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else(IFEM_COMMON_APP_BUILD)
  add_check_target()
endif()
