cmake_minimum_required(VERSION 3.5)

project(ShellEx)

get_filename_component(PARENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${PARENT_SOURCE_DIR})

# ANDES shell elements (3rd-party Fortran-90 code)
OPTION(IFEM_USE_ANDES "Compile with ANDES shell element?" OFF)
IF(IFEM_USE_ANDES)
  MESSAGE(STATUS "Compiling with ANDES shell element support")
  ENABLE_LANGUAGE(Fortran)
  FILE(GLOB ANDES_SOURCES ANDES/*.f90)
  STRING(APPEND CMAKE_Fortran_FLAGS " -cpp -frecursive -fwhole-file")
  STRING(APPEND CMAKE_CXX_FLAGS " -DHAS_ANDES")
ENDIF(IFEM_USE_ANDES)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
FILE(GLOB ShellEx_SOURCES *.C)
FILE(GLOB ShellEx_HEADERS *.h)
ADD_LIBRARY(ShellEx ${ShellEx_SOURCES} ${ANDES_SOURCES})
