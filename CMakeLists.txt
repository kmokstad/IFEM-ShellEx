CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

PROJECT(ShellEx)

# Add local modules
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${PROJECT_SOURCE_DIR}/../../../cmake/Modules)

# Required packages
if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  string(APPEND CMAKE_CXX_FLAGS " ${IFEM_CXX_FLAGS}")
endif()

INCLUDE_DIRECTORIES(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/..)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT TARGET Elasticity)
  add_subdirectory(../ Elasticity)
endif()

# ANDES shell elements (3rd-party Fortran-90 code)
OPTION(IFEM_USE_ANDES "Compile with ANDES shell element?" OFF)
IF(IFEM_USE_ANDES)
  MESSAGE(STATUS "Compiling with ANDES shell element support")
  ENABLE_LANGUAGE(Fortran)
  FILE(GLOB ANDES_SOURCES ANDES/*.f90)
  STRING(APPEND CMAKE_Fortran_FLAGS " -cpp -frecursive -fwhole-file")
  STRING(APPEND CMAKE_CXX_FLAGS " -DHAS_ANDES")
ENDIF(IFEM_USE_ANDES)

OPTION(IFEM_USE_FFLLIB "Compile with FFlLib?" OFF)
IF(IFEM_USE_FFLLIB)
  MESSAGE(STATUS "Compiling with use of FFlLib")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/FFlLib)
  add_subdirectory(FFlLib)
  STRING(APPEND CMAKE_CXX_FLAGS " -DHAS_FFLLIB")
ENDIF(IFEM_USE_FFLLIB)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
FILE(GLOB ShellEx_SOURCES *.C)
FILE(GLOB ShellEx_HEADERS *.h)
ADD_LIBRARY(ShellEx ${ShellEx_SOURCES} ${ANDES_SOURCES} )
target_link_libraries(ShellEx Elasticity ${IFEM_LIBRARIES})
if(TARGET FFlLib)
  target_link_libraries(ShellEx FFlLib)
endif(TARGET FFlLib)
