cmake_minimum_required(VERSION 3.22)

option(IFEM_USE_ANDES "Compile with ANDES shell element?" OFF)
option(IFEM_USE_FFLLIB "Compile with FFlLib?" OFF)
option(FFL_USE_MEMPOOL "Compile FFlLib with memory pool management?" OFF)

project(ShellEx LANGUAGES C CXX)

message(STATUS "Configuring project ${PROJECT_NAME}")

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

if(NOT TARGET Elasticity)
  add_subdirectory(.. Elasticity)
endif()

# ANDES shell elements (3rd-party Fortran-90 code)
if(IFEM_USE_ANDES)
  message(STATUS "Compiling with ANDES shell element support")
  enable_language(Fortran)
endif()

if(IFEM_USE_FFLLIB)
  message(STATUS "Compiling with use of FFlLib")
  add_subdirectory(FFlLib)
endif()

ifem_add_library(
  NAME
    ShellEx
  SOURCES
    AndesShell.C
    ASMu2DNastran.C
    SIMAndesShell.C
    SIMAndesSplit.C
    SIMShellModal.C
    WaveFunc.C
    ../Beam/ElasticBeam.C
    ../Beam/BeamProperty.C
    ../Linear/DynamicDriver.C
    ../Linear/ModalDriver.C
    ../Linear/MultiLoadCaseDriver.C
  HEADERS
    AndesShell.h
    ASMu2DNastran.h
    SIMAndesShell.h
    SIMAndesSplit.h
    SIMShellModal.h
    WaveFunc.h
  LIBRARIES
    Elasticity
)

target_include_directories(
    ShellEx
    PUBLIC
      ${PROJECT_SOURCE_DIR}/../Beam
      ${PROJECT_SOURCE_DIR}/../Linear
)

if(IFEM_USE_ANDES)
  target_sources(ShellEx PRIVATE
    ANDES/andes3ShellModule.f90
    ANDES/andes4ShellModule.f90
    ANDES/IFEM_wrapper.f90
    ANDES/kindModule.f90
    ANDES/manipMatrixModule.f90
    ANDES/pmatStiff.f90
    ANDES/strainAndStressUtils.f90
    ANDES/wavgmConstrEqn.f90
    ANDES/qbek35.f
    ANDES/qmrf35.f
    ANDES/teba35.f
    ANDES/tmrf35.f
    FFlLib/FFaLib/FFaAlgebra/FFaTensorTransformsInterface.f90
  )
  target_compile_options(ShellEx PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-cpp -frecursive -fwhole-file>)
  target_compile_definitions(ShellEx PUBLIC HAS_ANDES=1)
endif()

if(TARGET FFlLib)
  target_link_libraries(ShellEx PUBLIC FFlLib)
endif()

# Building executable
ifem_add_application(
  NAME
    ShellSim
  SOURCES
    main.C
  LIBRARIES
    ShellEx
)

# Regression tests
enable_testing()

ifem_add_regression_test(
  TARGET
    ShellSim
  CONDITIONS
    IFEM_USE_ANDES
    IFEM_USE_FFLLIB
  TEST_FILES
    BeamOnly-dyn.reg
    BeamOnly.reg
    DroppedBox.reg
    Fine-dyn.reg
    Fine-mlc.reg
    Fine-partial.reg
    Fine.reg
    L-beam10.reg
    LoadTest.reg
    Q4andBeam.reg
    Q4-dynp.reg
    Q4-dyn.reg
    Q4-generated.reg
    Q4noBeam.reg
    Q4.reg
    T3.reg
    Wheel.reg
)

ifem_add_regression_test(
  TARGET
    ShellSim
  CONDITIONS
    IFEM_USE_ANDES
    IFEM_USE_FFLLIB
    SPR_FOUND
  TEST_FILES
    Fine-dump.reg
    Fine-load.reg
)

ifem_add_vtf_test(
  TARGET
    ShellSim
  CONDITIONS
    IFEM_USE_ANDES
    IFEM_USE_FFLLIB
  TEST_FILES
    T3.vreg
    Q4.vreg
    DroppedBox.vreg
)

if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  if(IFEM_USE_ANDES AND IFEM_USE_FFLLIB AND SPR_FOUND)
    set_tests_properties(ShellSim+Fine-load.reg
                         PROPERTIES DEPENDS
                         ShellSim+Fine-dump.reg)
    set_tests_properties(ShellSim+Fine-dump.reg
                         PROPERTIES FIXTURES_SETUP Dump)
    set_tests_properties(ShellSim+Fine-load.reg
                         PROPERTIES FIXTURES_REQUIRED Dump)
  endif()
endif()

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else(IFEM_COMMON_APP_BUILD)
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(TARGET ShellSim)

# Installation
include(GNUInstallDirs)
install(TARGETS ShellSim DESTINATION ${CMAKE_INSTALL_BINDIR})
