cmake_minimum_required(VERSION 3.24)

project(MSTM VERSION 4.0.0 LANGUAGES Fortran)

option(MSTM_ENABLE_MPI "Build the distributed-memory MPI implementation" OFF)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR "MSTM currently requires GNU Fortran 15.2 or newer")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 15.2)
  message(FATAL_ERROR "GNU Fortran 15.2 or newer is required for Fortran 2023 mode")
endif()

set(mstm_mpi_source src/mpi_serial.f90)
if(MSTM_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  set(mstm_mpi_source src/mpi_parallel.f90)
endif()

set(mstm_sources
  src/intrinsics.f90
  ${mstm_mpi_source}
  src/core.f90
  src/fft_translation.f90
  src/scattering.f90
  src/solver.f90
  src/random_configuration.f90
  src/lmfit.f90
  src/input.f90
  src/main.f90
)

add_executable(mstm ${mstm_sources})

set_target_properties(mstm PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  Fortran_PREPROCESS OFF
  Fortran_FORMAT FREE
)
target_include_directories(mstm PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_compile_options(mstm PRIVATE -std=f2023)

if(MSTM_ENABLE_MPI)
  target_link_libraries(mstm PRIVATE MPI::MPI_Fortran)
endif()

find_program(FPRETTIFY_EXECUTABLE NAMES fprettify)
if(FPRETTIFY_EXECUTABLE)
  add_custom_target(format
    COMMAND "${FPRETTIFY_EXECUTABLE}"
      --config-file "${CMAKE_CURRENT_SOURCE_DIR}/.fprettify.rc"
      --silent
      ${mstm_sources}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Formatting Fortran sources"
    VERBATIM
  )
endif()

include(CTest)
if(BUILD_TESTING)
  function(add_mstm_example_test test_name input_name expected_outputs)
    set(test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/${test_name}")
    file(MAKE_DIRECTORY "${test_work_dir}")
    add_test(
      NAME "${test_name}"
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/examples/${input_name}"
        "-DMSTM_WORK_DIR=${test_work_dir}"
        "-DMSTM_EXPECTED_OUTPUTS=${expected_outputs}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSmokeTest.cmake"
    )
  endfunction()

  add_mstm_example_test(
    smoke
    mstm-2022b-fig1.inp
    "mstm-2022b-fig1.dat,nf-fig1.dat"
  )
  add_mstm_example_test(
    effective_medium_smoke
    effective-medium-smoke.inp
    "effective-medium-smoke.dat,anpeff.dat"
  )
endif()

install(TARGETS mstm RUNTIME DESTINATION bin)
