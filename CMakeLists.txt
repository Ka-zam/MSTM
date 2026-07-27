cmake_minimum_required(VERSION 3.24)

project(MSTM VERSION 4.0.0 LANGUAGES Fortran)

set(mstm_generated_directory "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${mstm_generated_directory}")
configure_file(
  cmake/mstm_version.f90.in
  "${mstm_generated_directory}/mstm_version.f90"
  @ONLY
)

option(MSTM_ENABLE_MPI "Build the distributed-memory MPI implementation" OFF)
option(MSTM_ENABLE_WARNINGS "Enable GNU Fortran diagnostics" ON)
option(MSTM_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(MSTM_ENABLE_RUNTIME_CHECKS "Enable GNU Fortran runtime checks" OFF)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR "MSTM currently requires GNU Fortran 15.2 or newer")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 15.2)
  message(FATAL_ERROR "GNU Fortran 15.2 or newer is required for Fortran 2023 mode")
endif()

add_library(mstm_fortran_options INTERFACE)
target_compile_options(mstm_fortran_options INTERFACE -std=f2023)
if(MSTM_ENABLE_WARNINGS)
  target_compile_options(mstm_fortran_options INTERFACE
    -Wall
    -Wextra
    -Wcharacter-truncation
    -Wimplicit-interface
    -Werror=implicit-interface
    -Wno-compare-reals
    -Wno-conversion
    -Wno-do-subscript
    -Wno-maybe-uninitialized
    -Wno-unused-dummy-argument
    -Wno-unused-variable
    -Wsurprising
  )
endif()
if(MSTM_WARNINGS_AS_ERRORS)
  target_compile_options(mstm_fortran_options INTERFACE -Werror)
endif()
if(MSTM_ENABLE_RUNTIME_CHECKS)
  target_compile_options(mstm_fortran_options INTERFACE
    -fbacktrace
    -fcheck=all
    -ffpe-trap=invalid,zero,overflow
    -finit-integer=-2147483647
    -finit-real=snan
  )
endif()

function(mstm_apply_fortran_options target_name)
  target_link_libraries("${target_name}" PRIVATE mstm_fortran_options)
endfunction()

set(mstm_mpi_source src/mpi_serial.f90)
if(MSTM_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  set(mstm_mpi_source src/mpi_parallel.f90)
endif()

set(mstm_gpfa_sources
  src/gpfa_radix2.f90
  src/gpfa_radix3.f90
  src/gpfa_radix5.f90
  src/gpfa_dispatch.f90
  src/gpfa_controller.f90
  src/gpfa_setup.f90
)

set(mstm_bessel_source src/bessel_functions.f90)

set(mstm_core_sources
  ${mstm_mpi_source}
  src/numerical_tables.f90
  src/runtime_support.f90
  src/coefficient_indexing.f90
  src/angular_functions.f90
  src/wave_functions.f90
  src/quadrature.f90
  src/surface.f90
  src/periodic_lattice.f90
  src/sphere_data.f90
  src/mie.f90
  src/translation_operator.f90
  src/translation_surface_interactions.f90
  src/translation_expansions.f90
  src/fft_translation.f90
  src/scattering_interactions.f90
  src/scattering_amplitudes.f90
  src/scattering_efficiencies.f90
  src/random_orientation_scattering.f90
  src/near_field.f90
  src/solver.f90
  src/random_configuration_state.f90
  src/random_configuration_sorting.f90
  src/random_configuration_geometry.f90
  src/random_configuration_dynamics.f90
  src/random_sphere_configuration.f90
  src/levenberg_marquardt.f90
  src/input_value_parsing.f90
  src/input_state.f90
  src/input_parser.f90
  src/effective_medium_analysis.f90
  src/input_reporting.f90
  src/configuration_data.f90
  src/scattering_matrix_driver.f90
  src/input_execution.f90
)

add_library(mstm_constants STATIC src/constants.f90)
set_target_properties(mstm_constants PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/constant-modules"
  Fortran_PREPROCESS OFF
  Fortran_FORMAT FREE
)
target_include_directories(mstm_constants PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/constant-modules")
mstm_apply_fortran_options(mstm_constants)

add_library(mstm_bessel STATIC ${mstm_bessel_source})
set_target_properties(mstm_bessel PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bessel-modules"
  Fortran_PREPROCESS OFF
  Fortran_FORMAT FREE
)
target_include_directories(mstm_bessel PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/bessel-modules")
mstm_apply_fortran_options(mstm_bessel)
target_link_libraries(mstm_bessel PUBLIC mstm_constants)

add_library(mstm_gpfa STATIC ${mstm_gpfa_sources})
set_target_properties(mstm_gpfa PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/gpfa-modules"
  Fortran_PREPROCESS OFF
  Fortran_FORMAT FREE
)
target_include_directories(mstm_gpfa PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/gpfa-modules")
mstm_apply_fortran_options(mstm_gpfa)
target_link_libraries(mstm_gpfa PUBLIC mstm_constants)

add_library(mstm_core STATIC ${mstm_core_sources})

set_target_properties(mstm_core PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  Fortran_PREPROCESS OFF
  Fortran_FORMAT FREE
)
target_include_directories(mstm_core PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/modules")
mstm_apply_fortran_options(mstm_core)
target_link_libraries(mstm_core PUBLIC mstm_bessel mstm_gpfa)

if(MSTM_ENABLE_MPI)
  target_link_libraries(mstm_core PUBLIC MPI::MPI_Fortran)
endif()

add_executable(mstm
  "${mstm_generated_directory}/mstm_version.f90"
  src/main.f90
)
set_target_properties(mstm PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/executable-modules"
)
mstm_apply_fortran_options(mstm)
target_link_libraries(mstm PRIVATE mstm_core)

find_program(FPRETTIFY_EXECUTABLE NAMES fprettify)
if(FPRETTIFY_EXECUTABLE)
  add_custom_target(format
    COMMAND "${FPRETTIFY_EXECUTABLE}"
      --config-file "${CMAKE_CURRENT_SOURCE_DIR}/.fprettify.rc"
      --silent
      src/constants.f90
      ${mstm_bessel_source}
      ${mstm_gpfa_sources}
      ${mstm_core_sources}
      src/main.f90
      tests/bessel_test.f90
      tests/concurrent_kernels_test.f90
      tests/gpfa_test.f90
      tests/output_regression_test.f90
      tests/translation_expansion_test.f90
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Formatting Fortran sources"
    VERBATIM
  )
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(bessel_test tests/bessel_test.f90)
  mstm_apply_fortran_options(bessel_test)
  target_link_libraries(bessel_test PRIVATE mstm_bessel)
  add_test(NAME bessel_unit COMMAND bessel_test)

  add_executable(gpfa_test tests/gpfa_test.f90)
  mstm_apply_fortran_options(gpfa_test)
  target_link_libraries(gpfa_test PRIVATE mstm_gpfa)
  add_test(NAME gpfa_unit COMMAND gpfa_test)

  add_executable(translation_expansion_test tests/translation_expansion_test.f90)
  mstm_apply_fortran_options(translation_expansion_test)
  target_link_libraries(translation_expansion_test PRIVATE mstm_core)
  if(MSTM_ENABLE_MPI)
    add_test(
      NAME translation_expansion_unit
      COMMAND
        "${MPIEXEC_EXECUTABLE}"
        "${MPIEXEC_NUMPROC_FLAG}"
        1
        $<TARGET_FILE:translation_expansion_test>
    )
  else()
    add_test(NAME translation_expansion_unit COMMAND translation_expansion_test)
  endif()

  add_executable(concurrent_kernels_test tests/concurrent_kernels_test.f90)
  mstm_apply_fortran_options(concurrent_kernels_test)
  target_link_libraries(concurrent_kernels_test PRIVATE mstm_core)
  add_test(NAME concurrent_kernels_unit COMMAND concurrent_kernels_test)

  add_executable(output_regression_test tests/output_regression_test.f90)
  mstm_apply_fortran_options(output_regression_test)

  function(add_mstm_cli_test test_name)
    if(MSTM_ENABLE_MPI)
      add_test(
        NAME "${test_name}"
        COMMAND
          "${MPIEXEC_EXECUTABLE}"
          "${MPIEXEC_NUMPROC_FLAG}"
          1
          $<TARGET_FILE:mstm>
          ${ARGN}
      )
    else()
      add_test(NAME "${test_name}" COMMAND $<TARGET_FILE:mstm> ${ARGN})
    endif()
  endfunction()

  add_mstm_cli_test(cli_help)
  set_tests_properties(cli_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")

  add_mstm_cli_test(cli_help_alias help)
  set_tests_properties(cli_help_alias PROPERTIES PASS_REGULAR_EXPRESSION "Usage:")

  add_mstm_cli_test(cli_version --version)
  set_tests_properties(cli_version PROPERTIES PASS_REGULAR_EXPRESSION "MSTM ${PROJECT_VERSION}")

  function(add_mstm_failure_test test_name input_name expected_error)
    set(test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/${test_name}")
    file(MAKE_DIRECTORY "${test_work_dir}")
    set(mpi_test_arguments)
    if(MSTM_ENABLE_MPI)
      list(APPEND mpi_test_arguments
        "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMSTM_MPI_RANKS=1"
      )
    endif()
    add_test(
      NAME "${test_name}"
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/${input_name}"
        "-DMSTM_WORK_DIR=${test_work_dir}"
        "-DMSTM_EXPECTED_ERROR=${expected_error}"
        ${mpi_test_arguments}
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunExpectedFailure.cmake"
    )
  endfunction()

  add_mstm_failure_test(invalid_input_error invalid-option.inp "Unknown input parameter")
  add_mstm_failure_test(missing_sphere_file_error missing-sphere-file.inp "Cannot open input file 'does-not-exist.pos'")

  function(add_mstm_regression_test test_name input_path expected_outputs regression_case)
    set(test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/${test_name}")
    file(MAKE_DIRECTORY "${test_work_dir}")
    set(mpi_test_arguments)
    if(MSTM_ENABLE_MPI)
      list(APPEND mpi_test_arguments
        "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMSTM_MPI_RANKS=1"
      )
    endif()
    add_test(
      NAME "${test_name}"
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${input_path}"
        "-DMSTM_WORK_DIR=${test_work_dir}"
        "-DMSTM_EXPECTED_OUTPUTS=${expected_outputs}"
        "-DMSTM_REGRESSION_EXECUTABLE=$<TARGET_FILE:output_regression_test>"
        "-DMSTM_REGRESSION_CASE=${regression_case}"
        ${mpi_test_arguments}
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSmokeTest.cmake"
    )
  endfunction()

  add_mstm_regression_test(
    smoke
    "${CMAKE_CURRENT_SOURCE_DIR}/examples/mstm-2022b-fig1.inp"
    "mstm-2022b-fig1.dat,nf-fig1.dat"
    figure1
  )
  add_mstm_regression_test(
    effective_medium_smoke
    "${CMAKE_CURRENT_SOURCE_DIR}/examples/effective-medium-smoke.inp"
    "effective-medium-smoke.dat,anpeff.dat"
    effective_medium
  )
  add_mstm_regression_test(
    two_sphere_broadside
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/two-sphere-broadside.inp"
    two-sphere-broadside.dat
    two_sphere_broadside
  )

  if(MSTM_ENABLE_MPI)
    set(mpi_test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_two_rank_effective_medium")
    file(MAKE_DIRECTORY "${mpi_test_work_dir}")
    add_test(
      NAME mpi_two_rank_effective_medium
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/examples/effective-medium-smoke.inp"
        "-DMSTM_WORK_DIR=${mpi_test_work_dir}"
        "-DMSTM_EXPECTED_OUTPUTS=effective-medium-smoke.dat,anpeff.dat"
        "-DMSTM_REGRESSION_EXECUTABLE=$<TARGET_FILE:output_regression_test>"
        "-DMSTM_REGRESSION_CASE=effective_medium"
        "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMSTM_MPI_RANKS=2"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSmokeTest.cmake"
    )
  endif()
endif()

install(TARGETS mstm RUNTIME DESTINATION bin)
