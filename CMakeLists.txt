cmake_minimum_required(VERSION 3.24)

project(MSTM VERSION 4.0.0 LANGUAGES Fortran)

find_package(LAPACK REQUIRED)

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
option(MSTM_ENABLE_COVERAGE "Instrument GNU Fortran targets for gcov coverage" OFF)
option(MSTM_ENABLE_NATIVE_ARCH "Generate code for the build machine's CPU" ON)
option(MSTM_ENABLE_VECTORIZATION_REPORT "Report optimized and missed vectorization opportunities" OFF)
set(MSTM_SERIAL_REFERENCE_EXECUTABLE "" CACHE FILEPATH
  "Serial MSTM executable used by MPI equivalence tests"
)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR "MSTM currently requires GNU Fortran 15.2 or newer")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 15.2)
  message(FATAL_ERROR "GNU Fortran 15.2 or newer is required for Fortran 2023 mode")
endif()

add_library(mstm_fortran_options INTERFACE)
target_compile_options(mstm_fortran_options INTERFACE -std=f2023)
if(MSTM_ENABLE_NATIVE_ARCH)
  target_compile_options(mstm_fortran_options INTERFACE -march=native)
endif()
if(MSTM_ENABLE_VECTORIZATION_REPORT)
  target_compile_options(mstm_fortran_options INTERFACE
    -fopt-info-vec-optimized
    -fopt-info-vec-missed
  )
endif()
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
if(MSTM_ENABLE_COVERAGE)
  find_program(GCOV_EXECUTABLE NAMES gcov REQUIRED)
  target_compile_options(mstm_fortran_options INTERFACE
    --coverage
    -O0
    -fprofile-abs-path
  )
  target_link_options(mstm_fortran_options INTERFACE --coverage)
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
  src/excitation.f90
  src/quadrature.f90
  src/surface.f90
  src/periodic_lattice.f90
  src/sphere_data.f90
  src/mie.f90
  src/translation_operator.f90
  src/translation_surface_interactions.f90
  src/translation_expansions.f90
  src/fft_translation.f90
  src/interaction_operators.f90
  src/scattering_interactions.f90
  src/scattering_amplitudes.f90
  src/scattering_efficiencies.f90
  src/random_orientation_scattering.f90
  src/near_field.f90
  src/direct_lu_solver.f90
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
  src/simulation_execution.f90
  src/configuration_averaging.f90
  src/random_orientation_averaging.f90
  src/incidence_averaging.f90
  src/simulation_averaging.f90
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
target_link_libraries(mstm_core PUBLIC mstm_bessel mstm_gpfa LAPACK::LAPACK)

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
      tests/direct_lu_solver_test.f90
      tests/excitation_test.f90
      tests/gpfa_test.f90
      tests/mie_test.f90
      tests/output_regression_test.f90
      tests/sphere_cluster_test.f90
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

  add_executable(mie_test tests/mie_test.f90)
  mstm_apply_fortran_options(mie_test)
  target_link_libraries(mie_test PRIVATE mstm_core)
  add_test(NAME mie_unit COMMAND mie_test)

  add_executable(direct_lu_solver_test tests/direct_lu_solver_test.f90)
  mstm_apply_fortran_options(direct_lu_solver_test)
  target_link_libraries(direct_lu_solver_test PRIVATE mstm_core)
  add_test(NAME direct_lu_solver_unit COMMAND direct_lu_solver_test)

  add_executable(excitation_test tests/excitation_test.f90)
  mstm_apply_fortran_options(excitation_test)
  target_link_libraries(excitation_test PRIVATE mstm_core)
  add_test(NAME excitation_unit COMMAND excitation_test)

  add_executable(sphere_cluster_test tests/sphere_cluster_test.f90)
  mstm_apply_fortran_options(sphere_cluster_test)
  target_link_libraries(sphere_cluster_test PRIVATE mstm_core)
  add_test(NAME sphere_cluster_unit COMMAND sphere_cluster_test)

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

  add_mstm_cli_test(cli_dry_run_alias --dry-run "${CMAKE_CURRENT_SOURCE_DIR}/tests/pec-sphere.inp")
  set_tests_properties(cli_dry_run_alias PROPERTIES PASS_REGULAR_EXPRESSION "Input validation successful")
  add_mstm_cli_test(pec_mixed_example_check --check "${CMAKE_CURRENT_SOURCE_DIR}/examples/pec-dielectric-cluster.inp")
  set_tests_properties(pec_mixed_example_check PROPERTIES PASS_REGULAR_EXPRESSION "PEC spheres: 1")

  set(validation_mode_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/validation_mode")
  file(MAKE_DIRECTORY "${validation_mode_work_dir}")
  set(validation_mode_mpi_arguments)
  if(MSTM_ENABLE_MPI)
    list(APPEND validation_mode_mpi_arguments
      "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
      "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
      "-DMSTM_MPI_RANKS=1"
    )
  endif()
  add_test(
    NAME cli_check_mode
    COMMAND "${CMAKE_COMMAND}"
      "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
      "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/validation-mode.inp"
      "-DMSTM_WORK_DIR=${validation_mode_work_dir}"
      "-DMSTM_FORBIDDEN_OUTPUTS=validation-run.log,validation-output.dat,validation-near-field.dat"
      ${validation_mode_mpi_arguments}
      -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunValidationMode.cmake"
  )

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
  add_mstm_failure_test(
    electric_dipole_inside_sphere_error
    electric-dipole-inside-sphere.inp
    "electric_dipole_position must lie outside every sphere"
  )
  add_mstm_failure_test(
    magnetic_current_intersects_sphere_error
    magnetic-current-intersects-sphere.inp
    "Magnetic-current segment 1 intersects sphere 1"
  )
  add_mstm_failure_test(missing_sphere_file_error missing-sphere-file.inp "Cannot open input file 'does-not-exist.pos'")
  add_mstm_failure_test(
    solver_nonconvergence_error
    solver-nonconvergence.inp
    "Fixed-orientation solver failed: iteration limit reached"
  )
  add_mstm_failure_test(
    nested_pec_error
    nested-pec.inp
    "only top-level solid PEC spheres are supported"
  )
  add_mstm_failure_test(
    malformed_pec_error
    malformed-pec.inp
    "Sphere record 1 in '${CMAKE_CURRENT_SOURCE_DIR}/tests/malformed-pec.inp' at line 6: PEC sphere records require x, y, z, radius, and PEC. Offending record: 0.d0,0.d0,0.d0,PEC"
  )

  add_mstm_cli_test(pec_sphere_smoke "${CMAKE_CURRENT_SOURCE_DIR}/tests/pec-sphere.inp")
  set_tests_properties(pec_sphere_smoke PROPERTIES PASS_REGULAR_EXPRESSION "PEC")

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
  add_mstm_regression_test(
    rigid_transform_invariants
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/rigid-transform-invariants.inp"
    rigid-transform-invariants.dat
    rigid_transform_invariants
  )
  add_mstm_regression_test(
    solver_validation
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/solver-validation.inp"
    solver-validation.dat
    solver_validation
  )
  add_mstm_regression_test(
    solver_reciprocity
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/solver-reciprocity.inp"
    solver-reciprocity.dat
    solver_reciprocity
  )
  add_mstm_regression_test(
    fft_translation_validation
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/fft-translation-validation.inp"
    fft-translation-validation.dat
    fft_translation_validation
  )
  add_mstm_regression_test(
    random_orientation_average
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/random-orientation-average.inp"
    random-orientation-average.dat
    random_orientation_average
  )
  add_mstm_regression_test(
    incidence_average
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/incidence-average.inp"
    incidence-average.dat
    incidence_average
  )
  add_mstm_regression_test(
    periodic_surface
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/periodic-surface.inp"
    periodic-surface.dat
    periodic_surface
  )
  add_mstm_regression_test(
    nested_sphere
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/nested-sphere.inp"
    nested-sphere.dat
    nested_sphere
  )
  add_mstm_regression_test(
    near_field_modes
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/near-field-modes.inp"
    "near-field-modes.dat,near-field-total.dat,near-field-scattered.dat"
    near_field_modes
  )
  add_mstm_regression_test(
    electric_dipole
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/electric-dipole.inp"
    "electric-dipole.dat,electric-dipole-near-field.dat"
    electric_dipole
  )
  add_mstm_regression_test(
    magnetic_current
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/magnetic-current.inp"
    "magnetic-current.dat,magnetic-current-near-field.dat"
    magnetic_current
  )
  add_mstm_regression_test(
    state_lifecycle
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/state-lifecycle.inp"
    state-lifecycle.dat
    state_lifecycle
  )
  add_mstm_regression_test(
    pec_sphere_regression
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/pec-sphere.inp"
    "pec-sphere.dat,pec-near-field.dat"
    pec_sphere
  )
  add_mstm_regression_test(
    pec_mixed_invariants
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/pec-mixed-invariants.inp"
    pec-mixed-invariants.dat
    pec_mixed_invariants
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

    set(mpi_solver_test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_two_rank_solver_validation")
    file(MAKE_DIRECTORY "${mpi_solver_test_work_dir}")
    add_test(
      NAME mpi_two_rank_solver_validation
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/solver-validation.inp"
        "-DMSTM_WORK_DIR=${mpi_solver_test_work_dir}"
        "-DMSTM_EXPECTED_OUTPUTS=solver-validation.dat"
        "-DMSTM_REGRESSION_EXECUTABLE=$<TARGET_FILE:output_regression_test>"
        "-DMSTM_REGRESSION_CASE=solver_validation"
        "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMSTM_MPI_RANKS=2"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSmokeTest.cmake"
    )

    set(mpi_fft_test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_two_rank_fft_validation")
    file(MAKE_DIRECTORY "${mpi_fft_test_work_dir}")
    add_test(
      NAME mpi_two_rank_fft_validation
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/fft-translation-validation.inp"
        "-DMSTM_WORK_DIR=${mpi_fft_test_work_dir}"
        "-DMSTM_EXPECTED_OUTPUTS=fft-translation-validation.dat"
        "-DMSTM_REGRESSION_EXECUTABLE=$<TARGET_FILE:output_regression_test>"
        "-DMSTM_REGRESSION_CASE=fft_translation_validation"
        "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMSTM_MPI_RANKS=2"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSmokeTest.cmake"
    )

    set(mpi_dipole_test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_two_rank_electric_dipole")
    file(MAKE_DIRECTORY "${mpi_dipole_test_work_dir}")
    add_test(
      NAME mpi_two_rank_electric_dipole
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/electric-dipole.inp"
        "-DMSTM_WORK_DIR=${mpi_dipole_test_work_dir}"
        "-DMSTM_EXPECTED_OUTPUTS=electric-dipole.dat,electric-dipole-near-field.dat"
        "-DMSTM_REGRESSION_EXECUTABLE=$<TARGET_FILE:output_regression_test>"
        "-DMSTM_REGRESSION_CASE=electric_dipole"
        "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMSTM_MPI_RANKS=2"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSmokeTest.cmake"
    )

    set(mpi_magnetic_current_test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_two_rank_magnetic_current")
    file(MAKE_DIRECTORY "${mpi_magnetic_current_test_work_dir}")
    add_test(
      NAME mpi_two_rank_magnetic_current
      COMMAND "${CMAKE_COMMAND}"
        "-DMSTM_EXECUTABLE=$<TARGET_FILE:mstm>"
        "-DMSTM_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/tests/magnetic-current.inp"
        "-DMSTM_WORK_DIR=${mpi_magnetic_current_test_work_dir}"
        "-DMSTM_EXPECTED_OUTPUTS=magnetic-current.dat,magnetic-current-near-field.dat"
        "-DMSTM_REGRESSION_EXECUTABLE=$<TARGET_FILE:output_regression_test>"
        "-DMSTM_REGRESSION_CASE=magnetic_current"
        "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
        "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        "-DMSTM_MPI_RANKS=2"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunSmokeTest.cmake"
    )

    if(MSTM_SERIAL_REFERENCE_EXECUTABLE)
      if(NOT EXISTS "${MSTM_SERIAL_REFERENCE_EXECUTABLE}")
        message(FATAL_ERROR
          "MSTM_SERIAL_REFERENCE_EXECUTABLE does not exist: ${MSTM_SERIAL_REFERENCE_EXECUTABLE}"
        )
      endif()

      function(add_mstm_parallel_equivalence_test test_name input_path output_name regression_case)
        set(test_work_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/${test_name}")
        file(MAKE_DIRECTORY "${test_work_dir}")
        add_test(
          NAME "${test_name}"
          COMMAND "${CMAKE_COMMAND}"
            "-DMSTM_SERIAL_EXECUTABLE=${MSTM_SERIAL_REFERENCE_EXECUTABLE}"
            "-DMSTM_MPI_EXECUTABLE=$<TARGET_FILE:mstm>"
            "-DMSTM_MPIEXEC=${MPIEXEC_EXECUTABLE}"
            "-DMSTM_MPI_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
            "-DMSTM_INPUT=${input_path}"
            "-DMSTM_OUTPUT_NAME=${output_name}"
            "-DMSTM_WORK_DIR=${test_work_dir}"
            "-DMSTM_REGRESSION_EXECUTABLE=$<TARGET_FILE:output_regression_test>"
            "-DMSTM_REGRESSION_CASE=${regression_case}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RunParallelEquivalence.cmake"
        )
        set_tests_properties("${test_name}" PROPERTIES PROCESSORS 4 TIMEOUT 300)
      endfunction()

      add_mstm_parallel_equivalence_test(
        serial_two_four_rank_solver_equivalence
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/solver-validation.inp"
        solver-validation.dat
        parallel_solver_equivalence
      )
      add_mstm_parallel_equivalence_test(
        serial_two_four_rank_fft_equivalence
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/fft-translation-validation.inp"
        fft-translation-validation.dat
        parallel_fft_equivalence
      )
      add_mstm_parallel_equivalence_test(
        serial_two_four_rank_dipole_equivalence
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/electric-dipole.inp"
        electric-dipole.dat
        parallel_dipole_equivalence
      )
      add_mstm_parallel_equivalence_test(
        serial_two_four_rank_magnetic_current_equivalence
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/magnetic-current.inp"
        magnetic-current.dat
        parallel_magnetic_current_equivalence
      )
    endif()
  endif()
endif()

if(MSTM_ENABLE_COVERAGE)
  if(NOT BUILD_TESTING)
    message(FATAL_ERROR "MSTM_ENABLE_COVERAGE requires BUILD_TESTING=ON")
  endif()
  add_custom_target(coverage
    COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure
    COMMAND "${CMAKE_COMMAND}"
      "-DMSTM_BUILD_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}"
      "-DMSTM_SOURCE_DIRECTORY=${CMAKE_CURRENT_SOURCE_DIR}"
      "-DMSTM_GCOV_EXECUTABLE=${GCOV_EXECUTABLE}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GenerateCoverage.cmake"
    DEPENDS
      mstm
      bessel_test
      concurrent_kernels_test
      direct_lu_solver_test
      excitation_test
      gpfa_test
      mie_test
      output_regression_test
      sphere_cluster_test
      translation_expansion_test
    COMMENT "Running tests and generating GNU coverage report"
    VERBATIM
  )
endif()

install(TARGETS mstm RUNTIME DESTINATION bin)
