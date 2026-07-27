if(NOT DEFINED MSTM_SERIAL_EXECUTABLE
    OR NOT DEFINED MSTM_MPI_EXECUTABLE
    OR NOT DEFINED MSTM_MPIEXEC
    OR NOT DEFINED MSTM_MPI_NUMPROC_FLAG
    OR NOT DEFINED MSTM_INPUT
    OR NOT DEFINED MSTM_OUTPUT_NAME
    OR NOT DEFINED MSTM_WORK_DIR
    OR NOT DEFINED MSTM_REGRESSION_EXECUTABLE
    OR NOT DEFINED MSTM_REGRESSION_CASE)
  message(FATAL_ERROR "Parallel equivalence test settings are incomplete")
endif()

set(serial_directory "${MSTM_WORK_DIR}/serial")
set(two_rank_directory "${MSTM_WORK_DIR}/mpi-2")
set(four_rank_directory "${MSTM_WORK_DIR}/mpi-4")
file(REMOVE_RECURSE "${serial_directory}" "${two_rank_directory}" "${four_rank_directory}")
file(MAKE_DIRECTORY "${serial_directory}" "${two_rank_directory}" "${four_rank_directory}")

function(run_mstm_case label working_directory)
  execute_process(
    COMMAND ${ARGN}
    WORKING_DIRECTORY "${working_directory}"
    RESULT_VARIABLE run_result
    OUTPUT_VARIABLE run_stdout
    ERROR_VARIABLE run_stderr
    TIMEOUT 180
  )
  if(NOT run_result EQUAL 0)
    message(FATAL_ERROR
      "${label} exited with ${run_result}\nstdout:\n${run_stdout}\nstderr:\n${run_stderr}"
    )
  endif()
  if(NOT EXISTS "${working_directory}/${MSTM_OUTPUT_NAME}")
    message(FATAL_ERROR "${label} did not create ${MSTM_OUTPUT_NAME}")
  endif()
endfunction()

run_mstm_case(
  "Serial reference"
  "${serial_directory}"
  "${MSTM_SERIAL_EXECUTABLE}"
  "${MSTM_INPUT}"
)
run_mstm_case(
  "Two-rank MPI run"
  "${two_rank_directory}"
  "${MSTM_MPIEXEC}"
  "${MSTM_MPI_NUMPROC_FLAG}"
  2
  "${MSTM_MPI_EXECUTABLE}"
  "${MSTM_INPUT}"
)
run_mstm_case(
  "Four-rank MPI run"
  "${four_rank_directory}"
  "${MSTM_MPIEXEC}"
  "${MSTM_MPI_NUMPROC_FLAG}"
  4
  "${MSTM_MPI_EXECUTABLE}"
  "${MSTM_INPUT}"
)

execute_process(
  COMMAND "${MSTM_REGRESSION_EXECUTABLE}"
    "${MSTM_REGRESSION_CASE}"
    "${serial_directory}"
    "${two_rank_directory}"
    "${four_rank_directory}"
  RESULT_VARIABLE comparison_result
  OUTPUT_VARIABLE comparison_stdout
  ERROR_VARIABLE comparison_stderr
)
if(NOT comparison_result EQUAL 0)
  message(FATAL_ERROR
    "Parallel result comparison exited with ${comparison_result}"
    "\nstdout:\n${comparison_stdout}\nstderr:\n${comparison_stderr}"
  )
endif()
