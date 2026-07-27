if(NOT DEFINED MSTM_EXECUTABLE
    OR NOT DEFINED MSTM_INPUT
    OR NOT DEFINED MSTM_WORK_DIR
    OR NOT DEFINED MSTM_EXPECTED_ERROR)
  message(FATAL_ERROR "The failure test requires executable, input, work directory, and expected error settings")
endif()

set(mstm_command "${MSTM_EXECUTABLE}" "${MSTM_INPUT}")
if(DEFINED MSTM_MPIEXEC)
  set(mstm_command
    "${MSTM_MPIEXEC}"
    "${MSTM_MPI_NUMPROC_FLAG}"
    "${MSTM_MPI_RANKS}"
    ${mstm_command}
  )
endif()

execute_process(
  COMMAND ${mstm_command}
  WORKING_DIRECTORY "${MSTM_WORK_DIR}"
  RESULT_VARIABLE mstm_result
  OUTPUT_VARIABLE mstm_stdout
  ERROR_VARIABLE mstm_stderr
  TIMEOUT 30
)

if(mstm_result EQUAL 0)
  message(FATAL_ERROR "MSTM unexpectedly accepted invalid input")
endif()

set(mstm_output "${mstm_stdout}\n${mstm_stderr}")
string(FIND "${mstm_output}" "${MSTM_EXPECTED_ERROR}" error_position)
if(error_position EQUAL -1)
  message(FATAL_ERROR
    "Expected '${MSTM_EXPECTED_ERROR}' in output\nstdout:\n${mstm_stdout}\nstderr:\n${mstm_stderr}"
  )
endif()
