if(NOT DEFINED MSTM_EXECUTABLE
    OR NOT DEFINED MSTM_INPUT
    OR NOT DEFINED MSTM_WORK_DIR
    OR NOT DEFINED MSTM_EXPECTED_OUTPUTS
    OR NOT DEFINED MSTM_REGRESSION_EXECUTABLE
    OR NOT DEFINED MSTM_REGRESSION_CASE)
  message(FATAL_ERROR
    "The smoke test requires executable, input, output, and regression-check settings"
  )
endif()

string(REPLACE "," ";" expected_outputs "${MSTM_EXPECTED_OUTPUTS}")
foreach(output_file IN LISTS expected_outputs)
  file(REMOVE "${MSTM_WORK_DIR}/${output_file}")
endforeach()

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
  TIMEOUT 180
)

if(NOT mstm_result EQUAL 0)
  message(FATAL_ERROR
    "MSTM exited with ${mstm_result}\nstdout:\n${mstm_stdout}\nstderr:\n${mstm_stderr}"
  )
endif()

foreach(output_file IN LISTS expected_outputs)
  if(NOT EXISTS "${MSTM_WORK_DIR}/${output_file}")
    message(FATAL_ERROR "MSTM did not create ${output_file}")
  endif()
endforeach()

execute_process(
  COMMAND "${MSTM_REGRESSION_EXECUTABLE}" "${MSTM_REGRESSION_CASE}" "${MSTM_WORK_DIR}"
  RESULT_VARIABLE regression_result
  OUTPUT_VARIABLE regression_stdout
  ERROR_VARIABLE regression_stderr
)

if(NOT regression_result EQUAL 0)
  message(FATAL_ERROR
    "Numerical regression check exited with ${regression_result}\nstdout:\n${regression_stdout}\nstderr:\n${regression_stderr}"
  )
endif()
