if(NOT DEFINED MSTM_EXECUTABLE
    OR NOT DEFINED MSTM_INPUT
    OR NOT DEFINED MSTM_WORK_DIR
    OR NOT DEFINED MSTM_EXPECTED_OUTPUTS)
  message(FATAL_ERROR
    "The smoke test requires an executable, input, working directory, and expected outputs"
  )
endif()

string(REPLACE "," ";" expected_outputs "${MSTM_EXPECTED_OUTPUTS}")
foreach(output_file IN LISTS expected_outputs)
  file(REMOVE "${MSTM_WORK_DIR}/${output_file}")
endforeach()

execute_process(
  COMMAND "${MSTM_EXECUTABLE}" "${MSTM_INPUT}"
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
