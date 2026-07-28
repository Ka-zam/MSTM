if(NOT DEFINED MSTM_EXECUTABLE
    OR NOT DEFINED MSTM_INPUT
    OR NOT DEFINED MSTM_WORK_DIR
    OR NOT DEFINED MSTM_FORBIDDEN_OUTPUTS)
  message(FATAL_ERROR "The validation-mode test requires executable, input, work directory, and forbidden outputs")
endif()

string(REPLACE "," ";" forbidden_outputs "${MSTM_FORBIDDEN_OUTPUTS}")
foreach(output_file IN LISTS forbidden_outputs)
  file(REMOVE "${MSTM_WORK_DIR}/${output_file}")
endforeach()
file(REMOVE "${MSTM_WORK_DIR}/temp_pos.dat")

set(mstm_command "${MSTM_EXECUTABLE}" --check "${MSTM_INPUT}")
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

if(NOT mstm_result EQUAL 0)
  message(FATAL_ERROR
    "Validation mode exited with ${mstm_result}\nstdout:\n${mstm_stdout}\nstderr:\n${mstm_stderr}"
  )
endif()

set(mstm_output "${mstm_stdout}\n${mstm_stderr}")
foreach(expected_text IN ITEMS
    "Input validation successful"
    "PEC spheres: 1"
    "Parsed sphere geometry")
  string(FIND "${mstm_output}" "${expected_text}" text_position)
  if(text_position EQUAL -1)
    message(FATAL_ERROR "Validation output omitted '${expected_text}'\n${mstm_output}")
  endif()
endforeach()

foreach(output_file IN LISTS forbidden_outputs)
  if(EXISTS "${MSTM_WORK_DIR}/${output_file}")
    message(FATAL_ERROR "Validation mode unexpectedly created ${output_file}")
  endif()
endforeach()
if(EXISTS "${MSTM_WORK_DIR}/temp_pos.dat")
  message(FATAL_ERROR "Embedded sphere data unexpectedly created temp_pos.dat")
endif()
