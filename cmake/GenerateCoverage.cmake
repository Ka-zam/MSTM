if(NOT DEFINED MSTM_BUILD_DIRECTORY
    OR NOT DEFINED MSTM_SOURCE_DIRECTORY
    OR NOT DEFINED MSTM_GCOV_EXECUTABLE)
  message(FATAL_ERROR "Coverage reporting requires build, source, and gcov paths")
endif()

set(coverage_directory "${MSTM_BUILD_DIRECTORY}/coverage")
set(coverage_report "${coverage_directory}/coverage.txt")
file(REMOVE_RECURSE "${coverage_directory}")
file(MAKE_DIRECTORY "${coverage_directory}")

file(GLOB_RECURSE coverage_data "${MSTM_BUILD_DIRECTORY}/CMakeFiles/*.gcda")
list(FILTER coverage_data INCLUDE REGEX
  "/mstm_(core|bessel|gpfa|constants)\\.dir/"
)
if(NOT coverage_data)
  message(FATAL_ERROR "No coverage data found; build and run the tests first")
endif()

file(WRITE "${coverage_report}"
  "MSTM GNU coverage report\n"
  "Source: ${MSTM_SOURCE_DIRECTORY}/src\n\n"
)

foreach(data_file IN LISTS coverage_data)
  execute_process(
    COMMAND "${MSTM_GCOV_EXECUTABLE}"
      --branch-probabilities
      --branch-counts
      --function-summaries
      --no-output
      "${data_file}"
    WORKING_DIRECTORY "${coverage_directory}"
    RESULT_VARIABLE gcov_result
    OUTPUT_VARIABLE gcov_output
    ERROR_VARIABLE gcov_error
  )
  if(NOT gcov_result EQUAL 0)
    message(FATAL_ERROR
      "gcov failed for ${data_file}\n${gcov_output}\n${gcov_error}"
    )
  endif()
  file(APPEND "${coverage_report}" "${gcov_output}")
endforeach()

message(STATUS "Coverage report: ${coverage_report}")
