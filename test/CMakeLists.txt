cmake_minimum_required(VERSION 3.5)

# Find Python for the test script
find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_subdirectory(unit)

# Define a single test that runs the Python compare script
add_test(
  NAME test_ves_only
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TOK/CAT/test_ves_only
  COMMAND ${Python3_EXECUTABLE} compare_outputs.py
)

add_test(
  NAME test_ves_and_parts
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TOK/CAT/test_ves_and_parts
  COMMAND ${Python3_EXECUTABLE} compare_outputs.py
)

add_test(
  NAME test_poincare_tok
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TOK/CAT/poincare
  COMMAND ${Python3_EXECUTABLE} compare_outputs.py
)

add_test(
  NAME test_switch_cases
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TOK/CAT
  COMMAND ${Python3_EXECUTABLE} test_switch_cases.py
)
