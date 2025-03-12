cmake_minimum_required(VERSION 3.5)

# Find Python for the test script
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Define a single test that runs the Python compare script
add_test(
  NAME test_ves_only
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TOK/CAT/test_ves_only
  COMMAND ${Python3_EXECUTABLE} compare_outputs.py
)
