set(UNIT_TEST_HELPERS
  test_assert_mod.f90
)

file(COPY fixtures DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

function(add_div3d_unit_test test_name test_source)
  add_executable(${test_name}
    ${UNIT_TEST_HELPERS}
    ${test_source}
  )
  set_target_properties(${test_name} PROPERTIES OUTPUT_NAME "${test_name}.exe")
  target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} PRIVATE libdiv3d_core libcommon libfxdr libstell ${NetCDF_Fortran_LIBRARY} ${NetCDF_C_LIBRARY} z stdc++ ${LAPACK_LIBRARIES})
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endfunction()

add_div3d_unit_test(test_math_routines test_math_routines.f90)
add_div3d_unit_test(test_inside_vessel test_inside_vessel.f90)
add_div3d_unit_test(test_read_parts test_read_parts.f90)
add_div3d_unit_test(test_intersections test_intersections.f90)
