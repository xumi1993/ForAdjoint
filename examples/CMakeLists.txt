file(GLOB APP_SOURCE "${PROJECT_SOURCE_DIR}/examples/ex*.f90")

foreach(exec_source ${APP_SOURCE})
  get_filename_component(exec_name ${exec_source} NAME_WE)
  add_executable(${exec_name} ${exec_source})
  set_target_properties(${exec_name} PROPERTIES PREFIX "x")
  target_link_libraries(${exec_name} foradjoint)
  set_target_properties(${exec_name} PROPERTIES LINKER_LANGUAGE Fortran)

endforeach()
