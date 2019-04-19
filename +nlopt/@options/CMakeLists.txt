# compile back-end MEX function for mexClass class
matlab_add_mex(NAME mexfcn SRC options_mexfcn.cpp mexObjectiveFunction.cpp mexConstraintFunction.cpp LINK_TO ${NLOPT_LIBRARY})

# if the compiler is MS VisualC++, add 'mexFunction' to the exported symbols
if(WIN32 AND MSVC)
  set_target_properties(mexfcn PROPERTIES LINK_FLAGS /EXPORT:mexFunction)
endif()

# install
file(RELATIVE_PATH DstRelativePath ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS mexfcn RUNTIME DESTINATION "${MATLAB_NLOPT_TOOLBOX_DIR}/${DstRelativePath}")
if (WIN32 AND BUILD_NLOPT_LIBS) 
  install(FILES "${NLOPT_ROOT}/bin/nlopt.dll" DESTINATION "${MATLAB_NLOPT_TOOLBOX_DIR}/${DstRelativePath}")
endif()

if(${Matlab_VERSION_STRING_INTERNAL} VERSION_LESS "8.4")
  set(OPTIONS_M "options.m.R2014a")
else()
  set(OPTIONS_M "options.m")
endif()
install(FILES ${OPTIONS_M}
  DESTINATION "${MATLAB_NLOPT_TOOLBOX_DIR}/${DstRelativePath}"
  RENAME "options.m")
