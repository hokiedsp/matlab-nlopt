cmake_minimum_required(VERSION 2.8.11) 
# for FindMatlab support (may require a later version to detect the latest Matlab release)

project (matlab-nlopt)

list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

option (BUILD_NLOPT_LIBS "Build NLopt as a shared library and place it along the MEX file" ON)
set(NLOPT_LIB_DIR "" CACHE PATH "If NLopt not built together, specify the location of NLopt library directory")
set(NLOPT_INCLUDE_DIR "" CACHE PATH "If NLopt not built together, specify the lcoation of NLopt include directory")

# Set the default installation directory
if (WIN32)
  set(MATLAB_USER_DIR "$ENV{USERPROFILE}/Documents/MATLAB")
else()
  set(MATLAB_USER_DIR "$ENV{home}/Documents/MATLAB")
endif()
if (NOT EXISTS ${MATLAB_USER_DIR})
  # if not found, show message and use the default CMake install directory
  message("Matlab user directory not found at ${MATLAB_USER_DIR}. It will be created if intallation directory not changed.")
endif()
set(MATLAB_NLOPT_INSTALL_DIR "${MATLAB_USER_DIR}/nlopt")
set(MATLAB_NLOPT_INSTALL_DIR ${MATLAB_NLOPT_INSTALL_DIR} CACHE PATH "Installation directory for Matlab-NLopt")

# set global C++ compilation/linking options

# Look for MATLAB API library paths
find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY)

# Build/find NLopt
if (BUILD_NLOPT_LIBS) 
  if (WIN32)
    include(ExternalProject)
    set(NLOPT_ROOT "${CMAKE_BINARY_DIR}/nlopt")
    ExternalProject_Add(nlopt
                        PREFIX ${NLOPT_ROOT}
                        GIT_REPOSITORY "https://github.com/stevengj/nlopt.git"
                        BINARY_DIR "${NLOPT_ROOT}/build"
                        CMAKE_ARGS
                          -DCMAKE_BUILD_TYPE=Release
                          -DNLOPT_CXX:BOOL=OFF
                          -DBUILD_SHARED_LIBS:BOOL=ON
                          -DNLOPT_PYTHON:BOOL=OFF
                          -DNLOPT_OCTAVE:BOOL=OFF
                          -DNLOPT_MATLAB:BOOL=OFF
                          -DNLOPT_GUILE:BOOL=OFF
                          -DNLOPT_SWIG:BOOL=OFF
                          -DCMAKE_INSTALL_PREFIX=${NLOPT_ROOT}
                       )
    # set include & library path
    set(NLOPT_INCLUDE_DIR,"${NLOPT_ROOT}/include")
    set(NLOPT_LIBRARY,"${NLOPT_ROOT}/lib/nlopt.lib")
  else() # 
    find_package(NLopt) # if installed sets NLOPT_FOUND, NLOPT_INCLUDE_DIR, & NLOPT_LIBRARY
    if (NOT NLOPT_FOUND) # install the full package with default options (sans )

      include(ExternalProject)
      ExternalProject_Add(nlopt
                          PREFIX ${NLOPT_ROOT}
                          GIT_REPOSITORY "https://github.com/stevengj/nlopt.git"
                          BINARY_DIR "${NLOPT_ROOT}/build"
                          CMAKE_ARGS
                            -DCMAKE_BUILD_TYPE=Release
                            -DNLOPT_CXX:BOOL=OFF
                            -DBUILD_SHARED_LIBS:BOOL=ON
                         )
    endif()
  endif()
else()
  find_package(NLopt REQUIRED) # sets NLOPT_INCLUDE_DIR & NLOPT_LIBRARY
endif()

include_directories(${NLOPT_INCLUDE_DIR})
include_directories("mexutils/include")

# set globally to ignore MSVC warning on STD library calls
if(WIN32 AND MSVC)
  # add_compile_options(/wd4996)
endif()


# force to change the CMake install folder
set(CMAKE_INSTALL_PREFIX ${MEXCPP_DEMO_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)

# add the install targets
install(DIRECTORY +nlopt DESTINATION .)

# first we can indicate the documentation build as an option and set it to ON by default
# option(BUILD_DOC "Build documentation" ON)

# documentation build
# find_package(Doxygen) # check if Doxygen is installed
# if (DOXYGEN_FOUND AND CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
#     # set input and output files
#     set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
#     set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#     # request to configure the file
#     configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#     message("Doxygen build started")

#     # note the option ALL which allows to build the docs together with the application
#     add_custom_target( doc_doxygen ALL
#         COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#         COMMENT "Generating API documentation with Doxygen"
#         VERBATIM )
# else ()
#   message("Doxygen need to be installed to generate the doxygen documentation")
# endif ()


# traverse subfolders (note mexutils not included as it is a header only library)
# add_subdirectory("+nlopt")
