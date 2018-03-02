# Find the NLopt include and library
# This module defines
# NLOPT_INCLUDE_DIR, where to find nlopt.h
# NLOPT_LIBRARY, the library to link against to use NLopt.
# NLOPT_SHAREDLIBRARY
# NLOPT_FOUND, If false, do not try to use NLopt.
# NLOPT_ROOT, if this module use this path to find NLopt header
# and libraries.
#
# In Windows, it looks for NLOPT_DIR environment variable if defined

include(FindPackageHandleStandardArgs)

# ###################################
# Exploring the possible NLOPT_ROOT
if(WIN32)
    if(NOT DEFINED NLOPT_ROOT)
        #look in the system environmental variable named "NLOPT_DIR"
        set(NLOPT_ROOT $ENV{NLOPT_DIR} CACHE PATH "NLOPT installation root path")
    endif()
    if(DEFINED NLOPT_ROOT AND NOT EXISTS ${NLOPT_ROOT})
        # if NLOPT_ROOT specified but erroneous
        message(WARNING "[NLOPT] the specified path for NLOPT_ROOT does not exist (${NLOPT_ROOT})")
    endif()
else()
  # Linux & Mac should not need one
  set(NLOPT_ROOT "" CACHE PATH "NLOPT installation root path")
endif()

include (GNUInstallDirs) # defines CMAKE_INSTALL_LIBDIR & CMAKE_INSTALL_INCLUDEDIR

# Find header and lib directories
find_path(NLOPT_INCLUDE_DIR nlopt.h
    PATHS
    ${NLOPT_ROOT}
    ${CMAKE_INSTALL_FULL_INCLUDEDIR}
    ${CMAKE_INSTALL_FULL_OLDINCLUDEDIR}
    "$ENV{ProgramFiles}/NLopt"
    "$ENV{USERPROFILE}/AppData/Local"
    "$ENV{USERPROFILE}/AppData/Local/Programs"
    "$ENV{USERPROFILE}/AppData/Local/include"
    "$ENV{SystemDrive}"
    PATH_SUFFIXES NLOPT
    DOC "Location of NLopt header file"
)

find_library(NLOPT_LIBRARY
    NAMES nlopt libnlopt
    PATHS
    ${NLOPT_ROOT}
    ${CMAKE_INSTALL_FULL_LIBDIR}
    "$ENV{ProgramFiles}/NLopt"
    "$ENV{USERPROFILE}/AppData/Local"
    "$ENV{USERPROFILE}/AppData/Local/Programs"
    "$ENV{USERPROFILE}/AppData/Local/lib"
    "$ENV{SystemDrive}"
    DOC "Location of NLopt library"
)

find_package_handle_standard_args(
  NLOPT
  REQUIRED_VARS NLOPT_INCLUDE_DIR NLOPT_LIBRARY)
