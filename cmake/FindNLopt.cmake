# Find the NLopt include and library
# This module defines
# NLOPT_INCLUDE_DIR, where to find nlopt.h
# NLOPT_LIBRARY, the library to link against to use NLopt.
# NLOPT_SHAREDLIBRARY
# NLOPT_FOUND, If false, do not try to use NLopt.
# NLOPT_ROOT_DIR, if this module use this path to find NLopt header
# and libraries.
#
# In Windows, it looks for NLOPT_DIR environment variable if defined

include(FindPackageHandleStandardArgs)

# ###################################
# Exploring the possible NLOPT_ROOT_DIR
if (NOT NLOPT_ROOT_DIR)
    set(NLOPT_ROOT_DIR $ENV{NLOPT_DIR} CACHE PATH "NLOPT installation root path")
endif (NOT NLOPT_ROOT_DIR)

include (GNUInstallDirs) # defines CMAKE_INSTALL_LIBDIR & CMAKE_INSTALL_INCLUDEDIR

# Find header and lib directories
find_path(NLOPT_INCLUDE_DIR nlopt.h
    PATHS
    ${NLOPT_ROOT_DIR}
    ${CMAKE_INSTALL_FULL_INCLUDEDIR}
    ${CMAKE_INSTALL_FULL_OLDINCLUDEDIR}
    "$ENV{ProgramFiles}"
    "$ENV{USERPROFILE}/AppData/Local"
    "$ENV{USERPROFILE}/AppData/Local/Programs"
    "$ENV{USERPROFILE}/AppData/Local/include"
    "$ENV{SystemDrive}"
    PATH_SUFFIXES nlopt/include
    DOC "Location of NLopt header file"
)

find_library(NLOPT_LIBRARY
    NAMES nlopt
    PATHS
        ${NLOPT_ROOT_DIR}
        ${CMAKE_INSTALL_FULL_LIBDIR}
        "$ENV{ProgramFiles}"
        "$ENV{USERPROFILE}/AppData/Local"
        "$ENV{USERPROFILE}/AppData/Local/Programs"
        "$ENV{USERPROFILE}/AppData/Local/lib"
        "$ENV{SystemDrive}"
    PATH_SUFFIXES nlopt/lib
    DOC "Location of NLopt library"
)

find_package_handle_standard_args(
  NLOPT
  REQUIRED_VARS NLOPT_INCLUDE_DIR NLOPT_LIBRARY)
