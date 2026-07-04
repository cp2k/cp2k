#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

if(NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "OUTPUT_FILE is not defined")
endif()

if(NOT DEFINED SOURCE_DIR)
  set(SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/..")
endif()

string(TIMESTAMP CP2K_BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%S")

# Get the latest abbreviated commit hash of the working branch. As a fall back,
# e.g. in a container without a .git directory, try reading a file named
# "REVISION". Running this here (rather than once at configure time) means a new
# commit is picked up on the next build, without having to re-run cmake.
execute_process(
  COMMAND bash -c "git log -1 --format=%h || cat REVISION"
  WORKING_DIRECTORY "${SOURCE_DIR}"
  ERROR_QUIET
  OUTPUT_VARIABLE CP2K_GIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT CP2K_GIT_HASH)
  set(CP2K_GIT_HASH "unknown")
endif()

set(CP2K_BUILD_INFO_CONTENT
    "#define __COMPILE_DATE \"${CP2K_BUILD_TIMESTAMP}\"\n#define __COMPILE_REVISION \"${CP2K_GIT_HASH}\"\n"
)

if(EXISTS "${OUTPUT_FILE}")
  file(READ "${OUTPUT_FILE}" CP2K_OLD_BUILD_INFO_CONTENT)
else()
  set(CP2K_OLD_BUILD_INFO_CONTENT "")
endif()

if(NOT CP2K_BUILD_INFO_CONTENT STREQUAL CP2K_OLD_BUILD_INFO_CONTENT)
  file(WRITE "${OUTPUT_FILE}" "${CP2K_BUILD_INFO_CONTENT}")
endif()
