#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

if(NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "OUTPUT_FILE is not defined")
endif()

string(TIMESTAMP CP2K_BUILD_TIMESTAMP "%Y-%m-%d %H:%M:%S")

set(CP2K_BUILD_INFO_CONTENT
    "#define __COMPILE_DATE \"${CP2K_BUILD_TIMESTAMP}\"\n")

if(EXISTS "${OUTPUT_FILE}")
  file(READ "${OUTPUT_FILE}" CP2K_OLD_BUILD_INFO_CONTENT)
else()
  set(CP2K_OLD_BUILD_INFO_CONTENT "")
endif()

if(NOT CP2K_BUILD_INFO_CONTENT STREQUAL CP2K_OLD_BUILD_INFO_CONTENT)
  file(WRITE "${OUTPUT_FILE}" "${CP2K_BUILD_INFO_CONTENT}")
endif()
