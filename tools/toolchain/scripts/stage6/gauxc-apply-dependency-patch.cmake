#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# FetchContent may repeat its patch step during CMake regeneration.
find_program(PATCH_EXECUTABLE patch REQUIRED)
execute_process(
  COMMAND "${PATCH_EXECUTABLE}" --batch --forward --dry-run -p1 -i
          "${PATCH_FILE}"
  RESULT_VARIABLE forward_status
  OUTPUT_QUIET ERROR_QUIET)
if(forward_status EQUAL 0)
  execute_process(COMMAND "${PATCH_EXECUTABLE}" --batch --forward -p1 -i
                          "${PATCH_FILE}" COMMAND_ERROR_IS_FATAL ANY)
else()
  execute_process(
    COMMAND "${PATCH_EXECUTABLE}" --batch --reverse --dry-run -p1 -i
            "${PATCH_FILE}"
    RESULT_VARIABLE reverse_status
    OUTPUT_QUIET ERROR_QUIET)
  if(NOT reverse_status EQUAL 0)
    message(
      FATAL_ERROR
        "Dependency patch neither applies cleanly nor is already applied: ${PATCH_FILE}"
    )
  endif()
endif()
