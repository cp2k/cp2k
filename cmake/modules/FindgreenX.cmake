#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(GREENX "greenx")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_GREENX IMPORTED_TARGET GLOBAL greenx)
endif()

#if(NOT CP2K_GREENX_FOUND)
#  cp2k_find_libraries(GREENX GXCommon)
# message("${CP2K_GREENX_INCLUDE_DIRS} ${CP2K_GREENX_LINK_LIBRARIES}")
#  cp2k_find_libraries(GREENX LibGXAC)
# message("${CP2K_GREENX_INCLUDE_DIRS} ${CP2K_GREENX_LINK_LIBRARIES}")
#  cp2k_find_libraries(GREENX LibGXMiniMax)
# message("${CP2K_GREENX_INCLUDE_DIRS} ${CP2K_GREENX_LINK_LIBRARIES}")

#endif()

if(NOT CP2K_GREENX_FOUND)
  set(CP2K_GREENX_LINK_LIBRARIES "")

  cp2k_find_libraries(GREENX GXCommon)
  set(GXCommon_LIBRARIES ${CP2K_GREENX_LINK_LIBRARIES})

  cp2k_find_libraries(GREENX LibGXAC)
  set(LibGXAC_LIBRARIES ${CP2K_GREENX_LINK_LIBRARIES})

  cp2k_find_libraries(GREENX LibGXMiniMax)
  set(LibGXMiniMax_LIBRARIES ${CP2K_GREENX_LINK_LIBRARIES})

  set(CP2K_GREENX_LINK_LIBRARIES "${GXCommon_LIBRARIES};${LibGXAC_LIBRARIES};${LibGXMiniMax_LIBRARIES}")

  message(STATUS "Final CP2K_GREENX_LINK_LIBRARIES: ${CP2K_GREENX_LINK_LIBRARIES}")

#  cp2k_include_dirs(GREENX "modules/gx_minimax.mod")
  message(STATUS "Final CP2K_GREENX_INCLUDE_DIRS: ${CP2K_GREENX_INCLUDE_DIRS}")
endif()

find_package_handle_standard_args(
  greenX DEFAULT_MSG CP2K_GREENX_INCLUDE_DIRS CP2K_GREENX_LINK_LIBRARIES)

 message("${CP2K_GREENX_INCLUDE_DIRS} ${CP2K_GREENX_LINK_LIBRARIES}")

if(CP2K_GREENX_FOUND)
  message(STATUS "Here" )
  if(NOT TARGET cp2k::greenx::greenx)
    add_library(cp2k::greenx::greenx INTERFACE IMPORTED)
  endif()

  set_target_properties(
    cp2k::greenx::greenx PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                    "${CP2K_GREENX_INCLUDE_DIRS}")

  target_link_libraries(cp2k::greenx::greenx
                        INTERFACE ${CP2K_GREENX_LINK_LIBRARIES})
endif()

mark_as_advanced(CP2K_GREENX_FOUND CP2K_GREENX_LINK_LIBRARIES
                 CP2K_GREENX_INCLUDE_DIRS)
