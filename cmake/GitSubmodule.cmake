#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Call to ensure that the git submodule in location `path` is loaded. If the
# submodule is not loaded, an error message that describes how to update the
# submodules is printed. Sets the variable name_avail to `ON` if the submodule
# is available, or `OFF` otherwise. copyright github.com/arbor-sim

function(check_git_submodule name path)
  set(success_var "${name}_avail")
  set(${success_var}
      ON
      PARENT_SCOPE)

  get_filename_component(dotgit "${path}/.git" ABSOLUTE)
  if(NOT EXISTS ${dotgit})
    message(
      FATAL_ERROR
        "\nThe git submodule for ${name} is not available.\n"
        "To check out all submodules use the following commands:\n"
        "    git submodule init\n"
        "    git submodule update\n"
        "Or download submodules recursively when checking out:\n"
        "    git clone --recursive https://github.com/cp2k/cp2k.git\n")
  endif()
endfunction()
