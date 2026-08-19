!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM twopt
   USE trajana_command_line,            ONLY: has_flag
   USE trajana_twopt_analysis,          ONLY: print_twopt_help,&
                                              run_twopt

   IMPLICIT NONE

   IF (COMMAND_ARGUMENT_COUNT() == 0 .OR. has_flag("--help") .OR. has_flag("-h")) THEN
      CALL print_twopt_help()
   ELSE
      CALL run_twopt()
   END IF
END PROGRAM twopt
