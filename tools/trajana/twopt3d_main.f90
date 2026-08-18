PROGRAM twopt3d
   USE trajana_command_line,            ONLY: has_flag
   USE trajana_twopt3d_analysis,        ONLY: print_twopt3d_help,&
                                              run_twopt3d

   IMPLICIT NONE

   IF (COMMAND_ARGUMENT_COUNT() == 0 .OR. has_flag("--help") .OR. has_flag("-h")) THEN
      CALL print_twopt3d_help()
   ELSE
      CALL run_twopt3d()
   END IF
END PROGRAM twopt3d
