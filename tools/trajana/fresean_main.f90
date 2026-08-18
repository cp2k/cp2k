PROGRAM fresean
   USE trajana_command_line,            ONLY: has_flag
   USE trajana_fresean_analysis,        ONLY: print_fresean_help,&
                                              run_fresean

   IMPLICIT NONE

   IF (has_flag("--help") .OR. has_flag("-h")) THEN
      CALL print_fresean_help()
   ELSE
      CALL run_fresean()
   END IF
END PROGRAM fresean
