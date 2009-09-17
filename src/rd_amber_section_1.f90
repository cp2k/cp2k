    CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
    i = 1
    DO WHILE ((i<=dim).AND.(.NOT.my_end))
       IF (parser_test_next_token(parser,error=error)=="EOL")&
            CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
       IF (my_end) EXIT
#if defined (LOWER_TO_UPPER)
       CALL parser_get_object  (parser,array1(i),lower_to_upper=.TRUE.,error=error)
#else
       CALL parser_get_object  (parser,array1(i),error=error)
#endif
       i = i + 1
    END DO
    ! Trigger end of file aborting
    CALL cp_assert(.NOT.my_end.OR.(i>dim), cp_fatal_level, cp_assertion_failed, routineP,&
         "End of file while reading section "//TRIM(section)//" in amber topology file!"//&
CPSourceFileRef,&
         error=error,failure=failure)
