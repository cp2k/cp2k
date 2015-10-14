! Common use statements and preprocessor macros
! should be included in the use statements

  USE cp_log_handling,                 ONLY: cp_logger_type,&
                                             cp_logger_p_type,&
                                             cp_logger_create,&
                                             cp_logger_release,&
                                             cp_logger_set,&
                                             cp_warning_level,&
                                             cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_logger_get_default_io_unit,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_default_logger_stack_size,&
                                             cp_logger_get_unit_nr,&
                                             cp_logger_generate_filename,&
                                             cp_get_default_logger,&
                                             cp_logger_would_log,&
                                             cp_add_default_logger,&
                                             cp_rm_default_logger,&
                                             cp_to_string
  USE cp_error_handling,               ONLY: cp_abort,&
                                             cp_warn,&
                                             cp_assertion_failed,&
                                             cp_internal_error,&
                                             cp_assert,&
                                             cp_precondition_failed,&
                                             cp_caller_error,&
                                             cp_wrong_args_error,&
                                             cp__a,cp__b,cp__w,cp__l

! Dangerous: Full path can be arbitrarily long and might overflow Fortran line.
#if !defined(__SHORT_FILE__)
#define __SHORT_FILE__ __FILE__
#endif

#define __LOCATION__ cp__l(__SHORT_FILE__,__LINE__)
#define CPWARN(msg) CALL cp__w(__SHORT_FILE__,__LINE__,msg)
#define CPABORT(msg) CALL cp__b(__SHORT_FILE__,__LINE__,msg)
#define CPASSERT(cond) IF(.NOT.(cond))CALL cp__a(__SHORT_FILE__,__LINE__)

! The MARK_USED macro can be used to mark an argument/variable as used.
! It is intended to make it possible to switch on -Werror=unused-dummy-argument,
! but deal elegantly with e.g. library wrapper routines that take arguments only used if the library is linked in. 
! This code should be valid for any Fortran variable, is always standard conforming,
! and will be optimized away completely by the compiler
#define MARK_USED(foo) IF(.FALSE.)THEN; DO ; IF(SIZE(SHAPE(foo))==-1) EXIT ;  END DO ; ENDIF
