! Basic use statements and preprocessor macros
! should be included in the use statements

  USE base_hooks,                      ONLY: cp__a,&
                                             cp__b,&
                                             cp__w,&
                                             cp__l,&
                                             cp_abort,&
                                             cp_warn,&
                                             timeset,&
                                             timestop


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
