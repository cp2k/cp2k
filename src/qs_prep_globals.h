#ifndef QS_PREP_GLOBALS_H
#define QS_PREP_GLOBALS_H

! The following macros are here to facilitate the use of error handling 
! proposed in qs_error_handling.
! they assume at least 'use qs_error_handling, only: qs_assert, qs_error_type'
! and 'use qs_log_handling, only: qs_to_string'
! They ere useful because they give a reference to the file and line
! number in the error message.

! this macro expands to a string that contains the filename.
! if the path is long the filename could make that some lines of code
! become too long and overlow (fortran compilers have a maximum line length)
! in this case substitute __FILE__ with "file" down here.
! obviously then the error messages will not give the filename.
! (otherwise make the file reference in the makefile relative vs. absolute)
#define QSSourceFileRef() __FILE__//' line '//qs_to_string(__LINE__)

! #define FD_ASSERT_WITH_FILE_NAME

#ifdef FD_ASSERT_WITH_FILE_NAME
! try to also print the file path (can give problems if the path is 
! too long)

! (qs_precondition_failed = 200)
#define QSPrecondition(cond,level,routineN,error,failure) \
call qs_assert(cond,level,200,\
routineN,'PRECONDITION(cond) failed in '//QSSourceFileRef(),error,failure)

! (qs_postcondition_failed = -200)
#define QSPostcondition(cond,level,routineN,error,failure) \
call qs_assert(cond,level,-200, \
routineN,'POSTCONDITION(cond) failed in '//QSSourceFileRef(),error,failure)

! (qs_invariant_failed = -100)
#define QSInvariant(cond,level,routineN,error,failure) \
call qs_assert(cond, level, -100, \
routineN,'INVARIANT (cond) failed in '//QSSourceFileRef(),error,failure)

! (qs_assert_failed = -300)
#define QSAssert(cond,level,routineN,error,failure) \
call qs_assert(cond,level,-300,\
routineN,'ASSERTION (cond) failed in '//QSSourceFileRef(),error,failure)

#define QSErrorMessage(level,routineN,msg,error) \
call qs_error_message(level,routineN,msg//' in file'//QSSourceFileRef(),error)

#else
! inlines the test but does not write the file name, but that should
! be easily recovered from the routineN variable (that should contain
! the module name).

! (qs_precondition_failed = 200)
#define QSPrecondition(cond,level,routineN,error,failure) \
if(cond) call qs_assert(.true.,level,200,routineN,\
'PRECONDITION(cond) failed at line '//qs_to_string(__LINE__),error,failure)

! (qs_postcondition_failed = -200)
#define QSPostcondition(cond,level,routineN,error,failure) \
if(cond) call qs_assert(.true.,level,-200,routineN,\
'POSTCONDITION(cond) failed at line '//qs_to_string(__LINE__),error,failure)

! (qs_invariant_failed = -100)
#define QSInvariant(cond,level,routineN,error,failure) \
if(cond) call qs_assert(.true., level, -100,routineN,\
'INVARIANT (cond) failed at line '//qs_to_string(__LINE__),error,failure)

! (qs_assert_failed = -300)
#define QSAssert(cond,level,routineN,error,failure) \
if(cond) call qs_assert(.true.,level,-300,routineN,\
'ASSERTION (cond) failed at line '//qs_to_string(__LINE__),error,failure)

#define QSErrorMessage(level,routineN,msg,error) \
call qs_error_message(level,routineN,msg//' at line '//qs_to_string(__LINE__),error)

#endif


#endif
