#ifndef CP_PREP_GLOBALS_H
#define CP_PREP_GLOBALS_H

! The following macros are here to facilitate the use of error handling 
! proposed in cp_error_handling.
! they assume at least 'use cp_error_handling, only: cp_assert, cp_error_type'
! and 'use cp_log_handling, only: cp_to_string'
! They ere useful because they give a reference to the file and line
! number in the error message.

! this macro expands to a string that contains the filename.
! if the path is long the filename could make that some lines of code
! become too long and overlow (fortran compilers have a maximum line length)
! in this case substitute __FILE__ with "file" down here.
! obviously then the error messages will not give the filename.
! (otherwise make the file reference in the makefile relative vs. absolute)
#define CPSourceFileRef() __FILE__//' line '//cp_to_string(__LINE__)

! #define FD_ASSERT_WITH_FILE_NAME

#ifdef FD_ASSERT_WITH_FILE_NAME
! try to also print the file path (can give problems if the path is 
! too long)

! (cp_precondition_failed = 200)
#define CPPrecondition(cond,level,routineP,error,failure) \
call cp_assert(cond,level,200,\
routineP,'PRECONDITION(cond) failed in '//CPSourceFileRef(),error,failure)

! (cp_postcondition_failed = -200)
#define CPPostcondition(cond,level,routineP,error,failure) \
call cp_assert(cond,level,-200, \
routineP,'POSTCONDITION(cond) failed in '//CPSourceFileRef(),error,failure)

! (cp_invariant_failed = -100)
#define CPInvariant(cond,level,routineP,error,failure) \
call cp_assert(cond, level, -100, \
routineP,'INVARIANT (cond) failed in '//CPSourceFileRef(),error,failure)

! (cp_assert_failed = -300)
#define CPAssert(cond,level,routineP,error,failure) \
call cp_assert(cond,level,-300,\
routineP,'ASSERTION (cond) failed in '//CPSourceFileRef(),error,failure)

#define CPErrorMessage(level,routineP,msg,error) \
call cp_error_message(level,routineP,msg//' in file'//CPSourceFileRef(),error)

#else
! inlines the test but does not write the file name, but that should
! be easily recovered from the routineP variable (that should contain
! the module name).

! (cp_precondition_failed = 200)
#define CPPrecondition(cond,level,routineP,error,failure) \
if(.not.(cond)) call cp_assert(.false.,level,200,routineP,\
'PRECONDITION(cond) failed at line '//cp_to_string(__LINE__),error,failure)

! (cp_postcondition_failed = -200)
#define CPPostcondition(cond,level,routineP,error,failure) \
if(.not.(cond)) call cp_assert(.false.,level,-200,routineP,\
'POSTCONDITION(cond) failed at line '//cp_to_string(__LINE__),error,failure)

! (cp_invariant_failed = -100)
#define CPInvariant(cond,level,routineP,error,failure) \
if(.not.(cond)) call cp_assert(.false., level, -100,routineP,\
'INVARIANT (cond) failed at line '//cp_to_string(__LINE__),error,failure)

! (cp_assert_failed = -300)
#define CPAssert(cond,level,routineP,error,failure) \
if(.not.(cond)) call cp_assert(.false.,level,-300,routineP,\
'ASSERTION (cond) failed at line '//cp_to_string(__LINE__),error,failure)

#define CPErrorMessage(level,routineP,msg,error) \
call cp_error_message(level,routineP,msg//' at line '//cp_to_string(__LINE__),error)

#endif

#define CPPreconditionNoFail(cond,level,routineP,error) \
if(.not.(cond)) call cp_assert(.false.,level,200,routineP,\
'PRECONDITION(cond) failed at line '//cp_to_string(__LINE__),error)

! (cp_postcondition_failed = -200)
#define CPPostconditionNoFail(cond,level,routineP,error) \
if(.not.(cond)) call cp_assert(.false.,level,-200,routineP,\
'POSTCONDITION(cond) failed at line '//cp_to_string(__LINE__),error)

! (cp_invariant_failed = -100)
#define CPInvariantNoFail(cond,level,routineP,error) \
if(.not.(cond)) call cp_assert(.false., level, -100,routineP,\
'INVARIANT (cond) failed at line '//cp_to_string(__LINE__),error)

! (cp_assert_failed = -300)
#define CPAssertNoFail(cond,level,routineP,error) \
if(.not.(cond)) call cp_assert(.false.,level,-300,routineP,\
'ASSERTION (cond) failed at line '//cp_to_string(__LINE__),error)


#define CPPreconditionNoErr(cond, level, routine_name) \
if(.not.(cond)) call cp_assert(.false., level, 200, routine_name ,\
"PRECONDITION failed on line " //cp_to_string(__LINE__))

#define CPPostconditionNoErr(cond, level, routine_name) \
if(.not.(cond)) call cp_assert(.false., level, -200, routine_name ,\
"POSTCONDITION failed on line " //cp_to_string(__LINE__))

#define CPInvariantNoErr(cond,level,routineP) \
if(.not.(cond)) call cp_assert(.false., level, -100,routineP,\
'INVARIANT (cond) failed at line '//cp_to_string(__LINE__))

! (cp_assert_failed = -300)
#define CPAssertNoErr(cond,level,routineP) \
if(.not.(cond)) call cp_assert(.false.,level,-300,routineP,\
'ASSERTION (cond) failed at line '//cp_to_string(__LINE__))


#endif
