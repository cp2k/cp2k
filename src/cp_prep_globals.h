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
#define CPSourceFileRef __FILE__//' line '//cp_to_string(__LINE__)

! if the following macro is defined the shortest possible form of macro 
! expansions is used (but then the error messages are not so meaningful)
#ifdef  __INTEL

#define FD_SHORT_EXPANSIONS

#endif

#ifndef FD_SHORT_EXPANSIONS

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

#else
! same as before but trying to use a small amount of characters
! the test is not inlined (you have always a function call)

#define CPPrecondition(cond,level,routineP,error,failure) \
call cp_assert(cond,level,200,routineP,'PRECOND',error,failure)
#define CPPostcondition(cond,level,routineP,error,failure) \
call cp_assert(cond,level,-200,routineP,'POSTCOND',error,failure)
#define CPInvariant(cond,level,routineP,error,failure) \
call cp_assert(cond, level, -100,routineP,'INVARIANT',error,failure)
#define CPAssert(cond,level,routineP,error,failure) \
call cp_assert(cond,level,-300,routineP,'ASSERT',error,failure)
#define CPErrorMessage(level,routineP,msg,error) \
call cp_error_message(level,routineP,msg,error)
#define CPPreconditionNoFail(cond,level,routineP,error) \
call cp_assert(cond,level,200,routineP,'PRECOND',error)
#define CPPostconditionNoFail(cond,level,routineP,error) \
call cp_assert(cond,level,-200,routineP,'POSTCOND',error)
#define CPInvariantNoFail(cond,level,routineP,error) \
call cp_assert(cond, level, -100,routineP,'INVAR',error)
#define CPAssertNoFail(cond,level,routineP,error) \
call cp_assert(cond,level,-300,routineP,'ASSERT',error)
#define CPPreconditionNoErr(cond, level, routine_name) \
call cp_assert(cond, level, 200, routine_name ,"PRECOND")
#define CPPostconditionNoErr(cond, level, routine_name) \
call cp_assert(cond, level, -200, routine_name ,"POSTCOND")
#define CPInvariantNoErr(cond,level,routineP) \
call cp_assert(cond, level, -100,routineP,'INVAR')
#define CPAssertNoErr(cond,level,routineP) \
call cp_assert(.false.,level,-300,routineP,'ASSERT')

#endif


#endif
