!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_rho_set_types [1.0] *
!!
!!   NAME
!!     xc_rho_set_types
!!
!!   FUNCTION
!!     contains the structure 
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE xc_rho_set_types
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l, cp_error_check
  USE pw_pool_types, ONLY: pw_pool_type,pw_pool_create_pw, &
       pw_pool_give_back_cr3d, pw_pool_create_cr3d, pw_pool_give_back_pw
  USE dft_types, ONLY: dft_control_type
  USE simulation_cell, only: cell_type
  USE cp_array_r_utils, only: cp_3d_r_p_type
  use kinds, only: dp
  use timings, only: timeset, timestop
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='xc_rho_set_types'
  INTEGER, SAVE :: last_rho_set_id=0

!!***
!****************************************************************************

!!****s* xc_rho_set_types/xc_rho_cflags_type [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_type
!!
!!   FUNCTION
!!     contains a flag for each component of xc_rho_set, so that you can
!!     use it to tell which components you need, which ones you need,....
!!
!!   NOTES
!!     low_level type without retain/release
!!
!!   ATTRIBUTES
!!     - rho: flags for rho (the total density)
!!     - rho_spin: flag for rhoa and rhob (the spin density with LSD)
!!     - drho: flag for drho (the gradient of rho)
!!     - drho_spin: flag for drhoa and drhob (the gradient of the spin
!!       density)
!!     - norm_drho: flag for norm_drho (the norm of the gradient of rho)
!!     - norm_drho_spin: flag for norm_drhoa, norm_drhob (the norm of the
!!       gradient of the spin density)
!!     - drhoa_drhob: flag for drhoa_drhob (the scalar product of the
!!       gradient of the two spin densities)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE xc_rho_cflags_type
     LOGICAL :: rho, rho_spin, drho, drho_spin,&
          norm_drho, norm_drho_spin, drhoa_drhob
  END TYPE xc_rho_cflags_type
!!***
!****************************************************************************


!!****s* xc_rho_set_types/xc_rho_set_type [1.0] *
!!
!!   NAME
!!     xc_rho_set_type
!!
!!   FUNCTION
!!     represent a density, with all the representation and data needed
!!     to perform a functional evaluation
!!
!!   NOTES
!!     the use of 3d arrays is the result of trying to use only basic
!!     types (to be generic and independent from the method), and
!!     avoiding copies using the actual structure.
!!     only the part defined by local bounds is guaranteed to be present,
!!     and it is the only meaningful part.
!!
!!     add rho_1_3 (rho^(1./3.))?
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - id_nr: identification number (unique)
!!     - local_bounds: the part of the 3d array on which the functional is 
!!       computed
!!     - owns: which components are owned by this structure (and should
!!       be deallocated
!!     - has: which components are present and up to date
!!     - rho: the density
!!     - drho: the gradient of the density (x,y and z direction)
!!     - norm_drho: the norm of the gradient of the density
!!     - rhoa, rhob: spin alpha and beta parts of the density in the LSD case
!!     - drhoa, drhob: gradient of the spin alpha and beta parts of the density
!!       in the LSD case (x,y and z direction)
!!     - norm_drhoa, norm_drhob: norm of the gradient of rhoa and rhob
!!     - drhoa_drhob: the scalar product between the gradient of rhoa and the
!!       one of rhob
!!
!!   AUTHOR
!!     fawzi & thomas
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi & thomas]
!!
!!   SOURCE
!***************************************************************************
  TYPE xc_rho_set_type
     INTEGER :: ref_count, id_nr
     INTEGER, DIMENSION(2,3) :: local_bounds
     TYPE(xc_rho_cflags_type) :: owns, has
! for spin restricted systems
     real(KIND = dp), dimension(:,:,:), pointer :: rho
     TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drho
     real(KIND = dp), dimension(:,:,:), pointer :: norm_drho
! for UNrestricted systems
     real(KIND = dp), dimension(:,:,:), pointer :: rhoa, rhob
     TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drhoa, drhob
     real(KIND = dp), dimension(:,:,:), pointer :: norm_drhoa, norm_drhob
     real(KIND = dp), dimension(:,:,:), pointer :: drhoa_drhob
  END TYPE xc_rho_set_type
!!***
!****************************************************************************

contains

!!****f* xc_rho_set_types/xc_rho_cflags_setall [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_setall
!!
!!   FUNCTION
!!     sets all the flags to the given value
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cflags: the flags to set
!!     - value: the value to set (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_cflags_setall(cflags, value, error)
    TYPE(xc_rho_cflags_type), intent(out) :: cflags
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure, my_value
    character(len=*), parameter :: routineN='xc_rho_cflags_setall',&
         routineP=moduleN//':'//routineN

    failure=.false.
    my_value=.false.
    if (present(value)) my_value=value
    cflags%rho=my_value
    cflags%rho_spin=my_value
    cflags%drho=my_value
    cflags%drho_spin=my_value
    cflags%norm_drho=my_value
    cflags%norm_drho_spin=my_value
    cflags%drhoa_drhob=my_value
  END SUBROUTINE xc_rho_cflags_setall
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_or [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_or
!!
!!   FUNCTION
!!     performs an or on all the flags (cflags1=cflag1.or.cflags2)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cflags1: the flags to change
!!     - cflags2: the flags to or with
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_cflags_or(cflags1, cflags2, error)
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), intent(in) :: cflags2
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure, my_value
    character(len=*), parameter :: routineN='xc_rho_cflags_or',&
         routineP=moduleN//':'//routineN

    failure=.false.

    cflags1%rho=cflags1%rhoc.or.flags2%rho
    cflags1%rho_spin=cflags1%rho_spin.or.cflags2%rho_spin
    cflags1%drho=cflags1%drho.or.cflags2%drho
    cflags1%drho_spin=cflags1%drho_spin.or.cflags2%drho_spin
    cflags1%norm_drho=cflags1%norm_drho.or.cflags2%norm_drho
    cflags1%norm_drho_spin=cflags1%norm_drho_spin.or.cflags2%norm_drho_spin
    cflags1%drhoa_drhob=cflags1%drhoa_drhob.or.cflags2%drhoa_drhob

  END SUBROUTINE xc_rho_cflags_or
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_and [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_and
!!
!!   FUNCTION
!!     performs an and on all the flags (cflags1=cflag1.and.cflags2)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cflags1: the flags to change
!!     - cflags2: the flags to and with
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_cflags_and(cflags1, cflags2, error)
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), intent(in) :: cflags2
    type(cp_error_type), optional, intent(inout) :: error

    LOGICAL :: failure, my_value
    character(len=*), parameter :: routineN='xc_rho_cflags_and',&
         routineP=moduleN//':'//routineN

    failure=.false.

    cflags1%rho=cflags1%rhoc.and.flags2%rho
    cflags1%rho_spin=cflags1%rho_spin.and.cflags2%rho_spin
    cflags1%drho=cflags1%drho.and.cflags2%drho
    cflags1%drho_spin=cflags1%drho_spin.and.cflags2%drho_spin
    cflags1%norm_drho=cflags1%norm_drho.and.cflags2%norm_drho
    cflags1%norm_drho_spin=cflags1%norm_drho_spin.and.cflags2%norm_drho_spin
    cflags1%drhoa_drhob=cflags1%drhoa_drhob.and.cflags2%drhoa_drhob

  END SUBROUTINE xc_rho_cflags_and
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_equal [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_equal
!!
!!   FUNCTION
!!     return true if the two cflags are equal
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - cflags1, cflags2: the flags to compare
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION xc_rho_cflags_equal(cflags1, cflags2, error) RESULT(equal)
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in) :: cflags2
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    logical :: equal

    LOGICAL :: failure
    character(len=*), parameter :: routineN='xc_rho_cflags_equal',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    equal=(cflags1%rho==flags2%rho.AND.&
         cflags1%rho_spin==cflags2%rho_spin.AND.&
         cflags1%drho==cflags2%drho.AND.&
         cflags1%drho_spin==cflags2%drho_spin.AND.&
         cflags1%norm_drho==cflags2%norm_drho.AND.&
         cflags1%norm_drho_spin==cflags2%norm_drho_spin.AND.&
         cflags1%drhoa_drhob==cflags2%drhoa_drhob)

  END FUNCTION xc_rho_cflags_equal
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_cflags_contain [1.0] *
!!
!!   NAME
!!     xc_rho_cflags_contain
!!
!!   FUNCTION
!!     return true if all the true flags in cflags2 are also true in cflags1
!!
!!   NOTES
!!     defines a partial ordering on the cflags (cflags1>cflags2)
!!
!!   ARGUMENTS
!!     - cflags1, cflags2: the flags to compare
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION xc_rho_cflags_contain(cflags1, cflags2, error) RESULT(contain)
    TYPE(xc_rho_cflags_type), INTENT(inout) :: cflags1
    TYPE(xc_rho_cflags_type), INTENT(in) :: cflags2
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    logical :: contain

    LOGICAL :: failure
    character(len=*), parameter :: routineN='xc_rho_cflags_contain',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    contain=((cflags1%rho.or..not.flags2%rho).AND.&
         (cflags1%rho_spin.or..not.cflags2%rho_spin).AND.&
         (cflags1%drho.OR..not.cflags2%drho).AND.&
         (cflags1%drho_spin.or..not.cflags2%drho_spin).AND.&
         (cflags1%norm_drho.or..not.cflags2%norm_drho).AND.&
         (cflags1%norm_drho_spin.or..not.cflags2%norm_drho_spin).AND.&
         (cflags1%drhoa_drhob.or..not.cflags2%drhoa_drhob))

  END FUNCTION xc_rho_cflags_contain
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_create [1.0] *
!!
!!   NAME
!!     xc_rho_set_create
!!
!!   FUNCTION
!!     allocates and does (minimal) initialization of a rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the structure to allocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_create(rho_set, error)
    TYPE(xc_rho_set_type), pointer :: rho_set
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='xc_rho_set_create',&
         routineP=moduleN//':'//routineN
    integer :: stat

    failure=.false.

    ALLOCATE(rho_set, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       rho_set%ref_count=1
       last_rho_set_id=last_rho_set_id+1
       rho_set%id_nr=last_rho_set_id
       rho_set%local_bounds(1,:)=-HUGE(0) ! we want to crash...
       rho_set%local_bounds(1,:)=HUGE(0)
       CALL xc_rho_cflags_setall(rho_set%owns,.TRUE.,error=error)
       CALL xc_rho_cflags_setall(rho_set%has,.FALSE.,error=error)
       NULLIFY(rho_set%rho)
       DO i=1,3
          NULLIFY(rho_set%drho(i)%array)
       END DO
       NULLIFY(rho_set%norm_drho,rho_set%rhoa,rho_set%rhob)
       DO i=1,3
          nullify(rho_set%drhoa(i)%array, rho_set%drhob(i)%array)
       END DO
       NULLIFY(rho_set%norm_drhoa, rho_set%norm_drhob, &
            rho_set%drhoa_drhob)
    END IF
  END SUBROUTINE xc_rho_set_create
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_retain [1.0] *
!!
!!   NAME
!!     xc_rho_set_retain
!!
!!   FUNCTION
!!     retains the given rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_retain(rho_set, error)
    TYPE(xc_rho_set_type), pointer :: rho_set
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='xc_rho_set_retain',&
         routineP=moduleN//':'//routineN

    failure=.false.
    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(rho_set%ref_count>0,cp_failure_level,routineP,error)
       rho_set%ref_count=rho_set%ref_count+1
    END IF
  END SUBROUTINE xc_rho_set_retain
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_release [1.0] *
!!
!!   NAME
!!     xc_rho_set_release
!!
!!   FUNCTION
!!     releases the given rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the structure to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_release(rho_set, local_bounds, rho, drho, norm_drho,&
       error)
    TYPE(xc_rho_set_type), pointer :: rho_set
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='xc_rho_set_release',&
         routineP=moduleN//':'//routineN
    integer :: stat

    failure=.false.

    IF (ASSOCIATED(rho_set)) THEN
       CPPreconditionNoFail(rho_set%ref_count>0,cp_failure_level,routineP,error,failure)
       rho_set%ref_count=rho_set%ref_count-1
       IF (rho_set%ref_count==0) THEN
          rho_set%local_bounds(1,:)=-HUGE(0) ! we want to crash...
          rho_set%local_bounds(1,:)=HUGE(0)
          IF (rho_set%owns%rho .AND. ASSOCIATED(rho_set%rho)) THEN
             DEALLOCATE(rho_set%rho, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%rho_spin) THEN
             IF (ASSOCIATED(rho_set%rhoa)) THEN
                DEALLOCATE(rho_set%rhoa, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%rhob)) THEN
                DEALLOCATE(rho_set%rhob, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
          IF (rho_set%owns%drho.AND.ASSOCIATED(rho_set%drho)) THEN
             DO i=1,3
                IF (ASSOCIATED(rho_set%drho(i)%array)) THEN
                   DEALLOCATE(rho_set%drho(i)%array, stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                END IF
             END DO
          END IF
          IF (rho_set%owns%drho_spin) THEN
             DO i=1,3
                IF (ASSOCIATED(rho_set%drhoa(i)%array)) THEN
                   DEALLOCATE(rho_set%drhoa(i)%array, stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                END IF
                IF (ASSOCIATED(rho_set%drhob(i)%array)) THEN
                   DEALLOCATE(rho_set%drhob(i)%array, stat=stat)
                   CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                END IF
             END DO
          END IF
          IF (rho_set%owns%norm_drho.and.associated(rho_set%norm_drho)) THEN
             DEALLOCATE(rho_set%norm_drho, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (rho_set%owns%norm_drho_spin) THEN
             IF (ASSOCIATED(rho_set%norm_drhoa)) THEN
                DEALLOCATE(rho_set%norm_drhoa, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
             IF (ASSOCIATED(rho_set%norm_drhob)) THEN
                DEALLOCATE(rho_set%norm_drhob, stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
          IF (rho_set%owns%drhoa_drhob.AND.ASSOCIATED(rho_set%drhoa_drhob)) THEN
             DEALLOCATE(rho_set%drhoa_drhob, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(rho_set, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    nullify(rho_set)
  END SUBROUTINE xc_rho_set_release
!***************************************************************************

!!****f* xc_rho_set_types/xc_rho_set_get [1.0] *
!!
!!   NAME
!!     xc_rho_set_get
!!
!!   FUNCTION
!!     returns the various attributes of rho_set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rho_set: the object you whant info about
!!     - can_return_null: if true the object returned can be null,
!!       if false (the default) it stops with an error if a requested
!!       component is not associated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     for the other components see the attributes of xc_rho_set
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_get(rho_set, can_return_null, rho, drho, norm_drho,&
       rhoa, rhob, norm_drhoa, norm_drhob, drhoa_drhob,error)
    TYPE(xc_rho_set_type), POINTER :: rho
    LOGICAL, INTENT(in), optional :: can_return_null
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: rho
    TYPE(cp_3d_r_p_type), DIMENSION(:), OPTIONAL:: drho
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: norm_drho
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: rhoa, rhob
    TYPE(cp_3d_r_p_type), DIMENSION(:), OPTIONAL :: drhoa, drhob
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: norm_drhoa, norm_drhob
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: drhoa_drhob
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, my_can_return_null, valid
    CHARACTER(len=*), PARAMETER :: routineN='xc_rho_set_get',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.
    my_can_return_null=.FALSE.
    IF (PRESENT(can_return_null)) my_can_return_null=can_return_null

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_cont>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(rho)) THEN
          rho => rho_set%rho
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(drho)) THEN
          drho => rho_set%drho
          DO i=1,3
             IF (.NOT.ASSOCIATED(rho_set%drho(i)%array)) THEN
                nullify(drho)
             END IF
          END DO
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(drho),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(norm_drho)) THEN
          norm_drho => rho_set%norm_drho
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(norm_drho),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rhoa)) THEN
          rhoa => rho_set%rhoa
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rhoa),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rhob)) THEN
          rhob => rho_set%rhob
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rhob),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(drhoa)) THEN
          drhoa => rho_set%drho
          DO i=1,3
             IF (.NOT.ASSOCIATED(rho_set%drhoa(i)%array)) THEN
                nullify(drhoa)
             END IF
          END DO
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(drhoa),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(drhob)) THEN
          drho => rho_set%drhob
          DO i=1,3
             IF (.NOT.ASSOCIATED(rho_set%drhob(i)%array)) THEN
                nullify(drhob)
             END IF
          END DO
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(drhob),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(norm_drhoa)) THEN
          norm_drhoa => rho_set%norm_drhoa
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(norm_drhoa),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(norm_drhob)) THEN
          norm_drhob => rho_set%norm_drhob
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(norm_drhob),cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(rhoa_rhob)) THEN
          rhoa_rhob => rho_set%rhoa_rhob
          CPPrecondition(my_can_return_null.OR.ASSOCIATED(rhoa_rhob),cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE xc_rho_set_get
!***************************************************************************


!!****f* xc_rho_set_types/xc_rho_set_update [1.0] *
!!
!!   NAME
!!     xc_rho_set_update
!!
!!   FUNCTION
!!     updates the given rho set with the density given by
!!     rho_r (and rho_g). The rho set will contain the components specified
!!     in needs
!!
!!   NOTES
!!     when 
!!
!!   ARGUMENTS
!!     - rho_set: the rho_set to update
!!     - rho_r: the new density (in r space)
!!     - rho_g: the new density (in g space, needed for smoothing and
!!       derivatives)
!!     - needs: the components of rho that are needed
!!     - dft_control: parameters (LSD/LDA, kind of derivative requested,...)
!!     - cell: cell parameters (to rescale the derivative)
!!     - pw_pool: pool for the allocation of pw and cr3d
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_rho_set_update(rho_set, rho_r, rho_g, needs,&
       dft_control,cell,pw_pool,error)
    TYPE(xc_rho_set_type), pointer :: rho_set
    TYPE(pw_p_type),DIMENSION(:), INTENT(inout), TARGET :: rho_r
    TYPE(pw_p_type),DIMENSION(:), POINTER        :: rho_g
    TYPE(xc_rho_cflags_type), INTENT(in) :: needs
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(cell_type), pointer :: cell
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, my_rho_r_local, my_rho_g_local
    character(len=*), parameter :: routineN='xc_rho_set_update',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, ispin, idir, nspins
    TYPE(pw_p_type), DIMENSION(2) :: my_rho_r
    TYPE(pw_p_type), DIMENSION(3,2) :: drho_r
    TYPE(pw_type), POINTER :: tmp_g, my_rho_g

    failure=.false.
    DO ispin=1,2
       nullify(my_rho_r(ispin)%pw
       DO idir=1,3
          nullify(drho_r(idir,ispin)%pw)
       END DO
    END DO
    NULLIFY(tmp_g,my_rho_g)

    CPPrecondition(ASSOCIATED(rho_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_set%ref_cont>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)

    nspins=dft_control%nspins
    rho_set%local_bounds=rho_r(1)%pw%pw_grid%bounds_local

! some checks
    SELECT CASE(dft_control%nspins)
    CASE(1)
       CPPrecondition(SIZE(rho_r)==1,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho_r(1)%pw),cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_r(1)%pw%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%rho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%norm_drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drhoa_drhob,cp_failure_level,routineP,error,failure)
    CASE(2)
       CPPrecondition(SIZE(rho_r)==2,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho_r(1)%pw),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(rho_r(2)%pw),cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_r(1)%pw%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
       CPPrecondition(rho_r(2)%pw%in_use==REALDATA3D,cp_failure_level,routineP,error,failure)
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

! preclean
    IF (rho_set%owns%rho) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rho, &
            accept_non_compatible=.TRUE., error=error)
    ELSE
       NULLIFY(rho_set%rho)
    END IF
    IF (rho_set%owns%drho) THEN
       DO idir=1,3
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%drho(idir)%array,&
               accept_non_compatible=.TRUE.,error=error)
       END DO
    ELSE
       NULLIFY(rho_set%drho)
    END IF
    IF (rho_set%owns%norm_drho) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%norm_drho,&
            accept_non_compatible=.TRUE.,error=error)
    ELSE
       NULLIFY(rho_set%norm_rho)
    END IF
    IF (rho_set%owns%rho_spin) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhoa,&
            accept_non_compatible=.TRUE., error=error)
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%rhob,&
            accept_non_compatible=.TRUE., error=error)
    ELSE
       NULLIFY(rho_set%rhoa,rho_set%rhob)
    END IF
    IF (rho_set%owns%drho_spin) THEN
       DO idir=1,3
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%drhoa(idir)%array,&
               accept_non_compatible=.TRUE.,error=error)
          CALL pw_pool_give_back_cr3d(pw_pool,rho_set%drhob(idir)%array,&
               accept_non_compatible=.TRUE.,error=error)
       END DO
    ELSE
       NULLIFY(rho_set%drhoa,rho_set%drhob)
    END IF
    IF (rho_set%owns%norm_drho_spin) THEN
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%norm_drhoa,&
            accept_non_compatible=.TRUE.,error=error)
       CALL pw_pool_give_back_cr3d(pw_pool,rho_set%norm_drhob,&
            accept_non_compatible=.TRUE.,error=error)
    ELSE
       NULLIFY(rho_set%norm_drhoa, rho_set%norm_drhob)
    END IF
    IF (rho_set%owns%drhoa_drhob) THEN
       CALL pw_pool_give_back_cr3d(pw_pool, rho_set%drhoa_drhob,&
            accept_non_compatible=.TRUE., error=error)
    ELSE
       NULLIFY(rho_set%drhoa_drhob)
    END IF

    CALL xc_rho_cflags_setall(rho_set%has,.FALSE.,error=error)

    gradient_f=(needs%drho_spin.OR.needs%norm_drho_spin.OR.&
         needs%drhoa_drhob.OR.needs%drho.OR.needs%norm_drho)

    DO ispin=1,nspins
! introduce a smoothing kernel on the density
       IF (dft_control%xc_rho_smooth_id==xc_rho_no_smooth) THEN
          my_rho_r_local=.FALSE.
          my_rho_g_local=.FALSE.
          my_rho_r(ispin)%pw => rho_r(ispin)%pw
          IF (gradient_f) THEN
             my_rho_g => rho_g(ispin)%pw
          END IF
       ELSE
          my_rho_r_local=.TRUE.
          my_rho_g_local=.TRUE.
          CALL pw_pool_create_pw(pw_pool, my_rho_r,&
               use_data=REALDATA3D, in_space=REALSPACE, &
               error=error)
          CALL pw_pool_create_pw(pw_pool, my_rho_g,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
               error=error)
          CALL pw_copy(rho_g(ispin)%pw,my_rho_g)
          SELECT CASE(dft_control%xc_rho_smooth_id)
          CASE (xc_rho_spline2_smooth)
             CALL pw_spline2_evaluate_values_g(my_rho_g,error=error)
          CASE (xc_rho_spline3_smooth)
             CALL pw_spline3_evaluate_values_g(my_rho_g,error=error)
          CASE (xc_rho_nn10)
             CALL pw_nn_smear_g(my_rho_g,&
                  nn10_coeffs, error=error)
          END SELECT
          CALL pw_transfer(my_rho_g,my_rho_r)
       END IF

       IF (gradient_f) THEN ! calculate the grad of rho
! normally when you need the gradient you need the whole gradient
! (for the partial integration)
! deriv rho
          DO idir=1,3
             NULLIFY(drho_r(idir,ispin)%pw)
             CALL pw_pool_create_pw(pw_pool,drho_r(idir,ispin)%pw, &
                  use_data=REALDATA3D, in_space=REALSPACE, &
                  error=error)
          END DO
          CALL pw_pool_create_pw(pw_pool, tmp_g,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
               error=error)
          SELECT CASE(dft_control%xc_deriv_method_id)
          CASE (xc_deriv_pw)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_derive ( tmp_g, nd(:,idir) )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline2)
             IF (.NOT.my_rho_g_local) THEN
                CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                     use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                     error=error)
                my_rho_g_local=.TRUE.
                CALL pw_copy(rho_g(ispin)%pw, my_rho_g)
             END IF
             CALL pw_spline2_interpolate_values_g(my_rho_g,error=error)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline2_smooth)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline2_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline3)
             IF (.NOT.my_rho_g_local) THEN
                CALL pw_pool_create_pw(pw_pool, my_rho_g,&
                     use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, &
                     error=error)
                CALL pw_copy(rho_g(ispin)%pw, my_rho_g)
                my_rho_g_local=.TRUE.
             END IF
             CALL pw_spline3_interpolate_values_g(my_rho_g,error=error)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE (xc_deriv_spline3_smooth)
             DO idir=1,3
                CALL pw_copy ( my_rho_g, tmp_g )
                CALL pw_spline3_deriv_g ( tmp_g, idir=idir, error=error )
                CALL pw_transfer ( tmp_g, drho_r(idir,ispin)%pw )
             END DO
          CASE default
             CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

          IF (dft_control%xc_deriv_method_id/=xc_deriv_pw) THEN
             DO idir=1,3
                drho_r_att(idir)%pw => drho_r(idir,ispin)%pw
             END DO
             CALL pw_spline_scale_deriv(drho_r_att, cell=cell,&
                  error=error)
          END IF

! my_rho_g might be either rho_g or the weights of the spline
! my_rho_r contains the function that has been derived
! drgrid contains the function that should have been derived:
!        the (maybe smoothed) rho

          IF (my_rho_g_local) THEN
             my_rho_g_local=.FALSE.
             CALL pw_pool_give_back_pw(pw_pool, my_rho_g ,error=error)
          END IF
          NULLIFY(my_rho_g)
          CALL pw_pool_give_back_pw(pw_pool, tmp_g ,error=error)

       END IF
    END DO

    SELECT CASE(dft_control%nspins)
    CASE(1)
       IF (needs%rho) THEN
          rho_set%rho => ma_rho_r(1)%pw%cr3d
          if (my_rho_r_local) nullify(my_rho_r(1)%pw%cr3d)
          rho_set%owns%rho=my_rho_r_local
          rho_set%has%rho=.true.
       END IF
       IF (needs%norm_drho) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drho, &
               error=error)
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%norm_drho(i,j,k)=sqrt(&
                        drho_r(1,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(2,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(3,1)%pw%cr3d(i,j,k)**2)
                END DO
             END DO
          END DO
          rho_set%owns%norm_drho=.true.
          rho_set%has%norm_drho=.true.
       END IF
       IF (needs%drho) THEN
          DO idir=1,3
             rho_set%drho(idir)%array => drho_r(1)%pw%cr3d
             nullify(drho_r(1)%pw%cr3d)
          END DO
          rho_set%owns%norm_drho=.true.
          rho_set%has%drho=.true.
       END IF
    CASE(2)
       IF (needs%rho) THEN
! this should basically never be the case unless you use LDA functionals
! with LSD

          CALL pw_pool_create_cr3d(pw_pool,rho_set%rho,error=error)
!assume that the bounds are the same?
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%rho(i,j,k)=rho_r(1)%pw%cr3d(i,j,k)+&
                        rho_r(2)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
          rho_set%owns%rho=.true.
          rho_set%has%rho=.true.
       END IF
       IF (needs%rho_spin) THEN

          rho_set%rhoa => my_rho_r(1)%pw%cr3d
          IF (my_rho_r_local) NULLIFY(my_rho_r(1)%pw%cr3d)

          rho_set%rhob => my_rho_r(2)%pw%cr3d
          if (my_rho_r_local) nullify(my_rho_r(2)%pw%cr3d)

          rho_set%owns%rho_spin=my_rho_r_local
          rho_set%has%rho_spin=.true.
       END IF
       IF (needs%norm_drho_spin) THEN

          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drhoa, &
               error=error)
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%norm_drhoa(i,j,k)=sqrt(&
                        drho_r(1,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(2,1)%pw%cr3d(i,j,k)**2+&
                        drho_r(3,1)%pw%cr3d(i,j,k)**2)
                END DO
             END DO
          END DO

          CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drhob, &
               error=error)
          rho_set%owns%norm_drho_spin=.true.
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%norm_drhob(i,j,k)=sqrt(&
                        drho_r(1,2)%pw%cr3d(i,j,k)**2+&
                        drho_r(2,2)%pw%cr3d(i,j,k)**2+&
                        drho_r(3,2)%pw%cr3d(i,j,k)**2)
                END DO
             END DO
          END DO

          rho_set%owns%norm_drho_spin=.true.
          rho_set%has%norm_drho_spin=.true.
       END IF
       IF (needs%drhoa_drhob) THEN
          CALL pw_pool_create_cr3d(pw_pool, rho_set%drhoa_drhob, &
               error=error)
!$omp parallel do private(i,j,k)
          DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
             DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                   rho_set%drhoa_drhob(i,j,k)=&
                        drho_r(1,1)%pw%cr3d(i,j,k)*drho_r(1,2)%pw%cr3d(i,j,k)+&
                        drho_r(2,1)%pw%cr3d(i,j,k)*drho_r(2,2)%pw%cr3d(i,j,k)+&
                        drho_r(3,1)%pw%cr3d(i,j,k)*drho_r(3,2)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
          rho_set%owns%drhoa_drhob=.true.
          rho_set%has%drhoa_drhob=.true.
       END IF
       IF (needs%drho) THEN
! this should basically never be the case unless you use LDA functionals
! with LSD
          DO idir=1,3
             CALL pw_pool_create_cr3d(pw_pool,rho_set%drho(idir)%array,&
                  error=error)
!assume that the bounds are the same?
!$omp parallel do private(i,j,k)
             DO k=rho_set%local_bounds(1,1),rho_set%local_bounds(2,1)
                DO j=rho_set%local_bounds(1,2),rho_set%local_bounds(2,2)
                   DO i=rho_set%local_bounds(1,3),rho_set%local_bounds(2,3)
                      rho_set%drho(idir)%array(i,j,k)=&
                           drho_r(idir,1)%pw%cr3d(i,j,k)+&
                           drho_r(idir,2)%pw%cr3d(i,j,k)
                   END DO
                END DO
             END DO
          END DO
          rho_set%owns%drho=.true.
          rho_set%has%drho=.true.
       END IF
       IF (needs%drho_spin) THEN
          DO idir=1,3
             rho_set%drhoa(idir)%array => drho_r(idir,1)%pw%cr3d
             NULLIFY(drho_r(idir,1)%pw%cr3d)
             rho_set%drhob(idir)%array => drho_r(idir,2)%pw%cr3d
             NULLIFY(drho_r(idir,2)%pw%cr3d)
          END DO
          rho_set%owns%drho_spin=.true.
          rho_set%has%drho_spin=.true.
       END IF
    END SELECT

! post cleanup
    DO ispin=1,nspins
       DO idir=1,3
          CALL pw_pool_give_back_pw(pw_pool, drho_r(idir,ispin), &
               accept_non_compatible=.TRUE., error=error)
       END DO
    END DO
    IF (my_rho_r_local) THEN
       DO ispin=1,nspins
          CALL pw_pool_give_back_pw(pw_pool, my_rho_r(ispin)%pw,&
               accept_non_compatible=.TRUE., error=error)
       END DO
    END IF

  END SUBROUTINE xc_rho_set_update
!***************************************************************************

END MODULE xc_rho_set_types
