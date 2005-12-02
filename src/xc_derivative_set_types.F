!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_derivative_set_types [1.0] *
!!
!!   NAME
!!     xc_derivative_set_types
!!
!!   FUNCTION
!!     represent a group ofunctional derivatives
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi & thomas
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE xc_derivative_set_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_assertion_failed
  USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_dealloc,&
                                             cp_sll_xc_deriv_insert_ordered,&
                                             cp_sll_xc_deriv_next,&
                                             cp_sll_xc_deriv_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_create_cr3d,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type,&
                                             pw_pool_create
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids, only: pw_grid_create,&
       pw_grid_release
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             standardize_derivative_desc
  USE xc_derivative_types,             ONLY: xc_derivative_create,&
                                             xc_derivative_release,&
                                             xc_derivative_type
  USE message_passing, only: MPI_COMM_SELF

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='xc_derivative_set_types'

  INTEGER, SAVE :: derivative_set_last_id_nr  = 0

  PUBLIC :: xc_derivative_set_type
  PUBLIC :: xc_dset_create, xc_dset_release, xc_dset_retain, &
            xc_dset_get_derivative, xc_dset_zero_all
!!***
!****************************************************************************

!!****t* xc_derivative_set_types/xc_derivative_set_type  *
!!
!!   NAME
!!     xc_derivative_set_type
!!
!!   FUNCTION
!!     A derivative set contains the different derivatives of a xc-functional
!!     in form of a linked list
!!
!!****
  TYPE xc_derivative_set_type
     INTEGER :: ref_count, id_nr
     TYPE(pw_pool_type), POINTER :: pw_pool
     TYPE(cp_sll_xc_deriv_type), POINTER :: derivs
  END TYPE xc_derivative_set_type

CONTAINS

!!****f* xc_derivative_set_types/xc_dset_get_derivative [1.0] *
!!
!!   NAME
!!     xc_dset_get_derivative
!!
!!   SYNOPSIS
!!     Function xc_dset_get_derivative(derivative_set, description,&
!!         allocate_deriv, error) Result(res)
!!       Type(xc_derivative_set_type), Pointer:: derivative_set
!!       Character(Len=*), Intent (IN):: description
!!       Logical, Intent (IN), Optional:: allocate_deriv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Type(xc_derivative_type), Pointer:: res
!!     End Function xc_dset_get_derivative
!!
!!   FUNCTION
!!     returns the requested xc_derivative
!!
!!   NOTES
!!     -
!!     
!!   INPUTS
!!     - derivative_set: the set where to search for the derivative
!!     - description   : the description of the derivative you want to have
!!     - allocate_deriv: if the derivative should be allocated when not present
!!                       Defaults to false.
!!     - error         : variable to control error logging, stopping,... 
!!                       see module cp_error_handling 
!!
!!*** **********************************************************************
  FUNCTION xc_dset_get_derivative(derivative_set, description, allocate_deriv, error) &
       RESULT(res)

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    CHARACTER(len=*), INTENT(in)             :: description
    LOGICAL, INTENT(in), OPTIONAL            :: allocate_deriv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(xc_derivative_type), POINTER        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_get_derivative', &
      routineP = moduleN//':'//routineN

    CHARACTER&
      (len=MAX_DERIVATIVE_DESC_LENGTH)       :: std_deriv_desc
    LOGICAL                                  :: failure, my_allocate_deriv
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d_ptr
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(xc_derivative_type), POINTER        :: deriv_att

    failure=.FALSE.
    NULLIFY(pos,deriv_att,cr3d_ptr)

    CPPrecondition(ASSOCIATED(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative_set%ref_count>0,cp_failure_level,routineP,error,failure)

    my_allocate_deriv=.FALSE.
    IF (PRESENT(allocate_deriv)) my_allocate_deriv=allocate_deriv
    NULLIFY(res)
    CALL standardize_derivative_desc(description,std_deriv_desc)
    IF (.NOT. failure) THEN
       pos => derivative_set%derivs
       DO WHILE(cp_sll_xc_deriv_next(pos,el_att=deriv_att,error=error))
          IF (deriv_att%desc == std_deriv_desc) THEN
             res => deriv_att
             EXIT
          END IF
       END DO
       IF (.NOT.ASSOCIATED(res).AND.my_allocate_deriv) THEN
          CALL pw_pool_create_cr3d(derivative_set%pw_pool,cr3d_ptr,error=error)
          cr3d_ptr=0.0_dp
          CALL xc_derivative_create(res, std_deriv_desc, & 
                                    cr3d_ptr=cr3d_ptr, error=error)
          CALL cp_sll_xc_deriv_insert_ordered(derivative_set%derivs,res,&
               error=error)
       END IF
    END IF
  END FUNCTION xc_dset_get_derivative
  !***************************************************************************

!!****f* xc_derivative_set_types/xc_dset_create [1.0] *
!!
!!   NAME
!!     xc_dset_create
!!
!!   SYNOPSIS
!!     Subroutine xc_dset_create(derivative_set, pw_pool, error)
!!       Type(xc_derivative_set_type), Pointer:: derivative_set
!!       Type(pw_pool_type), Pointer:: pw_pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_dset_create
!!
!!   FUNCTION
!!     creates a derivative set object
!!
!!   INPUTS
!!     - derivative_set : the set where to search for the derivative
!!     - pw_pool: pool where to get the cr3d arrays needed to store the
!!       derivatives
!!     - error   : variable to control error logging, stopping,... 
!!                 see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_dset_create(derivative_set, pw_pool, local_bounds, error)

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    TYPE(pw_pool_type), POINTER, optional    :: pw_pool
    INTEGER, DIMENSION(2,3), OPTIONAL        :: local_bounds
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(pw_grid_type), pointer              :: pw_grid

    failure=.FALSE.
    nullify(pw_grid)
    CPPrecondition(.not.ASSOCIATED(derivative_set),cp_failure_level,routineP,error,failure)

    ALLOCATE(derivative_set, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       NULLIFY(derivative_set%derivs)
       derivative_set%ref_count  = 1
       derivative_set_last_id_nr      = derivative_set_last_id_nr + 1
       derivative_set%id_nr      = derivative_set_last_id_nr
       IF (PRESENT(pw_pool)) THEN
          derivative_set%pw_pool => pw_pool
          CALL pw_pool_retain(pw_pool, error=error)
          IF (PRESENT(local_bounds)) THEN
             CALL cp_assert(ALL(pw_pool%pw_grid%bounds_local==local_bounds),&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "incompatible local_bounds and pw_pool",error,failure)
          END IF
       ELSE
          !FM ugly hack, should be replaced by a pool only for 3d arrays
          CPPrecondition(PRESENT(local_bounds),cp_failure_level,routineP,error,failure)
          CALL pw_grid_create(pw_grid,MPI_COMM_SELF)
          pw_grid%bounds_local=local_bounds
          NULLIFY(derivative_set%pw_pool)
          CALL pw_pool_create(derivative_set%pw_pool, pw_grid, error=error)
          CALL pw_grid_release(pw_grid,error=error)
       END IF
    END IF

  END SUBROUTINE xc_dset_create
  !***************************************************************************

!!****f* xc_derivative_set_types/xc_dset_release [1.0] *
!!
!!   NAME
!!     xc_dset_release
!!
!!   SYNOPSIS
!!     Subroutine xc_dset_release(derivative_set, error)
!!       Type(xc_derivative_set_type), Pointer:: derivative_set
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_dset_release
!!
!!   FUNCTION
!!     releases a derivative set
!!
!!   INPUTS
!!     - set     : the set to release
!!     - error   : variable to control error logging, stopping,... 
!!                 see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE xc_dset_release(derivative_set, error)

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(xc_derivative_type), POINTER        :: deriv_att

    failure=.FALSE.
    
    NULLIFY(deriv_att,pos)
    CPPrecondition(ASSOCIATED(derivative_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(derivative_set%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       derivative_set%ref_count = derivative_set%ref_count - 1
       IF (derivative_set%ref_count == 0) THEN
          pos => derivative_set%derivs
          DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att,error=error))
             CALL xc_derivative_release(deriv_att, pw_pool=derivative_set%pw_pool,&
                  error=error)
          END DO
          CALL cp_sll_xc_deriv_dealloc(derivative_set%derivs,error=error)
          CALL pw_pool_release(derivative_set%pw_pool,error=error)

          DEALLOCATE(derivative_set, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(derivative_set)

  END SUBROUTINE xc_dset_release
  !***************************************************************************

!!****f* xc_derivative_set_types/xc_dset_retain [1.0] *
!!
!!   NAME
!!     xc_dset_retain
!!
!!   SYNOPSIS
!!     Subroutine xc_dset_retain(deriv_set, error)
!!       Type(xc_derivative_set_type), Pointer:: deriv_set
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine xc_dset_retain
!!
!!   FUNCTION
!!     retains the given derivative set
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - deriv_set: the derivative set to retain
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
SUBROUTINE xc_dset_retain(deriv_set,error)
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(deriv_set%ref_count>0,cp_failure_level,routineP,error)
     deriv_set%ref_count=deriv_set%ref_count+1
  END IF
END SUBROUTINE xc_dset_retain
!***************************************************************************

SUBROUTINE xc_dset_zero_all(deriv_set, error)
  
  TYPE(xc_derivative_set_type), POINTER    :: deriv_set
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                            :: error
  
  LOGICAL                                  :: failure
  TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
  TYPE(xc_derivative_type), POINTER        :: deriv_att

  CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_zero_all', &
                                 routineP = moduleN//':'//routineN
  
  failure = .FALSE.
  NULLIFY(pos, deriv_att)
  
  CPPrecondition(ASSOCIATED(deriv_set),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     pos => deriv_set%derivs
     DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att, error=error)) 
        deriv_att%deriv_data = 0.0_dp
     END DO
  END IF

END SUBROUTINE xc_dset_zero_all

!***************************************************************************

  
END MODULE xc_derivative_set_types
