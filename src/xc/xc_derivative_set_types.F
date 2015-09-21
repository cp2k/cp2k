!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represent a group ofunctional derivatives
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi & thomas
! *****************************************************************************
MODULE xc_derivative_set_types
  USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_dealloc,&
                                             cp_sll_xc_deriv_insert_ordered,&
                                             cp_sll_xc_deriv_next,&
                                             cp_sll_xc_deriv_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: MPI_COMM_SELF
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_release
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_create_cr3d,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH,&
                                             standardize_derivative_desc
  USE xc_derivative_types,             ONLY: xc_derivative_create,&
                                             xc_derivative_release,&
                                             xc_derivative_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_derivative_set_types'

  INTEGER, SAVE :: derivative_set_last_id_nr  = 0

  PUBLIC :: xc_derivative_set_type
  PUBLIC :: xc_dset_create, xc_dset_release,&
            xc_dset_get_derivative, xc_dset_zero_all

! *****************************************************************************
!> \brief A derivative set contains the different derivatives of a xc-functional
!>      in form of a linked list
! *****************************************************************************
  TYPE xc_derivative_set_type
     INTEGER :: ref_count, id_nr
     TYPE(pw_pool_type), POINTER :: pw_pool
     TYPE(cp_sll_xc_deriv_type), POINTER :: derivs
  END TYPE xc_derivative_set_type

CONTAINS

! *****************************************************************************
!> \brief returns the requested xc_derivative
!> \param derivative_set the set where to search for the derivative
!> \param description the description of the derivative you want to have
!> \param allocate_deriv if the derivative should be allocated when not present
!>                        Defaults to false.
!> \retval res ...
! *****************************************************************************
  FUNCTION xc_dset_get_derivative(derivative_set, description, allocate_deriv) &
       RESULT(res)

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    CHARACTER(len=*), INTENT(in)             :: description
    LOGICAL, INTENT(in), OPTIONAL            :: allocate_deriv
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

    CPASSERT(ASSOCIATED(derivative_set))
    CPASSERT(derivative_set%ref_count>0)

    my_allocate_deriv=.FALSE.
    IF (PRESENT(allocate_deriv)) my_allocate_deriv=allocate_deriv
    NULLIFY(res)
    CALL standardize_derivative_desc(description,std_deriv_desc)
    pos => derivative_set%derivs
    DO WHILE(cp_sll_xc_deriv_next(pos,el_att=deriv_att))
       IF (deriv_att%desc == std_deriv_desc) THEN
          res => deriv_att
          EXIT
       END IF
    END DO
    IF (.NOT.ASSOCIATED(res).AND.my_allocate_deriv) THEN
       CALL pw_pool_create_cr3d(derivative_set%pw_pool,cr3d_ptr)
       cr3d_ptr=0.0_dp
       CALL xc_derivative_create(res, std_deriv_desc, &
                                 cr3d_ptr=cr3d_ptr)
       CALL cp_sll_xc_deriv_insert_ordered(derivative_set%derivs,res)
    END IF
  END FUNCTION xc_dset_get_derivative

! *****************************************************************************
!> \brief creates a derivative set object
!> \param derivative_set the set where to search for the derivative
!> \param pw_pool pool where to get the cr3d arrays needed to store the
!>        derivatives
!> \param local_bounds ...
! *****************************************************************************
  SUBROUTINE xc_dset_create(derivative_set, pw_pool, local_bounds)

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_pool
    INTEGER, DIMENSION(2, 3), OPTIONAL       :: local_bounds

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(pw_grid_type), POINTER              :: pw_grid

    failure=.FALSE.
    NULLIFY(pw_grid)
    CPASSERT(.not.ASSOCIATED(derivative_set))

    ALLOCATE(derivative_set)

    NULLIFY(derivative_set%derivs)
    derivative_set%ref_count  = 1
    derivative_set_last_id_nr      = derivative_set_last_id_nr + 1
    derivative_set%id_nr      = derivative_set_last_id_nr
    IF (PRESENT(pw_pool)) THEN
       derivative_set%pw_pool => pw_pool
       CALL pw_pool_retain(pw_pool)
       IF (PRESENT(local_bounds)) THEN
          CALL cp_assert(ALL(pw_pool%pw_grid%bounds_local==local_bounds),&
               cp_failure_level,cp_assertion_failed,routineP,&
               "incompatible local_bounds and pw_pool",failure)
       END IF
    ELSE
       !FM ugly hack, should be replaced by a pool only for 3d arrays
       CPASSERT(PRESENT(local_bounds))
       CALL pw_grid_create(pw_grid,MPI_COMM_SELF)
       pw_grid%bounds_local=local_bounds
       NULLIFY(derivative_set%pw_pool)
       CALL pw_pool_create(derivative_set%pw_pool, pw_grid)
       CALL pw_grid_release(pw_grid)
    END IF

  END SUBROUTINE xc_dset_create

! *****************************************************************************
!> \brief releases a derivative set
!> \param derivative_set the set to release
! *****************************************************************************
  SUBROUTINE xc_dset_release(derivative_set)

    TYPE(xc_derivative_set_type), POINTER    :: derivative_set

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(xc_derivative_type), POINTER        :: deriv_att

    failure=.FALSE.

    NULLIFY(deriv_att,pos)
    CPASSERT(ASSOCIATED(derivative_set))
    CPASSERT(derivative_set%ref_count>0)

    derivative_set%ref_count = derivative_set%ref_count - 1
    IF (derivative_set%ref_count == 0) THEN
       pos => derivative_set%derivs
       DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att))
          CALL xc_derivative_release(deriv_att, pw_pool=derivative_set%pw_pool)
       END DO
       CALL cp_sll_xc_deriv_dealloc(derivative_set%derivs)
       CALL pw_pool_release(derivative_set%pw_pool)

       DEALLOCATE(derivative_set)
    END IF
    NULLIFY(derivative_set)

  END SUBROUTINE xc_dset_release

! *****************************************************************************
!> \brief retains the given derivative set
!> \param deriv_set the derivative set to retain
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE xc_dset_retain(deriv_set)
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(deriv_set))
  CPASSERT(deriv_set%ref_count>0)
  deriv_set%ref_count=deriv_set%ref_count+1
END SUBROUTINE xc_dset_retain

! *****************************************************************************
!> \brief ...
!> \param deriv_set ...
! *****************************************************************************
SUBROUTINE xc_dset_zero_all(deriv_set)

    TYPE(xc_derivative_set_type), POINTER    :: deriv_set

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_dset_zero_all', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos
    TYPE(xc_derivative_type), POINTER        :: deriv_att

  failure = .FALSE.
  NULLIFY(pos, deriv_att)

  CPASSERT(ASSOCIATED(deriv_set))
  pos => deriv_set%derivs
  DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att))
     deriv_att%deriv_data = 0.0_dp
  END DO

END SUBROUTINE xc_dset_zero_all

END MODULE xc_derivative_set_types
