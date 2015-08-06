!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains type used for a Simulation Cell Optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
MODULE cell_opt_types
  USE cell_opt_utils,                  ONLY: get_ut_cell_matrix,&
                                             read_external_press_tensor
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_types'
  INTEGER, PRIVATE, SAVE :: last_cell_opt_env_id=0

  PUBLIC :: cell_opt_env_type,&
            cell_opt_env_create,&
            cell_opt_env_release

! *****************************************************************************
!> \brief Type containing all informations abour the simulation cell optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  TYPE cell_opt_env_type
     ! Simulation cell optimization parameters
     INTEGER                                    :: ref_count, id_nr
     LOGICAL                                    :: keep_angles,&
                                                   keep_symmetry
     REAL(KIND=dp)                              :: pres_ext, pres_int, pres_tol
     REAL(KIND=dp), DIMENSION(3,3)              :: mtrx
     REAL(KIND=dp), DIMENSION(3,3)              :: rot_matrix
     TYPE(cell_type), POINTER                   :: ref_cell
  END TYPE cell_opt_env_type


CONTAINS

! *****************************************************************************
!> \brief ...
!> \param cell_env ...
!> \param force_env ...
!> \param geo_section ...
!> \param error ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  SUBROUTINE cell_opt_env_create(cell_env, force_env, geo_section, error)
    TYPE(cell_opt_env_type), POINTER         :: cell_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: geo_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_opt_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: r
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    CPPostcondition(.NOT.ASSOCIATED(cell_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(cell_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(cell_env%ref_cell,cell,subsys,particles)
    cell_env%ref_count=1
    last_cell_opt_env_id=last_cell_opt_env_id+1
    cell_env%id_nr=last_cell_opt_env_id
    CALL force_env_get(force_env, cell=cell, subsys=subsys, error=error)
    CALL cell_create(cell_env%ref_cell,error=error)
    CALL cell_clone(cell, cell_env%ref_cell, error)
    CALL section_vals_val_get(geo_section,"KEEP_ANGLES",l_val=cell_env%keep_angles,error=error)
    CALL section_vals_val_get(geo_section,"KEEP_SYMMETRY",l_val=cell_env%keep_symmetry,error=error)
    CALL section_vals_val_get(geo_section,"PRESSURE_TOLERANCE",r_val=cell_env%pres_tol,error=error)

    ! First let's rotate the cell vectors in order to have an upper triangular matrix.
    CALL get_ut_cell_matrix(cell, error)

    ! Compute the rotation matrix that give the cell vectors in the "canonical" orientation
    cell_env%rot_matrix = MATMUL(cell_env%ref_cell%hmat,cell%h_inv)

    ! Get the external pressure
    CALL read_external_press_tensor(geo_section,cell,cell_env%pres_ext,cell_env%mtrx,&
         cell_env%rot_matrix, error)

    ! Rotate particles accordingly
    CALL cp_subsys_get(subsys, particles=particles, error=error)
    DO ip=1,particles%n_els
       r = MATMUL(TRANSPOSE(cell_env%rot_matrix),particles%els(ip)%r)
       particles%els(ip)%r = r
    END DO
  END SUBROUTINE cell_opt_env_create

! *****************************************************************************
!> \brief ...
!> \param cell_env ...
!> \param error ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  SUBROUTINE cell_opt_env_retain(cell_env, error)
    TYPE(cell_opt_env_type), POINTER         :: cell_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_opt_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(cell_env),cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(cell_env%ref_count>0,cp_failure_level,routineP,error)
    cell_env%ref_count=cell_env%ref_count+1
  END SUBROUTINE cell_opt_env_retain

! *****************************************************************************
!> \brief ...
!> \param cell_env ...
!> \param error ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! *****************************************************************************
  SUBROUTINE cell_opt_env_release(cell_env, error)
    TYPE(cell_opt_env_type), POINTER         :: cell_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cell_opt_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(cell_env)) THEN
       CPPreconditionNoFail(cell_env%ref_count>0,cp_failure_level,routineP,error)
       cell_env%ref_count=cell_env%ref_count-1
       IF (cell_env%ref_count==0) THEN
          CALL cell_release(cell_env%ref_cell,error)
          DEALLOCATE(cell_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE cell_opt_env_release

END MODULE cell_opt_types
