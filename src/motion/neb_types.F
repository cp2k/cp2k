!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Typo for Nudged Elastic Band Calculation
!> \note
!>      Numerical accuracy for parallel runs:
!>       Each replica starts the SCF run from the one optimized
!>       in a previous run. It may happen then energies and derivatives
!>       of a serial run and a parallel run could be slightly different
!>       'cause of a different starting density matrix.
!>       Exact results are obtained using:
!>          EXTRAPOLATION USE_GUESS in QS section (Teo 09.2006)
!> \author Teodoro Laino 10.2006
! *****************************************************************************
MODULE neb_types
  
  USE force_env_types,                 ONLY: force_env_type
  USE input_constants,                 ONLY: do_band_cartesian,&
                                             do_band_collective
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_types'

  PUBLIC :: neb_type,&
            neb_var_type,&
            neb_var_create,&
            neb_var_release

! *****************************************************************************
  TYPE neb_type
     ! NEB parameters
     INTEGER       :: id_type
     INTEGER       :: opt_type
     INTEGER       :: pot_type
     INTEGER       :: number_of_replica, nsize_xyz, nsize_int
     INTEGER       :: nsteps_it, istep
     INTEGER       :: nr_HE_image
     LOGICAL       :: rotate_frames, align_frames
     LOGICAL       :: optimize_end_points
     LOGICAL       :: use_colvar
     LOGICAL       :: reparametrize_frames
     INTEGER       :: spline_order
     REAL(KIND=dp) :: K, spring_energy, avg_distance
     REAL(KIND=dp) :: smoothing
     CHARACTER(LEN=default_string_length) :: opt_type_label
     ! Section used for restart
     TYPE(force_env_type), POINTER        :: force_env
     TYPE(section_vals_type), POINTER     :: root_section
     TYPE(section_vals_type), POINTER     :: motion_print_section
     TYPE(section_vals_type), POINTER     :: force_env_section
     TYPE(section_vals_type), POINTER     :: neb_section
  END TYPE neb_type

! *****************************************************************************
  TYPE neb_var_type
     INTEGER :: in_use, size_wrk(2)
     REAL(KIND=dp), DIMENSION(:, :),  POINTER  :: xyz, int, wrk
     REAL(KIND=dp), DIMENSION(:,:,:), POINTER  :: Mmatrix
  END TYPE neb_var_type

CONTAINS
! *****************************************************************************
!> \brief Creates a variable type for BAND calculation
!> \param neb_var ...
!> \param neb_env ...
!> \param full_allocation ...
!> \date   05.2007
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE neb_var_create(neb_var, neb_env, full_allocation)
    TYPE(neb_var_type), POINTER              :: neb_var
    TYPE(neb_type), POINTER                  :: neb_env
    LOGICAL, INTENT(IN), OPTIONAL            :: full_allocation

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_var_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: neb_nr_replica
    LOGICAL                                  :: allocate_all, failure

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(neb_var))
    allocate_all = .FALSE.
    IF (PRESENT(full_allocation)) allocate_all = full_allocation
    neb_nr_replica = neb_env%number_of_replica
    ALLOCATE(neb_var)
    NULLIFY(neb_var%xyz, neb_var%int, neb_var%wrk, neb_var%Mmatrix)
    IF (allocate_all) THEN
       ALLOCATE(neb_var%xyz(neb_env%nsize_xyz, neb_nr_replica))
       neb_var%xyz = 0.0_dp
    END IF
    IF (neb_env%use_colvar) THEN
       neb_var%in_use = do_band_collective
       CPASSERT(neb_env%nsize_int>0)
       ALLOCATE(neb_var%int(neb_env%nsize_int, neb_nr_replica))
       neb_var%int = 0.0_dp
       neb_var%wrk => neb_var%int
    ELSE
       neb_var%in_use = do_band_cartesian
       IF (.NOT.allocate_all) THEN
          ALLOCATE(neb_var%xyz(neb_env%nsize_xyz, neb_nr_replica))
          neb_var%xyz = 0.0_dp
       END IF
       neb_var%wrk => neb_var%xyz
    END IF
    neb_var%size_wrk(1) = SIZE(neb_var%wrk,1)
    neb_var%size_wrk(2) = SIZE(neb_var%wrk,2)

  END SUBROUTINE neb_var_create

! *****************************************************************************
!> \brief Releases a variable type for BAND calculation
!> \param neb_var ...
!> \date   05.2007
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE neb_var_release(neb_var)
    TYPE(neb_var_type), POINTER              :: neb_var

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_var_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(neb_var))
    IF (ASSOCIATED(neb_var%xyz)) THEN
       DEALLOCATE(neb_var%xyz)
    END IF
    IF (neb_var%in_use==do_band_collective) THEN
       DEALLOCATE(neb_var%int)
    END IF
    NULLIFY(neb_var%wrk)
    DEALLOCATE(neb_var)

  END SUBROUTINE neb_var_release

END MODULE neb_types
