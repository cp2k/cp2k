!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/neb_md_utils [1.0] *
!!
!!   NAME
!!     neb_md_utils
!!
!!   FUNCTION
!!     Module with utility to perform MD Nudged Elastic Band Calculation
!!
!!   NOTES
!!     Numerical accuracy for parallel runs:
!!      Each replica starts the SCF run from the one optimized
!!      in a previous run. It may happen then energies and derivatives
!!      of a serial run and a parallel run could be slightly different
!!      'cause of a different starting density matrix.
!!      Exact results are obtained using:
!!         EXTRAPOLATION USE_GUESS in QS section (Teo 09.2006)
!!
!!   AUTHOR
!!     Teodoro Laino 10.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
MODULE neb_md_utils
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_md_utils'

  PUBLIC :: neb_initialize_velocity,&
            control_vels_a,&
            control_vels_b,&
            get_temperatures

CONTAINS

!!****f* neb_md_utils/neb_initialize_velocity [1.0] *
!!
!!   NAME
!!     neb_initialize_velocity
!!
!!   FUNCTION
!!     Initialize velocities of replica in an MD optimized algorithm within NEB
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE neb_initialize_velocity(vels, neb_section, particle_set, i, iw, globenv,&
       para_env, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vels
    TYPE(section_vals_type), POINTER         :: neb_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: i, iw
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_initialize_velocity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, kind
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: akin, mass, mass_tot, sc, &
                                                temp, temp_ext, vcom(3)
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(unit_convert_type), POINTER         :: unit

    failure = .FALSE.
    md_section => section_vals_get_subs_vals(neb_section,"OPTIMIZE_BAND%MD",error=error)
    CALL section_vals_get(md_section,explicit=explicit,error=error)
    IF ((.NOT.failure).AND.(explicit)) THEN
       CALL section_vals_val_get(md_section,"TEMPERATURE",r_val=temp_ext,error=error)
       NULLIFY(unit)
       CALL get_cp2k_units ( unit )
       ! Initialize velocity according external temperature
       mass_tot = 0.0_dp
       vcom     = 0.0_dp
       vels(:,i)= 0.0_dp
       DO k = 1, SIZE(particle_set)
          kind = (k-1)*3
          vels(kind+1,i) = next_random_number(globenv%gaussian_rng_stream,error=error)
          vels(kind+2,i) = next_random_number(globenv%gaussian_rng_stream,error=error)
          vels(kind+3,i) = next_random_number(globenv%gaussian_rng_stream,error=error)
       END DO
       DO k = 1, SIZE(particle_set)
          kind = (k-1)*3
          mass = particle_set(k)%atomic_kind%mass
          mass_tot = mass_tot + mass
          vcom(1) = vcom(1) + mass * vels(kind+1,i)
          vcom(2) = vcom(2) + mass * vels(kind+2,i)
          vcom(3) = vcom(3) + mass * vels(kind+3,i)
       END DO
       vcom = vcom / mass_tot
       akin = 0.0_dp
       DO k = 1, SIZE(particle_set)
          kind = (k-1)*3
          mass = particle_set(k)%atomic_kind%mass
          vels(kind+1:kind+3,i) = vels(kind+1:kind+3,i) - vcom(:)
          akin = akin + 0.5_dp*(mass*DOT_PRODUCT(vels(kind+1:kind+3,i),vels(kind+1:kind+3,i)))
       END DO
       temp = 2.0_dp * akin / REAL ( 3*SIZE(particle_set),KIND=dp)
       ! Scale velocities to get the correct initial temperature and
       sc        = SQRT ( temp_ext / temp )
       vels(:,i) = vels(:,i) * sc
       ! Moving the velocities so that the COM velocity is zero
       akin = 0.0_dp
       vcom = 0.0_dp
       DO k = 1, SIZE(particle_set)
          kind = (k-1)*3
          mass = particle_set(k)%atomic_kind%mass
          vcom(1) = vcom(1) + mass * vels(kind+1,i)
          vcom(2) = vcom(2) + mass * vels(kind+2,i)
          vcom(3) = vcom(3) + mass * vels(kind+3,i)
          akin = akin + 0.5_dp*(mass*DOT_PRODUCT(vels(kind+1:kind+3,i),vels(kind+1:kind+3,i)))
       END DO
       vcom = vcom / mass_tot
       IF (iw>0) THEN
          temp = 2.0_dp * akin / REAL ( 3*SIZE(particle_set),KIND=dp)
          WRITE ( iw, '( A, T61, F18.2, A2 )' ) &
               ' NEB| Initial Temperature ', temp / unit % kelvin, " K"
          WRITE ( iw, '( A, T61, F20.12 )' ) &
               ' NEB| Centre of mass velocity in direction x:', vcom ( 1 )
          WRITE ( iw, '( A, T61, F20.12 )' ) &
               ' NEB| Centre of mass velocity in direction y:', vcom ( 2 )
          WRITE ( iw, '( A, T61, F20.12 )' ) &
               ' NEB| Centre of mass velocity in direction z:', vcom ( 3 )
          WRITE(iw,'(T2,"NEB|",75("*"))')
       END IF
    ELSE IF (.NOT.explicit) THEN
       vels(:,i) = 0.0_dp
    END IF
  END SUBROUTINE neb_initialize_velocity

!!****f* neb_md_utils/control_vels_a [1.0] *
!!
!!   NAME
!!     control_vels_a
!!
!!   FUNCTION
!!     Control on  velocities - I part
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE control_vels_a(vels, particle_set, tc_section, vc_section, &
       output_unit, istep, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: tc_section, vc_section
    INTEGER, INTENT(IN)                      :: output_unit, istep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'control_vels_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat, temp_tol_steps
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: ext_temp, f_annealing, &
                                                factor, scale, temp_tol, &
                                                temploc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: temperatures
    TYPE(unit_convert_type), POINTER         :: unit

    failure = .FALSE.
    NULLIFY(unit)
    CALL get_cp2k_units ( unit )
    ! Temperature control
    CALL section_vals_get(tc_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(tc_section,"TEMP_TOL_STEPS",i_val=temp_tol_steps, error=error)
       CALL section_vals_val_get(tc_section,"TEMPERATURE",r_val=ext_temp, error=error)
       CALL section_vals_val_get(tc_section,"TEMP_TOL",r_val=temp_tol, error=error)
       ALLOCATE(temperatures(SIZE(vels,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Computes temperatures
       CALL get_temperatures(vels,particle_set,temperatures,factor=1.0_dp,error=error)
       ! Possibly rescale
       IF (istep<=temp_tol_steps) THEN
          DO i = 2, SIZE(vels,2)-1
             temploc = temperatures(i)
             IF (ABS(temploc-ext_temp) > temp_tol) THEN
                IF (output_unit>0) THEN
                   WRITE(output_unit,'(T2,"NEB| Replica Nr.",I5,'//&
                        '"  - Velocity rescaled from: ",F12.6," to: ",F12.6,".")')&
                        i, temploc/ unit % kelvin, ext_temp/ unit % kelvin
                END IF
                scale = SQRT(ext_temp/temploc)
                vels(:,i) = scale * vels(:,i)
             END IF
          END DO
       END IF
       DEALLOCATE(temperatures,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    ! Annealing
    CALL section_vals_get(vc_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(vc_section,"ANNEALING",r_val=f_annealing, error=error)
       DO i = 2, SIZE(vels,2)-1
          vels(:,i) = f_annealing * vels(:,i)
       END DO
    END IF
  END SUBROUTINE control_vels_a

!!****f* neb_md_utils/control_vels_b [1.0] *
!!
!!   NAME
!!     control_vels_b
!!
!!   FUNCTION
!!     Control on velocities - II part
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE control_vels_b(vels, forces, vc_section, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vels, forces
    TYPE(section_vals_type), POINTER         :: vc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'control_vels_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: explicit, failure, lval
    REAL(KIND=dp)                            :: factor, norm

    failure = .FALSE.
    ! Check the sign of V.dot.F
    CALL section_vals_get(vc_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(vc_section,"PROJ_VELOCITY_VERLET",l_val=lval, error=error)
       IF (lval) THEN
          DO i = 2, SIZE(vels,2)-1
             norm   = DOT_PRODUCT(forces(:,i),forces(:,i))
             factor = DOT_PRODUCT(vels(:,i),forces(:,i))
             IF (factor>0.AND.(norm >= EPSILON(0.0_dp))) THEN
                vels(:,i) = factor/norm*forces(:,i)
             ELSE
                vels(:,i) = 0.0_dp
             END IF
          END DO
       END IF
       CALL section_vals_val_get(vc_section,"SD_LIKE",l_val=lval, error=error)
       IF (lval) THEN
          DO i = 2, SIZE(vels,2)-1
             vels(:,i) = 0.0_dp
          END DO
       END IF
    END IF
  END SUBROUTINE control_vels_b

!!****f* neb_md_utils/get_temperatures [1.0] *
!!
!!   NAME
!!     get_temperatures
!!
!!   FUNCTION
!!     Computes temperatures
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE get_temperatures(vels, particle_set, temperatures, ekin, factor, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: temperatures
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(OUT), OPTIONAL                  :: ekin
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_temperatures', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, kind
    REAL(KIND=dp)                            :: akin, mass, myfactor, temploc
    TYPE(unit_convert_type), POINTER         :: unit

    NULLIFY(unit)
    CALL get_cp2k_units ( unit )
    myfactor     = unit % kelvin
    IF (PRESENT(factor)) myfactor = factor
    IF (PRESENT(ekin))   ekin = 0.0_dp
    temperatures = 0.0_dp
    DO i = 2, SIZE(vels,2)-1
       akin = 0.0_dp
       DO k = 1, SIZE(particle_set)
          kind = (k-1)*3
          mass = particle_set(k)%atomic_kind%mass
          akin = akin + 0.5_dp*(mass*DOT_PRODUCT(vels(kind+1:kind+3,i),vels(kind+1:kind+3,i)))
       END DO
       IF (PRESENT(ekin)) ekin(i) = akin
       temploc = 2.0_dp * akin / REAL ( 3*SIZE(particle_set),KIND=dp)
       temperatures(i) = temploc / myfactor
    END DO
  END SUBROUTINE get_temperatures

END MODULE neb_md_utils
