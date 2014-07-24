!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module with utility to perform MD Nudged Elastic Band Calculation
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
MODULE neb_md_utils
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: band_md_opt,&
                                             do_band_collective
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE neb_types,                       ONLY: neb_type,&
                                             neb_var_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_types,                  ONLY: get_particle_pos_or_vel,&
                                             particle_type,&
                                             update_particle_pos_or_vel
  USE physcon,                         ONLY: kelvin,&
                                             massunit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_md_utils'

  PUBLIC :: neb_initialize_velocity,&
            control_vels_a,&
            control_vels_b,&
            get_temperatures

CONTAINS

! *****************************************************************************
!> \brief Initialize velocities of replica in an MD optimized algorithm within
!>        NEB
!> \param vels ...
!> \param neb_section ...
!> \param particle_set ...
!> \param i_rep ...
!> \param iw ...
!> \param globenv ...
!> \param neb_env ...
!> \param error ...
!> \par History
!>      25.11.2010 Consider core-shell model (MK)
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE neb_initialize_velocity(vels,neb_section,particle_set,i_rep,iw,&
                                     globenv,neb_env,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vels
    TYPE(section_vals_type), POINTER         :: neb_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: i_rep, iw
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'neb_initialize_velocity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ivar, k, natom, &
                                                nparticle, nvar
    REAL(KIND=dp)                            :: akin, mass, mass_tot, sc, &
                                                temp, temp_ext, tmp_r1
    REAL(KIND=dp), DIMENSION(3)              :: v, vcom
    TYPE(section_vals_type), POINTER         :: md_section

    IF (neb_env%opt_type == band_md_opt) THEN
       md_section => section_vals_get_subs_vals(neb_section,"OPTIMIZE_BAND%MD",error=error)
       CALL section_vals_val_get(md_section,"TEMPERATURE",r_val=temp_ext,error=error)
       ! Initialize velocity according to external temperature
       nparticle = SIZE(vels,1)
       natom = SIZE(particle_set)
       mass_tot = 0.0_dp
       vcom(1:3) = 0.0_dp
       vels(:,i_rep)  = 0.0_dp
       DO k=1,nparticle
          vels(k,i_rep) = next_random_number(globenv%gaussian_rng_stream,error=error)
       END DO
       ! Check always if BAND is working in Cartesian or in internal coordinates
       ! If working in cartesian coordinates let's get rid of the COM
       ! Compute also the total mass (both in Cartesian and internal)
       IF (neb_env%use_colvar) THEN
          nvar = nparticle
          mass_tot = REAL(nvar,KIND=dp)*massunit
       ELSE
          DO iatom=1,natom
             mass = particle_set(iatom)%atomic_kind%mass
             mass_tot = mass_tot + mass
             v(1:3) = get_particle_pos_or_vel(iatom,particle_set,vels(:,i_rep))
             vcom(1:3) = vcom(1:3) + mass*v(1:3)
          END DO
          vcom(1:3) = vcom(1:3)/mass_tot
       END IF
       ! Compute kinetic energy and temperature
       akin = 0.0_dp
       IF (neb_env%use_colvar) THEN
          DO ivar=1,nvar
             akin = akin + 0.5_dp*massunit*vels(ivar,i_rep)*vels(ivar,i_rep)
          END DO
       ELSE
          DO iatom=1,natom
             mass = particle_set(iatom)%atomic_kind%mass
             v(1:3) = -vcom(1:3)
             CALL update_particle_pos_or_vel(iatom,particle_set,v(1:3),vels(:,i_rep))
             akin = akin + 0.5_dp*mass*DOT_PRODUCT(v(1:3),v(1:3))
          END DO
          nvar = 3*natom
       END IF
       temp = 2.0_dp*akin/REAL(nvar,KIND=dp)
       ! Scale velocities to get the correct initial temperature and
       sc = SQRT(temp_ext/temp)
       vels(:,i_rep) = vels(:,i_rep)*sc
       ! Re-compute kinetic energya and temperature and velocity of COM
       akin = 0.0_dp
       vcom = 0.0_dp
       IF (neb_env%use_colvar) THEN
          DO ivar=1,nvar
             akin = akin + 0.5_dp*massunit*vels(ivar,i_rep)*vels(ivar,i_rep)
          END DO
       ELSE
          DO iatom=1,natom
             mass = particle_set(iatom)%atomic_kind%mass
             v(1:3) = get_particle_pos_or_vel(iatom,particle_set,vels(:,i_rep))
             vcom(1:3) = vcom(1:3) + mass*v(1:3)
             akin = akin + 0.5_dp*mass*DOT_PRODUCT(v(1:3),v(1:3))
          END DO
       END IF
       vcom(1:3) = vcom(1:3)/mass_tot
       ! Dump information
       IF (iw > 0) THEN
          temp = 2.0_dp*akin/REAL(nvar,KIND=dp)
          tmp_r1 = cp_unit_from_cp2k(temp,"K",error=error)
          WRITE (iw,'(A,T61,F18.2,A2)')&
             ' NEB| Initial Temperature ',tmp_r1," K"
          WRITE (iw,'(A,T61,F20.12)')&
             ' NEB| Centre of mass velocity in direction x:',vcom(1),&
             ' NEB| Centre of mass velocity in direction y:',vcom(2),&
             ' NEB| Centre of mass velocity in direction z:',vcom(3)
          WRITE (iw,'(T2,"NEB|",75("*"))')
       END IF
    ELSE
       vels(:,i_rep) = 0.0_dp
    END IF

  END SUBROUTINE neb_initialize_velocity

! *****************************************************************************
!> \brief Control on  velocities - I part
!> \param vels ...
!> \param particle_set ...
!> \param tc_section ...
!> \param vc_section ...
!> \param output_unit ...
!> \param istep ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE control_vels_a(vels, particle_set, tc_section, vc_section, &
       output_unit, istep, error)
    TYPE(neb_var_type), POINTER              :: vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: tc_section, vc_section
    INTEGER, INTENT(IN)                      :: output_unit, istep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'control_vels_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat, temp_tol_steps
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: ext_temp, f_annealing, scale, &
                                                temp_tol, temploc, tmp_r1, &
                                                tmp_r2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: temperatures

    failure = .FALSE.
    ! Temperature control
    CALL section_vals_get(tc_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(tc_section,"TEMP_TOL_STEPS",i_val=temp_tol_steps, error=error)
       CALL section_vals_val_get(tc_section,"TEMPERATURE",r_val=ext_temp, error=error)
       CALL section_vals_val_get(tc_section,"TEMP_TOL",r_val=temp_tol, error=error)
       ALLOCATE(temperatures(SIZE(vels%wrk,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Computes temperatures
       CALL get_temperatures(vels,particle_set,temperatures,factor=1.0_dp,error=error)
       ! Possibly rescale
       IF (istep<=temp_tol_steps) THEN
          DO i = 2, SIZE(vels%wrk,2)-1
             temploc = temperatures(i)
             IF (ABS(temploc-ext_temp) > temp_tol) THEN
                IF (output_unit>0) THEN
                   tmp_r1 =cp_unit_from_cp2k(temploc,"K",error=error)
                   tmp_r2 =cp_unit_from_cp2k(ext_temp,"K",error=error)
                   WRITE(output_unit,'(T2,"NEB| Replica Nr.",I5,'//&
                        '"  - Velocity rescaled from: ",F12.6," to: ",F12.6,".")')&
                        i, tmp_r1, tmp_r2

                END IF
                scale = SQRT(ext_temp/temploc)
                vels%wrk(:,i) = scale * vels%wrk(:,i)
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
       DO i = 2, SIZE(vels%wrk,2)-1
          vels%wrk(:,i) = f_annealing * vels%wrk(:,i)
       END DO
    END IF
  END SUBROUTINE control_vels_a

! *****************************************************************************
!> \brief Control on velocities - II part
!> \param vels ...
!> \param forces ...
!> \param vc_section ...
!> \param error ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE control_vels_b(vels, forces, vc_section, error)
    TYPE(neb_var_type), POINTER              :: vels, forces
    TYPE(section_vals_type), POINTER         :: vc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'control_vels_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: explicit, lval
    REAL(KIND=dp)                            :: factor, norm

! Check the sign of V.dot.F

    CALL section_vals_get(vc_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(vc_section,"PROJ_VELOCITY_VERLET",l_val=lval, error=error)
       IF (lval) THEN
          DO i = 2, SIZE(vels%wrk,2)-1
             norm   = DOT_PRODUCT(forces%wrk(:,i),forces%wrk(:,i))
             factor = DOT_PRODUCT(vels%wrk(:,i),forces%wrk(:,i))
             IF (factor>0.AND.(norm >= EPSILON(0.0_dp))) THEN
                vels%wrk(:,i) = factor/norm*forces%wrk(:,i)
             ELSE
                vels%wrk(:,i) = 0.0_dp
             END IF
          END DO
       END IF
       CALL section_vals_val_get(vc_section,"SD_LIKE",l_val=lval, error=error)
       IF (lval) THEN
          DO i = 2, SIZE(vels%wrk,2)-1
             vels%wrk(:,i) = 0.0_dp
          END DO
       END IF
    END IF
  END SUBROUTINE control_vels_b

! *****************************************************************************
!> \brief Computes temperatures
!> \param vels ...
!> \param particle_set ...
!> \param temperatures ...
!> \param ekin ...
!> \param factor ...
!> \param error ...
!> \par History
!>      24.11.2010 rewritten to include core-shell model (MK)
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE get_temperatures(vels,particle_set,temperatures,ekin,factor,error)

    TYPE(neb_var_type), POINTER              :: vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: temperatures
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(OUT), OPTIONAL                  :: ekin
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_temperatures', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_rep, iatom, ivar, n_rep, &
                                                natom, nvar
    REAL(KIND=dp)                            :: akin, mass, myfactor, temploc
    REAL(KIND=dp), DIMENSION(3)              :: v

    myfactor = kelvin
    IF (PRESENT(factor)) myfactor = factor
    IF (PRESENT(ekin)) ekin(:) = 0.0_dp
    temperatures(:) = 0.0_dp
    nvar = SIZE(vels%wrk,1)
    n_rep = SIZE(vels%wrk,2)
    natom = SIZE(particle_set)
    DO i_rep=2,n_rep-1
      akin = 0.0_dp
      IF (vels%in_use == do_band_collective) THEN
        DO ivar=1,nvar
          akin = akin + 0.5_dp*massunit*vels%wrk(ivar,i_rep)*vels%wrk(ivar,i_rep)
        END DO
      ELSE
        DO iatom=1,natom
          mass = particle_set(iatom)%atomic_kind%mass
          v(1:3) = get_particle_pos_or_vel(iatom,particle_set,vels%wrk(:,i_rep))
          akin = akin + 0.5_dp*mass*DOT_PRODUCT(v(1:3),v(1:3))
        END DO
        nvar = 3*natom
      END IF
      IF (PRESENT(ekin)) ekin(i_rep) = akin
      temploc = 2.0_dp*akin/REAL(nvar,KIND=dp)
      temperatures(i_rep) = temploc*myfactor
    END DO

  END SUBROUTINE get_temperatures

END MODULE neb_md_utils
