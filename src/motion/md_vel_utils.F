!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Collection of utilities for setting-up and handle velocities in MD
!>         runs
!> \author CJM
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
!>      Teodoro Laino [tlaino] 12.2008 - Preparing for VIRTUAL SITE constraints
!>                                       (patch by Marcel Baer)
! *****************************************************************************
MODULE md_vel_utils
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: &
       cell_type, use_perd_none, use_perd_x, use_perd_xy, use_perd_xyz, &
       use_perd_xz, use_perd_y, use_perd_yz, use_perd_z
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE force_env_utils,                 ONLY: force_env_rattle,&
                                             force_env_shake
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: npe_f_ensemble,&
                                             npe_i_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble,&
                                             reftraj_ensemble
  USE input_cp2k_binary_restarts,      ONLY: read_binary_velocities
  USE input_restart_force_eval,        ONLY: update_subsys
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: diamat_all
  USE md_ener_types,                   ONLY: md_ener_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             get_molecule_kind_set,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: kelvin
  USE shell_opt,                       ONLY: optimize_shell_core
  USE shell_potential_types,           ONLY: shell_kind_type
  USE simpar_types,                    ONLY: simpar_type
  USE thermal_region_types,            ONLY: thermal_region_type,&
                                             thermal_regions_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_vel_utils'


  PUBLIC :: temperature_control,&
            comvel_control,&
            angvel_control,&
            setup_velocities

CONTAINS

! *****************************************************************************
!> \brief compute center of mass position
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param is_fixed ...
!> \param rcom ...
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE compute_rcom(part,is_fixed,rcom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: rcom

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_rcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: denom, mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    rcom(:)  = 0.0_dp
    denom = 0.0_dp
    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_xz,use_perd_yz,use_perd_none)
          rcom(1) = rcom(1) + part(i)%r(1) * mass
          rcom(2) = rcom(2) + part(i)%r(2) * mass
          rcom(3) = rcom(3) + part(i)%r(3) * mass
          denom = denom + mass
       END SELECT
    END DO
    rcom = rcom/denom

  END SUBROUTINE compute_rcom

! *****************************************************************************
!> \brief compute center of mass velocity
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param is_fixed ...
!> \param vcom ...
!> \param ecom ...
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE compute_vcom(part,is_fixed,vcom,ecom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: vcom
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: ecom

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_vcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: denom, mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    vcom  = 0.0_dp
    denom = 0.0_dp
    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
       IF (mass.NE.0.0) THEN
          SELECT CASE(is_fixed(i))
          CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_xz,use_perd_yz,use_perd_none)
             vcom(1) = vcom(1) + part(i)%v(1) * mass
             vcom(2) = vcom(2) + part(i)%v(2) * mass
             vcom(3) = vcom(3) + part(i)%v(3) * mass
             denom = denom + mass
          END SELECT
       END IF
    END DO
    vcom = vcom/denom
    IF (PRESENT(ecom)) THEN
       ecom = 0.5_dp*denom*SUM(vcom*vcom)
    END IF

  END SUBROUTINE compute_vcom

! *****************************************************************************
!> \brief Copy atom velocities into core and shell velocities
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param shell_part ...
!> \param core_part ...
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE clone_core_shell_vel(part,shell_part,core_part)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part, shell_part, core_part

    CHARACTER(len=*), PARAMETER :: routineN = 'clone_core_shell_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: is_shell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, shell_active=is_shell)
       IF (is_shell) THEN
          shell_part( part(i)%shell_index )%v(:) =  part(i)%v(:)
          core_part( part(i)%shell_index )%v(:)  =  part(i)%v(:)
       END IF
    END DO

  END SUBROUTINE clone_core_shell_vel

! *****************************************************************************
!> \brief Compute the kinetic energy. Does not subtract the center of mass kinetic
!>      energy.
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param ireg ...
!> \retval ekin ...
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  FUNCTION compute_ekin(part,ireg) RESULT(ekin)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, INTENT(IN), OPTIONAL            :: ireg
    REAL(KIND=dp)                            :: ekin

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_ekin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    NULLIFY(atomic_kind)
    ekin = 0.0_dp
    IF(PRESENT(ireg)) THEN
      DO i = 1, SIZE(part)
       IF(part(i)%t_region_index==ireg) THEN
         atomic_kind => part(i)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
         ekin = ekin + 0.5_dp * mass * SUM(part(i)%v(:) * part(i)%v(:))
       END IF
      END DO
    ELSE
      DO i = 1, SIZE(part)
         atomic_kind => part(i)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
         ekin = ekin + 0.5_dp * mass * SUM(part(i)%v(:) * part(i)%v(:))
      END DO
    END IF

  END FUNCTION compute_ekin

! *****************************************************************************
!> \brief Rescale the velocity to mimic the given external kinetic temperature.
!>      Optionally also rescale vcom.
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param simpar ...
!> \param ekin ...
!> \param vcom ...
!> \param ireg ...
!> \param nfree ...
!> \param temp ...
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE rescale_vel(part,simpar,ekin,vcom,ireg,nfree,temp)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(simpar_type), POINTER               :: simpar
    REAL(KIND=dp), INTENT(INOUT)             :: ekin
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT), OPTIONAL                :: vcom
    INTEGER, INTENT(IN), OPTIONAL            :: ireg, nfree
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: temp

    CHARACTER(len=*), PARAMETER :: routineN = 'rescale_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_ireg, my_nfree
    REAL(KIND=dp)                            :: factor, my_temp

    IF(PRESENT(ireg).AND.PRESENT(nfree).AND.PRESENT(temp)) THEN
      my_ireg = ireg
      my_nfree = nfree
      my_temp = temp
    ELSEIF(PRESENT(nfree)) THEN
      my_ireg = 0
      my_nfree = nfree
      my_temp = simpar%temp_ext
    ELSE
      my_ireg = 0
      my_nfree = simpar%nfree
      my_temp = simpar%temp_ext
    END IF
    IF (my_nfree/=0) THEN
       factor = my_temp / ( 2.0_dp * ekin ) * REAL(my_nfree,KIND=dp)
    ELSE
       factor = 0.0_dp
    ENDIF
    ! Note:
    ! this rescaling is still wrong, it should take the masses into account
    ! rescaling is generally not correct, so needs fixing
    ekin = ekin * factor
    factor = SQRT(factor)
    IF(PRESENT(ireg)) THEN
      DO i = 1, SIZE(part)
        IF(part(i)%t_region_index == my_ireg) part(i)%v(:) = factor*part(i)%v(:)
      END DO
    ELSE
      DO i = 1, SIZE(part)
        part(i)%v(:) = factor*part(i)%v(:)
      END DO
      IF (PRESENT(vcom)) THEN
         vcom = factor*vcom
      END IF
    END IF

  END SUBROUTINE rescale_vel

! *****************************************************************************
!> \brief Rescale the velocity of separated regions independently
!> \param part ...
!> \param md_env ...
!> \param simpar ...
!> \param error ...
!> \par History
!>      2008-11
!> \author  MI
! *****************************************************************************
  SUBROUTINE rescale_vel_region(part,md_env,simpar,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rescale_vel_region', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ireg, nfree, nfree0, &
                                                nfree_done
    REAL(KIND=dp)                            :: ekin, temp
    TYPE(thermal_region_type), POINTER       :: t_region
    TYPE(thermal_regions_type), POINTER      :: thermal_regions

    NULLIFY(thermal_regions, t_region)

    CALL get_md_env(md_env,thermal_regions=thermal_regions,error=error)
    nfree_done = 0
    DO ireg = 1,thermal_regions%nregions
      NULLIFY(t_region)
      t_region => thermal_regions%thermal_region(ireg)
      nfree = t_region%npart*3
      ekin = compute_ekin(part,ireg)
      temp =  t_region%temp_expected
      CALL rescale_vel(part,simpar,ekin,ireg=ireg,nfree=nfree,temp=temp)
      nfree_done = nfree_done+nfree
      ekin = compute_ekin(part,ireg)
      temp = 2.0_dp* ekin/REAL(nfree,dp)*kelvin
      t_region%temperature = temp
    END DO
    nfree0 = simpar%nfree - nfree_done
    IF(nfree0>0) THEN
      ekin = compute_ekin(part,0)
      CALL rescale_vel(part,simpar,ekin,ireg=0,nfree=nfree0,temp=simpar%temp_ext)
      ekin = compute_ekin(part,0)
      temp = 2.0_dp* ekin/REAL(nfree0,dp)*kelvin
      thermal_regions%temp_reg0 = temp
    END IF
  END SUBROUTINE  rescale_vel_region

! *****************************************************************************
!> \brief subtract center of mass velocity
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param is_fixed ...
!> \param vcom ...
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE subtract_vcom(part,is_fixed,vcom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vcom

    CHARACTER(len=*), PARAMETER :: routineN = 'subtract_vcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    DO i = 1, SIZE(part)
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x)
          part(i)%v(2) = part(i)%v(2) - vcom(2)
          part(i)%v(3) = part(i)%v(3) - vcom(3)
       CASE(use_perd_y)
          part(i)%v(1) = part(i)%v(1) - vcom(1)
          part(i)%v(3) = part(i)%v(3) - vcom(3)
       CASE(use_perd_z)
          part(i)%v(1) = part(i)%v(1) - vcom(1)
          part(i)%v(2) = part(i)%v(2) - vcom(2)
       CASE(use_perd_xy)
          part(i)%v(3) = part(i)%v(3) - vcom(3)
       CASE(use_perd_xz)
          part(i)%v(2) = part(i)%v(2) - vcom(2)
       CASE(use_perd_yz)
          part(i)%v(1) = part(i)%v(1) - vcom(1)
       CASE(use_perd_none)
          part(i)%v(:) = part(i)%v(:) - vcom(:)
       END SELECT
    END DO
  END SUBROUTINE subtract_vcom

! *****************************************************************************
!> \brief compute the angular velocity
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param is_fixed ...
!> \param rcom ...
!> \param vang ...
!> \param error ...
!> \par History
!>      2007-11-9: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE compute_vang(part,is_fixed,rcom,vang,error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rcom
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: vang
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_vang', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: mass, proj
    REAL(KIND=dp), DIMENSION(3)              :: evals, mang, r
    REAL(KIND=dp), DIMENSION(3, 3)           :: iner
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    NULLIFY(atomic_kind)
    mang(:)   = 0.0_dp
    iner(:,:) = 0.0_dp
    DO i=1,SIZE(part)
       ! compute angular momentum and inertia tensor
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_xz,use_perd_yz,use_perd_none)
          r(:) = part(i)%r(:) - rcom(:)
          atomic_kind => part(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
          mang(1) = mang(1) + mass*(r(2)*part(i)%v(3) - r(3)*part(i)%v(2))
          mang(2) = mang(2) + mass*(r(3)*part(i)%v(1) - r(1)*part(i)%v(3))
          mang(3) = mang(3) + mass*(r(1)*part(i)%v(2) - r(2)*part(i)%v(1))

          iner(1,1) = iner(1,1) + mass*(r(2)*r(2) + r(3)*r(3))
          iner(2,2) = iner(2,2) + mass*(r(3)*r(3) + r(1)*r(1))
          iner(3,3) = iner(3,3) + mass*(r(1)*r(1) + r(2)*r(2))

          iner(1,2) = iner(1,2) - mass*r(1)*r(2)
          iner(2,3) = iner(2,3) - mass*r(2)*r(3)
          iner(3,1) = iner(3,1) - mass*r(3)*r(1)
       END SELECT
    END DO
    iner(2,1) = iner(1,2)
    iner(3,2) = iner(2,3)
    iner(1,3) = iner(3,1)

    ! Take the safest route, i.e. diagonalize the inertia tensor and solve
    ! the angular velocity only with the non-zero eigenvalues. A plain inversion
    ! would fail for linear molecules.
    CALL diamat_all(iner, evals,error=error)

    vang(:) = 0.0_dp
    DO i=1,3
       IF (evals(i) > 0.0_dp) THEN
         proj = SUM(iner(:,i)*mang)/evals(i)
         vang(1) = vang(1) + proj*iner(1,i)
         vang(2) = vang(2) + proj*iner(2,i)
         vang(3) = vang(3) + proj*iner(3,i)
       END IF
    END DO

  END SUBROUTINE compute_vang

! *****************************************************************************
!> \brief subtract the angular velocity
!>      *** is only used by initialize_velocities below ***
!> \param part ...
!> \param is_fixed ...
!> \param rcom ...
!> \param vang ...
!> \par History
!>      2007-11-9: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE subtract_vang(part,is_fixed,rcom,vang)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    INTEGER, DIMENSION(:), INTENT(IN)        :: is_fixed
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rcom, vang

    CHARACTER(len=*), PARAMETER :: routineN = 'subtract_vang', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(3)              :: r

    DO i=1,SIZE(part)
       r(:) = part(i)%r(:) - rcom(:)
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x)
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       CASE(use_perd_y)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       CASE(use_perd_z)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
       CASE(use_perd_xy)
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       CASE(use_perd_xz)
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
       CASE(use_perd_yz)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
       CASE(use_perd_none)
          part(i)%v(1) = part(i)%v(1) - (vang(2)*r(3) - vang(3)*r(2))
          part(i)%v(2) = part(i)%v(2) - (vang(3)*r(1) - vang(1)*r(3))
          part(i)%v(3) = part(i)%v(3) - (vang(1)*r(2) - vang(2)*r(1))
       END SELECT
    END DO

  END SUBROUTINE subtract_vang

! *****************************************************************************
!> \brief Initializes the velocities to the Maxwellian distribution
!> \param simpar ...
!> \param part ...
!> \param force_env ...
!> \param globenv ...
!> \param md_env ...
!> \param molecule_kinds ...
!> \param label ...
!> \param print_section ...
!> \param subsys_section ...
!> \param shell_present ...
!> \param shell_part ...
!> \param core_part ...
!> \param force_rescaling ...
!> \param para_env ...
!> \param error ...
!> \par History
!>      - is_fixed removed from particle_type
!>      - 2007-11-07: Cleanup (TV)
!>      - 2007-11-09: Added angvel_zero feature
!> \author CJM,MK,Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE initialize_velocities(simpar, part, force_env, globenv, md_env, &
       molecule_kinds, label, print_section, subsys_section, shell_present,&
       shell_part, core_part, force_rescaling, para_env, error)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    CHARACTER(LEN=*), INTENT(IN)             :: label
    TYPE(section_vals_type), POINTER         :: print_section, subsys_section
    LOGICAL, INTENT(IN)                      :: shell_present
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_part, core_part
    LOGICAL, INTENT(IN)                      :: force_rescaling
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_velocities', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_format
    INTEGER                                  :: handle, i, ifixd, &
                                                imolecule_kind, iw, natoms, &
                                                num_x, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: is_fixed
    LOGICAL                                  :: failure, success
    REAL(KIND=dp)                            :: ecom, ekin, mass, mass_tot, &
                                                temp, tmp_r1
    REAL(KIND=dp), DIMENSION(3)              :: rcom, vang, vcom
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    CALL timeset(routineN,handle)

    ! Initializing parameters
    failure = .FALSE.
    natoms = SIZE(part)
    NULLIFY (atomic_kind, fixd_list, logger, molecule_kind)
    NULLIFY (molecule_kind_set)

    ! Logging
    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",extension=".log",error=error)
    IF (iw>0) THEN
       num_x = (79-LEN_TRIM(ADJUSTL(label))-2)/2
       WRITE(my_format,'(A,I0,A,I0,A)')'(1X,',num_x,'("*"),1X,A,1X,',num_x,'("*"))'
       WRITE ( iw, TRIM(my_format))TRIM(ADJUSTL(label))
    END IF

    ! Build a list of all fixed atoms (if any)
    ALLOCATE (is_fixed(natoms),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    is_fixed = use_perd_none
    molecule_kind_set => molecule_kinds%els
    DO imolecule_kind=1,molecule_kinds%n_els
       molecule_kind => molecule_kind_set(imolecule_kind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,fixd_list=fixd_list)
       IF (ASSOCIATED(fixd_list)) THEN
          DO ifixd=1,SIZE(fixd_list)
             IF (.NOT.fixd_list(ifixd)%restraint%active) is_fixed(fixd_list(ifixd)%fixd) = fixd_list(ifixd)%itype
          END DO
       END IF
    END DO

    ! Compute the total mass when needed
    IF ( simpar%ensemble == nph_uniaxial_ensemble .OR.&
         simpar%ensemble == nph_uniaxial_damped_ensemble ) THEN
       mass_tot = 0.0_dp
       DO i = 1, natoms
          atomic_kind => part(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
          mass_tot = mass_tot + mass
       END DO
       simpar%v_shock = simpar%v_shock * SQRT ( mass_tot )
    END IF


    CALL read_input_velocities(simpar, part, force_env, md_env, subsys_section,&
             shell_present, shell_part, core_part, force_rescaling, para_env, is_fixed, success, error)
    IF(.NOT.success) &
       CALL generate_velocities(simpar, part, force_env, globenv, md_env, shell_present,&
                shell_part, core_part, is_fixed, iw, error)


    IF (iw>0) THEN
       ! Recompute vcom, ecom and ekin for IO
       CALL compute_vcom(part,is_fixed,vcom,ecom)
       ekin = compute_ekin(part) - ecom
       IF (simpar%nfree == 0) THEN
          CPPostcondition(ekin==0.0_dp,cp_failure_level,routineP,error,failure)
          temp = 0.0_dp
       ELSE
          temp = 2.0_dp * ekin / REAL ( simpar%nfree,KIND=dp)
       END IF
       tmp_r1 = cp_unit_from_cp2k(temp,"K",error=error)
       WRITE (iw, '( A, T61, F18.2, A2 )' ) ' Initial Temperature ', tmp_r1, " K"
       WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' COM velocity:', vcom ( 1 ), vcom ( 2 ), vcom ( 3 )

       ! compute and log rcom and vang if not periodic
       CALL force_env_get(force_env, cell=cell, error=error)
       IF (SUM(cell%perd(1:3)) == 0) THEN
          CALL compute_rcom(part,is_fixed,rcom)
          CALL compute_vang(part,is_fixed,rcom,vang,error)
          WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' COM position:', rcom ( 1 ), rcom ( 2 ), rcom ( 3 )
          WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' Angular velocity:', vang ( 1 ), vang ( 2 ), vang ( 3 )
       END IF
       WRITE ( iw, '( 1X,79("*"),/)' )
    END IF

    DEALLOCATE (is_fixed,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_print_key_finished_output(iw,logger,print_section,"PROGRAM_RUN_INFO", error=error)
    CALL timestop(handle)

  END SUBROUTINE initialize_velocities


! *****************************************************************************
!> \brief Read velocities from binary restart file if available
!> \param simpar ...
!> \param part ...
!> \param force_env ...
!> \param md_env ...
!> \param subsys_section ...
!> \param shell_present ...
!> \param shell_part ...
!> \param core_part ...
!> \param force_rescaling ...
!> \param para_env ...
!> \param is_fixed ...
!> \param success ...
!> \param error ...
!> \author CJM,MK,Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE read_input_velocities(simpar, part, force_env, md_env, subsys_section,&
       shell_present, shell_part, core_part, force_rescaling, para_env, is_fixed, success, error)
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    LOGICAL, INTENT(IN)                      :: shell_present
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_part, core_part
    LOGICAL, INTENT(IN)                      :: force_rescaling
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: is_fixed
    LOGICAL, INTENT(OUT)                     :: success
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_input_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, natoms, nshell, shell_index
    LOGICAL :: atomvel_explicit, atomvel_read, corevel_explicit, &
      corevel_read, failure, is_ok, rescale_regions, shellvel_explicit, &
      shellvel_read
    REAL(KIND=dp)                            :: ecom, ekin, fac_massc, &
                                                fac_masss, mass
    REAL(KIND=dp), DIMENSION(3)              :: vc, vcom, vs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_sll_val_type), POINTER           :: atom_list, core_list, &
                                                shell_list
    TYPE(section_vals_type), POINTER         :: atomvel_section, &
                                                corevel_section, &
                                                shellvel_section
    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(val_type), POINTER                  :: val

! Initializing parameters

    success = .FALSE.
    natoms = SIZE(part)
    atomvel_read = .FALSE.
    corevel_read = .FALSE.
    shellvel_read = .FALSE.
    NULLIFY (vel, atomic_kind, atom_list, core_list, shell_list)
    NULLIFY (atomvel_section,shellvel_section, corevel_section)
    NULLIFY (shell, thermal_regions, val)


    ! Core-Shell Model
    nshell = 0
    IF (shell_present) THEN
       CPPostcondition(ASSOCIATED(core_part),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(shell_part),cp_failure_level,routineP,error,failure)
       nshell = SIZE(shell_part)
    END IF

    atomvel_section  => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
    shellvel_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
    corevel_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)

    ! Read or initialize the particle velocities
    CALL section_vals_get(atomvel_section,explicit=atomvel_explicit,error=error)
    CALL section_vals_get(shellvel_section,explicit=shellvel_explicit,error=error)
    CALL section_vals_get(corevel_section,explicit=corevel_explicit,error=error)
    CPPostcondition(shellvel_explicit.EQV.corevel_explicit,cp_failure_level,routineP,error,failure)

    CALL read_binary_velocities("",part,force_env%root_section,para_env,&
                                subsys_section,atomvel_read,error)
    CALL read_binary_velocities("SHELL",shell_part,force_env%root_section,para_env,&
                                subsys_section,shellvel_read,error)
    CALL read_binary_velocities("CORE",core_part,force_env%root_section,para_env,&
                                subsys_section,corevel_read,error)

    IF (.NOT.(atomvel_explicit.OR.atomvel_read)) RETURN
    success = .TRUE.

     IF (.NOT.atomvel_read) THEN
        ! Read the atom velocities if explicitly given in the input file
        CALL section_vals_list_get(atomvel_section,"_DEFAULT_KEYWORD_",list=atom_list,error=error)
        DO i = 1, natoms
           is_ok = cp_sll_val_next(atom_list,val,error=error)
           CALL val_get(val,r_vals=vel,error=error)
           part(i)%v = vel
        END DO
     END IF
     DO i = 1, natoms
       SELECT CASE(is_fixed(i))
       CASE (use_perd_x)
          part(i)%v(1) = 0.0_dp
       CASE (use_perd_y)
          part(i)%v(2) = 0.0_dp
       CASE (use_perd_z)
          part(i)%v(3) = 0.0_dp
       CASE (use_perd_xy)
          part(i)%v(1) = 0.0_dp
          part(i)%v(2) = 0.0_dp
       CASE (use_perd_xz)
          part(i)%v(1) = 0.0_dp
          part(i)%v(3) = 0.0_dp
       CASE (use_perd_yz)
          part(i)%v(2) = 0.0_dp
          part(i)%v(3) = 0.0_dp
       CASE (use_perd_xyz)
          part(i)%v    = 0.0_dp
       END SELECT
     END DO
     IF (shell_present) THEN
        IF (shellvel_explicit) THEN
           ! If the atoms positions are given (?) and core and shell velocities are
           ! present in the input, read the latter.
           CALL section_vals_list_get(shellvel_section,"_DEFAULT_KEYWORD_",list=shell_list,error=error)
           CALL section_vals_list_get(corevel_section,"_DEFAULT_KEYWORD_",list=core_list,error=error)
           DO i=1,nshell
              is_ok = cp_sll_val_next(shell_list,val,error=error)
              CALL val_get(val,r_vals=vel,error=error)
              shell_part(i)%v = vel
              is_ok = cp_sll_val_next(core_list,val,error=error)
              CALL val_get(val,r_vals=vel,error=error)
              core_part(i)%v = vel
           END DO
        ELSE
           IF (.NOT.(shellvel_read.AND.corevel_read)) THEN
              ! Otherwise, just copy atom velocties into shell and core velocities.
              CALL clone_core_shell_vel(part,shell_part,core_part)
           END IF
        END IF
     END IF

     ! compute vcom, ecom and ekin
     CALL compute_vcom(part,is_fixed,vcom,ecom)
     ekin = compute_ekin(part) - ecom

     IF(simpar%do_thermal_region) THEN
       CALL get_md_env (md_env, thermal_regions=thermal_regions, error=error)
       IF(ASSOCIATED(thermal_regions)) THEN
         rescale_regions = thermal_regions%force_rescaling
       END IF
     ELSE
       rescale_regions = .FALSE.
     END IF
     IF (simpar%nfree /= 0 .AND. (force_rescaling .OR. rescale_regions)) THEN
        IF(simpar%do_thermal_region) THEN
          CALL rescale_vel_region(part,md_env,simpar,error=error)
        ELSE
          CALL rescale_vel(part,simpar,ekin,vcom=vcom)
        END IF

        ! After rescaling, the core and shell velocities must also adapt.
        DO i = 1, natoms
          shell_index = part(i)%shell_index
          IF(shell_present .AND. shell_index/=0) THEN
              atomic_kind => part(i)%atomic_kind
              CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass, shell=shell)
              fac_masss = shell%mass_shell/mass
              fac_massc = shell%mass_core/mass
              vs = shell_part(shell_index)%v
              vc = core_part(shell_index)%v

              shell_part(shell_index)%v(1) = part(i)%v(1) + fac_massc*(vs(1)-vc(1))
              shell_part(shell_index)%v(2) = part(i)%v(2) + fac_massc*(vs(2)-vc(2))
              shell_part(shell_index)%v(3) = part(i)%v(3) + fac_massc*(vs(3)-vc(3))
              core_part(shell_index)%v(1)  = part(i)%v(1) + fac_masss*(vc(1)-vs(1))
              core_part(shell_index)%v(2)  = part(i)%v(2) + fac_masss*(vc(2)-vs(2))
              core_part(shell_index)%v(3)  = part(i)%v(3) + fac_masss*(vc(3)-vs(3))
          END IF
        END DO
     END IF
  END SUBROUTINE read_input_velocities


! *****************************************************************************
!> \brief Initializing velocities deterministically on all processors, if not given in input
!> \param simpar ...
!> \param part ...
!> \param force_env ...
!> \param globenv ...
!> \param md_env ...
!> \param shell_present ...
!> \param shell_part ...
!> \param core_part ...
!> \param is_fixed ...
!> \param iw ...
!> \param error ...
!> \author CJM,MK,Toon Verstraelen <Toon.Verstraelen@gmail.com>, Ole Schuett
! *****************************************************************************
  SUBROUTINE generate_velocities(simpar, part, force_env, globenv, md_env,&
                          shell_present, shell_part, core_part, is_fixed, iw, error)
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(md_environment_type), POINTER       :: md_env
    LOGICAL, INTENT(IN)                      :: shell_present
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_part, core_part
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: is_fixed
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, natoms
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    NULLIFY (atomic_kind)
    natoms = SIZE(part)

     DO i = 1, natoms
        atomic_kind => part(i)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
        part(i)%v(1) = 0.0_dp
        part(i)%v(2) = 0.0_dp
        part(i)%v(3) = 0.0_dp
        IF (mass.NE.0.0) THEN
           SELECT CASE(is_fixed(i))
           CASE (use_perd_x)
              part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
              part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
           CASE (use_perd_y)
              part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
              part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
           CASE (use_perd_z)
              part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
              part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
           CASE (use_perd_xy)
              part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
           CASE (use_perd_xz)
              part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
           CASE (use_perd_yz)
              part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
           CASE (use_perd_none)
              part(i)%v(1) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
              part(i)%v(2) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
              part(i)%v(3) = next_random_number(globenv%gaussian_rng_stream,error=error) / SQRT(mass)
           END SELECT
        END IF
     END DO

     CALL normalize_velocities(simpar, part, force_env, md_env, is_fixed, error)
     CALL soften_velocities(simpar, part, force_env, md_env, is_fixed, iw, error)

     ! Initialize the core and the shell velocity. Atom velocities are just
     ! copied so that the initial relative core-shell velocity is zero.
     IF (shell_present) THEN
       CALL optimize_shell_core(force_env,part,shell_part, core_part, globenv, error=error)
     ENDIF
  END SUBROUTINE generate_velocities

! *****************************************************************************
!> \brief Direct velocities along a low-curvature direction in order to 
!>        favors MD trajectories to cross rapidly over small energy barriers
!>        into neighboring basins.
!> \param simpar ...
!> \param part ...
!> \param force_env ...
!> \param md_env ...
!> \param is_fixed ...
!> \param iw ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE soften_velocities(simpar, part, force_env, md_env, is_fixed, iw, error)
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: is_fixed
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'soften_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k
    REAL(KIND=dp), DIMENSION(SIZE(part), 3)  :: F, F_t, N, x0

    IF(simpar%soften_nsteps <= 0) RETURN !nothing todo

    CALL cp_assert(ALL(is_fixed==use_perd_none), cp_fatal_level, cp_assertion_failed,&
                 routineP,"Velocitiy softening with constraints is not supported.")

    !backup positions
    DO i = 1, SIZE(part)
          x0(i,:) = part(i)%r
    END DO

    DO k=1, simpar%soften_nsteps

       !use normalized velocities as displace direction
       DO i = 1, SIZE(part)
          N(i,:) = part(i)%v
       END DO
       N = N / SQRT(SUM(N**2))

       ! displace system temporarly to calculate forces
       DO i = 1, SIZE(part)
          part(i)%r = part(i)%r + simpar%soften_delta * N(i,:)
       END DO
       CALL force_env_calc_energy_force(force_env,error=error)

       ! calculate velocity update direction F_t
       DO i = 1, SIZE(part)
          F(i,:) = part(i)%f
       END DO
       F_t = F - N * SUM(N*F)

       ! restore positions and update velocities
       DO i = 1, SIZE(part)
          part(i)%r = x0(i,:)
          part(i)%v = part(i)%v + simpar%soften_alpha * F_t(i,:)
       END DO

       CALL normalize_velocities(simpar, part, force_env, md_env, is_fixed, error)
    END DO

    IF(iw>0) THEN
      WRITE (iw, "(A,T71, I10)") " Velocities softening Steps: ",simpar%soften_nsteps
      WRITE (iw, "(A,T71, E10.3)") " Velocities softening NORM(F_t): ",SQRT(SUM(F_t**2))
    END IF
  END SUBROUTINE soften_velocities


! *****************************************************************************
!> \brief Scale velocities according to temperature and remove rigid body motion.
!> \param simpar ...
!> \param part ...
!> \param force_env ...
!> \param md_env ...
!> \param is_fixed ...
!> \param error ...
!> \author CJM,MK,Toon Verstraelen <Toon.Verstraelen@gmail.com>, Ole Schuett
! *****************************************************************************
  SUBROUTINE normalize_velocities(simpar, part, force_env, md_env, is_fixed, error)
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_environment_type), POINTER       :: md_env
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: is_fixed
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'normalize_velocities', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: ekin
    REAL(KIND=dp), DIMENSION(3)              :: rcom, vang, vcom
    TYPE(cell_type), POINTER                 :: cell

     NULLIFY(cell)

     ! Subtract the vcom
     CALL compute_vcom(part,is_fixed,vcom)
     CALL subtract_vcom(part,is_fixed,vcom)
     ! If requested and the system is not periodic, subtract the angular velocity
     CALL force_env_get(force_env, cell=cell, error=error)
     IF (SUM(cell%perd(1:3)) == 0 .AND. simpar%angvel_zero) THEN
        CALL compute_rcom(part,is_fixed,rcom)
        CALL compute_vang(part,is_fixed,rcom,vang,error)
        CALL subtract_vang(part,is_fixed,rcom,vang)
     END IF
     ! Rescale the velocities
     IF(simpar%do_thermal_region) THEN
       CALL rescale_vel_region(part,md_env,simpar,error=error)
     ELSE
       ekin = compute_ekin(part)
       CALL rescale_vel(part,simpar,ekin)
     END IF
  END SUBROUTINE normalize_velocities


! *****************************************************************************
!> \brief Computes Ekin, VCOM and Temp for particles
!> \param subsys ...
!> \param md_ener ...
!> \param vsubtract ...
!> \param error ...
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE reset_vcom(subsys, md_ener, vsubtract, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: vsubtract
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reset_vcom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, handle, iatom, ikind, &
                                                natom, shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: ekin_old, imass_c, imass_s, &
                                                mass, v2
    REAL(KIND=dp), DIMENSION(3)              :: tmp, v
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(particles, atomic_kind, atomic_kinds, atom_list, shell)
    CALL timeset(routineN,handle)

    CALL cp_subsys_get(subsys,&
                       atomic_kinds=atomic_kinds,&
                       particles=particles,&
                       shell_particles=shell_particles,&
                       core_particles=core_particles,&
                       error=error)

    ekin_old = md_ener%ekin
    ! Possibly subtract a quantity from all velocities
       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list,&
             natom=natom, mass=mass, shell_active=is_shell, shell=shell)
          IF (is_shell) THEN
             tmp = 0.5_dp*vsubtract*mass
             imass_s = 1.0_dp/shell%mass_shell
             imass_c = 1.0_dp/shell%mass_core
             DO iatom = 1, natom
                atom = atom_list(iatom)
                shell_index = particles%els(atom)%shell_index
                shell_particles%els(shell_index)%v = shell_particles%els(shell_index)%v - tmp*imass_s
                core_particles%els(shell_index)%v  = core_particles%els(shell_index)%v  - tmp*imass_c
                particles%els(atom)%v = particles%els(atom)%v - vsubtract
             END DO
          ELSE
             DO iatom = 1, natom
                atom = atom_list(iatom)
                particles%els(atom)%v = particles%els(atom)%v - vsubtract
             END DO
          END IF
       END DO
    ! Compute Kinetic Energy and COM Velocity
    md_ener%vcom       = 0.0_dp
    md_ener%total_mass = 0.0_dp
    md_ener%ekin       = 0.0_dp
    DO ikind=1,atomic_kinds%n_els
       atomic_kind => atomic_kinds%els(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list, mass=mass, natom=natom)
          v2   = 0.0_dp
          v    = 0.0_dp
          DO iatom = 1, natom
             atom = atom_list(iatom)
             v2 = v2 + SUM(particles%els(atom)%v**2)
             v(1) = v(1) + particles%els(atom)%v(1)
             v(2) = v(2) + particles%els(atom)%v(2)
             v(3) = v(3) + particles%els(atom)%v(3)
          END DO
          md_ener%ekin  = md_ener%ekin + 0.5_dp*mass*v2
          md_ener%vcom(1)    = md_ener%vcom(1) + mass*v(1)
          md_ener%vcom(2)    = md_ener%vcom(2) + mass*v(2)
          md_ener%vcom(3)    = md_ener%vcom(3) + mass*v(3)
          md_ener%total_mass = md_ener%total_mass + REAL(natom,KIND=dp)*mass
    END DO
    md_ener%vcom = md_ener%vcom / md_ener%total_mass
    md_ener%constant = md_ener%constant - ekin_old + md_ener%ekin
    IF (md_ener%nfree /=0) THEN
         md_ener%temp_part = 2.0_dp*md_ener%ekin/REAL(md_ener%nfree,KIND=dp)*kelvin
    END IF
    CALL timestop(handle)

  END SUBROUTINE reset_vcom

! *****************************************************************************
!> \brief Scale velocities to get the correct temperature
!> \param subsys ...
!> \param md_ener ...
!> \param temp_expected ...
!> \param temp_tol ...
!> \param iw ...
!> \param error ...
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE scale_velocity(subsys, md_ener, temp_expected, temp_tol, iw, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: ekin_old, scale, temp_old

    IF (ABS(temp_expected - md_ener%temp_part/kelvin) > temp_tol) THEN
       scale = 0.0_dp
       IF (md_ener%temp_part>0.0_dp) scale = SQRT((temp_expected/md_ener%temp_part)*kelvin)
       ekin_old = md_ener%ekin
       temp_old = md_ener%temp_part
       md_ener%ekin       = 0.0_dp
       md_ener%temp_part  = 0.0_dp
       md_ener%vcom       = 0.0_dp
       md_ener%total_mass = 0.0_dp

       CALL scale_velocity_low(subsys,scale,ireg=0,ekin=md_ener%ekin,vcom=md_ener%vcom,error=error)
       IF(md_ener%nfree /=0) THEN
           md_ener%temp_part = 2.0_dp*md_ener%ekin/REAL(md_ener%nfree,KIND=dp)*kelvin
       END IF
       md_ener%constant = md_ener%constant - ekin_old + md_ener%ekin

       IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")"Temperature scaled to requested temperature:",&
             temp_old," K  ->",md_ener%temp_part," K"
       END IF
    END IF
  END SUBROUTINE scale_velocity

! *****************************************************************************
!> \brief Scale velocities of set of regions
!> \param md_env ...
!> \param subsys ...
!> \param md_ener ...
!> \param simpar ...
!> \param iw ...
!> \param error ...
!> \par author MI
! *****************************************************************************
  SUBROUTINE scale_velocity_region(md_env, subsys, md_ener, simpar, iw, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(simpar_type), POINTER               :: simpar
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity_region', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ireg, istat, nfree, &
                                                nfree_done, nregions
    REAL(KIND=dp)                            :: ekin, ekin_old, &
                                                ekin_total_new, fscale, &
                                                vcom(3), vcom_total(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: temp_new, temp_old
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(thermal_region_type), POINTER       :: t_region
    TYPE(thermal_regions_type), POINTER      :: thermal_regions

     NULLIFY( particles, part, thermal_regions, t_region)
     CALL cp_subsys_get(subsys, particles=particles, error=error)
     part => particles%els
     CALL get_md_env(md_env, thermal_regions=thermal_regions, error=error)

     nregions = thermal_regions%nregions
     nfree_done = 0
     ekin_total_new = 0.0_dp
     ekin_old = md_ener%ekin
     vcom_total = 0.0_dp
     ALLOCATE(temp_new(0:nregions),temp_old(0:nregions), STAT=istat)
     temp_new = 0.0_dp
     temp_old = 0.0_dp
     !loop regions
     DO ireg = 1,nregions
       NULLIFY(t_region)
       t_region => thermal_regions%thermal_region(ireg)
       nfree = 3*t_region%npart
       ekin = compute_ekin(part,ireg)
       IF(nfree > 0) t_region%temperature = 2.0_dp*ekin/REAL(nfree,KIND=dp)*kelvin
       temp_old(ireg) =  t_region%temperature
       IF(t_region%temp_tol > 0.0_dp .AND. &
         ABS(t_region%temp_expected - t_region%temperature/kelvin) > t_region%temp_tol) THEN
           fscale =  SQRT((t_region%temp_expected/t_region%temperature)*kelvin)
           CALL scale_velocity_low(subsys,fscale,ireg,ekin,vcom,error)
           t_region%temperature = 2.0_dp*ekin/REAL(nfree,KIND=dp)*kelvin
           temp_new(ireg) = t_region%temperature
       END IF
       nfree_done = nfree_done + nfree
       ekin_total_new = ekin_total_new + ekin
     END DO
     nfree = simpar%nfree - nfree_done
     ekin = compute_ekin(part,ireg=0)
     IF(nfree>0) thermal_regions%temp_reg0 = 2.0_dp*ekin/REAL(nfree,KIND=dp)*kelvin
     temp_old(0) =  thermal_regions%temp_reg0
     IF(simpar%temp_tol > 0.0_dp .AND. nfree>0) THEN
        IF (ABS(simpar%temp_ext - thermal_regions%temp_reg0/kelvin) > simpar%temp_tol) THEN
           fscale =  SQRT((simpar%temp_ext/thermal_regions%temp_reg0)*kelvin)
           CALL scale_velocity_low(subsys,fscale,0,ekin,vcom,error)
           thermal_regions%temp_reg0 = 2.0_dp*ekin/REAL(nfree,KIND=dp)*kelvin
           temp_new(0) = thermal_regions%temp_reg0
        END IF
     END IF
     ekin_total_new = ekin_total_new + ekin

     md_ener%ekin = ekin_total_new
     IF(md_ener%nfree /=0) THEN
         md_ener%temp_part = 2.0_dp*md_ener%ekin/REAL(md_ener%nfree,KIND=dp)*kelvin
     END IF
     md_ener%constant = md_ener%constant - ekin_old + md_ener%ekin

     DO ireg = 0,nregions
       IF (iw>0 .AND. temp_new(ireg)>0.0_dp) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,I5, A,F10.2,A,F10.2,A)")"Temperature region ", ireg, &
             " rescaled from:",temp_old(ireg)," K  to ",temp_new(ireg)," K"
       END IF
     END DO
     DEALLOCATE(temp_new,temp_old,STAT=istat)

  END SUBROUTINE scale_velocity_region

! *****************************************************************************
!> \brief Scale velocities  for a specific region
!> \param subsys ...
!> \param fscale ...
!> \param ireg ...
!> \param ekin ...
!> \param vcom ...
!> \param error ...
!> \par author MI
! *****************************************************************************
  SUBROUTINE scale_velocity_low(subsys,fscale,ireg,ekin,vcom,error)

    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(KIND=dp), INTENT(IN)                :: fscale
    INTEGER, INTENT(IN)                      :: ireg
    REAL(KIND=dp), INTENT(OUT)               :: ekin, vcom(3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, iatom, ikind, my_ireg, &
                                                natom, shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: imass, mass, tmass, v2
    REAL(KIND=dp), DIMENSION(3)              :: tmp, v, vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(atomic_kinds,particles,shell_particles,core_particles,shell,atom_list)

    my_ireg = ireg
    ekin = 0.0_dp
    tmass = 0.0_dp
    vcom = 0.0_dp

    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles,&
          shell_particles=shell_particles, core_particles=core_particles, error=error)

    DO ikind=1,atomic_kinds%n_els
       atomic_kind => atomic_kinds%els(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list, mass=mass,&
            natom=natom, shell_active=is_shell, shell=shell)
       IF(is_shell) THEN
          imass = 1.0_dp / mass
          v2 = 0.0_dp
          v  = 0.0_dp
          DO iatom=1,natom
             atom = atom_list(iatom)
             !check region
             IF( particles%els(atom)%t_region_index/=my_ireg) CYCLE

             particles%els(atom)%v(:) = fscale*particles%els(atom)%v
             shell_index = particles%els(atom)%shell_index
             vs = shell_particles%els(shell_index)%v
             vc = core_particles %els(shell_index)%v
             tmp(1) = imass*(vs(1)-vc(1))
             tmp(2) = imass*(vs(2)-vc(2))
             tmp(3) = imass*(vs(3)-vc(3))

             shell_particles%els(shell_index)%v(1) = particles%els(atom)%v(1) + tmp(1)*shell%mass_core
             shell_particles%els(shell_index)%v(2) = particles%els(atom)%v(2) + tmp(2)*shell%mass_core
             shell_particles%els(shell_index)%v(3) = particles%els(atom)%v(3) + tmp(3)*shell%mass_core

             core_particles%els(shell_index)%v(1) =  particles%els(atom)%v(1) - tmp(1)*shell%mass_shell
             core_particles%els(shell_index)%v(2) =  particles%els(atom)%v(2) - tmp(2)*shell%mass_shell
             core_particles%els(shell_index)%v(3) =  particles%els(atom)%v(3) - tmp(3)*shell%mass_shell

             ! kinetic energy and velocity of COM
             v2 = v2 + SUM(particles%els(atom)%v**2)
             v(1) = v(1) + particles%els(atom)%v(1)
             v(2) = v(2) + particles%els(atom)%v(2)
             v(3) = v(3) + particles%els(atom)%v(3)
             tmass = tmass + mass
          END DO
       ELSE
          v2 = 0.0_dp
          v = 0.0_dp
          DO iatom=1,natom
             atom = atom_list(iatom)
             !check region
             IF( particles%els(atom)%t_region_index/=my_ireg) CYCLE

             particles%els(atom)%v(:) = fscale*particles%els(atom)%v
            ! kinetic energy and velocity of COM
             v2 = v2 + SUM(particles%els(atom)%v**2)
             v(1) = v(1) + particles%els(atom)%v(1)
             v(2) = v(2) + particles%els(atom)%v(2)
             v(3) = v(3) + particles%els(atom)%v(3)
             tmass = tmass + mass
          END DO
       END IF
       ekin  = ekin + 0.5_dp*mass*v2
       vcom(1)    = vcom(1) + mass*v(1)
       vcom(2)    = vcom(2) + mass*v(2)
       vcom(3)    = vcom(3) + mass*v(3)

     END DO
     vcom = vcom / tmass

  END SUBROUTINE scale_velocity_low

! *****************************************************************************
!> \brief Scale internal motion of CORE-SHELL model to the correct temperature
!> \param subsys ...
!> \param md_ener ...
!> \param temp_expected ...
!> \param temp_tol ...
!> \param iw ...
!> \param error ...
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE scale_velocity_internal(subsys, md_ener,  temp_expected, temp_tol, iw, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity_internal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, iatom, ikind, natom, &
                                                shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: ekin_shell_old, fac_mass, &
                                                mass, scale, temp_shell_old, &
                                                v2
    REAL(KIND=dp), DIMENSION(3)              :: tmp, v, vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(atom_list,atomic_kinds,atomic_kind,core_particles,particles,shell_particles,shell)
    IF (ABS(temp_expected - md_ener%temp_shell/kelvin) > temp_tol) THEN
       scale = 0.0_dp
       IF (md_ener%temp_shell>EPSILON(0.0_dp)) scale = SQRT((temp_expected/md_ener%temp_shell)*kelvin)
       ekin_shell_old = md_ener%ekin_shell
       temp_shell_old = md_ener%temp_shell
       md_ener%ekin_shell = 0.0_dp
       md_ener%temp_shell = 0.0_dp

       CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, particles=particles, shell_particles=shell_particles,&
          core_particles=core_particles, error=error)

       DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom_list, mass=mass, natom=natom,&
             shell_active=is_shell, shell=shell)
          IF(is_shell) THEN
             fac_mass = 1.0_dp/mass
             v2 = 0.0_dp
             DO iatom= 1, natom
                atom = atom_list(iatom)
                shell_index = particles%els(atom)%shell_index
                vs  = shell_particles%els(shell_index)%v
                vc  = core_particles%els(shell_index)%v
                v   = particles%els(atom)%v
                tmp(1) = fac_mass*(vc(1)-vs(1))
                tmp(2) = fac_mass*(vc(2)-vs(2))
                tmp(3) = fac_mass*(vc(3)-vs(3))

                shell_particles%els(shell_index)%v(1) = v(1) - shell%mass_core*scale*tmp(1)
                shell_particles%els(shell_index)%v(2) = v(2) - shell%mass_core*scale*tmp(2)
                shell_particles%els(shell_index)%v(3) = v(3) - shell%mass_core*scale*tmp(3)

                core_particles%els(shell_index)%v(1) =  v(1) + shell%mass_shell*scale*tmp(1)
                core_particles%els(shell_index)%v(2) =  v(2) + shell%mass_shell*scale*tmp(2)
                core_particles%els(shell_index)%v(3) =  v(3) + shell%mass_shell*scale*tmp(3)

                vs  = shell_particles%els(shell_index)%v
                vc  = core_particles%els(shell_index)%v
                tmp(1) = vc(1) - vs(1)
                tmp(2) = vc(2) - vs(2)
                tmp(3) = vc(3) - vs(3)
                v2 = v2 + SUM(tmp**2)
             END DO
             md_ener%ekin_shell = md_ener%ekin_shell + 0.5_dp*shell%mass_core*shell%mass_shell*fac_mass*v2
          END IF
       END DO
       IF(md_ener%nfree_shell>0)THEN
           md_ener%temp_shell = 2.0_dp*md_ener%ekin_shell/REAL(md_ener%nfree_shell,KIND=dp)*kelvin
       END IF
       md_ener%constant = md_ener%constant - ekin_shell_old + md_ener%ekin_shell

       IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
               "Temperature  shell internal motion scaled to requested temperature:",&
               temp_shell_old," K  ->",md_ener%temp_shell," K"
       END IF
    ENDIF
  END SUBROUTINE scale_velocity_internal

! *****************************************************************************
!> \brief Scale barostat velocities to get the desired temperature
!> \param md_env ...
!> \param md_ener ...
!> \param temp_expected ...
!> \param temp_tol ...
!> \param iw ...
!> \param error ...
!> \par History
!>     MI 02.2008
! *****************************************************************************
  SUBROUTINE scale_velocity_baro(md_env, md_ener, temp_expected, temp_tol, iw, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), POINTER              :: md_ener
    REAL(KIND=dp), INTENT(IN)                :: temp_expected, temp_tol
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_velocity_baro', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nfree
    REAL(KIND=dp)                            :: ekin_old, scale, temp_old
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(simpar_type), POINTER               :: simpar

    NULLIFY( npt, simpar)
    CALL get_md_env ( md_env, simpar = simpar, npt = npt,  error=error)
    IF (ABS(temp_expected - md_ener%temp_baro/kelvin) > temp_tol) THEN
       scale = 0.0_dp
       IF (md_ener%temp_baro>0.0_dp) scale = SQRT((temp_expected/md_ener%temp_baro)*kelvin)
       ekin_old = md_ener%baro_kin
       temp_old = md_ener%temp_baro
       md_ener%baro_kin       = 0.0_dp
       md_ener%temp_baro  = 0.0_dp
      IF ( simpar%ensemble==npt_i_ensemble .OR. simpar%ensemble==npe_i_ensemble) THEN
         npt ( 1, 1 )%v = npt ( 1, 1 )%v*scale
         md_ener%baro_kin = 0.5_dp * npt ( 1, 1 )%v**2 * npt ( 1, 1 )%mass
       ELSE IF (simpar%ensemble==npt_f_ensemble .OR. simpar%ensemble==npe_f_ensemble) THEN
         md_ener%baro_kin = 0.0_dp
         DO i = 1, 3
            DO j = 1, 3
             npt(i,j)%v = npt(i,j)%v*scale
             md_ener%baro_kin = md_ener%baro_kin + 0.5_dp * npt(i,j)%v**2  * npt ( i, j )%mass
            END DO
         END DO
      END IF

      nfree = SIZE ( npt, 1 ) * SIZE ( npt, 2 )
      md_ener%temp_baro = 2.0_dp * md_ener%baro_kin / REAL(nfree,dp)*kelvin
      IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
               "Temperature  of barostat motion scaled to requested temperature:",&
               temp_old," K  ->",md_ener%temp_baro," K"
      END IF
   END IF
  END SUBROUTINE scale_velocity_baro

! *****************************************************************************
!> \brief Perform all temperature manipulations during a QS MD run.
!> \param simpar ...
!> \param md_env ...
!> \param md_ener ...
!> \param force_env ...
!> \param logger ...
!> \param error ...
!> \par History
!>     Creation (15.09.2003,MK)
!>     adapted to force_env (05.10.2003,fawzi)
!>     Cleaned (09.2007) Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE temperature_control(simpar, md_env, md_ener,force_env,logger, error)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'temperature_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,handle)
    NULLIFY(subsys, para_env)
    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(md_ener),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)
    iw = cp_print_key_unit_nr(logger,force_env%root_section,"MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
         extension=".mdLog",error=error)

    ! Control the particle motion
    IF(simpar%do_thermal_region) THEN
       CALL scale_velocity_region(md_env, subsys, md_ener, simpar, iw, error)
    ELSE
      IF (simpar%temp_tol > 0.0_dp ) THEN
         CALL scale_velocity(subsys, md_ener, simpar%temp_ext, simpar%temp_tol, iw, error)
      END IF
    END IF
    ! Control the internal core-shell motion
    IF(simpar%temp_sh_tol > 0.0_dp) THEN
       CALL scale_velocity_internal(subsys, md_ener, simpar%temp_sh_ext, simpar%temp_sh_tol, iw, error)
    END IF
    ! Control cell motion
    SELECT CASE (simpar%ensemble)
    CASE( nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, &
            npt_f_ensemble, npt_i_ensemble, npe_f_ensemble, npe_i_ensemble)
       IF(simpar%temp_baro_tol > 0.0_dp) THEN
          CALL scale_velocity_baro(md_env, md_ener, simpar%temp_baro_ext, simpar%temp_baro_tol, iw, error)
       END IF
    END SELECT

    CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
         "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    CALL timestop(handle)
  END SUBROUTINE temperature_control

! *****************************************************************************
!> \brief Set to 0 the velocity of the COM along MD runs, if required.
!> \param md_ener ...
!> \param force_env ...
!> \param md_section ...
!> \param logger ...
!> \param error ...
!> \par History
!>      Creation (29.04.2007,MI)
!>      Cleaned (09.2007) Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE comvel_control(md_ener,force_env, md_section, logger, error)

    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comvel_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: comvel_tol, temp_old, vel_com
    REAL(KIND=dp), DIMENSION(3)              :: vcom_old
    TYPE(cp_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,handle)
    NULLIFY(subsys)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,error=error)

    ! Print COMVEL and COM Position
    iw = cp_print_key_unit_nr(logger,md_section,"PRINT%CENTER_OF_MASS",extension=".mdLog",error=error)
    IF (iw>0) THEN
       WRITE (UNIT=iw,FMT="(/,T2,A,(T58,A3,F20.10))")&
            "Centre of mass motion (COM):","x =",md_ener%vcom(1),"y =",md_ener%vcom(2),"z =",md_ener%vcom(3)
    END IF
    CALL cp_print_key_finished_output(iw,logger,md_section,"PRINT%CENTER_OF_MASS", error=error)

    ! If requested rescale COMVEL
    CALL section_vals_val_get(md_section,"COMVEL_TOL",explicit=explicit,error=error)
    IF ( explicit ) THEN
       CALL section_vals_val_get(md_section,"COMVEL_TOL",r_val=comvel_tol,error=error)
       iw = cp_print_key_unit_nr(logger,md_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".mdLog",error=error)
       vel_com = SQRT(md_ener%vcom(1)**2+md_ener%vcom(2)**2+md_ener%vcom(3)**2)

       ! Subtract the velocity of the COM, if requested
       IF (vel_com > comvel_tol) THEN
          temp_old = md_ener%temp_part/kelvin
          vcom_old = md_ener%vcom
          CALL reset_vcom( subsys, md_ener, vsubtract=vcom_old, error=error)
          CALL scale_velocity(subsys, md_ener, temp_old, 0.0_dp, iw, error)
          IF (iw>0) THEN
             WRITE (UNIT=iw,FMT="(T2,'MD| ',A,3F16.10,A)") &
                  "Old VCOM = ",vcom_old(1:3)," a.u.",&
                  "New VCOM = ",md_ener%vcom(1:3)," a.u"
          END IF
       END IF
       CALL cp_print_key_finished_output(iw,logger,md_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

    CALL timestop(handle)
  END SUBROUTINE comvel_control

! *****************************************************************************
!> \brief Set to 0 the angular velocity along MD runs, if required.
!> \param md_ener ...
!> \param force_env ...
!> \param md_section ...
!> \param logger ...
!> \param error ...
!> \par History
!>      Creation (10.2009) Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE angvel_control(md_ener, force_env, md_section, logger, error)

    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'angvel_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ifixd, &
                                                imolecule_kind, iw, natoms, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: is_fixed
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: angvel_tol, rcom(3), &
                                                temp_old, vang(3), vang_new(3)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_list_type), POINTER        :: particles

    CALL timeset(routineN,handle)
    ! If requested rescale ANGVEL
    CALL section_vals_val_get(md_section,"ANGVEL_TOL",explicit=explicit,error=error)
    IF ( explicit ) THEN
       NULLIFY(subsys, cell)
       CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
       CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)

       IF (SUM(cell%perd(1:3)) == 0) THEN
          CALL section_vals_val_get(md_section,"ANGVEL_TOL",r_val=angvel_tol,error=error)
          iw = cp_print_key_unit_nr(logger,md_section,"PRINT%PROGRAM_RUN_INFO",&
               extension=".mdLog",error=error)

          CALL cp_subsys_get(subsys,molecule_kinds_new=molecule_kinds,&
               particles=particles, error=error)

          natoms = SIZE(particles%els)
          ! Build a list of all fixed atoms (if any)
          ALLOCATE (is_fixed(natoms),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          is_fixed = use_perd_none
          molecule_kind_set => molecule_kinds%els
          DO imolecule_kind=1,molecule_kinds%n_els
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind,fixd_list=fixd_list)
             IF (ASSOCIATED(fixd_list)) THEN
                DO ifixd=1,SIZE(fixd_list)
                   IF (.NOT.fixd_list(ifixd)%restraint%active) &
                        is_fixed(fixd_list(ifixd)%fixd) = fixd_list(ifixd)%itype
                END DO
             END IF
          END DO

          ! If requested and the system is not periodic, subtract the angular velocity
          CALL compute_rcom(particles%els,is_fixed,rcom)
          CALL compute_vang(particles%els,is_fixed,rcom,vang,error)
          IF (SQRT(DOT_PRODUCT(vang,vang))>angvel_tol) THEN
             CALL subtract_vang(particles%els,is_fixed,rcom,vang)

             ! Rescale velocities after removal
             temp_old = md_ener%temp_part/kelvin
             CALL scale_velocity(subsys, md_ener, temp_old, 0.0_dp, iw, error)
             CALL compute_vang(particles%els,is_fixed,rcom,vang_new,error)
             IF (iw>0) THEN
                WRITE (UNIT=iw,FMT="(T2,'MD| ',A,3F16.10,A)") &
                     "Old VANG = ",vang(1:3)," a.u.",&
                     "New VANG = ",vang_new(1:3)," a.u"
             END IF
          END IF

          DEALLOCATE (is_fixed,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          CALL cp_print_key_finished_output(iw,logger,md_section,&
               "PRINT%PROGRAM_RUN_INFO", error=error)
       END IF
    END IF

    CALL timestop(handle)
  END SUBROUTINE angvel_control

! *****************************************************************************
!> \brief Initialize Velocities for MD runs
!> \param force_env ...
!> \param simpar ...
!> \param globenv ...
!> \param md_env ...
!> \param md_section ...
!> \param constraint_section ...
!> \param write_binary_restart_file ...
!> \param error ...
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE setup_velocities(force_env, simpar, globenv, md_env, md_section, &
                              constraint_section, write_binary_restart_file, &
                              error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: md_section, constraint_section
    LOGICAL, INTENT(IN)                      :: write_binary_restart_file
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nconstraint, &
                                                nconstraint_fixd
    LOGICAL                                  :: apply_cns0, failure, &
                                                shell_adiabatic, shell_present
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                print_section, subsys_section

    failure = .FALSE.
    
    CALL timeset(routineN,handle)

    NULLIFY (atomic_kinds,cell,para_env,subsys,molecule_kinds,core_particles,particles)
    NULLIFY (shell_particles,core_particle_set,particle_set,shell_particle_set)
    NULLIFY (force_env_section,print_section,subsys_section)

    print_section => section_vals_get_subs_vals(md_section,"PRINT",error=error)
    apply_cns0 = .FALSE.
    IF (simpar%constraint) THEN
       CALL section_vals_val_get(constraint_section,"CONSTRAINT_INIT",l_val=apply_cns0,error=error)
    END IF
    ! Always initialize velocities and possibly restart them
    CALL force_env_get(force_env, subsys=subsys, cell=cell, para_env=para_env,&
         force_env_section=force_env_section, error=error )
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)

    CALL cp_subsys_get(subsys,&
                       atomic_kinds=atomic_kinds,&
                       core_particles=core_particles,&
                       molecule_kinds_new=molecule_kinds,&
                       particles=particles,&
                       shell_particles=shell_particles,&
                       error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
                             shell_present=shell_present,&
                             shell_adiabatic=shell_adiabatic)

    NULLIFY (core_particle_set)
    NULLIFY (particle_set)
    NULLIFY (shell_particle_set)
    particle_set => particles%els

    IF (shell_present.AND.shell_adiabatic) THEN
       ! Constraints are not yet implemented for core-shell models generally
       CALL get_molecule_kind_set(molecule_kind_set=molecule_kinds%els,&
                                  nconstraint=nconstraint,&
                                  nconstraint_fixd=nconstraint_fixd)
       CALL cp_assert((nconstraint - nconstraint_fixd == 0),&
                      cp_fatal_level,cp_assertion_failed,routineP,&
                      "Only the fixed atom constraint is implemented for core-shell models",&
                      only_ionode=.TRUE.)
!MK    CPPostcondition(.NOT.simpar%constraint,cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(shell_particles),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(core_particles),cp_failure_level,routineP,error,failure)
       shell_particle_set => shell_particles%els
       core_particle_set  =>  core_particles%els
    END IF

    CALL initialize_velocities(simpar,particle_set, molecule_kinds=molecule_kinds,&
         force_env=force_env, globenv=globenv, md_env=md_env,&
         label="Velocities initialization",print_section=print_section, &
         subsys_section=subsys_section, shell_present=(shell_present.AND.shell_adiabatic), &
         shell_part=shell_particle_set, core_part=core_particle_set, force_rescaling=.FALSE., &
         para_env=para_env,error=error)

    ! Apply constraints if required and rescale velocities..
    IF (simpar%ensemble /= reftraj_ensemble) THEN
       IF (apply_cns0) THEN
          CALL force_env_calc_energy_force ( force_env, calc_force=.TRUE.,error=error)
          CALL force_env_shake(force_env,shake_tol=simpar%shake_tol,&
               log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
               dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
          CALL force_env_rattle(force_env,shake_tol=simpar%shake_tol,&
               log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
               dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)
          IF (simpar%do_respa)THEN
             CALL force_env_calc_energy_force (force_env%sub_force_env(1)%force_env,&
                  calc_force=.TRUE.,error=error)
             CALL force_env_shake(force_env%sub_force_env(1)%force_env,&
                  shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                  lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
             CALL force_env_rattle(force_env%sub_force_env(1)%force_env,&
                  shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                  lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)
          END IF
          ! Reinitialize velocities rescaling properly after rattle
          subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
          CALL update_subsys(subsys_section,force_env,.FALSE.,write_binary_restart_file,error)

          CALL initialize_velocities(simpar,particle_set, molecule_kinds=molecule_kinds,&
               force_env=force_env, globenv=globenv, md_env=md_env,&
               label="Re-Initializing velocities after applying constraints",print_section=print_section, &
               subsys_section=subsys_section, shell_present=(shell_present.AND.shell_adiabatic), &
               shell_part=shell_particle_set, core_part=core_particle_set, force_rescaling=.TRUE., &
               para_env=para_env,error=error)
       END IF
    END IF

    ! Perform setup for a cascade run
    CALL initialize_cascade(simpar,particle_set,molecule_kinds,md_section,&
                            error=error)

    CALL timestop(handle)

  END SUBROUTINE setup_velocities

! *****************************************************************************
!> \brief   Perform the initialization for a cascade run
!> \param simpar ...
!> \param particle_set ...
!> \param molecule_kinds ...
!> \param md_section ...
!> \param error ...
!> \date    05.02.2012
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE initialize_cascade(simpar,particle_set,molecule_kinds,md_section,&
                                error)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_cascade', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2*default_string_length)   :: line
    INTEGER                                  :: handle, iatom, ifixd, &
                                                imolecule_kind, iparticle, &
                                                iw, natom, nparticle, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_index, is_fixed
    LOGICAL                                  :: failure, init_cascade, is_ok, &
                                                no_read_error
    REAL(KIND=dp)                            :: ecom, ekin, energy, norm, &
                                                temp, temperature
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: matom, weight
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vatom
    REAL(KIND=dp), DIMENSION(3)              :: vcom
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: atom_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(section_vals_type), POINTER         :: atom_list_section, &
                                                cascade_section, print_section
    TYPE(val_type), POINTER                  :: val

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atom_list)
    NULLIFY (atom_list_section)
    NULLIFY (atomic_kind)
    NULLIFY (cascade_section)
    NULLIFY (fixd_list)
    NULLIFY (molecule_kind)
    NULLIFY (molecule_kind_set)
    NULLIFY (logger)
    NULLIFY (val)

    logger => cp_error_get_logger(error)
    print_section => section_vals_get_subs_vals(md_section,"PRINT",error=error)
    iw = cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",extension=".log",error=error)

    cascade_section => section_vals_get_subs_vals(md_section,"CASCADE",error=error)
    CALL section_vals_val_get(cascade_section,"_SECTION_PARAMETERS_",l_val=init_cascade,error=error)

    nparticle = SIZE(particle_set)

    IF (init_cascade) THEN

       CALL section_vals_val_get(cascade_section,"ENERGY",r_val=energy,error=error)
       CALL cp_assert((energy >= 0.0_dp),cp_fatal_level,cp_assertion_failed,routineP,&
                      "Error occurred reading &CASCADE section: Negative energy found",&
                      only_ionode=.TRUE.)

       IF (iw > 0) THEN
          ekin = cp_unit_from_cp2k(energy,"keV",error=error)
          WRITE (UNIT=iw,FMT="(T2,A,T61,F20.6)")&
           "CASCADE| Energy [keV]",ekin
          WRITE (UNIT=iw,FMT="(T2,A)")&
           "CASCADE|"
       END IF

       ! Read the atomic velocities given in the input file
       atom_list_section => section_vals_get_subs_vals(cascade_section,"ATOM_LIST",error=error)
       CALL section_vals_val_get(atom_list_section,"_DEFAULT_KEYWORD_",n_rep_val=natom,error=error)
       CALL section_vals_list_get(atom_list_section,"_DEFAULT_KEYWORD_",list=atom_list,error=error)
       CALL cp_assert((natom > 0),cp_fatal_level,cp_assertion_failed,routineP,&
                      "Error occurred reading &CASCADE section: No atom list found",&
                      only_ionode=.TRUE.)

       IF (iw > 0) THEN
          WRITE (UNIT=iw,FMT="(T2,A,T11,A,3(11X,A),9X,A)")&
           "CASCADE| ","Atom index","v(x)","v(y)","v(z)","weight"
       END IF

       ALLOCATE (atom_index(natom),STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
       ALLOCATE (matom(natom),STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
       ALLOCATE (vatom(3,natom),STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
       ALLOCATE (weight(natom),STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)

       DO iatom=1,natom
          is_ok = cp_sll_val_next(atom_list,val,error=error)
          CALL val_get(val,c_val=line,error=error)
          ! Read atomic index, velocity vector, and weight
          no_read_error = .FALSE.
          READ (UNIT=line,FMT=*,ERR=999) atom_index(iatom),vatom(1:3,iatom),weight(iatom)
          no_read_error = .TRUE.
999       CALL cp_assert(no_read_error,cp_fatal_level,cp_assertion_failed,routineP,&
                         "Error occurred reading &CASCADE section. Last line read <"//&
                         TRIM(line)//">",only_ionode=.TRUE.)
          CALL cp_assert((atom_index(iatom) > 0).AND.((atom_index(iatom) <= nparticle)),&
                         cp_fatal_level,cp_assertion_failed,routineP,&
                         "Error occurred reading &CASCADE section: Invalid atom index found",&
                         only_ionode=.TRUE.)
          CALL cp_assert((weight(iatom) >= 0.0_dp),cp_fatal_level,cp_assertion_failed,routineP,&
                         "Error occurred reading &CASCADE section: Negative weight found",&
                         only_ionode=.TRUE.)
          IF (iw > 0) THEN
             WRITE (UNIT=iw,FMT="(T2,A,I10,4(1X,F14.6))")&
              "CASCADE| ",atom_index(iatom),vatom(1:3,iatom),weight(iatom)
          END IF
       END DO

       ! Normalise velocities and weights
       norm = 0.0_dp
       DO iatom=1,natom
          iparticle = atom_index(iatom)
          CALL cp_assert((particle_set(iparticle)%shell_index == 0),&
                         cp_warning_level,cp_assertion_failed,routineP,&
                         "Warning: The primary knock-on atom is a core-shell atom",&
                         only_ionode=.TRUE.)
          atomic_kind => particle_set(iparticle)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=matom(iatom))
          norm = norm + matom(iatom)*weight(iatom)
       END DO
       weight(:) = matom(:)*weight(:)*energy/norm
       DO iatom=1,natom
          norm = SQRT(DOT_PRODUCT(vatom(1:3,iatom),vatom(1:3,iatom)))
          vatom(1:3,iatom) = vatom(1:3,iatom)/norm
       END DO

       IF (iw > 0) THEN
          WRITE (UNIT=iw,FMT="(T2,A)")&
           "CASCADE|",&
           "CASCADE| Normalised velocities and additional kinetic energy [keV]",&
           "CASCADE|"
          WRITE (UNIT=iw,FMT="(T2,A,T11,A,3(11X,A),9X,A)")&
           "CASCADE| ","Atom index","v(x)","v(y)","v(z)","E(kin)"
          DO iatom=1,natom
             ekin = cp_unit_from_cp2k(weight(iatom),"keV",error=error)
             WRITE (UNIT=iw,FMT="(T2,A,I10,4(1X,F14.6))")&
              "CASCADE| ",atom_index(iatom),vatom(1:3,iatom),ekin
          END DO
       END IF

       ! Apply velocity modifications
       DO iatom=1,natom
          iparticle = atom_index(iatom)
          particle_set(iparticle)%v(:) = particle_set(iparticle)%v(:) +&
                                         SQRT(2.0_dp*weight(iatom)/matom(iatom))*vatom(1:3,iatom)
       END DO

       DEALLOCATE (atom_index,STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
       DEALLOCATE (matom,STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
       DEALLOCATE (vatom,STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)
       DEALLOCATE (weight,STAT=stat)
       CPPostcondition((stat == 0),cp_failure_level,routineP,error,failure)

       IF (iw > 0) THEN
          ! Build a list of all fixed atoms (if any)
          ALLOCATE (is_fixed(nparticle),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          is_fixed = use_perd_none
          molecule_kind_set => molecule_kinds%els
          DO imolecule_kind=1,molecule_kinds%n_els
             molecule_kind => molecule_kind_set(imolecule_kind)
             CALL get_molecule_kind(molecule_kind=molecule_kind,fixd_list=fixd_list)
             IF (ASSOCIATED(fixd_list)) THEN
                DO ifixd=1,SIZE(fixd_list)
                   IF (.NOT.fixd_list(ifixd)%restraint%active) is_fixed(fixd_list(ifixd)%fixd) = fixd_list(ifixd)%itype
                END DO
             END IF
          END DO
          ! Compute vcom, ecom and ekin for printout
          CALL compute_vcom(particle_set,is_fixed,vcom,ecom)
          ekin = compute_ekin(particle_set) - ecom
          IF (simpar%nfree == 0) THEN
             CPPostcondition((ekin == 0.0_dp),cp_failure_level,routineP,error,failure)
             temp = 0.0_dp
          ELSE
             temp = 2.0_dp*ekin/REAL(simpar%nfree,KIND=dp)
          END IF
          temperature = cp_unit_from_cp2k(temp,"K",error=error)
          WRITE (UNIT=iw,FMT="(T2,A)")&
           "CASCADE|"
          WRITE (UNIT=iw,FMT="(T2,A,T61,F18.2,A2)")&
           "CASCADE| Temperature after cascade initialization",temperature," K"
          WRITE (UNIT=iw,FMT="(T2,A,T30,3(1X,ES16.8),/)")&
           "CASCADE| COM velocity: ",vcom(1:3)
!MK          ! compute and log rcom and vang if not periodic
!MK          CALL force_env_get(force_env,cell=cell,error=error)
!MK          IF (SUM(cell%perd(1:3)) == 0) THEN
!MK             CALL compute_rcom(particle_set,is_fixed,rcom)
!MK             CALL compute_vang(particle_set,is_fixed,rcom,vang,error)
!MK             WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' COM position:',rcom(1:3)
!MK             WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' Angular velocity:',vang(1:3)
!MK          END IF
          DEALLOCATE (is_fixed,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

    END IF

    CALL cp_print_key_finished_output(iw,logger,print_section,"PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE initialize_cascade

END MODULE md_vel_utils
