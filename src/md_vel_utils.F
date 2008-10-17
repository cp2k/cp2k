!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Collection of utilities for setting-up and handle velocities in MD
!>         runs
!> \author CJM
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! *****************************************************************************
MODULE md_vel_utils
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_rattle,&
                                             force_env_shake
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       npe_f_ensemble, npe_i_ensemble, nph_uniaxial_damped_ensemble, &
       nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, &
       reftraj_ensemble, use_perd_none, use_perd_x, use_perd_xy, &
       use_perd_xyz, use_perd_xz, use_perd_y, use_perd_yz, use_perd_z
  USE input_cp2k_restarts,             ONLY: update_subsys
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE md_ener_types,                   ONLY: md_ener_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: kelvin
  USE shell_potential_types,           ONLY: shell_kind_type
  USE simpar_types,                    ONLY: simpar_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_vel_utils'

 
  PUBLIC :: temperature_control,&
            comvel_control,&
            setup_velocities

CONTAINS

! *****************************************************************************
!> \brief compute center of mass position
!>      *** is only used by initialize_velocities below ***
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
       SELECT CASE(is_fixed(i))
       CASE(use_perd_x,use_perd_y,use_perd_z,use_perd_xy,use_perd_xz,use_perd_yz,use_perd_none)
          vcom(1) = vcom(1) + part(i)%v(1) * mass
          vcom(2) = vcom(2) + part(i)%v(2) * mass
          vcom(3) = vcom(3) + part(i)%v(3) * mass
          denom = denom + mass
       END SELECT
    END DO
    vcom = vcom/denom
    IF (PRESENT(ecom)) THEN
       ecom = 0.5_dp*denom*SUM(vcom*vcom)
    END IF

  END SUBROUTINE compute_vcom

! *****************************************************************************
!> \brief Copy atom velocities into core and shell velocities
!>      *** is only used by initialize_velocities below ***
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
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  FUNCTION compute_ekin(part) RESULT(ekin)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    REAL(KIND=dp)                            :: ekin

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_ekin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    ekin = 0.0_dp
    DO i = 1, SIZE(part)
       atomic_kind => part(i)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
       ekin = ekin + 0.5_dp * mass * SUM(part(i)%v(:) * part(i)%v(:))
    END DO

  END FUNCTION compute_ekin

! *****************************************************************************
!> \brief Rescale the velocity to mimic the given external kinetic temperature.
!>      Optionally also rescale vcom.
!>      *** is only used by initialize_velocities below ***
!> \par History
!>      2007-11-6: created
!> \author Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE rescale_vel(part,simpar,ekin,vcom)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(simpar_type), POINTER               :: simpar
    REAL(KIND=dp), INTENT(INOUT)             :: ekin
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT), OPTIONAL                :: vcom

    CHARACTER(len=*), PARAMETER :: routineN = 'rescale_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: factor

    IF (simpar%nfree/=0) THEN
       factor = simpar%temp_ext / ( 2.0_dp * ekin ) * REAL(simpar%nfree,KIND=dp)
    ELSE
       factor = 0.0_dp
    ENDIF
    ! Note:
    ! this rescaling is still wrong, it should take the masses into account
    ! rescaling is generally not correct, so needs fixing
    ekin = ekin * factor
    factor = SQRT(factor)
    DO i = 1, SIZE(part)
       part(i)%v(:) = factor*part(i)%v(:)
    END DO
    IF (PRESENT(vcom)) THEN
       vcom = factor*vcom
    END IF

  END SUBROUTINE rescale_vel

! *****************************************************************************
!> \brief subtract center of mass velocity
!>      *** is only used by initialize_velocities below ***
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
!> \par History
!>      - is_fixed removed from particle_type
!>      - 2007-11-07: Cleanup (TV)
!>      - 2007-11-09: Added angvel_zero feature
!> \author CJM,MK,Toon Verstraelen <Toon.Verstraelen@gmail.com>
! *****************************************************************************
  SUBROUTINE initialize_velocities(simpar, part, force_env, para_env, globenv, &
       molecule_kinds, print_section, subsys_section, shell_present, shell_part, &
       core_part, force_rescaling, error)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: part
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(section_vals_type), POINTER         :: print_section, subsys_section
    LOGICAL, INTENT(IN)                      :: shell_present
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_part, core_part
    LOGICAL, INTENT(IN)                      :: force_rescaling
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ifixd, &
                                                imolecule_kind, iw, natoms, &
                                                nshell, shell_index, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: is_fixed
    LOGICAL                                  :: atomvel_explicit, &
                                                corevel_explicit, failure, &
                                                is_ok, shellvel_explicit
    REAL(KIND=dp)                            :: ecom, ekin, fac_massc, &
                                                fac_masss, mass, mass_tot, &
                                                temp
    REAL(KIND=dp), DIMENSION(3)              :: rcom, vang, vc, vcom, vs
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: atom_list, core_list, &
                                                shell_list
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(section_vals_type), POINTER         :: atomvel_section, &
                                                corevel_section, &
                                                shellvel_section
    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(unit_convert_type), POINTER         :: unit
    TYPE(val_type), POINTER                  :: val

    CALL timeset(routineN,handle)

    ! Initializing parameters
    failure = .FALSE.
    NULLIFY (atomic_kind, fixd_list, atom_list, shell_list, core_list, logger, molecule_kind)
    NULLIFY (molecule_kind_set, shell, unit, val)

    NULLIFY (atomvel_section,shellvel_section, corevel_section)
    atomvel_section  => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
    shellvel_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
    corevel_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
    
    ! Logging
    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",extension=".log",error=error)
    IF (iw>0) THEN
       WRITE ( iw, '( A, A, A )' ) ' **********************', &
            ' begin of velocity initialization ', '***********************'
    END IF

    ! Initializing parameters
    natoms = SIZE(part)
    nshell = 0

    ! Core-Shell Model
    IF (shell_present) THEN
       CPPostcondition(ASSOCIATED(core_part),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(shell_part),cp_failure_level,routineP,error,failure)
       nshell = SIZE(shell_part)
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

    ! Read or initialize the particle velocities
    CALL get_cp2k_units ( unit )
    CALL section_vals_get(atomvel_section,explicit=atomvel_explicit,error=error)
    CALL section_vals_get(shellvel_section,explicit=shellvel_explicit,error=error)
    CALL section_vals_get(corevel_section,explicit=corevel_explicit,error=error)
    CPPostcondition(shellvel_explicit.EQV.corevel_explicit,cp_failure_level,routineP,error,failure)

    IF (atomvel_explicit) THEN
       ! Read the atom velocities if explicitely given in the input file
       CALL section_vals_list_get(atomvel_section,"_DEFAULT_KEYWORD_",list=atom_list,error=error)
       DO i = 1, natoms
         is_ok=cp_sll_val_next(atom_list,val,error=error)
         CALL val_get(val,r_vals=vel,error=error)
         part(i)%v = vel
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
       IF(shell_present) THEN 
          IF (shellvel_explicit) THEN
             ! If the atoms positions are given (?) and core and shell velocities are
             ! present in the input, read the latter.
             CALL section_vals_list_get(shellvel_section,"_DEFAULT_KEYWORD_",list=shell_list,error=error)
             CALL section_vals_list_get(corevel_section,"_DEFAULT_KEYWORD_",list=core_list,error=error)
             DO i = 1 , nshell
                is_ok=cp_sll_val_next(shell_list,val,error=error)
                CALL val_get(val,r_vals=vel,error=error)
                shell_part(i)%v  = vel
                is_ok=cp_sll_val_next(core_list,val,error=error)
                CALL val_get(val,r_vals=vel,error=error)
                core_part(i)%v  = vel
             END DO
          ELSE
             ! Otherwise, just copy atom velocties into shell and core velocities.
             CALL clone_core_shell_vel(part,shell_part,core_part)
          END IF
       END IF
       
       ! compute vcom, ecom and ekin
       CALL compute_vcom(part,is_fixed,vcom,ecom)
       ekin = compute_ekin(part) - ecom

       IF (simpar%nfree /= 0 .AND. force_rescaling) THEN
          CALL rescale_vel(part,simpar,ekin,vcom)
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
    ELSE
       ! Initializing velocities deterministically on all processors, if not given in input
       DO i = 1, natoms
          atomic_kind => part(i)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
          part(i)%v(1) = 0.0_dp
          part(i)%v(2) = 0.0_dp
          part(i)%v(3) = 0.0_dp
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
       END DO
       
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
       ekin = compute_ekin(part)
       CALL rescale_vel(part,simpar,ekin)
       
       ! Initialize the core and the shell velocity. Atom velocities are just
       ! copied so that the initial relative core-shell velocity is zero.
       IF (shell_present) THEN
         CALL clone_core_shell_vel(part, shell_part, core_part)
       ENDIF 
    END IF

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
       WRITE (iw, '( A, T61, F18.2, A2 )' ) ' Initial Temperature ', temp / unit%kelvin, " K"
       WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' COM velocity:', vcom ( 1 ), vcom ( 2 ), vcom ( 3 )

       ! compute and log rcom and vang if not periodic
       CALL force_env_get(force_env, cell=cell, error=error)
       IF (SUM(cell%perd(1:3)) == 0) THEN
          CALL compute_rcom(part,is_fixed,rcom)
          CALL compute_vang(part,is_fixed,rcom,vang,error)
          WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' COM position:', rcom ( 1 ), rcom ( 2 ), rcom ( 3 )
          WRITE (iw, '( A, T21, F20.12 , F20.12 , F20.12 )' ) ' Angular velocity:', vang ( 1 ), vang ( 2 ), vang ( 3 )
       END IF
       WRITE ( iw, '( A, A, A, / )' ) ' ***********************', ' end of velocity initialization ', '************************'
    END IF
    
    DEALLOCATE (is_fixed,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_print_key_finished_output(iw,logger,print_section,"PROGRAM_RUN_INFO", error=error)
    CALL timestop(handle)    
  END SUBROUTINE initialize_velocities

! *****************************************************************************
!> \brief Computes Ekin, VCOM and Temp for particles
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

    INTEGER                                  :: atom, iatom, ikind, natom, &
                                                shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: ekin_old, imass, mass, scale, &
                                                temp_old, v2
    REAL(KIND=dp), DIMENSION(3)              :: tmp, v, vc, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

    IF (ABS(temp_expected - md_ener%temp_part/kelvin) > temp_tol) THEN
       scale = 0.0_dp
       IF (md_ener%temp_part>0.0_dp) scale = SQRT((temp_expected/md_ener%temp_part)*kelvin)
       ekin_old = md_ener%ekin
       temp_old = md_ener%temp_part
       md_ener%ekin       = 0.0_dp
       md_ener%temp_part  = 0.0_dp
       md_ener%vcom       = 0.0_dp
       md_ener%total_mass = 0.0_dp

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
                particles%els(atom)%v(:) = scale*particles%els(atom)%v
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
             END DO
          ELSE
             v2 = 0.0_dp
             v = 0.0_dp
             DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(1) = scale*particles%els(atom)%v(1)
                particles%els(atom)%v(2) = scale*particles%els(atom)%v(2)
                particles%els(atom)%v(3) = scale*particles%els(atom)%v(3)
               ! kinetic energy and velocity of COM
                v2 = v2 + SUM(particles%els(atom)%v**2)
                v(1) = v(1) + particles%els(atom)%v(1)
                v(2) = v(2) + particles%els(atom)%v(2)
                v(3) = v(3) + particles%els(atom)%v(3)     
             END DO 
          END IF
          md_ener%ekin  = md_ener%ekin + 0.5_dp*mass*v2
          md_ener%vcom(1)    = md_ener%vcom(1) + mass*v(1)
          md_ener%vcom(2)    = md_ener%vcom(2) + mass*v(2)
          md_ener%vcom(3)    = md_ener%vcom(3) + mass*v(3)
          md_ener%total_mass = md_ener%total_mass + REAL(natom,KIND=dp)*mass
       END DO 
       md_ener%vcom = md_ener%vcom / md_ener%total_mass
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
!> \brief Scale internal motion of CORE-SHELL model to the correct temperature
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
!> \par History
!>     MI 02.200s87
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
    IF (simpar%temp_tol > 0.0_dp ) THEN
       CALL scale_velocity(subsys, md_ener, simpar%temp_ext, simpar%temp_tol, iw, error)
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
!> \par History
!>      Creation (29.04.2007,MI)
!>      Cleaned (09.2007) Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE comvel_control(simpar,md_ener,force_env,logger, error)

    TYPE(simpar_type), POINTER               :: simpar
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comvel_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: temp_old, vel_com
    REAL(KIND=dp), DIMENSION(3)              :: vcom_old
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,handle)
    NULLIFY(subsys, para_env)
    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)

    ! Print COMVEL and COM Position
    iw = cp_print_key_unit_nr(logger,force_env%root_section,"MOTION%MD%PRINT%CENTER_OF_MASS",&
         extension=".mdLog",error=error)
    IF (iw>0) THEN
       WRITE (UNIT=iw,FMT="(/,T2,A,(T58,A3,F20.10))")&
            "Centre of mass motion (COM):","x =",md_ener%vcom(1),"y =",md_ener%vcom(2),"z =",md_ener%vcom(3)
    END IF
    CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
         "MOTION%MD%PRINT%CENTER_OF_MASS", error=error)

    ! If requested rescale COMVEL
    IF (simpar%comvel_tol > 0.0_dp ) THEN
       iw = cp_print_key_unit_nr(logger,force_env%root_section,"MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
            extension=".mdLog",error=error)
       vel_com = SQRT(md_ener%vcom(1)**2+md_ener%vcom(2)**2+md_ener%vcom(3)**2)

       ! Subtract the velocity of the COM, if requested
       IF (vel_com > simpar%comvel_tol) THEN
          temp_old = md_ener%temp_part/kelvin
          vcom_old = md_ener%vcom
          CALL reset_vcom( subsys, md_ener, vsubtract=vcom_old, error=error)
          CALL scale_velocity(subsys, md_ener, temp_old, 0.0_dp, iw, error)
          IF (iw>0) THEN
             WRITE (UNIT=iw,FMT="(T2,'MD| ',A,3F16.10,A)") &
                  "Old VCOM = ",vcom_old(1:3)," a.u.",&
                  "New VCOM =", md_ener%vcom(1:3)," a.u"
          END IF
       END IF
       CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
            "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

    CALL timestop(handle)
  END SUBROUTINE comvel_control

! *****************************************************************************
!> \brief Initialize Velocities for MD runs
!> \par History
!>     Teodoro Laino - University of Zurich - 09.2007 [tlaino]
! *****************************************************************************
  SUBROUTINE setup_velocities(force_env, simpar, globenv, md_section, &
       constraint_section, error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: md_section, constraint_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_velocities', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: apply_cns0, failure, &
                                                shell_adiabatic, shell_present
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
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
    print_section => section_vals_get_subs_vals(md_section,"PRINT",error=error)
    apply_cns0 = .FALSE.
    IF (simpar%constraint) THEN
       CALL section_vals_val_get(constraint_section,"CONSTRAINT_INIT",l_val=apply_cns0,error=error)
    END IF
    ! Always initialize velocities and possibly restart them
    CALL force_env_get(force_env, subsys=subsys, cell=cell, cell_ref=cell_ref, para_env=para_env,&
         force_env_section=force_env_section, error=error )
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)

    CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,molecule_kinds_new=molecule_kinds_new,&
         particles=particles,shell_particles=shell_particles,core_particles=core_particles,&
         error=error)
    
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,shell_present=shell_present,&
         shell_adiabatic=shell_adiabatic)

    NULLIFY(particle_set, shell_particle_set, core_particle_set)
    particle_set => particles%els
    IF(shell_present .AND. shell_adiabatic) THEN
       ! shell-model with constraint NYI
       CPPostcondition(.NOT.simpar%constraint,cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(shell_particles),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(core_particles),cp_failure_level,routineP,error,failure)
       shell_particle_set => shell_particles%els
       core_particle_set  =>  core_particles%els
    END IF

    CALL initialize_velocities(simpar,particle_set, molecule_kinds=molecule_kinds_new,&
       force_env=force_env, para_env=para_env, globenv=globenv, print_section=print_section,&
       subsys_section=subsys_section, shell_present=(shell_present.AND.shell_adiabatic), &
       shell_part=shell_particle_set, core_part=core_particle_set, force_rescaling=.FALSE., &
       error=error)
    
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
          CALL update_subsys(subsys_section,force_env,.FALSE.,error)

          CALL initialize_velocities(simpar,particle_set, molecule_kinds=molecule_kinds_new,&
             force_env=force_env, para_env=para_env, globenv=globenv, print_section=print_section, &
             subsys_section=subsys_section, shell_present=(shell_present.AND.shell_adiabatic), &
             shell_part=shell_particle_set, core_part=core_particle_set, force_rescaling=.TRUE., &
             error=error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE setup_velocities

END MODULE md_vel_utils
