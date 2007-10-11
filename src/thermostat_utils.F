!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!***** cp2k/thermostat_utils [1.0] *
!!
!!   NAME
!!     thermostat_utils
!!
!!   FUNCTION
!!     Utilities for thermostats
!!
!!   AUTHOR
!!     teo [tlaino] - University of Zurich - 10.2007
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE thermostat_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE csvr_system_types,               ONLY: csvr_system_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             map_info_type,&
                                             npt_info_type
  USE input_constants,                 ONLY: &
       do_region_global, do_region_massive, do_region_molecule, &
       do_thermo_communication, do_thermo_csvr, do_thermo_no_communication, &
       do_thermo_nose, dyn_coeff_mass_global, isokin_ensemble, &
       langevin_ensemble, nph_uniaxial_damped_ensemble, &
       nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, nve_ensemble, &
       nvt_ensemble, reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             get_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE termination,                     ONLY: stop_program
  USE thermostat_types,                ONLY: thermostat_info_type,&
                                             thermostat_type,&
                                             thermostats_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: compute_degrees_of_freedom,&
            setup_thermostat_info,&
            ke_region_coeffs,&
            ke_region_baro,&
            ke_region_particles,&
            ke_region_shells,&
            vel_rescale_coeffs,&
            vel_rescale_baro,&
            vel_rescale_particles,&
            vel_rescale_shells,&
            get_thermostat_energies


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_utils'

CONTAINS

!******************************************************************************
!!****** thermostat_utils/compute_degrees_of_freedom [1.0] *
!!
!!   NAME
!!    compute_degrees_of_freedom 
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE compute_degrees_of_freedom ( thermostats, cell, simpar, molecule_kind_set,&
       local_molecules, print_section, gci, region, error)

    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(cell_type), POINTER                 :: cell
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(global_constraint_type), POINTER    :: gci
    INTEGER, INTENT(IN)                      :: region
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_degrees_of_freedom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw, natom, nconstraint_ext, &
                                                nconstraint_int, nfree, &
                                                nrestraints_int
    TYPE(cp_logger_type), POINTER            :: logger

    simpar % nfree = 6 - SUM(cell%perd(1:3))
    ! Collect info about thermostats
    CALL setup_thermostat_info(thermostats%thermostat_info_part, molecule_kind_set,&
         local_molecules, region, simpar%ensemble, simpar%nfree, &
         error=error)

    ! Saving this value of simpar preliminar to the real count of constraints..
    simpar % nfree_rot_transl = simpar % nfree

    ! compute the total number of degrees of freedom for temperature
    CALL get_molecule_kind_set(molecule_kind_set=molecule_kind_set,&
         natom=natom, nconstraint=nconstraint_int, nrestraints=nrestraints_int)
    nconstraint_ext = gci%ntot - gci%nrestraint
    nfree = 3*natom - nconstraint_int - nconstraint_ext - simpar%nfree

    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)
    IF ( iw > 0 ) THEN
       WRITE ( iw, '( /, A )' ) ' Calculation of degrees of freedom'
       WRITE ( iw, '( T48, A, T71, I10 )' ) '       Number of atoms:', natom
       WRITE ( iw, '( T33, A, T71, I10 )' ) &
            ' Number of Intramolecular constraints:', nconstraint_int
       WRITE ( iw, '( T33, A, T71, I10 )' ) &
            ' Number of Intermolecular constraints:', nconstraint_ext
       WRITE ( iw, '( T34, A, T71, I10 )' ) &
            ' Invariants(translation + rotations):', simpar % nfree
       WRITE ( iw, '( T48, A, T71, I10 )' ) '    Degrees of freedom:', nfree
       WRITE ( iw, '()' )
       WRITE ( iw, '( /, A )' ) ' Restraints Information'
       WRITE ( iw, '( T33, A, T71, I10 )' ) &
            '  Number of Intramolecular restraints:', nrestraints_int
       WRITE ( iw, '( T33, A, T71, I10 )' ) &
            '  Number of Intermolecular restraints:', gci%nrestraint
    END IF
    CALL cp_print_key_finished_output(iw,logger,print_section,&
               "PROGRAM_RUN_INFO", error=error)
    simpar%nfree = nfree

  END SUBROUTINE compute_degrees_of_freedom

!******************************************************************************
!!****** thermostat_utils/setup_thermostat_info [1.0] *
!!
!!   NAME
!!     setup_thermostat_info
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE setup_thermostat_info(thermostat_info, molecule_kind_set,&
       local_molecules, region, ensemble, nfree, shell, error)
    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, INTENT(IN)                      :: region, ensemble
    INTEGER, INTENT(INOUT), OPTIONAL         :: nfree
    LOGICAL, INTENT(IN), OPTIONAL            :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_thermostat_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dis_type, i, ikind, natom, &
                                                nkind, nmol_local, nmolecule, &
                                                nshell, number, &
                                                sum_of_thermostats
    LOGICAL                                  :: do_shell, nointer
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(molecule_kind)
    nkind = SIZE ( molecule_kind_set )
    do_shell = .FALSE.
    IF (PRESENT(shell)) do_shell = shell
    ! Counting the global number of thermostats
    sum_of_thermostats = 0
    ! variable to denote independent thermostats (no communication necessary)
    nointer = .TRUE.
    number  = 0

    SELECT CASE(ensemble)
    CASE DEFAULT
       CALL stop_program (routineN, 'unknown ensemble ')
    CASE(nve_ensemble, isokin_ensemble, nph_uniaxial_ensemble, nph_uniaxial_damped_ensemble,&
         reftraj_ensemble, langevin_ensemble)
       ! Do Nothing
    CASE(nvt_ensemble, npt_i_ensemble, npt_f_ensemble)

       SELECT CASE(region)
       CASE(do_region_global)
          ! Global Thermostat
          nointer = .FALSE.
          sum_of_thermostats = 1
       CASE(do_region_molecule)
          ! Molecular Thermostat
          DO i = 1, nkind
             molecule_kind => molecule_kind_set ( i )
             CALL get_molecule_kind ( molecule_kind, nmolecule=nmolecule, nshell=nshell )
             IF ((do_shell).AND.(nshell==0)) nmolecule=0
             sum_of_thermostats = sum_of_thermostats + nmolecule
          END DO
          ! If we have ONE kind and ONE molecule, then effectively we have a GLOBAL thermostat
          ! and the degrees of freedom will be computed correctly for this special case
          IF ( ( nmolecule == 1 ) .AND. ( nkind == 1 ) ) nointer = .FALSE.
       CASE(do_region_massive)
          ! Massive Thermostat
          DO i = 1, nkind
             molecule_kind => molecule_kind_set ( i )
             CALL get_molecule_kind ( molecule_kind, nmolecule=nmolecule, &
                  natom=natom, nshell=nshell )
             IF (do_shell) natom = nshell
             sum_of_thermostats = sum_of_thermostats + 3 * natom * nmolecule
          END DO
       END SELECT
       number = sum_of_thermostats

       ! Here we decide on what parallel algorithm to use.
       ! if there are only massive and molecule type thermostats we can use
       ! a local scheme, in cases involving any combination with a
       ! global thermostat we assume a coupling of  degrees of freedom
       ! from different processors
       IF ( nointer ) THEN
          ! Distributed thermostats, no interaction
          dis_type = do_thermo_no_communication
          ! we only count thermostats on this processor
          number = 0
          DO ikind = 1, nkind
             nmol_local = local_molecules % n_el ( ikind )
             molecule_kind => molecule_kind_set ( ikind )
             CALL get_molecule_kind ( molecule_kind, natom=natom, nshell=nshell)
             IF (do_shell) THEN 
                natom = nshell
                IF (nshell==0) nmol_local = 0
             END IF
             IF (region==do_region_molecule) THEN
                number = number + nmol_local
             ELSE IF (region==do_region_massive) THEN
                number = number + 3*nmol_local*natom
             ELSE
                CALL stop_program (routineP, 'Invalid region setup' )
             END IF
          END DO
       ELSE
          ! REPlicated thermostats, INTERacting via communication
          dis_type = do_thermo_communication
       END IF
       
       IF (PRESENT(nfree)) THEN
          IF ((sum_of_thermostats>1).OR.(dis_type==do_thermo_no_communication)) THEN
             ! re-initializing simpar%nfree to zero because of multiple thermostats
             nfree = 0.0_dp
          END IF
       END IF
    END SELECT

    ! Saving information about thermostats
    thermostat_info%sum_of_thermostats    = sum_of_thermostats
    thermostat_info%number_of_thermostats = number
    thermostat_info%dis_type              = dis_type

  END SUBROUTINE setup_thermostat_info

!******************************************************************************
!!****** thermostat_utils/ke_region_coeffs [1.0] *
!!
!!   NAME
!!    ke_region_coeffs
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE ke_region_coeffs(map_info, dyn_coeff_set, coeffs, group, flops, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ke_region_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, icoef, iel, ikind, nkind
    LOGICAL                                  :: failure, gmass

    failure = .FALSE.
    gmass=(dyn_coeff_set%mass_storage==dyn_coeff_mass_global)
    map_info%v_scale = 1.0_dp
    map_info%s_kin   = 0.0_dp
    i  = 0               
    nkind = SIZE (dyn_coeff_set%coeffs_of_kind)
    DO ikind=1,nkind
       IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
          coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
          DO iel = 1, coeffs%n_els
             DO icoef = 1, coeffs%ncoef_atom
                i = i + 1
                IF(gmass) THEN
                   map_info%p_kin(1,i)%point  = map_info%p_kin(1,i)%point + &
                        dyn_coeff_set%global_mass*coeffs%vel(iel,icoef)**2
                ELSE
                   map_info%p_kin(1,i)%point  = map_info%p_kin(1,i)%point + &
                        coeffs%masses(iel,icoef)*coeffs%vel(iel,icoef)**2
                END IF
             END DO
          END DO
       END IF
    END DO
    IF ( map_info%dis_type == do_thermo_communication )  CALL mp_sum (map_info%s_kin, group )
    flops = flops + REAL(i*3, KIND=dp)

  END SUBROUTINE ke_region_coeffs

!******************************************************************************
!!****** thermostat_utils/vel_rescale_coeffs [1.0] *
!!
!!   NAME
!!     vel_rescale_coeffs
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE vel_rescale_coeffs(map_info, dyn_coeff_set, coeffs, flops, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vel_rescale_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, icoef, iel, ikind, nkind
    LOGICAL                                  :: failure

    failure = .FALSE.
    i=0
    nkind = SIZE (dyn_coeff_set%coeffs_of_kind)
    ! Now scale the coefficients velocities
    DO ikind=1,nkind
       IF(ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
          coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
          DO iel = 1, coeffs%n_els
             DO icoef = 1, coeffs%ncoef_atom
                i = i + 1
                coeffs%vel(iel,icoef)=coeffs%vel(iel,icoef)*map_info%p_scale(1,i)%point
             END DO
          END DO
       END IF
    END DO
    flops = flops + REAL (i,KIND=dp)

  END SUBROUTINE vel_rescale_coeffs

!******************************************************************************
!!****** thermostat_utils/ke_region_baro [1.0] *
!!
!!   NAME
!!    ke_region_baro
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE ke_region_baro(map_info, npt, group, flops, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ke_region_baro', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncoef
    LOGICAL                                  :: failure

    failure = .FALSE.
    map_info%v_scale = 1.0_dp
    map_info%s_kin = 0.0_dp
    ncoef = 0
    DO i = 1, SIZE ( npt, 1 )
       DO j = 1, SIZE ( npt, 2 )
          ncoef = ncoef + 1
          map_info%p_kin( 1, ncoef ) % point = map_info%p_kin ( 1, ncoef ) % point &
               + npt ( i, j ) % mass * npt ( i, j )%v**2
       END DO
    END DO

    IF (map_info%dis_type == do_thermo_communication ) CALL mp_sum (map_info%s_kin, group )
    flops = flops + REAL(ncoef*3, KIND=dp)

  END SUBROUTINE ke_region_baro

!******************************************************************************
!!****** thermostat_utils/vel_rescale_baro [1.0] *
!!
!!   NAME
!!     vel_rescale_baro
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE vel_rescale_baro(map_info, npt, flops, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vel_rescale_baro', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncoef

    ncoef = 0
    DO i = 1,  SIZE ( npt, 1 )
       DO j = 1,  SIZE ( npt, 2 )
          ncoef = ncoef + 1
          npt ( i, j ) % v = npt ( i, j ) % v * map_info%p_scale( 1, ncoef ) % point
       END DO
    END DO
    flops = flops + REAL (ncoef,KIND=dp)

  END SUBROUTINE vel_rescale_baro

!******************************************************************************
!!****** thermostat_utils/ke_region_particles [1.0] *
!!
!!   NAME
!!     ke_region_particles
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, flops, vel, error)
    
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ke_region_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, ii, ikind, imol, &
                                                imol_local, ipart, last_atom, &
                                                nmol_local
    LOGICAL                                  :: failure, present_vel
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    map_info%v_scale = 1.0_dp
    map_info%s_kin   = 0.0_dp
    present_vel = PRESENT(vel)
    ii = 0
    DO ikind = 1, SIZE ( molecule_kind_set )
       nmol_local = local_molecules % n_el ( ikind )
       DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          molecule => molecule_set ( imol )
          CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom ) 
          DO ipart = first_atom, last_atom
             ii = ii + 1
             atomic_kind => particle_set ( ipart ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
             IF (present_vel) THEN
                map_info%p_kin(1,ii)%point = map_info%p_kin(1,ii)%point+mass*vel(1,ipart)**2
                map_info%p_kin(2,ii)%point = map_info%p_kin(2,ii)%point+mass*vel(2,ipart)**2
                map_info%p_kin(3,ii)%point = map_info%p_kin(3,ii)%point+mass*vel(3,ipart)**2
             ELSE
                map_info%p_kin(1,ii)%point=map_info%p_kin(1,ii)%point+mass*particle_set(ipart)%v(1)**2
                map_info%p_kin(2,ii)%point=map_info%p_kin(2,ii)%point+mass*particle_set(ipart)%v(2)**2
                map_info%p_kin(3,ii)%point=map_info%p_kin(3,ii)%point+mass*particle_set(ipart)%v(3)**2
             END IF
          END DO
       END DO
    END DO
   
    IF (map_info%dis_type == do_thermo_communication )  CALL mp_sum (map_info%s_kin, group)
    flops = flops + REAL ( ii * 9,KIND=dp)

  END SUBROUTINE ke_region_particles
  
!******************************************************************************
!!****** thermostat_utils/vel_rescale_particles [1.0] *
!!
!!   NAME
!!     vel_rescale_particles
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE vel_rescale_particles(map_info, molecule_kind_set, molecule_set,&
       particle_set, local_molecules, shell_adiabatic, flops, shell_particle_set,&
       core_particle_set, vel, shell_vel, core_vel, error)
    
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vel_rescale_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, ii, ikind, imol, &
                                                imol_local, ipart, jj, &
                                                last_atom, nmol_local, &
                                                shell_index
    LOGICAL                                  :: failure, present_vel
    REAL(KIND=dp)                            :: fac_massc, fac_masss, mass, &
                                                vc(3), vs(3)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(shell_kind_type), POINTER           :: shell

    ii = 0
    jj = 0
    failure = .FALSE.
    present_vel = PRESENT(vel)
    ! Just few checks for consistency
    IF (present_vel) THEN
       IF (shell_adiabatic) THEN
          CPPostcondition(PRESENT(shell_vel),cp_failure_level,routineP,error,failure)
          CPPostcondition(PRESENT(core_vel),cp_failure_level,routineP,error,failure)
       END IF
    ELSE
       IF (shell_adiabatic) THEN
          CPPostcondition(PRESENT(shell_particle_set),cp_failure_level,routineP,error,failure)
          CPPostcondition(PRESENT(core_particle_set),cp_failure_level,routineP,error,failure)
       END IF
    END IF
    Kind: DO ikind = 1, SIZE ( molecule_kind_set )
       nmol_local = local_molecules % n_el ( ikind )
       Mol_local: DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          molecule => molecule_set ( imol )
          CALL get_molecule ( molecule, first_atom=first_atom, last_atom=last_atom ) 
          Particle: DO ipart = first_atom, last_atom
             ii = ii + 1
             IF (present_vel) THEN
                vel(1,ipart) =  vel(1,ipart) * map_info%p_scale(1,ii)%point
                vel(2,ipart) =  vel(2,ipart) * map_info%p_scale(2,ii)%point
                vel(3,ipart) =  vel(3,ipart) * map_info%p_scale(3,ii)%point
             ELSE
                particle_set(ipart)%v(1) =  particle_set(ipart)%v(1)*map_info%p_scale(1,ii)%point
                particle_set(ipart)%v(2) =  particle_set(ipart)%v(2)*map_info%p_scale(2,ii)%point
                particle_set(ipart)%v(3) =  particle_set(ipart)%v(3)*map_info%p_scale(3,ii)%point
             END IF
                ! If Shell Adiabatic then apply the NHC thermostat also to the Shells
                IF(shell_adiabatic) THEN
                   shell_index=particle_set ( ipart ) %shell_index
                   IF(shell_index/=0) THEN
                      jj = jj + 2
                      atomic_kind => particle_set ( ipart ) % atomic_kind
                      CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass,shell=shell)
                      fac_masss = shell%mass_shell/mass 
                      fac_massc = shell%mass_core/mass
                      IF (present_vel) THEN
                         vs(1:3) = shell_vel( 1:3,shell_index )
                         vc(1:3) = core_vel( 1:3,shell_index )
                         shell_vel( 1,shell_index ) = vel( 1, ipart ) + fac_massc*(vs(1)-vc(1))
                         shell_vel( 2,shell_index ) = vel( 2, ipart ) + fac_massc*(vs(2)-vc(2))
                         shell_vel( 3,shell_index ) = vel( 3, ipart ) + fac_massc*(vs(3)-vc(3))
                         core_vel ( 1,shell_index ) = vel( 1, ipart ) + fac_masss*(vc(1)-vs(1))
                         core_vel ( 2,shell_index ) = vel( 2, ipart ) + fac_masss*(vc(2)-vs(2))
                         core_vel ( 3,shell_index ) = vel( 3, ipart ) + fac_masss*(vc(3)-vs(3))                        
                      ELSE
                         vs(1:3) = shell_particle_set(shell_index)%v(1:3)
                         vc(1:3) = core_particle_set(shell_index)%v(1:3)
                         shell_particle_set(shell_index)%v(1)= particle_set(ipart)%v(1) + fac_massc*(vs(1)-vc(1))
                         shell_particle_set(shell_index)%v(2)= particle_set(ipart)%v(2) + fac_massc*(vs(2)-vc(2))
                         shell_particle_set(shell_index)%v(3)= particle_set(ipart)%v(3) + fac_massc*(vs(3)-vc(3))
                         core_particle_set (shell_index)%v(1)= particle_set(ipart)%v(1) + fac_masss*(vc(1)-vs(1))
                         core_particle_set (shell_index)%v(2)= particle_set(ipart)%v(2) + fac_masss*(vc(2)-vs(2))
                         core_particle_set (shell_index)%v(3)= particle_set(ipart)%v(3) + fac_masss*(vc(3)-vs(3))
                      END IF
                   END IF
             END IF
          END DO Particle
       END DO Mol_local
    END DO Kind
    flops = flops + REAL ((ii+jj*2)*3,KIND=dp)

  END SUBROUTINE vel_rescale_particles

!******************************************************************************
!!****** thermostat_utils/ke_region_shells [1.0] *
!!
!!   NAME
!!     ke_region_shells
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE ke_region_shells(map_info, particle_set, atomic_kind_set,&
       local_particles, group, flops, core_particle_set, shell_particle_set,&
       core_vel, shell_vel, error)
    
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    TYPE(particle_type), OPTIONAL, POINTER   :: core_particle_set(:), &
                                                shell_particle_set(:)
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: core_vel(:,:), shell_vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ke_region_shells', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii, iparticle, iparticle_kind, iparticle_local, &
      nparticle_kind, nparticle_local, shell_index
    LOGICAL                                  :: failure, is_shell, present_vel
    REAL(dp)                                 :: mass, mu_mass, v_sc(3)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    failure = .FALSE.
    present_vel = PRESENT(shell_vel)
    ! Preliminary checks for consistency usage
    IF (present_vel) THEN
       CPPostcondition(PRESENT(core_vel),cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(PRESENT(shell_particle_set),cp_failure_level,routineP,error,failure)
       CPPostcondition(PRESENT(core_particle_set),cp_failure_level,routineP,error,failure)
    END IF
    ! get force on first thermostat for all the chains in the system.
    map_info%v_scale = 1.0_dp
    map_info%s_kin = 0.0_dp
    ii = 0
    
    nparticle_kind = SIZE(atomic_kind_set)
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass,shell_active=is_shell,shell=shell)
       IF(is_shell) THEN
          mu_mass = shell%mass_shell*shell%mass_core/mass 
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             shell_index = particle_set(iparticle)%shell_index
             ii = ii + 1
             IF (present_vel) THEN
                v_sc(1) = core_vel(1,shell_index)-shell_vel(1,shell_index)
                v_sc(2) = core_vel(2,shell_index)-shell_vel(2,shell_index)
                v_sc(3) = core_vel(3,shell_index)-shell_vel(3,shell_index)
                map_info%p_kin(1,ii)%point = map_info%p_kin(1,ii)%point + mu_mass * v_sc(1)**2
                map_info%p_kin(2,ii)%point = map_info%p_kin(2,ii)%point + mu_mass * v_sc(2)**2
                map_info%p_kin(3,ii)%point = map_info%p_kin(3,ii)%point + mu_mass * v_sc(3)**2               
             ELSE
                v_sc(1) = core_particle_set(shell_index)%v(1)-shell_particle_set(shell_index)%v(1)
                v_sc(2) = core_particle_set(shell_index)%v(2)-shell_particle_set(shell_index)%v(2)
                v_sc(3) = core_particle_set(shell_index)%v(3)-shell_particle_set(shell_index)%v(3)
                map_info%p_kin(1,ii)%point = map_info%p_kin(1,ii)%point + mu_mass * v_sc(1)**2
                map_info%p_kin(2,ii)%point = map_info%p_kin(2,ii)%point + mu_mass * v_sc(2)**2
                map_info%p_kin(3,ii)%point = map_info%p_kin(3,ii)%point + mu_mass * v_sc(3)**2
             END IF
          END DO
       END IF
    END DO
    IF ( map_info%dis_type == do_thermo_communication )  CALL mp_sum (map_info%s_kin, group )
    flops = flops + REAL(ii* 9,KIND=dp)

  END SUBROUTINE ke_region_shells

!******************************************************************************
!!****** thermostat_utils/vel_rescale_shells [1.0] *
!!
!!   NAME
!!     vel_rescale_shells 
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     10.2007 [tlaino] - Teodoro Laino - University of Zurich
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE vel_rescale_shells(map_info, atomic_kind_set, particle_set, local_particles, &
       flops, shell_particle_set, core_particle_set, shell_vel, core_vel, vel, error)

    TYPE(map_info_type), POINTER             :: map_info
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(INOUT)             :: flops
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set(:), &
                                                core_particle_set(:)
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: shell_vel(:,:), &
                                                core_vel(:,:), vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'vel_rescale_shells', &
      routineP = moduleN//':'//routineN

    INTEGER :: ii, iparticle, iparticle_kind, iparticle_local, &
      nparticle_kind, nparticle_local, shell_index
    LOGICAL                                  :: failure, is_shell, present_vel
    REAL(dp)                                 :: mass, massc, masss, umass, &
                                                v(3), vc(3), vs(3)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    failure = .FALSE.
    present_vel = PRESENT(vel)
    ! Preliminary checks for consistency usage
    IF (present_vel) THEN
       CPPostcondition(PRESENT(shell_vel),cp_failure_level,routineP,error,failure)
       CPPostcondition(PRESENT(core_vel),cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(PRESENT(shell_particle_set),cp_failure_level,routineP,error,failure)
       CPPostcondition(PRESENT(core_particle_set),cp_failure_level,routineP,error,failure)
    END IF
    ii = 0
    nparticle_kind = SIZE(atomic_kind_set)
    ! now scale the core-shell velocities
    Kind: DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass,shell_active=is_shell,shell=shell)
       IF(is_shell) THEN
          umass = 1.0_dp/mass
          masss = shell%mass_shell*umass
          massc = shell%mass_core*umass
          
          nparticle_local = local_particles%n_el(iparticle_kind)
          Particles: DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             shell_index = particle_set(iparticle)%shell_index
             ii = ii + 1
             IF (present_vel) THEN
                vc(1:3) = core_vel(1:3,shell_index)
                vs(1:3) = shell_vel(1:3,shell_index)
                v(1:3) =  vel(1:3,iparticle)                
                shell_vel(1,shell_index) = v(1) + map_info%p_scale ( 1, ii ) % point * massc*(vs(1)-vc(1))
                shell_vel(2,shell_index) = v(2) + map_info%p_scale ( 2, ii ) % point * massc*(vs(2)-vc(2))
                shell_vel(3,shell_index) = v(3) + map_info%p_scale ( 3, ii ) % point * massc*(vs(3)-vc(3))
                core_vel(1,shell_index)  = v(1) + map_info%p_scale ( 1, ii ) % point * masss*(vc(1)-vs(1))
                core_vel(2,shell_index)  = v(2) + map_info%p_scale ( 2, ii ) % point * masss*(vc(2)-vs(2))
                core_vel(3,shell_index)  = v(3) + map_info%p_scale ( 3, ii ) % point * masss*(vc(3)-vs(3))               
             ELSE
                vc(1:3) = core_particle_set(shell_index)%v(1:3)
                vs(1:3) = shell_particle_set(shell_index)%v(1:3)
                v(1:3)  = particle_set(iparticle)%v(1:3)
                shell_particle_set(shell_index)%v(1) = v(1) + map_info%p_scale(1,ii)%point * massc*(vs(1)-vc(1))
                shell_particle_set(shell_index)%v(2) = v(2) + map_info%p_scale(2,ii)%point * massc*(vs(2)-vc(2))
                shell_particle_set(shell_index)%v(3) = v(3) + map_info%p_scale(3,ii)%point * massc*(vs(3)-vc(3))
                core_particle_set(shell_index)%v(1)  = v(1) + map_info%p_scale(1,ii)%point * masss*(vc(1)-vs(1))
                core_particle_set(shell_index)%v(2)  = v(2) + map_info%p_scale(2,ii)%point * masss*(vc(2)-vs(2))
                core_particle_set(shell_index)%v(3)  = v(3) + map_info%p_scale(3,ii)%point * masss*(vc(3)-vs(3))
             END IF
          END DO Particles
       END IF
    END DO Kind
    flops = flops + REAL (ii*12,KIND=dp)

  END SUBROUTINE vel_rescale_shells

  !******************************************************************************
  !!****** thermostat_utils/get_nhc_energies [1.0] *
  !!
  !!   NAME
  !!     get_nhc_energies
  !!
  !!   FUNCTION
  !!     Calculates kinetic energy and potential energy
  !!     of the nhc variables
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!***8*************************************************************************
  SUBROUTINE get_nhc_energies ( nhc, nhc_pot, nhc_kin, int_grp, error )
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    REAL(KIND=dp), INTENT(OUT)               :: nhc_pot, nhc_kin
    INTEGER                                  :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_nhc_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l, n

    nhc_kin = 0.0_dp
    nhc_pot = 0.0_dp
    DO n = 1, nhc%loc_num_nhc
       DO l = 1, nhc % nhc_len
          nhc_kin = nhc_kin +  &
               0.5_dp * nhc % nvt ( l, n ) % mass * &
               nhc % nvt ( l, n ) % v * nhc % nvt ( l, n ) % v

          nhc_pot = nhc_pot + nhc % nvt ( l, n ) % nkt *  &
               nhc % nvt ( l, n ) % eta
       END DO
    END DO
    IF (nhc%map_info%dis_type==do_thermo_no_communication) THEN
       CALL mp_sum(nhc_kin,int_grp)
       CALL mp_sum(nhc_pot,int_grp)
    END IF
  END SUBROUTINE get_nhc_energies

  !******************************************************************************
  !!****** thermostat_utils/get_csvr_energies [1.0] *
  !!
  !!   NAME
  !!     get_csvr_energies
  !!
  !!   FUNCTION
  !!     Calculates kinetic energy and potential energy
  !!     of the csvr thermostat
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] - University of Zurich 10.2007
  !!
  !!   MODIFICATION HISTORY
  !!    
  !!
  !!*****************************************************************************
  SUBROUTINE get_csvr_energies ( csvr, csvr_kin, int_grp, error )
    TYPE(csvr_system_type), POINTER          :: csvr
    REAL(KIND=dp), INTENT(OUT)               :: csvr_kin
    INTEGER                                  :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_csvr_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n

    csvr_kin = 0.0_dp
    DO n = 1, csvr%loc_num_csvr
       csvr_kin = csvr_kin + csvr%nvt(n)%thermostat_energy
    END DO
    IF (csvr%map_info%dis_type==do_thermo_no_communication) THEN
       CALL mp_sum(csvr_kin,int_grp)
    END IF
  END SUBROUTINE get_csvr_energies

  !******************************************************************************
  !!****** thermostat_utils/get_thermostat_energies [1.0] *
  !!
  !!   NAME
  !!     get_thermostat_energies
  !!
  !!   FUNCTION
  !!     Calculates energy associated with a thermostat
  !!
  !!   AUTHOR
  !!     Teodoro Laino [tlaino] - 10.2007 - University of Zurich
  !!
  !!   MODIFICATION HISTORY
  !!     
  !!
  !!*****************************************************************************
  SUBROUTINE get_thermostat_energies ( thermostat, thermostat_pot, thermostat_kin, int_grp, error )
    TYPE(thermostat_type), POINTER           :: thermostat
    REAL(KIND=dp), INTENT(OUT)               :: thermostat_pot, thermostat_kin
    INTEGER, INTENT(IN)                      :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_thermostat_energies', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    
    thermostat_pot = 0.0_dp
    thermostat_kin = 0.0_dp
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Energy associated with the Nose-Hoover thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL get_nhc_energies (thermostat%nhc, thermostat_pot, thermostat_kin, int_grp, error )
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN
          ! Energy associated with the CSVR thermostat
          CALL get_csvr_energies(thermostat%csvr, thermostat_kin, int_grp, error)
       END IF
    END IF

  END SUBROUTINE get_thermostat_energies

END MODULE thermostat_utils
!******************************************************************************
