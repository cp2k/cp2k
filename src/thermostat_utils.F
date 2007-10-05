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
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE input_constants,                 ONLY: &
       do_region_global, do_region_massive, do_region_molecule, &
       do_thermo_communication, do_thermo_no_communication, isokin_ensemble, &
       langevin_ensemble, nph_uniaxial_damped_ensemble, &
       nph_uniaxial_ensemble, npt_f_ensemble, npt_i_ensemble, nve_ensemble, &
       nvt_ensemble, reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             get_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type
  USE termination,                     ONLY: stop_program
  USE thermostat_types,                ONLY: thermostat_info_type,&
                                             thermostats_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: compute_degrees_of_freedom,&
            setup_thermostat_info    

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
       local_molecules, print_section, gci, error)

    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(cell_type), POINTER                 :: cell
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(global_constraint_type), POINTER    :: gci
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
         local_molecules, simpar%part_region, simpar%ensemble, simpar%nfree, &
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
          ! DIStributed thermostats, no interaction
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

END MODULE thermostat_utils
!******************************************************************************
