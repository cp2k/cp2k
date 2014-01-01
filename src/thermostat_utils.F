!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for thermostats
!> \author teo [tlaino] - University of Zurich - 10.2007
! *****************************************************************************
MODULE thermostat_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             map_info_type,&
                                             npt_info_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_constr_atomic, do_constr_molec, do_region_defined, &
       do_region_global, do_region_massive, do_region_molecule, do_thermo_al, &
       do_thermo_communication, do_thermo_csvr, do_thermo_gle, &
       do_thermo_no_communication, do_thermo_nose, isokin_ensemble, &
       langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_adiabatic_ensemble, nvt_ensemble, &
       reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             get_molecule_kind_set,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             molecule_type
  USE motion_utils,                    ONLY: rot_ana
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: femtoseconds
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_program
  USE thermostat_types,                ONLY: thermostat_info_type,&
                                             thermostat_type,&
                                             thermostats_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: compute_degrees_of_freedom,&
            compute_nfree,&
            setup_thermostat_info,&
            setup_adiabatic_thermostat_info,&
            ke_region_baro,&
            ke_region_particles,&
            ke_region_shells,&
            vel_rescale_baro,&
            vel_rescale_particles,&
            vel_rescale_shells,&
            get_thermostat_energies,&
            get_nhc_energies,&
            get_kin_energies,&
            communication_thermo_low2,&
            print_thermostats_status,&
            momentum_region_particles

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_utils'

CONTAINS

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE compute_nfree (  cell, simpar, molecule_kind_set,&
       local_molecules, molecules, print_section,particles, gci, error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_nfree', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natom, nconstraint_ext, &
                                                nconstraint_int, &
                                                nrestraints_int, rot_dof, &
                                                roto_trasl_dof
    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Retrieve information on number of atoms, constraints (external and internal)
    CALL get_molecule_kind_set(molecule_kind_set=molecule_kind_set,&
         natom=natom, nconstraint=nconstraint_int, nrestraints=nrestraints_int)

    ! Compute degrees of freedom
    CALL rot_ana(particles%els, dof=roto_trasl_dof, rot_dof=rot_dof, &
                 print_section=print_section, keep_rotations=.FALSE.,&
                 mass_weighted=.TRUE., natoms=natom, error=error)

    roto_trasl_dof = roto_trasl_dof - MIN(SUM(cell%perd(1:3)),rot_dof)

    ! Saving this value of simpar preliminar to the real count of constraints..
    simpar%nfree_rot_transl = roto_trasl_dof

    ! compute the total number of degrees of freedom for temperature
    nconstraint_ext = gci%ntot - gci%nrestraint
    simpar%nfree = 3*natom - nconstraint_int - nconstraint_ext - roto_trasl_dof

  END SUBROUTINE compute_nfree

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE compute_degrees_of_freedom ( thermostats, cell, simpar, molecule_kind_set,&
       local_molecules, molecules, particles, print_section, region_sections, gci,&
       region, qmmm_env, error)

    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(cell_type), POINTER                 :: cell
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: print_section, region_sections
    TYPE(global_constraint_type), POINTER    :: gci
    INTEGER, INTENT(IN)                      :: region
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_degrees_of_freedom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iw, natom, nconstraint_ext, &
                                                nconstraint_int, &
                                                nrestraints_int, rot_dof, &
                                                roto_trasl_dof
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    ! Retrieve information on number of atoms, constraints (external and internal)
    CALL get_molecule_kind_set(molecule_kind_set=molecule_kind_set,&
         natom=natom, nconstraint=nconstraint_int, nrestraints=nrestraints_int)

    ! Compute degrees of freedom
    CALL rot_ana(particles%els, dof=roto_trasl_dof, rot_dof=rot_dof, &
                 print_section=print_section, keep_rotations=.FALSE.,&
                 mass_weighted=.TRUE., natoms=natom, error=error)

    roto_trasl_dof = roto_trasl_dof - MIN(SUM(cell%perd(1:3)),rot_dof)

    ! Collect info about thermostats
    CALL setup_thermostat_info(thermostats%thermostat_info_part, molecule_kind_set,&
         local_molecules, molecules, particles, region, simpar%ensemble, roto_trasl_dof, &
         region_sections=region_sections, qmmm_env=qmmm_env, error=error)

    ! Saving this value of simpar preliminar to the real count of constraints..
    simpar%nfree_rot_transl = roto_trasl_dof

    ! compute the total number of degrees of freedom for temperature
    nconstraint_ext = gci%ntot - gci%nrestraint
    simpar%nfree = 3*natom - nconstraint_int - nconstraint_ext - roto_trasl_dof

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
            ' Invariants(translation + rotations):', roto_trasl_dof
       WRITE ( iw, '( T48, A, T71, I10 )' ) '    Degrees of freedom:', simpar%nfree
       WRITE ( iw, '()' )
       WRITE ( iw, '( /, A )' ) ' Restraints Information'
       WRITE ( iw, '( T33, A, T71, I10 )' ) &
            '  Number of Intramolecular restraints:', nrestraints_int
       WRITE ( iw, '( T33, A, T71, I10 )' ) &
            '  Number of Intermolecular restraints:', gci%nrestraint
    END IF
    CALL cp_print_key_finished_output(iw,logger,print_section,&
               "PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE compute_degrees_of_freedom

! *****************************************************************************
!> \author 10.2011  CJM - PNNL
! *****************************************************************************
  SUBROUTINE setup_adiabatic_thermostat_info(thermostat_info, molecule_kind_set, local_molecules,&
       molecules, particles, region, ensemble, nfree, shell, region_sections, qmmm_env,&
       error)
    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: region, ensemble
    INTEGER, INTENT(INOUT), OPTIONAL         :: nfree
    LOGICAL, INTENT(IN), OPTIONAL            :: shell
    TYPE(section_vals_type), POINTER         :: region_sections
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'setup_adiabatic_thermostat_info', &
      routineP = moduleN//':'//routineN

    INTEGER :: dis_type, first_atom, i, ikind, imol, imol_global, ipart, &
      itherm, katom, last_atom, natom, natom_local, nkind, nmol_local, &
      nmol_per_kind, nmolecule, nshell, number, stat, sum_of_thermostats
    INTEGER, POINTER                         :: molecule_list(:), &
                                                thermolist(:)
    LOGICAL                                  :: check, do_shell, failure, &
                                                nointer, on_therm
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule, molecule_set( : )

    failure = .FALSE.
    NULLIFY(molecule_kind,molecule,thermostat_info%map_loc_thermo_gen,thermolist)
    nkind = SIZE ( molecule_kind_set )
    do_shell = .FALSE.
    IF (PRESENT(shell)) do_shell = shell
    ! Counting the global number of thermostats
    sum_of_thermostats = 0
    ! Variable to denote independent thermostats (no communication necessary)
    nointer  = .TRUE.
    check    = .TRUE.
    number   = 0
    dis_type = do_thermo_no_communication

    CALL get_adiabatic_region_info(region_sections, sum_of_thermostats,&
         thermolist=thermolist,&
         local_molecules=local_molecules, molecule_kind_set=molecule_kind_set,&
         molecules=molecules, particles=particles, qmmm_env=qmmm_env, error=error)

!    map_loc_thermo_gen=>thermostat_info%map_loc_thermo_gen
    molecule_set => molecules%els
    SELECT CASE(ensemble)
    CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,'Unknown ensemble')
    CASE(nvt_adiabatic_ensemble)
        SELECT CASE(region)
        CASE(do_region_global)
           ! Global Thermostat
           nointer = .FALSE.
           sum_of_thermostats = 1
        CASE(do_region_molecule)
           ! Molecular Thermostat
           itherm=0
           DO ikind = 1, nkind
             molecule_kind => molecule_kind_set ( ikind )
             nmol_per_kind =  local_molecules % n_el ( ikind )
             CALL get_molecule_kind ( molecule_kind, natom = natom, &
                                      molecule_list=molecule_list )
! use thermolist ( ipart ) to get global indexing correct
             DO imol_global = 1, SIZE ( molecule_list )
               molecule => molecule_set ( molecule_list ( imol_global ) )
               CALL get_molecule ( molecule, first_atom=first_atom, &
                                   last_atom=last_atom )
               on_therm=.TRUE.
               DO katom = first_atom, last_atom
                 IF ( thermolist ( katom ) == HUGE ( 0 ) ) THEN
                    on_therm=.FALSE.
                    EXIT
                 END IF
               END DO
               IF ( on_therm ) THEN
                 itherm=itherm + 1
                 DO katom = first_atom, last_atom
                   thermolist ( katom ) = itherm
                 END DO
               END IF
             END DO
           END DO
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

    natom_local = 0
    DO ikind = 1, SIZE(molecule_kind_set)
       nmol_per_kind = local_molecules%n_el(ikind)
       DO imol = 1, nmol_per_kind
          i = local_molecules%list(ikind)%array(imol)
          molecule => molecule_set(i)
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
          DO ipart = first_atom, last_atom
             natom_local = natom_local + 1
          END DO
       END DO
    END DO

    ! Now map the local atoms with the corresponding thermostat
    ALLOCATE(thermostat_info%map_loc_thermo_gen(natom_local),stat=stat)
    thermostat_info%map_loc_thermo_gen  = HUGE ( 0 )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    natom_local = 0
    DO ikind = 1, SIZE(molecule_kind_set)
       nmol_per_kind = local_molecules%n_el(ikind)
       DO imol = 1, nmol_per_kind
          i = local_molecules%list(ikind)%array(imol)
          molecule => molecule_set(i)
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
          DO ipart = first_atom, last_atom
             natom_local = natom_local + 1
! only map the correct region to the thermostat
             IF ( thermolist (ipart ) /= HUGE ( 0 ) ) &
             thermostat_info%map_loc_thermo_gen(natom_local) = thermolist(ipart)
          END DO
       END DO
    END DO
        ! Here we decide which parallel algorithm to use.
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
                 CALL stop_program(routineN,moduleN,__LINE__,'Invalid region setup')
              END IF
           END DO
        ELSE
           ! REPlicated thermostats, INTERacting via communication
           dis_type = do_thermo_communication
           IF ((region==do_region_global).OR.(region==do_region_molecule)) number = 1
        END IF

        IF (PRESENT(nfree)) THEN
          ! re-initializing simpar%nfree to zero because of multiple thermostats in the adiabatic sampling
           nfree = 0.0_dp
        END IF
    END SELECT

    ! Saving information about thermostats
    thermostat_info%sum_of_thermostats    = sum_of_thermostats
    thermostat_info%number_of_thermostats = number
    thermostat_info%dis_type              = dis_type

    DEALLOCATE(thermolist, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE setup_adiabatic_thermostat_info

! *****************************************************************************
!> \author 10.2011 CJM -PNNL
! *****************************************************************************
  SUBROUTINE get_adiabatic_region_info(region_sections, sum_of_thermostats,&
       thermolist, local_molecules, molecule_kind_set, molecules, particles,&
       qmmm_env, error)
    TYPE(section_vals_type), POINTER         :: region_sections
    INTEGER, INTENT(INOUT), OPTIONAL         :: sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: thermolist( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_adiabatic_region_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER :: first_atom, i, ig, ikind, ilist, imol, ipart, itherm, jg, &
      last_atom, mregions, n_rep, nregions, output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    NULLIFY(tmplist, tmpstringlist, thermolist, molecule_kind, molecule, molecule_set)
    NULLIFY(logger)
    logger      => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
!    CPPostcondition(.NOT.(ASSOCIATED(map_loc_thermo_gen)),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(region_sections, n_repetition=nregions, error=error)
    ALLOCATE(thermolist(particles%n_els), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    thermolist = HUGE(0)
    molecule_set => molecules%els
    mregions = nregions
    itherm = 0
    DO ig = 1, mregions
       CALL section_vals_val_get(region_sections,"LIST", i_rep_section=ig, n_rep_val=n_rep, error=error)
       DO jg = 1, n_rep
          CALL section_vals_val_get(region_sections,"LIST", i_rep_section=ig,i_rep_val=jg,i_vals=tmplist, error=error)
          DO i = 1, SIZE(tmplist)
             ipart = tmplist(i)
             CPPostcondition(((ipart>0).AND.(ipart<=particles%n_els)),cp_failure_level,routineP,error,failure)
             IF (thermolist(ipart)==HUGE(0)) THEN
                itherm = itherm + 1
                thermolist(ipart) = itherm
             ELSE
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
          END DO
       END DO
       CALL section_vals_val_get(region_sections,"MOLNAME", i_rep_section=ig,n_rep_val=n_rep, error=error)
       DO jg = 1,  n_rep
          CALL section_vals_val_get(region_sections,"MOLNAME", i_rep_section=ig,i_rep_val=jg,c_vals=tmpstringlist, error=error)
          DO ilist = 1, SIZE(tmpstringlist)
             DO ikind = 1, SIZE(molecule_kind_set)
                molecule_kind => molecule_kind_set(ikind)
                IF (molecule_kind%name==tmpstringlist(ilist))  THEN
                   DO imol=1,SIZE(molecule_kind%molecule_list)
                      molecule => molecule_set(molecule_kind%molecule_list(imol))
                      CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
                      DO ipart = first_atom, last_atom
                         IF (thermolist(ipart)==HUGE(0)) THEN
                            itherm = itherm + 1
                            thermolist(ipart) = itherm
                         ELSE
                            CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         END IF
                      END DO
                   END DO
                END IF
             END DO
          END DO
       END DO
       CALL setup_thermostat_subsys(region_sections, qmmm_env, thermolist, molecule_set,&
            subsys_qm=.FALSE., ig=ig, sum_of_thermostats=sum_of_thermostats, nregions=nregions,&
            error=error)
       CALL setup_thermostat_subsys(region_sections, qmmm_env, thermolist, molecule_set,&
            subsys_qm=.TRUE., ig=ig, sum_of_thermostats=sum_of_thermostats, nregions=nregions,&
            error=error)
    END DO

    CPPostcondition(.NOT.ALL(thermolist==HUGE(0)),cp_failure_level,routineP,error,failure)

!    natom_local = 0
!    DO ikind = 1, SIZE(molecule_kind_set)
!       nmol_per_kind = local_molecules%n_el(ikind)
!       DO imol = 1, nmol_per_kind
!          i = local_molecules%list(ikind)%array(imol)
!          molecule => molecule_set(i)
!          CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
!          DO ipart = first_atom, last_atom
!             natom_local = natom_local + 1
!          END DO
!       END DO
!    END DO

    ! Now map the local atoms with the corresponding thermostat
!    ALLOCATE(map_loc_thermo_gen(natom_local),stat=stat)
!    map_loc_thermo_gen  = HUGE ( 0 )
!    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!    natom_local = 0
!    DO ikind = 1, SIZE(molecule_kind_set)
!       nmol_per_kind = local_molecules%n_el(ikind)
!       DO imol = 1, nmol_per_kind
!          i = local_molecules%list(ikind)%array(imol)
!          molecule => molecule_set(i)
!          CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
!          DO ipart = first_atom, last_atom
!             natom_local = natom_local + 1
! only map the correct region to the thermostat
!             IF ( thermolist (ipart ) /= HUGE ( 0 ) ) &
!             map_loc_thermo_gen(natom_local) = thermolist(ipart)
!          END DO
!       END DO
!    END DO

!    DEALLOCATE(thermolist, stat=stat)
!    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE get_adiabatic_region_info
! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE setup_thermostat_info(thermostat_info, molecule_kind_set, local_molecules,&
       molecules, particles, region, ensemble, nfree, shell, region_sections, qmmm_env,&
       error)
    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: region, ensemble
    INTEGER, INTENT(INOUT), OPTIONAL         :: nfree
    LOGICAL, INTENT(IN), OPTIONAL            :: shell
    TYPE(section_vals_type), POINTER         :: region_sections
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_thermostat_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dis_type, i, ikind, natom, &
                                                nkind, nmol_local, nmolecule, &
                                                nshell, number, &
                                                sum_of_thermostats
    LOGICAL                                  :: check, do_shell, failure, &
                                                nointer
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    failure = .FALSE.
    NULLIFY(molecule_kind)
    nkind = SIZE ( molecule_kind_set )
    do_shell = .FALSE.
    IF (PRESENT(shell)) do_shell = shell
    ! Counting the global number of thermostats
    sum_of_thermostats = 0
    ! Variable to denote independent thermostats (no communication necessary)
    nointer  = .TRUE.
    check    = .TRUE.
    number   = 0
    dis_type = do_thermo_no_communication

    SELECT CASE(ensemble)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,'Unknown ensemble')
    CASE(isokin_ensemble, nph_uniaxial_ensemble, nph_uniaxial_damped_ensemble,&
         reftraj_ensemble, langevin_ensemble)
       ! Do Nothing
    CASE(nve_ensemble, nvt_ensemble, nvt_adiabatic_ensemble, npt_i_ensemble, &
        npt_f_ensemble, npe_i_ensemble, npe_f_ensemble)
       IF (ensemble==nve_ensemble) check = do_shell
       IF (check) THEN
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
          CASE(do_region_defined)
             ! User defined region to thermostat..
             nointer = .FALSE.
             ! Determine the number of thermostats defined in the input
             CALL section_vals_get(region_sections, n_repetition=sum_of_thermostats, error=error)
             CALL cp_assert((sum_of_thermostats>=1),cp_failure_level,cp_assertion_failed,routineP,&
                  "Provide at least 1 region (&DEFINE_REGION) when using the thermostat type DEFINED"//&
 CPSourceFileRef,&
                  error=error,failure=failure)
          END SELECT

          ! Here we decide which parallel algorithm to use.
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
                   CALL stop_program(routineN,moduleN,__LINE__,'Invalid region setup')
                END IF
             END DO
          ELSE
             ! REPlicated thermostats, INTERacting via communication
             dis_type = do_thermo_communication
             IF     ((region==do_region_global).OR.(region==do_region_molecule)) THEN
                number = 1
             ELSE IF (region==do_region_defined) THEN
                CALL get_defined_region_info(region_sections, number, sum_of_thermostats,&
                     map_loc_thermo_gen=thermostat_info%map_loc_thermo_gen,&
                     local_molecules=local_molecules, molecule_kind_set=molecule_kind_set,&
                     molecules=molecules, particles=particles, qmmm_env=qmmm_env, error=error)
             END IF
          END IF

          IF (PRESENT(nfree)) THEN
             IF ((sum_of_thermostats>1).OR.(dis_type==do_thermo_no_communication)) THEN
                ! re-initializing simpar%nfree to zero because of multiple thermostats
                nfree = 0.0_dp
             END IF
          END IF
       END IF
    END SELECT

    ! Saving information about thermostats
    thermostat_info%sum_of_thermostats    = sum_of_thermostats
    thermostat_info%number_of_thermostats = number
    thermostat_info%dis_type              = dis_type
  END SUBROUTINE setup_thermostat_info

! *****************************************************************************
!> \author 11.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE get_defined_region_info(region_sections, number, sum_of_thermostats,&
       map_loc_thermo_gen, local_molecules, molecule_kind_set, molecules, particles,&
       qmmm_env, error)
    TYPE(section_vals_type), POINTER         :: region_sections
    INTEGER, INTENT(OUT), OPTIONAL           :: number
    INTEGER, INTENT(INOUT), OPTIONAL         :: sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: map_loc_thermo_gen
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set(:)
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_defined_region_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER :: first_atom, i, ig, ikind, ilist, imol, ipart, jg, last_atom, &
      mregions, n_rep, natom_local, nmol_per_kind, nregions, output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: thermolist, tmp, tmplist
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    NULLIFY(tmplist, tmpstringlist, thermolist, molecule_kind, molecule, molecule_set)
    NULLIFY(logger)
    logger      => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    CPPostcondition(.NOT.(ASSOCIATED(map_loc_thermo_gen)),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(region_sections, n_repetition=nregions, error=error)
    ALLOCATE(thermolist(particles%n_els), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    thermolist = HUGE(0)
    molecule_set => molecules%els
    mregions = nregions
    DO ig = 1, mregions
       CALL section_vals_val_get(region_sections,"LIST", i_rep_section=ig, n_rep_val=n_rep, error=error)
       DO jg = 1, n_rep
          CALL section_vals_val_get(region_sections,"LIST", i_rep_section=ig,i_rep_val=jg,i_vals=tmplist, error=error)
          DO i = 1, SIZE(tmplist)
             ipart = tmplist(i)
             CPPostcondition(((ipart>0).AND.(ipart<=particles%n_els)),cp_failure_level,routineP,error,failure)
             IF (thermolist(ipart)==HUGE(0)) THEN
                thermolist(ipart) = ig
             ELSE
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
          END DO
       END DO
       CALL section_vals_val_get(region_sections,"MOLNAME", i_rep_section=ig,n_rep_val=n_rep, error=error)
       DO jg = 1,  n_rep
          CALL section_vals_val_get(region_sections,"MOLNAME", i_rep_section=ig,i_rep_val=jg,c_vals=tmpstringlist, error=error)
          DO ilist = 1, SIZE(tmpstringlist)
             DO ikind = 1, SIZE(molecule_kind_set)
                molecule_kind => molecule_kind_set(ikind)
                IF (molecule_kind%name==tmpstringlist(ilist))  THEN
                   DO imol=1,SIZE(molecule_kind%molecule_list)
                      molecule => molecule_set(molecule_kind%molecule_list(imol))
                      CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
                      DO ipart = first_atom, last_atom
                         IF (thermolist(ipart)==HUGE(0)) THEN
                            thermolist(ipart) = ig
                         ELSE
                            CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         END IF
                      END DO
                   END DO
                END IF
             END DO
          END DO
       END DO
       CALL setup_thermostat_subsys(region_sections, qmmm_env, thermolist, molecule_set,&
            subsys_qm=.FALSE., ig=ig, sum_of_thermostats=sum_of_thermostats, nregions=nregions,&
            error=error)
       CALL setup_thermostat_subsys(region_sections, qmmm_env, thermolist, molecule_set,&
            subsys_qm=.TRUE., ig=ig, sum_of_thermostats=sum_of_thermostats, nregions=nregions,&
            error=error)
    END DO

    ! Dump IO warning for not thermalized particles
    IF (ANY(thermolist==HUGE(0))) THEN
       nregions = nregions + 1
       sum_of_thermostats = sum_of_thermostats + 1
       ALLOCATE(tmp(COUNT(thermolist==HUGE(0))),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ilist = 0
       DO i = 1, SIZE(thermolist)
          IF (thermolist(i)==HUGE(0)) THEN
             ilist = ilist+1
             tmp(ilist) = i
             thermolist(i)=nregions
          END IF
       END DO
       IF (output_unit>0) THEN
          WRITE (output_unit,'(A)')" WARNING| No thermostats defined for the following atoms:"
          IF (ilist>0) WRITE (output_unit,'(" WARNING|",9I8)')tmp
          WRITE (output_unit,'(A)')" WARNING| They will be included in a further unique thermostat!"
       END IF
       DEALLOCATE(tmp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CPPostcondition(ALL(thermolist/=HUGE(0)),cp_failure_level,routineP,error,failure)

    ! Now identify the local number of thermostats
    ALLOCATE(tmp(nregions),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tmp = 0
    natom_local = 0
    DO ikind = 1, SIZE(molecule_kind_set)
       nmol_per_kind = local_molecules%n_el(ikind)
       DO imol = 1, nmol_per_kind
          i = local_molecules%list(ikind)%array(imol)
          molecule => molecule_set(i)
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
          DO ipart = first_atom, last_atom
             natom_local = natom_local + 1
             tmp(thermolist(ipart)) = 1
          END DO
       END DO
    END DO
    number = SUM(tmp)
    DEALLOCATE(tmp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Now map the local atoms with the corresponding thermostat
    ALLOCATE(map_loc_thermo_gen(natom_local),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    natom_local = 0
    DO ikind = 1, SIZE(molecule_kind_set)
       nmol_per_kind = local_molecules%n_el(ikind)
       DO imol = 1, nmol_per_kind
          i = local_molecules%list(ikind)%array(imol)
          molecule => molecule_set(i)
          CALL get_molecule ( molecule, first_atom = first_atom, last_atom = last_atom )
          DO ipart = first_atom, last_atom
             natom_local = natom_local + 1
             map_loc_thermo_gen(natom_local) = thermolist(ipart)
          END DO
       END DO
    END DO

    DEALLOCATE(thermolist, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE get_defined_region_info

! *****************************************************************************
!> \author 11.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE setup_thermostat_subsys(region_sections, qmmm_env, thermolist,&
       molecule_set, subsys_qm, ig, sum_of_thermostats, nregions, error)
    TYPE(section_vals_type), POINTER         :: region_sections
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    INTEGER, DIMENSION(:), POINTER           :: thermolist
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    LOGICAL, INTENT(IN)                      :: subsys_qm
    INTEGER, INTENT(IN)                      :: ig
    INTEGER, INTENT(INOUT)                   :: sum_of_thermostats, nregions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_thermostat_subsys', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: label1, label2
    INTEGER                                  :: first_atom, i, imolecule, &
                                                ipart, last_atom, nrep, &
                                                thermo1
    INTEGER, DIMENSION(:), POINTER           :: atom_index1
    LOGICAL                                  :: failure
    TYPE(molecule_type), POINTER             :: molecule

    failure = .FALSE.
    label1               = "MM_SUBSYS"
    label2               = "QM_SUBSYS"
    IF (subsys_qm) THEN
       label1            = "QM_SUBSYS"
       label2            = "MM_SUBSYS"
    END IF
    CALL section_vals_val_get(region_sections,TRIM(label1),n_rep_val=nrep,i_rep_section=ig,error=error)
    IF (nrep==1) THEN
       IF (ASSOCIATED(qmmm_env)) THEN
          atom_index1          => qmmm_env%mm_atom_index
          IF (subsys_qm) THEN
             atom_index1       => qmmm_env%qm_atom_index
          END IF
          CALL section_vals_val_get(region_sections,TRIM(label1),i_val=thermo1,i_rep_section=ig,error=error)
          SELECT CASE(thermo1)
          CASE(do_constr_atomic)
             DO i = 1, SIZE(atom_index1)
                ipart = atom_index1(i)
                IF (subsys_qm.AND.qmmm_env%qmmm_link.AND.ASSOCIATED(qmmm_env%mm_link_atoms)) THEN
                   IF (ANY(ipart==qmmm_env%mm_link_atoms)) CYCLE
                END IF
                IF (thermolist(ipart)==HUGE(0)) THEN
                   thermolist(ipart) = ig
                ELSE
                   CALL stop_program(routineN,moduleN,__LINE__,&
                        'One atom ('//cp_to_string(ipart)//') of the '//&
                        TRIM(label1)//' was already assigned to'//&
                        ' the thermostatting region Nr.'//cp_to_string(thermolist(ipart))//&
                        '. Please check the input for inconsistencies!')
                END IF
             END DO
          CASE(do_constr_molec)
             DO imolecule=1,SIZE(molecule_set)
                molecule => molecule_set(imolecule)
                CALL get_molecule(molecule,first_atom=first_atom,last_atom=last_atom)
                IF (ANY(atom_index1 >= first_atom .AND. atom_index1 <= last_atom)) THEN
                   DO ipart = first_atom, last_atom
                      IF (thermolist(ipart)==HUGE(0)) THEN
                         thermolist(ipart) = ig
                      ELSE
                         CALL stop_program(routineN,moduleN,__LINE__,&
                              'One atom ('//cp_to_string(ipart)//') of the '//&
                              TRIM(label1)//' was already assigned to'//&
                              ' the thermostatting region Nr.'//cp_to_string(thermolist(ipart))//&
                              '. Please check the input for inconsistencies!')
                      END IF
                   END DO
                END IF
             END DO
          END SELECT
       ELSE
          sum_of_thermostats = sum_of_thermostats - 1
          nregions = nregions - 1
       END IF
    END IF
  END SUBROUTINE setup_thermostat_subsys

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE ke_region_baro(map_info, npt, group, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ke_region_baro', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ncoef
    LOGICAL                                  :: failure

    failure = .FALSE.
    map_info%v_scale = 1.0_dp
    map_info%s_kin   = 0.0_dp
    ncoef = 0
    DO i = 1, SIZE ( npt, 1 )
       DO j = 1, SIZE ( npt, 2 )
          ncoef = ncoef + 1
          map_info%p_kin( 1, ncoef ) % point = map_info%p_kin ( 1, ncoef ) % point &
               + npt ( i, j ) % mass * npt ( i, j )%v**2
       END DO
    END DO

    IF (map_info%dis_type == do_thermo_communication ) CALL mp_sum (map_info%s_kin, group )

  END SUBROUTINE ke_region_baro

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE vel_rescale_baro(map_info, npt, error)
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
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

  END SUBROUTINE vel_rescale_baro

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, vel, error)

    TYPE(map_info_type), POINTER             :: map_info
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    INTEGER, INTENT(IN)                      :: group
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
                IF ( ASSOCIATED ( map_info%p_kin ( 1, ii ) % point ) ) &
                map_info%p_kin(1,ii)%point = map_info%p_kin(1,ii)%point+mass*vel(1,ipart)**2
                IF ( ASSOCIATED ( map_info%p_kin ( 2, ii ) % point ) ) &
                map_info%p_kin(2,ii)%point = map_info%p_kin(2,ii)%point+mass*vel(2,ipart)**2
                IF ( ASSOCIATED ( map_info%p_kin ( 3, ii ) % point ) ) &
                map_info%p_kin(3,ii)%point = map_info%p_kin(3,ii)%point+mass*vel(3,ipart)**2
             ELSE
                IF ( ASSOCIATED ( map_info%p_kin ( 1, ii ) % point ) ) &
                map_info%p_kin(1,ii)%point=map_info%p_kin(1,ii)%point+mass*particle_set(ipart)%v(1)**2
                IF ( ASSOCIATED ( map_info%p_kin ( 2, ii ) % point ) ) &
                map_info%p_kin(2,ii)%point=map_info%p_kin(2,ii)%point+mass*particle_set(ipart)%v(2)**2
                IF ( ASSOCIATED ( map_info%p_kin ( 3, ii ) % point ) ) &
                map_info%p_kin(3,ii)%point=map_info%p_kin(3,ii)%point+mass*particle_set(ipart)%v(3)**2
             END IF
          END DO
       END DO
    END DO

    IF (map_info%dis_type == do_thermo_communication )  CALL mp_sum (map_info%s_kin, group)

  END SUBROUTINE ke_region_particles

! *****************************************************************************
!> \author 07.2009 MI
! *****************************************************************************
  SUBROUTINE momentum_region_particles (map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, vel, error)

    TYPE(map_info_type), POINTER             :: map_info
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'momentum_region_particles', &
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
                map_info%p_kin(1,ii)%point = map_info%p_kin(1,ii)%point+SQRT(mass)*vel(1,ipart)
                map_info%p_kin(2,ii)%point = map_info%p_kin(2,ii)%point+SQRT(mass)*vel(2,ipart)
                map_info%p_kin(3,ii)%point = map_info%p_kin(3,ii)%point+SQRT(mass)*vel(3,ipart)
             ELSE
                map_info%p_kin(1,ii)%point=map_info%p_kin(1,ii)%point+SQRT(mass)*particle_set(ipart)%v(1)
                map_info%p_kin(2,ii)%point=map_info%p_kin(2,ii)%point+SQRT(mass)*particle_set(ipart)%v(2)
                map_info%p_kin(3,ii)%point=map_info%p_kin(3,ii)%point+SQRT(mass)*particle_set(ipart)%v(3)
             END IF
          END DO
       END DO
    END DO

    IF (map_info%dis_type == do_thermo_communication )  CALL mp_sum (map_info%s_kin, group)

  END SUBROUTINE momentum_region_particles

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE vel_rescale_particles(map_info, molecule_kind_set, molecule_set,&
       particle_set, local_molecules, shell_adiabatic, shell_particle_set,&
       core_particle_set, vel, shell_vel, core_vel, error)

    TYPE(map_info_type), POINTER             :: map_info
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
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
                particle_set(ipart)%v(1) =  particle_set(ipart)%v(1) * map_info%p_scale(1,ii)%point
                particle_set(ipart)%v(2) =  particle_set(ipart)%v(2) * map_info%p_scale(2,ii)%point
                particle_set(ipart)%v(3) =  particle_set(ipart)%v(3) * map_info%p_scale(3,ii)%point
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

  END SUBROUTINE vel_rescale_particles

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE ke_region_shells(map_info, particle_set, atomic_kind_set,&
       local_particles, group, core_particle_set, shell_particle_set,&
       core_vel, shell_vel, error)

    TYPE(map_info_type), POINTER             :: map_info
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    INTEGER, INTENT(IN)                      :: group
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
    map_info%s_kin   = 0.0_dp
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

  END SUBROUTINE ke_region_shells

! *****************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE vel_rescale_shells(map_info, atomic_kind_set, particle_set, local_particles, &
       shell_particle_set, core_particle_set, shell_vel, core_vel, vel, error)

    TYPE(map_info_type), POINTER             :: map_info
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
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

  END SUBROUTINE vel_rescale_shells

! *****************************************************************************
!> \brief Calculates kinetic energy and potential energy of the nhc variables
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_nhc_energies ( nhc, nhc_pot, nhc_kin, para_env, array_kin, array_pot, error )
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    REAL(KIND=dp), INTENT(OUT)               :: nhc_pot, nhc_kin
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: array_kin, array_pot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_nhc_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imap, l, n, number, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: akin, vpot

    failure = .FALSE.
    number = nhc%glob_num_nhc
    ALLOCATE(akin(number),stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ALLOCATE(vpot(number),stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    akin = 0.0_dp
    vpot = 0.0_dp
    DO n = 1, nhc%loc_num_nhc
       imap = nhc%map_info%index(n)
       DO l = 1, nhc % nhc_len
          akin(imap) = akin(imap) + 0.5_dp * nhc%nvt(l,n)%mass * nhc%nvt(l,n)%v**2
          vpot(imap) = vpot(imap) + nhc%nvt(l,n)%nkt * nhc%nvt(l,n)%eta
       END DO
    END DO

    ! Handle the thermostat distribution
    IF      (nhc%map_info%dis_type==do_thermo_no_communication) THEN
       CALL mp_sum(akin,para_env%group)
       CALL mp_sum(vpot,para_env%group)
    ELSE IF (nhc%map_info%dis_type==do_thermo_communication) THEN
       CALL communication_thermo_low1(akin, number, para_env, error)
       CALL communication_thermo_low1(vpot, number, para_env, error)
    END IF
    nhc_kin = SUM(akin)
    nhc_pot = SUM(vpot)

    ! Possibly give back kinetic or potential energy arrays
    IF (PRESENT(array_pot)) THEN
       IF (ASSOCIATED(array_pot)) THEN
          CPPrecondition(SIZE(array_pot)==number,cp_fatal_level,routineP,error,failure)
       ELSE
          ALLOCATE(array_pot(number),stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       array_pot=vpot
    END IF
    IF (PRESENT(array_kin)) THEN
       IF (ASSOCIATED(array_kin)) THEN
          CPPrecondition(SIZE(array_kin)==number,cp_fatal_level,routineP,error,failure)
       ELSE
          ALLOCATE(array_kin(number),stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       array_kin=akin
    END IF
    DEALLOCATE(akin,stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    DEALLOCATE(vpot,stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
  END SUBROUTINE get_nhc_energies

! *****************************************************************************
!> \brief Calculates kinetic energy and potential energy
!>      of the csvr  and gle thermostats
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
!> \history generalized MI [07.2009]
! *****************************************************************************
  SUBROUTINE get_kin_energies(map_info,loc_num,glob_num,thermo_energy,thermostat_kin, &
                para_env, array_pot, array_kin, error)

    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, INTENT(IN)                      :: loc_num, glob_num
    REAL(dp), DIMENSION(:), INTENT(IN)       :: thermo_energy
    REAL(KIND=dp), INTENT(OUT)               :: thermostat_kin
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: array_pot, array_kin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kin_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imap, n, number, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: akin

    failure = .FALSE.
    number = glob_num
    ALLOCATE(akin(number),stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    akin     = 0.0_dp
    DO n = 1, loc_num
       imap = map_info%index(n)
       akin(imap) = thermo_energy(n)
    END DO

    ! Handle the thermostat distribution
    IF      (map_info%dis_type==do_thermo_no_communication) THEN
       CALL mp_sum(akin,para_env%group)
    ELSE IF (map_info%dis_type==do_thermo_communication) THEN
       CALL communication_thermo_low1(akin, number, para_env, error)
    END IF
    thermostat_kin = SUM(akin)

    ! Possibly give back kinetic or potential energy arrays
    IF (PRESENT(array_pot)) THEN
       IF (ASSOCIATED(array_pot)) THEN
          CPPrecondition(SIZE(array_pot)==number,cp_fatal_level,routineP,error,failure)
       ELSE
          ALLOCATE(array_pot(number),stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       array_pot=0.0_dp
    END IF
    IF (PRESENT(array_kin)) THEN
       IF (ASSOCIATED(array_kin)) THEN
          CPPrecondition(SIZE(array_kin)==number,cp_fatal_level,routineP,error,failure)
       ELSE
          ALLOCATE(array_kin(number),stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       array_kin=akin
    END IF
    DEALLOCATE(akin,stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
  END SUBROUTINE get_kin_energies

! *****************************************************************************
!> \brief Calculates the temperatures of the regions when a thermostat is
!>        applied
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
!> \history generalized MI [07.2009]
! *****************************************************************************
  SUBROUTINE get_temperatures ( map_info, loc_num, glob_num, nkt, dof, para_env,&
                                temp_tot, array_temp, error )
    TYPE(map_info_type), POINTER             :: map_info
    INTEGER, INTENT(IN)                      :: loc_num, glob_num
    REAL(dp), DIMENSION(:), INTENT(IN)       :: nkt, dof
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(OUT)               :: temp_tot
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: array_temp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_temperatures', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imap, imap2, n, number, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fdeg_of_free
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: akin, deg_of_free

    failure = .FALSE.
    number = glob_num
    ALLOCATE(akin(number),stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ALLOCATE(deg_of_free(number),stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    akin         = 0.0_dp
    deg_of_free  = 0.0_dp
    DO n = 1, loc_num
       imap  = map_info%index(n)
       imap2 = map_info%map_index(n)
       IF (nkt(n)==0.0_dp) CYCLE
       deg_of_free(imap) = REAL(dof(n),KIND=dp)
       akin(imap)        = map_info%s_kin(imap2)
    END DO

    ! Handle the thermostat distribution
    IF      (map_info%dis_type==do_thermo_no_communication) THEN
       CALL mp_sum(akin,para_env%group)
       CALL mp_sum(deg_of_free,para_env%group)
    ELSE IF (map_info%dis_type==do_thermo_communication) THEN
       CALL communication_thermo_low1(akin, number, para_env, error)
       CALL communication_thermo_low1(deg_of_free, number, para_env, error)
    END IF
    temp_tot     = SUM(akin)
    fdeg_of_free = SUM(deg_of_free)

    temp_tot = temp_tot/fdeg_of_free
    temp_tot = cp_unit_from_cp2k(temp_tot, "K_temp", error=error)
    ! Possibly give back temperatures of the full set of regions
    IF (PRESENT(array_temp)) THEN
       IF (ASSOCIATED(array_temp)) THEN
          CPPrecondition(SIZE(array_temp)==number,cp_fatal_level,routineP,error,failure)
       ELSE
          ALLOCATE(array_temp(number),stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       DO i = 1, number
          array_temp(i) = akin(i)/deg_of_free(i)
          array_temp(i) = cp_unit_from_cp2k(array_temp(i), "K_temp", error=error)
       END DO
    END IF
    DEALLOCATE(akin,stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    DEALLOCATE(deg_of_free,stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
  END SUBROUTINE get_temperatures

! *****************************************************************************
!> \brief Calculates energy associated with a thermostat
!> \author Teodoro Laino [tlaino] - 10.2007 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_thermostat_energies ( thermostat, thermostat_pot, thermostat_kin, para_env,&
        array_pot, array_kin, error )
    TYPE(thermostat_type), POINTER           :: thermostat
    REAL(KIND=dp), INTENT(OUT)               :: thermostat_pot, thermostat_kin
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: array_pot, array_kin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_thermostat_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: thermo_energy

    failure = .FALSE.

    thermostat_pot = 0.0_dp
    thermostat_kin = 0.0_dp
    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Energy associated with the Nose-Hoover thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          CALL get_nhc_energies (thermostat%nhc, thermostat_pot, thermostat_kin, para_env,&
                                 array_pot, array_kin, error )
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN
          ! Energy associated with the CSVR thermostat
          CPPostcondition(ASSOCIATED(thermostat%csvr),cp_failure_level,routineP,error,failure)
          ALLOCATE(thermo_energy(thermostat%csvr%loc_num_csvr),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1,thermostat%csvr%loc_num_csvr
             thermo_energy(i) = thermostat%csvr%nvt(i)%thermostat_energy
          END DO
          CALL get_kin_energies(thermostat%csvr%map_info,thermostat%csvr%loc_num_csvr,&
               thermostat%csvr%glob_num_csvr,thermo_energy,&
               thermostat_kin, para_env, array_pot, array_kin, error=error)
          DEALLOCATE(thermo_energy,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ELSE IF (thermostat%type_of_thermostat==do_thermo_gle) THEN
          ! Energy associated with the GLE thermostat
          CPPostcondition(ASSOCIATED(thermostat%gle),cp_failure_level,routineP,error,failure)
          ALLOCATE(thermo_energy(thermostat%gle%loc_num_gle),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1,thermostat%gle%loc_num_gle
             thermo_energy(i) = thermostat%gle%nvt(i)%thermostat_energy
          END DO
          CALL get_kin_energies(thermostat%gle%map_info,thermostat%gle%loc_num_gle,&
               thermostat%gle%glob_num_gle,thermo_energy,&
               thermostat_kin, para_env, array_pot, array_kin, error=error)
          DEALLOCATE(thermo_energy,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ![NB] nothing to do for Ad-Langevin?

       END IF
    END IF

  END SUBROUTINE get_thermostat_energies

! *****************************************************************************
!> \brief Calculates the temperatures for each region associated to a thermostat
!> \author Teodoro Laino [tlaino] - 02.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_region_temperatures ( thermostat, tot_temperature, para_env, array_temp, error )
    TYPE(thermostat_type), POINTER           :: thermostat
    REAL(KIND=dp), INTENT(OUT)               :: tot_temperature
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: array_temp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_region_temperatures', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: dof, nkt

    failure = .FALSE.

    IF (ASSOCIATED(thermostat)) THEN
       IF      (thermostat%type_of_thermostat==do_thermo_nose) THEN
          ! Energy associated with the Nose-Hoover thermostat
          CPPostcondition(ASSOCIATED(thermostat%nhc),cp_failure_level,routineP,error,failure)
          ALLOCATE(nkt(thermostat%nhc%loc_num_nhc), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dof(thermostat%nhc%loc_num_nhc), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1,thermostat%nhc%loc_num_nhc
             nkt(i) = thermostat%nhc%nvt(1,i)%nkt
             dof(i) = REAL(thermostat%nhc%nvt(1,i)%degrees_of_freedom,KIND=dp)
          END DO
          CALL get_temperatures ( thermostat%nhc%map_info, thermostat%nhc%loc_num_nhc, &
          thermostat%nhc%glob_num_nhc, nkt, dof, para_env, tot_temperature, array_temp, error )
          DEALLOCATE(nkt, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dof, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_csvr) THEN
          ! Energy associated with the CSVR thermostat
          CPPostcondition(ASSOCIATED(thermostat%csvr),cp_failure_level,routineP,error,failure)

          ALLOCATE(nkt(thermostat%csvr%loc_num_csvr), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dof(thermostat%csvr%loc_num_csvr), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1,thermostat%csvr%loc_num_csvr
             nkt(i) = thermostat%csvr%nvt(i)%nkt
             dof(i) = REAL(thermostat%csvr%nvt(i)%degrees_of_freedom,KIND=dp)
          END DO
          CALL get_temperatures ( thermostat%csvr%map_info, thermostat%csvr%loc_num_csvr, &
          thermostat%csvr%glob_num_csvr, nkt, dof, para_env, tot_temperature, array_temp, error )
          DEALLOCATE(nkt, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dof, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_al) THEN
          ! Energy associated with the AD_LANGEVIN thermostat
          CPPostcondition(ASSOCIATED(thermostat%al),cp_failure_level,routineP,error,failure)

          ALLOCATE(nkt(thermostat%al%loc_num_al), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dof(thermostat%al%loc_num_al), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1,thermostat%al%loc_num_al
             nkt(i) = thermostat%al%nvt(i)%nkt
             dof(i) = REAL(thermostat%al%nvt(i)%degrees_of_freedom,KIND=dp)
          END DO
          CALL get_temperatures ( thermostat%al%map_info, thermostat%al%loc_num_al, &
          thermostat%al%glob_num_al, nkt, dof, para_env, tot_temperature, array_temp, error )
          DEALLOCATE(nkt, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dof, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ELSE IF (thermostat%type_of_thermostat==do_thermo_gle) THEN
          ! Energy associated with the GLE thermostat
          CPPostcondition(ASSOCIATED(thermostat%gle),cp_failure_level,routineP,error,failure)

          ALLOCATE(nkt(thermostat%gle%loc_num_gle), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dof(thermostat%gle%loc_num_gle), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1,thermostat%gle%loc_num_gle
             nkt(i) = thermostat%gle%nvt(i)%nkt
             dof(i) = REAL(thermostat%gle%nvt(i)%degrees_of_freedom,KIND=dp)
          END DO
          CALL get_temperatures ( thermostat%gle%map_info, thermostat%gle%loc_num_gle, &
          thermostat%gle%glob_num_gle, nkt, dof, para_env, tot_temperature, array_temp, error )
          DEALLOCATE(nkt, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dof, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE get_region_temperatures

! *****************************************************************************
!> \brief Prints status of all thermostats during an MD run
!> \author Teodoro Laino [tlaino] - 02.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE print_thermostats_status(thermostats, para_env, my_pos, my_act, cell, itimes, time, error)
    TYPE(thermostats_type), POINTER          :: thermostats
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=default_string_length)     :: my_pos, my_act
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: itimes
    REAL(KIND=dp), INTENT(IN)                :: time
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_thermostats_status', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(thermostats)) THEN
       IF (ASSOCIATED(thermostats%thermostat_part)) THEN
          CALL print_thermostat_status(thermostats%thermostat_part, para_env, my_pos, my_act, cell, itimes, time, error)
       END IF
       IF (ASSOCIATED(thermostats%thermostat_shell)) THEN
          CALL print_thermostat_status(thermostats%thermostat_shell, para_env, my_pos, my_act, cell, itimes, time, error)
       END IF
       IF (ASSOCIATED(thermostats%thermostat_coef)) THEN
          CALL print_thermostat_status(thermostats%thermostat_coef, para_env, my_pos, my_act, cell, itimes, time, error)
       END IF
       IF (ASSOCIATED(thermostats%thermostat_baro)) THEN
          CALL print_thermostat_status(thermostats%thermostat_baro, para_env, my_pos, my_act, cell, itimes, time, error)
       END IF
    END IF
  END SUBROUTINE print_thermostats_status

! *****************************************************************************
!> \brief Prints status of a specific thermostat
!> \author Teodoro Laino [tlaino] - 02.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE print_thermostat_status ( thermostat, para_env, my_pos, my_act, cell, itimes, time, error )
    TYPE(thermostat_type), POINTER           :: thermostat
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(LEN=default_string_length)     :: my_pos, my_act
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: itimes
    REAL(KIND=dp), INTENT(IN)                :: time
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_thermostat_status', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat, unit
    LOGICAL                                  :: failure, new_file
    REAL(KIND=dp)                            :: thermo_kin, thermo_pot, &
                                                tot_temperature
    REAL(KIND=dp), DIMENSION(:), POINTER     :: array_kin, array_pot, &
                                                array_temp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    failure = .FALSE.
    NULLIFY(logger, print_key, array_pot, array_kin, array_temp)
    logger => cp_error_get_logger(error)

    IF (ASSOCIATED(thermostat)) THEN
       ! Print Energies
       print_key => section_vals_get_subs_vals(thermostat%section,"PRINT%ENERGY",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
          CALL get_thermostat_energies(thermostat, thermo_pot, thermo_kin, para_env, array_pot, array_kin, error)
          unit = cp_print_key_unit_nr(logger,thermostat%section,"PRINT%ENERGY",&
                 extension="."//TRIM(thermostat%label)//".tener",file_position=my_pos, &
                 file_action=my_act, is_new_file=new_file, error=error)
          IF(unit > 0) THEN
            IF (new_file) THEN
               WRITE(unit,'(A)')"# Thermostat Potential and Kinetic Energies - Total and per Region"
               WRITE(unit,'("#",3X,A,2X,A,13X,A,10X,A)')"Step Nr.","Time[fs]","Kin.[a.u.]","Pot.[a.u.]"
            END IF
            WRITE(UNIT=unit,FMT="(I8, F12.3,6X,2F20.10)") itimes,time*femtoseconds, thermo_kin, thermo_pot
            WRITE(unit, '(A,4F20.10)')"# KINETIC ENERGY REGIONS: ",array_kin(1:MIN(4,SIZE(array_kin)))
            DO i = 5, SIZE(array_kin), 4
              WRITE(UNIT=unit,FMT='("#",25X,4F20.10)')array_kin(i:MIN(i+3,SIZE(array_kin)))
            END DO
            WRITE(unit, '(A,4F20.10)')"# POTENT. ENERGY REGIONS: ",array_pot(1:MIN(4,SIZE(array_pot)))
            DO i = 5, SIZE(array_pot), 4
              WRITE(UNIT=unit,FMT='("#",25X,4F20.10)')array_pot(i:MIN(i+3,SIZE(array_pot)))
            END DO
            CALL m_flush(unit)
          END IF
          DEALLOCATE(array_kin, stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          DEALLOCATE(array_pot, stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          CALL cp_print_key_finished_output(unit,logger,thermostat%section,"PRINT%ENERGY", error=error)
       END IF
       ! Print Temperatures of the regions
       print_key => section_vals_get_subs_vals(thermostat%section,"PRINT%TEMPERATURE",error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
          CALL get_region_temperatures(thermostat, tot_temperature, para_env, array_temp, error)
          unit = cp_print_key_unit_nr(logger,thermostat%section,"PRINT%TEMPERATURE",&
                 extension="."//TRIM(thermostat%label)//".temp",file_position=my_pos,&
                 file_action=my_act, is_new_file=new_file, error=error)
          IF(unit > 0) THEN
            IF (new_file) THEN
               WRITE(unit,'(A)')"# Temperature Total and per Region"
               WRITE(unit,'("#",3X,A,2X,A,10X,A)')"Step Nr.","Time[fs]","Temp.[K]"
            END IF
            WRITE(UNIT=unit,FMT="(I8, F12.3,3X,F20.10)") itimes,time*femtoseconds, tot_temperature
            WRITE(unit, '(A,I10)')"# TEMPERATURE REGIONS: ",SIZE(array_temp)
            DO i = 1, SIZE(array_temp), 4
              WRITE(UNIT=unit,FMT='("#",22X,4F20.10)')array_temp(i:MIN(i+3,SIZE(array_temp)))
            END DO
            CALL m_flush(unit)
          END IF
          DEALLOCATE(array_temp, stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          CALL cp_print_key_finished_output(unit,logger,thermostat%section,"PRINT%TEMPERATURE", error=error)
       END IF
    END IF
  END SUBROUTINE print_thermostat_status

! *****************************************************************************
!> \brief Handles the communication for thermostats (1D array)
!> \author Teodoro Laino [tlaino] - University of Zurich 11.2007
! *****************************************************************************
  SUBROUTINE communication_thermo_low1(array, number, para_env, error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: array
    INTEGER, INTENT(IN)                      :: number
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'communication_thermo_low1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, icheck, ncheck, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work, work2

    failure = .FALSE.
    ALLOCATE(work(para_env%num_pe),stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    DO i = 1, number
       work = 0.0_dp
       work(para_env%mepos+1) = array(i)
       CALL mp_sum(work,para_env%group)
       ncheck   = COUNT(work/=0.0_dp)
       array(i) = 0.0_dp
       IF (ncheck/=0) THEN
          ALLOCATE(work2(ncheck),stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          ncheck = 0
          DO icheck = 1, para_env%num_pe
             IF (work(icheck)/=0.0_dp) THEN
                ncheck = ncheck + 1
                work2(ncheck) = work(icheck)
             END IF
          END DO
          CPPrecondition(ncheck==SIZE(work2),cp_fatal_level,routineP,error,failure)
          CPPrecondition(ALL(work2==work2(1)),cp_fatal_level,routineP,error,failure)

          array(i) = work2(1)
          DEALLOCATE(work2,stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
    END DO
    DEALLOCATE(work,stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
  END SUBROUTINE communication_thermo_low1

! *****************************************************************************
!> \brief Handles the communication for thermostats (2D array)
!> \author Teodoro Laino [tlaino] - University of Zurich 11.2007
! *****************************************************************************
  SUBROUTINE communication_thermo_low2(array, number1, number2, para_env, error)
    INTEGER, DIMENSION(:, :), INTENT(INOUT)  :: array
    INTEGER, INTENT(IN)                      :: number1, number2
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'communication_thermo_low2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, icheck, j, ncheck, stat
    INTEGER, DIMENSION(:, :), POINTER        :: work, work2
    LOGICAL                                  :: failure

    failure = .FALSE.
    ALLOCATE(work(number1,para_env%num_pe),stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    DO i = 1, number2
       work = 0
       work(:,para_env%mepos+1) = array(:,i)
       CALL mp_sum(work,para_env%group)
       ncheck = 0
       DO j = 1, para_env%num_pe
          IF (ANY(work(:,j)/=0)) THEN
             ncheck   = ncheck + 1
          END IF
       END DO
       array(:,i) = 0
       IF (ncheck/=0) THEN
          ALLOCATE(work2(number1,ncheck),stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          ncheck = 0
          DO icheck = 1, para_env%num_pe
             IF (ANY(work(:,icheck)/=0)) THEN
                ncheck = ncheck + 1
                work2(:,ncheck) = work(:,icheck)
             END IF
          END DO
          CPPrecondition(ncheck==SIZE(work2,2),cp_fatal_level,routineP,error,failure)
          DO j = 1, ncheck
             CPPrecondition(ALL(work2(:,j)==work2(:,1)),cp_fatal_level,routineP,error,failure)
          END DO
          array(:,i) = work2(:,1)
          DEALLOCATE(work2,stat=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
    END DO
    DEALLOCATE(work,stat=stat)
    CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
  END SUBROUTINE communication_thermo_low2

END MODULE thermostat_utils
