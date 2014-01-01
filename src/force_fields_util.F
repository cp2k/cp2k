!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Subroutine input_torsions changed (DG) 05-Dec-2000
!>      Output formats changed (DG) 05-Dec-2000
!>      JGH (26-01-2002) : force field parameters stored in tables, not in
!>        matrices. Input changed to have parameters labeled by the position
!>        and not atom pairs (triples etc)
!>      Teo (11.2005) : Moved all information on force field  pair_potential to
!>                      a much lighter memory structure
!> \author CJM
! *****************************************************************************
MODULE force_fields_util
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE colvar_types,                    ONLY: dist_colvar_id
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE f77_blas
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_create,&
                                             fist_nonbond_env_set,&
                                             fist_nonbond_env_type
  USE force_field_types,               ONLY: &
       allocate_bend_kind_set, allocate_bond_kind_set, &
       allocate_impr_kind_set, allocate_opbend_kind_set, &
       allocate_torsion_kind_set, allocate_ub_kind_set, amber_info_type, &
       bend_kind_type, bond_kind_type, charmm_info_type, &
       deallocate_bend_kind_set, deallocate_bond_kind_set, force_field_type, &
       gromos_info_type, impr_kind_dealloc_ref, impr_kind_type, &
       input_info_type, opbend_kind_type, torsion_kind_dealloc_ref, &
       torsion_kind_type, ub_kind_dealloc_ref, ub_kind_type
  USE force_fields_all,                ONLY: &
       force_field_pack_bend, force_field_pack_bond, force_field_pack_charge, &
       force_field_pack_charges, force_field_pack_damp, &
       force_field_pack_eicut, force_field_pack_impr, &
       force_field_pack_nonbond, force_field_pack_nonbond14, &
       force_field_pack_opbend, force_field_pack_pol, &
       force_field_pack_radius, force_field_pack_shell, &
       force_field_pack_splines, force_field_pack_tors, force_field_pack_ub, &
       force_field_unique_bend, force_field_unique_bond, &
       force_field_unique_impr, force_field_unique_opbend, &
       force_field_unique_tors, force_field_unique_ub
  USE input_constants,                 ONLY: do_ff_undef
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: &
       atom_type, bend_type, bond_type, colvar_constraint_type, &
       g3x3_constraint_type, g4x6_constraint_type, get_molecule_kind, &
       impr_type, molecule_kind_type, opbend_type, set_molecule_kind, &
       torsion_type, ub_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE pair_potential_types,            ONLY: pair_potential_pp_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields_util'

  PRIVATE
  PUBLIC :: force_field_pack,&
            force_field_qeff_output,&
            clean_intra_force_kind,&
            get_generic_info

CONTAINS

! *****************************************************************************
!> \brief Pack in all the information needed for the force_field
! *****************************************************************************
  SUBROUTINE force_field_pack (particle_set,atomic_kind_set,molecule_kind_set, &
       molecule_set, ewald_env,fist_nonbond_env,ff_type,root_section, qmmm, &
       qmmm_env, mm_section, subsys_section, shell_particle_set,core_particle_set,&
       cell, error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(section_vals_type), POINTER         :: root_section
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: mm_section, subsys_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_particle_set, &
                                                core_particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: Ainfo
    INTEGER                                  :: handle, iw, iw2, iw3, iw4, &
                                                output_unit
    LOGICAL                                  :: do_zbl, explicit, failure, &
                                                fatal, ignore_fatal, my_qmmm
    REAL(KIND=dp)                            :: ewald_rcut, verlet_skin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    TYPE(amber_info_type), POINTER           :: amb_info
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gromos_info_type), POINTER          :: gro_info
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(pair_potential_pp_type), POINTER    :: potparm_nonbond, &
                                                potparm_nonbond14
    TYPE(section_vals_type), POINTER         :: charges_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    fatal   = .FALSE.
    ignore_fatal = ff_type%ignore_missing_critical
    NULLIFY(logger, Ainfo, charges_section, charges)
    logger => cp_error_get_logger(error)
    ! Error unit
    output_unit= cp_logger_get_default_io_unit(logger)

    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)
    iw2= cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO/SPLINE_INFO",&
         extension=".mmLog",error=error)
    iw3= cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO/SPLINE_DATA",&
         extension=".mmLog",error=error)
    iw4= cp_print_key_unit_nr(logger,mm_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".mmLog",error=error)
    NULLIFY(potparm_nonbond14, potparm_nonbond)
    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm
    inp_info => ff_type%inp_info
    chm_info => ff_type%chm_info
    gro_info => ff_type%gro_info
    amb_info => ff_type%amb_info
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Determine the number of unique bond kind and allocate bond_kind_set
    !-----------------------------------------------------------------------------
    CALL  force_field_unique_bond (particle_set, molecule_kind_set, molecule_set, &
         ff_type, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Determine the number of unique bend kind and allocate bend_kind_set
    !-----------------------------------------------------------------------------
    CALL force_field_unique_bend (particle_set, molecule_kind_set, molecule_set, &
         ff_type, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Determine the number of unique Urey-Bradley kind and allocate ub_kind_set
    !-----------------------------------------------------------------------------
    CALL force_field_unique_ub (particle_set, molecule_kind_set, molecule_set,   &
         error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Determine the number of unique torsion kind and allocate torsion_kind_set
    !-----------------------------------------------------------------------------
    CALL force_field_unique_tors (particle_set, molecule_kind_set, molecule_set, &
         ff_type, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Determine the number of unique impr kind and allocate impr_kind_set
    !-----------------------------------------------------------------------------
    CALL force_field_unique_impr (particle_set, molecule_kind_set, molecule_set, &
         ff_type, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Determine the number of unique opbend kind and allocate opbend_kind_set
    !-----------------------------------------------------------------------------
    CALL force_field_unique_opbend (particle_set, molecule_kind_set, molecule_set, &
         ff_type, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Bonds
    !-----------------------------------------------------------------------------
    CALL force_field_pack_bond (particle_set, molecule_kind_set, molecule_set,   &
         fatal, Ainfo, chm_info, inp_info, gro_info, amb_info, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Bends
    !-----------------------------------------------------------------------------
    CALL force_field_pack_bend (particle_set, molecule_kind_set, molecule_set,   &
         fatal, Ainfo, chm_info, inp_info, gro_info, amb_info, error)
    ! Give information and abort if any bond or angle parameter is missing..
    CALL release_FF_missing_par(fatal,ignore_fatal,AInfo,output_unit,iw,error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 9. Urey-Bradley
    !-----------------------------------------------------------------------------
    CALL force_field_pack_ub (particle_set, molecule_kind_set, molecule_set,     &
         Ainfo, chm_info, inp_info, iw, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. Torsions
    !-----------------------------------------------------------------------------
    CALL force_field_pack_tors (particle_set, molecule_kind_set, molecule_set,   &
         Ainfo, chm_info, inp_info, gro_info, amb_info, iw, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Impropers
    !-----------------------------------------------------------------------------
    CALL force_field_pack_impr (particle_set, molecule_kind_set, molecule_set,   &
         Ainfo, chm_info, inp_info, gro_info, amb_info, error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 12. Out of plane bends
    !-----------------------------------------------------------------------------
    CALL force_field_pack_opbend (particle_set, molecule_kind_set, molecule_set, &
       Ainfo, inp_info, error)
    ! Give information only if any Urey-Bradley, Torsion, improper or opbend is missing
    ! continue calculation..
    CALL release_FF_missing_par(fatal,ignore_fatal,AInfo,output_unit,iw,error)

    charges_section => section_vals_get_subs_vals(mm_section,"FORCEFIELD%CHARGES",error=error)
    CALL section_vals_get(charges_section, explicit=explicit, error=error)
    IF (.NOT.explicit) THEN
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 13.a Set up atomic_kind_set()%fist_potentail%[qeff] and shell
       !      potential parameters
       !-----------------------------------------------------------------------------
       CALL force_field_pack_charge(atomic_kind_set, qmmm_env, fatal, iw, iw4,&
            Ainfo, my_qmmm, inp_info, error)
       ! Give information only if charge is missing and abort..
       CALL release_FF_missing_par(fatal,ignore_fatal,AInfo,output_unit,iw,error)
    ELSE
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 13.b Setup a global array of classical charges - this avoids the packing and
       !      allows the usage of different charges for same atomic types
       !-----------------------------------------------------------------------------
       CALL force_field_pack_charges(charges, charges_section, particle_set, my_qmmm,&
            qmmm_env, inp_info, iw4, error)
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 14. Set up the radius of the electrostatic multipole in Fist
    !-----------------------------------------------------------------------------
    CALL force_field_pack_radius(atomic_kind_set, iw, subsys_section, error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 15. Set up the polarizable FF parameters
    !-----------------------------------------------------------------------------
    CALL force_field_pack_pol(atomic_kind_set, iw, inp_info, error)
    CALL force_field_pack_damp(atomic_kind_set, iw, inp_info, error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 16. Set up  Shell potential
    !-----------------------------------------------------------------------------
    CALL force_field_pack_shell ( particle_set, atomic_kind_set,&
         molecule_kind_set, molecule_set, root_section, subsys_section,&
         shell_particle_set, core_particle_set, cell, iw, inp_info, error)
    IF (ff_type%do_nonbonded) THEN
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 17. Set up potparm_nonbond14
       !-----------------------------------------------------------------------------
       ! Move the data from the info structures to potparm_nonbond
       CALL force_field_pack_nonbond14(atomic_kind_set, ff_type, qmmm_env, iw, Ainfo,&
            chm_info, inp_info, gro_info, amb_info, potparm_nonbond14, ewald_env, &
            iw2, iw3, error)
       ! Give information if any 1-4 is missing.. continue calculation..
       CALL release_FF_missing_par(fatal,ignore_fatal,AInfo,output_unit,iw,error)
       ! Create the spline data
       CALL section_vals_val_get(mm_section,"FORCEFIELD%ZBL_SCATTERING",l_val=do_zbl,error=error)
       CALL force_field_pack_splines(atomic_kind_set, ff_type, iw2, iw3, iw4, &
            potparm_nonbond14, do_zbl, nonbonded_type="NONBONDED14", error=error)
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 18. Set up potparm_nonbond
       !-----------------------------------------------------------------------------
       ! Move the data from the info structures to potparm_nonbond
       CALL force_field_pack_nonbond(atomic_kind_set, ff_type, qmmm_env, &
            fatal, iw, Ainfo, chm_info, inp_info, gro_info, amb_info, &
            potparm_nonbond, ewald_env, error)
       ! Give information and abort if any pair potential spline is missing..
       CALL release_FF_missing_par(fatal,ignore_fatal,AInfo,output_unit,iw,error)
       ! Create the spline data
       CALL section_vals_val_get(mm_section,"FORCEFIELD%ZBL_SCATTERING",l_val=do_zbl,error=error)
       CALL force_field_pack_splines(atomic_kind_set, ff_type, iw2, iw3, iw4, &
            potparm_nonbond, do_zbl, nonbonded_type="NONBONDED", error=error)
    END IF
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 19. Create nonbond environment
    !-----------------------------------------------------------------------------
    CALL ewald_env_get(ewald_env, rcut=ewald_rcut, error=error)
    CALL section_vals_val_get(mm_section,"NEIGHBOR_LISTS%VERLET_SKIN",&
         r_val=verlet_skin,error=error)
    CALL fist_nonbond_env_create (fist_nonbond_env, atomic_kind_set, &
         potparm_nonbond14, potparm_nonbond, ff_type%do_nonbonded, &
         verlet_skin, ewald_rcut, ff_type%ei_scale14, &
         ff_type%vdw_scale14, ff_type%shift_cutoff, error)
    CALL fist_nonbond_env_set(fist_nonbond_env, charges=charges, error=error)
    ! Compute the electrostatic interaction cutoffs.
    CALL force_field_pack_eicut(atomic_kind_set, ff_type, potparm_nonbond, &
         ewald_env, error)

    CALL cp_print_key_finished_output(iw4,logger,mm_section,"PRINT%PROGRAM_RUN_INFO",error=error)
    CALL cp_print_key_finished_output(iw3,logger,mm_section,"PRINT%FF_INFO/SPLINE_DATA",error=error)
    CALL cp_print_key_finished_output(iw2,logger,mm_section,"PRINT%FF_INFO/SPLINE_INFO",error=error)
    CALL cp_print_key_finished_output(iw,logger,mm_section,"PRINT%FF_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE force_field_pack

! *****************************************************************************
!> \brief Store informations on possible missing ForceFields parameters
! *****************************************************************************
  SUBROUTINE release_FF_missing_par(fatal, ignore_fatal, array, output_unit, iw, error)
    LOGICAL, INTENT(INOUT)                   :: fatal, ignore_fatal
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: array
    INTEGER, INTENT(IN)                      :: output_unit, iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_FF_missing_par', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(array)) THEN
       IF (output_unit>0) THEN
          WRITE(output_unit,*)
          WRITE(output_unit,'(T2,"FORCEFIELD|",A)')&
               " WARNING: A non Critical ForceField parameter is missing! CP2K GOES!",&
               " Non critical parameters are:Urey-Bradley,Torsions,Impropers, Opbends and 1-4",&
               " All missing parameters will not contribute to the potential energy!"
          IF (fatal.OR.iw>0) THEN
             WRITE(output_unit,*)
             DO i = 1,SIZE(array)
                WRITE(output_unit,'(A)')array(i)
             END DO
          END IF
          IF (.NOT.fatal.AND.iw<0) THEN
             WRITE(output_unit,'(T2,"FORCEFIELD|",A)')&
                  " Activate the print key FF_INFO to have a list of missing parameters"
             WRITE(output_unit,*)
          END IF
       END IF
       DEALLOCATE(array,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (fatal) THEN
      IF (ignore_fatal) THEN
        IF (output_unit>0) THEN
          WRITE(output_unit,*)
          WRITE(output_unit,'(T2,"FORCEFIELD|",A)')&
               " WARNING: Ignoring missing critical FF parameters! CP2K GOES!",&
               " Critical parameters are: Bonds, Bends and Charges",&
               " All missing parameters will not contribute to the potential energy!"
        END IF
      ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Missing critical ForceField parameters!")
      END IF
    END IF
  END SUBROUTINE release_FF_missing_par

! *****************************************************************************
!> \brief Compute the total qeff charges for each molecule kind and total system
! *****************************************************************************
  SUBROUTINE force_field_qeff_output (particle_set,molecule_kind_set,&
       molecule_set,mm_section,charges,error)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(section_vals_type), POINTER         :: mm_section
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_qeff_output', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atmname, molname
    INTEGER                                  :: first, handle, iatom, imol, &
                                                iw, j, jatom
    LOGICAL                                  :: failure, shell_active
    REAL(KIND=dp)                            :: mass, mass_mol, mass_sum, &
                                                qeff, qeff_mol, qeff_sum
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(shell_kind_type), POINTER           :: shell

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    qeff     = 0.0_dp
    qeff_mol = 0.0_dp
    qeff_sum = 0.0_dp
    mass_sum = 0.0_dp

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Sum of [qeff,mass] for each molecule_kind
    !-----------------------------------------------------------------------------
    DO imol=1,SIZE(molecule_kind_set)
       qeff_mol=0.0_dp
       mass_mol=0.0_dp
       molecule_kind => molecule_kind_set(imol)

       j=molecule_kind_set(imol)%molecule_list(1)
       molecule => molecule_set(j)
       CALL get_molecule(molecule=molecule, first_atom=first)

       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            name=molname,atom_list=atom_list)
       DO iatom=1,SIZE(atom_list)
          atomic_kind => atom_list(iatom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               name=atmname,qeff=qeff,mass=mass, shell_active=shell_active, shell=shell)
          IF(shell_active) qeff = shell%charge_core+shell%charge_shell
          IF (ASSOCIATED(charges)) THEN
             jatom=first-1+iatom
             qeff=charges(jatom)
          END IF
          IF(iw>0) WRITE(iw,*) "      atom ",iatom," ",TRIM(atmname)," charge = ",qeff
          qeff_mol = qeff_mol + qeff
          mass_mol = mass_mol + mass
       END DO
       CALL set_molecule_kind(molecule_kind=molecule_kind,charge=qeff_mol,mass=mass_mol)
       IF(iw>0) WRITE(iw,*) "    Mol Kind ",TRIM(molname)," charge = ",qeff_mol
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Sum of [qeff,mass] for particle_set
    !-----------------------------------------------------------------------------
    DO iatom=1,SIZE(particle_set)
       atomic_kind => particle_set(iatom)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
            name=atmname,qeff=qeff,mass=mass,shell_active=shell_active, shell=shell)
       IF(shell_active) qeff = shell%charge_core+shell%charge_shell
       IF (ASSOCIATED(charges)) THEN
          qeff=charges(iatom)
       END IF
       IF(iw>0) WRITE(iw,*) "      atom ",iatom," ",TRIM(atmname),&
            " charge = ",qeff
       qeff_sum = qeff_sum + qeff
       mass_sum = mass_sum + mass
    END DO
    IF(iw>0) WRITE(iw,'(A,F20.10)') "    Total system charge = ",qeff_sum
    IF(iw>0) WRITE(iw,'(A,F20.10)') "    Total system mass   = ",mass_sum

    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(handle)
  END SUBROUTINE force_field_qeff_output

! *****************************************************************************
!> \brief Removes UNSET force field types
! *****************************************************************************
  SUBROUTINE clean_intra_force_kind (molecule_kind_set,mm_section,error)

    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'clean_intra_force_kind', &
      routineP = moduleN//':'//routineN

    INTEGER :: atm2_a, atm2_b, atm2_c, atm_a, atm_b, atm_c, atm_d, counter, &
      handle, i, ibend, ibond, icolv, ig3x3, ig4x6, iimpr, ikind, iopbend, &
      itorsion, iub, iw, j, k, nbend, nbond, newkind, ng3x3, ng4x6, nimpr, &
      nopbend, ntorsion, nub, stat
    INTEGER, POINTER                         :: bad1(:), bad2(:)
    LOGICAL                                  :: failure, unsetme, valid_kind
    TYPE(bend_kind_type), DIMENSION(:), &
      POINTER                                :: bend_kind_set, &
                                                new_bend_kind_set
    TYPE(bend_type), DIMENSION(:), POINTER   :: bend_list, new_bend_list
    TYPE(bond_kind_type), DIMENSION(:), &
      POINTER                                :: bond_kind_set, &
                                                new_bond_kind_set
    TYPE(bond_type), DIMENSION(:), POINTER   :: bond_list, new_bond_list
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(g3x3_constraint_type), &
      DIMENSION(:), POINTER                  :: g3x3_list
    TYPE(g4x6_constraint_type), &
      DIMENSION(:), POINTER                  :: g4x6_list
    TYPE(impr_kind_type), DIMENSION(:), &
      POINTER                                :: impr_kind_set, &
                                                new_impr_kind_set
    TYPE(impr_type), DIMENSION(:), POINTER   :: impr_list, new_impr_list
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(opbend_kind_type), DIMENSION(:), &
      POINTER                                :: new_opbend_kind_set, &
                                                opbend_kind_set
    TYPE(opbend_type), DIMENSION(:), POINTER :: new_opbend_list, opbend_list
    TYPE(torsion_kind_type), DIMENSION(:), &
      POINTER                                :: new_torsion_kind_set, &
                                                torsion_kind_set
    TYPE(torsion_type), DIMENSION(:), &
      POINTER                                :: new_torsion_list, torsion_list
    TYPE(ub_kind_type), DIMENSION(:), &
      POINTER                                :: new_ub_kind_set, ub_kind_set
    TYPE(ub_type), DIMENSION(:), POINTER     :: new_ub_list, ub_list

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Lets Tag the unwanted bonds due to the use of distance constraint
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            colv_list=colv_list,&
            nbond=nbond,&
            bond_list=bond_list)
       IF (ASSOCIATED(colv_list)) THEN
          DO icolv=1,SIZE(colv_list)
             IF ((colv_list(icolv)%type_id == dist_colvar_id).AND.&
                  ((.NOT.colv_list(icolv)%use_points).OR.(SIZE(colv_list(icolv)%i_atoms)==2))) THEN
                atm_a = colv_list(icolv)%i_atoms(1)
                atm_b = colv_list(icolv)%i_atoms(2)
                DO ibond=1,nbond
                   unsetme = .FALSE.
                   atm2_a = bond_list(ibond)%a
                   atm2_b = bond_list(ibond)%b
                   IF(atm2_a==atm_a .AND. atm2_b==atm_b) unsetme=.TRUE.
                   IF(atm2_a==atm_b .AND. atm2_b==atm_a) unsetme=.TRUE.
                   IF(unsetme) bond_list(ibond)%id_type = do_ff_undef
                END DO
             END IF
          END DO
       END IF
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Lets Tag the unwanted bends due to the use of distance constraint
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            colv_list=colv_list,&
            nbend=nbend,&
            bend_list=bend_list)
       IF (ASSOCIATED(colv_list)) THEN
          DO ibend=1,nbend
             unsetme = .FALSE.
             atm_a = bend_list(ibend)%a
             atm_b = bend_list(ibend)%b
             atm_c = bend_list(ibend)%c
             DO icolv=1,SIZE(colv_list)
                IF ((colv_list(icolv)%type_id == dist_colvar_id).AND.&
                     ((.NOT.colv_list(icolv)%use_points).OR.(SIZE(colv_list(icolv)%i_atoms)==2))) THEN
                   atm2_a = colv_list(icolv)%i_atoms(1)
                   atm2_b = colv_list(icolv)%i_atoms(2)
                   ! Check that the bonds we're constraining does not involve atoms defining
                   ! a bend..
                   IF  ((( atm2_a == atm_a ).AND.( atm2_b == atm_c )).OR.&
                        (( atm2_a == atm_c ).AND.( atm2_b == atm_a ))) THEN
                      unsetme=.TRUE.
                      EXIT
                   END IF
                END IF
             END DO
             IF(unsetme) bend_list(ibend)%id_type = do_ff_undef
          END DO
       END IF
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Lets Tag the unwanted bonds due to the use of 3x3
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng3x3=ng3x3,&
            g3x3_list=g3x3_list,&
            nbond=nbond,&
            bond_list=bond_list)
       DO ig3x3=1,ng3x3
          atm_a = g3x3_list(ig3x3)%a
          atm_b = g3x3_list(ig3x3)%b
          atm_c = g3x3_list(ig3x3)%c
          DO ibond=1,nbond
             unsetme = .FALSE.
             atm2_a = bond_list(ibond)%a
             atm2_b = bond_list(ibond)%b
             IF(atm2_a==atm_a .AND. atm2_b==atm_b) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_c .AND. atm2_b==atm_c) unsetme=.TRUE.
             IF(unsetme) bond_list(ibond)%id_type = do_ff_undef
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Lets Tag the unwanted bends due to the use of 3x3
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng3x3=ng3x3,&
            g3x3_list=g3x3_list,&
            nbend=nbend,&
            bend_list=bend_list)
       DO ig3x3=1,ng3x3
          atm_a = g3x3_list(ig3x3)%a
          atm_b = g3x3_list(ig3x3)%b
          atm_c = g3x3_list(ig3x3)%c
          DO ibend=1,nbend
             unsetme = .FALSE.
             atm2_a = bend_list(ibend)%a
             atm2_b = bend_list(ibend)%b
             atm2_c = bend_list(ibend)%c
             IF(atm2_a==atm_a .AND. atm2_b==atm_b .AND. atm2_c==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_c .AND. atm2_c==atm_b) unsetme=.TRUE.
             IF(atm2_a==atm_b .AND. atm2_b==atm_a .AND. atm2_c==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_b .AND. atm2_b==atm_c .AND. atm2_c==atm_a) unsetme=.TRUE.
             IF(unsetme) bend_list(ibend)%id_type = do_ff_undef
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Lets Tag the unwanted bonds due to the use of 4x6
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng4x6=ng4x6,&
            g4x6_list=g4x6_list,&
            nbond=nbond,&
            bond_list=bond_list)
       DO ig4x6=1,ng4x6
          atm_a = g4x6_list(ig4x6)%a
          atm_b = g4x6_list(ig4x6)%b
          atm_c = g4x6_list(ig4x6)%c
          atm_d = g4x6_list(ig4x6)%d
          DO ibond=1,nbond
             unsetme = .FALSE.
             atm2_a = bond_list(ibond)%a
             atm2_b = bond_list(ibond)%b
             IF(atm2_a==atm_a .AND. atm2_b==atm_b) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_a .AND. atm2_b==atm_d) unsetme=.TRUE.
             IF(unsetme) bond_list(ibond)%id_type = do_ff_undef
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Lets Tag the unwanted bends due to the use of 4x6
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ng4x6=ng4x6,&
            g4x6_list=g4x6_list,&
            nbend=nbend,&
            bend_list=bend_list)
       DO ig4x6=1,ng4x6
          atm_a = g4x6_list(ig4x6)%a
          atm_b = g4x6_list(ig4x6)%b
          atm_c = g4x6_list(ig4x6)%c
          atm_d = g4x6_list(ig4x6)%d
          DO ibend=1,nbend
             unsetme = .FALSE.
             atm2_a = bend_list(ibend)%a
             atm2_b = bend_list(ibend)%b
             atm2_c = bend_list(ibend)%c
             IF(atm2_a==atm_b .AND. atm2_b==atm_a .AND. atm2_c==atm_c) unsetme=.TRUE.
             IF(atm2_a==atm_b .AND. atm2_b==atm_a .AND. atm2_c==atm_d) unsetme=.TRUE.
             IF(atm2_a==atm_c .AND. atm2_b==atm_a .AND. atm2_c==atm_d) unsetme=.TRUE.
             IF(unsetme) bend_list(ibend)%id_type = do_ff_undef
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Count the number of UNSET bond kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nbond=nbond,&
            bond_kind_set=bond_kind_set,&
            bond_list=bond_list)
       IF(nbond>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old BOND Count: ",&
               SIZE(bond_list),SIZE(bond_kind_set)
          IF (iw>0) WRITE(iw,'(2I6)')(bond_list(ibond)%a,bond_list(ibond)%b,ibond=1,SIZE(bond_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(bond_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO ibond=1,SIZE(bond_kind_set)
             unsetme=.FALSE.
             IF(bond_kind_set(ibond)%id_type==do_ff_undef) unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(bond_list)
                IF(bond_list(i)%id_type/=do_ff_undef.AND.&
                     bond_list(i)%bond_kind%kind_number==ibond) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(ibond) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(bond_kind_set)-SUM(bad1)
             CALL allocate_bond_kind_set(new_bond_kind_set,counter,error)
             counter=0
             DO ibond=1,SIZE(bond_kind_set)
                IF(bad1(ibond)==0) THEN
                   counter=counter+1
                   new_bond_kind_set(counter)= bond_kind_set(ibond)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(bond_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO ibond=1,SIZE(bond_list)
                unsetme = .FALSE.
                IF(bond_list(ibond)%bond_kind%id_type==do_ff_undef) unsetme = .TRUE.
                IF(bond_list(ibond)%id_type==do_ff_undef)           unsetme = .TRUE.
                IF(unsetme) bad2(ibond) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(bond_list)-SUM(bad2)
                ALLOCATE(new_bond_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO ibond=1,SIZE(bond_list)
                   IF(bad2(ibond)==0) THEN
                      counter=counter+1
                      new_bond_list(counter) = bond_list(ibond)
                      newkind = bond_list(ibond)%bond_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_bond_list(counter)%bond_kind => new_bond_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nbond=SIZE(new_bond_list),&
                     bond_kind_set=new_bond_kind_set,&
                     bond_list=new_bond_list)
                DO ibond=1,SIZE(new_bond_kind_set)
                   new_bond_kind_set(ibond)%kind_number=ibond
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL deallocate_bond_kind_set(bond_kind_set,error)
             DEALLOCATE(bond_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New BOND Count: ",&
                  SIZE(new_bond_list),SIZE(new_bond_kind_set)
             IF (iw>0) WRITE(iw,'(2I6)')(new_bond_list(ibond)%a,new_bond_list(ibond)%b,&
                  ibond=1,SIZE(new_bond_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Count the number of UNSET bend kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nbend=nbend,&
            bend_kind_set=bend_kind_set,&
            bend_list=bend_list)
       IF(nbend>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old BEND Count: ",&
               SIZE(bend_list),SIZE(bend_kind_set)
          IF (iw>0) WRITE(iw,'(3I6)')(bend_list(ibend)%a,bend_list(ibend)%b,&
               bend_list(ibend)%c,ibend=1,SIZE(bend_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(bend_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO ibend=1,SIZE(bend_kind_set)
             unsetme=.FALSE.
             IF(bend_kind_set(ibend)%id_type==do_ff_undef) unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(bend_list)
                IF(bend_list(i)%id_type/=do_ff_undef.AND.&
                     bend_list(i)%bend_kind%kind_number==ibend) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(ibend) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(bend_kind_set)-SUM(bad1)
             CALL allocate_bend_kind_set(new_bend_kind_set,counter,error)
             counter=0
             DO ibend=1,SIZE(bend_kind_set)
                IF(bad1(ibend)==0) THEN
                   counter=counter+1
                   new_bend_kind_set(counter)= bend_kind_set(ibend)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(bend_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO ibend=1,SIZE(bend_list)
                unsetme = .FALSE.
                IF(bend_list(ibend)%bend_kind%id_type==do_ff_undef) unsetme = .TRUE.
                IF(bend_list(ibend)%id_type==do_ff_undef)           unsetme = .TRUE.
                IF(unsetme) bad2(ibend) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(bend_list)-SUM(bad2)
                ALLOCATE(new_bend_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO ibend=1,SIZE(bend_list)
                   IF(bad2(ibend)==0) THEN
                      counter=counter+1
                      new_bend_list(counter) = bend_list(ibend)
                      newkind = bend_list(ibend)%bend_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_bend_list(counter)%bend_kind => new_bend_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nbend=SIZE(new_bend_list),&
                     bend_kind_set=new_bend_kind_set,&
                     bend_list=new_bend_list)
                DO ibend=1,SIZE(new_bend_kind_set)
                   new_bend_kind_set(ibend)%kind_number=ibend
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL deallocate_bend_kind_set(bend_kind_set,error)
             DEALLOCATE(bend_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New BEND Count: ",&
                  SIZE(new_bend_list),SIZE(new_bend_kind_set)
             IF (iw>0) WRITE(iw,'(3I6)')(new_bend_list(ibend)%a,new_bend_list(ibend)%b,&
                  new_bend_list(ibend)%c,ibend=1,SIZE(new_bend_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 9. Count the number of UNSET Urey-Bradley kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nub=nub,&
            ub_kind_set=ub_kind_set,&
            ub_list=ub_list)
       IF(nub>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old UB Count: ",&
               SIZE(ub_list),SIZE(ub_kind_set)
          IF (iw>0) WRITE(iw,'(3I6)')(ub_list(iub)%a,ub_list(iub)%b,&
               ub_list(iub)%c,iub=1,SIZE(ub_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(ub_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO iub=1,SIZE(ub_kind_set)
             unsetme=.FALSE.
             IF(ub_kind_set(iub)%id_type==do_ff_undef) unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(ub_list)
                IF(ub_list(i)%id_type/=do_ff_undef.AND.&
                     ub_list(i)%ub_kind%kind_number==iub) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(iub) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(ub_kind_set)-SUM(bad1)
             CALL allocate_ub_kind_set(new_ub_kind_set,counter,error)
             counter=0
             DO iub=1,SIZE(ub_kind_set)
                IF(bad1(iub)==0) THEN
                   counter=counter+1
                   new_ub_kind_set(counter)= ub_kind_set(iub)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(ub_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO iub=1,SIZE(ub_list)
                unsetme = .FALSE.
                IF(ub_list(iub)%ub_kind%id_type==do_ff_undef) unsetme = .TRUE.
                IF(ub_list(iub)%id_type==do_ff_undef)         unsetme = .TRUE.
                IF(unsetme) bad2(iub) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(ub_list)-SUM(bad2)
                ALLOCATE(new_ub_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO iub=1,SIZE(ub_list)
                   IF(bad2(iub)==0) THEN
                      counter=counter+1
                      new_ub_list(counter) = ub_list(iub)
                      newkind = ub_list(iub)%ub_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_ub_list(counter)%ub_kind => new_ub_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nub=SIZE(new_ub_list),&
                     ub_kind_set=new_ub_kind_set,&
                     ub_list=new_ub_list)
                DO iub=1,SIZE(new_ub_kind_set)
                   new_ub_kind_set(iub)%kind_number=iub
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CALL ub_kind_dealloc_ref(ub_kind_set,error=error)
             DEALLOCATE(ub_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New UB Count: ",&
                  SIZE(new_ub_list),SIZE(new_ub_kind_set)
             IF (iw>0) WRITE(iw,'(3I6)')(new_ub_list(iub)%a,new_ub_list(iub)%b,&
                  new_ub_list(iub)%c,iub=1,SIZE(new_ub_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. Count the number of UNSET torsion kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            ntorsion=ntorsion,&
            torsion_kind_set=torsion_kind_set,&
            torsion_list=torsion_list)
       IF(ntorsion>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old TORSION Count: ",&
               SIZE(torsion_list),SIZE(torsion_kind_set)
          IF (iw>0) WRITE(iw,'(4I6)')(torsion_list(itorsion)%a,torsion_list(itorsion)%b,&
               torsion_list(itorsion)%c,torsion_list(itorsion)%d,itorsion=1,SIZE(torsion_list))
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(torsion_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO itorsion=1,SIZE(torsion_kind_set)
             unsetme=.FALSE.
             IF(torsion_kind_set(itorsion)%id_type==do_ff_undef) unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(torsion_list)
                IF(torsion_list(i)%id_type/=do_ff_undef.AND.&
                     torsion_list(i)%torsion_kind%kind_number==itorsion) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(itorsion) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(torsion_kind_set)-SUM(bad1)
             CALL allocate_torsion_kind_set(new_torsion_kind_set,counter,error)
             counter=0
             DO itorsion=1,SIZE(torsion_kind_set)
                IF(bad1(itorsion)==0) THEN
                   counter=counter+1
                   new_torsion_kind_set(counter)= torsion_kind_set(itorsion)
                   i = SIZE(torsion_kind_set(itorsion)%m)
                   j = SIZE(torsion_kind_set(itorsion)%k)
                   k = SIZE(torsion_kind_set(itorsion)%phi0)
                   ALLOCATE(new_torsion_kind_set(counter)%m(i),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE(new_torsion_kind_set(counter)%k(i),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   ALLOCATE(new_torsion_kind_set(counter)%phi0(i),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                   new_torsion_kind_set(counter)%m= torsion_kind_set(itorsion)%m
                   new_torsion_kind_set(counter)%k= torsion_kind_set(itorsion)%k
                   new_torsion_kind_set(counter)%phi0= torsion_kind_set(itorsion)%phi0
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(torsion_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO itorsion=1,SIZE(torsion_list)
                unsetme = .FALSE.
                IF(torsion_list(itorsion)%torsion_kind%id_type==do_ff_undef) unsetme = .TRUE.
                IF(torsion_list(itorsion)%id_type==do_ff_undef)              unsetme = .TRUE.
                IF(unsetme) bad2(itorsion) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(torsion_list)-SUM(bad2)
                ALLOCATE(new_torsion_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO itorsion=1,SIZE(torsion_list)
                   IF(bad2(itorsion)==0) THEN
                      counter=counter+1
                      new_torsion_list(counter) = torsion_list(itorsion)
                      newkind = torsion_list(itorsion)%torsion_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_torsion_list(counter)%torsion_kind => new_torsion_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     ntorsion=SIZE(new_torsion_list),&
                     torsion_kind_set=new_torsion_kind_set,&
                     torsion_list=new_torsion_list)
                DO itorsion=1,SIZE(new_torsion_kind_set)
                   new_torsion_kind_set(itorsion)%kind_number=itorsion
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO itorsion=1,SIZE(torsion_kind_set)
                CALL torsion_kind_dealloc_ref(torsion_kind_set(itorsion),error=error)
             END DO
             DEALLOCATE(torsion_kind_set,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(torsion_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New TORSION Count: ",&
                  SIZE(new_torsion_list),SIZE(new_torsion_kind_set)
             IF (iw>0) WRITE(iw,'(4I6)')(new_torsion_list(itorsion)%a,new_torsion_list(itorsion)%b,&
                  new_torsion_list(itorsion)%c,new_torsion_list(itorsion)%d,itorsion=1,&
                  SIZE(new_torsion_list))
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Count the number of UNSET improper kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nimpr=nimpr,&
            impr_kind_set=impr_kind_set,&
            impr_list=impr_list)
       IF(nimpr>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old IMPROPER Count: ",&
               SIZE(impr_list),SIZE(impr_kind_set)
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(impr_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO iimpr=1,SIZE(impr_kind_set)
             unsetme=.FALSE.
             IF(impr_kind_set(iimpr)%id_type==do_ff_undef) unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(impr_list)
                IF(impr_list(i)%id_type/=do_ff_undef.AND.&
                     impr_list(i)%impr_kind%kind_number==iimpr) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(iimpr) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(impr_kind_set)-SUM(bad1)
             CALL allocate_impr_kind_set(new_impr_kind_set,counter,error)
             counter=0
             DO iimpr=1,SIZE(impr_kind_set)
                IF(bad1(iimpr)==0) THEN
                   counter=counter+1
                   new_impr_kind_set(counter)= impr_kind_set(iimpr)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(impr_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO iimpr=1,SIZE(impr_list)
                unsetme = .FALSE.
                IF(impr_list(iimpr)%impr_kind%id_type==do_ff_undef) unsetme = .TRUE.
                IF(impr_list(iimpr)%id_type==do_ff_undef)           unsetme = .TRUE.
                IF(unsetme) bad2(iimpr) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(impr_list)-SUM(bad2)
                ALLOCATE(new_impr_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO iimpr=1,SIZE(impr_list)
                   IF(bad2(iimpr)==0) THEN
                      counter=counter+1
                      new_impr_list(counter) = impr_list(iimpr)
                      newkind = impr_list(iimpr)%impr_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_impr_list(counter)%impr_kind => new_impr_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nimpr=SIZE(new_impr_list),&
                     impr_kind_set=new_impr_kind_set,&
                     impr_list=new_impr_list)
                DO iimpr=1,SIZE(new_impr_kind_set)
                   new_impr_kind_set(iimpr)%kind_number=iimpr
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DO iimpr=1,SIZE(impr_kind_set)
                CALL impr_kind_dealloc_ref(impr_kind_set(iimpr),error=error)
             END DO
             DEALLOCATE(impr_kind_set,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(impr_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New IMPROPER Count: ",&
                  SIZE(new_impr_list),SIZE(new_impr_kind_set)
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Count the number of UNSET opbends kinds there are
    !-----------------------------------------------------------------------------
    DO ikind=1,SIZE(molecule_kind_set)
       molecule_kind => molecule_kind_set(ikind)
       CALL get_molecule_kind(molecule_kind=molecule_kind,&
            nopbend=nopbend,&
            opbend_kind_set=opbend_kind_set,&
            opbend_list=opbend_list)
       IF(nopbend>0) THEN
          IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") Old OPBEND Count: ",&
               SIZE(opbend_list),SIZE(opbend_kind_set)
          NULLIFY(bad1,bad2)
          ALLOCATE(bad1(SIZE(opbend_kind_set)),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          bad1(:)=0
          DO iopbend=1,SIZE(opbend_kind_set)
             unsetme=.FALSE.
             IF(opbend_kind_set(iopbend)%id_type==do_ff_undef) unsetme = .TRUE.
             valid_kind=.FALSE.
             DO i=1,SIZE(opbend_list)
                IF(opbend_list(i)%id_type/=do_ff_undef.AND.&
                     opbend_list(i)%opbend_kind%kind_number==iopbend) THEN
                   valid_kind=.TRUE.
                   EXIT
                END IF
             END DO
             IF(.NOT.valid_kind) unsetme = .TRUE.
             IF(unsetme) bad1(iopbend) = 1
          END DO
          IF(SUM(bad1)/=0) THEN
             counter = SIZE(opbend_kind_set)-SUM(bad1)
             CALL allocate_opbend_kind_set(new_opbend_kind_set,counter,error)
             counter=0
             DO iopbend=1,SIZE(opbend_kind_set)
                IF(bad1(iopbend)==0) THEN
                   counter=counter+1
                   new_opbend_kind_set(counter)= opbend_kind_set(iopbend)
                END IF
             END DO
             counter=0
             ALLOCATE(bad2(SIZE(opbend_list)),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             bad2(:)=0
             DO iopbend=1,SIZE(opbend_list)
                unsetme = .FALSE.
                IF(opbend_list(iopbend)%opbend_kind%id_type==do_ff_undef) unsetme = .TRUE.
                IF(opbend_list(iopbend)%id_type==do_ff_undef)           unsetme = .TRUE.
                IF(unsetme) bad2(iopbend) = 1
             END DO
             IF(SUM(bad2)/=0) THEN
                counter = SIZE(opbend_list)-SUM(bad2)
                ALLOCATE(new_opbend_list(counter),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                counter=0
                DO iopbend=1,SIZE(opbend_list)
                   IF(bad2(iopbend)==0) THEN
                      counter=counter+1
                      new_opbend_list(counter) = opbend_list(iopbend)
                      newkind = opbend_list(iopbend)%opbend_kind%kind_number
                      newkind = newkind - SUM(bad1(1:newkind))
                      new_opbend_list(counter)%opbend_kind => new_opbend_kind_set(newkind)
                   END IF
                END DO
                CALL set_molecule_kind(molecule_kind=molecule_kind,&
                     nopbend=SIZE(new_opbend_list),&
                     opbend_kind_set=new_opbend_kind_set,&
                     opbend_list=new_opbend_list)
                DO iopbend=1,SIZE(new_opbend_kind_set)
                   new_opbend_kind_set(iopbend)%kind_number=iopbend
                END DO
             END IF
             DEALLOCATE(bad2,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(opbend_kind_set,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(opbend_list,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF(iw>0) WRITE(iw,*) "    Mol(",ikind,") New OPBEND Count: ",&
                  SIZE(new_opbend_list),SIZE(new_opbend_kind_set)
          END IF
          DEALLOCATE(bad1,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! 12. Count the number of UNSET NONBOND14 kinds there are
    !-                NEED TO REMOVE EXTRAS HERE   - IKUO
    !---------------------------------------------------------------------------
    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(handle)
  END SUBROUTINE clean_intra_force_kind

! *****************************************************************************
!> \brief Reads from the input structure all information for generic functions
! *****************************************************************************
  SUBROUTINE get_generic_info(gen_section, func_name, xfunction, parameters, values,&
       var_values, size_variables, i_rep_sec, input_variables, error)
    TYPE(section_vals_type), POINTER         :: gen_section
    CHARACTER(LEN=*), INTENT(IN)             :: func_name
    CHARACTER(LEN=default_path_length), &
      INTENT(OUT)                            :: xfunction
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: parameters
    REAL(KIND=dp), DIMENSION(:), POINTER     :: values
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: var_values
    INTEGER, INTENT(IN), OPTIONAL            :: size_variables, i_rep_sec
    CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: input_variables
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_generic_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par, my_par_tmp, my_units, &
                                                my_units_tmp, my_var
    INTEGER                                  :: i, ind, irep, isize, j, &
                                                mydim, n_par, n_units, n_val, &
                                                nblank, stat
    LOGICAL                                  :: check, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val, my_val_tmp

    failure = .FALSE.
    NULLIFY(my_var, my_par, my_val, my_par_tmp, my_val_tmp)
    NULLIFY (my_units)
    NULLIFY (my_units_tmp)
    IF (ASSOCIATED(parameters)) THEN
       DEALLOCATE(parameters,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(values)) THEN
       DEALLOCATE(values,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    irep = 1
    IF (PRESENT(i_rep_sec)) irep = i_rep_sec
    mydim = 0
    CALL section_vals_val_get(gen_section,TRIM(func_name),i_rep_section=irep,c_val=xfunction,error=error)
    CALL compress(xfunction, full=.TRUE.)
    IF (PRESENT(input_variables)) THEN
       ALLOCATE(my_var(SIZE(input_variables)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       my_var = input_variables
    ELSE
       CALL section_vals_val_get(gen_section,"VARIABLES",i_rep_section=irep,c_vals=my_var,error=error)
    END IF
    IF (ASSOCIATED(my_var)) THEN
       mydim = SIZE(my_var)
    END IF
    IF (PRESENT(size_variables)) THEN
       CPPrecondition(mydim==size_variables,cp_failure_level,routineP,error,failure)
    END IF
    ! Check for presence of Parameters
    CALL section_vals_val_get(gen_section,"PARAMETERS",i_rep_section=irep,n_rep_val=n_par,error=error)
    CALL section_vals_val_get(gen_section,"VALUES",i_rep_section=irep,n_rep_val=n_val,error=error)
    check = (n_par>0).EQV.(n_val>0)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    CALL section_vals_val_get(gen_section,"UNITS",i_rep_section=irep,n_rep_val=n_units,error=error)
    IF (n_par>0) THEN
       ! Parameters
       ALLOCATE(my_par(0),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(my_val(0),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, n_par
          isize  = SIZE(my_par)
          CALL section_vals_val_get(gen_section,"PARAMETERS",i_rep_section=irep,i_rep_val=i,c_vals=my_par_tmp,error=error)
          nblank = COUNT(my_par_tmp=="")
          CALL reallocate(my_par, 1, isize+SIZE(my_par_tmp)-nblank)
          ind = 0
          DO j = 1, SIZE(my_par_tmp)
             IF (my_par_tmp(j)=="") CYCLE
             ind = ind + 1
             my_par(isize+ind)=my_par_tmp(j)
          END DO
       END DO
       DO i = 1, n_val
          isize = SIZE(my_val)
          CALL section_vals_val_get(gen_section,"VALUES",i_rep_section=irep,i_rep_val=i,r_vals=my_val_tmp,error=error)
          CALL reallocate(my_val,1, isize+SIZE(my_val_tmp))
          my_val(isize+1:isize+SIZE(my_val_tmp))=my_val_tmp
       END DO
       CPPrecondition(SIZE(my_par)==SIZE(my_val),cp_failure_level,routineP,error,failure)
       ! Optionally read the units for each parameter value
       ALLOCATE (my_units(0),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (n_units > 0) THEN
          DO i=1,n_units
             isize = SIZE(my_units)
             CALL section_vals_val_get(gen_section,"UNITS",i_rep_section=irep,i_rep_val=i,c_vals=my_units_tmp,error=error)
             nblank = COUNT(my_units_tmp == "")
             CALL reallocate(my_units,1,isize+SIZE(my_units_tmp)-nblank)
             ind = 0
             DO j=1,SIZE(my_units_tmp)
                IF (my_units_tmp(j) == "") CYCLE
                ind = ind + 1
                my_units(isize+ind) = my_units_tmp(j)
             END DO
          END DO
          CPPrecondition(SIZE(my_units)==SIZE(my_val),cp_failure_level,routineP,error,failure)
       END IF
       mydim=mydim+SIZE(my_val)
       IF (SIZE(my_val)==0) THEN
          DEALLOCATE(my_par,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(my_val,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(my_units,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    ! Handle trivial case of a null function defined
    ALLOCATE(parameters(mydim),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(values(mydim),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (mydim>0) THEN
       parameters(1:SIZE(my_var)) = my_var
       values(1:SIZE(my_var))     = 0.0_dp
       IF (PRESENT(var_values)) THEN
          CPPrecondition(SIZE(var_values)==SIZE(my_var),cp_failure_level,routineP,error,failure)
          values(1:SIZE(my_var)) = var_values
       END IF
       IF (ASSOCIATED(my_val)) THEN
          DO i=1,SIZE(my_val)
             parameters(SIZE(my_var)+i) = my_par(i)
             IF (n_units > 0) THEN
                values(SIZE(my_var)+i) = cp_unit_to_cp2k(my_val(i),TRIM(ADJUSTL(my_units(i))),error=error)
             ELSE
                values(SIZE(my_var)+i) = my_val(i)
             END IF
          END DO
       END IF
    END IF
    IF (ASSOCIATED(my_par)) THEN
       DEALLOCATE(my_par,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(my_val)) THEN
       DEALLOCATE(my_val,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(my_units)) THEN
       DEALLOCATE(my_units,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (PRESENT(input_variables)) THEN
       DEALLOCATE(my_var,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE get_generic_info

END MODULE force_fields_util
