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
MODULE force_fields
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE f77_blas
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
  USE force_field_kind_types,          ONLY: do_ff_amber,&
                                             do_ff_charmm,&
                                             do_ff_g87,&
                                             do_ff_g96,&
                                             do_ff_undef
  USE force_field_types,               ONLY: deallocate_ff_type,&
                                             force_field_type,&
                                             init_ff_type
  USE force_fields_ext,                ONLY: read_force_field_amber,&
                                             read_force_field_charmm,&
                                             read_force_field_gromos
  USE force_fields_input,              ONLY: read_force_field_section
  USE force_fields_util,               ONLY: clean_intra_force_kind,&
                                             force_field_pack,&
                                             force_field_qeff_output
  USE input_constants,                 ONLY: do_skip_14
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields'

  PRIVATE
  PUBLIC :: force_field_control

CONTAINS

! *****************************************************************************
!> \brief 1. If reading in from external file, make sure its there first
!>      2. Read in the force_field from the corresponding locations
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param ewald_env ...
!> \param fist_nonbond_env ...
!> \param root_section ...
!> \param para_env ...
!> \param qmmm ...
!> \param qmmm_env ...
!> \param subsys_section ...
!> \param mm_section ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param cell ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE force_field_control(atomic_kind_set, particle_set, &
       molecule_kind_set, molecule_set, ewald_env, fist_nonbond_env, &
       root_section, para_env, qmmm, qmmm_env, subsys_section, mm_section, &
       shell_particle_set, core_particle_set, cell, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section, mm_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_particle_set, &
                                                core_particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: exclude_ei, exclude_vdw, &
                                                handle, iw
    LOGICAL                                  :: failure, found
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(force_field_type)                   :: ff_type
    TYPE(section_vals_type), POINTER         :: topology_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    !-----------------------------------------------------------------------------
    ! 1. Initialize the ff_type structure type
    !-----------------------------------------------------------------------------
    CALL init_ff_type(ff_type,error)

    !-----------------------------------------------------------------------------
    ! 2. Read in the force field section in the input file if any
    !-----------------------------------------------------------------------------
    CALL read_force_field_section(ff_type,para_env,mm_section,error)

    !-----------------------------------------------------------------------------
    ! 2.1 In case exclusion 1-4 was requested, we need to modify the values of
    !     the scale factors setting them to zero..
    !-----------------------------------------------------------------------------
    topology_section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY",error=error)
    CALL section_vals_val_get(topology_section,"EXCLUDE_VDW",i_val=exclude_vdw,error=error)
    CALL section_vals_val_get(topology_section,"EXCLUDE_EI",i_val=exclude_ei,error=error)
    IF (exclude_vdw==do_skip_14) ff_type%vdw_scale14 = 0.0_dp
    IF (exclude_ei ==do_skip_14) ff_type%ei_scale14  = 0.0_dp

    !-----------------------------------------------------------------------------
    ! 3. If reading in from external file, make sure its there first
    !-----------------------------------------------------------------------------
    SELECT CASE (ff_type%ff_type)
    CASE (do_ff_charmm,do_ff_amber,do_ff_g96,do_ff_g87)
       INQUIRE(FILE=ff_type%ff_file_name,EXIST=found)
       IF (.NOT.found) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Force field file missing",para_env)
       END IF
    CASE (do_ff_undef)
       ! Do Nothing
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Force field type not implemented",para_env)
    END SELECT

    !-----------------------------------------------------------------------------
    ! 4. Read in the force field from the corresponding locations
    !-----------------------------------------------------------------------------
    SELECT CASE (ff_type%ff_type)
    CASE (do_ff_charmm)
       CALL read_force_field_charmm(ff_type,para_env,mm_section,error)
    CASE (do_ff_amber)
       CALL read_force_field_amber(ff_type,para_env,mm_section,particle_set,error)
    CASE (do_ff_g87,do_ff_g96)
       CALL read_force_field_gromos(ff_type,para_env,mm_section,error)
    CASE (do_ff_undef)
       ! Do Nothing
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Force field type not implemented",para_env)
    END SELECT

    !-----------------------------------------------------------------------------
    ! 5. Possibly print the top file
    !-----------------------------------------------------------------------------
    CALL print_pot_parameter_file(ff_type, mm_section, error)

    !-----------------------------------------------------------------------------
    ! 6. Pack all force field info into different structures
    !-----------------------------------------------------------------------------
    CALL force_field_pack (particle_set,atomic_kind_set, molecule_kind_set,molecule_set,&
         ewald_env,fist_nonbond_env,ff_type,root_section, qmmm, qmmm_env, mm_section,&
         subsys_section, shell_particle_set=shell_particle_set, &
         core_particle_set=core_particle_set, cell=cell, error=error)

    !-----------------------------------------------------------------------------
    ! 7. Output total system charge assigned to qeff
    !-----------------------------------------------------------------------------
    CALL force_field_qeff_output (particle_set, molecule_kind_set,&
         molecule_set, mm_section, fist_nonbond_env%charges, error)

    !-----------------------------------------------------------------------------
    ! 8. Clean up "UNSET" bond,bend,UB,TORSION,IMPR,ONFO kinds
    !-----------------------------------------------------------------------------
    CALL clean_intra_force_kind (molecule_kind_set,mm_section,error)

    !-----------------------------------------------------------------------------
    ! 9. Cleanup the ff_type structure type
    !-----------------------------------------------------------------------------
    CALL deallocate_ff_type(ff_type, error=error)

    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE force_field_control

! *****************************************************************************
!> \brief Prints force field information in a pot file
!> \param ff_type ...
!> \param mm_section ...
!> \param error ...
!> \author Teodoro Laino [tlaino, teodoro.laino-AT-gmail.com] - 11.2008
! *****************************************************************************
  SUBROUTINE print_pot_parameter_file(ff_type, mm_section, error)

    TYPE(force_field_type)                   :: ff_type
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_pot_parameter_file', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw, m
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eps, k, phi0, r0, sigma, &
                                                theta0
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,mm_section,"PRINT%FF_PARAMETER_FILE",&
         error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_PARAMETER_FILE",&
            middle_name="force_field",extension=".pot",error=error)
       IF (iw>0) THEN
          ! Header
          WRITE(iw,1000)"Force Field Parameter File dumped into CHARMM FF style"
       END IF
       SELECT CASE(ff_type%ff_type)
       CASE (do_ff_charmm)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Dumping FF parameter file for CHARMM FF  not implemented!", &
               error=error, error_level=cp_warning_level)
       CASE (do_ff_amber)
          IF (iw>0) THEN
             ! Bonds
             WRITE(iw,1001)
             DO i = 1, SIZE(ff_type%amb_info%bond_a)
                k = cp_unit_from_cp2k(ff_type%amb_info%bond_k(i),"kcalmol*angstrom^-2",error=error)
                r0= cp_unit_from_cp2k(ff_type%amb_info%bond_r0(i),"angstrom",error=error)
                WRITE(iw, 2001) ff_type%amb_info%bond_a(i),&
                                ff_type%amb_info%bond_b(i),&
                                k, r0
             END DO
             ! Angles
             WRITE(iw,1002)
             DO i = 1, SIZE(ff_type%amb_info%bend_a)
                k      = cp_unit_from_cp2k(ff_type%amb_info%bend_k(i),"kcalmol*rad^-2",error=error)
                theta0 = cp_unit_from_cp2k(ff_type%amb_info%bend_theta0(i),"deg",error=error)
                WRITE(iw, 2002) ff_type%amb_info%bend_a(i),&
                                ff_type%amb_info%bend_b(i),&
                                ff_type%amb_info%bend_c(i),&
                                k, theta0
             END DO
             ! Torsions
             WRITE(iw,1003)
             DO i = 1, SIZE(ff_type%amb_info%torsion_a)
                k    = cp_unit_from_cp2k(ff_type%amb_info%torsion_k(i),"kcalmol",error=error)
                m    = ff_type%amb_info%torsion_m(i)
                phi0 = cp_unit_from_cp2k(ff_type%amb_info%torsion_phi0(i),"deg",error=error)
                WRITE(iw, 2003) ff_type%amb_info%torsion_a(i),&
                                ff_type%amb_info%torsion_b(i),&
                                ff_type%amb_info%torsion_c(i),&
                                ff_type%amb_info%torsion_d(i),&
                                k, m, phi0
             END DO
             ! Lennard-Jones
             WRITE(iw,1005)
             DO i = 1, SIZE(ff_type%amb_info%nonbond_a)
                eps   = cp_unit_from_cp2k(ff_type%amb_info%nonbond_eps(i),"kcalmol",error=error)
                sigma = cp_unit_from_cp2k(ff_type%amb_info%nonbond_rmin2(i),"angstrom", error=error)
                WRITE(iw, 2005) ff_type%amb_info%nonbond_a(i),&
                                eps, sigma
             END DO
          END IF
       CASE (do_ff_g87,do_ff_g96)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Dumping FF parameter file for GROMOS FF not implemented!", &
               error=error, error_level=cp_warning_level)
       CASE (do_ff_undef)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Dumping FF parameter file for INPUT FF  not implemented!", &
               error=error, error_level=cp_warning_level)
       END SELECT
       IF (iw>0) THEN
          WRITE(iw,'(/,A)')"END"
       END IF
       CALL cp_print_key_finished_output(iw,logger,mm_section,&
            "PRINT%FF_PARAMETER_FILE",error=error)
    END IF
    CALL timestop(handle)
    RETURN
1000 FORMAT("*>>>>>>>",T12,A,T73,"<<<<<<<")
1001 FORMAT(/,"BONDS",/,"!",/,"!V(bond) = Kb(b - b0)**2",/,"!",/,"!Kb: kcal/mole/A**2",/, &
         "!b0: A",/,"!",/,"! atom type           Kb              b0",/,"!")
1002 FORMAT(/,"ANGLES",/,"!",/,"!V(angle) = Ktheta(Theta - Theta0)**2",/,"!",/,           &
         "!V(Urey-Bradley) = Kub(S - S0)**2",/,"!",/,"!Ktheta: kcal/mole/rad**2",/,  &
         "!Theta0: degrees",/,"!Kub: kcal/mole/A**2 (Urey-Bradley)",/,"!S0: A",/,    &
         "!",/,"!   atom types              Ktheta          Theta0       Kub        S0",/,"!")
1003 FORMAT(/,"DIHEDRALS",/,"!",/,"!V(dihedral) = Kchi(1 + cos(n(chi) - delta))",/,       &
         "!",/,"!Kchi: kcal/mole",/,"!n: multiplicity",/,"!delta: degrees",/,        &
         "!",/,"!     atom types                    Kchi       n       delta",/,"!")
1005 FORMAT(/,"NONBONDED",/,"!",/,                                                        &
         "!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]",/, &
         "!",/,"!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)",/,               &
         "!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j",/,"!",/,                       &
         "!atom         ignored        epsilon       Rmin/2      ignored   eps,1-4       "//&
         "Rmin/2,1-4",/,"!")

2001 FORMAT(A6,1X,A6,1X,2F15.9)                     ! bond
2002 FORMAT(A6,1X,A6,1X,A6,1X,2F15.9)               ! angle
2003 FORMAT(A6,1X,A6,1X,A6,1X,A6,1X,F15.9,I5,F15.9) ! torsion
2005 FORMAT(A6,1X,"    0.000000000",2F15.9)         ! nonbond
  END SUBROUTINE print_pot_parameter_file

END MODULE force_fields
