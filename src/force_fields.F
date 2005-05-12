!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/force_fields [1.0] *
!!
!!   NAME
!!     force_fields
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Subroutine input_torsions changed (DG) 05-Dec-2000
!!     Output formats changed (DG) 05-Dec-2000
!!     JGH (26-01-2002) : force field parameters stored in tables, not in
!!       matrices. Input changed to have parameters labeled by the position
!!       and not atom pairs (triples etc)
!!
!!   SOURCE
!******************************************************************************

MODULE force_fields

  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             get_atomic_kind,&
                                             init_atomic_kind_set,&
                                             set_atomic_kind
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_files,                        ONLY: open_file
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             allocate_potential,&
                                             get_potential,&
                                             set_potential
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type, &
                                             fist_nonbond_env_create, &
                                             fist_nonbond_env_release
  USE force_field_types,               ONLY: allocate_bond_kind_set,&
                                             allocate_bend_kind_set,&
                                             allocate_ub_kind_set,&
                                             ub_kind_dealloc_ref,&
                                             torsion_kind_dealloc_ref,&
                                             impr_kind_dealloc_ref,&
                                             allocate_onfo_kind_set,&
                                             deallocate_bond_kind_set,&
                                             deallocate_bend_kind_set,&
                                             allocate_impr_kind_set,&
                                             allocate_torsion_kind_set,&
                                             bend_kind_type,&
                                             bond_kind_type,&
                                             impr_kind_type,&
                                             onfo_kind_type,&
                                             torsion_kind_type,&
                                             ub_kind_type,&
                                             amber_info_type,&
                                             charmm_info_type,&
                                             force_field_type,&
                                             input_info_type,&
                                             init_ff_type,&
                                             deallocate_ff_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE manybody_potential,              ONLY: spline_manybody_control
  USE mathconstants,                   ONLY: pi, twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set,&
                                             atom_type,&
                                             bend_type,&
                                             bond_type,&
                                             dist_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             impr_type,&
                                             molecule_kind_type,&
                                             onfo_type,&
                                             set_molecule_kind,&
                                             torsion_type,&
                                             ub_type
  USE molecule_types_new,              ONLY: molecule_type,&
                                             get_molecule
  USE pair_potential,                  ONLY: spline_nonbond_control
  USE pair_potential_types,            ONLY: pair_potential_type,&
                                             eam_pot_type,       &
                                             reallocate_potparm, &
                                             pot_parameters_create, &
                                             pot_parameters_release,&
                                             clean_pot_parameter
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE particle_types,                  ONLY: particle_type
  USE splines,                         ONLY: spline_environment_type, &
                                             spline_env_create, &
                                             spline_env_release, &
                                             spline_env_retain,&
                                             spline_data_release, &
                                             spline_data_retain, &
                                             spline_data_type
  USE string_utilities,                ONLY: str_search,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             constraint_info_type
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE input_constants,                 ONLY: do_qmmm_none
  USE physcon,                         ONLY: bohr

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "force_fields"

  PRIVATE
  PUBLIC :: force_field_control,&
            read_force_field_section,&
            read_force_field_charmm,&
            read_force_field_amber

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** force_fields/force_field_control [1.0] *
!!
!!   NAME
!!     force_field_control
!!
!!   FUNCTION
!!     1. If reading in from external file, make sure its there first
!!     2. Read in the force_field from the corresponding locations
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE force_field_control (atomic_kind_set,particle_set,&
                                molecule_kind_set,molecule_set,&
                                ewald_env,fist_nonbond_env,&
                                globenv, qmmm, qmmm_env, error)

  TYPE(atomic_kind_type), DIMENSION(:), POINTER       :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER          :: particle_set
  TYPE(molecule_kind_type), DIMENSION(:), POINTER     :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), POINTER          :: molecule_set
  TYPE ( ewald_environment_type ), POINTER :: ewald_env
  TYPE ( fist_nonbond_env_type ), POINTER ::  fist_nonbond_env
  TYPE(global_environment_type), INTENT(IN),OPTIONAL  :: globenv
  LOGICAL, INTENT(IN), OPTIONAL :: qmmm
  TYPE (qmmm_env_mm_type), POINTER, OPTIONAL :: qmmm_env
  TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error

  CHARACTER(len=*), PARAMETER :: routineN = 'force_field_control', &
      routineP = moduleN//':'//routineN

  INTEGER                                  :: handle, iw
  LOGICAL                                  :: found
  TYPE(force_field_type)                   :: ff_type

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "Entering force_field_control"

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Initialize the ff_type structure type
  !-----------------------------------------------------------------------------
  CALL init_ff_type(ff_type)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in the force field section in the input file if any
  !-----------------------------------------------------------------------------
  CALL read_force_field_section(ff_type,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. If reading in from external file, make sure its there first
  !-----------------------------------------------------------------------------
  IF(ff_type%ff_type/="OFF") THEN
    !------------------------------------------------------
    !Make sure the file really do exist before proceding
    IF(globenv%ionode) THEN
      INQUIRE(FILE=ff_type%ff_file_name,EXIST=found)
      IF(found) THEN
      ELSE
        CALL stop_program ("force_field_control","force field file missing")
      END IF
    END IF
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Read in the force field from the corresponding locations
  !-----------------------------------------------------------------------------
  IF(ff_type%ff_type/="OFF") THEN
    SELECT CASE (ff_type%ff_type)
    CASE ("CHM")
      CALL read_force_field_charmm(ff_type,globenv)
    CASE ("AMB")
      CALL stop_program ("force_field_control",&
                         "force field type not implemented")
      CALL read_force_field_amber(ff_type,globenv)
    CASE DEFAULT
      CALL stop_program ("force_field_control",&
                         "force field type not implemented")
    END SELECT
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Pack all force field info into different structures
  !-----------------------------------------------------------------------------
  CALL force_field_pack (particle_set,atomic_kind_set,&
                         molecule_kind_set,molecule_set,&
                         ewald_env,fist_nonbond_env,&
                         ff_type,globenv, qmmm, qmmm_env, error=error)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Output total system charge assigned to qeff
  !-----------------------------------------------------------------------------
  CALL force_field_qeff_output (particle_set,atomic_kind_set,&
                                molecule_kind_set,molecule_set,&
                                globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Clean up "UNSET" bond,bend,UB,TORSION,IMPR,ONFO kinds
  !-----------------------------------------------------------------------------
  CALL clean_intra_force_kind (molecule_kind_set,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Cleanup the ff_type structure type
  !-----------------------------------------------------------------------------
  CALL deallocate_ff_type(ff_type)


  IF(PM1) WRITE(iw,*) "Exiting  force_field_control"

  CALL timestop(0.0_dp,handle)
  CALL write_checkpoint_information("leaving "//routineN,globenv)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE force_field_control

!******************************************************************************
!!****** force_fields/read_force_field_section [1.0] *
!!
!!   NAME
!!     read_force_field_section
!!
!!   FUNCTION
!!     Reads the force_field input section
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (30.11.2001) : moved determination of setup variables to 
!!                        molecule_input 
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_force_field_section ( ff_type , globenv )


    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv


    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length)     :: string

    INTEGER                                  :: iw,ilen,ierror

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering read_force_field_section"


  CALL parser_init(globenv%input_file_name,globenv)


!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror /= 0 ) THEN
     IF (globenv%ionode .AND. globenv%print_level >= 2) THEN
        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 9
        CALL get_next(string,ilen)
        CALL uppercase(string)
        IF (string == 'PARMFILE') THEN
          ilen = 40
          CALL get_next(ff_type%ff_type,ilen)
          CALL uppercase(ff_type%ff_type)
          CALL get_next(ff_type%ff_file_name,ilen)
        ELSE IF(string == 'RCUT_NB') THEN
          CALL get_next(ff_type%rcut_nb)
          CALL convert_to_cp2k_units("ANGSTROM",length=ff_type%rcut_nb)
        END IF
        CALL read_line
     END DO
   END IF
    CALL parser_end()

  CALL parser_init(globenv%input_file_name,globenv)
!..read in non-bonded interaction information
    CALL input_nonbonded(ff_type,globenv)
    CALL input_nonbonded14(ff_type,globenv)

!..read the charges parameters
    CALL input_charges(ff_type,globenv)

!..read the bond parameters
    CALL input_bonds(ff_type,globenv)

!..read the bend parameters
    CALL input_bends(ff_type,globenv)

!..read the torsion parameters
    CALL input_torsions(ff_type,globenv)

!..read the distance constraint parameters
    CALL parser_end()
 
  IF(PM1) WRITE(iw,*) "  Exiting  read_force_field_section"

END SUBROUTINE read_force_field_section


!!*****
!******************************************************************************
!!****** force_fields/input_nonbonded14 [1.0] *
!!
!!   NAME
!!     input_nonbonded14
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I nonbonded14                                                               I
!!I    lennard-jones atom1 atom2 epsilon sigma rcut                           I
!!I    williams      atom1 atom2 a b c rcut                                   I
!!I    goodwin       atom1 atom2 vr0 d dc m mc rcut                           I
!!I end nonbonded14                                                           I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_nonbonded14(ff_type,globenv)

    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(input_info_type),POINTER            :: inp_info
    CHARACTER(LEN=default_string_length)     :: at1, at2
    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length)     :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw, &
                                                np
    REAL(KIND=dp)                                :: a, b, c, d, dc, epsilon, m, &
                                                mc, rcut, sigma, vr0

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

  inp_info => ff_type%inp_info

  np = 0

!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror /= 0 ) THEN
     IF (globenv%ionode .AND. globenv%print_level >= 2) THEN
        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 9
        CALL get_next(string,ilen)
        CALL uppercase(string)
        IF (index(string,'NONBONDED14') /= 0 ) THEN
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_nonbonded14','unknown option')
              CASE ( 'LENNARD-JONES')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded14,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( epsilon )
                 CALL get_next ( sigma )
                 CALL get_next ( rcut )
                 CALL convert_to_cp2k_units ( "KELVIN", energy = epsilon )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = sigma )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = rcut )
                 inp_info%nonbonded14(np) %type = string2
                 inp_info%nonbonded14(np) %at1 = at1
                 inp_info%nonbonded14(np) %at2 = at2
                 inp_info%nonbonded14(np) %lj%epsilon = epsilon
                 inp_info%nonbonded14(np) %lj%sigma6 = sigma**6
                 inp_info%nonbonded14(np) %lj%sigma12 = sigma**12
                 inp_info%nonbonded14(np) %rcutsq = rcut*rcut
              CASE ( 'WILLIAMS')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded14,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( a )
                 CALL get_next ( b )
                 CALL get_next ( c )
                 CALL get_next ( rcut )
                 CALL convert_to_cp2k_units ( "KELVIN", energy = a )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = b, l_power = -1  )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = c, l_power = 6  )
                 CALL convert_to_cp2k_units ( "KELVIN", energy = c )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = rcut )
                 inp_info%nonbonded14(np) %type = string2
                 inp_info%nonbonded14(np) %at1 = at1
                 inp_info%nonbonded14(np) %at2 = at2
                 inp_info%nonbonded14(np) %willis%a = a
                 inp_info%nonbonded14(np) %willis%b = b
                 inp_info%nonbonded14(np) %willis%c = c
                 inp_info%nonbonded14(np) %rcutsq = rcut*rcut
              CASE ( 'GOODWIN')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded14,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( vr0 )
                 CALL get_next ( d )
                 CALL get_next ( dc )
                 CALL get_next ( m )
                 CALL get_next ( mc )
                 CALL get_next ( rcut )
                 inp_info%nonbonded14(np) %type = string2
                 inp_info%nonbonded14(np) %at1 = at1
                 inp_info%nonbonded14(np) %at2 = at2
                 inp_info%nonbonded14(np) %goodwin%vr0 = vr0
                 inp_info%nonbonded14(np) %goodwin%d = d
                 inp_info%nonbonded14(np) %goodwin%dc = dc
                 inp_info%nonbonded14(np) %goodwin%m = m
                 inp_info%nonbonded14(np) %goodwin%mc = mc
                 inp_info%nonbonded14(np) %rcutsq = rcut*rcut
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
  END IF

END SUBROUTINE input_nonbonded14


!!*****
!******************************************************************************
!!****** force_fields/input_nonbonded [1.0] *
!!
!!   NAME
!!     input_nonbonded
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I nonbonded                                                                 I
!!I    lennard-jones atom1 atom2 epsilon sigma rcut                           I
!!I    williams      atom1 atom2 a b c rcut                                   I
!!I    goodwin       atom1 atom2 vr0 d dc m mc rcut                           I
!!I    ipbv          atom1 atom2  rcut                                        I
!!I end nonbonded                                                             I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_nonbonded(ff_type,globenv)

    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(input_info_type),POINTER            :: inp_info
    CHARACTER(LEN=default_string_length)     :: at1, at2
    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length)     :: string, string2, file_name
    INTEGER                                  :: i, ierror, ilen, iw, &
                                                np
    REAL(KIND=dp)                                :: a, b, c, d, dc, epsilon, m, &
                                                mc, rcut, sigma, vr0

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

  inp_info => ff_type%inp_info

  np = 0

!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror /= 0 ) THEN
     IF (globenv%ionode .AND. globenv%print_level >= 2) THEN
        WRITE ( iw, '( A )' ) ' No input section &FORCE_FIELD found on file '
     END IF
  ELSE
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 9
        CALL get_next(string,ilen)
        CALL uppercase(string)
        IF (index(string,'NONBONDED') /= 0 ) THEN
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_nonbonded','unknown option')
              CASE ( 'LENNARD-JONES')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( epsilon )
                 CALL get_next ( sigma )
                 CALL get_next ( rcut )
                 CALL convert_to_cp2k_units ( "KELVIN", energy = epsilon )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = sigma )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = rcut )
                 inp_info%nonbonded(np) %type = string2
                 inp_info%nonbonded(np) %at1 = at1
                 inp_info%nonbonded(np) %at2 = at2
                 inp_info%nonbonded(np) %lj%epsilon = epsilon
                 inp_info%nonbonded(np) %lj%sigma6 = sigma**6
                 inp_info%nonbonded(np) %lj%sigma12 = sigma**12
                 inp_info%nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'WILLIAMS')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( a )
                 CALL get_next ( b )
                 CALL get_next ( c )
                 CALL get_next ( rcut )
                 CALL convert_to_cp2k_units ( "KELVIN", energy = a )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = b, l_power = -1  )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = c, l_power = 6  )
                 CALL convert_to_cp2k_units ( "KELVIN", energy = c )
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = rcut )
                 inp_info%nonbonded(np) %type = string2
                 inp_info%nonbonded(np) %at1 = at1
                 inp_info%nonbonded(np) %at2 = at2
                 inp_info%nonbonded(np) %willis%a = a
                 inp_info%nonbonded(np) %willis%b = b
                 inp_info%nonbonded(np) %willis%c = c
                 inp_info%nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'EAM')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 inp_info%nonbonded(np) %type = string2
                 inp_info%nonbonded(np) %at1 = at1
                 inp_info%nonbonded(np) %at2 = at2
                 ilen = 80
                 CALL get_next(file_name,ilen)
                 inp_info%nonbonded(np) %eam%eam_file_name = TRIM ( file_name )
                 CALL read_eam_data ( inp_info%nonbonded(np)%eam, globenv )
                 inp_info%nonbonded(np) %rcutsq = inp_info%nonbonded(np)%eam%acutal* & 
                                                  inp_info%nonbonded(np)%eam%acutal
              CASE ( 'IPBV')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( rcut )
                 IF ( rcut /= 7.8_dp ) &
                 WRITE ( *, * ) '****IPBV cut-off changed to 7.8 Angstroms****'
                 rcut = 7.8_dp
                 CALL convert_to_cp2k_units ( "ANGSTROM", length = rcut )
                 inp_info%nonbonded(np) %type = string2
                 inp_info%nonbonded(np) %at1 = at1
                 inp_info%nonbonded(np) %at2 = at2
                 IF  ( ( at1(1:1) == 'O' ) .AND. ( at2 (1:1) == 'O' ) ) THEN
                   inp_info%nonbonded(np) %ipbv%rcore = 0.9_dp       ! a.u.
                   inp_info%nonbonded(np) %ipbv%m = 0.0_dp           ! Kelvin/a.u. 
                   inp_info%nonbonded(np) %ipbv%b = 1.0e23_dp        ! Hartree
                   ! Hartree*a.u.^2
                   inp_info%nonbonded(np) %ipbv% a ( 2 ) = 4.786380682394_dp 
                   inp_info%nonbonded(np) %ipbv% a ( 3 ) = -1543.407053545_dp
                   inp_info%nonbonded(np) %ipbv% a ( 4 ) = 88783.31188529_dp
                   inp_info%nonbonded(np) %ipbv% a ( 5 ) = -2361200.155376_dp
                   inp_info%nonbonded(np) %ipbv% a ( 6 ) = 35940504.84679_dp
                   inp_info%nonbonded(np) %ipbv% a ( 7 ) = -339762743.6358_dp
                   inp_info%nonbonded(np) %ipbv% a ( 8 ) = 2043874926.466_dp
                   inp_info%nonbonded(np) %ipbv% a ( 9 ) = -7654856796.383_dp
                   inp_info%nonbonded(np) %ipbv% a ( 10 ) = 16195251405.65_dp
                   inp_info%nonbonded(np) %ipbv% a ( 11 ) = -13140392992.18_dp
                   inp_info%nonbonded(np) %ipbv% a ( 12 ) = -9285572894.245_dp 
                   inp_info%nonbonded(np) %ipbv% a ( 13 ) = 8756947519.029_dp
                   inp_info%nonbonded(np) %ipbv% a ( 14 ) = 15793297761.67_dp
                   inp_info%nonbonded(np) %ipbv% a ( 15 ) = 12917180227.21_dp
                 ELSEIF ( ( at1(1:1) == 'O' ) .AND. ( at2 (1:1) == 'H' ) ) THEN
                   inp_info%nonbonded(np) %ipbv%rcore = 2.95_dp       ! a.u.
                   inp_info%nonbonded(np) %ipbv%m =  -0.00437054E0_dp ! Hartree/a.u. 
                   inp_info%nonbonded(np) %ipbv%b =  0.0313652_dp     ! Hartree
                   ! Hartree*a.u.^2
                   inp_info%nonbonded(np) %ipbv% a ( 2 ) = -195.7716013277_dp
                   inp_info%nonbonded(np) %ipbv% a ( 3 ) = 15343.78613395_dp
                   inp_info%nonbonded(np) %ipbv% a ( 4 ) = -530864.4586516_dp
                   inp_info%nonbonded(np) %ipbv% a ( 5 ) = 10707934.39058_dp
                   inp_info%nonbonded(np) %ipbv% a ( 6 ) = -140099704.7890_dp
                   inp_info%nonbonded(np) %ipbv% a ( 7 ) = 1250943273.785_dp
                   inp_info%nonbonded(np) %ipbv% a ( 8 ) = -7795458330.676_dp
                   inp_info%nonbonded(np) %ipbv% a ( 9 ) = 33955897217.31_dp
                   inp_info%nonbonded(np) %ipbv% a ( 10 ) = -101135640744.0_dp
                   inp_info%nonbonded(np) %ipbv% a ( 11 ) = 193107995718.7_dp
                   inp_info%nonbonded(np) %ipbv% a ( 12 ) = -193440560940.0_dp
                   inp_info%nonbonded(np) %ipbv% a ( 13 ) =-4224406093.918E0_dp
                   inp_info%nonbonded(np) %ipbv% a ( 14 ) = 217192386506.5E0_dp 
                   inp_info%nonbonded(np) %ipbv% a ( 15 ) = -157581228915.5_dp
                 ELSEIF ( ( at1(1:1) == 'H' ) .AND. ( at2 (1:1) == 'O' ) ) THEN
                   inp_info%nonbonded(np) %ipbv%rcore = 2.95_dp       ! a.u.
                   inp_info%nonbonded(np) %ipbv%m =  -0.00437054E0_dp ! Hartree/a.u. 
                   inp_info%nonbonded(np) %ipbv%b =  0.0313652_dp     ! Hartree
                   ! Hartree*a.u.^2
                   inp_info%nonbonded(np) %ipbv% a ( 2 ) = -195.7716013277_dp
                   inp_info%nonbonded(np) %ipbv% a ( 3 ) = 15343.78613395_dp
                   inp_info%nonbonded(np) %ipbv% a ( 4 ) = -530864.4586516_dp
                   inp_info%nonbonded(np) %ipbv% a ( 5 ) = 10707934.39058_dp
                   inp_info%nonbonded(np) %ipbv% a ( 6 ) = -140099704.7890_dp
                   inp_info%nonbonded(np) %ipbv% a ( 7 ) = 1250943273.785_dp
                   inp_info%nonbonded(np) %ipbv% a ( 8 ) = -7795458330.676_dp
                   inp_info%nonbonded(np) %ipbv% a ( 9 ) = 33955897217.31_dp
                   inp_info%nonbonded(np) %ipbv% a ( 10 ) = -101135640744.0_dp
                   inp_info%nonbonded(np) %ipbv% a ( 11 ) = 193107995718.7_dp
                   inp_info%nonbonded(np) %ipbv% a ( 12 ) = -193440560940.0_dp
                   inp_info%nonbonded(np) %ipbv% a ( 13 ) =-4224406093.918E0_dp
                   inp_info%nonbonded(np) %ipbv% a ( 14 ) = 217192386506.5E0_dp 
                   inp_info%nonbonded(np) %ipbv% a ( 15 ) = -157581228915.5_dp
                 ELSEIF ( ( at1(1:1) == 'H' ) .AND. ( at2 (1:1) == 'H' ) ) THEN
                   inp_info%nonbonded(np) %ipbv%rcore = 3.165_dp     ! a.u.
                   inp_info%nonbonded(np) %ipbv%m = 0.00377391E0_dp  ! Hartree/a.u. 
                   inp_info%nonbonded(np) %ipbv%b = -0.02062040_dp   ! Hartree
                   ! Hartree*a.u.^2
                   inp_info%nonbonded(np) %ipbv% a ( 2 ) = -26.29456010782_dp
                   inp_info%nonbonded(np) %ipbv% a ( 3 ) =  2373.352548248_dp
                   inp_info%nonbonded(np) %ipbv% a ( 4 ) = -93880.43551360_dp
                   inp_info%nonbonded(np) %ipbv% a ( 5 ) =  2154624.884809_dp
                   inp_info%nonbonded(np) %ipbv% a ( 6 ) = -31965151.34955_dp
                   inp_info%nonbonded(np) %ipbv% a ( 7 ) = 322781785.3278_dp
                   inp_info%nonbonded(np) %ipbv% a ( 8 ) = -2271097368.668_dp
                   inp_info%nonbonded(np) %ipbv% a ( 9 ) = 11169163192.90_dp
                   inp_info%nonbonded(np) %ipbv% a ( 10 ) =-37684457778.47_dp
                   inp_info%nonbonded(np) %ipbv% a ( 11 ) = 82562104256.03_dp
                   inp_info%nonbonded(np) %ipbv% a ( 12 ) = -100510435213.4_dp
                   inp_info%nonbonded(np) %ipbv% a ( 13 ) = 24570342714.65E0_dp 
                   inp_info%nonbonded(np) %ipbv% a ( 14 ) = 88766181532.94E0_dp
                   inp_info%nonbonded(np) %ipbv% a ( 15 ) = -79705131323.98_dp
                 ELSE
                   CALL stop_program ("input_nonbonded","IPBV only for WATER")
                 ENDIF
                 inp_info%nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'GOODWIN')
                 np = np + 1
                 CALL reallocate_potparm(inp_info%nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( vr0 )
                 CALL get_next ( d )
                 CALL get_next ( dc )
                 CALL get_next ( m )
                 CALL get_next ( mc )
                 CALL get_next ( rcut )
                 inp_info%nonbonded(np) %type = string2
                 inp_info%nonbonded(np) %at1 = at1
                 inp_info%nonbonded(np) %at2 = at2
                 inp_info%nonbonded(np) %goodwin%vr0 = vr0
                 inp_info%nonbonded(np) %goodwin%d = d
                 inp_info%nonbonded(np) %goodwin%dc = dc
                 inp_info%nonbonded(np) %goodwin%m = m
                 inp_info%nonbonded(np) %goodwin%mc = mc
                 inp_info%nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
  END IF

END SUBROUTINE input_nonbonded

!!*****
!******************************************************************************
!!****** force_fields/input_charges [1.0] *
!!
!!   NAME
!!     input_charges
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I charges
!!I    ATMNAME1 qeff
!!I end charges                                                               I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_charges(ff_type,globenv)

    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(input_info_type),POINTER            :: inp_info
    CHARACTER(LEN=default_string_length)              :: label
    CHARACTER(LEN=default_string_length)              :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

  inp_info => ff_type%inp_info

!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror == 0 ) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 7
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'CHARGES') /= 0 ) THEN
           i = 0
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              IF(string2=="END") EXIT

              i = i + 1
              CALL reallocate(inp_info%charge_atm,1,i)
              CALL reallocate(inp_info%charge,1,i)

              READ(string2(:),*) inp_info%charge_atm(i)
              CALL get_next ( inp_info%charge(i) )
           END DO
        END IF
        CALL read_line
     END DO
  END IF

END SUBROUTINE input_charges

!!*****
!******************************************************************************
!!****** force_fields/input_bonds [1.0] *
!!
!!   NAME
!!     input_bonds
!!
!!   FUNCTION
!!     Read all the bond parameters. Put them in the bond_parm table.
!!     While we read the bond_pot_type,
!!     The for we allow for is :
!!     k_2*(r-r0)**2/2+k_3*(r-r0)**3/3+k_4*(r-r0)**4/4
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Bond parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I bonds
!!I    harmonic ATMNAE1 ATMNAME2 k1 r0
!!I    quartic  ATMNAE1 ATMNAME2 k1 k2 k3 r0
!!I end bonds                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bonds(ff_type,globenv)

    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(input_info_type),POINTER            :: inp_info
    CHARACTER(LEN=default_string_length)              :: label
    CHARACTER(LEN=default_string_length)              :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

  inp_info => ff_type%inp_info

!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror == 0 ) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 5
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'BONDS') /= 0 ) THEN
           i = 0
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_bond','unknown option')
              CASE ( 'HARMONIC')
                 i = i + 1
                 CALL reallocate(inp_info%bond_a,1,i)
                 CALL reallocate(inp_info%bond_b,1,i)
                 CALL reallocate(inp_info%bond_k,1,3,1,i)
                 CALL reallocate(inp_info%bond_r0,1,i)
                 ilen = default_string_length
                 CALL get_next ( inp_info%bond_a(i), ilen )
                 CALL get_next ( inp_info%bond_b(i), ilen )
                 CALL get_next ( inp_info%bond_k(1,i) )
                 inp_info%bond_k(2,i)=0.0_dp
                 inp_info%bond_k(3,i)=0.0_dp
                 CALL get_next ( inp_info%bond_r0(i) )
              CASE ( 'QUARTIC' )
                 i = i + 1
                 CALL reallocate(inp_info%bond_a,1,i)
                 CALL reallocate(inp_info%bond_b,1,i)
                 CALL reallocate(inp_info%bond_k,1,3,1,i)
                 CALL reallocate(inp_info%bond_r0,1,i)
                 ilen = default_string_length
                 CALL get_next ( inp_info%bond_a(i), ilen )
                 CALL get_next ( inp_info%bond_b(i), ilen )
                 CALL get_next ( inp_info%bond_k(1,i) )
                 CALL get_next ( inp_info%bond_k(2,i) )
                 CALL get_next ( inp_info%bond_k(3,i) )
                 CALL get_next ( inp_info%bond_r0(i) )
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
  END IF

END SUBROUTINE input_bonds

!!*****
!******************************************************************************
!!****** force_fields/input_bends [1.0] *
!!
!!   NAME
!!     input_bends
!!
!!   FUNCTION
!!     Read all the bend parameters. Put them in the bend_parm table.
!!     While we read the bend_pot_type,
!!     we (so far) will do nothing with this information. It appears that
!!     all major parameters sets, such as CHARMM or AMBER, all have only
!!     only one type of bend interaction (harmonic).
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Bend parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I bends
!!I    harmonic ATMNAE1 ATMNAME2 ATMNAE3 k theta
!!I end bends                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bends(ff_type,globenv)

    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(input_info_type),POINTER            :: inp_info
    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw
    REAL(KIND=dp)                                :: todeg, torad

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

  inp_info => ff_type%inp_info

  torad = pi / 180.0_dp
  todeg = 180.0_dp/pi

!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror == 0 ) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 5
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'BENDS') /= 0 ) THEN
           i=0
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_bond','unknown option')
              CASE ( 'HARMONIC')
                 i = i + 1
                 CALL reallocate(inp_info%bend_a,1,i)
                 CALL reallocate(inp_info%bend_b,1,i)
                 CALL reallocate(inp_info%bend_c,1,i)
                 CALL reallocate(inp_info%bend_k,1,i)
                 CALL reallocate(inp_info%bend_theta0,1,i)

                 ilen = default_string_length
                 CALL get_next ( inp_info%bend_a(i), ilen )
                 CALL get_next ( inp_info%bend_b(i), ilen )
                 CALL get_next ( inp_info%bend_c(i), ilen )
                 CALL get_next ( inp_info%bend_k(i) )
                 CALL get_next ( inp_info%bend_theta0(i) )
               !  CALL convert_to_cp2k_units ( "KELVIN",  &
               !                 bend_k = inp_info%bend_k ( i ) )
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
  END IF

END SUBROUTINE input_bends

!!*****
!******************************************************************************
!!****** force_fields/input_torsions [1.0] *
!!
!!   NAME
!!     input_torsions
!!
!!   FUNCTION
!!     Read all the torsion parameters. Put them in the torsion_parm table.
!!     We have defined the
!!     torsion_parameter type with 7 coefficients for 0 to 6. Hopefully
!!     this will be generic enough to handle all torsion potentials.
!!     Since different types will be used, we need to store the type of torsion
!!     as well.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Torsion parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I torsions
!!I    jorgensen tp[1-4]                                                      I
!!I end torsions                                                              I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_torsions(ff_type,globenv)

    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(input_info_type),POINTER            :: inp_info

    CHARACTER(LEN=default_string_length)  :: label
    CHARACTER(LEN=default_string_length)  :: type
    CHARACTER(LEN=default_string_length)  :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw, j

  LOGICAL                                  :: PM1, PM2, PM3, PM4

!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF
!------------------------------------------------------------------------------

  inp_info => ff_type%inp_info

!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror == 0 ) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
        ilen = 8
        CALL get_next(string,ilen)
        CALL uppercase ( string )
        IF (index(string,'TORSIONS') /= 0 ) THEN
           i = 0
           DO
              CALL read_line
              ilen = 0
              CALL get_next(string2,ilen)
              CALL uppercase(string2)
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_torsion','unknown option')
              CASE ( 'JORGENSEN')
                 i = i + 1
                 CALL reallocate(inp_info%torsion_a,1,i)
                 CALL reallocate(inp_info%torsion_b,1,i)
                 CALL reallocate(inp_info%torsion_c,1,i)
                 CALL reallocate(inp_info%torsion_d,1,i)
                 CALL reallocate(inp_info%torsion_k,1,i)
                 CALL reallocate(inp_info%torsion_cosphi0,1,i)

                 ilen = default_string_length
                 CALL get_next ( inp_info%torsion_a(i), ilen )
                 CALL get_next ( inp_info%torsion_b(i), ilen )
                 CALL get_next ( inp_info%torsion_c(i), ilen )
                 CALL get_next ( inp_info%torsion_d(i), ilen )
                 CALL get_next ( inp_info%torsion_k(i) )
                 CALL get_next ( inp_info%torsion_cosphi0(i) )
              !   CALL convert_to_cp2k_units ( "KELVIN",  &
              !                     torsion_a0 = inp_info%torsion_a ( i ), &
              !                     torsion_a1 = inp_info%torsion_b ( i ), &
              !                     torsion_a2 = inp_info%torsion_c ( i ), &
              !                     torsion_a3 = inp_info%torsion_d ( i ), &
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line()
     END DO
  END IF

END SUBROUTINE input_torsions


!!*****
!******************************************************************************

!******************************************************************************
!!****** force_fields/read_force_field_charmm [1.0] *
!!
!!   NAME
!!     read_force_field_charmm
!!
!!   FUNCTION
!!     Reads the charmm force_field
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_force_field_charmm ( ff_type , globenv )


    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(charmm_info_type),POINTER           :: chm_info

    LOGICAL                                  :: PM1, PM2, PM3, PM4
    INTEGER                                  :: i, istat, ilen, iw
    INTEGER                                  :: nbond,nbend,nub,ntorsion,&
                                                nimpr,nonfo,nnonbond,dummy

    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length)     :: string,string2,string3,string4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "    Entering read_force_field_charmm"

  chm_info => ff_type%chm_info
  CALL parser_init(ff_type%ff_file_name,globenv)

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Read in all the Bonds info from the param file here
  !      Vbond = Kb(b-b0)^2
  !      UNITS for Kb: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
  !      FACTOR of "2" rolled into Kb
  !-----------------------------------------------------------------------------
  nbond = 0
  label = 'BONDS'
  CALL search_label ( label, istat, ignore_case=.TRUE. )
  IF (istat == 0) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 7
       CALL get_next(string,ilen)
       CALL uppercase ( string )
       IF(string == "ANGLES") EXIT
       CALL get_next(string2,ilen)
       nbond = nbond + 1
       CALL reallocate(chm_info%bond_a,1,nbond)
       CALL reallocate(chm_info%bond_b,1,nbond)
       CALL reallocate(chm_info%bond_k,1,3,1,nbond)
       CALL reallocate(chm_info%bond_r0,1,nbond)
       chm_info%bond_a(nbond) = string
       chm_info%bond_b(nbond) = string2
       CALL get_next(chm_info%bond_k(1,nbond))
       CALL get_next(chm_info%bond_r0(nbond))
       IF(PM2) WRITE(iw,*) "    CHM BOND ",nbond,&
                           TRIM(chm_info%bond_a(nbond))," ",&
                           TRIM(chm_info%bond_b(nbond))," ",&
                           chm_info%bond_k(1,nbond),&
                           chm_info%bond_r0(nbond)
       ! Do some units conversion into internal atomic units
       CALL convert_to_cp2k_units("ANGSTROM",length=chm_info%bond_r0(nbond))
       CALL convert_to_cp2k_units("KCAL",energy=chm_info%bond_k(1,nbond))
       CALL convert_to_cp2k_units("BOHR",length=chm_info%bond_k(1,nbond))
       CALL convert_to_cp2k_units("BOHR",length=chm_info%bond_k(1,nbond))
       chm_info%bond_k(1,nbond) = 2.0_dp * chm_info%bond_k(1,nbond)
       chm_info%bond_k(2,nbond) =0.0_dp
       chm_info%bond_k(3,nbond) =0.0_dp
       CALL read_line
     END DO
  END IF

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in all the Bends and UB info from the param file here
  !      Vangle = Ktheta(theta-theta0)^2
  !      UNITS for Ktheta: [(kcal/mol)/(rad^2)] to [Eh/(rad^2)]
  !      FACTOR of "2" rolled into Ktheta
  !      Vub = Kub(S-S0)^2
  !      UNITS for Kub: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
  !      FACTOR of "2" rolled into Kub
  !-----------------------------------------------------------------------------
  nbend = 0
  nub = 0
  label = 'ANGLES'
  CALL search_label ( label, istat, ignore_case=.TRUE. )
  IF (istat == 0) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 7
       CALL get_next(string,ilen)
       CALL uppercase ( string )
       IF(string == "DIHEDRA") EXIT
       CALL get_next(string2,ilen)
       CALL get_next(string3,ilen)
       nbend = nbend + 1
       CALL reallocate(chm_info%bend_a,1,nbend)
       CALL reallocate(chm_info%bend_b,1,nbend)
       CALL reallocate(chm_info%bend_c,1,nbend)
       CALL reallocate(chm_info%bend_k,1,nbend)
       CALL reallocate(chm_info%bend_theta0,1,nbend)
       chm_info%bend_a(nbend) = string
       chm_info%bend_b(nbend) = string2
       chm_info%bend_c(nbend) = string3
       CALL get_next(chm_info%bend_k(nbend))
       CALL get_next(chm_info%bend_theta0(nbend))
       IF(PM2) WRITE(iw,*) "    CHM BEND ",nbend,&
                           TRIM(chm_info%bend_a(nbend))," ",&
                           TRIM(chm_info%bend_b(nbend))," ",&
                           TRIM(chm_info%bend_c(nbend))," ",&
                           chm_info%bend_k(nbend),&
                           chm_info%bend_theta0(nbend)
       ! Do some units conversion into internal atomic units
       CALL convert_to_cp2k_units("DEGREE",angle=chm_info%bend_theta0(nbend))
       CALL convert_to_cp2k_units("KCAL",energy=chm_info%bend_k(nbend))
       chm_info%bend_k(nbend) = 2.0_dp * chm_info%bend_k(nbend)
       IF (test_next()=='N') THEN
         nub = nub + 1
         CALL reallocate(chm_info%ub_a,1,nub)
         CALL reallocate(chm_info%ub_b,1,nub)
         CALL reallocate(chm_info%ub_c,1,nub)
         CALL reallocate(chm_info%ub_k,1,nub)
         CALL reallocate(chm_info%ub_r0,1,nub)
         chm_info%ub_a(nub) = string
         chm_info%ub_b(nub) = string2
         chm_info%ub_c(nub) = string3
         CALL get_next(chm_info%ub_k(nub))
         CALL get_next(chm_info%ub_r0(nub))
         IF(PM2) WRITE(iw,*) "    CHM UB ",nub,&
                             TRIM(chm_info%ub_a(nub))," ",&
                             TRIM(chm_info%ub_b(nub))," ",&
                             TRIM(chm_info%ub_c(nub))," ",&
                             chm_info%ub_k(nub),&
                             chm_info%ub_r0(nub)
         ! Do some units conversion into internal atomic units
         CALL convert_to_cp2k_units("ANGSTROM",length=chm_info%ub_r0(nub))
         CALL convert_to_cp2k_units("KCAL",energy=chm_info%ub_k(nub))
         CALL convert_to_cp2k_units("BOHR",length=chm_info%ub_k(nub))
         CALL convert_to_cp2k_units("BOHR",length=chm_info%ub_k(nub))
         chm_info%ub_k(nub) = 2.0_dp * chm_info%ub_k(nub)
       END IF
       CALL read_line
     END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Read in all the Dihedrals info from the param file here
  !      Vtorsion = Kphi(1+COS(n(phi)-delta))
  !      UNITS for Kphi: [(kcal/mol)] to [Eh]
  !-----------------------------------------------------------------------------
  ntorsion = 0
  label = 'DIHEDRALS'
  CALL search_label ( label, istat, ignore_case=.TRUE. )
  IF (istat == 0) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 7
       CALL get_next(string,ilen)
       CALL uppercase ( string )
       IF(string == "IMPROPE") EXIT
       CALL get_next(string2,ilen)
       CALL get_next(string3,ilen)
       CALL get_next(string4,ilen)
       ntorsion = ntorsion + 1
       CALL reallocate(chm_info%torsion_a,1,ntorsion)
       CALL reallocate(chm_info%torsion_b,1,ntorsion)
       CALL reallocate(chm_info%torsion_c,1,ntorsion)
       CALL reallocate(chm_info%torsion_d,1,ntorsion)
       CALL reallocate(chm_info%torsion_k,1,ntorsion)
       CALL reallocate(chm_info%torsion_m,1,ntorsion)
       CALL reallocate(chm_info%torsion_cosphi0,1,ntorsion)
       chm_info%torsion_a(ntorsion) = string
       chm_info%torsion_b(ntorsion) = string2
       chm_info%torsion_c(ntorsion) = string3
       chm_info%torsion_d(ntorsion) = string4
       CALL get_next(chm_info%torsion_k(ntorsion))
       CALL get_next(chm_info%torsion_m(ntorsion))
       CALL get_next(chm_info%torsion_cosphi0(ntorsion))
       IF(PM2) WRITE(iw,*) "    CHM TORSION ",ntorsion,&
                           TRIM(chm_info%torsion_a(ntorsion))," ",&
                           TRIM(chm_info%torsion_b(ntorsion))," ",&
                           TRIM(chm_info%torsion_c(ntorsion))," ",&
                           TRIM(chm_info%torsion_d(ntorsion))," ",&
                           chm_info%torsion_k(ntorsion),&
                           chm_info%torsion_m(ntorsion),&
                           chm_info%torsion_cosphi0(ntorsion)
       ! Do some units conversion into internal atomic units
       CALL convert_to_cp2k_units("DEGREE",&
                                  angle=chm_info%torsion_cosphi0(ntorsion))
       CALL convert_to_cp2k_units("KCAL",energy=chm_info%torsion_k(ntorsion))
       chm_info%torsion_cosphi0(ntorsion)=&
                                  COS(chm_info%torsion_cosphi0(ntorsion))
       CALL read_line
     END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Read in all the Improper info from the param file here
  !      Vimpr = Kpsi(psi-psi0)^2
  !      UNITS for Kpsi: [(kcal/mol)/(rad^2)] to [Eh/(rad^2)]
  !      FACTOR of "2" rolled into Kpsi
  !-----------------------------------------------------------------------------
  nimpr = 0
  label = 'IMPROPER'
  CALL search_label ( label, istat, ignore_case=.TRUE. )
  IF (istat == 0) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 7
       CALL get_next(string,ilen)
       CALL uppercase ( string )
       IF(string == "NONBOND") EXIT
       CALL get_next(string2,ilen)
       CALL get_next(string3,ilen)
       CALL get_next(string4,ilen)
       nimpr = nimpr + 1
       CALL reallocate(chm_info%impr_a,1,nimpr)
       CALL reallocate(chm_info%impr_b,1,nimpr)
       CALL reallocate(chm_info%impr_c,1,nimpr)
       CALL reallocate(chm_info%impr_d,1,nimpr)
       CALL reallocate(chm_info%impr_k,1,nimpr)
       CALL reallocate(chm_info%impr_phi0,1,nimpr)
       chm_info%impr_a(nimpr) = string
       chm_info%impr_b(nimpr) = string2
       chm_info%impr_c(nimpr) = string3
       chm_info%impr_d(nimpr) = string4
       CALL get_next(chm_info%impr_k(nimpr))
       CALL get_next(dummy)
       CALL get_next(chm_info%impr_phi0(nimpr))
       IF(PM2) WRITE(iw,*) "    CHM IMPROPERS ",nimpr,&
                           TRIM(chm_info%impr_a(nimpr))," ",&
                           TRIM(chm_info%impr_b(nimpr))," ",&
                           TRIM(chm_info%impr_c(nimpr))," ",&
                           TRIM(chm_info%impr_d(nimpr))," ",&
                           chm_info%impr_k(nimpr),&
                           chm_info%impr_phi0(nimpr)
       ! Do some units conversion into internal atomic units
       CALL convert_to_cp2k_units("DEGREE",angle=chm_info%impr_phi0(nimpr))
       CALL convert_to_cp2k_units("KCAL",energy=chm_info%impr_k(nimpr))
       chm_info%impr_k(nimpr) = 2.0_dp * chm_info%impr_k(nimpr)
       CALL read_line
     END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Read in all the Nonbonded info from the param file here
  !-----------------------------------------------------------------------------
  nnonbond = 0
  nonfo = 0
  label = 'NONBONDED'
  CALL search_label ( label, istat )
  IF (istat == 0) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 7
       CALL get_next(string,ilen)
       CALL uppercase ( string )
       IF(string == "HBOND".OR.string == "END") EXIT
       nnonbond = nnonbond + 1
       CALL reallocate(chm_info%nonbond_a,1,nnonbond)
       CALL reallocate(chm_info%nonbond_eps,1,nnonbond)
       CALL reallocate(chm_info%nonbond_rmin2,1,nnonbond)
       chm_info%nonbond_a(nnonbond) = string
       CALL get_next(chm_info%nonbond_eps(nnonbond))
       CALL get_next(chm_info%nonbond_eps(nnonbond))
       CALL get_next(chm_info%nonbond_rmin2(nnonbond))
       IF(PM2) WRITE(iw,*) "    CHM NONBOND ",nnonbond,&
                           TRIM(chm_info%nonbond_a(nnonbond))," ",&
                           chm_info%nonbond_eps(nnonbond),&
                           chm_info%nonbond_rmin2(nnonbond)
       CALL convert_to_cp2k_units("ANGSTROM",&
                                  length=chm_info%nonbond_rmin2(nnonbond))
       CALL convert_to_cp2k_units("KCAL",&
                                  energy=chm_info%nonbond_eps(nnonbond))
       IF (test_next()=='N') THEN
         nonfo = nonfo + 1
         CALL reallocate(chm_info%nonbond_a_14,1,nonfo)
         CALL reallocate(chm_info%nonbond_eps_14,1,nonfo)
         CALL reallocate(chm_info%nonbond_rmin2_14,1,nonfo)
         chm_info%nonbond_a_14(nonfo) = chm_info%nonbond_a(nnonbond)
         CALL get_next(chm_info%nonbond_eps_14(nonfo))
         CALL get_next(chm_info%nonbond_eps_14(nonfo))
         CALL get_next(chm_info%nonbond_rmin2_14(nonfo))
         IF(PM2) WRITE(iw,*) "    CHM ONFO ",nonfo,&
                             TRIM(chm_info%nonbond_a_14(nonfo))," ",&
                             chm_info%nonbond_eps_14(nonfo),&
                             chm_info%nonbond_rmin2_14(nonfo)
         CALL convert_to_cp2k_units("ANGSTROM",&
                                  length=chm_info%nonbond_rmin2_14(nonfo))
         CALL convert_to_cp2k_units("KCAL",&
                                  energy=chm_info%nonbond_eps_14(nonfo))
       END IF
       CALL read_line
     END DO
  END IF

  CALL parser_end

  IF(PM1) WRITE(iw,*) "    Exiting  read_force_field_charmm"

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE read_force_field_charmm


!******************************************************************************
!!****** force_fields/read_force_field_amber [1.0] *
!!
!!   NAME
!!     read_force_field_amber
!!
!!   FUNCTION
!!     Reads the amber force_field
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!
!!
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE read_force_field_amber ( ff_type , globenv )


    TYPE(force_field_type), &
      INTENT(INOUT)                          :: ff_type
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(input_info_type),POINTER            :: inp_info

  LOGICAL                                  :: PM1, PM2, PM3, PM4
    INTEGER                                  :: i, istat, ilen, iw

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Do something
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Do something else
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Hopefully Do the last thing
  !-----------------------------------------------------------------------------


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE read_force_field_amber



!******************************************************************************
!!****** force_field_util/force_field_pack [1.0] *
!!
!!   NAME
!!     force_field_pack
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Pack in all the information needed for the force_field
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

SUBROUTINE force_field_pack (particle_set,atomic_kind_set,&
                             molecule_kind_set,molecule_set,&
                             ewald_env,fist_nonbond_env,ff_type,globenv,&
                             qmmm, qmmm_env, error)

  TYPE ( ewald_environment_type ), POINTER :: ewald_env
  TYPE ( fist_nonbond_env_type ), POINTER ::  fist_nonbond_env
  TYPE(particle_type), DIMENSION(:), POINTER             :: particle_set
  TYPE(atomic_kind_type), DIMENSION(:), POINTER          :: atomic_kind_set
  TYPE(molecule_kind_type), DIMENSION(:), POINTER        :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), POINTER             :: molecule_set
  TYPE(force_field_type), INTENT(INOUT)                  :: ff_type
  TYPE(global_environment_type), INTENT(IN), OPTIONAL    :: globenv
  TYPE(pair_potential_type), POINTER                     :: potparm_bond(:,:)
  TYPE(pair_potential_type), POINTER                     :: potparm_nonbond(:,:)
  TYPE(cp_error_type), INTENT (INOUT), &
    OPTIONAL                                    :: error
  LOGICAL, INTENT(IN), OPTIONAL :: qmmm
  TYPE (qmmm_env_mm_type), POINTER, OPTIONAL :: qmmm_env

  TYPE(spline_environment_type), POINTER                 :: spline_bond_env
  TYPE(spline_environment_type), POINTER                 :: spline_nonbond_env
  TYPE(spline_environment_type), POINTER                 :: spline_manybody_env

  CHARACTER(len=*), PARAMETER :: routineN = 'force_field_pack', &
    routineP = moduleN//':'//routineN

  TYPE(input_info_type),POINTER     :: inp_info
  TYPE(charmm_info_type),POINTER    :: chm_info

  INTEGER                                  :: istat, iw, i, j, k, ikind, jkind, n
  LOGICAL                                  :: PM1, PM2, PM3, PM4

  CHARACTER(LEN=default_string_length)   :: name
  TYPE(atomic_kind_type), POINTER        :: atomic_kind
  TYPE(molecule_kind_type), POINTER      :: molecule_kind

  INTEGER                                :: natom,nbond,nbend,nub,ntorsion,&
                                            nonfo,nimpr,nkinds
  INTEGER                                :: atm_a,atm_b,atm_c,atm_d
  CHARACTER(LEN=default_string_length)   :: name_atm_a,name_atm_b, name_atm_a_local, name_atm_b_local
  CHARACTER(LEN=default_string_length)   :: name_atm_c,name_atm_d
  CHARACTER(LEN=default_string_length)   :: name_atm_a2,name_atm_b2
  CHARACTER(LEN=default_string_length)   :: name_atm_c2,name_atm_d2
  TYPE(atom_type), DIMENSION(:), POINTER            :: atom_list
  TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
  TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
  TYPE(torsion_type), DIMENSION(:), POINTER         :: torsion_list
  TYPE(onfo_type), DIMENSION(:), POINTER            :: onfo_list
  TYPE(ub_type), DIMENSION(:), POINTER              :: ub_list
  TYPE(impr_type), DIMENSION(:), POINTER            :: impr_list
  LOGICAL                                  :: found,found_a,found_b
  INTEGER                                  :: counter
  INTEGER,POINTER                          :: map_bond_kind(:)
  INTEGER,POINTER                          :: map_bend_kind(:)
  INTEGER,POINTER                          :: map_ub_kind(:)
  INTEGER,POINTER                          :: map_torsion_kind(:)
  INTEGER,POINTER                          :: map_impr_kind(:)
  INTEGER,POINTER                          :: map_onfo_kind(:)
  TYPE(bond_kind_type), DIMENSION(:), POINTER       :: bond_kind_set
  TYPE(bend_kind_type), DIMENSION(:), POINTER       :: bend_kind_set
  TYPE(ub_kind_type), DIMENSION(:), POINTER         :: ub_kind_set
  TYPE(onfo_kind_type), DIMENSION(:), POINTER       :: onfo_kind_set
  TYPE(impr_kind_type), DIMENSION(:), POINTER       :: impr_kind_set
  TYPE(torsion_kind_type), DIMENSION(:), POINTER    :: torsion_kind_set

  REAL(KIND=dp)                                :: rmin,epsilon, default_rcutsq
  REAL(KIND=dp)                                :: rmin2_a,rmin2_b
  REAL(KIND=dp)                                :: epsilon_a,epsilon_b,scale_factor
  INTEGER                                  :: ii, jj, imul, ilink
  LOGICAL :: only_qm, my_qmmm, is_link_atom, failure
  INTEGER, DIMENSION(:), POINTER :: my_atom_list

  CHARACTER(LEN=default_string_length)     :: atmname
  REAL(KIND=dp), DIMENSION(:), POINTER         :: charge
  TYPE(elp_potential_type), POINTER        :: elp_potential

  INTEGER, DIMENSION(:), POINTER           :: molecule_list
  INTEGER                                  :: first,last
  TYPE(molecule_type), POINTER        :: molecule
  TYPE ( spline_data_type ), POINTER :: spline_data

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  NULLIFY(atomic_kind,molecule_kind,atom_list,bond_list,bend_list,&
       torsion_list,onfo_list,map_bond_kind,map_bend_kind,map_ub_kind,&
       map_torsion_kind,bond_kind_set,bend_kind_set,ub_kind_set,&
       onfo_kind_set,impr_kind_set,torsion_kind_set,my_atom_list)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  only_qm = .FALSE.
  my_qmmm = .FALSE.
  failure = .FALSE.
  default_rcutsq = (10.0_dp*bohr)*(10.0_dp*bohr)
  IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering force_field_pack"

  inp_info => ff_type%inp_info
  chm_info => ff_type%chm_info


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Determine the number of unique bond kind and allocate bond_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nbond=nbond,bond_list=bond_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    IF(nbond>0) THEN
      ALLOCATE(map_bond_kind(nbond),STAT=istat)
      IF(istat/=0) CALL stop_memory('force_field_pack','map_bond_kind')
      counter=0
      DO j=1,nbond
        atm_a = bond_list(j)%a
        !atomic_kind => atom_list(atm_a)%atomic_kind
        atomic_kind => particle_set(atm_a+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = bond_list(j)%b
        !atomic_kind => atom_list(atm_b)%atomic_kind
        atomic_kind => particle_set(atm_b+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        found = .FALSE.
        DO k=1,j-1
          atm_a = bond_list(k)%a
          !atomic_kind => atom_list(atm_a)%atomic_kind
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = bond_list(k)%b
          !atomic_kind => atom_list(atm_b)%atomic_kind
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_bond_kind(j) = map_bond_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_bond_kind(j) = counter
        END IF
      END DO
      NULLIFY(bond_kind_set)
      CALL allocate_bond_kind_set(bond_kind_set,counter)
      DO j=1,nbond
        bond_list(j)%bond_kind => bond_kind_set(map_bond_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             bond_kind_set=bond_kind_set,bond_list=bond_list)
      DEALLOCATE(map_bond_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_bond_kind')
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Determine the number of unique bend kind and allocate bend_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nbend=nbend,bend_list=bend_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    IF(nbend>0) THEN
      ALLOCATE(map_bend_kind(nbend),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_bend_kind')
      counter=0
      DO j=1,nbend
        atm_a = bend_list(j)%a
        !atomic_kind => atom_list(atm_a)%atomic_kind
        atomic_kind => particle_set(atm_a+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = bend_list(j)%b
        !atomic_kind => atom_list(atm_b)%atomic_kind
        atomic_kind => particle_set(atm_b+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        atm_c = bend_list(j)%c
        !atomic_kind => atom_list(atm_c)%atomic_kind
        atomic_kind => particle_set(atm_c+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_c)
        found = .FALSE.
        DO k=1,j-1
          atm_a = bend_list(k)%a
          !atomic_kind => atom_list(atm_a)%atomic_kind
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = bend_list(k)%b
          !atomic_kind => atom_list(atm_b)%atomic_kind
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          atm_c = bend_list(k)%c
          !atomic_kind => atom_list(atm_c)%atomic_kind
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_c2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_c2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_bend_kind(j) = map_bend_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_bend_kind(j) = counter
        END IF
      END DO
      CALL allocate_bend_kind_set(bend_kind_set,counter)
      DO j=1,nbend
        bend_list(j)%bend_kind => bend_kind_set(map_bend_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             bend_kind_set=bend_kind_set,bend_list=bend_list)
      DEALLOCATE(map_bend_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_bend_kind')
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Determine the number of unique Urey-Bradley kind and allocate ub_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nub=nub,ub_list=ub_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    IF(nub>0) THEN
      ALLOCATE(map_ub_kind(nub),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_ub_kind')
      counter=0
      DO j=1,nub
        atm_a = ub_list(j)%a
        !atomic_kind => atom_list(atm_a)%atomic_kind
        atomic_kind => particle_set(atm_a+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = ub_list(j)%b
        !atomic_kind => atom_list(atm_b)%atomic_kind
        atomic_kind => particle_set(atm_b+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        atm_c = ub_list(j)%c
        !atomic_kind => atom_list(atm_c)%atomic_kind
        atomic_kind => particle_set(atm_c+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_c)
        found = .FALSE.
        DO k=1,j-1
          atm_a = ub_list(k)%a
          !atomic_kind => atom_list(atm_a)%atomic_kind
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = ub_list(k)%b
          !atomic_kind => atom_list(atm_b)%atomic_kind
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          atm_c = ub_list(k)%c
          !atomic_kind => atom_list(atm_c)%atomic_kind
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_c2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_c2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_ub_kind(j) = map_ub_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_ub_kind(j) = counter
        END IF
      END DO
      CALL allocate_ub_kind_set(ub_kind_set,counter)
      DO j=1,nub
        ub_list(j)%ub_kind => ub_kind_set(map_ub_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             ub_kind_set=ub_kind_set,ub_list=ub_list)
      DEALLOCATE(map_ub_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_ub_kind')
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Determine the number of unique torsion kind and allocate torsion_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           ntorsion=ntorsion,torsion_list=torsion_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    IF(ntorsion>0) THEN
      ALLOCATE(map_torsion_kind(ntorsion),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_torsion_kind')
      counter=0
      DO j=1,ntorsion
        name_atm_a = ''
        name_atm_b = ''
        name_atm_c = ''
        name_atm_d = ''
        name_atm_a2 = ''
        name_atm_b2 = ''
        name_atm_c2 = ''
        name_atm_d2 = ''
        atm_a = torsion_list(j)%a
        !atomic_kind => atom_list(atm_a)%atomic_kind
        atomic_kind => particle_set(atm_a+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = torsion_list(j)%b
        !atomic_kind => atom_list(atm_b)%atomic_kind
        atomic_kind => particle_set(atm_b+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        atm_c = torsion_list(j)%c
        !atomic_kind => atom_list(atm_c)%atomic_kind
        atomic_kind => particle_set(atm_c+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_c)
        atm_d = torsion_list(j)%d
        !atomic_kind => atom_list(atm_d)%atomic_kind
        atomic_kind => particle_set(atm_d+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_d)
        found = .FALSE.
        DO k=1,j-1
          atm_a = torsion_list(k)%a
          !atomic_kind => atom_list(atm_a)%atomic_kind
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = torsion_list(k)%b
          !atomic_kind => atom_list(atm_b)%atomic_kind
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          atm_c = torsion_list(k)%c
          !atomic_kind => atom_list(atm_c)%atomic_kind
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_c2)
          atm_d = torsion_list(k)%d
          !atomic_kind => atom_list(atm_d)%atomic_kind
          atomic_kind => particle_set(atm_d+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_d2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_d)==TRIM(name_atm_d2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_d2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_d)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_torsion_kind(j) = map_torsion_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_torsion_kind(j) = counter
        END IF
      END DO
      CALL allocate_torsion_kind_set(torsion_kind_set,counter)
      DO j=1,ntorsion
        torsion_list(j)%torsion_kind => torsion_kind_set(map_torsion_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             torsion_kind_set=torsion_kind_set,torsion_list=torsion_list)
      DEALLOCATE(map_torsion_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_torsion_kind')
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Determine the number of unique impr kind and allocate impr_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nimpr=nimpr,impr_list=impr_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    IF(nimpr>0) THEN
      ALLOCATE(map_impr_kind(nimpr),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_impr_kind')
      counter=0
      DO j=1,nimpr
        name_atm_a = ''
        name_atm_b = ''
        name_atm_c = ''
        name_atm_d = ''
        name_atm_a2 = ''
        name_atm_b2 = ''
        name_atm_c2 = ''
        name_atm_d2 = ''
        atm_a = impr_list(j)%a
        !atomic_kind => atom_list(atm_a)%atomic_kind
        atomic_kind => particle_set(atm_a+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = impr_list(j)%b
        !atomic_kind => atom_list(atm_b)%atomic_kind
        atomic_kind => particle_set(atm_b+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        atm_c = impr_list(j)%c
        !atomic_kind => atom_list(atm_c)%atomic_kind
        atomic_kind => particle_set(atm_c+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_c)
        atm_d = impr_list(j)%d
        !atomic_kind => atom_list(atm_d)%atomic_kind
        atomic_kind => particle_set(atm_d+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_d)
        found = .FALSE.
        DO k=1,j-1
          atm_a = impr_list(k)%a
          !atomic_kind => atom_list(atm_a)%atomic_kind
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = impr_list(k)%b
          !atomic_kind => atom_list(atm_b)%atomic_kind
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          atm_c = impr_list(k)%c
          !atomic_kind => atom_list(atm_c)%atomic_kind
          atomic_kind => particle_set(atm_c+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_c2)
          atm_d = impr_list(k)%d
          !atomic_kind => atom_list(atm_d)%atomic_kind
          atomic_kind => particle_set(atm_d+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_d2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_d)==TRIM(name_atm_d2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_d2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_c2)) .AND. &
               (TRIM(name_atm_c)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_d)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_impr_kind(j) = map_impr_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_impr_kind(j) = counter
        END IF
      END DO
      CALL allocate_impr_kind_set(impr_kind_set,counter)
      DO j=1,nimpr
        impr_list(j)%impr_kind => impr_kind_set(map_impr_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             impr_kind_set=impr_kind_set,impr_list=impr_list)
      DEALLOCATE(map_impr_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_impr_kind')
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Determine the number of unique onfo kind and allocate onfo_kind_set
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nonfo=nonfo,onfo_list=onfo_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    IF(nonfo>0) THEN
      ALLOCATE(map_onfo_kind(nonfo),STAT=istat)
      IF (istat /= 0) CALL stop_memory('force_field_pack','map_onfo_kind')
      counter=0
      DO j=1,nonfo
        name_atm_a = ''
        name_atm_b = ''
        name_atm_a2 = ''
        name_atm_b2 = ''
        atm_a = onfo_list(j)%a
        !atomic_kind => atom_list(atm_a)%atomic_kind
        atomic_kind => particle_set(atm_a+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_a)
        atm_b = onfo_list(j)%b
        !atomic_kind => atom_list(atm_b)%atomic_kind
        atomic_kind => particle_set(atm_b+first-1)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=name_atm_b)
        found = .FALSE.
        DO k=1,j-1
          atm_a = onfo_list(k)%a
          !atomic_kind => atom_list(atm_a)%atomic_kind
          atomic_kind => particle_set(atm_a+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_a2)
          atm_b = onfo_list(k)%b
          !atomic_kind => atom_list(atm_b)%atomic_kind
          atomic_kind => particle_set(atm_b+first-1)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               name=name_atm_b2)
          IF( ((TRIM(name_atm_a)==TRIM(name_atm_a2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_b2))) .OR. &
              ((TRIM(name_atm_a)==TRIM(name_atm_b2)) .AND. &
               (TRIM(name_atm_b)==TRIM(name_atm_a2))) ) THEN
            found = .TRUE.
            map_onfo_kind(j) = map_onfo_kind(k)
            EXIT
          END IF
        END DO
        IF(.NOT.found) THEN
          counter=counter+1
          map_onfo_kind(j) = counter
        END IF
      END DO
      CALL allocate_onfo_kind_set(onfo_kind_set,counter)
      DO j=1,nonfo
        onfo_list(j)%onfo_kind => onfo_kind_set(map_onfo_kind(j))
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             onfo_kind_set=onfo_kind_set,onfo_list=onfo_list)
      DEALLOCATE(map_onfo_kind,STAT=istat)
      IF(istat/=0) CALL stop_memory('force_field_pack','map_onfo_kind')
    END IF
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. BONDS
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nbond=nbond,bond_list=bond_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    DO j=1,nbond
      atm_a = bond_list(j)%a
      !atomic_kind => atom_list(atm_a)%atomic_kind
      atomic_kind => particle_set(atm_a+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = bond_list(j)%b
      !atomic_kind => atom_list(atm_b)%atomic_kind
      atomic_kind => particle_set(atm_b+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      found = .FALSE.
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)
      ! loop over params from amber
!     IF(ASSOCIATED(ffam_info%ff_bond_a)) THEN
!       DO k=1,SIZE(ffam_info%ff_bond_a)
!         IF( ((TRIM(ffam_info%ff_bond_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(ffam_info%ff_bond_b(k))==TRIM(name_atm_b))) .OR. &
!             ((TRIM(ffam_info%ff_bond_a(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_bond_b(k))==TRIM(name_atm_a))) ) THEN
!           bond_list(j)%bond_kind%k  = ffam_info%ff_bond_k(k)
!           bond_list(j)%bond_kind%r0 = ffam_info%ff_bond_r0(k)
!         IF(found .AND. PM1) &
!            WRITE(iw,*) "    Multiple Bond declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

      ! loop over params from charmm
      IF(ASSOCIATED(chm_info%bond_a)) THEN
        DO k=1,SIZE(chm_info%bond_a)
          IF( ((TRIM(chm_info%bond_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(chm_info%bond_b(k))==TRIM(name_atm_b))) .OR. &
              ((TRIM(chm_info%bond_a(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%bond_b(k))==TRIM(name_atm_a))) ) THEN
            bond_list(j)%bond_kind%name = "CHARMM"
            bond_list(j)%bond_kind%k    = chm_info%bond_k(:,k)
            bond_list(j)%bond_kind%r0   = chm_info%bond_r0(k)
          IF(found .AND. PM1) &
             WRITE(iw,*) "    Multiple Bond declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(inp_info%bond_a)) THEN
        DO k=1,SIZE(inp_info%bond_a)
          IF( ((TRIM(inp_info%bond_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(inp_info%bond_b(k))==TRIM(name_atm_b))) .OR. &
              ((TRIM(inp_info%bond_a(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%bond_b(k))==TRIM(name_atm_a))) ) THEN
            bond_list(j)%bond_kind%name = "INPUT"
            bond_list(j)%bond_kind%k    = inp_info%bond_k(:,k)
            bond_list(j)%bond_kind%r0   = inp_info%bond_r0(k)
          IF(found .AND. PM1) &
             WRITE(iw,*) "    Multiple Bond declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) THEN
        IF(PM1) WRITE(iw,*) "    Bond ",j," : ",&
                            TRIM(name_atm_a)," ",&
                            TRIM(name_atm_b)
        CALL stop_program("FF_pack","missing bond parm")
      END IF
      !
      ! QM/MM modifications
      !
      IF (only_qm) THEN
         bond_list(j)%name           = "UNSET"
         bond_list(j)%bond_kind%name = "UNSET"
      END IF
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           bond_list=bond_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. BENDS
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nbend=nbend,bend_list=bend_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    DO j=1,nbend
      atm_a = bend_list(j)%a
      !atomic_kind => atom_list(atm_a)%atomic_kind
      atomic_kind => particle_set(atm_a+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = bend_list(j)%b
      !atomic_kind => atom_list(atm_b)%atomic_kind
      atomic_kind => particle_set(atm_b+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      atm_c = bend_list(j)%c
      !atomic_kind => atom_list(atm_c)%atomic_kind
      atomic_kind => particle_set(atm_c+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_c)
      found = .FALSE.
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c)
      ! loop over params from amber
!     IF(ASSOCIATED(ffam_info%ff_bend_a)) THEN
!       DO k=1,SIZE(ffam_info%ff_bend_a)
!         IF( ((TRIM(ffam_info%ff_bend_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(ffam_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_bend_c(k))==TRIM(name_atm_c))) .OR. &
!             ((TRIM(ffam_info%ff_bend_a(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(ffam_info%ff_bend_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_bend_c(k))==TRIM(name_atm_a))) ) THEN
!           bend_list(j)%bend_kind%k      = ffam_info%ff_bend_k(k)
!           bend_list(j)%bend_kind%theta0 = ffam_info%ff_bend_theta0(k)
!           IF(found .AND. PM1) &
!              WRITE(iw,*) "    Multiple Bend declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

      ! loop over params from charmm
      IF(ASSOCIATED(chm_info%bend_a)) THEN
        DO k=1,SIZE(chm_info%bend_a)
          IF( ((TRIM(chm_info%bend_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(chm_info%bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%bend_c(k))==TRIM(name_atm_c))) .OR. &
              ((TRIM(chm_info%bend_a(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(chm_info%bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%bend_c(k))==TRIM(name_atm_a))) ) THEN
            bend_list(j)%bend_kind%name   = "CHARMM"
            bend_list(j)%bend_kind%k      = chm_info%bend_k(k)
            bend_list(j)%bend_kind%theta0 = chm_info%bend_theta0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Bend declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(inp_info%bend_a)) THEN
        DO k=1,SIZE(inp_info%bend_a)
          IF( ((TRIM(inp_info%bend_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(inp_info%bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%bend_c(k))==TRIM(name_atm_c))) .OR. &
              ((TRIM(inp_info%bend_a(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(inp_info%bend_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%bend_c(k))==TRIM(name_atm_a))) ) THEN
            bend_list(j)%bend_kind%name   = "INPUT"
            bend_list(j)%bend_kind%k      = inp_info%bend_k(k)
            bend_list(j)%bend_kind%theta0 = inp_info%bend_theta0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Bend declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) THEN
        IF(PM1) WRITE(iw,*) "    Bend ",j," : ",&
                            TRIM(name_atm_a)," ",&
                            TRIM(name_atm_b)," ",&
                            TRIM(name_atm_c)
        CALL stop_program("FF_pack","missing bend parm")
      END IF
      !
      ! QM/MM modifications
      !
      IF (only_qm) THEN
         bend_list(j)%name           = "UNSET"
         bend_list(j)%bend_kind%name = "UNSET"
      END IF
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           bend_list=bend_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 9. Urey-Bradley
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nub=nub,ub_list=ub_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    DO j=1,nub
      atm_a = ub_list(j)%a
      !atomic_kind => atom_list(atm_a)%atomic_kind
      atomic_kind => particle_set(atm_a+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = ub_list(j)%b
      !atomic_kind => atom_list(atm_b)%atomic_kind
      atomic_kind => particle_set(atm_b+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      atm_c = ub_list(j)%c
      !atomic_kind => atom_list(atm_c)%atomic_kind
      atomic_kind => particle_set(atm_c+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_c)
      found = .FALSE.
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c)
      ! loop over params from amber
!     IF(ASSOCIATED(ffam_info%ff_ub_a)) THEN
!       DO k=1,SIZE(ffam_info%ff_ub_a)
!         IF( ((TRIM(ffam_info%ff_ub_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(ffam_info%ff_ub_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_ub_c(k))==TRIM(name_atm_c))) .OR. &
!             ((TRIM(ffam_info%ff_ub_a(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(ffam_info%ff_ub_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(ffam_info%ff_ub_c(k))==TRIM(name_atm_a))) ) THEN
!           ub_list(j)%ub_kind%k      = ffam_info%ff_ub_k(k)
!           ub_list(j)%ub_kind%theta0 = ffam_info%ff_ub_theta0(k)
!           IF(found .AND. PM1) &
!              WRITE(iw,*) "    Multiple Bend declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

      ! loop over params from charmm
      IF(ASSOCIATED(chm_info%ub_a)) THEN
        DO k=1,SIZE(chm_info%ub_a)
          IF( ((TRIM(chm_info%ub_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(chm_info%ub_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%ub_c(k))==TRIM(name_atm_c))) .OR. &
              ((TRIM(chm_info%ub_a(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(chm_info%ub_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%ub_c(k))==TRIM(name_atm_a))) ) THEN
            ub_list(j)%ub_kind%name   = "CHARMM"
            ub_list(j)%ub_kind%k      = chm_info%ub_k(k)
            ub_list(j)%ub_kind%r0     = chm_info%ub_r0(k)
            IF(PM2) WRITE(iw,*) "    Found UB ",TRIM(name_atm_a)," ",&
                                TRIM(name_atm_b)," ",TRIM(name_atm_c)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple UB declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(inp_info%ub_a)) THEN
        DO k=1,SIZE(inp_info%ub_a)
          IF( ((TRIM(inp_info%ub_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(inp_info%ub_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%ub_c(k))==TRIM(name_atm_c))) .OR. &
              ((TRIM(inp_info%ub_a(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(inp_info%ub_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%ub_c(k))==TRIM(name_atm_a))) ) THEN
            ub_list(j)%ub_kind%name   = "INPUT"
            ub_list(j)%ub_kind%k      = inp_info%ub_k(k)
            ub_list(j)%ub_kind%r0     = inp_info%ub_r0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple UB declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) THEN
        IF(PM2) WRITE(iw,*) "    Not Found UB ",TRIM(name_atm_a)," ",&
                                TRIM(name_atm_b)," ",TRIM(name_atm_c)
        ub_list(j)%ub_kind%k      = 0.0_dp
        ub_list(j)%ub_kind%r0     = 0.0_dp
      END IF
      !
      ! QM/MM modifications
      !
      IF (only_qm) THEN
         ub_list(j)%name         = "UNSET"
         ub_list(j)%ub_kind%name = "UNSET"
      END IF
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           ub_list=ub_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 10. TORSION
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           ntorsion=ntorsion,torsion_list=torsion_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    DO j=1,ntorsion
     IF(torsion_list(j)%torsion_kind%name == "UNSET") THEN
      atm_a = torsion_list(j)%a
      !atomic_kind => atom_list(atm_a)%atomic_kind
      atomic_kind => particle_set(atm_a+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = torsion_list(j)%b
      !atomic_kind => atom_list(atm_b)%atomic_kind
      atomic_kind => particle_set(atm_b+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      atm_c = torsion_list(j)%c
      !atomic_kind => atom_list(atm_c)%atomic_kind
      atomic_kind => particle_set(atm_c+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_c)
      atm_d = torsion_list(j)%d
      !atomic_kind => atom_list(atm_d)%atomic_kind
      atomic_kind => particle_set(atm_d+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_d)
      found = .FALSE.
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c, id4=name_atm_d)
      ! loop over params from amber
!     IF(ASSOCIATED(chm_info%torsion_a)) THEN
!       DO k=1,SIZE(chm_info%torsion_a)
!         IF( ((TRIM(chm_info%torsion_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(chm_info%torsion_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(chm_info%torsion_c(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(chm_info%torsion_d(k))==TRIM(name_atm_d))) .OR. &
!             ((TRIM(chm_info%torsion_a(k))==TRIM(name_atm_d)) .AND. &
!              (TRIM(chm_info%torsion_b(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(chm_info%torsion_c(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(chm_info%torsion_d(k))==TRIM(name_atm_a))) ) THEN
!           torsion_list(j)%torsion_kind%k      = chm_info%torsion_k(k)
!           torsion_list(j)%torsion_kind%phi0 = chm_info%torsion_cosphi0(k)
!           IF(found .AND. PM1) &
!              WRITE(iw,*) "    Multiple Torsion declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

      ! loop over params from charmm
      IF(ASSOCIATED(chm_info%torsion_a)) THEN
        DO k=1,SIZE(chm_info%torsion_a)
          IF( ((TRIM(chm_info%torsion_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(chm_info%torsion_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%torsion_c(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(chm_info%torsion_d(k))==TRIM(name_atm_d))) .OR. &
              ((TRIM(chm_info%torsion_a(k))==TRIM(name_atm_d)) .AND. &
               (TRIM(chm_info%torsion_b(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(chm_info%torsion_c(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%torsion_d(k))==TRIM(name_atm_a))) ) THEN
            imul = torsion_list(j)%torsion_kind%nmul + 1
            CALL reallocate(torsion_list(j)%torsion_kind%k,1,imul)
            CALL reallocate(torsion_list(j)%torsion_kind%m,1,imul)
            CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,imul)
            torsion_list(j)%torsion_kind%name    = "CHARMM"
            torsion_list(j)%torsion_kind%k(imul) = chm_info%torsion_k(k)
            torsion_list(j)%torsion_kind%m(imul) = chm_info%torsion_m(k)
            torsion_list(j)%torsion_kind%cosphi0(imul) = chm_info%torsion_cosphi0(k)
            torsion_list(j)%torsion_kind%nmul = imul
            IF(found .AND. PM1) &
               WRITE(iw,'(A,4A10)') "     Multiple Torsion for: ",&
                                    name_atm_a,name_atm_b,&
                                    name_atm_c,name_atm_d
            found = .TRUE.
          END IF
        END DO
        IF(.NOT.found) THEN
          DO k=1,SIZE(chm_info%torsion_a)
            IF( ((TRIM(chm_info%torsion_a(k))==TRIM("X")) .AND. &
                 (TRIM(chm_info%torsion_b(k))==TRIM(name_atm_b)) .AND. &
                 (TRIM(chm_info%torsion_c(k))==TRIM(name_atm_c)) .AND. &
                 (TRIM(chm_info%torsion_d(k))==TRIM("X"))) .OR. &
                ((TRIM(chm_info%torsion_a(k))==TRIM("X")) .AND. &
                 (TRIM(chm_info%torsion_b(k))==TRIM(name_atm_c)) .AND. &
                 (TRIM(chm_info%torsion_c(k))==TRIM(name_atm_b)) .AND. &
                 (TRIM(chm_info%torsion_d(k))==TRIM("X"))) ) THEN
              imul = torsion_list(j)%torsion_kind%nmul + 1
              CALL reallocate(torsion_list(j)%torsion_kind%k,1,imul)
              CALL reallocate(torsion_list(j)%torsion_kind%m,1,imul)
              CALL reallocate(torsion_list(j)%torsion_kind%cosphi0,1,imul)
              torsion_list(j)%torsion_kind%name    = "CHARMM"
              torsion_list(j)%torsion_kind%k(imul) = chm_info%torsion_k(k)
              torsion_list(j)%torsion_kind%m(imul) = chm_info%torsion_m(k)
              torsion_list(j)%torsion_kind%cosphi0(imul) = chm_info%torsion_cosphi0(k)
              torsion_list(j)%torsion_kind%nmul = imul
              IF(PM1) &
                 WRITE(iw,'(A,4A10)') "     Wildcard Torsion for: ",&
                                      name_atm_a,name_atm_b,&
                                      name_atm_c,name_atm_d
              found = .TRUE.
            END IF
          END DO
        END IF
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(inp_info%torsion_a)) THEN
        DO k=1,SIZE(inp_info%torsion_a)
          IF( ((TRIM(inp_info%torsion_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(inp_info%torsion_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%torsion_c(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(inp_info%torsion_d(k))==TRIM(name_atm_d))) .OR. &
              ((TRIM(inp_info%torsion_a(k))==TRIM(name_atm_d)) .AND. &
               (TRIM(inp_info%torsion_b(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(inp_info%torsion_c(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%torsion_d(k))==TRIM(name_atm_a))) ) THEN
            torsion_list(j)%torsion_kind%name    = "INPUT"
            torsion_list(j)%torsion_kind%k       = inp_info%torsion_k(k)
            torsion_list(j)%torsion_kind%cosphi0 = inp_info%torsion_cosphi0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Torsion declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) THEN
        IF(PM1) WRITE(iw,*) "    Torsion PARAM not found ",j," : ",&
                            TRIM(name_atm_a)," ",&
                            TRIM(name_atm_b)," ",&
                            TRIM(name_atm_c)," ",&
                            TRIM(name_atm_d)
      END IF
      !
      ! QM/MM modifications
      !
      IF (only_qm) THEN
         IF (PM1) WRITE(iw,*)"    Torsion PARAM between QM atoms ",j," : ",&
                            TRIM(name_atm_a)," ",&
                            TRIM(name_atm_b)," ",&
                            TRIM(name_atm_c)," ",&
                            TRIM(name_atm_d)," ",&
                            torsion_list(j)%a,&
                            torsion_list(j)%b,&
                            torsion_list(j)%c,&
                            torsion_list(j)%d

         torsion_list(j)%torsion_kind%name = "UNSET"
         torsion_list(j)%name              = "UNSET"
      END IF
     END IF

    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           torsion_list=torsion_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 11. IMPROPERS
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nimpr=nimpr,impr_list=impr_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    DO j=1,nimpr
      atm_a = impr_list(j)%a
      !atomic_kind => atom_list(atm_a)%atomic_kind
      atomic_kind => particle_set(atm_a+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a)
      atm_b = impr_list(j)%b
      !atomic_kind => atom_list(atm_b)%atomic_kind
      atomic_kind => particle_set(atm_b+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b)
      atm_c = impr_list(j)%c
      !atomic_kind => atom_list(atm_c)%atomic_kind
      atomic_kind => particle_set(atm_c+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_c)
      atm_d = impr_list(j)%d
      !atomic_kind => atom_list(atm_d)%atomic_kind
      atomic_kind => particle_set(atm_d+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_d)
      found = .FALSE.
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b, id3=name_atm_c, id4=name_atm_d)
      ! loop over params from amber
!     IF(ASSOCIATED(chm_info%impr_a)) THEN
!       DO k=1,SIZE(chm_info%impr_a)
!         IF( ((TRIM(chm_info%impr_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(chm_info%impr_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(chm_info%impr_c(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(chm_info%impr_d(k))==TRIM(name_atm_d))) .OR. &
!             ((TRIM(chm_info%impr_a(k))==TRIM(name_atm_d)) .AND. &
!              (TRIM(chm_info%impr_b(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(chm_info%impr_c(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(chm_info%impr_d(k))==TRIM(name_atm_a))) ) THEN
!           impr_list(j)%impr_kind%k      = chm_info%impr_k(k)
!           impr_list(j)%impr_kind%phi0 = chm_info%impr_phi0(k)
!           IF(found .AND. PM1) &
!              WRITE(iw,*) "    Multiple Impropers declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF
!
      ! loop over params from charmm
      IF(ASSOCIATED(chm_info%impr_a)) THEN
        DO k=1,SIZE(chm_info%impr_a)
          IF( ((TRIM(chm_info%impr_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(chm_info%impr_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%impr_c(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(chm_info%impr_d(k))==TRIM(name_atm_d))) .OR. &
              ((TRIM(chm_info%impr_a(k))==TRIM(name_atm_d)) .AND. &
               (TRIM(chm_info%impr_b(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(chm_info%impr_c(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(chm_info%impr_d(k))==TRIM(name_atm_a))) ) THEN
            impr_list(j)%impr_kind%name = "CHARMM"
            impr_list(j)%impr_kind%k    = chm_info%impr_k(k)
            impr_list(j)%impr_kind%phi0 = chm_info%impr_phi0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Improper declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
       END DO
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(inp_info%impr_a)) THEN
        DO k=1,SIZE(inp_info%impr_a)
          IF( ((TRIM(inp_info%impr_a(k))==TRIM(name_atm_a)) .AND. &
               (TRIM(inp_info%impr_b(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%impr_c(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(inp_info%impr_d(k))==TRIM(name_atm_d))) .OR. &
              ((TRIM(inp_info%impr_a(k))==TRIM(name_atm_d)) .AND. &
               (TRIM(inp_info%impr_b(k))==TRIM(name_atm_c)) .AND. &
               (TRIM(inp_info%impr_c(k))==TRIM(name_atm_b)) .AND. &
               (TRIM(inp_info%impr_d(k))==TRIM(name_atm_a))) ) THEN
            impr_list(j)%impr_kind%name = "INPUT"
            impr_list(j)%impr_kind%k    = inp_info%impr_k(k)
            impr_list(j)%impr_kind%phi0 = inp_info%impr_phi0(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple Improper declarations: overwriting"
            found = .TRUE.
            EXIT
          END IF
       END DO
      END IF

      IF(.NOT.found) THEN
        IF(PM2) WRITE(iw,*) "    IMPR ",j," : ",&
                            TRIM(name_atm_a)," ",&
                            TRIM(name_atm_b)," ",&
                            TRIM(name_atm_c)," ",&
                            TRIM(name_atm_d)
        !FIXME FIXME FIXME FIXME
        impr_list(j)%impr_kind%k      = 0.0_dp
        impr_list(j)%impr_kind%phi0   = 0.0_dp
        !CALL stop_program("FF_pack","missing impr parm")

      END IF
     !
     ! QM/MM modifications
     !
     IF (only_qm) THEN
        impr_list(j)%impr_kind%name = "UNSET"           
        impr_list(j)%name           = "UNSET"           
     END IF

    END DO
   CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           impr_list=impr_list)
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 12. Set spline_bond_env, use for ONFO interaction
  !-----------------------------------------------------------------------------
  CALL pot_parameters_create ( potparm_bond,  &
                               SIZE ( atomic_kind_set ), &
                               error )
  DO i=1,SIZE(atomic_kind_set)
    atomic_kind => atomic_kind_set(i)
    CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_a_local)
    DO j=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(j)
      CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_b_local)
      found = .FALSE.
      found_a = .FALSE.
      found_b = .FALSE.
      name_atm_a = name_atm_a_local
      name_atm_b = name_atm_b_local
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)

      ! loop over params from amber
      ! loop over params from amber
      ! loop over params from amber
      ! loop over params from amber

      ! loop over params from charmm
      IF(ASSOCIATED(chm_info%nonbond_a_14)) THEN
        ii = 0 
        jj = 0
        DO k=1,SIZE(chm_info%nonbond_a_14)
          IF(TRIM(name_atm_a)==TRIM(chm_info%nonbond_a_14(k))) THEN
            ii = k
            rmin2_a = chm_info%nonbond_rmin2_14(k)
            epsilon_a = chm_info%nonbond_eps_14(k)
            found_a = .TRUE.
          END IF
        END DO
        DO k=1,SIZE(chm_info%nonbond_a_14)
          IF(TRIM(name_atm_b)==TRIM(chm_info%nonbond_a_14(k))) THEN
            jj = k
            rmin2_b = chm_info%nonbond_rmin2_14(k)
            epsilon_b = chm_info%nonbond_eps_14(k)
            found_b = .TRUE.
          END IF
        END DO

        IF(.NOT.found_a) THEN
          DO k=1,SIZE(chm_info%nonbond_a)
            IF(TRIM(name_atm_a)==TRIM(chm_info%nonbond_a(k))) THEN
              ii = k
              rmin2_a = chm_info%nonbond_rmin2(k)
              epsilon_a = chm_info%nonbond_eps(k)
            END IF
          END DO
        END IF
        IF(.NOT.found_b) THEN
          DO k=1,SIZE(chm_info%nonbond_a)
            IF(TRIM(name_atm_b)==TRIM(chm_info%nonbond_a(k))) THEN
              jj = k
              rmin2_b = chm_info%nonbond_rmin2(k)
              epsilon_b = chm_info%nonbond_eps(k)
            END IF
          END DO
        END IF

        IF(ii/=0 .AND. jj/=0) THEN
          !rmin = chm_info%nonbond_rmin2_14(ii)+chm_info%nonbond_rmin2_14(jj)
          !epsilon = sqrt(chm_info%nonbond_eps_14(ii)*&
          !               chm_info%nonbond_eps_14(jj))
          rmin = rmin2_a + rmin2_b
          epsilon = sqrt(epsilon_a*epsilon_b)
          potparm_bond(i,j)%type       = 'LENNARD-JONES'
          potparm_bond(i,j)%at1        = name_atm_a
          potparm_bond(i,j)%at2        = name_atm_b
          potparm_bond(i,j)%lj%epsilon = epsilon
          potparm_bond(i,j)%lj%sigma6  = 0.5_dp*rmin**6
          potparm_bond(i,j)%lj%sigma12 = 0.25_dp*rmin**12
          IF(ff_type%rcut_nb>0.0_dp) THEN
            potparm_bond(i,j)%rcutsq     = ff_type%rcut_nb*ff_type%rcut_nb
          ELSE
            !Corresponds to an cutoff of 10A
            potparm_bond(i,j)%rcutsq     = default_rcutsq
          END IF

          IF(found .AND. PM1) &
             WRITE(iw,*) "    Multiple LJ-POT declarations: overwriting"
          found = .TRUE.
        END IF
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(inp_info%nonbonded14)) THEN
        DO k=1,SIZE(inp_info%nonbonded14)
          IF(PM2) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
              " with ",TRIM(inp_info%nonbonded14(k)%at1),&
              TRIM(inp_info%nonbonded14(k)%at2)
          IF(((TRIM(name_atm_a)==TRIM(inp_info%nonbonded14(k)%at1))  .AND. &
              (TRIM(name_atm_b)==TRIM(inp_info%nonbonded14(k)%at2))) .OR.  &
             ((TRIM(name_atm_b)==TRIM(inp_info%nonbonded14(k)%at1))  .AND. &
              (TRIM(name_atm_a)==TRIM(inp_info%nonbonded14(k)%at2))) ) THEN
            potparm_bond(i,j) = inp_info%nonbonded14(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple LJ-POT declarations: overwriting"
            IF(PM2) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a),&
                    " ",TRIM(name_atm_b)
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) THEN
        IF(PM2) WRITE(iw,*) "    spline_bond_env ",j," : ",&
                            TRIM(name_atm_a)," ",&
                            TRIM(name_atm_b)
        potparm_bond(i,j)%type       = 'UNSET'
        !Temporary Fix until 1-4 clean up is done
        potparm_bond(i,j)%type       = 'LENNARD-JONES'
        potparm_bond(i,j)%at1        = ''
        potparm_bond(i,j)%at2        = ''
        potparm_bond(i,j)%lj%epsilon = 0.0_dp
        potparm_bond(i,j)%lj%sigma6  = 0.0_dp
        potparm_bond(i,j)%lj%sigma12 = 0.0_dp
        potparm_bond(i,j)%rcutsq     = default_rcutsq
      END IF

      IF (only_qm) THEN
         CALL clean_pot_parameter( potparm_bond(i,j) )
         potparm_bond(i,j)%type       = 'LENNARD-JONES'
         potparm_bond(i,j)%rcutsq     = default_rcutsq
      END IF

    END DO
  END DO
  NULLIFY(spline_bond_env)
  nkinds = SIZE ( potparm_bond, 1 )
  CALL spline_env_create ( spline_bond_env, nkinds, error )
  CALL spline_nonbond_control ( spline_bond_env, potparm_bond, atomic_kind_set )
! point potparm ( ikind, jkind ) % spline_data => spline_data_p ( n )
  DO ikind = 1, SIZE ( potparm_bond, 1 )
    DO jkind = ikind, SIZE ( potparm_bond, 2 )
      n = spline_bond_env % spltab ( ikind, jkind )
      spline_data => spline_bond_env % spline_data_p ( n ) % spline_data
! do i,j
      CALL spline_data_retain ( spline_data, error )
      CALL spline_data_release (  &
         potparm_bond ( ikind, jkind ) % pair_spline_data, error )
      potparm_bond ( ikind, jkind ) % pair_spline_data => spline_data
! do j,i
      CALL spline_data_retain ( spline_data, error )
      CALL spline_data_release (  &
         potparm_bond ( jkind, ikind ) % pair_spline_data, error )
      potparm_bond ( jkind, ikind ) % pair_spline_data => spline_data
    END DO
  END DO

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 13. ONFO
  !-----------------------------------------------------------------------------
  DO i=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(i)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           molecule_list=molecule_list,&
                           natom=natom,atom_list=atom_list,&
                           nonfo=nonfo,onfo_list=onfo_list)
    molecule=>molecule_set(molecule_list(1))
    CALL get_molecule(molecule=molecule,first_atom=first,last_atom=last)
    DO j=1,nonfo
      atm_a = onfo_list(j)%a
      !atomic_kind => atom_list(atm_a)%atomic_kind
      atomic_kind => particle_set(atm_a+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_a, kind_number = ikind )
      atm_b = onfo_list(j)%b
      !atomic_kind => atom_list(atm_b)%atomic_kind
      atomic_kind => particle_set(atm_b+first-1)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=name_atm_b, kind_number = jkind )
      found = .FALSE.
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)
      ! loop over params from amber
!     IF(ASSOCIATED(chm_info%onfo_a)) THEN
!       DO k=1,SIZE(chm_info%onfo_a)
!         IF( ((TRIM(chm_info%onfo_a(k))==TRIM(name_atm_a)) .AND. &
!              (TRIM(chm_info%onfo_b(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(chm_info%onfo_c(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(chm_info%onfo_d(k))==TRIM(name_atm_d))) .OR. &
!             ((TRIM(chm_info%onfo_a(k))==TRIM(name_atm_d)) .AND. &
!              (TRIM(chm_info%onfo_b(k))==TRIM(name_atm_c)) .AND. &
!              (TRIM(chm_info%onfo_c(k))==TRIM(name_atm_b)) .AND. &
!              (TRIM(chm_info%onfo_d(k))==TRIM(name_atm_a))) ) THEN
!           onfo_list(j)%onfo_kind%k      = chm_info%onfo_k(k)
!           onfo_list(j)%onfo_kind%phi0 = chm_info%onfo_phi0(k)
!           IF(found .AND. PM1) &
!              WRITE(iw,*) "    Multiple Impropers declarations: overwriting"
!           found = .TRUE.
!           EXIT
!         END IF
!       END DO
!     END IF

!      IF(ASSOCIATED(chm_info%nonbond_a_14)) THEN
      onfo_list(j)%onfo_kind%name   = "CHARMM"
! Point to the correct potparm % ( ikind, jkind ) % spline_data
      CALL spline_data_retain (potparm_bond (ikind, jkind ) % pair_spline_data,error=error)
      CALL spline_data_release(onfo_list(j)%onfo_kind%spline_data, error=error)
      onfo_list(j)%onfo_kind%spline_data => &
      potparm_bond ( ikind, jkind ) % pair_spline_data
      
!      END IF

!      ! always have the input param last to overwrite all the other ones
!      IF(ASSOCIATED(inp_info%nonbonded14)) THEN
!        DO k=1,SIZE(inp_info%nonbonded14)
!          IF(PM2) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
!              " with ",TRIM(inp_info%nonbonded14(k)%at1),&
!              TRIM(inp_info%nonbonded14(k)%at2)
!          IF(((TRIM(name_atm_a)==TRIM(inp_info%nonbonded14(k)%at1))  .AND. &
!              (TRIM(name_atm_b)==TRIM(inp_info%nonbonded14(k)%at2))) .OR.  &
!             ((TRIM(name_atm_b)==TRIM(inp_info%nonbonded14(k)%at1))  .AND. &
!              (TRIM(name_atm_a)==TRIM(inp_info%nonbonded14(k)%at2))) ) THEN
!            potparm_bond(i,j) = inp_info%nonbonded14(k)
!            IF(found .AND. PM1) &
!               WRITE(iw,*) "    Multiple LJ-POT declarations: overwriting"
!            IF(PM2) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a),&
!                    " ",TRIM(name_atm_b)
!            found = .TRUE.
!            EXIT
!          END IF
!        END DO
!      END IF
!
!      IF(.NOT.found) THEN
      IF(PM2) WRITE(iw,*) "    ONE-FOUR ",j," : ",&
                          TRIM(name_atm_a)," ",&
                          TRIM(name_atm_b)
!      END IF
    END DO
    CALL set_molecule_kind(molecule_kind=molecule_kind,&
                           onfo_list=onfo_list)
  END DO
  CALL pot_parameters_release ( potparm_bond, error )


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 14. Set up atomic_kind_set()%elp_potentail%qeff (PART 2)
  !-----------------------------------------------------------------------------
  ALLOCATE(charge(SIZE(atomic_kind_set)),STAT=istat)
  IF (istat /= 0) CALL stop_memory('force_field_pack','charge')
  DO i=1,SIZE(atomic_kind_set)
    atomic_kind => atomic_kind_set(i)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
                         elp_potential=elp_potential,&
                         atom_list=my_atom_list,&
                         name=atmname) 
    CALL get_potential(potential=elp_potential,&
                       qeff=charge(i))

    found = .TRUE.
    only_qm = qmmm_ff_precond_only_qm(id1=atmname,is_link=is_link_atom)
    IF(charge(i)<=-100.0) found = .FALSE.

!   ! loop over params from amber
!   IF(ASSOCIATED(ff_type%amb_info%charge_atm)) THEN
!     DO j=1,SIZE(ff_type%amb_info%charge_atm)
!       IF(TRIM(ff_type%amb_info%charge_atm(j))==atmname) THEN
!         charge(i) = ff_type%amb_info%charge(j)
!         IF(found .AND. PM1) &
!           WRITE(iw,*) "    Multiple Charge declarations: overwriting"
!         found = .TRUE.
!       END IF
!       IF(found) EXIT
!     END DO
!   END IF
  
    ! loop over params from charmm
    IF(ASSOCIATED(ff_type%chm_info%charge_atm)) THEN
      DO j=1,SIZE(ff_type%chm_info%charge_atm)
        IF(TRIM(ff_type%chm_info%charge_atm(j))==atmname) THEN
          charge(i) = ff_type%chm_info%charge(j)
          IF(found .AND. PM1) &
            WRITE(iw,*) "    Multiple Charge declarations: overwriting"
          found = .TRUE.
        END IF
        !IF(found) EXIT
      END DO
    END IF
  
    ! always have the input param last to overwrite all the other ones
    IF(ASSOCIATED(ff_type%inp_info%charge_atm)) THEN
      DO j=1,SIZE(ff_type%inp_info%charge_atm)
        IF(TRIM(ff_type%inp_info%charge_atm(j))==atmname) THEN
          charge(i) = ff_type%inp_info%charge(j)
          IF(found .AND. PM1) &
            WRITE(iw,*) "    Multiple Charge declarations: overwriting"
          found = .TRUE.
        END IF
        !IF(found) EXIT
      END DO
    END IF

    IF(.NOT.found) THEN
      CALL stop_program("FF_coordinate_pack","missing qeff parm "//atmname)
    END IF
    !
    ! QM/MM modifications
    !
    IF (only_qm.AND.my_qmmm) THEN
       IF (qmmm_env%qmmm_coupl_type /= do_qmmm_none) THEN 
          scale_factor = 0.0_dp
          IF (is_link_atom) THEN
             !
             ! Find the scaling factor...
             !
             DO ilink = 1, SIZE(qmmm_env%mm_link_atoms)
                IF (ANY(my_atom_list == qmmm_env%mm_link_atoms(ilink))) EXIT
             END DO
             CPPostcondition(ilink <= SIZE(qmmm_env%mm_link_atoms),cp_failure_level,routineP,error,failure)
             scale_factor = qmmm_env%fist_scale_charge_link(ilink)
          END IF
          charge(i) = charge(i) * scale_factor
       END IF
    END IF
    
    CALL set_potential(potential=elp_potential,&
                       qeff=charge(i))
    CALL set_atomic_kind(atomic_kind=atomic_kind,&
                         elp_potential=elp_potential)
  END DO
  DEALLOCATE (charge,STAT=istat)
  IF (istat /= 0) THEN
    CALL stop_memory('force_field_pack','release work storage')
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 15. Set spline_nonbond_env
  !-----------------------------------------------------------------------------
  CALL pot_parameters_create ( potparm_nonbond,  &
                               SIZE ( atomic_kind_set ), &
                               error )
  DO i=1,SIZE(atomic_kind_set)
    atomic_kind => atomic_kind_set(i)
    CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_a_local)
    DO j=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(j)
      CALL get_atomic_kind(atomic_kind=atomic_kind,name=name_atm_b_local)
      found = .FALSE.
      name_atm_a = name_atm_a_local
      name_atm_b = name_atm_b_local
      only_qm = qmmm_ff_precond_only_qm(id1=name_atm_a, id2=name_atm_b)

      ! loop over params from amber
      ! loop over params from amber
      ! loop over params from amber
      ! loop over params from amber

      ! loop over params from charmm
      IF(ASSOCIATED(chm_info%nonbond_a)) THEN
        ii = 0 
        jj = 0
        DO k=1,SIZE(chm_info%nonbond_a)
          IF(TRIM(name_atm_a)==TRIM(chm_info%nonbond_a(k))) THEN
            ii = k
          END IF
        END DO
        DO k=1,SIZE(chm_info%nonbond_a)
          IF(TRIM(name_atm_b)==TRIM(chm_info%nonbond_a(k))) THEN
            jj = k
          END IF
        END DO

        IF(ii/=0 .AND. jj/=0) THEN
          rmin = chm_info%nonbond_rmin2(ii)+chm_info%nonbond_rmin2(jj)
          epsilon = sqrt(chm_info%nonbond_eps(ii)*&
                         chm_info%nonbond_eps(jj))
          potparm_nonbond(i,j)%type       = 'LENNARD-JONES'
          potparm_nonbond(i,j)%at1        = name_atm_a
          potparm_nonbond(i,j)%at2        = name_atm_b
          potparm_nonbond(i,j)%lj%epsilon = epsilon
          potparm_nonbond(i,j)%lj%sigma6  = 0.5_dp*rmin**6
          potparm_nonbond(i,j)%lj%sigma12 = 0.25_dp*rmin**12
          IF(ff_type%rcut_nb>0.0_dp) THEN
            potparm_nonbond(i,j)%rcutsq     = ff_type%rcut_nb*ff_type%rcut_nb
          ELSE
            !Corresponds to an cutoff of 10A
            potparm_nonbond(i,j)%rcutsq     = default_rcutsq
          END IF

          IF(found .AND. PM1) &
             WRITE(iw,*) "    Multiple LJ-POT declarations: overwriting"
          found = .TRUE.
        END IF
      END IF

      ! always have the input param last to overwrite all the other ones
      IF(ASSOCIATED(inp_info%nonbonded)) THEN
        DO k=1,SIZE(inp_info%nonbonded)
          IF(PM2) WRITE(iw,*) "    TESTING ",TRIM(name_atm_a),TRIM(name_atm_b),&
              " with ",TRIM(inp_info%nonbonded(k)%at1),&
              TRIM(inp_info%nonbonded(k)%at2)
          IF(((TRIM(name_atm_a)==TRIM(inp_info%nonbonded(k)%at1))  .AND. &
              (TRIM(name_atm_b)==TRIM(inp_info%nonbonded(k)%at2))) .OR.  &
             ((TRIM(name_atm_b)==TRIM(inp_info%nonbonded(k)%at1))  .AND. &
              (TRIM(name_atm_a)==TRIM(inp_info%nonbonded(k)%at2))) ) THEN
            potparm_nonbond(i,j) = inp_info%nonbonded(k)
            IF(found .AND. PM1) &
               WRITE(iw,*) "    Multiple LJ-POT declarations: overwriting"
            IF(PM2) WRITE(iw,*) "    FOUND ",TRIM(name_atm_a),&
                    " ",TRIM(name_atm_b)
            found = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      IF(.NOT.found) THEN
        IF(PM1) WRITE(iw,*) "    spline_nonbond_env ",j," : ",&
                            TRIM(name_atm_a)," ",&
                            TRIM(name_atm_b)
        CALL stop_program("FF_pack","missing pair pot")
      END IF

      IF (only_qm) THEN
         CALL clean_pot_parameter( potparm_nonbond(i,j) )
         potparm_nonbond(i,j)%type       = 'LENNARD-JONES'         
         potparm_nonbond(i,j)%rcutsq     = default_rcutsq
      END IF
    END DO
  END DO
  nkinds = SIZE ( potparm_nonbond, 1 )
!---------------------------------------
! create the pair potential spline environment
!---------------------------------------
  CALL spline_env_create(spline_nonbond_env, nkinds, error)
  CALL spline_nonbond_control(spline_nonbond_env, potparm_nonbond, &
                              atomic_kind_set,ewald_env)
!---------------------------------------
! create the manybody spline environment
!---------------------------------------
  CALL spline_env_create(spline_manybody_env, nkinds, error)
  CALL spline_manybody_control(spline_manybody_env, potparm_nonbond)
                              
  DO ikind = 1, SIZE ( potparm_nonbond, 1 )
    DO jkind = ikind, SIZE ( potparm_nonbond, 2)
!---------------------------------------
! 1st: The pair-potential contribution
!---------------------------------------
      n = spline_nonbond_env % spltab ( ikind, jkind )
      spline_data => spline_nonbond_env % spline_data_p ( n ) % spline_data
! do i,j
      CALL spline_data_retain ( spline_data, error )
      CALL spline_data_release (  &
         potparm_nonbond ( ikind, jkind ) % pair_spline_data, error )
      potparm_nonbond ( ikind, jkind ) % pair_spline_data => spline_data
! do j,i
      CALL spline_data_retain ( spline_data, error )
      CALL spline_data_release (  &
         potparm_nonbond ( jkind, ikind ) % pair_spline_data, error )
      potparm_nonbond ( jkind, ikind ) % pair_spline_data => spline_data
!---------------------------------------
! 2nd: The many-body contribution
!---------------------------------------
      n = spline_manybody_env % spltab ( ikind, jkind )
      spline_data => spline_manybody_env % spline_data_p ( n ) % spline_data
! do i,j
      IF ( ASSOCIATED ( spline_data ) ) THEN
        CALL spline_data_retain ( spline_data, error )
        CALL spline_data_release (  &
           potparm_nonbond ( ikind, jkind ) % rho_spline_data, error )
        potparm_nonbond ( ikind, jkind ) % rho_spline_data => spline_data
! do j,i
        CALL spline_data_retain ( spline_data, error )
        CALL spline_data_release (  &
           potparm_nonbond ( jkind, ikind ) % rho_spline_data, error )
        potparm_nonbond ( jkind, ikind ) % rho_spline_data => spline_data
      ENDIF 
    END DO
  END DO
  CALL fist_nonbond_env_create (  fist_nonbond_env, potparm_nonbond, error )

! Release all spline_environments
  CALL spline_env_release (spline_manybody_env,error)
  CALL spline_env_release (spline_nonbond_env,error)
  CALL spline_env_release (spline_bond_env,error)


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  force_field_pack"
  CALL write_checkpoint_information("leaving "//routineN,globenv)

END SUBROUTINE force_field_pack



!******************************************************************************
!!****** force_field_util/force_field_qeff_output [1.0] *
!!
!!   NAME
!!     force_field_qeff_output
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     Compute the total qeff charges for each molecule kind and total system 
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************

SUBROUTINE force_field_qeff_output (particle_set,atomic_kind_set,&
                             molecule_kind_set,molecule_set,globenv)

  TYPE(particle_type), DIMENSION(:), POINTER             :: particle_set
  TYPE(atomic_kind_type), DIMENSION(:), POINTER          :: atomic_kind_set
  TYPE(molecule_kind_type), DIMENSION(:), POINTER        :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), POINTER             :: molecule_set
  TYPE(global_environment_type), INTENT(IN), OPTIONAL    :: globenv

  CHARACTER(len=*), PARAMETER :: routineN = 'force_field_qeff_output', &
    routineP = moduleN//':'//routineN

  INTEGER                                  :: istat, iw, i, j, k
  LOGICAL                                  :: PM1, PM2, PM3, PM4

  CHARACTER(LEN=default_string_length)     :: atmname,molname
  TYPE(atomic_kind_type), POINTER          :: atomic_kind
  TYPE(molecule_kind_type), POINTER        :: molecule_kind

  INTEGER                                  :: iatom,imol
  REAL(KIND=dp)                                :: qeff,qeff_mol,qeff_sum

  TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering force_field_qeff_output"

  qeff     = 0.0_dp
  qeff_mol = 0.0_dp
  qeff_sum = 0.0_dp

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Sum of qeff for each molecule_kind
  !-----------------------------------------------------------------------------
  DO imol=1,SIZE(molecule_kind_set)
    qeff_mol=0.0_dp
    molecule_kind => molecule_kind_set(imol)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           name=molname,atom_list=atom_list)
    DO iatom=1,SIZE(atom_list)
      atomic_kind => atom_list(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           name=atmname,qeff=qeff)
      IF(PM2) WRITE(iw,*) "      atom ",iatom," ",TRIM(atmname)," charge = ",qeff
      qeff_mol = qeff_mol + qeff
    END DO
    IF(PM1) WRITE(iw,*) "    Mol Kind ",TRIM(molname)," charge = ",qeff_mol
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Sum of qeff for particle_set
  !-----------------------------------------------------------------------------
  DO iatom=1,SIZE(particle_set)
    atomic_kind => particle_set(iatom)%atomic_kind
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
                         name=atmname,qeff=qeff)
    IF(PM2) WRITE(iw,*) "      atom ",iatom," ",TRIM(atmname),&
                        " charge = ",qeff
    qeff_sum = qeff_sum + qeff
  END DO
  IF(PM1) WRITE(iw,*) "    Total system charge = ",qeff_sum

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  force_field_qeff_output"
  CALL write_checkpoint_information("leaving "//routineN,globenv)

END SUBROUTINE force_field_qeff_output


!******************************************************************************
!!****** force_field/clean_intra_forc_kind [1.0] *
!!
!!   NAME
!!     clean_intra_forc_kind
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     blah blah blah blah blah
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************
SUBROUTINE clean_intra_force_kind (molecule_kind_set,globenv)

  TYPE(molecule_kind_type), DIMENSION(:), POINTER        :: molecule_kind_set
  TYPE(global_environment_type), INTENT(IN), OPTIONAL    :: globenv

  CHARACTER(len=*), PARAMETER :: routineN = 'clean_intra_force_kind', &
    routineP = moduleN//':'//routineN

  INTEGER                                  :: istat, iw, i, j, k
  LOGICAL                                  :: PM1, PM2, PM3, PM4

  TYPE(bond_kind_type), DIMENSION(:), POINTER       :: bond_kind_set
  TYPE(bond_kind_type), DIMENSION(:), POINTER       :: new_bond_kind_set
  TYPE(bond_type), DIMENSION(:), POINTER            :: bond_list
  TYPE(bond_type), DIMENSION(:), POINTER            :: new_bond_list
  TYPE(bend_kind_type), DIMENSION(:), POINTER       :: bend_kind_set
  TYPE(bend_kind_type), DIMENSION(:), POINTER       :: new_bend_kind_set
  TYPE(bend_type), DIMENSION(:), POINTER            :: bend_list
  TYPE(bend_type), DIMENSION(:), POINTER            :: new_bend_list
  TYPE(ub_kind_type), DIMENSION(:), POINTER         :: ub_kind_set
  TYPE(ub_kind_type), DIMENSION(:), POINTER         :: new_ub_kind_set
  TYPE(ub_type), DIMENSION(:), POINTER              :: ub_list
  TYPE(ub_type), DIMENSION(:), POINTER              :: new_ub_list
  TYPE(torsion_kind_type), DIMENSION(:), POINTER    :: torsion_kind_set
  TYPE(torsion_kind_type), DIMENSION(:), POINTER    :: new_torsion_kind_set
  TYPE(torsion_type), DIMENSION(:), POINTER         :: torsion_list
  TYPE(torsion_type), DIMENSION(:), POINTER         :: new_torsion_list
  TYPE(impr_kind_type), DIMENSION(:), POINTER       :: impr_kind_set
  TYPE(impr_kind_type), DIMENSION(:), POINTER       :: new_impr_kind_set
  TYPE(impr_type), DIMENSION(:), POINTER            :: impr_list
  TYPE(impr_type), DIMENSION(:), POINTER            :: new_impr_list
  INTEGER                                  :: ikind,newkind,ibond,ibend,iub
  INTEGER                                  :: itorsion,iimpr
  TYPE(molecule_kind_type), POINTER        :: molecule_kind
  INTEGER, POINTER                         :: bad1(:),bad2(:)
  INTEGER                                  :: counter
  INTEGER                                  :: nbond,nbend,nub,ntorsion,nimpr
  INTEGER                                  :: idist,ndist
  INTEGER                                  :: ig3x3,ng3x3
  INTEGER                                  :: ig4x6,ng4x6
  INTEGER                                  :: atm_a,atm_b,atm_c,atm_d
  INTEGER                                  :: atm2_a,atm2_b,atm2_c,atm2_d
  INTEGER                                  :: central
  TYPE(dist_constraint_type), DIMENSION(:), POINTER :: dist_list
  TYPE(g3x3_constraint_type), DIMENSION(:), POINTER :: g3x3_list
  TYPE(g4x6_constraint_type), DIMENSION(:), POINTER :: g4x6_list
  LOGICAL                                  :: unsetme,valid_kind
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  PM1 = .FALSE.
  PM2 = .FALSE.
  PM3 = .FALSE.
  PM4 = .FALSE.
  IF(PRESENT(globenv)) THEN
    CALL write_checkpoint_information("entering "//routineN,globenv)
    iw = globenv%scr
    IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
    IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
  END IF

  IF(PM1) WRITE(iw,*) "  Entering clean_intra_force_kind"
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Lets Tag the unwanted bonds due to the use of distance constraint
  !-----------------------------------------------------------------------------
  DO ikind=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(ikind)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           ndist=ndist,&
                           dist_list=dist_list,&
                           nbond=nbond,&
                           bond_kind_set=bond_kind_set,&
                           bond_list=bond_list)
    DO idist=1,ndist
      atm_a = dist_list(idist)%a
      atm_b = dist_list(idist)%b
      DO ibond=1,nbond
        unsetme = .FALSE.
        atm2_a = bond_list(ibond)%a
        atm2_b = bond_list(ibond)%b
        IF(atm2_a==atm_a .AND. atm2_b==atm_b) unsetme=.TRUE.
        IF(atm2_a==atm_b .AND. atm2_b==atm_a) unsetme=.TRUE.
        IF(unsetme) bond_list(ibond)%name = "UNSET"
      END DO
    END DO
  END DO


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Lets Tag the unwanted bends due to the use of distance constraint
  !-----------------------------------------------------------------------------
  DO ikind=1,SIZE(molecule_kind_set)
    molecule_kind => molecule_kind_set(ikind)
    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                           ndist=ndist,&
                           dist_list=dist_list,&
                           nbend=nbend,&
                           bend_kind_set=bend_kind_set,&
                           bend_list=bend_list)
    DO ibend=1,nbend
      counter=0
      unsetme = .FALSE.
      atm_a = bend_list(ibend)%a
      atm_b = bend_list(ibend)%b
      atm_c = bend_list(ibend)%c
      central = 0
      DO idist=1,ndist
        atm2_a = dist_list(idist)%a
        atm2_b = dist_list(idist)%b
! Check to see if we have a central bend atom in dist list
        IF ( central == 0 ) THEN
          IF ( atm2_a == atm_b ) central = atm2_a
          IF ( atm2_b == atm_b ) central = atm2_b
        ENDIF
! If so, check the other two atoms for a match
        IF ( atm2_a == central ) THEN 
          IF(atm2_b==atm_a.OR.atm2_b==atm_c) counter=counter+1
! If so, check the other two atoms for a match
        ELSEIF ( atm2_b == central ) THEN
          IF(atm2_a==atm_a.OR.atm2_a==atm_c) counter=counter+1
! If not, check the other two atoms for a match
        ELSE
          IF(atm2_a==atm_a.OR.atm2_a==atm_c) counter=counter+1
          IF(atm2_b==atm_a.OR.atm2_b==atm_c) counter=counter+1
        ENDIF
      END DO
! must have central atom to be involved in bend
      IF ( central /= 0 ) THEN
        IF(counter>=4) unsetme=.TRUE.
        IF(unsetme) bend_list(ibend)%name = "UNSET"
      ENDIF
    END DO
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
                           bond_kind_set=bond_kind_set,&
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
        IF(unsetme) bond_list(ibond)%name = "UNSET"
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
                           bend_kind_set=bend_kind_set,&
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
        IF(unsetme) bend_list(ibend)%name = "UNSET"
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
                           bond_kind_set=bond_kind_set,&
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
        IF(unsetme) bond_list(ibond)%name = "UNSET"
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
                           bend_kind_set=bend_kind_set,&
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
        IF(unsetme) bend_list(ibend)%name = "UNSET"
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
    IF(PM1) WRITE(iw,*) "    Mol(",ikind,") Old BOND Count: ",&
                        SIZE(bond_list),SIZE(bond_kind_set)
    IF (PM2) WRITE(iw,'(2I6)')(bond_list(ibond)%a,bond_list(ibond)%b,ibond=1,size(bond_list))
    NULLIFY(bad1,bad2)
    ALLOCATE(bad1(SIZE(bond_kind_set)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
    bad1(:)=0
    DO ibond=1,SIZE(bond_kind_set)
      unsetme=.FALSE.
      IF(bond_kind_set(ibond)%name=="UNSET") unsetme = .TRUE.
      valid_kind=.FALSE.
      DO i=1,SIZE(bond_list)
        IF(bond_list(i)%name/="UNSET".AND.&
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
      CALL allocate_bond_kind_set(new_bond_kind_set,counter)
      counter=0
      DO ibond=1,SIZE(bond_kind_set)
        IF(bad1(ibond)==0) THEN
          counter=counter+1
          new_bond_kind_set(counter)= bond_kind_set(ibond)     
        END IF
      END DO
      counter=0
      ALLOCATE(bad2(SIZE(bond_list)),STAT=istat)
      IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad2')
      bad2(:)=0
      DO ibond=1,SIZE(bond_list)
        unsetme = .FALSE.
        IF(bond_list(ibond)%bond_kind%name=="UNSET") unsetme = .TRUE.
        IF(bond_list(ibond)%name=="UNSET")           unsetme = .TRUE.
        IF(unsetme) bad2(ibond) = 1
      END DO
      IF(SUM(bad2)/=0) THEN
        counter = SIZE(bond_list)-SUM(bad2)
        ALLOCATE(new_bond_list(counter),STAT=istat)
        IF(istat/=0) CALL stop_memory('clean_intra_force_kind','new_bond_list')
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
      DEALLOCATE(bad2,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','bad2')
      CALL deallocate_bond_kind_set(bond_kind_set)
      DEALLOCATE(bond_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','bond_list')
      IF(PM1) WRITE(iw,*) "    Mol(",ikind,") New BOND Count: ",&
                          SIZE(new_bond_list),SIZE(new_bond_kind_set)
      IF (PM2) WRITE(iw,'(2I6)')(new_bond_list(ibond)%a,new_bond_list(ibond)%b,ibond=1,SIZE(new_bond_list))
    END IF
    DEALLOCATE(bad1,STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
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
    IF(PM1) WRITE(iw,*) "    Mol(",ikind,") Old BEND Count: ",&
                        SIZE(bend_list),SIZE(bend_kind_set)
    IF (PM2) WRITE(iw,'(3I6)')(bend_list(ibend)%a,bend_list(ibend)%b,&
         bend_list(ibend)%c,ibend=1,SIZE(bend_list))
    NULLIFY(bad1,bad2)
    ALLOCATE(bad1(SIZE(bend_kind_set)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
    bad1(:)=0
    DO ibend=1,SIZE(bend_kind_set)
      unsetme=.FALSE.
      IF(bend_kind_set(ibend)%name=="UNSET") unsetme = .TRUE.
      valid_kind=.FALSE.
      DO i=1,SIZE(bend_list)
        IF(bend_list(i)%name/="UNSET".AND.&
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
      CALL allocate_bend_kind_set(new_bend_kind_set,counter)
      counter=0
      DO ibend=1,SIZE(bend_kind_set)
        IF(bad1(ibend)==0) THEN
          counter=counter+1
          new_bend_kind_set(counter)= bend_kind_set(ibend)     
        END IF
      END DO
      counter=0
      ALLOCATE(bad2(SIZE(bend_list)),STAT=istat)
      IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad2')
      bad2(:)=0
      DO ibend=1,SIZE(bend_list)
        unsetme = .FALSE.
        IF(bend_list(ibend)%bend_kind%name=="UNSET") unsetme = .TRUE.
        IF(bend_list(ibend)%name=="UNSET")           unsetme = .TRUE.
        IF(unsetme) bad2(ibend) = 1
      END DO
      IF(SUM(bad2)/=0) THEN
        counter = SIZE(bend_list)-SUM(bad2)
        ALLOCATE(new_bend_list(counter),STAT=istat)
        IF(istat/=0) CALL stop_memory('clean_intra_force_kind','new_bend_list')
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
      DEALLOCATE(bad2,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','bad2')
      CALL deallocate_bend_kind_set(bend_kind_set)
      DEALLOCATE(bend_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','bend_list')
      IF(PM1) WRITE(iw,*) "    Mol(",ikind,") New BEND Count: ",&
                          SIZE(new_bend_list),SIZE(new_bend_kind_set)
    IF (PM2) WRITE(iw,'(3I6)')(new_bend_list(ibend)%a,new_bend_list(ibend)%b,&
         new_bend_list(ibend)%c,ibend=1,SIZE(new_bend_list))
    END IF
    DEALLOCATE(bad1,STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
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
    IF(PM1) WRITE(iw,*) "    Mol(",ikind,") Old UB Count: ",&
                        SIZE(ub_list),SIZE(ub_kind_set)
    IF (PM2) WRITE(iw,'(3I6)')(ub_list(iub)%a,ub_list(iub)%b,&
         ub_list(iub)%c,iub=1,SIZE(ub_list))
    NULLIFY(bad1,bad2)
    ALLOCATE(bad1(SIZE(ub_kind_set)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
    bad1(:)=0
    DO iub=1,SIZE(ub_kind_set)
      unsetme=.FALSE.
      IF(ub_kind_set(iub)%name=="UNSET") unsetme = .TRUE.
      valid_kind=.FALSE.
      DO i=1,SIZE(ub_list)
        IF(ub_list(i)%name/="UNSET".AND.&
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
      CALL allocate_ub_kind_set(new_ub_kind_set,counter)
      counter=0
      DO iub=1,SIZE(ub_kind_set)
        IF(bad1(iub)==0) THEN
          counter=counter+1
          new_ub_kind_set(counter)= ub_kind_set(iub)     
        END IF
      END DO
      counter=0
      ALLOCATE(bad2(SIZE(ub_list)),STAT=istat)
      IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad2')
      bad2(:)=0
      DO iub=1,SIZE(ub_list)
        unsetme = .FALSE.
        IF(ub_list(iub)%ub_kind%name=="UNSET") unsetme = .TRUE.
        IF(ub_list(iub)%name=="UNSET")           unsetme = .TRUE.
        IF(unsetme) bad2(iub) = 1
      END DO
      IF(SUM(bad2)/=0) THEN
        counter = SIZE(ub_list)-SUM(bad2)
        ALLOCATE(new_ub_list(counter),STAT=istat)
        IF(istat/=0) &
          CALL stop_memory('clean_intra_force_kind','new_ub_list')
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
      DEALLOCATE(bad2,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','bad2')
      CALL ub_kind_dealloc_ref(ub_kind_set)
      DEALLOCATE(ub_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','ub_list')
      IF(PM1) WRITE(iw,*) "    Mol(",ikind,") New UB Count: ",&
                          SIZE(new_ub_list),SIZE(new_ub_kind_set)
    IF (PM2) WRITE(iw,'(3I6)')(new_ub_list(iub)%a,new_ub_list(iub)%b,&
         new_ub_list(iub)%c,iub=1,SIZE(new_ub_list))
    END IF
    DEALLOCATE(bad1,STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
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
    IF(PM1) WRITE(iw,*) "    Mol(",ikind,") Old TORSION Count: ",&
                        SIZE(torsion_list),SIZE(torsion_kind_set)
    IF (PM2) WRITE(iw,'(4I6)')(torsion_list(itorsion)%a,torsion_list(itorsion)%b,&
         torsion_list(itorsion)%c,torsion_list(itorsion)%d,itorsion=1,SIZE(torsion_list))
    NULLIFY(bad1,bad2)
    ALLOCATE(bad1(SIZE(torsion_kind_set)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
    bad1(:)=0
    DO itorsion=1,SIZE(torsion_kind_set)
      unsetme=.FALSE.
      IF(torsion_kind_set(itorsion)%name=="UNSET") unsetme = .TRUE.
      valid_kind=.FALSE.
      DO i=1,SIZE(torsion_list)
        IF(torsion_list(i)%name/="UNSET".AND.&
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
      CALL allocate_torsion_kind_set(new_torsion_kind_set,counter)
      counter=0
      DO itorsion=1,SIZE(torsion_kind_set)
        IF(bad1(itorsion)==0) THEN
          counter=counter+1
          new_torsion_kind_set(counter)= torsion_kind_set(itorsion)     
          i = SIZE(torsion_kind_set(itorsion)%m)
          j = SIZE(torsion_kind_set(itorsion)%k)
          k = SIZE(torsion_kind_set(itorsion)%cosphi0)
          ALLOCATE(new_torsion_kind_set(counter)%m(i),STAT=istat)
          IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','m')
          ALLOCATE(new_torsion_kind_set(counter)%k(i),STAT=istat)
          IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','k')
          ALLOCATE(new_torsion_kind_set(counter)%cosphi0(i),STAT=istat)
          IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','cosphi0')
          new_torsion_kind_set(counter)%m= torsion_kind_set(itorsion)%m
          new_torsion_kind_set(counter)%k= torsion_kind_set(itorsion)%k
          new_torsion_kind_set(counter)%cosphi0= torsion_kind_set(itorsion)%cosphi0     
        END IF
      END DO
      counter=0
      ALLOCATE(bad2(SIZE(torsion_list)),STAT=istat)
      IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad2')
      bad2(:)=0
      DO itorsion=1,SIZE(torsion_list)
        unsetme = .FALSE.
        IF(torsion_list(itorsion)%torsion_kind%name=="UNSET") unsetme = .TRUE.
        IF(torsion_list(itorsion)%name=="UNSET")           unsetme = .TRUE.
        IF(unsetme) bad2(itorsion) = 1
      END DO
      IF(SUM(bad2)/=0) THEN
        counter = SIZE(torsion_list)-SUM(bad2)
        ALLOCATE(new_torsion_list(counter),STAT=istat)
        IF(istat/=0) &
          CALL stop_memory('clean_intra_force_kind','new_torsion_list')
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
      DEALLOCATE(bad2,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','bad2')
      DO itorsion=1,SIZE(torsion_kind_set)
        CALL torsion_kind_dealloc_ref(torsion_kind_set(itorsion))
      END DO
      DEALLOCATE(torsion_kind_set,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','torsion_kind_set')
      DEALLOCATE(torsion_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','torsion_list')
      IF(PM1) WRITE(iw,*) "    Mol(",ikind,") New TORSION Count: ",&
                          SIZE(new_torsion_list),SIZE(new_torsion_kind_set)
    IF (PM2) WRITE(iw,'(4I6)')(new_torsion_list(itorsion)%a,new_torsion_list(itorsion)%b,&
         new_torsion_list(itorsion)%c,new_torsion_list(itorsion)%d,itorsion=1,&
         SIZE(new_torsion_list))
    END IF
    DEALLOCATE(bad1,STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
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
    IF(PM1) WRITE(iw,*) "    Mol(",ikind,") Old IMPROPER Count: ",&
                        SIZE(impr_list),SIZE(impr_kind_set)
    NULLIFY(bad1,bad2)
    ALLOCATE(bad1(SIZE(impr_kind_set)),STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
    bad1(:)=0
    DO iimpr=1,SIZE(impr_kind_set)
      unsetme=.FALSE.
      IF(impr_kind_set(iimpr)%name=="UNSET") unsetme = .TRUE.
      valid_kind=.FALSE.
      DO i=1,SIZE(impr_list)
        IF(impr_list(i)%name/="UNSET".AND.&
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
      CALL allocate_impr_kind_set(new_impr_kind_set,counter)
      counter=0
      DO iimpr=1,SIZE(impr_kind_set)
        IF(bad1(iimpr)==0) THEN
          counter=counter+1
          new_impr_kind_set(counter)= impr_kind_set(iimpr)     
        END IF
      END DO
      counter=0
      ALLOCATE(bad2(SIZE(impr_list)),STAT=istat)
      IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad2')
      bad2(:)=0
      DO iimpr=1,SIZE(impr_list)
        unsetme = .FALSE.
        IF(impr_list(iimpr)%impr_kind%name=="UNSET") unsetme = .TRUE.
        IF(impr_list(iimpr)%name=="UNSET")           unsetme = .TRUE.
        IF(unsetme) bad2(iimpr) = 1
      END DO
      IF(SUM(bad2)/=0) THEN
        counter = SIZE(impr_list)-SUM(bad2)
        ALLOCATE(new_impr_list(counter),STAT=istat)
        IF(istat/=0) CALL stop_memory('clean_intra_force_kind','new_impr_list')
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
      DEALLOCATE(bad2,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','bad2')
      DO iimpr=1,SIZE(impr_kind_set)
        CALL impr_kind_dealloc_ref(impr_kind_set(iimpr))
      END DO
      DEALLOCATE(impr_kind_set,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','impr_kind_set')
      DEALLOCATE(impr_list,STAT=istat)
      IF(istat/=0) CALL stop_memory('clean_intra_force_kind','impr_list')
      IF(PM1) WRITE(iw,*) "    Mol(",ikind,") New IMPROPER Count: ",&
                          SIZE(new_impr_list),SIZE(new_impr_kind_set)
    END IF
    DEALLOCATE(bad1,STAT=istat)
    IF (istat /= 0) CALL stop_memory('clean_intra_force_kind','bad1')
   END IF
  END DO
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  IF(PM1) WRITE(iw,*) "  Exiting  clean_intra_force_kind"
  CALL write_checkpoint_information("leaving "//routineN,globenv)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE clean_intra_force_kind

!******************************************************************************
!!****** force_field/read_eam_data [1.0] *
!!
!!   NAME
!!     read_eam_data
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     reads EAM potential from library
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*****************************************************************************
SUBROUTINE read_eam_data ( eam, globenv )
  TYPE ( eam_pot_type ), INTENT ( INOUT )                    :: eam
  TYPE ( global_environment_type ), INTENT(IN), OPTIONAL     :: globenv
!locals
  INTEGER :: unit_number, scr, i, ntyp, isos
  REAL ( dp ) :: rmass, rlatt
  CHARACTER ( len = default_string_length ) :: title3
  REAL ( dp ), POINTER, DIMENSION ( : ) :: rval, rhoval, phi, phip, frho, frhop, &
                                           rho, rhop
  REAL ( dp ) :: drar, drhoar, acutal

  NULLIFY ( rval, rhoval, phi, phip, frho, frhop, rho, rhop )
  IF ( globenv % ionode ) THEN
    ALLOCATE ( rho ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate rho' )
    rho=0.0_dp
    ALLOCATE ( rhop ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate rhop' )
    rhop=0.0_dp
    ALLOCATE ( rval ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate rval' )
    rval=0.0_dp
    ALLOCATE ( rhoval ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate rhoval' )
    rhoval=0.0_dp
    ALLOCATE ( phi ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate phi' )
    phi=0.0_dp
    ALLOCATE ( phip ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate phip' )
    phip=0.0_dp
    ALLOCATE ( frho ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate frho' )
    frho=0.0_dp
    ALLOCATE ( frhop ( 5000 ), STAT = isos )
    IF (isos/=0) CALL stop_program('read_eam_data', &
    'failed to allocate frhop' )
    frhop=0.0_dp

    IF ( PRESENT ( globenv ) ) scr = globenv % scr
    CALL open_file(file_name=TRIM(eam % eam_file_name),&
                 file_form="FORMATTED",&
                 file_position="REWIND",&
                 file_status="OLD",&
                 unit_number=unit_number)
    READ ( unit_number, '(A72)' ) title3
    WRITE ( scr, * ) title3

    READ ( unit_number, * ) ntyp,rmass,rlatt

!  WRITE ( scr, '(1x,72('*'),//1x,'input data from the EAM functions',  &
!         //1x,'title: ',a72,                                           &
!         /1x,'-> Z =',i5,8x,'A =',f10.3,                               &
!         /1x,'-> a =',f10.4)') title3,ntyp,rmass,rlatt

    READ ( unit_number, * ) drar, drhoar, acutal
        CALL convert_to_cp2k_units ( "ANGSTROM", length = drar )
        CALL convert_to_cp2k_units ( "ANGSTROM", length = acutal )

!  WRITE(scr,'('  dr, drho, rcut = ',3e12.4)') eam % drar, eam % drhoar, eam % acutal
!
    DO i = 1,5000
      READ ( unit_number, * ) rho(i), rhop(i)
      CALL convert_to_cp2k_units ( "ANGSTROM", length = rhop ( i ), l_power = -1 )
      rval(i) = REAL ( i-1 ) * drar
      rhoval(i) = REAL ( i-1 )* drhoar
    END DO
!
    DO i = 1,5000
       READ ( unit_number, * ) phi(i), phip(i)
       CALL convert_to_cp2k_units ( "EV", energy = phi ( i ) )
       CALL convert_to_cp2k_units ( "EV", energy = phip ( i ) )
       CALL convert_to_cp2k_units ( "ANGSTROM", length = phip ( i ), l_power = -1 )
    END DO
!
    DO i = 1,5000
       READ ( unit_number, * ) frho(i), frhop(i)
       CALL convert_to_cp2k_units ( "EV", energy = frho ( i ) )
       CALL convert_to_cp2k_units ( "EV", energy = frhop ( i ) )
    END DO
!
    WRITE ( scr, * ) 'Finished EAM data'
    CLOSE ( unit_number )
  ENDIF
  CALL mp_bcast ( drar, globenv % source, globenv % group )
  CALL mp_bcast ( drhoar, globenv % source, globenv % group )
  CALL mp_bcast ( acutal, globenv % source, globenv % group )

  eam % drar = drar
  eam % drhoar = drhoar
  eam % acutal = acutal
  
  CALL reallocate ( rho, 1, 5000 )
  CALL mp_bcast ( rho, globenv % source, globenv % group )

  CALL reallocate ( rhop, 1, 5000 )
  CALL mp_bcast ( rhop, globenv % source, globenv % group )

  CALL reallocate ( rval, 1, 5000 )
  CALL mp_bcast ( rval, globenv % source, globenv % group )

  CALL reallocate ( rhoval, 1, 5000 )
  CALL mp_bcast ( rhoval, globenv % source, globenv % group )

  CALL reallocate ( phi, 1, 5000 )
  CALL mp_bcast ( phi, globenv % source, globenv % group )

  CALL reallocate ( phip, 1, 5000 )
  CALL mp_bcast ( phip, globenv % source, globenv % group )

  CALL reallocate ( frho, 1, 5000 )
  CALL mp_bcast ( frho, globenv % source, globenv % group )

  CALL reallocate ( frhop, 1, 5000 )
  CALL mp_bcast ( frhop, globenv % source, globenv % group )
  
  DO i = 1, 5000
    eam % rho ( i ) = rho ( i )
    eam % rhop ( i ) = rhop ( i )
    eam % rval ( i ) = rval ( i )
    eam % rhoval ( i ) = rhoval ( i )
    eam % phi ( i ) = phi ( i )
    eam % phip ( i ) = phip ( i )
    eam % frho ( i ) = frho ( i )
    eam % frhop ( i ) = frhop ( i )
  END DO
  DEALLOCATE ( rho, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate rho' )

  DEALLOCATE ( rhop, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate rhop' )

  DEALLOCATE ( rval, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate rval' )

  DEALLOCATE ( rhoval, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate rhoval' )

  DEALLOCATE ( phi, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate phi' )

  DEALLOCATE ( phip, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate phip' )

  DEALLOCATE ( frho, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate frho' )

  DEALLOCATE ( frhop, STAT = isos )
  IF (isos/=0) CALL stop_program('read_eam_data', &
  'failed to deallocate frhop' )
END SUBROUTINE read_eam_data 
!******************************************************************************
END MODULE force_fields
!******************************************************************************
