!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
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

  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE pair_potential_types,            ONLY: reallocate_potparm
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE string_utilities,                ONLY: str_search,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE topology_types,                  ONLY: topology_parameters_type,&
                                             constraint_info_type,&
                                             force_field_amber_info_type,&
                                             force_field_charmm_info_type,&
                                             force_field_input_info_type


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_force_field_section,&
            read_force_field_charmm,&
            read_force_field_amber

!!*****
!******************************************************************************

CONTAINS

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

SUBROUTINE read_force_field_section ( topology , globenv )


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    INTEGER                                  :: iw

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

!..read in non-bonded interaction information
    CALL input_nonbonded(topology,globenv)

!..read the charges parameters
    CALL input_charges(topology,globenv)

!..read the bond parameters
    CALL input_bonds(topology,globenv)

!..read the bend parameters
    CALL input_bends(topology,globenv)

!..read the torsion parameters
    CALL input_torsions(topology,globenv)

!..read the distance constraint parameters
    CALL input_constraints(topology,globenv)

    CALL parser_end()
 
  IF(PM1) WRITE(iw,*) "  Exiting  read_force_field_section"

END SUBROUTINE read_force_field_section

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
!!I end nonbonded                                                             I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_nonbonded(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(force_field_input_info_type),POINTER       :: ffin_info
    CHARACTER(LEN=default_string_length)     :: at1, at2
    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length)     :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw, &
                                                np
    REAL(dbl)                                :: a, b, c, d, dc, epsilon, m, &
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

  ffin_info => topology%ffin_info

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
                 CALL reallocate_potparm(ffin_info%ff_nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( epsilon )
                 CALL get_next ( sigma )
                 CALL get_next ( rcut )
                 ffin_info%ff_nonbonded(np) %type = string2
                 ffin_info%ff_nonbonded(np) %at1 = at1
                 ffin_info%ff_nonbonded(np) %at2 = at2
                 ffin_info%ff_nonbonded(np) %lj%epsilon = epsilon
                 ffin_info%ff_nonbonded(np) %lj%sigma6 = sigma**6
                 ffin_info%ff_nonbonded(np) %lj%sigma12 = sigma**12
                 ffin_info%ff_nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'WILLIAMS')
                 np = np + 1
                 CALL reallocate_potparm(ffin_info%ff_nonbonded,1,np)
                 ilen = 0
                 CALL get_next(at1,ilen)
                 ilen = 0
                 CALL get_next(at2,ilen)
                 CALL get_next ( a )
                 CALL get_next ( b )
                 CALL get_next ( c )
                 CALL get_next ( rcut )
                 ffin_info%ff_nonbonded(np) %type = string2
                 ffin_info%ff_nonbonded(np) %at1 = at1
                 ffin_info%ff_nonbonded(np) %at2 = at2
                 ffin_info%ff_nonbonded(np) %willis%a = a
                 ffin_info%ff_nonbonded(np) %willis%b = b
                 ffin_info%ff_nonbonded(np) %willis%c = c
                 ffin_info%ff_nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'GOODWIN')
                 np = np + 1
                 CALL reallocate_potparm(ffin_info%ff_nonbonded,1,np)
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
                 ffin_info%ff_nonbonded(np) %type = string2
                 ffin_info%ff_nonbonded(np) %at1 = at1
                 ffin_info%ff_nonbonded(np) %at2 = at2
                 ffin_info%ff_nonbonded(np) %goodwin%vr0 = vr0
                 ffin_info%ff_nonbonded(np) %goodwin%d = d
                 ffin_info%ff_nonbonded(np) %goodwin%dc = dc
                 ffin_info%ff_nonbonded(np) %goodwin%m = m
                 ffin_info%ff_nonbonded(np) %goodwin%mc = mc
                 ffin_info%ff_nonbonded(np) %rcutsq = rcut*rcut
              CASE ( 'END')
                 ilen = 0
                 CALL get_next(string2,ilen)
                 EXIT
              END SELECT
           END DO
        END IF
        CALL read_line
     END DO
     CALL convert_to_cp2k_units ( "KELVIN", ff_nonbonded = ffin_info % ff_nonbonded )
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

SUBROUTINE input_charges(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(force_field_input_info_type),POINTER       :: ffin_info
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

  ffin_info => topology%ffin_info

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
              WRITE(*,*) "string2 ",string2
              IF(string2=="END") EXIT

              i = i + 1
              CALL reallocate(ffin_info%ff_charge_atm,1,i)
              CALL reallocate(ffin_info%ff_charge,1,i)

              READ(string2(:),*) ffin_info%ff_charge_atm(i)
              CALL get_next ( ffin_info%ff_charge(i) )

              WRITE(*,*) "CCCCHARGE ",ffin_info%ff_charge_atm(i)," ",ffin_info%ff_charge(i)

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
!!     we (so far) will do nothing with this information. It appears that
!!     all major parameters sets, such as CHARMM or AMBER, all have only
!!     one type of bonded interaction (harmonic).
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
!!I    harmonic ATMNAE1 ATMNAME2 k r0
!!I end bonds                                                                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_bonds(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(force_field_input_info_type),POINTER       :: ffin_info
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

  ffin_info => topology%ffin_info

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
              WRITE(*,*) "string2 ",string2
              SELECT CASE (string2)
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser( 'input_bond','unknown option')
              CASE ( 'HARMONIC')
                 i = i + 1
                 CALL reallocate(ffin_info%ff_bond_a,1,i)
                 CALL reallocate(ffin_info%ff_bond_b,1,i)
                 CALL reallocate(ffin_info%ff_bond_k,1,i)
                 CALL reallocate(ffin_info%ff_bond_r0,1,i)

                 ilen = default_string_length
                 CALL get_next ( ffin_info%ff_bond_a(i), ilen )
                 CALL get_next ( ffin_info%ff_bond_b(i), ilen )
                 CALL get_next ( ffin_info%ff_bond_k(i) )
                 CALL get_next ( ffin_info%ff_bond_r0(i) )
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

SUBROUTINE input_bends(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(force_field_input_info_type),POINTER       :: ffin_info
    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string, string2
    INTEGER                                  :: i, ierror, ilen, iw
    REAL(dbl)                                :: todeg, torad

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

  ffin_info => topology%ffin_info

  torad = pi / 180.0_dbl
  todeg = 180.0_dbl/pi

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
                 CALL reallocate(ffin_info%ff_bend_a,1,i)
                 CALL reallocate(ffin_info%ff_bend_b,1,i)
                 CALL reallocate(ffin_info%ff_bend_c,1,i)
                 CALL reallocate(ffin_info%ff_bend_k,1,i)
                 CALL reallocate(ffin_info%ff_bend_theta0,1,i)

                 ilen = default_string_length
                 CALL get_next ( ffin_info%ff_bend_a(i), ilen )
                 CALL get_next ( ffin_info%ff_bend_b(i), ilen )
                 CALL get_next ( ffin_info%ff_bend_c(i), ilen )
                 CALL get_next ( ffin_info%ff_bend_k(i) )
                 CALL get_next ( ffin_info%ff_bend_theta0(i) )
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

SUBROUTINE input_torsions(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(force_field_input_info_type),POINTER       :: ffin_info

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

  ffin_info => topology%ffin_info

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
                 CALL reallocate(ffin_info%ff_torsion_a,1,i)
                 CALL reallocate(ffin_info%ff_torsion_b,1,i)
                 CALL reallocate(ffin_info%ff_torsion_c,1,i)
                 CALL reallocate(ffin_info%ff_torsion_d,1,i)
                 CALL reallocate(ffin_info%ff_torsion_k,1,i)
                 CALL reallocate(ffin_info%ff_torsion_phi0,1,i)

                 ilen = default_string_length
                 CALL get_next ( ffin_info%ff_torsion_a(i), ilen )
                 CALL get_next ( ffin_info%ff_torsion_b(i), ilen )
                 CALL get_next ( ffin_info%ff_torsion_c(i), ilen )
                 CALL get_next ( ffin_info%ff_torsion_d(i), ilen )
                 CALL get_next ( ffin_info%ff_torsion_k(i) )
                 CALL get_next ( ffin_info%ff_torsion_phi0(i) )
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
!!****** force_fields/input_constraints [1.0] *
!!
!!   NAME
!!     input_constraints
!!
!!   FUNCTION
!!     Read all the distance parameters. Put them in the
!!     constraint_distance array.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) Distance parameters are now stored in tables. The position
!!        within the table is used as handle for the topology
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &force_field ... &end                                            I
!!I                                                                           I
!!I constraints                                                               I
!!I   HYDROGEN                                                                I
!!I   ATOM ??????                                                             I
!!I   DIST MOLKIND ATM_A ATM_B D_AB                                           I
!!I   G3x3 MOLKIND ATM_A ATM_B ATM_C D_AB D_AC D_BC                           I
!!I                                                                           I
!!I end constraints                                                           I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE input_constraints(topology,globenv)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(constraint_info_type),POINTER       :: cons_info
    CHARACTER(LEN=12)                        :: label
    CHARACTER(LEN=40)                        :: string
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

  cons_info => topology%cons_info


!..parse the input section
  label = '&FORCE_FIELD'
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF (ierror == 0 ) THEN
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 7
        CALL get_next ( string, ilen )
        CALL uppercase ( string )

        IF ( INDEX ( string, 'ATOM' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_atom=.TRUE.
           CALL stop_program ( "input_constraints", &
                               "CONATOM not done yet " )
        ELSE IF ( INDEX ( string, 'CONHYDR' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_hydr=.TRUE.
           CALL stop_program ( "input_constraints", &
               "CONHYDR will not work until force field is completed " )
        ELSE IF ( INDEX ( string, 'DIST' ) /= 0 ) THEN
           CALL read_line()
           topology%constraint=.TRUE.
           topology%const_dist=.TRUE.
           ilen = 40
           CALL get_next(string,ilen)
           CALL uppercase(string)
           i=cons_info%nconst_dist + 1
           CALL reallocate(cons_info%const_dist_mol,1,i)
           CALL reallocate(cons_info%const_dist_a,1,i)
           CALL reallocate(cons_info%const_dist_b,1,i)
           CALL reallocate(cons_info%const_dist_dab,1,i)
           CALL get_next(cons_info%const_dist_mol(i))
           CALL get_next(cons_info%const_dist_a(i))
           CALL get_next(cons_info%const_dist_b(i))
           CALL get_next(cons_info%const_dist_dab(i))
           cons_info%nconst_dist = i
        ELSE IF ( INDEX ( string, 'G3X3' ) /= 0 ) THEN
           topology%constraint=.TRUE.
           topology%const_33=.TRUE.
           ilen = 40
           CALL get_next(string,ilen)
           CALL uppercase(string)
           i=cons_info%nconst_g33 + 1
           CALL reallocate(cons_info%const_g33_mol,1,i)
           CALL reallocate(cons_info%const_g33_a,1,i)
           CALL reallocate(cons_info%const_g33_b,1,i)
           CALL reallocate(cons_info%const_g33_c,1,i)
           CALL reallocate(cons_info%const_g33_dab,1,i)
           CALL reallocate(cons_info%const_g33_dac,1,i)
           CALL reallocate(cons_info%const_g33_dbc,1,i)
           CALL get_next(cons_info%const_g33_mol(i))
           CALL get_next(cons_info%const_g33_a(i))
           CALL get_next(cons_info%const_g33_b(i))
           CALL get_next(cons_info%const_g33_c(i))
           CALL get_next(cons_info%const_g33_dab(i))
           CALL get_next(cons_info%const_g33_dac(i))
           CALL get_next(cons_info%const_g33_dbc(i))
           cons_info%nconst_g33 = i
        END IF
        CALL read_line
     END DO
  END IF

END SUBROUTINE input_constraints

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

SUBROUTINE read_force_field_charmm ( topology , globenv )


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    TYPE(force_field_charmm_info_type),POINTER       :: ffch_info

  LOGICAL                                  :: PM1, PM2, PM3, PM4
    INTEGER                                  :: i, istat, ilen, iw
    INTEGER                                  :: nbond,nbend,ntorsion,nnonbond

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

  ffch_info => topology%ffch_info
  CALL parser_init(topology%ff_file_name,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Read in all the Bonds info from the param file here
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
       CALL reallocate(ffch_info%ff_bond_a,1,nbond)
       CALL reallocate(ffch_info%ff_bond_b,1,nbond)
       CALL reallocate(ffch_info%ff_bond_k,1,nbond)
       CALL reallocate(ffch_info%ff_bond_r0,1,nbond)
       ffch_info%ff_bond_a(nbond) = string
       ffch_info%ff_bond_b(nbond) = string2
       CALL get_next(ffch_info%ff_bond_k(nbond))
       CALL get_next(ffch_info%ff_bond_r0(nbond))
       IF(PM1) WRITE(iw,*) "    CHM BOND ",nbond,&
                           TRIM(ffch_info%ff_bond_a(nbond))," ",&
                           TRIM(ffch_info%ff_bond_b(nbond))," ",&
                           ffch_info%ff_bond_k(nbond),&
                           ffch_info%ff_bond_r0(nbond)
       CALL read_line
     END DO
  END IF

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in all the Bends info from the param file here
  !-----------------------------------------------------------------------------
  nbend = 0
  label = 'ANGLES'
  CALL search_label ( label, istat, ignore_case=.TRUE. )
  IF (istat == 0) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 7
       CALL get_next(string,ilen)
       CALL uppercase ( string )
       IF(string == "DIHEDRALS") EXIT
       CALL get_next(string2,ilen)
       CALL get_next(string3,ilen)
       nbend = nbend + 1
       CALL reallocate(ffch_info%ff_bend_a,1,nbend)
       CALL reallocate(ffch_info%ff_bend_b,1,nbend)
       CALL reallocate(ffch_info%ff_bend_c,1,nbend)
       CALL reallocate(ffch_info%ff_bend_k,1,nbend)
       CALL reallocate(ffch_info%ff_bend_theta0,1,nbend)
       ffch_info%ff_bend_a(nbend) = string
       ffch_info%ff_bend_b(nbend) = string2
       ffch_info%ff_bend_c(nbend) = string3
       CALL get_next(ffch_info%ff_bend_k(nbend))
       CALL get_next(ffch_info%ff_bend_theta0(nbend))
       IF(PM1) WRITE(iw,*) "    CHM BEND ",nbend,&
                           TRIM(ffch_info%ff_bend_a(nbend))," ",&
                           TRIM(ffch_info%ff_bend_b(nbend))," ",&
                           TRIM(ffch_info%ff_bend_c(nbend))," ",&
                           ffch_info%ff_bend_k(nbend),&
                           ffch_info%ff_bend_theta0(nbend)
       CALL read_line
     END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Read in all the Dihedrals info from the param file here
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
       IF(string == "NONBONDED") EXIT
       CALL get_next(string2,ilen)
       CALL get_next(string3,ilen)
       CALL get_next(string4,ilen)
       ntorsion = ntorsion + 1
       CALL reallocate(ffch_info%ff_torsion_a,1,ntorsion)
       CALL reallocate(ffch_info%ff_torsion_b,1,ntorsion)
       CALL reallocate(ffch_info%ff_torsion_c,1,ntorsion)
       CALL reallocate(ffch_info%ff_torsion_c,1,ntorsion)
       CALL reallocate(ffch_info%ff_torsion_k,1,ntorsion)
       CALL reallocate(ffch_info%ff_torsion_phi0,1,ntorsion)
       ffch_info%ff_torsion_a(ntorsion) = string
       ffch_info%ff_torsion_b(ntorsion) = string2
       ffch_info%ff_torsion_c(ntorsion) = string3
       ffch_info%ff_torsion_d(ntorsion) = string4
       CALL get_next(ffch_info%ff_torsion_k(ntorsion))
       CALL get_next(ffch_info%ff_torsion_phi0(ntorsion))
       IF(PM1) WRITE(iw,*) "    CHM TORSION ",ntorsion,&
                           TRIM(ffch_info%ff_torsion_a(ntorsion))," ",&
                           TRIM(ffch_info%ff_torsion_b(ntorsion))," ",&
                           TRIM(ffch_info%ff_torsion_c(ntorsion))," ",&
                           TRIM(ffch_info%ff_torsion_d(ntorsion))," ",&
                           ffch_info%ff_torsion_k(ntorsion),&
                           ffch_info%ff_torsion_phi0(ntorsion)
       CALL read_line
     END DO
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Read in all the Nonbonded info from the param file here
  !-----------------------------------------------------------------------------
  nnonbond = 0
  label = 'NONBONDED'
  CALL search_label ( label, istat, ignore_case=.TRUE. )
  IF (istat == 0) THEN
     CALL read_line
     DO WHILE (test_next()/='X')
       ilen = 7
       CALL get_next(string,ilen)
       CALL uppercase ( string )
       IF(string == "HBOND") EXIT
       nnonbond = nnonbond + 1
       CALL reallocate(ffch_info%ff_nonbond_a,1,nnonbond)
       CALL reallocate(ffch_info%ff_nonbond_eps,1,nnonbond)
       CALL reallocate(ffch_info%ff_nonbond_rmin2,1,nnonbond)
       CALL reallocate(ffch_info%ff_nonbond_eps_14,1,nnonbond)
       CALL reallocate(ffch_info%ff_nonbond_rmin2_14,1,nnonbond)
       ffch_info%ff_nonbond_a(nnonbond) = string
       CALL get_next(ffch_info%ff_nonbond_eps(nnonbond))
       CALL get_next(ffch_info%ff_nonbond_rmin2(nnonbond))
       CALL get_next(ffch_info%ff_nonbond_eps_14(nnonbond))
       CALL get_next(ffch_info%ff_nonbond_rmin2_14(nnonbond))
       IF(PM1) WRITE(iw,*) "    CHM NONBOND ",nnonbond,&
                           TRIM(ffch_info%ff_nonbond_a(nnonbond))," ",&
                           ffch_info%ff_nonbond_eps(nnonbond),&
                           ffch_info%ff_nonbond_rmin2(nnonbond),&
                           ffch_info%ff_nonbond_eps_14(nnonbond),&
                           ffch_info%ff_nonbond_rmin2_14(nnonbond)
       CALL read_line
     END DO
  END IF


  CALL parser_end


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

SUBROUTINE read_force_field_amber ( topology , globenv )


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN),OPTIONAL                    :: globenv

    !TYPE(force_field_input_info_type),POINTER       :: ffin_info

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


END MODULE force_fields

!******************************************************************************
