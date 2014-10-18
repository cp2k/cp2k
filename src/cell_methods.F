!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles all functions related to the CELL
!> \par History
!>      11.2008 Teodoro Laino [tlaino] - deeply cleaning cell_type from units
!>      10.2014 Moved many routines to cell_types.F.
!> \author Matthias KracK (16.01.2002, based on a earlier version of CJM, JGH)
! *****************************************************************************
MODULE cell_methods
  USE cell_types,                      ONLY: &
       cell_clone, cell_create, cell_sym_none, cell_type, get_cell, &
       init_cell, set_cell_param, use_perd_none, use_perd_x, use_perd_xy, &
       use_perd_xyz, use_perd_xz, use_perd_y, use_perd_yz, use_perd_z
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                             cp_unit_to_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: do_cell_cp2k,&
                                             do_cell_xsc
  USE input_cp2k_subsys,               ONLY: create_cell_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: &
       section_get_keyword, section_release, section_type, section_vals_get, &
       section_vals_get_subs_vals, section_vals_type, section_vals_val_get, &
       section_vals_val_set, section_vals_val_unset
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: degree
  USE mathlib,                         ONLY: angle
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_methods'


  ! Public subroutines
  PUBLIC :: read_cell, write_cell

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param cell_section ...
!> \param check_for_ref ...
!> \param para_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2005 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE read_cell( cell, cell_ref, use_ref_cell, cell_section,&
    check_for_ref, para_env, error)
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(OUT), OPTIONAL           :: use_ref_cell
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: cell_section
    LOGICAL, INTENT(IN), OPTIONAL            :: check_for_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_cell', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_per
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    LOGICAL :: cell_read_a, cell_read_abc, cell_read_b, cell_read_c, &
      cell_read_file, check, failure, my_check
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_angles, cell_par
    TYPE(section_vals_type), POINTER         :: cell_ref_section

    failure  = .FALSE.
    my_check = .TRUE.
    NULLIFY(cell_ref_section, cell_par, multiple_unit_cell)
    IF (.NOT.failure) THEN
       IF (.NOT.ASSOCIATED(cell)) CALL cell_create(cell,error=error)
       IF (.NOT.ASSOCIATED(cell_ref)) CALL cell_create(cell_ref,error=error)
       IF (PRESENT(check_for_ref)) my_check = check_for_ref

       cell%deth         = 0.0_dp
       cell%orthorhombic = .FALSE.
       cell%perd(:)      = 1
       cell%symmetry_id  = cell_sym_none
       cell%hmat(:,:)    = 0.0_dp
       cell%h_inv(:,:)   = 0.0_dp
       cell_read_file    = .FALSE.
       cell_read_a       = .FALSE.
       cell_read_b       = .FALSE.
       cell_read_c       = .FALSE.
       ! Trying to read cell info from file
       CALL section_vals_val_get(cell_section,"CELL_FILE_NAME",explicit=cell_read_file,error=error)
       IF (cell_read_file) CALL read_cell_from_external_file(cell_section, para_env, error)

       ! Trying to read cell info from the separate A, B, C vectors
       ! If cell information is provided through file A,B,C contain the file information..
       ! a print warning is shown on screen..
       CALL section_vals_val_get(cell_section,"A",explicit=cell_read_a,error=error)
       IF (cell_read_a) THEN
          CALL section_vals_val_get(cell_section,"A",r_vals=cell_par,error=error)
          cell%hmat(:,1) = cell_par(:)
       END IF
       CALL section_vals_val_get(cell_section,"B",explicit=cell_read_b,error=error)
       IF (cell_read_b) THEN
          CALL section_vals_val_get(cell_section,"B",r_vals=cell_par,error=error)
          cell%hmat(:,2) = cell_par(:)
       END IF
       CALL section_vals_val_get(cell_section,"C",explicit=cell_read_c,error=error)
       IF (cell_read_c) THEN
          CALL section_vals_val_get(cell_section,"C",r_vals=cell_par,error=error)
          cell%hmat(:,3) = cell_par(:)
       END IF
       check = ((cell_read_a.EQV.cell_read_b).AND.(cell_read_b.EQV.cell_read_c))
       CALL cp_assert(check,cp_warning_level,cp_assertion_failed,routineP,&
            "Cell Information provided through vectors A, B and C. Not all three "//&
            "vectors were provided! Cell setup may be incomplete!"//&
CPSourceFileRef,&
            only_ionode=.TRUE.)

       ! Very last option.. Trying to read cell info from ABC keyword
       CALL section_vals_val_get(cell_section,"ABC",explicit=cell_read_abc,error=error)
       IF (cell_read_abc) THEN
          check = (cell_read_a.OR.cell_read_b.OR.cell_read_c)
          CALL cp_assert(.NOT.check,cp_warning_level,cp_assertion_failed,routineP,&
               "Cell Information provided through vectors A, B and C in conjunction with ABC."//&
               " The definition of the ABC keyword will override the one provided by A,B and C."//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          cell%hmat = 0.0_dp
          CALL section_vals_val_get(cell_section,"ABC",r_vals=cell_par,error=error)
          CALL section_vals_val_get(cell_section,"ALPHA_BETA_GAMMA",r_vals=cell_angles,error=error)
          CALL set_cell_param(cell,cell_par,cell_angles,do_init_cell=.FALSE.,error=error)
       END IF

       ! Multiple unit cell
       CALL section_vals_val_get(cell_section,"MULTIPLE_UNIT_CELL",i_vals=multiple_unit_cell,error=error)
       IF (ANY(multiple_unit_cell/=1)) CALL set_multiple_unit_cell(cell, multiple_unit_cell, error)

       CALL section_vals_val_get(cell_section,"PERIODIC",i_val=my_per,error=error)
       SELECT CASE(my_per)
       CASE(use_perd_x)
          cell%perd = (/1,0,0/)
       CASE(use_perd_y)
          cell%perd = (/0,1,0/)
       CASE(use_perd_z)
          cell%perd = (/0,0,1/)
       CASE(use_perd_xy)
          cell%perd = (/1,1,0/)
       CASE(use_perd_xz)
          cell%perd = (/1,0,1/)
       CASE(use_perd_yz)
          cell%perd = (/0,1,1/)
       CASE(use_perd_xyz)
          cell%perd = (/1,1,1/)
       CASE(use_perd_none)
          cell%perd = (/0,0,0/)
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       ! Load requested cell symmetry
       CALL section_vals_val_get(cell_section,"SYMMETRY",i_val=cell%symmetry_id,error=error)

       ! Initialize cell
       CALL init_cell(cell)

       IF (.NOT.my_check) RETURN
       cell_ref_section => section_vals_get_subs_vals(cell_section,&
            "CELL_REF",error=error)
       IF (parsed_cp2k_input(cell_ref_section,check_this_section=.TRUE.,error=error)) THEN
          IF(PRESENT(use_ref_cell) ) use_ref_cell = .TRUE.
          CALL read_cell(cell_ref, cell_ref, use_ref_cell, cell_section=cell_ref_section,&
               check_for_ref=.FALSE., para_env=para_env, error=error)
       ELSE
          CALL cell_clone (cell, cell_ref, error)
          IF ( PRESENT ( use_ref_cell ) ) use_ref_cell = .FALSE.
       END IF
    END IF

  END SUBROUTINE read_cell

! *****************************************************************************
!> \brief utility function to ease the transition to the new input.
!>      returns true if the new input was parsed
!> \param input_file the parsed input file
!> \param check_this_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION parsed_cp2k_input(input_file,check_this_section,error) RESULT(res)
    TYPE(section_vals_type), POINTER         :: input_file
    LOGICAL, INTENT(IN), OPTIONAL            :: check_this_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'parsed_cp2k_input', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_check
    TYPE(section_vals_type), POINTER         :: glob_section

    failure=.FALSE.
    my_check = .FALSE.
    IF (PRESENT(check_this_section)) my_check = check_this_section
    res=ASSOCIATED(input_file)
    IF (res) THEN
       CPPrecondition(input_file%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          IF (.NOT.my_check) THEN
             glob_section => section_vals_get_subs_vals(input_file,"GLOBAL",&
                  error=error)
             CALL section_vals_get(glob_section,explicit=res,error=error)
          ELSE
             CALL section_vals_get(input_file,explicit=res,error=error)
          END IF
       END IF
    END IF
  END FUNCTION parsed_cp2k_input

! *****************************************************************************
!> \brief   Setup of the multiple unit_cell
!> \param cell ...
!> \param multiple_unit_cell ...
!> \param error ...
!> \date    05.2009
!> \author  Teodoro Laino [tlaino]
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_multiple_unit_cell(cell, multiple_unit_cell, error)

    TYPE(cell_type), POINTER                 :: cell
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_multiple_unit_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.

    ! Fail is one of the value is set to zero..
    CALL cp_assert(ALL(multiple_unit_cell>0),cp_fatal_level,cp_assertion_failed,routineP,&
         "CELL%MULTIPLE_UNIT_CELL accepts only integer values larger than 0! "//&
         "A value of 0 or negative is meaningless!"//&
CPSourceFileRef)

    ! scale abc accordingly user request
    cell%hmat(:,1) = cell%hmat(:,1)*multiple_unit_cell(1)
    cell%hmat(:,2) = cell%hmat(:,2)*multiple_unit_cell(2)
    cell%hmat(:,3) = cell%hmat(:,3)*multiple_unit_cell(3)

  END SUBROUTINE set_multiple_unit_cell

! *****************************************************************************
!> \brief   Read cell information from an external file
!> \param cell_section ...
!> \param para_env ...
!> \param error ...
!> \date    02.2008
!> \author  Teodoro Laino [tlaino] - University of Zurich
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_cell_from_external_file(cell_section, para_env, error)

    TYPE(section_vals_type), POINTER         :: cell_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_cell_from_external_file', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: cell_file_name
    INTEGER                                  :: i, idum, j, my_format, n_rep, &
                                                stat
    LOGICAL                                  :: explicit, failure, my_end
    REAL(KIND=dp)                            :: xdum
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_par
    TYPE(cp_parser_type), POINTER            :: parser

    failure = .FALSE.
    NULLIFY(parser)
    CALL section_vals_val_get(cell_section,"CELL_FILE_NAME",c_val=cell_file_name,error=error)
    CALL section_vals_val_get(cell_section,"CELL_FILE_FORMAT",i_val=my_format, error=error)
    CALL parser_create(parser,cell_file_name, para_env=para_env,error=error)
    CALL parser_get_next_line(parser,1,error=error)
    SELECT CASE(my_format)
    CASE (do_cell_cp2k)
       my_end = .FALSE.
       DO WHILE (.NOT.my_end)
          READ(parser%input_line,*)idum,xdum,hmat(:,1),hmat(:,2),hmat(:,3)
          CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
       END DO
    CASE (do_cell_xsc)
       READ(parser%input_line,*)idum,hmat(:,1),hmat(:,2),hmat(:,3)
    END SELECT
    CALL parser_release(parser,error=error)
    CALL section_vals_val_unset(cell_section,"CELL_FILE_NAME",error=error)
    CALL section_vals_val_unset(cell_section,"CELL_FILE_FORMAT",error=error)
    ! Conver to CP2K units
    DO i = 1, 3
       DO j = 1, 3
          hmat(j,i) = cp_unit_to_cp2k(hmat(j,i), "angstrom", error=error)
       END DO
    END DO
    ! Check if the cell was already defined
    explicit = .FALSE.
    CALL section_vals_val_get(cell_section,"A",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    CALL section_vals_val_get(cell_section,"B",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    CALL section_vals_val_get(cell_section,"C",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    CALL section_vals_val_get(cell_section,"ABC",n_rep_val=n_rep,error=error)
    explicit = explicit .OR. (n_rep==1)
    ! Possibly print a warning
    CALL cp_assert(.NOT.explicit,cp_warning_level,cp_assertion_failed,routineP,&
         "Cell specification (A,B,C or ABC) provided together with the external "//&
         "cell setup! Ignoring (A,B,C or ABC) and proceeding with info read from the "//&
         "external file! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
    ! Copy cell information in the A, B, C fields..(we may need them later on..)
    ALLOCATE(cell_par(3), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    cell_par = hmat(:,1)
    CALL section_vals_val_set(cell_section,"A",r_vals_ptr=cell_par,error=error)
    ALLOCATE(cell_par(3), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    cell_par = hmat(:,2)
    CALL section_vals_val_set(cell_section,"B",r_vals_ptr=cell_par,error=error)
    ALLOCATE(cell_par(3), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    cell_par = hmat(:,3)
    CALL section_vals_val_set(cell_section,"C",r_vals_ptr=cell_par,error=error)
    ! Unset possible keywords
    CALL section_vals_val_unset(cell_section,"ABC",error=error)
    CALL section_vals_val_unset(cell_section,"ALPHA_BETA_GAMMA",error=error)

  END SUBROUTINE read_cell_from_external_file

! *****************************************************************************
!> \brief   Write the cell parameters to the output unit.
!> \param cell ...
!> \param subsys_section ...
!> \param cell_ref ...
!> \param label ...
!> \param error ...
!> \date    02.06.2000
!> \par     History
!>    - 11.2008 Teodoro Laino [tlaino] - rewrite and enabling user driven units
!> \author  Matthias Krack
!> \version 1.0
! *****************************************************************************
  RECURSIVE SUBROUTINE write_cell(cell,subsys_section,cell_ref,label,error)

    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cell_type), OPTIONAL, POINTER       :: cell_ref
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: label
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_cell', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_label, unit_str
    INTEGER                                  :: output_unit
    REAL(KIND=dp)                            :: alpha, beta, gamma, val
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section

    NULLIFY (enum)
    NULLIFY (keyword)
    NULLIFY (logger)
    NULLIFY (section)
    logger => cp_error_get_logger(error)
    my_label = "CELL|"
    IF (PRESENT(label)) my_label = TRIM(label)
    output_unit = cp_print_key_unit_nr(logger,subsys_section,"PRINT%CELL",&
                                       extension=".Log",error=error)
    CALL section_vals_val_get(subsys_section,"PRINT%CELL%UNIT",c_val=unit_str,error=error)
    IF (output_unit > 0) THEN
       CALL get_cell(cell=cell,abc=abc,alpha=alpha,beta=beta,gamma=gamma)
       WRITE (UNIT=output_unit, FMT='( )')
       val = cp_unit_from_cp2k(cell%deth,TRIM(unit_str)//"^3",error=error)
       WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.3)")&
        TRIM(my_label)//" Volume ["//TRIM(unit_str)//"^3]:",val
       val = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
       WRITE (UNIT=output_unit,FMT="(T2,A,T30,3F10.3,4X,A6,F11.3)")&
        TRIM(my_label)//" Vector a ["//TRIM(unit_str)//"]:",cell%hmat(:,1)*val,&
        "|a| = ",abc(1)*val,&
        TRIM(my_label)//" Vector b ["//TRIM(unit_str)//"]:",cell%hmat(:,2)*val,&
        "|b| = ",abc(2)*val,&
        TRIM(my_label)//" Vector c ["//TRIM(unit_str)//"]:",cell%hmat(:,3)*val,&
        "|c| = ",abc(3)*val
       WRITE (UNIT=output_unit,FMT="(T2,A,T70,F11.3)")&
        TRIM(my_label)//" Angle (b,c), alpha [degree]: ",alpha,&
        TRIM(my_label)//" Angle (a,c), beta  [degree]: ",beta,&
        TRIM(my_label)//" Angle (a,b), gamma [degree]: ",gamma
        IF (cell%symmetry_id /= cell_sym_none) THEN
           CALL create_cell_section(section,error=error)
           keyword => section_get_keyword(section,"SYMMETRY",error=error)
           CALL keyword_get(keyword,enum=enum,error=error)
           WRITE (UNIT=output_unit,FMT="(T2,A,T61,A20)")&
            TRIM(my_label)//" Requested initial symmetry: ",&
            ADJUSTR(TRIM(enum_i2c(enum,cell%symmetry_id,error=error)))
           CALL section_release(section,error=error)
        END IF
        IF (cell%orthorhombic) THEN
           WRITE (UNIT=output_unit,FMT="(T2,A,T78,A3)")&
            TRIM(my_label)//" Numerically orthorhombic: ","YES"
        ELSE
           WRITE (UNIT=output_unit,FMT="(T2,A,T78,A3)")&
            TRIM(my_label)//" Numerically orthorhombic: "," NO"
        END IF
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%CELL",error=error)

    IF (PRESENT(cell_ref)) THEN
       CALL write_cell(cell_ref, subsys_section, label="CELL_REF|", error=error)
    END IF

  END SUBROUTINE write_cell

END MODULE cell_methods
