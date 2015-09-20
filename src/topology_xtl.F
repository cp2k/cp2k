!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Handles XTL (Molecular Simulations, Inc (MSI)) files
!> \author Teodoro Laino [tlaino]
!> \date   05.2009
! *****************************************************************************
MODULE topology_xtl
  USE cell_methods,                    ONLY: write_cell
  USE cell_types,                      ONLY: cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             pbc,&
                                             scaled_to_real,&
                                             set_cell_param
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object,&
                                             parser_search_string
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE input_section_types,             ONLY: section_get_rval,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: matvec_3x3
  USE memory_utilities,                ONLY: reallocate
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_xtl'

  PRIVATE
  PUBLIC :: read_coordinate_xtl

CONTAINS

! *****************************************************************************
!> \brief  Performs the real task of reading the proper information from the XTL
!>         file
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \date   05.2009
!> \par    Format Information implemented:
!>            TITLE
!>            DIMENSION
!>            CELL
!>            SYMMETRY
!>            SYM MAT
!>            ATOMS
!>            EOF
!>
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE read_coordinate_xtl (topology,para_env,subsys_section)
    TYPE(topology_parameters_type)           :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_xtl', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000
    REAL(KIND=dp), PARAMETER                 :: threshold = 1.0E-6_dp

    CHARACTER(LEN=default_string_length)     :: strtmp
    INTEGER                                  :: dimensions, handle, icol, ii, &
                                                isym, iw, jj, natom, &
                                                natom_orig, newsize
    INTEGER, DIMENSION(3)                    :: periodic
    LOGICAL                                  :: check, failure, found, my_end
    REAL(KIND=dp)                            :: pfactor
    REAL(KIND=dp), DIMENSION(3)              :: cell_angles, cell_lengths, r, &
                                                r1, r2, s, transl_vec
    REAL(KIND=dp), DIMENSION(3, 3)           :: rot_mat
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    NULLIFY(parser, logger)
    failure = .FALSE.
    logger => cp_get_default_logger()
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XTL_INFO",&
         extension=".subsysLog")
    CALL timeset(routineN,handle)

    pfactor=section_get_rval(subsys_section,"TOPOLOGY%MEMORY_PROGRESSION_FACTOR")
    ! Element is assigned on the basis of the atm_name
    topology%aa_element = .TRUE.

    atom_info => topology%atom_info
    CALL reallocate(atom_info%id_molname,1,nblock)
    CALL reallocate(atom_info%id_resname,1,nblock)
    CALL reallocate(atom_info%resid,1,nblock)
    CALL reallocate(atom_info%id_atmname,1,nblock)
    CALL reallocate(atom_info%r,1,3,1,nblock)
    CALL reallocate(atom_info%atm_mass,1,nblock)
    CALL reallocate(atom_info%atm_charge,1,nblock)
    CALL reallocate(atom_info%occup,1,nblock)
    CALL reallocate(atom_info%beta,1,nblock)
    CALL reallocate(atom_info%id_element,1,nblock)

    IF (iw>0)  WRITE(iw,*) "    Reading in XTL file ",TRIM(topology%coord_file_name)
    CALL parser_create(parser,topology%coord_file_name,para_env=para_env)

    ! Check for TITLE
    CALL parser_search_string(parser,"TITLE",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.)
    IF (found) THEN
       IF (iw>0) WRITE(iw,'(/,A)')" XTL_INFO| TITLE :: "//TRIM(parser%input_line(parser%icol:))
    END IF

    ! Check for   _chemical_formula_sum
    CALL parser_search_string(parser,"DIMENSION",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.)
    IF (found) THEN
       IF (iw>0) WRITE(iw,'(A)')" XTL_INFO| DIMENSION :: "//TRIM(parser%input_line(parser%icol:))
       CALL parser_get_object(parser,dimensions)
       IF (dimensions/=3) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="XTL file with working DIMENSION different from 3 cannot be parsed! ", &
               error_level=cp_fatal_level)
       END IF
    ELSE
       ! Assuming by default we work in 3D-periodic systems
       dimensions = 3
    END IF

    ! Parsing cell infos
    periodic = 1
    ! Check for   _cell_length_a
    CALL parser_search_string(parser,"CELL",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field CELL was not found in XTL file! "//&
CPSourceFileRef)
    CALL parser_get_next_line(parser,1)
    ! CELL LENGTH A
    CALL parser_get_object(parser,cell_lengths(1))
    cell_lengths(1) = cp_unit_to_cp2k(cell_lengths(1),"angstrom")
    ! CELL LENGTH B
    CALL parser_get_object(parser,cell_lengths(2))
    cell_lengths(2) = cp_unit_to_cp2k(cell_lengths(2),"angstrom")
    ! CELL LENGTH C
    CALL parser_get_object(parser,cell_lengths(3))
    cell_lengths(3) = cp_unit_to_cp2k(cell_lengths(3),"angstrom")

    ! CELL ANGLE ALPHA
    CALL parser_get_object(parser,cell_angles(1))
    cell_angles(1) = cp_unit_to_cp2k(cell_angles(1),"deg")
    ! CELL ANGLE BETA
    CALL parser_get_object(parser,cell_angles(2))
    cell_angles(2) = cp_unit_to_cp2k(cell_angles(2),"deg")
    ! CELL ANGLE GAMMA
    CALL parser_get_object(parser,cell_angles(3))
    cell_angles(3) = cp_unit_to_cp2k(cell_angles(3),"deg")

    ! Create cell
    NULLIFY (cell)
    CALL cell_create(cell)
    CALL set_cell_param(cell,cell_lengths,cell_angles,periodic=periodic,&
                        do_init_cell=.TRUE.)
    CALL write_cell(cell,subsys_section,label="CELL_XTL|")

    ! Parse atoms info and fractional coordinates
    ! Check for   _atom_site_label
    CALL parser_search_string(parser,"ATOMS",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field ATOMS was not found in XTL file! "//&
CPSourceFileRef)
    CALL parser_get_next_line(parser,1)
    ! Paranoic syntax check.. if this fails one should improve the description of XTL files
    found = (INDEX(parser%input_line,"NAME       X          Y          Z")/=0)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field ATOMS in XTL file, is not followed by name and coordinates tags! "//&
CPSourceFileRef)
    CALL parser_get_next_line(parser,1)
    ! Parse real info
    natom = 0
    DO WHILE (INDEX(parser%input_line,"EOF")==0)
       natom = natom + 1
       ! Resize in case needed
       IF (natom > SIZE(atom_info%id_molname)) THEN
          newsize = INT(pfactor*natom)
          CALL reallocate(atom_info%id_molname,1,newsize)
          CALL reallocate(atom_info%id_resname,1,newsize)
          CALL reallocate(atom_info%resid,1,newsize)
          CALL reallocate(atom_info%id_atmname,1,newsize)
          CALL reallocate(atom_info%r,1,3,1,newsize)
          CALL reallocate(atom_info%atm_mass,1,newsize)
          CALL reallocate(atom_info%atm_charge,1,newsize)
          CALL reallocate(atom_info%occup,1,newsize)
          CALL reallocate(atom_info%beta,1,newsize)
          CALL reallocate(atom_info%id_element,1,newsize)
       END IF
       ! NAME
       CALL parser_get_object  (parser,strtmp)
       atom_info%id_atmname(natom) = str2id(strtmp)
       atom_info%id_molname(natom) = str2id(s2s("MOL"//TRIM(ADJUSTL(cp_to_string(natom)))))
       atom_info%id_resname(natom) = atom_info%id_molname(natom)
       atom_info%resid(natom)      = 1
       atom_info%id_element(natom) = atom_info%id_atmname(natom)
       ! X
       CALL parser_get_object(parser,atom_info%r(1,natom))
       ! Y
       CALL parser_get_object(parser,atom_info%r(2,natom))
       ! Z
       CALL parser_get_object(parser,atom_info%r(3,natom))
       s = atom_info%r(1:3,natom)
       CALL scaled_to_real(atom_info%r(1:3,natom),s,cell)
       CALL parser_get_next_line(parser,1,at_end=my_end)
       IF (my_end) EXIT
    END DO
    ! Preliminary check: check if atoms provided are really unique.. this is a paranoic
    ! check since they should be REALLY unique.. anyway..
    DO ii = 1, natom
       r1 = atom_info%r(1:3,ii)
       DO jj = ii+1,natom
          r2 = atom_info%r(1:3,jj)
          r  = pbc(r1-r2,cell)
          check = (SQRT(DOT_PRODUCT(r,r))>=threshold)
          CPPostcondition(check,cp_failure_level,routineP,failure)
       END DO
    END DO
    ! Parse Symmetry Group and generation elements..
    ! Check for SYMMETRY
    CALL parser_search_string(parser,"SYMMETRY",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.)
    IF (found) THEN
       IF (iw>0) WRITE(iw,'(A)')" XTL_INFO| Symmetry Infos :: "//TRIM(parser%input_line(parser%icol:))
    END IF

    ! Check for SYM MAT
    CALL parser_search_string(parser,"SYM MAT",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.)
    CALL cp_assert(found,cp_warning_level,cp_assertion_failed,routineP,&
         "The field SYM MAT was not found in XTL file! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
    IF (iw>0) WRITE(iw,'(A,I0)')" XTL_INFO| Number of atoms before applying symmetry operations :: ",natom
    IF (iw>0) WRITE(iw,'(A10,1X,3F12.6)')(TRIM(id2str(atom_info%id_atmname(ii))),atom_info%r(1:3,ii),ii=1,natom)
    IF (found) THEN
       ! Apply symmetry elements and generate the whole set of atoms in the unit cell
       isym       = 0
       natom_orig = natom
       DO WHILE (found)
          isym = isym + 1
          icol = INDEX(parser%input_line,"SYM MAT")+8
          READ(parser%input_line(icol:),*)((rot_mat(ii,jj),jj=1,3),ii=1,3),transl_vec(1:3)
          Loop_over_unique_atoms: DO ii = 1, natom_orig
             ! Rotate and apply translation
             CALL matvec_3x3(r1, rot_mat, atom_info%r(1:3,ii))
             r1 = r1 + transl_vec
             ! Verify if this atom is really unique..
             check = .TRUE.
             DO jj = 1, natom
                r2 = atom_info%r(1:3,jj)
                r  = pbc(r1-r2,cell)
                IF (SQRT(DOT_PRODUCT(r,r))<=threshold) THEN
                   check = .FALSE.
                   EXIT
                END IF
             END DO
             ! If the atom generated is unique let's add to the atom set..
             IF (check) THEN
                natom = natom + 1
                ! Resize in case needed
                IF (natom > SIZE(atom_info%id_molname)) THEN
                   newsize = INT(pfactor*natom)
                   CALL reallocate(atom_info%id_molname,1,newsize)
                   CALL reallocate(atom_info%id_resname,1,newsize)
                   CALL reallocate(atom_info%resid,1,newsize)
                   CALL reallocate(atom_info%id_atmname,1,newsize)
                   CALL reallocate(atom_info%r,1,3,1,newsize)
                   CALL reallocate(atom_info%atm_mass,1,newsize)
                   CALL reallocate(atom_info%atm_charge,1,newsize)
                   CALL reallocate(atom_info%occup,1,newsize)
                   CALL reallocate(atom_info%beta,1,newsize)
                   CALL reallocate(atom_info%id_element,1,newsize)
                END IF
                atom_info%id_atmname(natom) = atom_info%id_atmname(ii)
                atom_info%id_molname(natom) = atom_info%id_molname(ii)
                atom_info%id_resname(natom) = atom_info%id_resname(ii)
                atom_info%resid(natom)      = atom_info%resid(ii)
                atom_info%id_element(natom) = atom_info%id_element(ii)
                atom_info%r(1:3,natom)      = r1
             END IF
          END DO Loop_over_unique_atoms
          CALL parser_search_string(parser,"SYM MAT",ignore_case=.FALSE.,found=found,&
               begin_line=.FALSE.,search_from_begin_of_file=.FALSE.)
       END DO
    END IF
    IF (iw>0) WRITE(iw,'(A,I0)')" XTL_INFO| Number of symmetry operations :: ",isym
    IF (iw>0) WRITE(iw,'(A,I0)')" XTL_INFO| Number of total atoms :: ",natom
    IF (iw>0) WRITE(iw,'(A10,1X,3F12.6)')(TRIM(id2str(atom_info%id_atmname(ii))),atom_info%r(1:3,ii),ii=1,natom)

    ! Releasse local cell type and parser
    CALL cell_release(cell)
    CALL parser_release(parser)

    ! Reallocate all structures with the exact NATOM size
    CALL reallocate(atom_info%id_molname,1,natom)
    CALL reallocate(atom_info%id_resname,1,natom)
    CALL reallocate(atom_info%resid,1,natom)
    CALL reallocate(atom_info%id_atmname,1,natom)
    CALL reallocate(atom_info%r,1,3,1,natom)
    CALL reallocate(atom_info%atm_mass,1,natom)
    CALL reallocate(atom_info%atm_charge,1,natom)
    CALL reallocate(atom_info%occup,1,natom)
    CALL reallocate(atom_info%beta,1,natom)
    CALL reallocate(atom_info%id_element,1,natom)

    topology%natoms            = natom
    topology%molname_generated = .TRUE.
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/XTL_INFO")
    CALL timestop(handle)
  END SUBROUTINE read_coordinate_xtl


END MODULE topology_xtl
