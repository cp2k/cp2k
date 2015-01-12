!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Handles CIF (Crystallographic Information File) files
!> \author Teodoro Laino [tlaino]
!> \date   12.2008
! *****************************************************************************
MODULE topology_cif
  USE cell_methods,                    ONLY: write_cell
  USE cell_types,                      ONLY: cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             pbc,&
                                             real_to_scaled,&
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
  USE fparser,                         ONLY: evalf,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_section_types,             ONLY: section_get_rval,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE string_utilities,                ONLY: s2a
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_cif'

  PRIVATE
  PUBLIC :: read_coordinate_cif

CONTAINS

! *****************************************************************************
!> \brief  Performs the real task of reading the proper information from the CIF
!>         file
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \param error ...
!> \date   12.2008
!> \par    Format Information implemented:
!>            _chemical_name
!>            _chemical_formula_sum
!>            _cell_length_a
!>            _cell_length_b
!>            _cell_length_c
!>            _cell_angle_alpha
!>            _cell_angle_beta
!>            _cell_angle_gamma
!>            _symmetry_space_group_name_h-m
!>            _symmetry_equiv_pos_as_xyz
!>            _space_group_symop_operation_xyz
!>            _atom_site_label
!>            _atom_site_type_symbol
!>            _atom_site_fract_x
!>            _atom_site_fract_y
!>            _atom_site_fract_z
!>
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE read_coordinate_cif (topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type)           :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_cif', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock = 1000
    REAL(KIND=dp), PARAMETER                 :: threshold = 1.0E-3_dp

    CHARACTER(LEN=default_string_length)     :: s_tag, strtmp
    INTEGER                                  :: handle, ii, iln0, iln1, iln2, &
                                                iln3, isym, itype, iw, jj, &
                                                natom, natom_orig, newsize
    INTEGER, DIMENSION(20)                   :: ifield
    INTEGER, DIMENSION(3)                    :: periodic
    LOGICAL                                  :: check, failure, found, my_end
    REAL(KIND=dp)                            :: pfactor
    REAL(KIND=dp), DIMENSION(3)              :: cell_angles, cell_lengths, r, &
                                                r1, r2, s, s_tmp
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    NULLIFY(parser, logger)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/CIF_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)
    pfactor=section_get_rval(subsys_section,"TOPOLOGY%MEMORY_PROGRESSION_FACTOR",error)

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

    IF (iw>0)  WRITE(iw,*) "    Reading in CIF file ",TRIM(topology%coord_file_name)
    CALL parser_create(parser,topology%coord_file_name,&
         para_env=para_env,apply_preprocessing=.FALSE.,error=error)

    ! Check for   _chemical_name
    CALL parser_search_string(parser,"_chemical_name",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    IF (found) THEN
       IF (iw>0) WRITE(iw,'(/,A)')" CIF_INFO| _chemical_name :: "//TRIM(parser%input_line(parser%icol:))
    END IF

    ! Check for   _chemical_formula_sum
    CALL parser_search_string(parser,"_chemical_formula_sum",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    IF (found) THEN
       IF (iw>0) WRITE(iw,'(A)')" CIF_INFO| _chemical_formula_sum :: "//TRIM(parser%input_line(parser%icol:))
    END IF

    ! Parsing cell infos
    periodic = 1
    ! Check for   _cell_length_a
    CALL parser_search_string(parser,"_cell_length_a",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field (_cell_length_a) was not found in CIF file! "//&
CPSourceFileRef)
    CALL cif_get_real(parser,cell_lengths(1),error)
    cell_lengths(1) = cp_unit_to_cp2k(cell_lengths(1),"angstrom",error=error)

    ! Check for   _cell_length_b
    CALL parser_search_string(parser,"_cell_length_b",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field (_cell_length_b) was not found in CIF file! "//&
CPSourceFileRef)
    CALL cif_get_real(parser,cell_lengths(2),error)
    cell_lengths(2) = cp_unit_to_cp2k(cell_lengths(2),"angstrom",error=error)

    ! Check for   _cell_length_c
    CALL parser_search_string(parser,"_cell_length_c",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field (_cell_length_c) was not found in CIF file! "//&
CPSourceFileRef)
    CALL cif_get_real(parser,cell_lengths(3),error)
    cell_lengths(3) = cp_unit_to_cp2k(cell_lengths(3),"angstrom",error=error)

    ! Check for   _cell_angle_alpha
    CALL parser_search_string(parser,"_cell_angle_alpha",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field (_cell_angle_alpha) was not found in CIF file! "//&
CPSourceFileRef)
    CALL cif_get_real(parser,cell_angles(1),error)
    cell_angles(1) = cp_unit_to_cp2k(cell_angles(1),"deg",error=error)

    ! Check for   _cell_angle_beta
    CALL parser_search_string(parser,"_cell_angle_beta",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field (_cell_angle_beta) was not found in CIF file! "//&
CPSourceFileRef)
    CALL cif_get_real(parser,cell_angles(2),error)
    cell_angles(2) = cp_unit_to_cp2k(cell_angles(2),"deg",error=error)

    ! Check for   _cell_angle_gamma
    CALL parser_search_string(parser,"_cell_angle_gamma",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field (_cell_angle_gamma) was not found in CIF file! "//&
CPSourceFileRef)
    CALL cif_get_real(parser,cell_angles(3),error)
    cell_angles(3) = cp_unit_to_cp2k(cell_angles(3),"deg",error=error)

    ! Create cell
    NULLIFY (cell)
    CALL cell_create(cell,error=error)
    CALL set_cell_param(cell,cell_lengths,cell_angles,periodic=periodic,&
                        do_init_cell=.TRUE.,error=error)
    CALL write_cell(cell,subsys_section,label="CELL_CIF|",error=error)

    ! Parse atoms info and fractional coordinates
    ! Check for   _atom_site_label
    CALL parser_search_string(parser,"_atom_site_label",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    CALL cp_assert(found,cp_fatal_level,cp_assertion_failed,routineP,&
         "The field _atom_site_label was not found in CIF file! "//&
CPSourceFileRef)
    ifield     = 0
    ii         = 1
    ifield(ii) = 0
    CALL parser_get_next_line(parser,1,error=error)
    DO WHILE (INDEX(parser%input_line,"_atom_site_")/=0)
       ii          = ii + 1
       CPPostcondition(ii<=20,cp_failure_level,routineP,error,failure)

       ifield (ii) = -1
       ! Check for   _atom_site_type_symbol
       IF (INDEX(parser%input_line,"_atom_site_type_symbol")/=0) ifield(ii) = 0
       ! Check for   _atom_site_fract_x
       IF (INDEX(parser%input_line,"_atom_site_fract_x")/=0) ifield(ii) = 1
       ! Check for   _atom_site_fract_y
       IF (INDEX(parser%input_line,"_atom_site_fract_y")/=0) ifield(ii) = 2
       ! Check for   _atom_site_fract_z
       IF (INDEX(parser%input_line,"_atom_site_fract_z")/=0) ifield(ii) = 3
       CALL parser_get_next_line(parser,1,error=error)
    END DO
    ! Parse real info
    natom = 0
    DO WHILE ((INDEX(parser%input_line,"loop_")==0).AND.(parser%input_line(1:1)/="_"))
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
       ii    = 1
       itype = 0
       DO WHILE (ANY(ifield(ii:)>0))
          SELECT CASE (ifield(ii))
          CASE (-1)
             ! Skip this field
             CALL parser_get_object (parser,s_tag,error=error)
          CASE (0)
             CALL parser_get_object  (parser,strtmp,error=error)
             atom_info%id_atmname(natom) = str2id(strtmp)
             atom_info%id_molname(natom) = str2id(s2s("MOL"//TRIM(ADJUSTL(cp_to_string(natom)))))
             atom_info%id_resname(natom) = atom_info%id_molname(natom)
             atom_info%resid(natom)      = 1
             atom_info%id_element(natom) = atom_info%id_atmname(natom)
          CASE (1)
             CALL cif_get_real(parser,atom_info%r(1,natom),error)
          CASE (2)
             CALL cif_get_real(parser,atom_info%r(2,natom),error)
          CASE (3)
             CALL cif_get_real(parser,atom_info%r(3,natom),error)
          CASE DEFAULT
             ! Never reach this point..
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
          ii = ii + 1
       END DO
       s = atom_info%r(1:3,natom)
       CALL scaled_to_real(atom_info%r(1:3,natom),s,cell)
       CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
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
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
       END DO
    END DO
    ! Parse Symmetry Group and generation elements..
    ! Check for   _symmetry_space_group_name_h-m
    CALL parser_search_string(parser,"_symmetry_space_group_name_h-m",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    IF (found) THEN
       IF (iw>0) WRITE(iw,'(A)')" CIF_INFO| _symmetry_space_group_name_h-m :: "//TRIM(parser%input_line(parser%icol:))
    END IF

    ! Check for   _symmetry_equiv_pos_as_xyz
    ! Check for   _space_group_symop_operation_xyz
    CALL parser_search_string(parser,"_symmetry_equiv_pos_as_xyz",ignore_case=.FALSE.,found=found,&
         begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    IF (.NOT. found) THEN
       CALL parser_search_string(parser,"_space_group_symop_operation_xyz",ignore_case=.FALSE.,found=found,&
            begin_line=.FALSE.,search_from_begin_of_file=.TRUE.,error=error)
    END IF
    CALL cp_assert(found,cp_warning_level,cp_assertion_failed,routineP,&
         "The fields (_symmetry_equiv_pos_as_xyz) or (_space_group_symop_operation_xyz) were not found in CIF file! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)
    IF (iw>0) WRITE(iw,'(A,I0)')" CIF_INFO| Number of atoms before applying symmetry operations :: ",natom
    IF (iw>0) WRITE(iw,'(A10,1X,3F12.6)')(TRIM(id2str(atom_info%id_atmname(ii))),atom_info%r(1:3,ii),ii=1,natom)
    isym = 0
    IF (found) THEN
       ! Apply symmetry elements and generate the whole set of atoms in the unit cell
       CALL parser_get_next_line(parser,1,error=error)
       isym       = 0
       natom_orig = natom
       DO WHILE ((INDEX(parser%input_line,"loop_")==0).AND.(parser%input_line(1:1)/="_"))
          isym = isym + 1
          iln0 = INDEX(parser%input_line(1:)     ,"'")
          iln1 = INDEX(parser%input_line(iln0+1:)     ,",")+iln0
          iln2 = INDEX(parser%input_line(iln1+1:),",")+iln1
          IF (iln0==0) THEN
            iln3 = LEN_TRIM(parser%input_line)+1
          ELSE
            iln3 = INDEX(parser%input_line(iln2+1:),"'")+iln2
          END IF
          CPPostcondition(iln1/=0,cp_failure_level,routineP,error,failure)
          CPPostcondition(iln2/=iln1,cp_failure_level,routineP,error,failure)
          CPPostcondition(iln3/=iln2,cp_failure_level,routineP,error,failure)
          CALL initf(3)
          CALL parsef(1,TRIM(parser%input_line(iln0+1:iln1-1)),s2a("x","y","z"))
          CALL parsef(2,TRIM(parser%input_line(iln1+1:iln2-1)),s2a("x","y","z"))
          CALL parsef(3,TRIM(parser%input_line(iln2+1:iln3-1)),s2a("x","y","z"))
          Loop_over_unique_atoms: DO ii = 1, natom_orig
             CALL real_to_scaled(s_tmp, atom_info%r(1:3,ii), cell)
             s(1) = evalf(1,(/s_tmp(1),s_tmp(2),s_tmp(3)/))
             s(2) = evalf(2,(/s_tmp(1),s_tmp(2),s_tmp(3)/))
             s(3) = evalf(3,(/s_tmp(1),s_tmp(2),s_tmp(3)/))
             CALL scaled_to_real(r1, s, cell)
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
                atom_info%id_element(natom) = atom_info%id_element(ii)
                atom_info%resid(natom)      = atom_info%resid(ii)
                atom_info%r(1:3,natom)         = r1
             END IF
          END DO Loop_over_unique_atoms
          CALL finalizef()
          CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
          IF (my_end) EXIT
       END DO
    END IF
    IF (iw>0) WRITE(iw,'(A,I0)')" CIF_INFO| Number of symmetry operations :: ",isym
    IF (iw>0) WRITE(iw,'(A,I0)')" CIF_INFO| Number of total atoms :: ",natom
    IF (iw>0) WRITE(iw,'(A10,1X,3F12.6)')(TRIM(id2str(atom_info%id_atmname(ii))),atom_info%r(1:3,ii),ii=1,natom)

    ! Releasse local cell type and parser
    CALL cell_release(cell,error=error)
    CALL parser_release(parser,error=error)

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
         "PRINT%TOPOLOGY_INFO/CIF_INFO",error=error)
    CALL timestop(handle)
  END SUBROUTINE read_coordinate_cif

! *****************************************************************************
!> \brief  Reads REAL from the CIF file.. This wrapper is needed in order to
!>         treat properly the accuracy specified in the CIF file, i.e. 3.45(6)
!> \param parser ...
!> \param r ...
!> \param error ...
!> \date   12.2008
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE cif_get_real(parser, r, error)
    TYPE(cp_parser_type), POINTER            :: parser
    REAL(KIND=dp), INTENT(OUT)               :: r
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cif_get_real', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: s_tag
    INTEGER                                  :: iln

    CALL parser_get_object (parser,s_tag,error=error)
    iln = LEN_TRIM(s_tag)
    IF (INDEX(s_tag,"(")/=0) iln = INDEX(s_tag,"(")-1
    READ(s_tag(1:iln),*) r
  END SUBROUTINE cif_get_real

END MODULE topology_cif
