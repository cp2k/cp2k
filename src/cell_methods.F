!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Handles all functions related to the CELL
!> \par History
!>      11.2008 Teodoro Laino [tlaino] - deeply cleaning cell_type from units
!>      10.2014 Moved many routines to cell_types.F.
!> \author Matthias KracK (16.01.2002, based on a earlier version of CJM, JGH)
! **************************************************************************************************
MODULE cell_methods
   USE cell_types,                      ONLY: &
        cell_clone, cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, &
        cell_sym_hexagonal_gamma_60, cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, &
        cell_sym_none, cell_sym_orthorhombic, cell_sym_rhombohedral, cell_sym_tetragonal_ab, &
        cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type, get_cell, &
        plane_distance, use_perd_none, use_perd_x, use_perd_xy, use_perd_xyz, use_perd_xz, &
        use_perd_y, use_perd_yz, use_perd_z
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object,&
                                              parser_search_string
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                              cp_unit_to_cp2k
   USE input_constants,                 ONLY: do_cell_cif,&
                                              do_cell_cp2k,&
                                              do_cell_extxyz,&
                                              do_cell_pdb,&
                                              do_cell_xsc,&
                                              do_coord_cif,&
                                              do_coord_cp2k,&
                                              do_coord_pdb,&
                                              do_coord_xyz
   USE input_cp2k_subsys,               ONLY: create_cell_section
   USE input_enumeration_types,         ONLY: enum_i2c,&
                                              enumeration_type
   USE input_keyword_types,             ONLY: keyword_get,&
                                              keyword_type
   USE input_section_types,             ONLY: &
        section_get_keyword, section_release, section_type, section_vals_get, &
        section_vals_get_subs_vals, section_vals_type, section_vals_val_get, section_vals_val_set, &
        section_vals_val_unset
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: default_output_unit
   USE mathconstants,                   ONLY: degree,&
                                              sqrt3
   USE mathlib,                         ONLY: angle,&
                                              det_3x3,&
                                              inv_3x3
   USE message_passing,                 ONLY: mp_para_env_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_methods'

   PUBLIC :: cell_create, &
             init_cell, &
             read_cell, &
             read_cell_cif, &
             read_cell_extxyz, &
             read_cell_pdb, &
             set_cell_param, &
             write_cell, &
             write_cell_low

CONTAINS

! **************************************************************************************************
!> \brief allocates and initializes a cell
!> \param cell the cell to initialize
!> \param hmat the h matrix that defines the cell
!> \param periodic periodicity of the cell
!> \param tag ...
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cell_create(cell, hmat, periodic, tag)

      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: periodic
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: tag

      CPASSERT(.NOT. ASSOCIATED(cell))
      ALLOCATE (cell)
      cell%ref_count = 1
      IF (PRESENT(periodic)) THEN
         cell%perd = periodic
      ELSE
         cell%perd = 1
      END IF
      cell%orthorhombic = .FALSE.
      cell%symmetry_id = cell_sym_none
      IF (PRESENT(hmat)) CALL init_cell(cell, hmat)
      IF (PRESENT(tag)) cell%tag = tag

   END SUBROUTINE cell_create

! **************************************************************************************************
!> \brief   Initialise/readjust a simulation cell after hmat has been changed
!> \param cell ...
!> \param hmat ...
!> \param periodic ...
!> \date    16.01.2002
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_cell(cell, hmat, periodic)

      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: periodic

      REAL(KIND=dp), PARAMETER                           :: eps_hmat = 1.0E-14_dp

      INTEGER                                            :: dim
      REAL(KIND=dp)                                      :: a, acosa, acosah, acosg, alpha, asina, &
                                                            asinah, asing, beta, gamma, norm, &
                                                            norm_c
      REAL(KIND=dp), DIMENSION(3)                        :: abc

      CPASSERT(ASSOCIATED(cell))

      IF (PRESENT(hmat)) cell%hmat(:, :) = hmat(:, :)
      IF (PRESENT(periodic)) cell%perd(:) = periodic(:)

      cell%deth = ABS(det_3x3(cell%hmat))

      IF (cell%deth < 1.0E-10_dp) THEN
         CALL write_cell_low(cell, "angstrom", default_output_unit)
         CALL cp_abort(__LOCATION__, &
                       "An invalid set of cell vectors was specified. "// &
                       "The cell volume is too small")
      END IF

      SELECT CASE (cell%symmetry_id)
      CASE (cell_sym_cubic, &
            cell_sym_tetragonal_ab, &
            cell_sym_tetragonal_ac, &
            cell_sym_tetragonal_bc, &
            cell_sym_orthorhombic)
         CALL get_cell(cell=cell, abc=abc)
         abc(2) = plane_distance(0, 1, 0, cell=cell)
         abc(3) = plane_distance(0, 0, 1, cell=cell)
         SELECT CASE (cell%symmetry_id)
         CASE (cell_sym_cubic)
            abc(1:3) = SUM(abc(1:3))/3.0_dp
         CASE (cell_sym_tetragonal_ab, &
               cell_sym_tetragonal_ac, &
               cell_sym_tetragonal_bc)
            SELECT CASE (cell%symmetry_id)
            CASE (cell_sym_tetragonal_ab)
               a = 0.5_dp*(abc(1) + abc(2))
               abc(1) = a
               abc(2) = a
            CASE (cell_sym_tetragonal_ac)
               a = 0.5_dp*(abc(1) + abc(3))
               abc(1) = a
               abc(3) = a
            CASE (cell_sym_tetragonal_bc)
               a = 0.5_dp*(abc(2) + abc(3))
               abc(2) = a
               abc(3) = a
            END SELECT
         END SELECT
         cell%hmat(1, 1) = abc(1); cell%hmat(1, 2) = 0.0_dp; cell%hmat(1, 3) = 0.0_dp
         cell%hmat(2, 1) = 0.0_dp; cell%hmat(2, 2) = abc(2); cell%hmat(2, 3) = 0.0_dp
         cell%hmat(3, 1) = 0.0_dp; cell%hmat(3, 2) = 0.0_dp; cell%hmat(3, 3) = abc(3)
      CASE (cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120)
         CALL get_cell(cell=cell, abc=abc)
         a = 0.5_dp*(abc(1) + abc(2))
         acosg = 0.5_dp*a
         asing = sqrt3*acosg
         IF (cell%symmetry_id == cell_sym_hexagonal_gamma_120) acosg = -acosg
         cell%hmat(1, 1) = a; cell%hmat(1, 2) = acosg; cell%hmat(1, 3) = 0.0_dp
         cell%hmat(2, 1) = 0.0_dp; cell%hmat(2, 2) = asing; cell%hmat(2, 3) = 0.0_dp
         cell%hmat(3, 1) = 0.0_dp; cell%hmat(3, 2) = 0.0_dp; cell%hmat(3, 3) = abc(3)
      CASE (cell_sym_rhombohedral)
         CALL get_cell(cell=cell, abc=abc)
         a = SUM(abc(1:3))/3.0_dp
         alpha = (angle(cell%hmat(:, 3), cell%hmat(:, 2)) + &
                  angle(cell%hmat(:, 1), cell%hmat(:, 3)) + &
                  angle(cell%hmat(:, 1), cell%hmat(:, 2)))/3.0_dp
         acosa = a*COS(alpha)
         asina = a*SIN(alpha)
         acosah = a*COS(0.5_dp*alpha)
         asinah = a*SIN(0.5_dp*alpha)
         norm = acosa/acosah
         norm_c = SQRT(1.0_dp - norm*norm)
         cell%hmat(1, 1) = a; cell%hmat(1, 2) = acosa; cell%hmat(1, 3) = acosah*norm
         cell%hmat(2, 1) = 0.0_dp; cell%hmat(2, 2) = asina; cell%hmat(2, 3) = asinah*norm
         cell%hmat(3, 1) = 0.0_dp; cell%hmat(3, 2) = 0.0_dp; cell%hmat(3, 3) = a*norm_c
      CASE (cell_sym_monoclinic)
         CALL get_cell(cell=cell, abc=abc)
         beta = angle(cell%hmat(:, 1), cell%hmat(:, 3))
         cell%hmat(1, 1) = abc(1); cell%hmat(1, 2) = 0.0_dp; cell%hmat(1, 3) = abc(3)*COS(beta)
         cell%hmat(2, 1) = 0.0_dp; cell%hmat(2, 2) = abc(2); cell%hmat(2, 3) = 0.0_dp
         cell%hmat(3, 1) = 0.0_dp; cell%hmat(3, 2) = 0.0_dp; cell%hmat(3, 3) = abc(3)*SIN(beta)
      CASE (cell_sym_monoclinic_gamma_ab)
         ! Cell symmetry with a = b, alpha = beta = 90 degree and gammma not equal 90 degree
         CALL get_cell(cell=cell, abc=abc)
         a = 0.5_dp*(abc(1) + abc(2))
         gamma = angle(cell%hmat(:, 1), cell%hmat(:, 2))
         acosg = a*COS(gamma)
         asing = a*SIN(gamma)
         cell%hmat(1, 1) = a; cell%hmat(1, 2) = acosg; cell%hmat(1, 3) = 0.0_dp
         cell%hmat(2, 1) = 0.0_dp; cell%hmat(2, 2) = asing; cell%hmat(2, 3) = 0.0_dp
         cell%hmat(3, 1) = 0.0_dp; cell%hmat(3, 2) = 0.0_dp; cell%hmat(3, 3) = abc(3)
      CASE (cell_sym_triclinic)
         ! Nothing to do
      END SELECT

      ! Do we have an (almost) orthorhombic cell?
      IF ((ABS(cell%hmat(1, 2)) < eps_hmat) .AND. (ABS(cell%hmat(1, 3)) < eps_hmat) .AND. &
          (ABS(cell%hmat(2, 1)) < eps_hmat) .AND. (ABS(cell%hmat(2, 3)) < eps_hmat) .AND. &
          (ABS(cell%hmat(3, 1)) < eps_hmat) .AND. (ABS(cell%hmat(3, 2)) < eps_hmat)) THEN
         cell%orthorhombic = .TRUE.
      ELSE
         cell%orthorhombic = .FALSE.
      END IF

      ! Retain an exact orthorhombic cell
      ! (off-diagonal elements must remain zero identically to keep QS fast)
      IF (cell%orthorhombic) THEN
         cell%hmat(1, 2) = 0.0_dp
         cell%hmat(1, 3) = 0.0_dp
         cell%hmat(2, 1) = 0.0_dp
         cell%hmat(2, 3) = 0.0_dp
         cell%hmat(3, 1) = 0.0_dp
         cell%hmat(3, 2) = 0.0_dp
      END IF

      dim = COUNT(cell%perd == 1)
      IF ((dim == 1) .AND. (.NOT. cell%orthorhombic)) THEN
         CPABORT("Non-orthorhombic and not periodic")
      END IF

      ! Update deth and hmat_inv with enforced symmetry
      cell%deth = ABS(det_3x3(cell%hmat))
      IF (cell%deth < 1.0E-10_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "An invalid set of cell vectors was obtained after applying "// &
                       "the requested cell symmetry. The cell volume is too small")
      END IF
      cell%h_inv = inv_3x3(cell%hmat)

   END SUBROUTINE init_cell

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param cell_section ...
!> \param check_for_ref ...
!> \param para_env ...
!> \par History
!>      03.2005 created [teo]
!>      03.2026 revamped logic with pdb and extxyz parsers
!> \author Teodoro Laino
! **************************************************************************************************
   RECURSIVE SUBROUTINE read_cell(cell, cell_ref, use_ref_cell, cell_section, &
                                  topology_section, check_for_ref, para_env)

      TYPE(cell_type), POINTER                           :: cell, cell_ref
      LOGICAL, INTENT(INOUT), OPTIONAL                   :: use_ref_cell
      TYPE(section_vals_type), OPTIONAL, POINTER         :: cell_section, topology_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_for_ref
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=default_path_length)                 :: ang_unit,&
                                                            cell_file_name,&
                                                            coord_file_name,&
                                                            cryst,&
                                                            len_unit
      INTEGER                                            :: cell_file_format,&
                                                            coord_file_format,&
                                                            i, id1, id2, idum,&
                                                            j, my_per
      INTEGER, DIMENSION(:), POINTER                     :: multiple_unit_cell
      LOGICAL :: cell_read_a, cell_read_abc, cell_read_alpha_beta_gamma, cell_read_b, cell_read_c, &
         cell_read_file, cif_cell, my_check_ref, my_end, pdb_cryst, tmp_comb_abc, tmp_comb_cell, &
         tmp_comb_top, topo_read_coord
      REAL(KIND=dp)                                      :: xdum
      REAL(KIND=dp), DIMENSION(3, 3)                     :: read_mat
      REAL(KIND=dp), DIMENSION(3)                        :: read_ang, read_len
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cell_par
      REAL(KIND=dp), PARAMETER                           :: eps = EPSILON(0.0_dp)
      TYPE(cp_parser_type)                               :: parser
      TYPE(section_vals_type), POINTER                   :: cell_ref_section

      my_check_ref = .TRUE.
      NULLIFY (cell_ref_section, cell_par, multiple_unit_cell)
      IF (.NOT. ASSOCIATED(cell)) CALL cell_create(cell, tag="CELL")
      IF (.NOT. ASSOCIATED(cell_ref)) CALL cell_create(cell_ref, tag="CELL_REF")
      IF (PRESENT(check_for_ref)) my_check_ref = check_for_ref

      cell%deth = 0.0_dp
      cell%orthorhombic = .FALSE.
      cell%perd(:) = 1
      cell%symmetry_id = cell_sym_none
      cell%hmat(:, :) = 0.0_dp
      cell%h_inv(:, :) = 0.0_dp
      cell_read_file = .FALSE.
      cell_read_a = .FALSE.
      cell_read_b = .FALSE.
      cell_read_c = .FALSE.
      cell_read_abc = .FALSE.
      cell_read_alpha_beta_gamma = .FALSE.
      read_mat(:, :) = 0.0_dp
      read_ang(:) = 0.0_dp
      read_len(:) = 0.0_dp

      ! Precedence of retrieving cell information from input:
      ! 1. CELL/CELL_FILE_NAME
      ! 2. CELL/ABC and optionally CELL/ALPHA_BETA_GAMMA
      ! 3. CELL/A, CELL/B, CELL/C
      ! 4. TOPOLOGY/COORD_FILE_NAME, if topology_section is present
      ! The actual order of processing is 4 -> 1 -> 2 -> 3, with
      ! case 4 merged to case 1 (if file format permits) first.
      ! Store data into either read_mat or read_ang and read_len
      ! in CP2K units, which will be converted to cell%hamt and A, B, C.
      CALL section_vals_val_get(cell_section, "A", explicit=cell_read_a)
      CALL section_vals_val_get(cell_section, "B", explicit=cell_read_b)
      CALL section_vals_val_get(cell_section, "C", explicit=cell_read_c)
      CALL section_vals_val_get(cell_section, "ABC", explicit=cell_read_abc)
      CALL section_vals_val_get(cell_section, "ALPHA_BETA_GAMMA", explicit=cell_read_alpha_beta_gamma)
      CALL section_vals_val_get(cell_section, "CELL_FILE_NAME", explicit=cell_read_file)

      ! Case 4
      tmp_comb_top = (.NOT. (cell_read_file .OR. cell_read_abc))
      tmp_comb_top = (tmp_comb_top .AND. (.NOT. cell_read_a))
      tmp_comb_top = (tmp_comb_top .AND. (.NOT. cell_read_b))
      tmp_comb_top = (tmp_comb_top .AND. (.NOT. cell_read_c))
      IF (tmp_comb_top) THEN
         CALL cp_warn(__LOCATION__, &
                      "None of the keywords CELL_FILE_NAME, ABC, or A, B, C "// &
                      "are specified in CELL section. Will now attempt to read "// &
                      "TOPOLOGY/COORD_FILE_NAME if its format can be parsed for "// &
                      "cell information.")
         IF (ASSOCIATED(topology_section)) THEN
            CALL section_vals_val_get(topology_section, "COORD_FILE_NAME", explicit=topo_read_coord)
            IF (topo_read_coord) THEN
               CALL section_vals_val_get(topology_section, "COORD_FILE_NAME", c_val=coord_file_name)
               CALL section_vals_val_get(topology_section, "COORD_FILE_FORMAT", i_val=coord_file_format)
               SELECT CASE (coord_file_format) ! Add formats with both cell and coord parser manually
               CASE (do_coord_cif)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_NAME", c_val=coord_file_name)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_FORMAT", i_val=do_cell_cif)
               CASE (do_coord_cp2k)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_NAME", c_val=coord_file_name)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_FORMAT", i_val=do_cell_cp2k)
               CASE (do_coord_pdb)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_NAME", c_val=coord_file_name)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_FORMAT", i_val=do_cell_pdb)
               CASE (do_coord_xyz)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_NAME", c_val=coord_file_name)
                  CALL section_vals_val_set(cell_section, "CELL_FILE_FORMAT", i_val=do_cell_extxyz)
               CASE DEFAULT
                  CALL cp_abort(__LOCATION__, &
                               "COORD_FILE_FORMAT is not set to one of implemented CELL_FILE_FORMAT "// &
                               "and cannot be parsed for cell information!" )
               END SELECT
            ELSE
               CALL cp_abort(__LOCATION__, &
                            "COORD_FILE_NAME is not set, so no cell information available!")
            END IF
         ELSE
            CALL cp_warn(__LOCATION__, &
                         "TOPOLOGY section is not available, so COORD_FILE_NAME cannot "// &
                         "be parsed for cell information in lieu of missing CELL settings.")
         END IF
      END IF
      ! Former logic in SUBROUTINE read_cell_from_external_file is moved here
      CALL section_vals_val_get(cell_section, "CELL_FILE_NAME", explicit=cell_read_file)
      IF (cell_read_file) THEN ! Case 1
         len_unit = "angstrom" ! Initialize default units; when adding new cell format
         ang_unit = "deg"      ! parser, units should be set for specific case
         tmp_comb_cell = (cell_read_abc .OR. (cell_read_a .OR. (cell_read_b .OR. cell_read_c)))
         IF (tmp_comb_cell) &
            CALL cp_warn(__LOCATION__, &
                        "Cell Information provided through A, B, C, or ABC in conjunction "// &
                        "with CELL_FILE_NAME. The definition in external file will override "// &
                        "other ones.")
         CALL section_vals_val_get(cell_section, "CELL_FILE_NAME", c_val=cell_file_name)
         CALL section_vals_val_get(cell_section, "CELL_FILE_FORMAT", i_val=cell_file_format)
         SELECT CASE (cell_file_format)
         CASE (do_cell_cp2k)
            len_unit = "angstrom"
            CALL parser_create(parser, cell_file_name, para_env=para_env)
            CALL parser_get_next_line(parser, 1)
            my_end = .FALSE.
            DO WHILE (.NOT. my_end)
               READ (parser%input_line, *) idum, xdum, read_mat(:, 1), read_mat(:, 2), read_mat(:, 3)
               CALL parser_get_next_line(parser, 1, at_end=my_end)
            END DO
            CALL parser_release(parser)
         CASE (do_cell_xsc)
            len_unit = "angstrom"
            CALL parser_create(parser, cell_file_name, para_env=para_env)
            CALL parser_get_next_line(parser, 1)
            READ (parser%input_line, *) idum, read_mat(:, 1), read_mat(:, 2), read_mat(:, 3)
            CALL parser_release(parser)
         CASE (do_cell_extxyz)
            len_unit = "angstrom"
            CALL parser_create(parser, cell_file_name, para_env=para_env)
            CALL parser_get_next_line(parser, 2)
            CALL uppercase(parser%input_line)
            id1 = INDEX(parser%input_line, "LATTICE=")
            IF (id1 > 0) THEN
               id2 = INDEX(parser%input_line(id1 + 9:), '"')
               READ (parser%input_line(id1 + 9:id1 + id2 + 7), *) read_mat(:, 1), read_mat(:, 2), read_mat(:, 3)
               CALL parser_release(parser)
            ELSE
               CALL parser_release(parser)
               CALL cp_abort(__LOCATION__, &
                            "The field lattice= was not found on comment line "// &
                            "of XYZ file, so cell information cannot be set via "// &
                            "extended XYZ specification!")
            END IF
         CASE (do_cell_pdb)
            len_unit = "angstrom"
            ang_unit = "deg"
            CALL parser_create(parser, cell_file_name, para_env=para_env)
            CALL parser_search_string(parser, "CRYST1", ignore_case=.FALSE., found=pdb_cryst, &
                                begin_line=.TRUE., search_from_begin_of_file=.TRUE.)
            IF (pdb_cryst) THEN
               READ (parser%input_line, *) cryst, read_len(:), read_ang(:)
               CALL parser_release(parser)
            ELSE
               CALL parser_release(parser)
               CALL cp_abort(__LOCATION__, &
                            "The field CRYST1 was not found in PDB file, so cell "// &
                            "information cannot be set!")
            END IF
         CASE (do_cell_cif)
            len_unit = "angstrom"
            ang_unit = "deg"
            CALL parser_create(parser, cell_file_name, para_env=para_env)
            cif_cell = .FALSE.
            ! _cell_length_a
            CALL parser_search_string(parser, "_cell_length_a", ignore_case=.FALSE., &
                                      found=cif_cell, begin_line=.FALSE., &
                                      search_from_begin_of_file=.TRUE.)
            IF (.NOT. cif_cell) CALL cp_abort(__LOCATION__, "The field "// &
                                                            "_cell_length_a was not found in CIF file!")
            CALL cif_get_real(parser, read_len(1))
            ! _cell_length_b
            CALL parser_search_string(parser, "_cell_length_b", ignore_case=.FALSE., &
                                      found=cif_cell, begin_line=.FALSE., &
                                      search_from_begin_of_file=.TRUE.)
            IF (.NOT. cif_cell) CALL cp_abort(__LOCATION__, "The field "// &
                                                            "_cell_length_b was not found in CIF file!")
            CALL cif_get_real(parser, read_len(2))
            ! _cell_length_c
            CALL parser_search_string(parser, "_cell_length_c", ignore_case=.FALSE., &
                                      found=cif_cell, begin_line=.FALSE., &
                                      search_from_begin_of_file=.TRUE.)
            IF (.NOT. cif_cell) CALL cp_abort(__LOCATION__, "The field "// &
                                                            "_cell_length_c was not found in CIF file!")
            CALL cif_get_real(parser, read_len(3))
            ! _cell_angle_alpha
            CALL parser_search_string(parser, "_cell_angle_alpha", ignore_case=.FALSE., &
                                      found=cif_cell, begin_line=.FALSE., &
                                      search_from_begin_of_file=.TRUE.)
            IF (.NOT. cif_cell) CALL cp_abort(__LOCATION__, "The field "// &
                                                            "_cell_angle_alpha was not found in CIF file!")
            CALL cif_get_real(parser, read_ang(1))
            ! _cell_angle_beta
            CALL parser_search_string(parser, "_cell_angle_beta", ignore_case=.FALSE., &
                                      found=cif_cell, begin_line=.FALSE., &
                                      search_from_begin_of_file=.TRUE.)
            IF (.NOT. cif_cell) CALL cp_abort(__LOCATION__, "The field "// &
                                                            "_cell_angle_beta was not found in CIF file!")
            CALL cif_get_real(parser, read_ang(2))
            ! _cell_angle_gamma
            CALL parser_search_string(parser, "_cell_angle_gamma", ignore_case=.FALSE., &
                                      found=cif_cell, begin_line=.FALSE., &
                                      search_from_begin_of_file=.TRUE.)
            IF (.NOT. cif_cell) CALL cp_abort(__LOCATION__, "The field "// &
                                                            "_cell_angle_gamma was not found in CIF file!")
            CALL cif_get_real(parser, read_ang(3))
            CALL parser_release(parser)
         CASE DEFAULT
            CALL cp_abort(__LOCATION__, &
                          "CELL_FILE_FORMAT is not set to one of implemented options "// &
                          "and cannot be parsed for cell information!" )
         END SELECT
         ! Unit conversion
         DO i = 1, 3
            read_len(i) = cp_unit_to_cp2k(read_len(i), TRIM(len_unit))
            read_ang(i) = cp_unit_to_cp2k(read_ang(i), TRIM(ang_unit))
            DO j = 1, 3
               read_mat(j, i) = cp_unit_to_cp2k(read_mat(j, i), TRIM(len_unit))
            END DO
         END DO
      ELSE
         IF (cell_read_abc) THEN ! Case 2
            CALL section_vals_val_get(cell_section, "ABC", r_vals=cell_par)
            read_len = cell_par
            CALL section_vals_val_get(cell_section, "ALPHA_BETA_GAMMA", r_vals=cell_par)
            read_ang = cell_par
            IF (cell_read_a .OR. cell_read_b .OR. cell_read_c) &
               CALL cp_warn(__LOCATION__, &
                            "Cell information provided through vectors A, B or C in conjunction with ABC. "// &
                            "The definition of the ABC keyword will override the one provided by A, B and C.")
         ELSE ! Case 3
            tmp_comb_abc = ((cell_read_a .EQV. cell_read_b) .AND. (cell_read_b .EQV. cell_read_c))
            IF (tmp_comb_abc) THEN
               CALL section_vals_val_get(cell_section, "A", r_vals=cell_par)
               read_mat(:, 1) = cell_par(:)
               CALL section_vals_val_get(cell_section, "B", r_vals=cell_par)
               read_mat(:, 2) = cell_par(:)
               CALL section_vals_val_get(cell_section, "C", r_vals=cell_par)
               read_mat(:, 3) = cell_par(:)
               IF (cell_read_alpha_beta_gamma) &
                  CALL cp_warn(__LOCATION__, &
                               "The keyword ALPHA_BETA_GAMMA is ignored because it was used without the "// &
                               "keyword ABC.")
            ELSE
               CALL cp_abort(__LOCATION__, &
                             "Neither of the keywords CELL_FILE_NAME or ABC are specified, "// &
                             "and cell vector settings in A, B, C are incomplete!")
            END IF
         END IF
      END IF
      ! Unset other cell options
      IF (cell_read_file) &
         CALL section_vals_val_unset(cell_section, "CELL_FILE_NAME")
         CALL section_vals_val_unset(cell_section, "CELL_FILE_FORMAT")
      IF (cell_read_abc) &
         CALL section_vals_val_unset(cell_section, "ABC")
      IF (cell_read_alpha_beta_gamma) &
      CALL section_vals_val_unset(cell_section, "ALPHA_BETA_GAMMA")
      IF (cell_read_a) &
         CALL section_vals_val_unset(cell_section, "A")
      IF (cell_read_b) &
         CALL section_vals_val_unset(cell_section, "B")
      IF (cell_read_c) &
         CALL section_vals_val_unset(cell_section, "C")
      ! Convert read_mat or read_len and read_ang to actual cell%hmat
      IF (ANY(read_mat(:, :) > eps)) THEN
         ! Make a warning before storing cell vectors that
         ! do not form a triangular matrix.
         IF ((ABS(read_mat(2, 1)) > eps) .OR. &
             (ABS(read_mat(3, 1)) > eps) .OR. &
             (ABS(read_mat(3, 2)) > eps)) &
             CALL cp_warn(__LOCATION__, &
                         "Cell vectors are read but cell matrix is not "// &
                         "a lower triangle, not conforming to the program "// &
                         "convention that A lies along the X-axis and "// &
                         "B is in the XY plane.")
         cell%hmat = read_mat
      ELSE
         IF (ANY(read_ang(:) > eps) .AND. ANY(read_len(:) > eps)) THEN
            CALL set_cell_param(cell, cell_length=read_len, cell_angle=read_ang, &
                                do_init_cell=.FALSE.)
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "No meaningful cell information is read from parser!")
         END IF
      END IF
      ! Modify CELL section to use A, B, C keywords for future reuse.
      ALLOCATE (cell_par(3))
      cell_par = cell%hmat(:, 1)
      CALL section_vals_val_set(cell_section, "A", r_vals_ptr=cell_par)
      ALLOCATE (cell_par(3))
      cell_par = cell%hmat(:, 2)
      CALL section_vals_val_set(cell_section, "B", r_vals_ptr=cell_par)
      ALLOCATE (cell_par(3))
      cell_par = cell%hmat(:, 3)
      CALL section_vals_val_set(cell_section, "C", r_vals_ptr=cell_par)

      ! Multiple unit cell
      CALL section_vals_val_get(cell_section, "MULTIPLE_UNIT_CELL", i_vals=multiple_unit_cell)
      IF (ANY(multiple_unit_cell /= 1)) CALL set_multiple_unit_cell(cell, multiple_unit_cell)

      CALL section_vals_val_get(cell_section, "PERIODIC", i_val=my_per)
      SELECT CASE (my_per)
      CASE (use_perd_x)
         cell%perd = [1, 0, 0]
      CASE (use_perd_y)
         cell%perd = [0, 1, 0]
      CASE (use_perd_z)
         cell%perd = [0, 0, 1]
      CASE (use_perd_xy)
         cell%perd = [1, 1, 0]
      CASE (use_perd_xz)
         cell%perd = [1, 0, 1]
      CASE (use_perd_yz)
         cell%perd = [0, 1, 1]
      CASE (use_perd_xyz)
         cell%perd = [1, 1, 1]
      CASE (use_perd_none)
         cell%perd = [0, 0, 0]
      CASE DEFAULT
         CPABORT("")
      END SELECT

      ! Load requested cell symmetry
      CALL section_vals_val_get(cell_section, "SYMMETRY", i_val=cell%symmetry_id)

      ! Initialize cell
      CALL init_cell(cell)

      IF (my_check_ref) THEN
         ! Recursive check for reference cell requested
         cell_ref_section => section_vals_get_subs_vals(cell_section, "CELL_REF")
         IF (parsed_cp2k_input(cell_ref_section, check_this_section=.TRUE.)) THEN
            IF (PRESENT(use_ref_cell)) use_ref_cell = .TRUE.
            CALL read_cell(cell_ref, cell_ref, use_ref_cell=use_ref_cell, &
                           cell_section=cell_ref_section, check_for_ref=.FALSE., &
                           para_env=para_env)
         ELSE
            CALL cell_clone(cell, cell_ref, tag="CELL_REF")
            IF (PRESENT(use_ref_cell)) use_ref_cell = .FALSE.
         END IF
      END IF

   END SUBROUTINE read_cell

! **************************************************************************************************
!> \brief utility function to ease the transition to the new input.
!>      returns true if the new input was parsed
!> \param input_file the parsed input file
!> \param check_this_section ...
!> \return ...
!> \author fawzi
! **************************************************************************************************
   FUNCTION parsed_cp2k_input(input_file, check_this_section) RESULT(res)

      TYPE(section_vals_type), POINTER                   :: input_file
      LOGICAL, INTENT(IN), OPTIONAL                      :: check_this_section
      LOGICAL                                            :: res

      LOGICAL                                            :: my_check
      TYPE(section_vals_type), POINTER                   :: glob_section

      my_check = .FALSE.
      IF (PRESENT(check_this_section)) my_check = check_this_section
      res = ASSOCIATED(input_file)
      IF (res) THEN
         CPASSERT(input_file%ref_count > 0)
         IF (.NOT. my_check) THEN
            glob_section => section_vals_get_subs_vals(input_file, "GLOBAL")
            CALL section_vals_get(glob_section, explicit=res)
         ELSE
            CALL section_vals_get(input_file, explicit=res)
         END IF
      END IF

   END FUNCTION parsed_cp2k_input

! **************************************************************************************************
!> \brief   Sets the cell using the internal parameters (a,b,c) (alpha,beta,gamma)
!>          using the convention: a parallel to the x axis, b in the x-y plane and
!>          and c univoquely determined; gamma is the angle between a and b; beta
!>          is the angle between c and a and alpha is the angle between c and b
!> \param cell ...
!> \param cell_length ...
!> \param cell_angle ...
!> \param periodic ...
!> \param do_init_cell ...
!> \date    03.2008
!> \author  Teodoro Laino
! **************************************************************************************************
   SUBROUTINE set_cell_param(cell, cell_length, cell_angle, periodic, do_init_cell)

      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: cell_length, cell_angle
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: periodic
      LOGICAL, INTENT(IN)                                :: do_init_cell

      REAL(KIND=dp), PARAMETER                           :: eps = EPSILON(0.0_dp)

      REAL(KIND=dp)                                      :: cos_alpha, cos_beta, cos_gamma, sin_gamma

      CPASSERT(ASSOCIATED(cell))
      CPASSERT(ALL(cell_angle /= 0.0_dp))

      cos_gamma = COS(cell_angle(3)); IF (ABS(cos_gamma) < eps) cos_gamma = 0.0_dp
      IF (ABS(ABS(cos_gamma) - 1.0_dp) < eps) cos_gamma = SIGN(1.0_dp, cos_gamma)
      sin_gamma = SIN(cell_angle(3)); IF (ABS(sin_gamma) < eps) sin_gamma = 0.0_dp
      IF (ABS(ABS(sin_gamma) - 1.0_dp) < eps) sin_gamma = SIGN(1.0_dp, sin_gamma)
      cos_beta = COS(cell_angle(2)); IF (ABS(cos_beta) < eps) cos_beta = 0.0_dp
      IF (ABS(ABS(cos_beta) - 1.0_dp) < eps) cos_beta = SIGN(1.0_dp, cos_beta)
      cos_alpha = COS(cell_angle(1)); IF (ABS(cos_alpha) < eps) cos_alpha = 0.0_dp
      IF (ABS(ABS(cos_alpha) - 1.0_dp) < eps) cos_alpha = SIGN(1.0_dp, cos_alpha)

      cell%hmat(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
      cell%hmat(:, 2) = [cos_gamma, sin_gamma, 0.0_dp]
      cell%hmat(:, 3) = [cos_beta, (cos_alpha - cos_gamma*cos_beta)/sin_gamma, 0.0_dp]
      cell%hmat(3, 3) = SQRT(1.0_dp - cell%hmat(1, 3)**2 - cell%hmat(2, 3)**2)

      cell%hmat(:, 1) = cell%hmat(:, 1)*cell_length(1)
      cell%hmat(:, 2) = cell%hmat(:, 2)*cell_length(2)
      cell%hmat(:, 3) = cell%hmat(:, 3)*cell_length(3)

      IF (do_init_cell) THEN
         IF (PRESENT(periodic)) THEN
            CALL init_cell(cell=cell, periodic=periodic)
         ELSE
            CALL init_cell(cell=cell)
         END IF
      END IF

   END SUBROUTINE set_cell_param

! **************************************************************************************************
!> \brief   Setup of the multiple unit_cell
!> \param cell ...
!> \param multiple_unit_cell ...
!> \date    05.2009
!> \author  Teodoro Laino [tlaino]
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE set_multiple_unit_cell(cell, multiple_unit_cell)

      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(:), POINTER                     :: multiple_unit_cell

      CPASSERT(ASSOCIATED(cell))

      ! Abort, if one of the value is set to zero
      IF (ANY(multiple_unit_cell <= 0)) &
         CALL cp_abort(__LOCATION__, &
                       "CELL%MULTIPLE_UNIT_CELL accepts only integer values larger than 0! "// &
                       "A value of 0 or negative is meaningless!")

      ! Scale abc according to user request
      cell%hmat(:, 1) = cell%hmat(:, 1)*multiple_unit_cell(1)
      cell%hmat(:, 2) = cell%hmat(:, 2)*multiple_unit_cell(2)
      cell%hmat(:, 3) = cell%hmat(:, 3)*multiple_unit_cell(3)

   END SUBROUTINE set_multiple_unit_cell

! **************************************************************************************************
!> \brief  Reads cell information from CIF file
!> \param cif_file_name ...
!> \param cell ...
!> \param para_env ...
!> \date   12.2008
!> \par    Format Information implemented:
!>            _cell_length_a
!>            _cell_length_b
!>            _cell_length_c
!>            _cell_angle_alpha
!>            _cell_angle_beta
!>            _cell_angle_gamma
!>
!> \author Teodoro Laino [tlaino]
!>         moved from topology_cif (1/2019 JHU)
! **************************************************************************************************
   SUBROUTINE read_cell_cif(cif_file_name, cell, para_env)

      CHARACTER(len=*)                                   :: cif_file_name
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_cell_cif'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(3)                        :: cell_angles, cell_lengths
      TYPE(cp_parser_type)                               :: parser

      CALL timeset(routineN, handle)

      CALL parser_create(parser, cif_file_name, &
                         para_env=para_env, apply_preprocessing=.FALSE.)

      ! Parsing cell infos
      periodic = 1
      ! Check for   _cell_length_a
      CALL parser_search_string(parser, "_cell_length_a", ignore_case=.FALSE., found=found, &
                                begin_line=.FALSE., search_from_begin_of_file=.TRUE.)
      IF (.NOT. found) &
         CPABORT("The field (_cell_length_a) was not found in CIF file! ")
      CALL cif_get_real(parser, cell_lengths(1))
      cell_lengths(1) = cp_unit_to_cp2k(cell_lengths(1), "angstrom")

      ! Check for   _cell_length_b
      CALL parser_search_string(parser, "_cell_length_b", ignore_case=.FALSE., found=found, &
                                begin_line=.FALSE., search_from_begin_of_file=.TRUE.)
      IF (.NOT. found) &
         CPABORT("The field (_cell_length_b) was not found in CIF file! ")
      CALL cif_get_real(parser, cell_lengths(2))
      cell_lengths(2) = cp_unit_to_cp2k(cell_lengths(2), "angstrom")

      ! Check for   _cell_length_c
      CALL parser_search_string(parser, "_cell_length_c", ignore_case=.FALSE., found=found, &
                                begin_line=.FALSE., search_from_begin_of_file=.TRUE.)
      IF (.NOT. found) &
         CPABORT("The field (_cell_length_c) was not found in CIF file! ")
      CALL cif_get_real(parser, cell_lengths(3))
      cell_lengths(3) = cp_unit_to_cp2k(cell_lengths(3), "angstrom")

      ! Check for   _cell_angle_alpha
      CALL parser_search_string(parser, "_cell_angle_alpha", ignore_case=.FALSE., found=found, &
                                begin_line=.FALSE., search_from_begin_of_file=.TRUE.)
      IF (.NOT. found) &
         CPABORT("The field (_cell_angle_alpha) was not found in CIF file! ")
      CALL cif_get_real(parser, cell_angles(1))
      cell_angles(1) = cp_unit_to_cp2k(cell_angles(1), "deg")

      ! Check for   _cell_angle_beta
      CALL parser_search_string(parser, "_cell_angle_beta", ignore_case=.FALSE., found=found, &
                                begin_line=.FALSE., search_from_begin_of_file=.TRUE.)
      IF (.NOT. found) &
         CPABORT("The field (_cell_angle_beta) was not found in CIF file! ")
      CALL cif_get_real(parser, cell_angles(2))
      cell_angles(2) = cp_unit_to_cp2k(cell_angles(2), "deg")

      ! Check for   _cell_angle_gamma
      CALL parser_search_string(parser, "_cell_angle_gamma", ignore_case=.FALSE., found=found, &
                                begin_line=.FALSE., search_from_begin_of_file=.TRUE.)
      IF (.NOT. found) &
         CPABORT("The field (_cell_angle_gamma) was not found in CIF file! ")
      CALL cif_get_real(parser, cell_angles(3))
      cell_angles(3) = cp_unit_to_cp2k(cell_angles(3), "deg")

      ! Create cell
      CALL set_cell_param(cell, cell_lengths, cell_angles, periodic=periodic, &
                          do_init_cell=.TRUE.)

      CALL parser_release(parser)

      CALL timestop(handle)

   END SUBROUTINE read_cell_cif

! **************************************************************************************************
!> \brief  Reads REAL from the CIF file.. This wrapper is needed in order to
!>         treat properly the accuracy specified in the CIF file, i.e. 3.45(6)
!> \param parser ...
!> \param r ...
!> \date   12.2008
!> \author Teodoro Laino [tlaino]
! **************************************************************************************************
   SUBROUTINE cif_get_real(parser, r)

      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      REAL(KIND=dp), INTENT(OUT)                         :: r

      CHARACTER(LEN=default_string_length)               :: s_tag
      INTEGER                                            :: iln

      CALL parser_get_object(parser, s_tag)
      iln = LEN_TRIM(s_tag)
      IF (INDEX(s_tag, "(") /= 0) iln = INDEX(s_tag, "(") - 1
      READ (s_tag(1:iln), *) r

   END SUBROUTINE cif_get_real

! **************************************************************************************************
!> \brief  Reads cell information from comment line of extended xyz file
!> \param extxyz_file_name ...
!> \param cell ...
!> \param para_env ...
!> \date   03.2026
!> \par    Extended xyz format has comment on the second line in the form as
!>               Lattice="Ax Ay Az Bx By Bz Cx Cy Cz"
!>         where Ax, Ay, Az are three Cartesian components of cell vector A,
!>         Bx, By, Bz are components of B, Cx, Cy, Cz are components of C,
!>         all in the unit of angstrom.
! **************************************************************************************************
   SUBROUTINE read_cell_extxyz(extxyz_file_name, cell, para_env)

      CHARACTER(len=*)                                   :: extxyz_file_name
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_cell_extxyz'

      INTEGER                                            :: handle, i, id1, id2, j
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cp_parser_type)                               :: parser

      CALL timeset(routineN, handle)

      CALL parser_create(parser, extxyz_file_name, &
                         para_env=para_env, apply_preprocessing=.FALSE.)
      CALL parser_get_next_line(parser, 2)
      CALL uppercase(parser%input_line)
      id1 = INDEX(parser%input_line, "LATTICE=")
      IF (id1 > 0) THEN
         id2 = INDEX(parser%input_line(id1 + 9:), '"')
         READ (parser%input_line(id1 + 9:id1 + id2 + 7), *) hmat(:, 1), hmat(:, 2), hmat(:, 3)
         DO i = 1, 3
            DO j = 1, 3
               cell%hmat(j, i) = cp_unit_to_cp2k(hmat(j, i), "angstrom")
            END DO
         END DO

      ELSE
         CALL cp_abort(__LOCATION__, &
                       "The field (lattice=) was not found on comment line "// &
                       "of XYZ file, so cell information cannot be set via "// &
                       "extended XYZ specification! ")
      END IF

      CALL parser_release(parser)

      CALL timestop(handle)

   END SUBROUTINE read_cell_extxyz

! **************************************************************************************************
!> \brief  Reads cell information from CRYST1 record of PDB file
!> \param pdb_file_name ...
!> \param cell ...
!> \param para_env ...
!> \date   03.2026
!> \par    CRYST1 record may contain space group and Z value at the end,
!>         but here only the first entries are read:
!>         COLUMNS       DATA  TYPE    FIELD          DEFINITION
!>         -------------------------------------------------------------
!>          1 -  6       Record name   "CRYST1"
!>          7 - 15       Real(9.3)     a              a (Angstroms).
!>         16 - 24       Real(9.3)     b              b (Angstroms).
!>         25 - 33       Real(9.3)     c              c (Angstroms).
!>         34 - 40       Real(7.2)     alpha          alpha (degrees).
!>         41 - 47       Real(7.2)     beta           beta (degrees).
!>         48 - 54       Real(7.2)     gamma          gamma (degrees).
! **************************************************************************************************
   SUBROUTINE read_cell_pdb(pdb_file_name, cell, para_env)

      CHARACTER(len=*)                                   :: pdb_file_name
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_cell_pdb'

      CHARACTER(LEN=default_string_length)               :: cryst
      INTEGER                                            :: handle, i
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(3)                        :: cell_angles, cell_lengths
      TYPE(cp_parser_type)                               :: parser

      CALL timeset(routineN, handle)

      CALL parser_create(parser, pdb_file_name, &
                         para_env=para_env, apply_preprocessing=.FALSE.)

      CALL parser_search_string(parser, "CRYST1", ignore_case=.FALSE., found=found, &
                                begin_line=.TRUE., search_from_begin_of_file=.TRUE.)
      IF (.NOT. found) &
         CPABORT("The field (CRYST1) was not found in PDB file! ")

      periodic = 1
      READ (parser%input_line, *) cryst, cell_lengths(:), cell_angles(:)
      DO i = 1, 3
         cell_lengths(i) = cp_unit_to_cp2k(cell_lengths(i), "angstrom")
         cell_angles(i) = cp_unit_to_cp2k(cell_angles(i), "deg")
      END DO
      CALL set_cell_param(cell, cell_lengths, cell_angles, periodic=periodic, &
                          do_init_cell=.TRUE.)

      CALL parser_release(parser)

      CALL timestop(handle)

   END SUBROUTINE read_cell_pdb

! **************************************************************************************************
!> \brief   Write the cell parameters to the output unit.
!> \param cell ...
!> \param subsys_section ...
!> \param tag ...
!> \date    02.06.2000
!> \par     History
!>    - 11.2008 Teodoro Laino [tlaino] - rewrite and enabling user driven units
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_cell(cell, subsys_section, tag)

      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: subsys_section
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: tag

      CHARACTER(LEN=default_string_length)               :: label, unit_str
      INTEGER                                            :: output_unit
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (PRESENT(tag)) THEN
         label = TRIM(tag)//"|"
      ELSE
         label = TRIM(cell%tag)//"|"
      END IF

      output_unit = cp_print_key_unit_nr(logger, subsys_section, "PRINT%CELL", extension=".Log")
      CALL section_vals_val_get(subsys_section, "PRINT%CELL%UNIT", c_val=unit_str)
      CALL write_cell_low(cell, unit_str, output_unit, label)
      CALL cp_print_key_finished_output(output_unit, logger, subsys_section, "PRINT%CELL")

   END SUBROUTINE write_cell

! **************************************************************************************************
!> \brief  Write the cell parameters to the output unit
!> \param cell ...
!> \param unit_str ...
!> \param output_unit ...
!> \param label ...
!> \date  17.05.2023
!> \par   History
!>        - Extracted from write_cell (17.05.2023, MK)
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE write_cell_low(cell, unit_str, output_unit, label)

      TYPE(cell_type), POINTER                           :: cell
      CHARACTER(LEN=*), INTENT(IN)                       :: unit_str
      INTEGER, INTENT(IN)                                :: output_unit
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: label

      CHARACTER(LEN=12)                                  :: tag
      CHARACTER(LEN=3)                                   :: string
      CHARACTER(LEN=default_string_length)               :: my_label
      REAL(KIND=dp)                                      :: alpha, beta, gamma, val
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(enumeration_type), POINTER                    :: enum
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: section

      NULLIFY (enum)
      NULLIFY (keyword)
      NULLIFY (section)

      IF (output_unit > 0) THEN
         CALL get_cell(cell=cell, abc=abc, alpha=alpha, beta=beta, gamma=gamma, tag=tag)
         IF (PRESENT(label)) THEN
            my_label = label
         ELSE
            my_label = TRIM(tag)//"|"
         END IF
         val = cp_unit_from_cp2k(cell%deth, TRIM(unit_str)//"^3")
         WRITE (UNIT=output_unit, FMT="(/,T2,A,T61,F20.6)") &
            TRIM(my_label)//" Volume ["//TRIM(unit_str)//"^3]:", val
         val = cp_unit_from_cp2k(1.0_dp, TRIM(unit_str))
         WRITE (UNIT=output_unit, FMT="(T2,A,T30,3F10.3,3X,A6,F12.6)") &
            TRIM(my_label)//" Vector a ["//TRIM(unit_str)//"]:", cell%hmat(:, 1)*val, &
            "|a| = ", abc(1)*val, &
            TRIM(my_label)//" Vector b ["//TRIM(unit_str)//"]:", cell%hmat(:, 2)*val, &
            "|b| = ", abc(2)*val, &
            TRIM(my_label)//" Vector c ["//TRIM(unit_str)//"]:", cell%hmat(:, 3)*val, &
            "|c| = ", abc(3)*val
         WRITE (UNIT=output_unit, FMT="(T2,A,T69,F12.6)") &
            TRIM(my_label)//" Angle (b,c), alpha [degree]: ", alpha, &
            TRIM(my_label)//" Angle (a,c), beta  [degree]: ", beta, &
            TRIM(my_label)//" Angle (a,b), gamma [degree]: ", gamma
         IF (cell%symmetry_id /= cell_sym_none) THEN
            CALL create_cell_section(section)
            keyword => section_get_keyword(section, "SYMMETRY")
            CALL keyword_get(keyword, enum=enum)
            WRITE (UNIT=output_unit, FMT="(T2,A,T61,A20)") &
               TRIM(my_label)//" Requested initial symmetry: ", &
               ADJUSTR(TRIM(enum_i2c(enum, cell%symmetry_id)))
            CALL section_release(section)
         END IF
         IF (cell%orthorhombic) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               TRIM(my_label)//" Numerically orthorhombic: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               TRIM(my_label)//" Numerically orthorhombic: ", " NO"
         END IF
         IF (SUM(cell%perd(1:3)) == 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T77,A4)") &
               TRIM(my_label)//" Periodicity", "NONE"
         ELSE
            string = ""
            IF (cell%perd(1) == 1) string = TRIM(string)//"X"
            IF (cell%perd(2) == 1) string = TRIM(string)//"Y"
            IF (cell%perd(3) == 1) string = TRIM(string)//"Z"
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               TRIM(my_label)//" Periodicity", ADJUSTR(string)
         END IF
      END IF

   END SUBROUTINE write_cell_low

END MODULE cell_methods
