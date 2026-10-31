!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Parse CELL_OPT constraints, build symmetry constraint rows, and
!>        compute the cell enthalpy gradient.
!> \par History
!>      03.2008 - Teodoro Laino [tlaino] - University of Zurich - Cell Optimization
! **************************************************************************************************
MODULE cell_opt_utils
   USE cell_manifold_utils,             ONLY: &
        comp_ax, comp_bx, comp_by, comp_cx, comp_cy, comp_cz, geom_a, geom_alpha, geom_b, &
        geom_beta, geom_c, geom_gamma, gfix_alpha, gfix_beta, gfix_gamma, mask_a, mask_all, &
        mask_b, mask_c, mask_x, mask_xy, mask_xz, mask_y, mask_yz, mask_z, min_cell_length, &
        n_geom, n_slot, slot_ax, slot_bx, slot_by, slot_cx, slot_cy, slot_cz
   USE cell_types,                      ONLY: &
        cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_none, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE input_constants,                 ONLY: fix_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi,&
                                              sqrt3
   USE space_groups_types,              ONLY: spgr_type
   USE string_utilities,                ONLY: uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_utils'

   PUBLIC :: cell_opt_constraint_label, &
             comp_ax, comp_bx, comp_by, comp_cx, comp_cy, comp_cz, &
             get_dg_dh, &
             gfix_alpha, gfix_beta, gfix_gamma, &
             mask_a, mask_all, mask_b, mask_c, mask_x, mask_xy, mask_xz, mask_y, mask_yz, mask_z, &
             parse_cell_constraint, &
             parse_token, &
             read_external_press_tensor, &
             spgr_point_group_symmetry_id, &
             symmetry_constraint_rows

CONTAINS

! **************************************************************************************************
!> \brief Reads the external pressure tensor
!> \param geo_section ...
!> \param cell ...
!> \param pres_ext ...
!> \param mtrx ...
!> \param input_to_opt rotation from the input Cartesian frame to the optimization frame
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! **************************************************************************************************
   SUBROUTINE read_external_press_tensor(geo_section, cell, pres_ext, mtrx, input_to_opt)
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(OUT)                         :: pres_ext
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: mtrx
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: input_to_opt

      REAL(KIND=dp), PARAMETER                           :: eps_aniso_press = 1.0E-5_dp

      INTEGER                                            :: i, ind, j
      LOGICAL                                            :: check
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pres_ext_tens
      REAL(KIND=dp), DIMENSION(:), POINTER               :: pvals

      NULLIFY (pvals)
      mtrx = 0.0_dp
      pres_ext_tens = 0.0_dp
      pres_ext = 0.0_dp
      CALL section_vals_val_get(geo_section, "EXTERNAL_PRESSURE", r_vals=pvals)
      check = (SIZE(pvals) == 1) .OR. (SIZE(pvals) == 9)
      IF (.NOT. check) THEN
         CPABORT("EXTERNAL_PRESSURE can have 1 or 9 components only!")
      END IF

      IF (SIZE(pvals) == 9) THEN
         ind = 0
         DO i = 1, 3
            DO j = 1, 3
               ind = ind + 1
               pres_ext_tens(j, i) = pvals(ind)
            END DO
         END DO
         ! Express the pressure tensor in the canonical optimization frame.
         pres_ext_tens = MATMUL(input_to_opt, &
                                MATMUL(pres_ext_tens, TRANSPOSE(input_to_opt)))
         DO i = 1, 3
            pres_ext = pres_ext + pres_ext_tens(i, i)
         END DO
         pres_ext = pres_ext/3.0_dp
         DO i = 1, 3
            pres_ext_tens(i, i) = pres_ext_tens(i, i) - pres_ext
         END DO
      ELSE
         pres_ext = pvals(1)
      END IF

      IF (ANY(ABS(pres_ext_tens) > eps_aniso_press)) THEN
         mtrx = cell%deth*MATMUL(cell%h_inv, MATMUL(pres_ext_tens, TRANSPOSE(cell%h_inv)))
      END IF

   END SUBROUTINE read_external_press_tensor

! **************************************************************************************************
!> \brief Parses the tokens of the CONSTRAINT keyword into a cell component mask
!> \param constraint_tokens value tokens of the CONSTRAINT keyword (char_t list)
!> \param cell simulation cell (only the periodicity is used)
!> \param explicit whether the CONSTRAINT keyword was given explicitly
!> \param mask union of the fixed cell components (comp_* and mask_* bits)
! **************************************************************************************************
   SUBROUTINE parse_cell_constraint(constraint_tokens, cell, explicit, mask)
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)         :: constraint_tokens
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: explicit
      INTEGER, INTENT(OUT)                               :: mask

      CHARACTER(LEN=default_string_length)               :: reject_reason, token
      INTEGER                                            :: axis_index, nonperiodic_axis, overlap, &
                                                            required_mask, tok_index, tok_mask
      LOGICAL                                            :: accepted

      ! The input tokenizer delivers a bare keyword as one blank token.
      IF (SIZE(constraint_tokens) == 0 .OR. &
          (SIZE(constraint_tokens) == 1 .AND. LEN_TRIM(constraint_tokens(1)) == 0)) THEN
         CALL cp_abort(__LOCATION__, &
                       "CELL_OPT CONSTRAINT requires a value.")
      END IF

      mask = 0
      DO tok_index = 1, SIZE(constraint_tokens)
         token = ADJUSTL(constraint_tokens(tok_index))
         CALL uppercase(token)
         IF (token == "NONE") THEN
            ! NONE must be the first and only token: leftover tokens on the
            ! line or an already-accumulated mask mean it is being combined.
            IF (tok_index < SIZE(constraint_tokens) .OR. mask /= 0) THEN
               CALL cp_abort(__LOCATION__, &
                             "CELL_OPT CONSTRAINT NONE cannot be combined "// &
                             "with other tokens.")
            END IF
            mask = 0
            CYCLE
         END IF
         CALL parse_token(token, accepted, tok_mask, reject_reason)
         IF (.NOT. accepted) THEN
            CALL cp_abort(__LOCATION__, TRIM(reject_reason))
         END IF
         overlap = IAND(tok_mask, mask)
         IF (overlap == tok_mask) THEN
            CALL cp_warn(__LOCATION__, &
                         "Redundant CELL_OPT CONSTRAINT token: "//TRIM(token)// &
                         " (already constrained: "// &
                         TRIM(ADJUSTL(cell_opt_constraint_label(overlap)))//").")
         ELSE IF (overlap /= 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "CELL_OPT CONSTRAINT token "//TRIM(token)//" partially overlaps: "// &
                         "already constrained "// &
                         TRIM(ADJUSTL(cell_opt_constraint_label(overlap)))// &
                         ", newly fixed "// &
                         TRIM(ADJUSTL(cell_opt_constraint_label(IAND(NOT(mask), tok_mask))))//".")
         END IF
         mask = IOR(mask, tok_mask)
      END DO

      ! For 2D periodicity the non-periodic cell direction must be constrained.
      ! Complete the constraint automatically if the keyword was not given.
      IF (COUNT(cell%perd /= 0) == 2) THEN
         nonperiodic_axis = 0
         DO axis_index = 1, 3
            IF (cell%perd(axis_index) == 0) nonperiodic_axis = axis_index
         END DO
         required_mask = 0
         SELECT CASE (nonperiodic_axis)
         CASE (1)
            required_mask = mask_x
         CASE (2)
            required_mask = mask_y
         CASE (3)
            required_mask = mask_z
         END SELECT
         IF (.NOT. explicit .AND. mask == 0) THEN
            mask = required_mask
            CALL cp_warn(__LOCATION__, &
                         "2D CELL_OPT: constraining non-periodic cell direction.")
         ELSE IF (IAND(mask, required_mask) /= required_mask) THEN
            CALL cp_abort(__LOCATION__, &
                          "2D CELL_OPT needs constrained non-periodic cell direction.")
         END IF
      END IF

   END SUBROUTINE parse_cell_constraint

! **************************************************************************************************
!> \brief Returns the label of the cell constraint mask
!> \param mask bit mask of fixed cell components (comp_*, mask_*)
!> \return 4-character legacy label or space-separated component names
! **************************************************************************************************
   PURE FUNCTION cell_opt_constraint_label(mask) RESULT(label)
      INTEGER, INTENT(IN)                                :: mask
      CHARACTER(LEN=default_string_length)               :: label

      SELECT CASE (mask)
      CASE (fix_none)
         label = "NONE"
      CASE (mask_x)
         label = "   X"
      CASE (mask_y)
         label = "   Y"
      CASE (mask_z)
         label = "   Z"
      CASE (mask_xy)
         label = "  XY"
      CASE (mask_xz)
         label = "  XZ"
      CASE (mask_yz)
         label = "  YZ"
      CASE DEFAULT
         ! Composite constraint: spell out the fixed components
         label = ""
         IF (IAND(mask, comp_ax) /= 0) label = TRIM(label)//" AX"
         IF (IAND(mask, comp_bx) /= 0) label = TRIM(label)//" BX"
         IF (IAND(mask, comp_by) /= 0) label = TRIM(label)//" BY"
         IF (IAND(mask, comp_cx) /= 0) label = TRIM(label)//" CX"
         IF (IAND(mask, comp_cy) /= 0) label = TRIM(label)//" CY"
         IF (IAND(mask, comp_cz) /= 0) label = TRIM(label)//" CZ"
         IF (IAND(mask, gfix_alpha) /= 0) label = TRIM(label)//" ALPHA"
         IF (IAND(mask, gfix_beta) /= 0) label = TRIM(label)//" BETA"
         IF (IAND(mask, gfix_gamma) /= 0) label = TRIM(label)//" GAMMA"
         label = ADJUSTL(label)
      END SELECT

   END FUNCTION cell_opt_constraint_label

! **************************************************************************************************
!> \brief Translate a crystal system into cell constraint rows.
!> \param symmetry_id crystal system of the canonical cell
!> \param mask comp_* bits for canonical entries fixed to zero by symmetry
!> \param lin_w linear coupling normals, one normal per column
!> \param n_lin_rows number of active columns in lin_w
!> \param eq_pair index pairs [m, n] in 1..6 for geometric equalities, one pair per column
!> \param n_eq_rows number of active columns in eq_pair
!> \note Rhombohedral (a = b = c, alpha = beta = gamma) and monoclinic_gamma_ab
!>       (a = b with free gamma) produce curved geometric equalities. All other
!>       crystal systems produce constant linear rows with zero curvature.
! **************************************************************************************************
   SUBROUTINE symmetry_constraint_rows(symmetry_id, mask, lin_w, n_lin_rows, eq_pair, n_eq_rows)
      INTEGER, INTENT(IN)                                :: symmetry_id
      INTEGER, INTENT(OUT)                               :: mask
      REAL(KIND=dp), DIMENSION(n_slot, n_slot), &
         INTENT(OUT)                                     :: lin_w
      INTEGER, INTENT(OUT)                               :: n_lin_rows
      INTEGER, DIMENSION(2, n_geom), INTENT(OUT)         :: eq_pair
      INTEGER, INTENT(OUT)                               :: n_eq_rows

      REAL(KIND=dp), PARAMETER                           :: half_sqrt3 = 0.5_dp*sqrt3

      mask = 0
      lin_w = 0.0_dp
      n_lin_rows = 0
      eq_pair = 0
      n_eq_rows = 0

      SELECT CASE (symmetry_id)
      CASE (cell_sym_orthorhombic)
         mask = comp_bx + comp_cx + comp_cy
      CASE (cell_sym_tetragonal_ab)
         mask = comp_bx + comp_cx + comp_cy
         n_lin_rows = 1
         lin_w(slot_ax, 1) = 1.0_dp
         lin_w(slot_by, 1) = -1.0_dp
      CASE (cell_sym_tetragonal_ac)
         mask = comp_bx + comp_cx + comp_cy
         n_lin_rows = 1
         lin_w(slot_ax, 1) = 1.0_dp
         lin_w(slot_cz, 1) = -1.0_dp
      CASE (cell_sym_tetragonal_bc)
         mask = comp_bx + comp_cx + comp_cy
         n_lin_rows = 1
         lin_w(slot_by, 1) = 1.0_dp
         lin_w(slot_cz, 1) = -1.0_dp
      CASE (cell_sym_cubic)
         mask = comp_bx + comp_cx + comp_cy
         n_lin_rows = 2
         lin_w(slot_ax, 1) = 1.0_dp
         lin_w(slot_by, 1) = -1.0_dp
         lin_w(slot_by, 2) = 1.0_dp
         lin_w(slot_cz, 2) = -1.0_dp
      CASE (cell_sym_hexagonal_gamma_60)
         mask = comp_cx + comp_cy
         n_lin_rows = 2
         lin_w(slot_ax, 1) = -0.5_dp
         lin_w(slot_bx, 1) = 1.0_dp
         lin_w(slot_ax, 2) = -half_sqrt3
         lin_w(slot_by, 2) = 1.0_dp
      CASE (cell_sym_hexagonal_gamma_120)
         mask = comp_cx + comp_cy
         n_lin_rows = 2
         lin_w(slot_ax, 1) = 0.5_dp
         lin_w(slot_bx, 1) = 1.0_dp
         lin_w(slot_ax, 2) = -half_sqrt3
         lin_w(slot_by, 2) = 1.0_dp
      CASE (cell_sym_rhombohedral)
         n_eq_rows = 4
         eq_pair(:, 1) = [geom_a, geom_b]
         eq_pair(:, 2) = [geom_b, geom_c]
         eq_pair(:, 3) = [geom_alpha, geom_beta]
         eq_pair(:, 4) = [geom_beta, geom_gamma]
      CASE (cell_sym_monoclinic)
         mask = comp_bx + comp_cy
      CASE (cell_sym_monoclinic_gamma_ab)
         mask = comp_cx + comp_cy
         n_eq_rows = 1
         eq_pair(:, 1) = [geom_a, geom_b]
      CASE (cell_sym_none, cell_sym_triclinic)
         CONTINUE
      CASE DEFAULT
         CPABORT("Unknown cell symmetry_id")
      END SELECT

   END SUBROUTINE symmetry_constraint_rows

! **************************************************************************************************
!> \brief Identify the canonical crystal system preserved by a space group's point group.
!> \param spgr space group with integer lattice rotations in rotations_subset
!> \param cell symmetric start cell in the canonical frame
!> \return symmetry_id of the invariant metric form, or cell_sym_triclinic when
!>         the setting adds no canonical cell constraints
!> \note Averaging the metric G = hmat^T hmat over the point group produces the
!>       invariant metric form {G : R G R^T = G for all R}. Classifying the lengths
!>       and angles of this averaged metric maps the point group to the row sets of
!>       symmetry_constraint_rows. Non-canonical settings (such as a monoclinic
!>       unique axis along a or c) return cell_sym_triclinic and add zero cell rows.
! **************************************************************************************************
   FUNCTION spgr_point_group_symmetry_id(spgr, cell) RESULT(symmetry_id)
      TYPE(spgr_type), INTENT(IN), POINTER               :: spgr
      TYPE(cell_type), POINTER                           :: cell
      INTEGER                                            :: symmetry_id

      REAL(KIND=dp), PARAMETER :: ang_120 = 2.0_dp*pi/3.0_dp, ang_60 = pi/3.0_dp, &
         ang_90 = 0.5_dp*pi, eq_angle_tol = 1.0E-6_dp, rel_len_tol = 1.0E-6_dp, &
         special_ang_tol = 1.0E-5_dp

      INTEGER                                            :: ir, n90, nop
      REAL(KIND=dp)                                      :: ang(3), gsym(3, 3), len_tol, lens(3), &
                                                            metric(3, 3), rot(3, 3)

      symmetry_id = cell_sym_triclinic
      nop = spgr%n_operations_subset
      IF (nop == 0) RETURN

      metric = MATMUL(TRANSPOSE(cell%hmat), cell%hmat)
      gsym = 0.0_dp
      DO ir = 1, nop
         ! rotations_subset stores R_frac^T (r_out = rot^T r_in in spgr_apply_rotations_coord),
         ! so R_frac^T G R_frac takes the form rot * metric * rot^T.
         rot = REAL(spgr%rotations_subset(:, :, ir), dp)
         gsym = gsym + MATMUL(rot, MATMUL(metric, TRANSPOSE(rot)))
      END DO
      gsym = gsym/REAL(nop, dp)

      lens = [SQRT(gsym(1, 1)), SQRT(gsym(2, 2)), SQRT(gsym(3, 3))]
      IF (MINVAL(lens) < min_cell_length) RETURN
      len_tol = rel_len_tol*MAXVAL(lens)
      ang(1) = ACOS(MAX(-1.0_dp, MIN(1.0_dp, gsym(2, 3)/(lens(2)*lens(3)))))
      ang(2) = ACOS(MAX(-1.0_dp, MIN(1.0_dp, gsym(1, 3)/(lens(1)*lens(3)))))
      ang(3) = ACOS(MAX(-1.0_dp, MIN(1.0_dp, gsym(1, 2)/(lens(1)*lens(2)))))
      n90 = COUNT(ABS(ang - ang_90) < special_ang_tol)

      IF (n90 == 3) THEN
         IF (ABS(lens(1) - lens(2)) < len_tol .AND. ABS(lens(2) - lens(3)) < len_tol) THEN
            symmetry_id = cell_sym_cubic
         ELSE IF (ABS(lens(1) - lens(2)) < len_tol) THEN
            symmetry_id = cell_sym_tetragonal_ab
         ELSE IF (ABS(lens(1) - lens(3)) < len_tol) THEN
            symmetry_id = cell_sym_tetragonal_ac
         ELSE IF (ABS(lens(2) - lens(3)) < len_tol) THEN
            symmetry_id = cell_sym_tetragonal_bc
         ELSE
            symmetry_id = cell_sym_orthorhombic
         END IF
      ELSE IF (n90 == 2) THEN
         IF (ABS(ang(2) - ang_90) >= special_ang_tol) THEN
            ! CP2K monoclinic convention places the unique axis along b (beta /= 90 deg).
            symmetry_id = cell_sym_monoclinic
         ELSE IF (ABS(ang(3) - ang_90) >= special_ang_tol .AND. &
                  ABS(lens(1) - lens(2)) < len_tol) THEN
            IF (ABS(ang(3) - ang_120) < special_ang_tol) THEN
               symmetry_id = cell_sym_hexagonal_gamma_120
            ELSE IF (ABS(ang(3) - ang_60) < special_ang_tol) THEN
               symmetry_id = cell_sym_hexagonal_gamma_60
            ELSE
               symmetry_id = cell_sym_monoclinic_gamma_ab
            END IF
         END IF
      ELSE
         IF (ABS(ang(1) - ang(2)) < eq_angle_tol .AND. ABS(ang(2) - ang(3)) < eq_angle_tol .AND. &
             ABS(ang(1) - ang_90) >= special_ang_tol .AND. &
             ABS(lens(1) - lens(2)) < len_tol .AND. ABS(lens(2) - lens(3)) < len_tol) THEN
            symmetry_id = cell_sym_rhombohedral
         END IF
      END IF

   END FUNCTION spgr_point_group_symmetry_id

! **************************************************************************************************
!> \brief Classify one CONSTRAINT token without aborting (decision half of the parser)
!> \param token one upper-case token of the CONSTRAINT value line
!> \param accepted .TRUE. when the token contributes a mask
!> \param token_mask the constraint mask contributed by an accepted token
!> \param reject_reason empty for accepted tokens, or the abort message otherwise
!> \note parse_cell_constraint applies this classification and turns a rejection
!>       into a cp_abort. Keeping the decision here keeps it unit-testable:
!>       the unit test checks the classification table without triggering aborts.
! **************************************************************************************************
   PURE SUBROUTINE parse_token(token, accepted, token_mask, reject_reason)
      CHARACTER(LEN=*), INTENT(IN)                       :: token
      LOGICAL, INTENT(OUT)                               :: accepted
      INTEGER, INTENT(OUT)                               :: token_mask
      CHARACTER(LEN=default_string_length), INTENT(OUT)  :: reject_reason

      accepted = .TRUE.
      token_mask = 0
      reject_reason = ""
      SELECT CASE (token)
      CASE ("A", "AX")
         token_mask = mask_a
      CASE ("B")
         token_mask = mask_b
      CASE ("BX")
         token_mask = comp_bx
      CASE ("BY")
         token_mask = comp_by
      CASE ("C", "Z")
         token_mask = mask_c
      CASE ("CX")
         token_mask = comp_cx
      CASE ("CY")
         token_mask = comp_cy
      CASE ("CZ")
         token_mask = comp_cz
      CASE ("X")
         token_mask = mask_x
      CASE ("Y")
         token_mask = mask_y
      CASE ("XY")
         token_mask = mask_xy
      CASE ("XZ")
         token_mask = mask_xz
      CASE ("YZ")
         token_mask = mask_yz
      CASE ("ALPHA")
         token_mask = gfix_alpha
      CASE ("BETA")
         token_mask = gfix_beta
      CASE ("GAMMA")
         token_mask = gfix_gamma
      CASE ("AY", "AZ", "BZ")
         accepted = .FALSE.
         reject_reason = "CELL_OPT CONSTRAINT token "//TRIM(token)//" is always zero in "// &
                         "the cell convention (a lies along x, b lies in the xy plane) "// &
                         "and cannot be constrained."
      CASE DEFAULT
         accepted = .FALSE.
         reject_reason = "Unknown CELL_OPT CONSTRAINT token: "//TRIM(token)// &
                         ". Accepted tokens: NONE, the directions X, Y, Z,"// &
                         " XY, XZ, YZ, the vectors A, B, C, the components"// &
                         " AX, BX, BY, CX, CY, CZ, and the angles ALPHA, BETA, GAMMA."
      END SELECT
   END SUBROUTINE parse_token

! **************************************************************************************************
!> \brief Computes the derivatives for the cell
!> \param gradient enthalpy gradient +d(E + pres_ext*V)/dq in the six canonical components
!> \param av_ptens ...
!> \param pres_ext ...
!> \param cell ...
!> \param mtrx ...
!> \param pres_int ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! **************************************************************************************************
   SUBROUTINE get_dg_dh(gradient, av_ptens, pres_ext, cell, mtrx, pres_int)

      REAL(KIND=dp), DIMENSION(n_slot), INTENT(OUT)      :: gradient
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: av_ptens
      REAL(KIND=dp), INTENT(IN)                          :: pres_ext
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: mtrx
      REAL(KIND=dp), INTENT(OUT)                         :: pres_int

      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(3, 3)                     :: correction, pten_hinv_old, ptens

      gradient = 0.0_dp

      ptens = av_ptens

      ! Evaluating the internal pressure
      pres_int = 0.0_dp
      DO i = 1, 3
         pres_int = pres_int + ptens(i, i)
      END DO
      pres_int = pres_int/3.0_dp

      ptens(1, 1) = av_ptens(1, 1) - pres_ext
      ptens(2, 2) = av_ptens(2, 2) - pres_ext
      ptens(3, 3) = av_ptens(3, 3) - pres_ext

      pten_hinv_old = cell%deth*MATMUL(cell%h_inv, ptens)
      correction = MATMUL(mtrx, cell%hmat)

      gradient(slot_ax) = pten_hinv_old(1, 1) - correction(1, 1)
      gradient(slot_bx) = pten_hinv_old(2, 1) - correction(2, 1)
      gradient(slot_by) = pten_hinv_old(2, 2) - correction(2, 2)
      gradient(slot_cx) = pten_hinv_old(3, 1) - correction(3, 1)
      gradient(slot_cy) = pten_hinv_old(3, 2) - correction(3, 2)
      gradient(slot_cz) = pten_hinv_old(3, 3) - correction(3, 3)

      gradient = -gradient

   END SUBROUTINE get_dg_dh

END MODULE cell_opt_utils
