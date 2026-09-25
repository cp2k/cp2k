!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief contains a functional that calculates the energy and its derivatives
!>      for the geometry optimizer
!> \par History
!>      03.2008 - Teodoro Laino [tlaino] - University of Zurich - Cell Optimization
! **************************************************************************************************
MODULE cell_opt_utils
   USE cell_types,                      ONLY: &
        cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE input_constants,                 ONLY: fix_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: sqrt3
   USE mathlib,                         ONLY: angle,&
                                              det_3x3
   USE string_utilities,                ONLY: uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_utils'

   PUBLIC :: fixes_partially, get_dg_dh, keep_angles_conflict, &
             keep_symmetry_coupled_set, keep_symmetry_conflict, parse_token, &
             read_external_press_tensor, &
             apply_cell_constraints, cell_opt_constraint_label, parse_cell_constraint, &
             rescale_new_cell_volume, slot_is_fixed

   ! Masks of the six independent cell variables in the canonical frame
   ! (a along x, b in the xy plane): AX, BX, BY, CX, CY, CZ.
   ! A set bit means that a CONSTRAINT fixes the component.
   ! constraint_id = 0 (fix_none) leaves every component relaxable.
   INTEGER, PARAMETER, PUBLIC :: comp_ax = 2**0, & ! a(x)
                                 comp_bx = 2**1, & ! b(x)
                                 comp_by = 2**2, & ! b(y)
                                 comp_cx = 2**3, & ! c(x)
                                 comp_cy = 2**4, & ! c(y)
                                 comp_cz = 2**5  ! c(z)
   ! Shorthands for complete cell vectors, Cartesian directions and coupled
   ! component sets (unions of disjoint bits, hence representable as sums)
   INTEGER, PARAMETER, PUBLIC :: mask_a = comp_ax, &
                                 mask_b = comp_bx + comp_by, &
                                 mask_c = comp_cx + comp_cy + comp_cz, &
                                 mask_x = comp_ax + comp_bx + comp_cx, &
                                 mask_y = comp_bx + comp_by + comp_cy, &
                                 mask_z = mask_c, &
                                 mask_xy = comp_ax + comp_bx + comp_by + comp_cx + comp_cy, &
                                 mask_xz = comp_ax + comp_bx + comp_cx + comp_cy + comp_cz, &
                                 mask_yz = comp_bx + comp_by + comp_cx + comp_cy + comp_cz, &
                                 mask_diag = comp_ax + comp_by + comp_cz, &
                                 mask_diag_ab = comp_ax + comp_by, &
                                 mask_diag_ac = comp_ax + comp_cz, &
                                 mask_diag_bc = comp_by + comp_cz, &
                                 mask_ab_plane = comp_ax + comp_bx + comp_by, &
                                 mask_all = comp_ax + comp_bx + comp_by + comp_cx + comp_cy + comp_cz

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

      IF (ANY(ABS(pres_ext_tens) > 1.0E-5_dp)) THEN
         mtrx = cell%deth*MATMUL(cell%h_inv, MATMUL(pres_ext_tens, TRANSPOSE(cell%h_inv)))
      END IF

   END SUBROUTINE read_external_press_tensor

! **************************************************************************************************
!> \brief Parses the tokens of the CONSTRAINT keyword into a cell component mask
!> \param constraint_string value line of the CONSTRAINT keyword (lchar_t string)
!> \param cell simulation cell (only the periodicity is used)
!> \param explicit whether the CONSTRAINT keyword was given explicitly
!> \param mask union of the fixed cell components (comp_* and mask_* bits)
! **************************************************************************************************
   SUBROUTINE parse_cell_constraint(constraint_string, cell, explicit, mask)
      CHARACTER(LEN=*), INTENT(IN)                       :: constraint_string
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: explicit
      INTEGER, INTENT(OUT)                               :: mask

      CHARACTER(LEN=default_string_length)               :: reject_reason, token
      CHARACTER(LEN=LEN(constraint_string))              :: line
      INTEGER                                            :: axis_index, isep, n_tokens, &
                                                            nonperiodic_axis, overlap, &
                                                            required_mask, tok_mask
      LOGICAL                                            :: accepted

      ! The lchar_t tokenizer cuts the value at the first quote character,
      ! so a quoted value arrives as an empty string. Give the value
      ! unquoted. Restart echoes are bare and round-trip correctly.
      line = ADJUSTL(constraint_string)

      mask = 0
      n_tokens = 0
      DO WHILE (LEN_TRIM(line) > 0)
         isep = SCAN(line, " "//ACHAR(9))
         IF (isep == 0) THEN
            token = TRIM(line)
            line = ""
         ELSE
            token = TRIM(line(:isep - 1))
            line = ADJUSTL(line(isep + 1:))
         END IF
         n_tokens = n_tokens + 1
         IF (n_tokens > 16) THEN
            CALL cp_abort(__LOCATION__, &
                          "CELL_OPT CONSTRAINT accepts at most 16 tokens.")
         END IF
         CALL uppercase(token)
         IF (token == "NONE") THEN
            ! NONE must be the first and only token: leftover tokens on the
            ! line or an already-accumulated mask mean it is being combined.
            IF (LEN_TRIM(line) > 0 .OR. mask /= 0) THEN
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
      IF (n_tokens == 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "CELL_OPT CONSTRAINT requires a value. "// &
                       "Give the value unquoted: the input tokenizer cuts "// &
                       "a quoted value at the first quote character.")
      END IF

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
!> \brief Returns the label of the consolidated cell constraint mask
!> \param mask bit mask of fixed cell components (comp_* / mask_* bits)
!> \return the historical four-character labels ("   X", "  XY", "NONE") for
!>         the legacy constraint values (padded for backward compatible
!>         printing), the fixed components spelled out ("AX BX BY") otherwise
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
         label = ADJUSTL(label)
      END SELECT

   END FUNCTION cell_opt_constraint_label

! **************************************************************************************************
!> \brief Whether the mask fixes some but not all components of a coupled slot set
!> \param mask bit mask of fixed cell components (comp_* / mask_* bits)
!> \param slots the coupled set of components (a mask_* union)
!> \return ...
!> \note How it compares: IAND(mask, slots) is a set intersection. It holds
!>       the components that the mask fixes inside the coupled set. The
!>       predicate is true when the intersection is non-empty and is not the
!>       full set. In words: the constraint fixes at least one component of
!>       the set, and it leaves at least one component free.
!>       Why the guards need it: a complete freeze and no freeze both keep
!>       the KEEP_ANGLES and KEEP_SYMMETRY projections consistent. A partial
!>       fixation breaks them, so the guards reject exactly this case.
! **************************************************************************************************
   PURE FUNCTION fixes_partially(mask, slots) RESULT(partial)
      INTEGER, INTENT(IN)                                :: mask, slots
      LOGICAL                                            :: partial

      partial = IAND(mask, slots) /= 0 .AND. IAND(mask, slots) /= slots
   END FUNCTION fixes_partially

! **************************************************************************************************
!> \brief Whether the constraint mask fixes the given component slot
!> \param mask bit mask of fixed cell components
!> \param slot component index in the canonical order AX..CZ (1..6)
!> \return ...
!> \note How it compares: BTEST(mask, slot - 1) reads one bit of the mask.
!>       Each bit position corresponds to one component slot, in the order
!>       comp_ax ... comp_cz = 2**0 ... 2**5. The slot argument counts from
!>       one so that it indexes the gradient components directly:
!>       gradient(1) = AX ... gradient(6) = CZ. BTEST counts bits from zero.
! **************************************************************************************************
   PURE FUNCTION slot_is_fixed(mask, slot) RESULT(fixed)
      INTEGER, INTENT(IN)                                :: mask, slot
      LOGICAL                                            :: fixed

      fixed = BTEST(mask, slot - 1)
   END FUNCTION slot_is_fixed

! **************************************************************************************************
!> \brief Whether a CONSTRAINT mask conflicts with the KEEP_ANGLES projection
!> \param mask bit mask of fixed cell components
!> \param orthorhombic whether the cell is orthorhombic
!> \return conflict flag; the caller aborts when it is true
!> \note KEEP_ANGLES projects the gradient slots {BX,BY} onto b and {CX,CY,CZ}
!>       onto c. A conflict arises only on a non-orthorhombic cell, and only
!>       when the mask fixes part of the b vector or part of the c vector.
! **************************************************************************************************
   PURE FUNCTION keep_angles_conflict(mask, orthorhombic) RESULT(conflict)
      INTEGER, INTENT(IN)                                :: mask
      LOGICAL, INTENT(IN)                                :: orthorhombic
      LOGICAL                                            :: conflict

      conflict = .NOT. orthorhombic .AND. &
                 (fixes_partially(mask, mask_b) .OR. fixes_partially(mask, mask_c))
   END FUNCTION keep_angles_conflict

! **************************************************************************************************
!> \brief Whether a CONSTRAINT mask conflicts with the KEEP_SYMMETRY projection
!> \param mask bit mask of fixed cell components
!> \param symmetry_id crystal system of the cell (cell_sym_* constants)
!> \return conflict flag; the caller aborts when it is true
!> \note Each crystal system couples a set of cell components. A conflict
!>       arises when the mask fixes part of such a set. Systems without
!>       coupling never conflict.
! **************************************************************************************************
   PURE FUNCTION keep_symmetry_conflict(mask, symmetry_id) RESULT(conflict)
      INTEGER, INTENT(IN)                                :: mask, symmetry_id
      LOGICAL                                            :: conflict

      INTEGER                                            :: coupled_set

      coupled_set = keep_symmetry_coupled_set(symmetry_id)
      conflict = coupled_set /= 0 .AND. fixes_partially(mask, coupled_set)
   END FUNCTION keep_symmetry_conflict

! **************************************************************************************************
!> \brief The set of cell components coupled by a crystal system (for diagnostics)
!> \param symmetry_id crystal system of the cell (cell_sym_* constants)
!> \return mask of the coupled set; zero for systems without coupling
!> \note Mirrors the SELECT CASE in keep_symmetry_conflict. The abort messages
!>       print this label so the user can see which components the symmetry
!>       couples against their constraint.
! **************************************************************************************************
   PURE FUNCTION keep_symmetry_coupled_set(symmetry_id) RESULT(coupled_set)
      INTEGER, INTENT(IN)                                :: symmetry_id
      INTEGER                                            :: coupled_set

      SELECT CASE (symmetry_id)
      CASE (cell_sym_cubic)
         coupled_set = mask_diag
      CASE (cell_sym_rhombohedral)
         coupled_set = mask_all
      CASE (cell_sym_hexagonal_gamma_60, &
            cell_sym_hexagonal_gamma_120, &
            cell_sym_monoclinic_gamma_ab)
         coupled_set = mask_ab_plane
      CASE (cell_sym_tetragonal_ab)
         coupled_set = mask_diag_ab
      CASE (cell_sym_tetragonal_ac)
         coupled_set = mask_diag_ac
      CASE (cell_sym_tetragonal_bc)
         coupled_set = mask_diag_bc
      CASE DEFAULT
         coupled_set = 0
      END SELECT
   END FUNCTION keep_symmetry_coupled_set

! **************************************************************************************************
!> \brief Classify one CONSTRAINT token without aborting (decision half of the parser)
!> \param token one upper-case token of the CONSTRAINT value line
!> \param accepted .TRUE. when the token contributes a mask
!> \param token_mask the component mask contributed by an accepted token
!> \param reject_reason empty for accepted tokens; the abort message otherwise
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
      CASE ("ALPHA", "BETA", "GAMMA")
         accepted = .FALSE.
         reject_reason = "CELL_OPT CONSTRAINT on angles (ALPHA, BETA, GAMMA) is not supported yet."
      CASE ("AY", "AZ", "BZ")
         accepted = .FALSE.
         reject_reason = "CELL_OPT CONSTRAINT token "//TRIM(token)//" is always zero in "// &
                         "the cell convention (a lies along x, b lies in the xy plane) "// &
                         "and cannot be constrained."
      CASE DEFAULT
         accepted = .FALSE.
         reject_reason = "Unknown CELL_OPT CONSTRAINT token: "//TRIM(token)// &
                         ". Accepted tokens: NONE, the directions X, Y, Z,"// &
                         " XY, XZ, YZ, the vectors A, B, C, and the components"// &
                         " AX, BX, BY, CX, CY, CZ. Several tokens may be combined."
      END SELECT
   END SUBROUTINE parse_token

! **************************************************************************************************
!> \brief Computes the derivatives for the cell
!> \param gradient ...
!> \param av_ptens ...
!> \param pres_ext ...
!> \param cell ...
!> \param mtrx ...
!> \param keep_angles ...
!> \param keep_symmetry ...
!> \param pres_int ...
!> \param pres_constr ...
!> \param constraint_id ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 03.2008
! **************************************************************************************************
   SUBROUTINE get_dg_dh(gradient, av_ptens, pres_ext, cell, mtrx, keep_angles, &
                        keep_symmetry, pres_int, pres_constr, constraint_id)

      REAL(KIND=dp), DIMENSION(:), POINTER               :: gradient
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: av_ptens
      REAL(KIND=dp), INTENT(IN)                          :: pres_ext
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: mtrx
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_angles, keep_symmetry
      REAL(KIND=dp), INTENT(OUT)                         :: pres_int, pres_constr
      INTEGER, INTENT(IN), OPTIONAL                      :: constraint_id

      INTEGER                                            :: i, my_constraint_id
      LOGICAL                                            :: my_keep_angles, my_keep_symmetry
      LOGICAL, DIMENSION(3)                              :: diag_free
      REAL(KIND=dp), DIMENSION(3)                        :: diag_ptens
      REAL(KIND=dp), DIMENSION(3, 3)                     :: correction, pten_hinv_old, ptens

      my_keep_angles = .FALSE.
      IF (PRESENT(keep_angles)) my_keep_angles = keep_angles
      my_keep_symmetry = .FALSE.
      IF (PRESENT(keep_symmetry)) my_keep_symmetry = keep_symmetry
      gradient = 0.0_dp
      IF (PRESENT(constraint_id)) THEN
         my_constraint_id = constraint_id
      ELSE
         my_constraint_id = fix_none
      END IF

      gradient = 0.0_dp

      ptens = av_ptens

      ! Evaluating the internal pressure
      pres_int = 0.0_dp
      DO i = 1, 3
         pres_int = pres_int + ptens(i, i)
      END DO
      pres_int = pres_int/3.0_dp

      ! Nominal pressure averaged over the relaxed diagonal cell components
      ! (AX, BY, CZ couple to the diagonal of the pressure tensor). If all
      ! diagonal components are fixed, no pressure criterion can be relaxed
      ! and the external pressure trivially satisfies the tolerance.
      diag_free = [IAND(my_constraint_id, comp_ax) == 0, &
                   IAND(my_constraint_id, comp_by) == 0, &
                   IAND(my_constraint_id, comp_cz) == 0]
      diag_ptens = [ptens(1, 1), ptens(2, 2), ptens(3, 3)]
      IF (ANY(diag_free)) THEN
         pres_constr = SUM(diag_ptens, MASK=diag_free)/COUNT(diag_free)
      ELSE
         pres_constr = pres_ext
      END IF

      ptens(1, 1) = av_ptens(1, 1) - pres_ext
      ptens(2, 2) = av_ptens(2, 2) - pres_ext
      ptens(3, 3) = av_ptens(3, 3) - pres_ext

      pten_hinv_old = cell%deth*MATMUL(cell%h_inv, ptens)
      correction = MATMUL(mtrx, cell%hmat)

      gradient(1) = pten_hinv_old(1, 1) - correction(1, 1)
      gradient(2) = pten_hinv_old(2, 1) - correction(2, 1)
      gradient(3) = pten_hinv_old(2, 2) - correction(2, 2)
      gradient(4) = pten_hinv_old(3, 1) - correction(3, 1)
      gradient(5) = pten_hinv_old(3, 2) - correction(3, 2)
      gradient(6) = pten_hinv_old(3, 3) - correction(3, 3)

      CALL apply_cell_constraints(gradient, cell, my_keep_angles, my_keep_symmetry, my_constraint_id)

      gradient = -gradient

   END SUBROUTINE get_dg_dh

! **************************************************************************************************
!> \brief Apply cell constraints
!> \param gradient ...
!> \param cell ...
!> \param keep_angles ...
!> \param keep_symmetry ...
!> \param constraint_id ...
!> \author Matthias Krack (October 26, 2017, MK)
! **************************************************************************************************
   SUBROUTINE apply_cell_constraints(gradient, cell, keep_angles, keep_symmetry, constraint_id)

      REAL(KIND=dp), DIMENSION(:), POINTER               :: gradient
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: keep_angles, keep_symmetry
      INTEGER, INTENT(IN)                                :: constraint_id

      INTEGER                                            :: slot
      REAL(KIND=dp)                                      :: a, a_length, ab_length, b_length, cosa, &
                                                            cosah, cosg, deriv_gamma, g, gamma, &
                                                            norm, norm_b, norm_c, sina, sinah, sing

      IF (keep_angles) THEN
         ! If we want to keep the angles constant we have to project out the
         ! components of the cell angles
         norm_b = DOT_PRODUCT(cell%hmat(:, 2), cell%hmat(:, 2))
         norm = DOT_PRODUCT(cell%hmat(1:2, 2), gradient(2:3))
         gradient(2:3) = cell%hmat(1:2, 2)/norm_b*norm
         norm_c = DOT_PRODUCT(cell%hmat(:, 3), cell%hmat(:, 3))
         norm = DOT_PRODUCT(cell%hmat(1:3, 3), gradient(4:6))
         gradient(4:6) = cell%hmat(1:3, 3)/norm_c*norm
         ! Retain an exact orthorhombic cell
         ! (off-diagonal elements must remain zero identically to keep QS fast)
         IF (cell%orthorhombic) THEN
            gradient(2) = 0.0_dp
            gradient(4) = 0.0_dp
            gradient(5) = 0.0_dp
         END IF
      END IF

      IF (keep_symmetry) THEN
         SELECT CASE (cell%symmetry_id)
         CASE (cell_sym_cubic, &
               cell_sym_tetragonal_ab, &
               cell_sym_tetragonal_ac, &
               cell_sym_tetragonal_bc, &
               cell_sym_orthorhombic)
            SELECT CASE (cell%symmetry_id)
            CASE (cell_sym_cubic)
               g = (gradient(1) + gradient(3) + gradient(6))/3.0_dp
               gradient(1) = g
               gradient(3) = g
               gradient(6) = g
            CASE (cell_sym_tetragonal_ab, &
                  cell_sym_tetragonal_ac, &
                  cell_sym_tetragonal_bc)
               SELECT CASE (cell%symmetry_id)
               CASE (cell_sym_tetragonal_ab)
                  g = 0.5_dp*(gradient(1) + gradient(3))
                  gradient(1) = g
                  gradient(3) = g
               CASE (cell_sym_tetragonal_ac)
                  g = 0.5_dp*(gradient(1) + gradient(6))
                  gradient(1) = g
                  gradient(6) = g
               CASE (cell_sym_tetragonal_bc)
                  g = 0.5_dp*(gradient(3) + gradient(6))
                  gradient(3) = g
                  gradient(6) = g
               END SELECT
            CASE (cell_sym_orthorhombic)
               ! Nothing else to do
            END SELECT
            gradient(2) = 0.0_dp
            gradient(4) = 0.0_dp
            gradient(5) = 0.0_dp
         CASE (cell_sym_hexagonal_gamma_60)
            g = 0.5_dp*(gradient(1) + 0.5_dp*(gradient(2) + sqrt3*gradient(3)))
            gradient(1) = g
            gradient(2) = 0.5_dp*g
            gradient(3) = sqrt3*gradient(2)
            gradient(4) = 0.0_dp
            gradient(5) = 0.0_dp
         CASE (cell_sym_hexagonal_gamma_120)
            g = 0.5_dp*(gradient(1) - 0.5_dp*(gradient(2) - sqrt3*gradient(3)))
            gradient(1) = g
            gradient(2) = -0.5_dp*g
            gradient(3) = -sqrt3*gradient(2)
            gradient(4) = 0.0_dp
            gradient(5) = 0.0_dp
         CASE (cell_sym_rhombohedral)
            a = (angle(cell%hmat(:, 3), cell%hmat(:, 2)) + &
                 angle(cell%hmat(:, 1), cell%hmat(:, 3)) + &
                 angle(cell%hmat(:, 1), cell%hmat(:, 2)))/3.0_dp
            cosa = COS(a)
            sina = SIN(a)
            cosah = COS(0.5_dp*a)
            sinah = SIN(0.5_dp*a)
            norm = cosa/cosah
            norm_c = SQRT(1.0_dp - norm*norm)
            g = (gradient(1) + gradient(2)*cosa + gradient(3)*sina + &
                 gradient(4)*cosah*norm + gradient(5)*sinah*norm + gradient(6)*norm_c)/3.0_dp
            gradient(1) = g
            gradient(2) = g*cosa
            gradient(3) = g*sina
            gradient(4) = g*cosah*norm
            gradient(5) = g*sinah*norm
            gradient(6) = g*norm_c
         CASE (cell_sym_monoclinic)
            gradient(2) = 0.0_dp
            gradient(5) = 0.0_dp
         CASE (cell_sym_monoclinic_gamma_ab)
            ! Cell symmetry with a = b, alpha = beta = 90 degree and gammma not equal 90 degree
            a_length = SQRT(cell%hmat(1, 1)*cell%hmat(1, 1) + &
                            cell%hmat(2, 1)*cell%hmat(2, 1) + &
                            cell%hmat(3, 1)*cell%hmat(3, 1))
            b_length = SQRT(cell%hmat(1, 2)*cell%hmat(1, 2) + &
                            cell%hmat(2, 2)*cell%hmat(2, 2) + &
                            cell%hmat(3, 2)*cell%hmat(3, 2))
            ab_length = 0.5_dp*(a_length + b_length)
            gamma = angle(cell%hmat(:, 1), cell%hmat(:, 2))
            cosg = COS(gamma)
            sing = SIN(gamma)
            ! Here, g is the average derivative of the cell vector length ab_length,
            ! and deriv_gamma is the derivative of the angle gamma
            g = 0.5_dp*(gradient(1) + cosg*gradient(2) + sing*gradient(3))
            deriv_gamma = (gradient(3)*cosg - gradient(2)*sing)/b_length
            gradient(1) = g
            gradient(2) = g*cosg - ab_length*sing*deriv_gamma
            gradient(3) = g*sing + ab_length*cosg*deriv_gamma
            gradient(4) = 0.0_dp
            gradient(5) = 0.0_dp
         CASE (cell_sym_triclinic)
            ! Nothing to do
         END SELECT
      END IF

      ! Invariant: the projections above run first, the hard constraint
      ! zeroing runs last. KEEP_VOLUME is handled separately by
      ! rescale_new_cell_volume. No projected gradient survives into a fixed
      ! slot. The cell block is always the last six components of the
      ! optimization vector.
      DO slot = 1, 6
         IF (slot_is_fixed(constraint_id, slot)) gradient(slot) = 0.0_dp
      END DO

   END SUBROUTINE apply_cell_constraints

! **************************************************************************************************
!> \brief Rescale x(idg+1:idg+6) according to vol
!> \param vol ...
!> \param x ...
!> \param idg ...
! **************************************************************************************************
   SUBROUTINE rescale_new_cell_volume(vol, x, idg)

      REAL(KIND=dp), INTENT(IN)                          :: vol
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x
      INTEGER, INTENT(IN)                                :: idg

      INTEGER                                            :: i, j, my_idg
      REAL(KIND=dp)                                      :: ratio
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat_tmp

      CPASSERT((SIZE(x) == idg + 6))
      my_idg = idg

      hmat_tmp(:, :) = 0.0_dp
      DO i = 1, 3
         DO j = 1, i
            my_idg = my_idg + 1
            hmat_tmp(j, i) = x(my_idg)
         END DO
      END DO

      ratio = vol/ABS(det_3x3(hmat_tmp))
      ratio = ratio**(1.0_dp/3.0_dp) ! no need to use EXP(LOG(ratio)/3) for cube root

      my_idg = idg
      DO i = 1, 3
         DO j = 1, i
            my_idg = my_idg + 1
            x(my_idg) = x(my_idg)*ratio
         END DO
      END DO

   END SUBROUTINE rescale_new_cell_volume

END MODULE cell_opt_utils
