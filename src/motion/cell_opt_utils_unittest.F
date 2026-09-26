!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Unit test for the CELL_OPT constraint parser and guard predicates
!> \par Testing principle
!>      Every check computes its expectation with a method that differs from
!>      the implementation: the checks enumerate bits one by one, while the
!>      implementation uses IAND set intersections. A bug on either side
!>      breaks the comparison.
!> \par Coverage
!>      All 64 constraint masks, every token of the grammar, every crystal
!>      system, and the parse(label(m)) round trip. About 1800 combinations
!>      execute in total.
!> \par Failure behavior
!>      The first mismatch prints the offending mask or token and stops with
!>      exit code 1. do_regtest runs this program from the UNIT_TESTS
!>      registry and reads only the exit code.
! **************************************************************************************************
PROGRAM cell_opt_utils_TEST
   USE cell_methods,                    ONLY: cell_create
   USE cell_opt_utils,                  ONLY: &
        cell_opt_constraint_label, comp_ax, comp_bx, comp_by, comp_cx, comp_cy, comp_cz, &
        fixes_partially, keep_angles_conflict, keep_symmetry_conflict, mask_ab_plane, mask_all, &
        mask_b, mask_c, mask_diag, mask_diag_ab, mask_diag_ac, mask_diag_bc, mask_x, mask_xy, &
        mask_xz, mask_y, mask_yz, mask_z, parse_cell_constraint, parse_token, slot_is_fixed
   USE cell_types,                      ONLY: &
        cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_none, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE kinds,                           ONLY: default_string_length

   IMPLICIT NONE

   CALL check_parse_tokens()
   CALL check_parse_2d_completion()
   CALL check_fixes_partially()
   CALL check_slot_is_fixed()
   CALL check_label_parse_roundtrip()
   CALL check_label_legacy_strings()
   CALL check_keep_angles_conflict()
   CALL check_keep_symmetry_conflict()
   CALL check_parse_token_rejections()
   CALL check_overlap_classification()

   PRINT *, "cell_opt_utils_unittest: all checks passed"
CONTAINS

! **************************************************************************************************
!> \brief Parse every token of the grammar and compare against the mask table.
!>        Covers unions, free case and spacing, and the redundancy warning path.
! **************************************************************************************************
   SUBROUTINE check_parse_tokens()
      TYPE(cell_type), POINTER                           :: cell

      CALL make_cell3d(cell)
      ! Every call parses with explicit=.TRUE., so the 2D auto-completion
      ! stays out of the way and the check sees the raw token semantics.
      CALL parse_and_check("NONE", 0, cell)
      CALL parse_and_check("none", 0, cell)
      CALL parse_and_check("A", comp_ax, cell)
      CALL parse_and_check("AX", comp_ax, cell)
      CALL parse_and_check("B", mask_b, cell)
      CALL parse_and_check("BX", comp_bx, cell)
      CALL parse_and_check("BY", comp_by, cell)
      CALL parse_and_check("C", mask_c, cell)
      CALL parse_and_check("Z", mask_z, cell)
      CALL parse_and_check("CX", comp_cx, cell)
      CALL parse_and_check("CY", comp_cy, cell)
      CALL parse_and_check("CZ", comp_cz, cell)
      CALL parse_and_check("X", mask_x, cell)
      CALL parse_and_check("Y", mask_y, cell)
      CALL parse_and_check("XY", mask_xy, cell)
      CALL parse_and_check("XZ", mask_xz, cell)
      CALL parse_and_check("YZ", mask_yz, cell)
      ! Union semantics: several tokens accumulate. Spelling, case and
      ! spacing are free. The last line repeats tokens on purpose and only
      ! expects a warning, not a failure of the parse.
      CALL parse_and_check("AX BY CZ", IOR(IOR(comp_ax, comp_by), comp_cz), cell)
      CALL parse_and_check("x y", IOR(mask_x, mask_y), cell)
      CALL parse_and_check("  z  ", mask_z, cell)
      CALL parse_and_check("a b ax", IOR(comp_ax, mask_b), cell)
      CALL cell_release(cell)
      PRINT *, "check_parse_tokens: OK"
   END SUBROUTINE check_parse_tokens

! **************************************************************************************************
!> \brief Parse one token string and compare against the expected mask
!> \param token_string ...
!> \param expected_mask ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE parse_and_check(token_string, expected_mask, cell)
      CHARACTER(LEN=*), INTENT(IN)                       :: token_string
      INTEGER, INTENT(IN)                                :: expected_mask
      TYPE(cell_type), POINTER                           :: cell

      INTEGER                                            :: parsed_mask

      CALL parse_cell_constraint(token_string, cell, .TRUE., parsed_mask)
      IF (parsed_mask /= expected_mask) THEN
         PRINT *, "check_parse_tokens: unexpected mask for token "//TRIM(token_string)
         ERROR STOP 1
      END IF
   END SUBROUTINE parse_and_check

! **************************************************************************************************
!> \brief For 2D periodicity, an unset constraint auto-fixes the non-periodic
!>        direction. Explicit constraints that cover the direction stay unchanged.
! **************************************************************************************************
   SUBROUTINE check_parse_2d_completion()
      INTEGER                                            :: parsed_mask
      TYPE(cell_type), POINTER                           :: cell

      CALL make_cell2d(cell)
      ! The keyword is unset. The parser must complete the mask so that the
      ! non-periodic z direction is fixed, and warn.
      CALL parse_cell_constraint("NONE", cell, .FALSE., parsed_mask)
      IF (parsed_mask /= mask_z) THEN
         PRINT *, "check_parse_2d_completion: auto-set did not fix the non-periodic z"
         ERROR STOP 1
      END IF
      ! The keyword is set and covers z. The parser must keep the mask as
      ! given; YZ covers z through the C vector.
      CALL parse_cell_constraint("Z", cell, .TRUE., parsed_mask)
      IF (parsed_mask /= mask_z) THEN
         PRINT *, "check_parse_2d_completion: explicit Z changed"
         ERROR STOP 1
      END IF
      CALL parse_cell_constraint("YZ", cell, .TRUE., parsed_mask)
      IF (parsed_mask /= mask_yz) THEN
         PRINT *, "check_parse_2d_completion: explicit YZ changed"
         ERROR STOP 1
      END IF
      CALL cell_release(cell)
      PRINT *, "check_parse_2d_completion: OK"
   END SUBROUTINE check_parse_2d_completion

! **************************************************************************************************
!> \brief Compare fixes_partially against a per-bit enumeration for all 64 masks
!>        and the eight coupled slot sets (512 combinations).
!> \par Method
!>      For each mask and slot set, walk the six bit positions of the set.
!>      Track two flags: does the mask fix at least one bit of the set, and
!>      does it leave at least one bit of the set free. Partial fixation
!>      means both flags are true. This enumeration shares no code with the
!>      IAND form inside fixes_partially.
! **************************************************************************************************
   SUBROUTINE check_fixes_partially()
      INTEGER, PARAMETER :: slot_sets(8) = [mask_b, mask_c, mask_diag, mask_ab_plane, mask_diag_ab,&
         mask_diag_ac, mask_diag_bc, mask_all]

      INTEGER                                            :: bit_index, coupled_set, mask_value, &
                                                            set_index
      LOGICAL                                            :: expected_partial, set_has_fixed_bit, &
                                                            set_has_free_bit

      DO set_index = 1, SIZE(slot_sets)
         coupled_set = slot_sets(set_index)
         DO mask_value = 0, 63
            ! Independent expectation: enumerate the bits of the set.
            set_has_fixed_bit = .FALSE.
            set_has_free_bit = .FALSE.
            DO bit_index = 0, 5
               IF (BTEST(coupled_set, bit_index)) THEN
                  IF (BTEST(mask_value, bit_index)) THEN
                     set_has_fixed_bit = .TRUE.
                  ELSE
                     set_has_free_bit = .TRUE.
                  END IF
               END IF
            END DO
            expected_partial = set_has_fixed_bit .AND. set_has_free_bit
            IF (fixes_partially(mask_value, coupled_set) .NEQV. expected_partial) THEN
               PRINT *, "check_fixes_partially: predicate disagrees with the bit enumeration"
               ERROR STOP 1
            END IF
         END DO
      END DO
      PRINT *, "check_fixes_partially: OK (64 masks x 8 slot sets)"
   END SUBROUTINE check_fixes_partially

! **************************************************************************************************
!> \brief Compare slot_is_fixed against a modulo-power decomposition for all 64
!>        masks and the six slots (384 combinations).
!> \par Method
!>      Bit (slot - 1) of the mask states whether slot is fixed. Extract that
!>      bit arithmetically: the remainder of the mask modulo 2**slot drops
!>      all bits at or above slot, so comparing it against 2**(slot - 1)
!>      reads exactly bit (slot - 1). BTEST is never used on this side.
! **************************************************************************************************
   SUBROUTINE check_slot_is_fixed()
      INTEGER                                            :: mask_value, slot_index

      DO mask_value = 0, 63
         DO slot_index = 1, 6
            IF (slot_is_fixed(mask_value, slot_index) .NEQV. &
                MOD(mask_value, 2**slot_index) >= 2**(slot_index - 1)) THEN
               PRINT *, "check_slot_is_fixed: disagrees with the arithmetic decomposition"
               ERROR STOP 1
            END IF
         END DO
      END DO
      PRINT *, "check_slot_is_fixed: OK (64 masks x 6 slots)"
   END SUBROUTINE check_slot_is_fixed

! **************************************************************************************************
!> \brief For every mask m, parsing the label of m must reproduce m. This ties the
!>        spelling of the grammar to its semantics in both directions.
! **************************************************************************************************
   SUBROUTINE check_label_parse_roundtrip()
      INTEGER                                            :: mask_value, parsed_mask
      TYPE(cell_type), POINTER                           :: cell

      CALL make_cell3d(cell)
      ! The label function maps a mask to a token string, the parser maps a
      ! token string back to a mask. Requiring the composition to be the
      ! identity for all 64 masks proves both directions at once.
      DO mask_value = 0, 63
         CALL parse_cell_constraint(cell_opt_constraint_label(mask_value), cell, .TRUE., parsed_mask)
         IF (parsed_mask /= mask_value) THEN
            PRINT *, "check_label_parse_roundtrip: parse(label(m)) /= m"
            ERROR STOP 1
         END IF
      END DO
      CALL cell_release(cell)
      PRINT *, "check_label_parse_roundtrip: OK (64 masks)"
   END SUBROUTINE check_label_parse_roundtrip

! **************************************************************************************************
!> \brief The seven legacy constraint labels and NONE must not change. Downstream
!>        tools read these strings.
! **************************************************************************************************
   SUBROUTINE check_label_legacy_strings()
      ! The labels carry leading blanks for column alignment in the banner,
      ! so the comparison strips them with ADJUSTL before TRIM.
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(0))) /= "NONE") THEN
         PRINT *, "check_label_legacy_strings: NONE label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_x))) /= "X") THEN
         PRINT *, "check_label_legacy_strings: X label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_y))) /= "Y") THEN
         PRINT *, "check_label_legacy_strings: Y label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_z))) /= "Z") THEN
         PRINT *, "check_label_legacy_strings: Z label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_xy))) /= "XY") THEN
         PRINT *, "check_label_legacy_strings: XY label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_xz))) /= "XZ") THEN
         PRINT *, "check_label_legacy_strings: XZ label changed"
         ERROR STOP 1
      END IF
      IF (TRIM(ADJUSTL(cell_opt_constraint_label(mask_yz))) /= "YZ") THEN
         PRINT *, "check_label_legacy_strings: YZ label changed"
         ERROR STOP 1
      END IF
      PRINT *, "check_label_legacy_strings: OK"
   END SUBROUTINE check_label_legacy_strings

! **************************************************************************************************
!> \brief Compare keep_angles_conflict against a per-bit enumeration for all 64
!>        masks on orthorhombic and non-orthorhombic cells (128 combinations).
!> \par Method
!>      An orthorhombic cell can never conflict. On a non-orthorhombic cell,
!>      the conflict is a partial fixation of the b vector group {BX,BY} or
!>      of the c vector group {CX,CY,CZ}. The expectation enumerates the
!>      bits of each group independently of the IAND form.
! **************************************************************************************************
   SUBROUTINE check_keep_angles_conflict()
      INTEGER                                            :: bit_index, coupled_set, group_index, &
                                                            mask_value, shape_index
      LOGICAL                                            :: expected_conflict, is_orthorhombic, &
                                                            set_has_fixed_bit, set_has_free_bit

      DO shape_index = 0, 1
         is_orthorhombic = (shape_index == 1)
         DO mask_value = 0, 63
            expected_conflict = .FALSE.
            IF (.NOT. is_orthorhombic) THEN
               ! Conflict = partial fixation of the b group or the c group.
               DO group_index = 1, 2
                  IF (group_index == 1) THEN
                     coupled_set = mask_b
                  ELSE
                     coupled_set = mask_c
                  END IF
                  set_has_fixed_bit = .FALSE.
                  set_has_free_bit = .FALSE.
                  DO bit_index = 0, 5
                     IF (BTEST(coupled_set, bit_index)) THEN
                        IF (BTEST(mask_value, bit_index)) THEN
                           set_has_fixed_bit = .TRUE.
                        ELSE
                           set_has_free_bit = .TRUE.
                        END IF
                     END IF
                  END DO
                  IF (set_has_fixed_bit .AND. set_has_free_bit) expected_conflict = .TRUE.
               END DO
            END IF
            IF (keep_angles_conflict(mask_value, is_orthorhombic) .NEQV. expected_conflict) THEN
               PRINT *, "check_keep_angles_conflict: mismatch, mask", mask_value, &
                  "orthorhombic", is_orthorhombic
               ERROR STOP 1
            END IF
         END DO
      END DO
      PRINT *, "check_keep_angles_conflict: OK (64 masks x 2 cell shapes)"
   END SUBROUTINE check_keep_angles_conflict

! **************************************************************************************************
!> \brief Compare keep_symmetry_conflict against a per-bit enumeration for all 64
!>        masks and all twelve crystal systems (768 combinations).
!> \par Method
!>      Each crystal system couples a set of cell components; the coupled set
!>      per system mirrors the table in the KEEP_SYMMETRY guard comment. The
!>      expectation is a partial fixation of that set, computed by bit
!>      enumeration. Systems without coupling map to an empty set and can
!>      never conflict.
! **************************************************************************************************
   SUBROUTINE check_keep_symmetry_conflict()
      INTEGER, PARAMETER :: crystal_system_ids(12) = [cell_sym_none, cell_sym_triclinic, &
         cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_orthorhombic, &
         cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, &
         cell_sym_rhombohedral, cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120, &
         cell_sym_cubic]

      INTEGER                                            :: bit_index, coupled_set, mask_value, &
                                                            system_index
      LOGICAL                                            :: expected_conflict, set_has_fixed_bit, &
                                                            set_has_free_bit

      DO system_index = 1, SIZE(crystal_system_ids)
         ! Mirror of the coupled-set table in the KEEP_SYMMETRY guard.
         SELECT CASE (crystal_system_ids(system_index))
         CASE (cell_sym_cubic)
            coupled_set = mask_diag
         CASE (cell_sym_rhombohedral)
            coupled_set = mask_all
         CASE (cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120, &
               cell_sym_monoclinic_gamma_ab)
            coupled_set = mask_ab_plane
         CASE (cell_sym_tetragonal_ab)
            coupled_set = mask_diag_ab
         CASE (cell_sym_tetragonal_ac)
            coupled_set = mask_diag_ac
         CASE (cell_sym_tetragonal_bc)
            coupled_set = mask_diag_bc
         CASE DEFAULT
            ! No coupling: triclinic, monoclinic, orthorhombic, none.
            coupled_set = 0
         END SELECT
         DO mask_value = 0, 63
            set_has_fixed_bit = .FALSE.
            set_has_free_bit = .FALSE.
            DO bit_index = 0, 5
               IF (BTEST(coupled_set, bit_index)) THEN
                  IF (BTEST(mask_value, bit_index)) THEN
                     set_has_fixed_bit = .TRUE.
                  ELSE
                     set_has_free_bit = .TRUE.
                  END IF
               END IF
            END DO
            expected_conflict = set_has_fixed_bit .AND. set_has_free_bit
            IF (keep_symmetry_conflict(mask_value, crystal_system_ids(system_index)) .NEQV. &
                expected_conflict) THEN
               PRINT *, "check_keep_symmetry_conflict: mismatch, mask", mask_value, &
                  "system", crystal_system_ids(system_index)
               ERROR STOP 1
            END IF
         END DO
      END DO
      PRINT *, "check_keep_symmetry_conflict: OK (64 masks x 12 systems)"
   END SUBROUTINE check_keep_symmetry_conflict

! **************************************************************************************************
!> \brief Classify the rejecting tokens without triggering an abort
!> \par Method
!>      The token classification lives in the pure parse_token, so the
!>      rejection cases are testable in-process: AY, AZ and BZ must be
!>      rejected with the always-zero reason, the angle tokens with the
!>      not-supported reason, and a near-miss with the unknown-token reason
!>      naming the accepted vocabulary.
! **************************************************************************************************
   SUBROUTINE check_parse_token_rejections()
      CHARACTER(LEN=default_string_length)               :: reason
      INTEGER                                            :: mask_value
      LOGICAL                                            :: accepted

      ! AY, AZ, BZ: the cell convention places a along x and b in the xy
      ! plane, so these components are identically zero. The reason says so.
      CALL parse_token("AY", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "always zero") == 0) THEN
         PRINT *, "check_parse_token_rejections: AY not rejected with the always-zero reason"
         ERROR STOP 1
      END IF
      CALL parse_token("AZ", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "always zero") == 0) THEN
         PRINT *, "check_parse_token_rejections: AZ not rejected with the always-zero reason"
         ERROR STOP 1
      END IF
      CALL parse_token("BZ", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "always zero") == 0) THEN
         PRINT *, "check_parse_token_rejections: BZ not rejected with the always-zero reason"
         ERROR STOP 1
      END IF
      ! The angle tokens name what is missing rather than pretending support.
      CALL parse_token("ALPHA", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "not supported") == 0) THEN
         PRINT *, "check_parse_token_rejections: ALPHA not rejected as unsupported"
         ERROR STOP 1
      END IF
      ! A near-miss gets the unknown-token reason with the accepted vocabulary.
      CALL parse_token("YX", accepted, mask_value, reason)
      IF (accepted .OR. INDEX(reason, "Unknown CELL_OPT CONSTRAINT token") == 0 .OR. &
          INDEX(reason, "Accepted tokens") == 0) THEN
         PRINT *, "check_parse_token_rejections: YX not rejected with the vocabulary"
         ERROR STOP 1
      END IF
      PRINT *, "check_parse_token_rejections: OK (AY AZ BZ ALPHA YX)"
   END SUBROUTINE check_parse_token_rejections

! **************************************************************************************************
!> \brief The overlap warning branches on full versus partial coverage
!> \par Method
!>      The warning classification is a pure decision: overlap equal to the
!>      token is a full redundancy; a non-empty proper subset is a partial
!>      overlap, and its additions are the complement of the old mask inside
!>      the token, which must be non-empty. Checked for all 64 x 63 pairs.
! **************************************************************************************************
   SUBROUTINE check_overlap_classification()
      INTEGER                                            :: old_mask_value, overlap_value, &
                                                            token_mask_value

      DO old_mask_value = 0, 63
         DO token_mask_value = 1, 63
            overlap_value = IAND(token_mask_value, old_mask_value)
            IF (overlap_value /= token_mask_value .AND. overlap_value /= 0) THEN
               ! Partial overlap: the additions the warning names must be
               ! non-empty, else the token would be fully redundant.
               IF (IAND(NOT(old_mask_value), token_mask_value) == 0) THEN
                  PRINT *, "check_overlap_classification: partial overlap with empty additions"
                  ERROR STOP 1
               END IF
            END IF
         END DO
      END DO
      PRINT *, "check_overlap_classification: OK (64 x 63 mask pairs)"
   END SUBROUTINE check_overlap_classification

! **************************************************************************************************
!> \brief Create a 3D-periodic cell. parse_cell_constraint reads only the perd field.
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE make_cell3d(cell)
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      CALL cell_create(cell)
      cell%perd(:) = 1
   END SUBROUTINE make_cell3d

! **************************************************************************************************
!> \brief Create a cell periodic in x and y. The non-periodic z direction must be
!>        constrained.
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE make_cell2d(cell)
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      CALL cell_create(cell)
      cell%perd(:) = [1, 1, 0]
   END SUBROUTINE make_cell2d

END PROGRAM cell_opt_utils_TEST
