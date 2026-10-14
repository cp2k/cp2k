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
!>      breaks the comparison. The projection checks take their expectation
!>      from the guard predicates and their behavior from
!>      apply_cell_constraints, so the hand-mirrored coupling tables and the
!>      executable projections verify each other.
!> \par Coverage
!>      All 64 constraint masks, every token of the grammar, every crystal
!>      system, and the parse(label(m)) round trip. About 1800 combinations
!>      execute in total, plus the guard decisions checked against the real
!>      projections: 64 masks on 14 cell shapes with four gradient draws
!>      each, every draw passed through apply_cell_constraints twice, and
!>      the pressure-criterion predicate over 64 masks, both cell shapes
!>      and both KEEP_ANGLES settings.
!> \par Failure behavior
!>      The first mismatch prints the offending mask or token and stops with
!>      exit code 1. do_regtest runs this program from the UNIT_TESTS
!>      registry and reads only the exit code.
! **************************************************************************************************
PROGRAM cell_opt_utils_TEST
   USE cell_methods,                    ONLY: cell_create
   USE cell_opt_utils,                  ONLY: &
        apply_cell_constraints, cell_opt_constraint_label, comp_ax, comp_bx, comp_by, comp_cx, &
        comp_cy, comp_cz, fixes_partially, keep_angles_conflict, keep_symmetry_conflict, &
        mask_ab_plane, mask_all, mask_b, mask_c, mask_diag, mask_diag_ab, mask_diag_ac, &
        mask_diag_bc, mask_x, mask_xy, mask_xz, mask_y, mask_yz, mask_z, parse_cell_constraint, &
        parse_token, pres_crit_shear_pinned, slot_is_fixed
   USE cell_types,                      ONLY: &
        cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, cell_sym_hexagonal_gamma_60, &
        cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_none, cell_sym_orthorhombic, &
        cell_sym_rhombohedral, cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, &
        cell_sym_tetragonal_bc, cell_sym_triclinic, cell_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp

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
   CALL check_projection_guard_consistency()
   CALL check_pres_crit_shear_pinned()

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

      CHARACTER(LEN=default_string_length)               :: tokens(16)
      INTEGER                                            :: ntok, parsed_mask

      CALL split_tokens(token_string, tokens, ntok)
      CALL parse_cell_constraint(tokens(:ntok), cell, .TRUE., parsed_mask)
      IF (parsed_mask /= expected_mask) THEN
         PRINT *, "check_parse_tokens: unexpected mask for token "//TRIM(token_string)
         ERROR STOP 1
      END IF
   END SUBROUTINE parse_and_check

! **************************************************************************************************
!> \brief Splits a blank-separated string into tokens, mimicking how the input
!>        tokenizer delivers the value of a char_t list keyword
!> \param string ...
!> \param tokens ...
!> \param ntok ...
! **************************************************************************************************
   SUBROUTINE split_tokens(string, tokens, ntok)
      CHARACTER(LEN=*), INTENT(IN)                       :: string
      CHARACTER(LEN=default_string_length), INTENT(OUT)  :: tokens(:)
      INTEGER, INTENT(OUT)                               :: ntok

      CHARACTER(LEN=LEN(string))                         :: rest
      INTEGER                                            :: isep

      rest = ADJUSTL(string)
      ntok = 0
      DO WHILE (LEN_TRIM(rest) > 0)
         isep = INDEX(rest, " ")
         ntok = ntok + 1
         IF (isep == 0) THEN
            tokens(ntok) = TRIM(rest)
            EXIT
         END IF
         tokens(ntok) = rest(:isep - 1)
         rest = ADJUSTL(rest(isep + 1:))
      END DO
   END SUBROUTINE split_tokens

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
      CALL parse_cell_constraint(["NONE"], cell, .FALSE., parsed_mask)
      IF (parsed_mask /= mask_z) THEN
         PRINT *, "check_parse_2d_completion: auto-set did not fix the non-periodic z"
         ERROR STOP 1
      END IF
      ! The keyword is set and covers z. The parser must keep the mask as
      ! given; YZ covers z through the C vector.
      CALL parse_cell_constraint(["Z"], cell, .TRUE., parsed_mask)
      IF (parsed_mask /= mask_z) THEN
         PRINT *, "check_parse_2d_completion: explicit Z changed"
         ERROR STOP 1
      END IF
      CALL parse_cell_constraint(["YZ"], cell, .TRUE., parsed_mask)
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
      CHARACTER(LEN=default_string_length)               :: tokens(16)
      INTEGER                                            :: mask_value, ntok, parsed_mask
      TYPE(cell_type), POINTER                           :: cell

      CALL make_cell3d(cell)
      ! The label function maps a mask to a token string, the parser maps a
      ! token string back to a mask. Requiring the composition to be the
      ! identity for all 64 masks proves both directions at once.
      DO mask_value = 0, 63
         CALL split_tokens(cell_opt_constraint_label(mask_value), tokens, ntok)
         CALL parse_cell_constraint(tokens(:ntok), cell, .TRUE., parsed_mask)
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
!> \brief The guard decisions must agree with the projections they protect
!> \par Method
!>      A guard rejects a mask because the projection and the hard zeroing
!>      would fight over the coupled slots. The signature of that fight is
!>      non-idempotence: with a consistent mask, applying
!>      apply_cell_constraints to its own output reproduces it, and the
!>      output lies in the projected subspace with every fixed slot exactly
!>      zero. With a conflicting mask, the projection pushes components back
!>      into slots the zeroing has emptied, so a second application moves
!>      the result for a generic gradient. The expectation comes from the
!>      guard predicate, the behavior from the projection code: the coupled
!>      set table and the SELECT CASE in apply_cell_constraints verify each
!>      other instead of a third hand-written copy.
!> \note The monoclinic_gamma_ab projection is idempotent only on a cell
!>       with equal a and b lengths, and the KEEP_ANGLES projections divide
!>       by the norms of b and c. The shaped cells below satisfy both.
! **************************************************************************************************
   SUBROUTINE check_projection_guard_consistency()
      INTEGER, PARAMETER :: crystal_system_ids(12) = [cell_sym_none, cell_sym_triclinic, &
         cell_sym_monoclinic, cell_sym_monoclinic_gamma_ab, cell_sym_orthorhombic, &
         cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, &
         cell_sym_rhombohedral, cell_sym_hexagonal_gamma_60, cell_sym_hexagonal_gamma_120, &
         cell_sym_cubic]

      INTEGER                                            :: mask_value, shape_index, system_index
      LOGICAL                                            :: expected_conflict, is_orthorhombic
      TYPE(cell_type), POINTER                           :: cell

      ! Part one: KEEP_ANGLES against keep_angles_conflict, on an
      ! orthorhombic and a triclinic cell. The triclinic cell has nonzero
      ! in-plane components in b and c, so a partial fixation of either
      ! vector is visible in the projections.
      DO shape_index = 1, 2
         is_orthorhombic = (shape_index == 2)
         IF (is_orthorhombic) THEN
            CALL make_projection_cell(cell, .TRUE., cell_sym_orthorhombic)
         ELSE
            CALL make_projection_cell(cell, .FALSE., cell_sym_triclinic)
         END IF
         DO mask_value = 0, 63
            expected_conflict = keep_angles_conflict(mask_value, is_orthorhombic)
            CALL projection_case(cell, .TRUE., .FALSE., mask_value, expected_conflict, &
                                 "check_projection_guard_consistency: KEEP_ANGLES")
         END DO
         CALL cell_release(cell)
      END DO

      ! Part two: KEEP_SYMMETRY against keep_symmetry_conflict, on one
      ! shaped cell per crystal system.
      DO system_index = 1, SIZE(crystal_system_ids)
         CALL make_projection_cell(cell, .FALSE., crystal_system_ids(system_index))
         DO mask_value = 0, 63
            expected_conflict = keep_symmetry_conflict(mask_value, crystal_system_ids(system_index))
            CALL projection_case(cell, .FALSE., .TRUE., mask_value, expected_conflict, &
                                 "check_projection_guard_consistency: KEEP_SYMMETRY")
         END DO
         CALL cell_release(cell)
      END DO
      PRINT *, "check_projection_guard_consistency: OK (64 masks x 14 cell shapes x 4 draws)"
   END SUBROUTINE check_projection_guard_consistency

! **************************************************************************************************
!> \brief One cell and mask: a stationary projection iff the guard sees no conflict
!> \param cell shaped cell carrying orthorhombic flag and symmetry id
!> \param use_keep_angles run the KEEP_ANGLES projection
!> \param use_keep_symmetry run the KEEP_SYMMETRY projection
!> \param mask_value constraint mask under test
!> \param expect_conflict guard decision for this mask and cell
!> \param label check name prefix for failure messages
! **************************************************************************************************
   SUBROUTINE projection_case(cell, use_keep_angles, use_keep_symmetry, mask_value, &
                              expect_conflict, label)
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: use_keep_angles, use_keep_symmetry
      INTEGER, INTENT(IN)                                :: mask_value
      LOGICAL, INTENT(IN)                                :: expect_conflict
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: draw, slot
      LOGICAL                                            :: any_draw_moved
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), DIMENSION(:), POINTER               :: g1, g2

      ALLOCATE (g1(6), g2(6))
      any_draw_moved = .FALSE.
      DO draw = 1, 4
         DO slot = 1, 6
            ! Deterministic entries of order one, varied per draw so an
            ! accidental symmetry of one gradient cannot mask a moved
            ! result. The fractional part keeps entries away from zero.
            g1(slot) = MOD(draw*7 + slot*13, 19) - 9 + 0.37_dp*slot + 0.11_dp*draw
         END DO
         CALL apply_cell_constraints(g1, cell, use_keep_angles, use_keep_symmetry, mask_value)
         g2 = g1
         CALL apply_cell_constraints(g2, cell, use_keep_angles, use_keep_symmetry, mask_value)
         scale = MAX(1.0_dp, MAXVAL(ABS(g1)))
         ! The hard zeroing runs last, so fixed slots are exactly zero in
         ! the output of the first application, whatever the projections
         ! did before it.
         DO slot = 1, 6
            IF (slot_is_fixed(mask_value, slot) .AND. g1(slot) /= 0.0_dp) THEN
               PRINT *, TRIM(label)//": fixed slot nonzero, mask", mask_value, "slot", slot
               ERROR STOP 1
            END IF
         END DO
         IF (expect_conflict) THEN
            IF (MAXVAL(ABS(g2 - g1)) > 1.0E-6_dp*scale) any_draw_moved = .TRUE.
         ELSE
            ! Consistent mask: the second application reproduces the first
            ! up to rounding, and the output lies in the projected
            ! subspace the projections define.
            IF (MAXVAL(ABS(g2 - g1)) > 1.0E-12_dp*scale) THEN
               PRINT *, TRIM(label)//": projection not stationary, mask", mask_value
               ERROR STOP 1
            END IF
            IF (use_keep_angles) CALL check_angles_output(g1, cell, scale, label, mask_value)
            IF (use_keep_symmetry) CALL check_symmetry_output(g1, cell, scale, label, mask_value)
         END IF
      END DO
      IF (expect_conflict .AND. .NOT. any_draw_moved) THEN
         PRINT *, TRIM(label)//": conflict mask never moved the projection, mask", mask_value
         ERROR STOP 1
      END IF
      DEALLOCATE (g1, g2)
   END SUBROUTINE projection_case

! **************************************************************************************************
!> \brief A KEEP_ANGLES output must follow the cell vectors it projects onto
!> \param g projected gradient (six cell components)
!> \param cell shaped cell (hmat and orthorhombic flag are read)
!> \param scale magnitude scale for the comparison tolerances
!> \param label check name prefix for failure messages
!> \param mask_value constraint mask under test
! **************************************************************************************************
   SUBROUTINE check_angles_output(g, cell, scale, label, mask_value)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: g
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: scale
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      INTEGER, INTENT(IN)                                :: mask_value

      IF (cell%orthorhombic) THEN
         ! An orthorhombic cell keeps its off-axis slots exactly zero.
         IF (g(2) /= 0.0_dp .OR. g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": orthorhombic off-axis slot nonzero, mask", mask_value
            ERROR STOP 1
         END IF
      ELSE
         ! b lies in the xy plane: the projected pair is parallel to b when
         ! the two-dimensional cross product vanishes.
         IF (ABS(g(2)*cell%hmat(2, 2) - g(3)*cell%hmat(1, 2)) > 1.0E-12_dp*scale) THEN
            PRINT *, TRIM(label)//": b pair not parallel to b, mask", mask_value
            ERROR STOP 1
         END IF
         ! c spans all three directions: every cross-product component
         ! with c must vanish.
         IF (ABS(g(5)*cell%hmat(3, 3) - g(6)*cell%hmat(2, 3)) > 1.0E-12_dp*scale .OR. &
             ABS(g(6)*cell%hmat(1, 3) - g(4)*cell%hmat(3, 3)) > 1.0E-12_dp*scale .OR. &
             ABS(g(4)*cell%hmat(2, 3) - g(5)*cell%hmat(1, 3)) > 1.0E-12_dp*scale) THEN
            PRINT *, TRIM(label)//": c triple not parallel to c, mask", mask_value
            ERROR STOP 1
         END IF
      END IF
   END SUBROUTINE check_angles_output

! **************************************************************************************************
!> \brief A KEEP_SYMMETRY output must satisfy the coupling of its crystal system
!> \param g projected gradient (six cell components)
!> \param cell shaped cell (symmetry id is read)
!> \param scale magnitude scale for the comparison tolerances
!> \param label check name prefix for failure messages
!> \param mask_value constraint mask under test
!> \note Equalities are asserted where the projection enforces equal slots
!>       (cubic, tetragonal) or fixed ratios (hexagonal). Rhombohedral
!>       couples all six slots with angle-dependent ratios and is carried
!>       by the stationarity check alone. Systems without coupling
!>       (none, triclinic) impose nothing.
! **************************************************************************************************
   SUBROUTINE check_symmetry_output(g, cell, scale, label, mask_value)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: g
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: scale
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      INTEGER, INTENT(IN)                                :: mask_value

      SELECT CASE (cell%symmetry_id)
      CASE (cell_sym_cubic)
         IF (ABS(g(1) - g(3)) > 1.0E-12_dp*scale .OR. ABS(g(1) - g(6)) > 1.0E-12_dp*scale .OR. &
             g(2) /= 0.0_dp .OR. g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": cubic coupling violated, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_tetragonal_ab)
         IF (ABS(g(1) - g(3)) > 1.0E-12_dp*scale .OR. g(2) /= 0.0_dp .OR. &
             g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": tetragonal_ab coupling violated, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_tetragonal_ac)
         IF (ABS(g(1) - g(6)) > 1.0E-12_dp*scale .OR. g(2) /= 0.0_dp .OR. &
             g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": tetragonal_ac coupling violated, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_tetragonal_bc)
         IF (ABS(g(3) - g(6)) > 1.0E-12_dp*scale .OR. g(2) /= 0.0_dp .OR. &
             g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": tetragonal_bc coupling violated, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_orthorhombic)
         IF (g(2) /= 0.0_dp .OR. g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": orthorhombic off-axis slots nonzero, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_monoclinic)
         IF (g(2) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": monoclinic off-axis slots nonzero, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_hexagonal_gamma_60)
         IF (ABS(g(2) - 0.5_dp*g(1)) > 1.0E-12_dp*scale .OR. &
             ABS(g(3) - SQRT(3.0_dp)*g(2)) > 1.0E-12_dp*scale .OR. &
             g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": hexagonal_gamma_60 coupling violated, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_hexagonal_gamma_120)
         IF (ABS(g(2) + 0.5_dp*g(1)) > 1.0E-12_dp*scale .OR. &
             ABS(g(3) + SQRT(3.0_dp)*g(2)) > 1.0E-12_dp*scale .OR. &
             g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": hexagonal_gamma_120 coupling violated, mask", mask_value
            ERROR STOP 1
         END IF
      CASE (cell_sym_monoclinic_gamma_ab)
         IF (g(4) /= 0.0_dp .OR. g(5) /= 0.0_dp) THEN
            PRINT *, TRIM(label)//": monoclinic_gamma_ab off-axis slots nonzero, mask", mask_value
            ERROR STOP 1
         END IF
      CASE DEFAULT
         ! none, triclinic, rhombohedral: carried by the stationarity
         ! check.
      END SELECT
   END SUBROUTINE check_symmetry_output

! **************************************************************************************************
!> \brief The pressure-criterion predicate must fire exactly on the pinned-shear domain
!> \par Method
!>      The oracle enumerates slot bits by position: the diagonal slots
!>      AX, BY, CZ are bits 0, 2 and 5, the shear slots BX, CX, CY are
!>      bits 1, 3 and 4. The criterion is expected to fire on a
!>      non-orthorhombic cell when some but not all diagonal slots are
!>      fixed and either the mask pins a shear slot or KEEP_ANGLES pins
!>      the angles together with a fixed diagonal slot. Zero or all
!>      diagonal slots fixed never fire (the hydrostatic trace, resp. the
!>      external-pressure bypass). The firing counts per cell shape and
!>      KEEP_ANGLES setting are asserted as well, and mask NONE never
!>      fires, so existing KEEP_ANGLE regtests without a constraint stay
!>      warning-free.
! **************************************************************************************************
   SUBROUTINE check_pres_crit_shear_pinned()
      INTEGER                                            :: count_fired, ka_index, mask_value, &
                                                            shape_index
      LOGICAL                                            :: all_diag_fixed, diag_hit, expected, &
                                                            is_orthorhombic, shear_hit, &
                                                            use_keep_angles

      DO shape_index = 0, 1
         is_orthorhombic = (shape_index == 1)
         DO ka_index = 0, 1
            use_keep_angles = (ka_index == 1)
            count_fired = 0
            DO mask_value = 0, 63
               ! Independent expectation from the bit positions.
               all_diag_fixed = BTEST(mask_value, 0) .AND. BTEST(mask_value, 2) .AND. &
                                BTEST(mask_value, 5)
               shear_hit = BTEST(mask_value, 1) .OR. BTEST(mask_value, 3) .OR. &
                           BTEST(mask_value, 4)
               diag_hit = BTEST(mask_value, 0) .OR. BTEST(mask_value, 2) .OR. &
                          BTEST(mask_value, 5)
               expected = .NOT. is_orthorhombic .AND. .NOT. all_diag_fixed .AND. &
                          (shear_hit .OR. (use_keep_angles .AND. diag_hit))
               IF (pres_crit_shear_pinned(mask_value, is_orthorhombic, use_keep_angles) &
                   .NEQV. expected) THEN
                  PRINT *, "check_pres_crit_shear_pinned: mismatch, mask", mask_value, &
                     "orthorhombic", is_orthorhombic, "keep_angles", use_keep_angles
                  ERROR STOP 1
               END IF
               IF (expected) count_fired = count_fired + 1
            END DO
            ! Firing counts of the pinned-shear domain: an orthorhombic
            ! cell never fires; on an oblique cell the mask domain has 56
            ! members of which the all-diagonals-fixed ones are excluded,
            ! and KEEP_ANGLES adds every mask with a fixed diagonal slot.
            IF (is_orthorhombic) THEN
               IF (count_fired /= 0) THEN
                  PRINT *, "check_pres_crit_shear_pinned: fired on an orthorhombic cell"
                  ERROR STOP 1
               END IF
               CYCLE
            END IF
            IF (use_keep_angles) THEN
               IF (count_fired /= 55) THEN
                  PRINT *, "check_pres_crit_shear_pinned: expected 55 masks, got", count_fired
                  ERROR STOP 1
               END IF
            ELSE
               IF (count_fired /= 49) THEN
                  PRINT *, "check_pres_crit_shear_pinned: expected 49 masks, got", count_fired
                  ERROR STOP 1
               END IF
            END IF
         END DO
      END DO
      ! Effective warned set under KEEP_ANGLES: the merged abort guard
      ! (keep_angles_conflict) rejects every mask that fixes part of b or c,
      ! so only whole-vector and a-only masks reach the warning; of those,
      ! NONE stays outside the domain and ABC takes the bypass.
      count_fired = 0
      DO mask_value = 0, 63
         IF (pres_crit_shear_pinned(mask_value, .FALSE., .TRUE.) .AND. &
             .NOT. keep_angles_conflict(mask_value, .FALSE.)) THEN
            count_fired = count_fired + 1
            IF (.NOT. ANY(mask_value == [1, 6, 7, 56, 57, 62])) THEN
               PRINT *, "check_pres_crit_shear_pinned: unexpected survivor", mask_value
               ERROR STOP 1
            END IF
         END IF
      END DO
      IF (count_fired /= 6) THEN
         PRINT *, "check_pres_crit_shear_pinned: expected 6 effective masks, got", count_fired
         ERROR STOP 1
      END IF
      PRINT *, "check_pres_crit_shear_pinned: OK (64 masks x 2 cell shapes x 2 flags, 6 effective)"
   END SUBROUTINE check_pres_crit_shear_pinned

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

! **************************************************************************************************
!> \brief Create a 3D-periodic cell with a lattice consistent with symmetry_id
!> \param cell ...
!> \param is_orthorhombic sets the orthorhombic flag read by KEEP_ANGLES
!> \param symmetry_id crystal system read by KEEP_SYMMETRY
!> \note apply_cell_constraints reads hmat, the orthorhombic flag and
!>       symmetry_id. Each system whose projection reads hmat gets a
!>       consistent lattice: monoclinic_gamma_ab needs equal a and b
!>       lengths, hexagonal needs the named gamma, rhombohedral needs equal
!>       pairwise angles. The other projections are lattice independent and
!>       share the generic triclinic lattice.
! **************************************************************************************************
   SUBROUTINE make_projection_cell(cell, is_orthorhombic, symmetry_id)
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, INTENT(IN)                                :: is_orthorhombic
      INTEGER, INTENT(IN)                                :: symmetry_id

      NULLIFY (cell)
      CALL cell_create(cell)
      cell%perd(:) = 1
      cell%orthorhombic = is_orthorhombic
      cell%symmetry_id = symmetry_id
      SELECT CASE (symmetry_id)
      CASE (cell_sym_orthorhombic)
         cell%hmat = RESHAPE([10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 9.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                              8.0_dp], [3, 3])
      CASE (cell_sym_monoclinic_gamma_ab, cell_sym_hexagonal_gamma_60)
         ! a and b of length 2 with gamma = 60 degree, c orthogonal.
         cell%hmat = RESHAPE([2.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.7320508075688772_dp, 0.0_dp, &
                              0.0_dp, 0.0_dp, 8.0_dp], [3, 3])
      CASE (cell_sym_hexagonal_gamma_120)
         cell%hmat = RESHAPE([2.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 1.7320508075688772_dp, 0.0_dp, &
                              0.0_dp, 0.0_dp, 8.0_dp], [3, 3])
      CASE (cell_sym_rhombohedral)
         ! Three vectors of length 1 with pairwise angles of 60 degree.
         cell%hmat = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.8660254037844386_dp, 0.0_dp, &
                              0.5_dp, 0.2886751345948129_dp, 0.8164965809277260_dp], [3, 3])
      CASE DEFAULT
         ! Generic triclinic lattice: b and c carry nonzero in-plane
         ! components, so partial fixations show up in the projections.
         cell%hmat = RESHAPE([10.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 9.0_dp, 0.0_dp, 0.5_dp, 1.5_dp, &
                              8.0_dp], [3, 3])
      END SELECT
   END SUBROUTINE make_projection_cell

END PROGRAM cell_opt_utils_TEST
