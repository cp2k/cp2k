!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains type used for a Simulation Cell Optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
MODULE cell_opt_types

   USE cell_manifold_utils,             ONLY: &
        cell_manifold_init, cell_manifold_is_redundant, cell_manifold_pack_hmat, &
        cell_manifold_pin_semidefinite, cell_manifold_status_label, cell_manifold_type, comp_bx, &
        comp_cx, comp_cy, gfix_all, gfix_alpha, gfix_beta, gfix_gamma, gfix_v, info_licq, info_ok, &
        mask_all, n_constraint_bits, n_geom, n_slot
   USE cell_methods,                    ONLY: canonicalize_cell_matrix,&
                                              cell_create
   USE cell_opt_utils,                  ONLY: cell_opt_constraint_label,&
                                              parse_cell_constraint,&
                                              read_external_press_tensor,&
                                              symmetry_constraint_rows
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_release,&
                                              cell_sym_none,&
                                              cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_constants,                 ONLY: fix_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE string_utilities,                ONLY: uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_types'
   REAL(KIND=dp), PARAMETER, PRIVATE :: canonical_frame_warn_tol = 1.0E-6_dp

   PUBLIC :: cell_opt_env_type, &
             cell_opt_env_create, &
             cell_opt_env_release, &
             cell_opt_setup_manifold, &
             cell_opt_fold_symmetry_rows

! **************************************************************************************************
!> \brief Type containing all informations abour the simulation cell optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   TYPE cell_opt_env_type
      ! Simulation cell optimization parameters
      ! Bit mask of cell components and angles from CONSTRAINT (comp_* and gfix_*, 0 = fix_none).
      INTEGER                                    :: constraint_id = fix_none
      LOGICAL                                    :: keep_angles = .FALSE., &
                                                    keep_symmetry = .FALSE., &
                                                    keep_volume = .FALSE.
      REAL(KIND=dp)                              :: pres_ext = 0.0_dp, pres_int = 0.0_dp, pres_tol = 0.0_dp, pres_constr = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)             :: mtrx = 0.0_dp
      ! Fixed rotations between the user's Cartesian frame and the canonical
      ! frame used for the six independent cell optimization variables.
      REAL(KIND=dp), DIMENSION(3, 3)             :: input_to_opt = 0.0_dp, &
                                                    opt_to_input = 0.0_dp
      TYPE(cell_type), POINTER                   :: opt_cell => NULL()
      ! Active cell constraints and reduced-coordinate chart.
      TYPE(cell_manifold_type)                   :: manifold
   END TYPE cell_opt_env_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cell_env ...
!> \param force_env ...
!> \param geo_section ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   SUBROUTINE cell_opt_env_create(cell_env, force_env, geo_section)
      TYPE(cell_opt_env_type), INTENT(OUT)               :: cell_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: geo_section

      CHARACTER(LEN=default_string_length)               :: label
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: canonical_tokens, constraint_tokens
      INTEGER                                            :: i, output_unit
      LOGICAL                                            :: constraint_explicit, &
                                                            has_frame_dep_constraint
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (cell_env%opt_cell, cell, canonical_tokens, constraint_tokens)
      CALL force_env_get(force_env, cell=cell)
      CALL cell_create(cell_env%opt_cell)
      CALL cell_clone(cell, cell_env%opt_cell, tag="CELL_OPT_CANONICAL")
      CALL canonicalize_cell_matrix(cell_env%opt_cell)
      ! Keep the physical cell in the input orientation while the optimizer
      ! works with the canonical cell and its six non-redundant components.
      cell_env%input_to_opt = MATMUL(cell_env%opt_cell%hmat, cell%h_inv)
      cell_env%opt_to_input = MATMUL(cell%hmat, cell_env%opt_cell%h_inv)
      CALL section_vals_val_get(geo_section, "KEEP_VOLUME", l_val=cell_env%keep_volume)
      CALL section_vals_val_get(geo_section, "KEEP_ANGLES", l_val=cell_env%keep_angles)
      CALL section_vals_val_get(geo_section, "KEEP_SYMMETRY", l_val=cell_env%keep_symmetry)
      CALL section_vals_val_get(geo_section, "PRESSURE_TOLERANCE", r_val=cell_env%pres_tol)
      CALL section_vals_val_get(geo_section, "CONSTRAINT", &
                                c_vals=constraint_tokens, explicit=constraint_explicit)
      CALL parse_cell_constraint(constraint_tokens, cell, constraint_explicit, &
                                 cell_env%constraint_id)

      IF (cell_env%keep_symmetry .AND. cell_env%opt_cell%symmetry_id == cell_sym_none) THEN
         ! &SUBSYS/&CELL SYMMETRY defaults to NONE and is never auto-detected.
         CALL cp_warn(__LOCATION__, &
                      "CELL_OPT: KEEP_SYMMETRY is active, but &SUBSYS/&CELL SYMMETRY "// &
                      "is NONE. CELL_OPT applies no symmetry constraint. Set "// &
                      "&SUBSYS/&CELL SYMMETRY to specify the crystal system.")
      END IF

      has_frame_dep_constraint = (.NOT. constraint_explicit .AND. &
                                  cell_env%constraint_id /= fix_none)
      ! Only for explicitly set values: an unset keyword must stay unparsed so
      ! the restart echo keeps omitting the line, as it always has.
      IF (constraint_explicit) THEN
         ALLOCATE (canonical_tokens(SIZE(constraint_tokens)))
         canonical_tokens = constraint_tokens
         CALL uppercase(canonical_tokens)
         DO i = 1, SIZE(canonical_tokens)
            SELECT CASE (TRIM(ADJUSTL(canonical_tokens(i))))
            CASE ("NONE", "A", "B", "C", "ALPHA", "BETA", "GAMMA")
               CONTINUE
            CASE DEFAULT
               has_frame_dep_constraint = .TRUE.
            END SELECT
         END DO
         CALL section_vals_val_set(geo_section, "CONSTRAINT", c_vals_ptr=canonical_tokens)
      END IF

      IF (has_frame_dep_constraint .AND. &
          MAXVAL(ABS(cell%hmat - cell_env%opt_cell%hmat)) > &
          canonical_frame_warn_tol*MAX(1.0_dp, MAXVAL(ABS(cell%hmat)))) THEN
         CALL cp_warn(__LOCATION__, &
                      "CELL_OPT: Axis and component constraints act on the canonical "// &
                      "upper-triangular cell, not on the rotated input axes.")
      END IF

      ! Get the external pressure
      CALL read_external_press_tensor(geo_section, cell_env%opt_cell, cell_env%pres_ext, &
                                      cell_env%mtrx, cell_env%input_to_opt)

      ! Print cell optimisation setup
      NULLIFY (logger)
      logger => cp_get_default_logger()
      output_unit = cp_print_key_unit_nr(logger, geo_section, "PRINT%CELL", extension=".Log")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A,T61,F20.1)") &
            "CELL_OPT| Pressure tolerance [bar]: ", cp_unit_from_cp2k(cell_env%pres_tol, "bar")
         IF (cell_env%keep_volume) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep volume of cell: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep volume of cell: ", " NO"
         END IF
         IF (cell_env%keep_angles) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep angles between the cell vectors: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep angles between the cell vectors: ", " NO"
         END IF
         IF (cell_env%keep_symmetry) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep cell symmetry: ", "YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T78,A3)") &
               "CELL_OPT| Keep cell symmetry: ", " NO"
         END IF
         label = cell_opt_constraint_label(cell_env%constraint_id)
         IF (COUNT(cell%perd /= 0) == 2 .AND. .NOT. constraint_explicit) THEN
            label = TRIM(label)//" (auto)"
         END IF
         WRITE (UNIT=output_unit, FMT="(T2,A,T77,A)") &
            "CELL_OPT| Constraint: ", TRIM(label)
      END IF
      CALL cp_print_key_finished_output(output_unit, logger, geo_section, "PRINT%CELL")

   END SUBROUTINE cell_opt_env_create

! **************************************************************************************************
!> \brief Initialize cell constraints from CONSTRAINT, KEEP_SYMMETRY, KEEP_ANGLES, and KEEP_VOLUME.
!> \param cell_env cell optimization environment to populate
!> \param emit_warnings optional flag to emit user warnings (default .TRUE.)
! **************************************************************************************************
   SUBROUTINE cell_opt_setup_manifold(cell_env, emit_warnings)
      TYPE(cell_opt_env_type), INTENT(INOUT)             :: cell_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: emit_warnings

      INTEGER                                            :: angle_bits, info
      LOGICAL                                            :: do_warn
      REAL(KIND=dp)                                      :: q0(n_slot)

      do_warn = .TRUE.
      IF (PRESENT(emit_warnings)) do_warn = emit_warnings
      q0 = cell_manifold_pack_hmat(cell_env%opt_cell%hmat)

      ! Anchor q0 with explicit component rows, then fold explicit angle rows.
      CALL cell_manifold_init(cell_env%manifold, q0, IAND(cell_env%constraint_id, mask_all), info)
      IF (info /= info_ok) THEN
         CALL cp_abort(__LOCATION__, &
                       "CELL_OPT: initial cell or CONSTRAINT is invalid ("// &
                       TRIM(cell_manifold_status_label(info))//").")
      END IF
      CALL try_add_mask_bits(cell_env%manifold, IAND(cell_env%constraint_id, gfix_all), &
                             "CONSTRAINT", do_warn)

      IF (cell_env%keep_symmetry) THEN
         CALL cell_opt_fold_symmetry_rows(cell_env, cell_env%opt_cell%symmetry_id, &
                                          emit_warnings=do_warn)
      END IF

      IF (cell_env%keep_angles) THEN
         ! On an orthorhombic cell, BX = CX = CY = 0 fixes all angles to 90 deg.
         ! Off-diagonal hmat entries stay 0.0_dp so cell%orthorhombic stays TRUE.
         IF (cell_env%opt_cell%orthorhombic) THEN
            angle_bits = comp_bx + comp_cx + comp_cy
         ELSE
            angle_bits = gfix_alpha + gfix_beta + gfix_gamma
         END IF
         CALL try_add_mask_bits(cell_env%manifold, angle_bits, "KEEP_ANGLES", do_warn)
      END IF

      IF (cell_env%keep_volume) THEN
         CALL try_add_mask_bits(cell_env%manifold, gfix_v, "KEEP_VOLUME", do_warn)
      END IF

   END SUBROUTINE cell_opt_setup_manifold

! **************************************************************************************************
!> \brief Add independent constraint rows for crystal system symmetry_id.
!> \param cell_env cell optimization environment
!> \param symmetry_id crystal system identifier (cell_sym_*)
!> \param emit_warnings optional flag to emit user warnings (default .TRUE.)
!> \note This routine skips redundant rows and aborts on unresolvable singularities.
!>       When a row is second-order semidefinite at q0, it pins the implied
!>       directions and emits a warning. KEEP_SPACE_GROUP passes the point-group
!>       crystal system here because the point group fixes the same cell metric.
! **************************************************************************************************
   SUBROUTINE cell_opt_fold_symmetry_rows(cell_env, symmetry_id, emit_warnings)
      TYPE(cell_opt_env_type), INTENT(INOUT)             :: cell_env
      INTEGER, INTENT(IN)                                :: symmetry_id
      LOGICAL, INTENT(IN), OPTIONAL                      :: emit_warnings

      INTEGER                                            :: cand_eq_pair(2, n_geom), cand_n_eq, &
                                                            cand_n_lin, i, info, n_next, &
                                                            probe_eq_pair(2, n_geom), sym_bits
      LOGICAL                                            :: do_warn
      REAL(KIND=dp)                                      :: cand_lin_w(n_slot, n_slot), &
                                                            probe_lin_w(n_slot, n_slot)
      TYPE(cell_manifold_type)                           :: m_probe

      do_warn = .TRUE.
      IF (PRESENT(emit_warnings)) do_warn = emit_warnings

      CALL symmetry_constraint_rows(symmetry_id, sym_bits, cand_lin_w, cand_n_lin, &
                                    cand_eq_pair, cand_n_eq)
      CALL try_add_mask_bits(cell_env%manifold, sym_bits, "symmetry component constraint", do_warn)

      DO i = 1, cand_n_lin
         n_next = cell_env%manifold%n_lin_rows + 1
         IF (n_next > n_slot) EXIT
         probe_lin_w = cell_env%manifold%lin_w
         probe_lin_w(:, n_next) = cand_lin_w(:, i)
         CALL cell_manifold_init( &
            m_probe, cell_env%manifold%q0, cell_env%manifold%mask, info, &
            lin_w=probe_lin_w(:, 1:n_next), &
            eq_pair=cell_env%manifold%eq_pair(:, 1:cell_env%manifold%n_eq_rows))
         IF (info == info_ok) THEN
            cell_env%manifold = m_probe
         ELSE
            CALL resolve_degenerate_constraint(cell_env%manifold, info, &
                                               "symmetry linear constraint", 0, do_warn, &
                                               lin_w=cand_lin_w(:, i:i))
         END IF
      END DO

      DO i = 1, cand_n_eq
         n_next = cell_env%manifold%n_eq_rows + 1
         IF (n_next > n_geom - 1) EXIT
         probe_eq_pair = cell_env%manifold%eq_pair
         probe_eq_pair(:, n_next) = cand_eq_pair(:, i)
         CALL cell_manifold_init( &
            m_probe, cell_env%manifold%q0, cell_env%manifold%mask, info, &
            lin_w=cell_env%manifold%lin_w(:, 1:cell_env%manifold%n_lin_rows), &
            eq_pair=probe_eq_pair(:, 1:n_next))
         IF (info == info_ok) THEN
            cell_env%manifold = m_probe
         ELSE
            CALL resolve_degenerate_constraint(cell_env%manifold, info, &
                                               "symmetry equality constraint", 0, do_warn, &
                                               eq_pair=cand_eq_pair(:, i:i))
         END IF
      END DO

   END SUBROUTINE cell_opt_fold_symmetry_rows

! **************************************************************************************************
!> \brief Add each independent bit from candidate_bits into m%mask.
!> \param m active constraint state
!> \param candidate_bits comp_* or gfix_* flags to test one bit at a time
!> \param source_label constraint source label for warning and abort messages
!> \param do_warn .TRUE. to emit a warning when this routine pins a row
! **************************************************************************************************
   SUBROUTINE try_add_mask_bits(m, candidate_bits, source_label, do_warn)
      TYPE(cell_manifold_type), INTENT(INOUT)            :: m
      INTEGER, INTENT(IN)                                :: candidate_bits
      CHARACTER(LEN=*), INTENT(IN)                       :: source_label
      LOGICAL, INTENT(IN)                                :: do_warn

      INTEGER                                            :: b, flag, info
      TYPE(cell_manifold_type)                           :: m_probe

      DO b = 0, n_constraint_bits - 1
         flag = ISHFT(1, b)
         IF (IAND(candidate_bits, flag) /= 0 .AND. IAND(m%mask, flag) == 0) THEN
            CALL cell_manifold_init(m_probe, m%q0, IOR(m%mask, flag), info, &
                                    lin_w=m%lin_w(:, 1:m%n_lin_rows), &
                                    eq_pair=m%eq_pair(:, 1:m%n_eq_rows))
            IF (info == info_ok) THEN
               m = m_probe
            ELSE
               CALL resolve_degenerate_constraint(m, info, source_label, flag, do_warn)
            END IF
         END IF
      END DO

   END SUBROUTINE try_add_mask_bits

! **************************************************************************************************
!> \brief Resolve or reject a candidate constraint row that failed LICQ at q0.
!> \param m active constraint state. This routine updates m when it pins a row.
!> \param info failure code from cell_manifold_init
!> \param source_label constraint source label for warning and abort messages
!> \param mask candidate comp_* or gfix_* bit, or 0 for a lin_w or eq_pair test
!> \param do_warn .TRUE. to emit a warning when this routine pins a row
!> \param lin_w optional 6 x 1 column of candidate linear weights
!> \param eq_pair optional 2 x 1 index pair for a candidate equality
! **************************************************************************************************
   SUBROUTINE resolve_degenerate_constraint(m, info, source_label, mask, do_warn, lin_w, eq_pair)
      TYPE(cell_manifold_type), INTENT(INOUT)            :: m
      INTEGER, INTENT(IN)                                :: info
      CHARACTER(LEN=*), INTENT(IN)                       :: source_label
      INTEGER, INTENT(IN)                                :: mask
      LOGICAL, INTENT(IN)                                :: do_warn
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: lin_w(:, :)
      INTEGER, INTENT(IN), OPTIONAL                      :: eq_pair(:, :)

      INTEGER                                            :: old_mask, old_n_lin, pinned_bits
      LOGICAL                                            :: is_pinned

      IF (info == info_licq) THEN
         old_mask = m%mask
         old_n_lin = m%n_lin_rows
         CALL cell_manifold_pin_semidefinite(m, mask, is_pinned, lin_w=lin_w, eq_pair=eq_pair)
         IF (is_pinned) THEN
            IF (do_warn) THEN
               pinned_bits = IAND(IAND(NOT(old_mask), m%mask), mask_all)
               IF (pinned_bits /= 0 .AND. m%n_lin_rows == old_n_lin) THEN
                  CALL cp_warn(__LOCATION__, &
                               "CELL_OPT: "//TRIM(source_label)//" is locally degenerate at "// &
                               "the initial cell. CELL_OPT pins implied cell component(s): "// &
                               TRIM(ADJUSTL(cell_opt_constraint_label(pinned_bits)))//".")
               ELSE
                  CALL cp_warn(__LOCATION__, &
                               "CELL_OPT: "//TRIM(source_label)//" is locally degenerate at "// &
                               "the initial cell. CELL_OPT pins implied tangent direction(s).")
               END IF
            END IF
            RETURN
         END IF
         IF (cell_manifold_is_redundant(m, mask, lin_w=lin_w, eq_pair=eq_pair)) RETURN
      END IF

      CALL cp_abort(__LOCATION__, &
                    "CELL_OPT: "//TRIM(source_label)//" is locally degenerate at "// &
                    "the initial cell ("// &
                    TRIM(cell_manifold_status_label(info))//").")

   END SUBROUTINE resolve_degenerate_constraint

! **************************************************************************************************
!> \brief ...
!> \param cell_env ...
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   SUBROUTINE cell_opt_env_release(cell_env)
      TYPE(cell_opt_env_type), INTENT(INOUT)             :: cell_env

      CALL cell_release(cell_env%opt_cell)

   END SUBROUTINE cell_opt_env_release

END MODULE cell_opt_types
