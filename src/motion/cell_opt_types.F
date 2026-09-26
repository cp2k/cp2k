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

   USE cell_methods,                    ONLY: canonicalize_cell_matrix,&
                                              cell_create
   USE cell_opt_utils,                  ONLY: &
        cell_opt_constraint_label, fixes_partially, mask_ab_plane, mask_all, mask_b, mask_c, &
        mask_diag, mask_diag_ab, mask_diag_ac, mask_diag_bc, parse_cell_constraint, &
        read_external_press_tensor
   USE cell_types,                      ONLY: &
        cell_clone, cell_release, cell_sym_cubic, cell_sym_hexagonal_gamma_120, &
        cell_sym_hexagonal_gamma_60, cell_sym_monoclinic_gamma_ab, cell_sym_rhombohedral, &
        cell_sym_tetragonal_ab, cell_sym_tetragonal_ac, cell_sym_tetragonal_bc, cell_type
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
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_types'

   PUBLIC :: cell_opt_env_type, &
             cell_opt_env_create, &
             cell_opt_env_release

! **************************************************************************************************
!> \brief Type containing all informations abour the simulation cell optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
   TYPE cell_opt_env_type
      ! Simulation cell optimization parameters
      ! constraint_id is a mask of fixed cell components (comp_* / mask_*):
      ! a set bit means fixed, 0 (= fix_none) leaves all components relaxable
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

      CHARACTER(LEN=default_string_length)               :: constraint_string, label
      INTEGER                                            :: output_unit
      LOGICAL                                            :: bad_constraint, constraint_explicit
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (cell_env%opt_cell, cell)
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
                                c_val=constraint_string, explicit=constraint_explicit)
      CALL parse_cell_constraint(constraint_string, cell, constraint_explicit, &
                                 cell_env%constraint_id)
      ! Canonicalize the echoed value so restart files record the effective
      ! constraint instead of the user's spelling (case and token order).
      ! Only for explicitly set values: an unset keyword must stay unparsed so
      ! the restart echo keeps omitting the line, as it always has.
      IF (constraint_explicit) THEN
         CALL section_vals_val_set(geo_section, "CONSTRAINT", &
                                   c_val=TRIM(ADJUSTL( &
                                              cell_opt_constraint_label(cell_env%constraint_id))))
      END IF

      ! Reject CONSTRAINT combinations that would silently corrupt the
      ! KEEP_ANGLES projection or the KEEP_SYMMETRY symmetrization.
      ! Slot order of apply_cell_constraints: AX, BX, BY, CX, CY, CZ.
      ! KEEP_ANGLES projects the gradient slots {BX,BY} onto b and {CX,CY,CZ} onto c.
      ! Gradient slots that a legacy value zeroes:
      !    x  = {AX,BX,CX}          B partial,  C partial   -> breaks parallelism
      !    y  = {BX,BY,CY}          B complete, C partial   -> breaks parallelism
      !    z  = {CX,CY,CZ}          B untouched, C complete -> safe
      !    xy = {AX..CY}            B complete, C partial   -> breaks parallelism
      !    xz = {AX,BX,CX,CY,CZ}    B partial,  C complete  -> breaks parallelism
      !    yz = {BX..CZ}            B complete, C complete  -> safe
      ! A partial group breaks the parallelism. A complete or untouched group keeps it.
      IF (cell_env%keep_angles .AND. .NOT. cell_env%opt_cell%orthorhombic) THEN
         ! KEEP_ANGLES makes the b and c gradient components parallel to the
         ! current b and c cell vectors. A CONSTRAINT that fixes only part of
         ! the b or c vector destroys this parallelism. A constraint that
         ! fixes the vector completely leaves it intact.
         IF (fixes_partially(cell_env%constraint_id, mask_b) .OR. &
             fixes_partially(cell_env%constraint_id, mask_c)) THEN
            CALL cp_abort(__LOCATION__, &
                          "CELL_OPT: KEEP_ANGLES with a CONSTRAINT that fixes "// &
                          "only part of the b or c cell vector would destroy "// &
                          "the fixed angles.")
         END IF
      END IF

      IF (cell_env%keep_symmetry .AND. cell_env%constraint_id /= fix_none) THEN
         ! KEEP_SYMMETRY couples sets of cell components. A CONSTRAINT is
         ! compatible only if it leaves such a set untouched, or fixes the
         ! set completely. A partial fixation breaks the symmetry group
         ! after the gradient symmetrization.
         bad_constraint = .FALSE.
         SELECT CASE (cell_env%opt_cell%symmetry_id)
         CASE (cell_sym_cubic)
            ! All three cell vector lengths are coupled
            bad_constraint = fixes_partially(cell_env%constraint_id, mask_diag)
         CASE (cell_sym_rhombohedral)
            ! Even all six cell components are coupled
            bad_constraint = fixes_partially(cell_env%constraint_id, mask_all)
         CASE (cell_sym_hexagonal_gamma_60, &
               cell_sym_hexagonal_gamma_120, &
               cell_sym_monoclinic_gamma_ab)
            ! The in-plane (a,b) components are coupled
            bad_constraint = fixes_partially(cell_env%constraint_id, mask_ab_plane)
         CASE (cell_sym_tetragonal_ab)
            ! The lengths of a and b are coupled
            bad_constraint = fixes_partially(cell_env%constraint_id, mask_diag_ab)
         CASE (cell_sym_tetragonal_ac)
            ! The lengths of a and c are coupled
            bad_constraint = fixes_partially(cell_env%constraint_id, mask_diag_ac)
         CASE (cell_sym_tetragonal_bc)
            ! The lengths of b and c are coupled
            bad_constraint = fixes_partially(cell_env%constraint_id, mask_diag_bc)
         CASE DEFAULT
            ! orthorhombic, monoclinic, triclinic and none have no coupling
         END SELECT
         IF (bad_constraint) THEN
            CALL cp_abort(__LOCATION__, &
                          "CELL_OPT: CONSTRAINT fixes only part of the cell "// &
                          "components coupled by KEEP_SYMMETRY and would "// &
                          "break the cell symmetry.")
         END IF
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
         WRITE (UNIT=output_unit, FMT="(T2,A,T77,A)") &
            "CELL_OPT| Constraint: ", TRIM(label)
      END IF
      CALL cp_print_key_finished_output(output_unit, logger, geo_section, "PRINT%CELL")

   END SUBROUTINE cell_opt_env_create

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
