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
   USE cell_opt_utils,                  ONLY: read_external_press_tensor
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
   USE input_constants,                 ONLY: fix_none,&
                                              fix_x,&
                                              fix_xy,&
                                              fix_xz,&
                                              fix_y,&
                                              fix_yz,&
                                              fix_z
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
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

      CHARACTER(LEN=4)                                   :: label
      INTEGER                                            :: output_unit, required_constraint_id
      LOGICAL                                            :: bad_constraint, constraint_explicit, &
                                                            valid_constraint
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
                                i_val=cell_env%constraint_id, explicit=constraint_explicit)
      IF (COUNT(cell%perd /= 0) == 2) THEN
         required_constraint_id = fix_none
         IF (cell%perd(1) == 0) required_constraint_id = fix_x
         IF (cell%perd(2) == 0) required_constraint_id = fix_y
         IF (cell%perd(3) == 0) required_constraint_id = fix_z

         SELECT CASE (required_constraint_id)
         CASE (fix_x)
            valid_constraint = cell_env%constraint_id == fix_x .OR. &
                               cell_env%constraint_id == fix_xy .OR. &
                               cell_env%constraint_id == fix_xz
         CASE (fix_y)
            valid_constraint = cell_env%constraint_id == fix_y .OR. &
                               cell_env%constraint_id == fix_xy .OR. &
                               cell_env%constraint_id == fix_yz
         CASE (fix_z)
            valid_constraint = cell_env%constraint_id == fix_z .OR. &
                               cell_env%constraint_id == fix_xz .OR. &
                               cell_env%constraint_id == fix_yz
         CASE DEFAULT
            valid_constraint = .FALSE.
         END SELECT

         IF (cell_env%constraint_id == fix_none .AND. .NOT. constraint_explicit) THEN
            cell_env%constraint_id = required_constraint_id
            CALL cp_warn(__LOCATION__, &
                         "2D CELL_OPT: constraining non-periodic cell direction.")
         ELSE IF (.NOT. valid_constraint) THEN
            CALL cp_abort(__LOCATION__, &
                          "2D CELL_OPT needs constrained non-periodic cell direction.")
         END IF
      END IF

      ! Reject CONSTRAINT combinations that would silently corrupt the
      ! KEEP_ANGLES projection or the KEEP_SYMMETRY symmetrization
      IF (cell_env%keep_angles .AND. .NOT. cell_env%opt_cell%orthorhombic) THEN
         ! KEEP_ANGLES makes the b and c gradient components parallel to the
         ! current b and c cell vectors. A CONSTRAINT fixing only part of the
         ! b or c vector destroys this parallelism, while constraints fixing
         ! a vector completely (fix_z, fix_yz) leave it intact.
         IF (cell_env%constraint_id == fix_x .OR. cell_env%constraint_id == fix_y .OR. &
             cell_env%constraint_id == fix_xy .OR. cell_env%constraint_id == fix_xz) THEN
            CALL cp_abort(__LOCATION__, &
                          "CELL_OPT: KEEP_ANGLES with a CONSTRAINT that fixes "// &
                          "only part of the b or c cell vector would destroy "// &
                          "the fixed angles.")
         END IF
      END IF

      IF (cell_env%keep_symmetry .AND. cell_env%constraint_id /= fix_none) THEN
         ! A CONSTRAINT is compatible with KEEP_SYMMETRY only if it leaves
         ! every set of cell components coupled by the symmetry untouched or
         ! fixes such a set completely; a partial fixation breaks the
         ! symmetry group after the gradient symmetrization.
         bad_constraint = .FALSE.
         SELECT CASE (cell_env%opt_cell%symmetry_id)
         CASE (cell_sym_cubic, cell_sym_rhombohedral)
            ! All three (cubic), respectively even all six (rhombohedral),
            ! cell components are coupled: any partial constraint breaks it
            bad_constraint = .TRUE.
         CASE (cell_sym_tetragonal_ab, &
               cell_sym_hexagonal_gamma_60, &
               cell_sym_hexagonal_gamma_120, &
               cell_sym_monoclinic_gamma_ab)
            ! The in-plane (a,b) components are coupled
            bad_constraint = cell_env%constraint_id == fix_x .OR. &
                             cell_env%constraint_id == fix_y .OR. &
                             cell_env%constraint_id == fix_xz .OR. &
                             cell_env%constraint_id == fix_yz
         CASE (cell_sym_tetragonal_ac)
            ! The (a,c) components are coupled
            bad_constraint = cell_env%constraint_id == fix_x .OR. &
                             cell_env%constraint_id == fix_z .OR. &
                             cell_env%constraint_id == fix_xy .OR. &
                             cell_env%constraint_id == fix_yz
         CASE (cell_sym_tetragonal_bc)
            ! The (b,c) components are coupled
            bad_constraint = cell_env%constraint_id == fix_y .OR. &
                             cell_env%constraint_id == fix_z .OR. &
                             cell_env%constraint_id == fix_xy .OR. &
                             cell_env%constraint_id == fix_xz
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
         SELECT CASE (cell_env%constraint_id)
         CASE (fix_x)
            label = "   X"
         CASE (fix_y)
            label = "   Y"
         CASE (fix_z)
            label = "   Z"
         CASE (fix_xy)
            label = "  XY"
         CASE (fix_xz)
            label = "  XZ"
         CASE (fix_yz)
            label = "  YZ"
         CASE (fix_none)
            label = "NONE"
         END SELECT
         WRITE (UNIT=output_unit, FMT="(T2,A,T77,A4)") &
            "CELL_OPT| Constraint: ", label
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
