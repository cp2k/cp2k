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

   USE cell_methods,                    ONLY: cell_create
   USE cell_opt_utils,                  ONLY: read_external_press_tensor
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_release,&
                                              cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_constants,                 ONLY: &
        canonicalize_cell_auto, canonicalize_cell_false, fix_none, fix_x, fix_xy, fix_xz, fix_y, &
        fix_yz, fix_z
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE particle_list_types,             ONLY: particle_list_type
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
      REAL(KIND=dp), DIMENSION(3, 3)             :: rot_matrix = 0.0_dp
      TYPE(cell_type), POINTER                   :: ref_cell => NULL()
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
      INTEGER                                            :: canonicalize_mode, ip, output_unit, &
                                                            required_constraint_id
      LOGICAL                                            :: canonical_cell, constraint_explicit, &
                                                            keep_space_group, valid_constraint
      REAL(KIND=dp), DIMENSION(3)                        :: r
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(section_vals_type), POINTER                   :: cell_section

      NULLIFY (cell_env%ref_cell, cell, cell_section, subsys, particles)
      CALL force_env_get(force_env, cell=cell, subsys=subsys)
      CALL cell_create(cell_env%ref_cell)
      CALL cell_clone(cell, cell_env%ref_cell, tag="REF_CELL_OPT")
      CALL section_vals_val_get(geo_section, "KEEP_VOLUME", l_val=cell_env%keep_volume)
      CALL section_vals_val_get(geo_section, "KEEP_ANGLES", l_val=cell_env%keep_angles)
      CALL section_vals_val_get(geo_section, "KEEP_SYMMETRY", l_val=cell_env%keep_symmetry)
      CALL section_vals_val_get(geo_section, "PRESSURE_TOLERANCE", r_val=cell_env%pres_tol)
      CALL section_vals_val_get(geo_section, "CONSTRAINT", &
                                i_val=cell_env%constraint_id, explicit=constraint_explicit)
      CALL section_vals_val_get(geo_section, "KEEP_SPACE_GROUP", l_val=keep_space_group)

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

      cell_section => section_vals_get_subs_vals(force_env%force_env_section, "SUBSYS%CELL")
      CALL section_vals_val_get(cell_section, "CANONICALIZE", i_val=canonicalize_mode)
      canonical_cell = MAXVAL(ABS([cell%hmat(2, 1), cell%hmat(3, 1), cell%hmat(3, 2)])) <= &
         1.0E-12_dp*MAX(1.0_dp, MAXVAL(ABS(cell%hmat)))

      IF (.NOT. canonical_cell) THEN
         IF (canonicalize_mode == canonicalize_cell_auto) THEN
            CALL cp_abort(__LOCATION__, &
                          "CELL_OPT requires a canonical cell matrix because it optimizes only "// &
                          "six independent matrix components. Set "// &
                          "FORCE_EVAL%SUBSYS%CELL%CANONICALIZE TRUE or specify the cell using "// &
                          "ABC and ALPHA_BETA_GAMMA. "// &
                          "To retain the legacy unconstrained behavior, explicitly set "// &
                          "CANONICALIZE FALSE.")
         ELSE IF (canonicalize_mode == canonicalize_cell_false) THEN
            IF (keep_space_group .OR. cell_env%keep_angles .OR. cell_env%keep_symmetry .OR. &
                cell_env%keep_volume .OR. cell_env%constraint_id /= fix_none) THEN
               CALL cp_abort(__LOCATION__, &
                             "A noncanonical cell cannot be combined with CELL_OPT constraints, "// &
                             "KEEP_ANGLES, KEEP_SPACE_GROUP, KEEP_SYMMETRY, or KEEP_VOLUME. Set "// &
                             "FORCE_EVAL%SUBSYS%CELL%CANONICALIZE TRUE.")
            END IF
            CALL cp_warn(__LOCATION__, &
                         "CELL_OPT with CANONICALIZE FALSE retains three general-orientation "// &
                         "cell components while optimizing the other six. This legacy expert "// &
                         "mode is coordinate-frame dependent.")
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "CELL_OPT received a noncanonical cell despite requested canonicalization.")
         END IF
      END IF

      ! Compute the rotation matrix that give the cell vectors in the "canonical" orientation
      cell_env%rot_matrix = MATMUL(cell_env%ref_cell%hmat, cell%h_inv)

      ! Get the external pressure
      CALL read_external_press_tensor(geo_section, cell, cell_env%pres_ext, cell_env%mtrx, &
                                      cell_env%rot_matrix)

      ! Rotate particles accordingly
      CALL cp_subsys_get(subsys, particles=particles)
      DO ip = 1, particles%n_els
         r = MATMUL(TRANSPOSE(cell_env%rot_matrix), particles%els(ip)%r)
         particles%els(ip)%r = r
      END DO

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

      CALL cell_release(cell_env%ref_cell)

   END SUBROUTINE cell_opt_env_release

END MODULE cell_opt_types
