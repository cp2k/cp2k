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
   USE cell_opt_utils,                  ONLY: parse_cell_constraint,&
                                              read_external_press_tensor
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_release,&
                                              cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
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
      CHARACTER(LEN=32)                          :: constraint_label = "NONE"
      INTEGER, DIMENSION(6)                      :: cell_relax = 1
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

      INTEGER                                            :: output_unit
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

      ! Parse and check cell constraints
      CALL parse_cell_constraint(geo_section, cell, cell_env%cell_relax, cell_env%constraint_label)

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
         WRITE (UNIT=output_unit, FMT="(T2,A,T49,A32)") &
            "CELL_OPT| Constraint: ", ADJUSTR(cell_env%constraint_label)
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
