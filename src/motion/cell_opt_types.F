!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains type used for a Simulation Cell Optimization
!> \par History
!>      none
!> \author Teodoro Laino - created [tlaino] - 03.2008 - Zurich University
! **************************************************************************************************
MODULE cell_opt_types
   USE cell_opt_utils,                  ONLY: get_ut_cell_matrix,&
                                              read_external_press_tensor
   USE cell_types,                      ONLY: cell_clone,&
                                              cell_create,&
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
   USE particle_list_types,             ONLY: particle_list_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_opt_types'
   INTEGER, PRIVATE, SAVE :: last_cell_opt_env_id = 0

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
      INTEGER                                    :: constraint_id, id_nr, ref_count
      LOGICAL                                    :: keep_angles, &
                                                    keep_symmetry
      REAL(KIND=dp)                              :: pres_ext, pres_int, pres_tol, pres_constr
      REAL(KIND=dp), DIMENSION(3, 3)             :: mtrx
      REAL(KIND=dp), DIMENSION(3, 3)             :: rot_matrix
      TYPE(cell_type), POINTER                   :: ref_cell
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
      TYPE(cell_opt_env_type), POINTER                   :: cell_env
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: geo_section

      CHARACTER(LEN=4)                                   :: label
      INTEGER                                            :: ip, output_unit
      REAL(KIND=dp), DIMENSION(3)                        :: r
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(particle_list_type), POINTER                  :: particles

      CPASSERT(.NOT. ASSOCIATED(cell_env))
      ALLOCATE (cell_env)
      NULLIFY (cell_env%ref_cell, cell, subsys, particles)
      cell_env%ref_count = 1
      last_cell_opt_env_id = last_cell_opt_env_id + 1
      cell_env%id_nr = last_cell_opt_env_id
      CALL force_env_get(force_env, cell=cell, subsys=subsys)
      CALL cell_create(cell_env%ref_cell)
      CALL cell_clone(cell, cell_env%ref_cell)
      CALL section_vals_val_get(geo_section, "KEEP_ANGLES", l_val=cell_env%keep_angles)
      CALL section_vals_val_get(geo_section, "KEEP_SYMMETRY", l_val=cell_env%keep_symmetry)
      CALL section_vals_val_get(geo_section, "PRESSURE_TOLERANCE", r_val=cell_env%pres_tol)
      CALL section_vals_val_get(geo_section, "CONSTRAINT", i_val=cell_env%constraint_id)

      ! First let's rotate the cell vectors in order to have an upper triangular matrix.
      CALL get_ut_cell_matrix(cell)

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
      TYPE(cell_opt_env_type), POINTER                   :: cell_env

      IF (ASSOCIATED(cell_env)) THEN
         CPASSERT(cell_env%ref_count > 0)
         cell_env%ref_count = cell_env%ref_count - 1
         IF (cell_env%ref_count == 0) THEN
            CALL cell_release(cell_env%ref_cell)
            DEALLOCATE (cell_env)
         END IF
      END IF
   END SUBROUTINE cell_opt_env_release

END MODULE cell_opt_types
