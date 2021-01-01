!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Control parameters for optimizers that work with CDFT constraints
!> \par   History
!>                 separated from scf_control_types [03.2018]
!> \author Nico Holmberg [03.2018]
! **************************************************************************************************
MODULE qs_cdft_opt_types

   USE input_constants,                 ONLY: &
        broyden_type_1, broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, &
        broyden_type_2, broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, &
        outer_scf_optimizer_broyden, outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_opt_types'

   ! Public data types

   PUBLIC :: cdft_opt_type

   ! Public subroutines

   PUBLIC :: cdft_opt_type_create, &
             cdft_opt_type_release, &
             cdft_opt_type_read, &
             cdft_opt_type_write, &
             cdft_opt_type_copy

! **************************************************************************************************
!> \brief contains the parameters needed by CDFT specific optimizers
!> \param build_jacobian logical which determines if the inverse Jacobian should be computed
!> \param jacobian_step the step size for calculating the finite difference Jacobian
!> \param newton_step the step size used by the Newton optimizer with values between 0 and 1
!> \param newton_step_save permanent copy of the above
!> \param jacobian_type the finite difference scheme to compute the Jacobian
!> \param broyden_type the variant of Broyden's method to use
!> \param jacobian_freq control parameters defining how often the Jacobian is built
!> \param ijacobian counter to track how many SCF iterations/energy evaluations have passed since
!>                  the last Jacobian rebuild
!> \param broyden_update logical which determines if a Broyden update is needed
!> \param max_ls the maximum number of backtracking line search steps to perform
!> \param continue_ls continue line search until max steps are reached or until the gradient
!>                    no longer decreases
!> \param factor_ls line search parameter used in generating a new step size
!> \par History
!>                 created [03.2018]
!> \author Nico Holmberg [03.2018]
! **************************************************************************************************

   TYPE cdft_opt_type
      LOGICAL                              :: build_jacobian
      LOGICAL                              :: broyden_update
      LOGICAL                              :: continue_ls
      LOGICAL                              :: jacobian_restart
      REAL(KIND=dp)                        :: newton_step
      REAL(KIND=dp)                        :: newton_step_save
      REAL(KIND=dp)                        :: factor_ls
      REAL(KIND=dp), DIMENSION(:), &
         ALLOCATABLE                       :: jacobian_step
      REAL(KIND=dp), DIMENSION(:), &
         POINTER                           :: jacobian_vector
      INTEGER                              :: jacobian_type
      INTEGER                              :: broyden_type
      INTEGER                              :: jacobian_freq(2)
      INTEGER                              :: ijacobian(2)
      INTEGER                              :: max_ls
   END TYPE cdft_opt_type

CONTAINS

! **************************************************************************************************
!> \brief allocates and initializes the CDFT optimizer control object with default values
!> \param cdft_opt_control the object to initialize
!> \par History
!>      03.2018 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_opt_type_create(cdft_opt_control)

      TYPE(cdft_opt_type), POINTER                       :: cdft_opt_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cdft_opt_type_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(cdft_opt_control))
      ALLOCATE (cdft_opt_control)

      ! Load the default values

      cdft_opt_control%jacobian_type = -1
      cdft_opt_control%broyden_type = -1
      cdft_opt_control%jacobian_freq(:) = 1
      cdft_opt_control%newton_step = 1.0_dp
      cdft_opt_control%newton_step_save = 1.0_dp
      cdft_opt_control%factor_ls = 0.5_dp
      cdft_opt_control%ijacobian(:) = 0
      cdft_opt_control%max_ls = 0
      cdft_opt_control%build_jacobian = .FALSE.
      cdft_opt_control%broyden_update = .FALSE.
      cdft_opt_control%continue_ls = .FALSE.
      cdft_opt_control%jacobian_restart = .FALSE.
      NULLIFY (cdft_opt_control%jacobian_vector)

      CALL timestop(handle)

   END SUBROUTINE cdft_opt_type_create

! **************************************************************************************************
!> \brief releases the CDFT optimizer control object
!> \param cdft_opt_control the object to release
!> \par History
!>      03.2018 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_opt_type_release(cdft_opt_control)

      TYPE(cdft_opt_type), POINTER                       :: cdft_opt_control

      IF (ASSOCIATED(cdft_opt_control)) THEN
         IF (ASSOCIATED(cdft_opt_control%jacobian_vector)) &
            DEALLOCATE (cdft_opt_control%jacobian_vector)
         IF (ALLOCATED(cdft_opt_control%jacobian_step)) &
            DEALLOCATE (cdft_opt_control%jacobian_step)

         DEALLOCATE (cdft_opt_control)
      END IF

      NULLIFY (cdft_opt_control)

   END SUBROUTINE cdft_opt_type_release

! **************************************************************************************************
!> \brief reads the parameters of the CDFT optimizer type
!> \param cdft_opt_control the object that will contain the values read
!> \param inp_section the input section that contains the values that are read
!> \par History
!>      03.2018 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_opt_type_read(cdft_opt_control, inp_section)

      TYPE(cdft_opt_type), POINTER                       :: cdft_opt_control
      TYPE(section_vals_type), POINTER                   :: inp_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cdft_opt_type_read'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: exists
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rtmplist
      TYPE(section_vals_type), POINTER                   :: cdft_opt_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(cdft_opt_control))
      cdft_opt_section => section_vals_get_subs_vals(inp_section, "CDFT_OPT")

      CALL section_vals_val_get(cdft_opt_section, "MAX_LS", &
                                i_val=cdft_opt_control%max_ls)
      CALL section_vals_val_get(cdft_opt_section, "JACOBIAN_TYPE", &
                                i_val=cdft_opt_control%jacobian_type)
      CALL section_vals_val_get(cdft_opt_section, "JACOBIAN_STEP", &
                                r_vals=rtmplist)
      ALLOCATE (cdft_opt_control%jacobian_step(SIZE(rtmplist)))
      cdft_opt_control%jacobian_step(:) = rtmplist
      CALL section_vals_val_get(cdft_opt_section, "BROYDEN_TYPE", &
                                i_val=cdft_opt_control%broyden_type)
      CALL section_vals_val_get(cdft_opt_section, "CONTINUE_LS", &
                                l_val=cdft_opt_control%continue_ls)
      CALL section_vals_val_get(cdft_opt_section, "FACTOR_LS", &
                                r_val=cdft_opt_control%factor_ls)
      IF (cdft_opt_control%factor_ls .LE. 0.0_dp .OR. &
          cdft_opt_control%factor_ls .GE. 1.0_dp) &
         CALL cp_abort(__LOCATION__, &
                       "Keyword FACTOR_LS must be between 0.0 and 1.0.")
      CALL section_vals_val_get(cdft_opt_section, "JACOBIAN_FREQ", explicit=exists)
      IF (exists) THEN
         CALL section_vals_val_get(cdft_opt_section, "JACOBIAN_FREQ", &
                                   i_vals=tmplist)
         IF (SIZE(tmplist) /= 2) &
            CALL cp_abort(__LOCATION__, &
                          "Keyword JACOBIAN_FREQ takes exactly two input values.")
         IF (ANY(tmplist .LT. 0)) &
            CALL cp_abort(__LOCATION__, &
                          "Keyword JACOBIAN_FREQ takes only positive values.")
         IF (ALL(tmplist .EQ. 0)) &
            CALL cp_abort(__LOCATION__, &
                          "Both values to keyword JACOBIAN_FREQ cannot be zero.")
         cdft_opt_control%jacobian_freq(:) = tmplist(1:2)
      END IF
      CALL section_vals_val_get(cdft_opt_section, "JACOBIAN_RESTART", &
                                l_val=cdft_opt_control%jacobian_restart)
      IF (cdft_opt_control%jacobian_restart) THEN
         CALL section_vals_val_get(cdft_opt_section, "JACOBIAN_VECTOR", &
                                   r_vals=rtmplist)
         ALLOCATE (cdft_opt_control%jacobian_vector(SIZE(rtmplist)))
         cdft_opt_control%jacobian_vector = rtmplist
      END IF

      CALL timestop(handle)

   END SUBROUTINE cdft_opt_type_read

! **************************************************************************************************
!> \brief writes information about the CDFT optimizer object
!> \param cdft_opt_control the CDFT optimizer object
!> \param optimizer the type of optimizer to use
!> \param output_unit the output unit handle
!> \par History
!>      03.2018 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_opt_type_write(cdft_opt_control, optimizer, output_unit)
      TYPE(cdft_opt_type), POINTER                       :: cdft_opt_control
      INTEGER                                            :: optimizer, output_unit

      CPASSERT(ASSOCIATED(cdft_opt_control))

      SELECT CASE (optimizer)
      CASE DEFAULT
         ! Do nothing
      CASE (outer_scf_optimizer_broyden)
         WRITE (output_unit, '(T3,A)') "Optimization with Broyden's method"
         SELECT CASE (cdft_opt_control%broyden_type)
         CASE (broyden_type_1)
            WRITE (output_unit, '(A)') "                  variant : 1st method"
         CASE (broyden_type_1_explicit)
            WRITE (output_unit, '(A)') "                  variant : 1st method with explicit initial Jacobian"
         CASE (broyden_type_1_ls)
            WRITE (output_unit, '(A)') "                  variant : 1st method with backtracking line search"
         CASE (broyden_type_1_explicit_ls)
            WRITE (output_unit, '(A)') &
               "                  variant : 1st method with explicit initial Jacobian"
            WRITE (output_unit, '(A)') &
               "                            and backtracking line search"
         CASE (broyden_type_2)
            WRITE (output_unit, '(A)') "                  variant : 2nd method"
         CASE (broyden_type_2_explicit)
            WRITE (output_unit, '(A)') "                  variant : 2nd method with explicit initial Jacobian"
         CASE (broyden_type_2_ls)
            WRITE (output_unit, '(A)') "                  variant : 2nd method with backtracking line search"
         CASE (broyden_type_2_explicit_ls)
            WRITE (output_unit, '(A)') &
               "                  variant : 2nd method with explicit initial Jacobian"
            WRITE (output_unit, '(A)') &
               "                            and backtracking line search"
         END SELECT
      CASE (outer_scf_optimizer_newton)
         WRITE (output_unit, '(T3,A)') "Optimization with Newton's method"
      CASE (outer_scf_optimizer_newton_ls)
         WRITE (output_unit, '(T3,A)') "Optimization with Newton's method using backtracking line search"
      END SELECT
      SELECT CASE (optimizer)
      CASE DEFAULT
         ! Do nothing
      CASE (outer_scf_optimizer_broyden, outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls)
         IF (cdft_opt_control%jacobian_freq(2) > 0) THEN
            WRITE (output_unit, '(T6,A,I4,A)') &
               "The Jacobian is restarted every ", cdft_opt_control%jacobian_freq(2), " energy evaluation"
            IF (cdft_opt_control%jacobian_freq(1) > 0) &
               WRITE (output_unit, '(T29,A,I4,A)') &
               "or every ", cdft_opt_control%jacobian_freq(1), " CDFT SCF iteration"
         ELSE
            WRITE (output_unit, '(T6,A,I4,A)') &
               "The Jacobian is restarted every ", cdft_opt_control%jacobian_freq(1), " CDFT SCF iteration"
         END IF
         WRITE (output_unit, '(T3,A,F8.4)') &
            "Optimizer step size: ", cdft_opt_control%newton_step_save
      END SELECT

   END SUBROUTINE cdft_opt_type_write

! **************************************************************************************************
!> \brief copies settings between two CDFT optimizer control objects retaining both
!> \param new the object where to copy the settings
!> \param old the object from where to copy the settings
!> \par History
!>      03.2018 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE cdft_opt_type_copy(new, old)

      TYPE(cdft_opt_type), POINTER                       :: new, old

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cdft_opt_type_copy'

      INTEGER                                            :: handle

      ! Do nothing if cdft_opt_type is not allocated
      ! this happens if CDFT is performed with an optimizer other than Broyden/Newton
      IF (.NOT. ASSOCIATED(old)) RETURN

      CALL timeset(routineN, handle)

      IF (.NOT. ASSOCIATED(new)) CALL cdft_opt_type_create(new)
      new%max_ls = old%max_ls
      new%continue_ls = old%continue_ls
      new%factor_ls = old%factor_ls
      new%jacobian_type = old%jacobian_type
      new%jacobian_freq(:) = old%jacobian_freq(:)
      new%newton_step = old%newton_step
      new%newton_step_save = old%newton_step_save
      new%ijacobian(:) = old%ijacobian(:)
      new%build_jacobian = old%build_jacobian
      new%broyden_type = old%broyden_type
      new%broyden_update = old%broyden_update
      IF (ALLOCATED(new%jacobian_step)) DEALLOCATE (new%jacobian_step)
      ALLOCATE (new%jacobian_step(SIZE(old%jacobian_step)))
      new%jacobian_step(:) = old%jacobian_step
      IF (old%jacobian_restart) THEN
         ! Transfer restart vector for inverse Jacobian matrix
         ! (qs_calculate_inverse_jacobian handles deallocation of transferred vector)
         new%jacobian_restart = .TRUE.
         ALLOCATE (new%jacobian_vector(SIZE(old%jacobian_vector)))
         new%jacobian_vector = old%jacobian_vector
         DEALLOCATE (old%jacobian_vector)
         old%jacobian_restart = .FALSE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE cdft_opt_type_copy

END MODULE qs_cdft_opt_types
