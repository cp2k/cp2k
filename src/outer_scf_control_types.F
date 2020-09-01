!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief parameters that control the outer loop of an SCF iteration
!> \par History
!>      09.2018 created by moving outer SCF types to separate module [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
MODULE outer_scf_control_types

   USE input_constants,                 ONLY: outer_scf_optimizer_broyden,&
                                              outer_scf_optimizer_newton,&
                                              outer_scf_optimizer_newton_ls
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_cdft_opt_types,               ONLY: cdft_opt_type,&
                                              cdft_opt_type_create,&
                                              cdft_opt_type_read
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'outer_scf_control_types'
   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   ! Public data types

   PUBLIC :: outer_scf_control_type, &
             qs_outer_scf_type

   ! Public subroutines

   PUBLIC :: outer_scf_read_parameters

! **************************************************************************************************
!> \brief contains the parameters needed by a scf run
!> \param density_guess how to choose the initial density
!>        (CORE,RANDOM,RESTART,ATOMIC,FROZEN)
!> \param eps_eigval wanted error on the eigenvalues
!> \param eps_scf whanted error on the whole scf
!> \param level_shift amount of level shift
!> \param p_mix how to mix the new and old densities in non diss iterations
!> \param eps_lumos error on the lumos calculated at the end of the scf
!> \param max_iter_lumus maxumum number of iterations used to calculate
!>        the lumos at the end of the scf
!> \param max_scf max scf iterations
!> \param id_nr unique number to identify an scf control
!> \param ref_count reference count (see cp2k/doc/ReferenceCounting.html)
!> \param added_mos additional number of MOs that might be used in the SCF
!> \param step_size the optimizer step size
!> \param cdft_opt_control settings for optimizers that work only together with CDFT constraints
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************

   TYPE outer_scf_control_type
      LOGICAL       :: have_scf
      INTEGER       :: max_scf
      REAL(KIND=dp) :: eps_scf, step_size
      INTEGER       :: TYPE
      INTEGER       :: optimizer
      INTEGER       :: diis_buffer_length
      INTEGER       :: extrapolation_order
      INTEGER       :: bisect_trust_count
      TYPE(cdft_opt_type), POINTER :: cdft_opt_control
   END TYPE outer_scf_control_type

   TYPE qs_outer_scf_type
      INTEGER :: iter_count
      LOGICAL :: deallocate_jacobian
      ! these are the variable of outer loop.
      ! right now, we assume that they can be easily written as
      ! small arrays, but we might want to go the cp_fm_types
      ! at a later stage
      ! also, we just store the full iteration history
      REAL(KIND=dp), DIMENSION(:), POINTER :: energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: variables
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: gradient
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: inv_jacobian
      INTEGER, DIMENSION(:), POINTER :: count
   END TYPE qs_outer_scf_type

CONTAINS

! **************************************************************************************************
!> \brief reads the parameters of the outer_scf section into the given outer_scf
!> \param outer_scf the object that wil contain the values read
!> \param outer_scf_section the input section
!> \par History
!>      09.2018 created by separating from scf_c_read_parameters [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE outer_scf_read_parameters(outer_scf, outer_scf_section)

      TYPE(outer_scf_control_type)                       :: outer_scf
      TYPE(section_vals_type), POINTER                   :: outer_scf_section

      LOGICAL                                            :: exists

      CALL section_vals_val_get(outer_scf_section, "_SECTION_PARAMETERS_", &
                                l_val=outer_scf%have_scf)
      IF (outer_scf%have_scf) THEN
         CALL section_vals_val_get(outer_scf_section, "EPS_SCF", &
                                   r_val=outer_scf%eps_scf)
         CALL section_vals_val_get(outer_scf_section, "STEP_SIZE", &
                                   r_val=outer_scf%step_size, explicit=exists)
         CALL section_vals_val_get(outer_scf_section, "DIIS_BUFFER_LENGTH", &
                                   i_val=outer_scf%diis_buffer_length)
         CALL section_vals_val_get(outer_scf_section, "BISECT_TRUST_COUNT", &
                                   i_val=outer_scf%bisect_trust_count)
         CALL section_vals_val_get(outer_scf_section, "TYPE", &
                                   i_val=outer_scf%type)
         CALL section_vals_val_get(outer_scf_section, "MAX_SCF", &
                                   i_val=outer_scf%max_scf)
         CALL section_vals_val_get(outer_scf_section, "EXTRAPOLATION_ORDER", &
                                   i_val=outer_scf%extrapolation_order)
         CALL section_vals_val_get(outer_scf_section, "OPTIMIZER", &
                                   i_val=outer_scf%optimizer)
         ! Optimizer specific initializations
         SELECT CASE (outer_scf%optimizer)
         CASE DEFAULT
            ! Do nothing
         CASE (outer_scf_optimizer_broyden, outer_scf_optimizer_newton, &
               outer_scf_optimizer_newton_ls)
            ! CDFT optimizer -> read CDFT_OPT section
            CALL cdft_opt_type_create(outer_scf%cdft_opt_control)
            CALL cdft_opt_type_read(outer_scf%cdft_opt_control, &
                                    outer_scf_section)
            IF (exists) THEN
               outer_scf%cdft_opt_control%newton_step = ABS(outer_scf%step_size)
               ! Permanent copy needed in case line search is performed
               outer_scf%cdft_opt_control%newton_step_save = &
                  outer_scf%cdft_opt_control%newton_step
            END IF
         END SELECT
      END IF

   END SUBROUTINE outer_scf_read_parameters

END MODULE outer_scf_control_types
