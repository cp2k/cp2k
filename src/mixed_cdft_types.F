!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for mixed CDFT calculations
!> \par   History
!>                 Separated CDFT routines from mixed_environment_types
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
MODULE mixed_cdft_types
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: becke_control_release,&
                                              becke_restraint_type
   USE cp_log_handling,                 ONLY: cp_logger_release,&
                                              cp_logger_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_release,&
                                              pw_env_type
   USE qs_kind_types,                   ONLY: deallocate_qs_kind_set,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

! **************************************************************************************************
!> \brief Buffers for load balancing
!> \param rank indices of the processors the data in this buffer should be sent to
!> \param tag mpi tags for the messages to send
!> \param cavity the cavity to send
!> \param weight the weight to send
!> \param gradients the gradients to send
! **************************************************************************************************
   TYPE buffers
      INTEGER                                          :: rank(2), tag(2)
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                            :: cavity, weight
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)                         :: gradients
   END TYPE buffers
! **************************************************************************************************
!> \brief To build array of buffers
!> \param buffs the pointer to the buffers type
! **************************************************************************************************
   TYPE p_buffers
      TYPE(buffers), DIMENSION(:), POINTER             :: buffs
   END TYPE p_buffers
! **************************************************************************************************
!> \brief Information about load balancing
!> \param matrix_info size of the target_list array to receive and grid point bounds of the data
!> \param target_list the target_list array of the processor that sends me data
! **************************************************************************************************
   TYPE repl_info
      INTEGER, DIMENSION(:), POINTER                   :: matrix_info
      INTEGER, DIMENSION(:, :), POINTER                :: target_list
   END TYPE repl_info
! **************************************************************************************************
!> \brief Load balancing control for mixed CDFT calculation
!> \param my_source index of the processor which will send this processor data
!> \param distributed bounds that determine which grid points this processor will compute after
!>                    applying load balancing (is_special = .FALSE.)
!> \param my_dest_repl the dest_list arrays of all processors which send additional work to this
!>                     processor (indices of the processors where the redistributed slices should be
!>                     returned)
!> \param dest_tags_repl tags for the send messages (is_special = .FALSE.)
!> \param more_work allow heavily overloaded processors to redistribute more_work slices
!> \param bo bounds of the data that this processor will send to other processors which tells the
!>           receivers how to rearrange the data correctly
!> \param expected_work a list of the estimated work per processor
!> \param prediction_error the difference between the estimated and actual work per processor
!> \param target_list a list of processors to send data and the size of data to send
!> \param recv_work flag that determines if this processor will receive data from others
!> \param send_work flag that determines if this processor will send data to others
!> \param recv_work_repl list of processor indices where this processor will send data during load
!>                       balancing
!> \param load_scale allow underloaded processors to accept load_scale additional work
!> \param very_overloaded value to determine which processors are heavily overloaded
!> \param cavity the cavity that this processor builds in addition to its own cavity defined
!>               on the grid points which were redistributed to this processor
!> \param weight the weight that this processor builds in addition to its own weight
!> \param gradients the gradients that this processor builds in addition to its own gradients
!> \param sendbuffer buffer to hold the data this processor will send
!> \param sendbuffer buffer to hold the data this processor will receive
!> \param recv_info additional information on the data this processor will receive
! **************************************************************************************************
   TYPE mixed_cdft_dlb_type
      INTEGER                                          :: my_source, distributed(2), &
                                                          my_dest_repl(2), dest_tags_repl(2), &
                                                          more_work
      INTEGER, DIMENSION(:), POINTER                   :: bo, expected_work, &
                                                          prediction_error
      INTEGER, DIMENSION(:, :), POINTER                 :: target_list
      LOGICAL                                          :: recv_work, send_work
      LOGICAL, DIMENSION(:), POINTER                   :: recv_work_repl
      REAL(KIND=dp)                                    :: load_scale, very_overloaded
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                              :: cavity, weight
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)                            :: gradients
      ! Should convert to TYPE(p_buffers), POINTER
      TYPE(buffers), DIMENSION(:), POINTER             :: sendbuff
      TYPE(p_buffers), DIMENSION(:), POINTER           :: recvbuff
      TYPE(repl_info), DIMENSION(:), POINTER           :: recv_info
   END TYPE mixed_cdft_dlb_type
! **************************************************************************************************
!> \brief Main mixed CDFT control type
!> \param sim_step counter to keep track of the simulation step for MD
!> \param multiplicity spin multiplicity
!> \param constraint_type determins what kind of constraint to use
!> \param combined_type for combined density+spin constraint, determines which atoms to apply the
!>                      spin constraint
!> \param source_list a list of processors which will send this processor data
!> \param dest_list a list of processors which this processor will send data to
!> \param recv_bo bounds of the data which this processor will receive (is_special = .FALSE.)
!> \param source_list_save permanent copy of source_list which might get reallocated during
!>                         load balancing
!> \param dest_list_save permanent copy of dest_list which might get reallocated during
!>                       load balancing
!> \param source_list_bo bounds of the data which this processor will receive (is_special = .TRUE.)
!> \param dest_list_bo bounds of the data this processor will send (is_special = .TRUE.)
!> \param source_bo_save permanent copy of source_list_bo
!> \param deset_bo_save permanent copy of dest_list_bo
!> \param is_pencil flag controlling which scheme to use for constraint replication
!> \param dlb flag to enable dynamic load balancing
!> \param is_special another flag controlling which scheme to use for constraint replication
!> \param first_iteration flag to mark the first iteration e.g. during MD to output information
!> \param calculate_metric flag which determines if the coupling reliablity metric should be computed
!> \param wnf_ovelap_method flag to enable the wavefunction overlap method for computing the coupling
!> \param has_unit_metric flag to determine if the basis set has unit metric
!> \param eps_rho_rspace threshold to determine when the realspace density can be considered zero
!> \param sim_dt timestep of the MD simulation
!> \param eps_svd value that controls which matrix inversion method to use
!> \param weight the constraint weight function
!> \param cavity the confinement cavity: the weight function is nonzero only within the cavity
!> \param becke_control container for becke_restraint_type
!> \param sendbuff buffer that holds the data to be replicated
!> \param blacs_env the blacs_env needed to redistribute arrays during a coupling calculation
!> \param dlb_control container for load balancing structures
!> \param qs_kind_set the qs_kind_set needed to setup a confinement cavity
!> \param pw_env the pw_env that holds the fully distributed realspace grid
! **************************************************************************************************
   TYPE mixed_cdft_type
      INTEGER                                          :: sim_step, multiplicity, &
                                                          constraint_type, combined_type
      INTEGER, POINTER, DIMENSION(:)                   :: source_list, dest_list, &
                                                          recv_bo, source_list_save, &
                                                          dest_list_save
      INTEGER, POINTER, DIMENSION(:, :)                :: source_list_bo, dest_list_bo, &
                                                          source_bo_save, dest_bo_save
      LOGICAL                                          :: is_pencil, dlb, &
                                                          is_special, first_iteration, &
                                                          calculate_metric, &
                                                          wfn_overlap_method, &
                                                          has_unit_metric, &
                                                          is_parallel
      REAL(KIND=dp)                                    :: eps_rho_rspace, sim_dt, &
                                                          eps_svd
      REAL(KIND=dp), POINTER, DIMENSION(:, :, :)       :: weight, cavity
      TYPE(becke_restraint_type), POINTER              :: becke_control
      TYPE(buffers), DIMENSION(:), POINTER             :: sendbuff
      TYPE(cp_blacs_env_type), POINTER                 :: blacs_env
      TYPE(mixed_cdft_dlb_type), POINTER               :: dlb_control
      TYPE(qs_kind_type), DIMENSION(:), &
         POINTER                                       :: qs_kind_set
      TYPE(pw_env_type), POINTER                       :: pw_env
      TYPE(cp_logger_type), POINTER                    :: sub_logger
   END TYPE mixed_cdft_type

! **************************************************************************************************
!> \brief Container for constraint settings to check consistency of force_evals
! **************************************************************************************************
   TYPE mixed_cdft_settings_type
      LOGICAL                                            :: is_spherical, &
                                                            is_odd
      INTEGER                                            :: nbecke, &
                                                            max_nkinds
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, DIMENSION(:), POINTER                     :: grid_span, &
                                                            spherical, &
                                                            odd
      INTEGER, DIMENSION(:, :), POINTER                  :: si, &
                                                            rs_dims, &
                                                            atoms, &
                                                            npts
      REAL(KIND=dp)                                      :: radius
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoff, &
                                                            rel_cutoff
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sr, &
                                                            coeffs, &
                                                            confine_bounds, &
                                                            cutoffs, &
                                                            radii
   END TYPE mixed_cdft_settings_type

! *** Public data types ***

   PUBLIC :: mixed_cdft_type, &
             mixed_cdft_settings_type

! *** Public subroutines ***

   PUBLIC :: mixed_cdft_type_create, &
             mixed_cdft_type_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_cdft_types'

CONTAINS

! **************************************************************************************************
!> \brief inits the given mixed_cdft_type
!> \param cdft_control the object to init
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_type_create(cdft_control)
      TYPE(mixed_cdft_type), POINTER                     :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_type_create', &
         routineP = moduleN//':'//routineN

      NULLIFY (cdft_control%pw_env, cdft_control%blacs_env, cdft_control%qs_kind_set)
      NULLIFY (cdft_control%dlb_control, cdft_control%dest_list_bo, cdft_control%dest_list)
      NULLIFY (cdft_control%dest_bo_save, cdft_control%dest_list_save, cdft_control%source_list)
      NULLIFY (cdft_control%source_list_save, cdft_control%source_bo_save, cdft_control%source_list_bo)
      NULLIFY (cdft_control%cavity, cdft_control%weight, cdft_control%sendbuff)
      NULLIFY (cdft_control%becke_control, cdft_control%recv_bo)
      NULLIFY (cdft_control%sub_logger)

   END SUBROUTINE mixed_cdft_type_create

! **************************************************************************************************
!> \brief releases the given mixed_cdft_type
!> \param cdft_control the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_type_release(cdft_control)
      TYPE(mixed_cdft_type), POINTER                     :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_type_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      CALL pw_env_release(cdft_control%pw_env)
      IF (ASSOCIATED(cdft_control%dest_list)) &
         DEALLOCATE (cdft_control%dest_list)
      IF (ASSOCIATED(cdft_control%dest_list_save)) &
         DEALLOCATE (cdft_control%dest_list_save)
      IF (ASSOCIATED(cdft_control%dest_list_bo)) &
         DEALLOCATE (cdft_control%dest_list_bo)
      IF (ASSOCIATED(cdft_control%dest_bo_save)) &
         DEALLOCATE (cdft_control%dest_bo_save)
      IF (ASSOCIATED(cdft_control%source_list)) &
         DEALLOCATE (cdft_control%source_list)
      IF (ASSOCIATED(cdft_control%source_list_save)) &
         DEALLOCATE (cdft_control%source_list_save)
      IF (ASSOCIATED(cdft_control%source_list_bo)) &
         DEALLOCATE (cdft_control%source_list_bo)
      IF (ASSOCIATED(cdft_control%source_bo_save)) &
         DEALLOCATE (cdft_control%source_bo_save)
      IF (ASSOCIATED(cdft_control%recv_bo)) &
         DEALLOCATE (cdft_control%recv_bo)
      IF (ASSOCIATED(cdft_control%weight)) &
         DEALLOCATE (cdft_control%weight)
      IF (ASSOCIATED(cdft_control%cavity)) &
         DEALLOCATE (cdft_control%cavity)
      IF (ASSOCIATED(cdft_control%dlb_control)) &
         CALL mixed_cdft_dlb_release(cdft_control%dlb_control)
      IF (ASSOCIATED(cdft_control%sendbuff)) THEN
         DO i = 1, SIZE(cdft_control%sendbuff)
            CALL mixed_cdft_buffers_release(cdft_control%sendbuff(i))
         END DO
         DEALLOCATE (cdft_control%sendbuff)
      END IF
      IF (ASSOCIATED(cdft_control%becke_control)) &
         CALL becke_control_release(cdft_control%becke_control)
      IF (ASSOCIATED(cdft_control%blacs_env)) &
         CALL cp_blacs_env_release(cdft_control%blacs_env)
      IF (ASSOCIATED(cdft_control%qs_kind_set)) &
         CALL deallocate_qs_kind_set(cdft_control%qs_kind_set)
      IF (ASSOCIATED(cdft_control%sub_logger)) &
         CALL cp_logger_release(cdft_control%sub_logger)
      DEALLOCATE (cdft_control)

   END SUBROUTINE mixed_cdft_type_release

! **************************************************************************************************
!> \brief releases the given load balancing control
!> \param dlb_control the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_dlb_release(dlb_control)
      TYPE(mixed_cdft_dlb_type), POINTER                 :: dlb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_dlb_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(dlb_control%recv_work_repl)) &
         DEALLOCATE (dlb_control%recv_work_repl)
      IF (ASSOCIATED(dlb_control%sendbuff)) THEN
         DO i = 1, SIZE(dlb_control%sendbuff)
            CALL mixed_cdft_buffers_release(dlb_control%sendbuff(i))
         END DO
         DEALLOCATE (dlb_control%sendbuff)
      END IF
      IF (ASSOCIATED(dlb_control%recvbuff)) THEN
         DO i = 1, SIZE(dlb_control%recvbuff)
            CALL mixed_cdft_p_buffers_release(dlb_control%recvbuff(i))
         END DO
         DEALLOCATE (dlb_control%recvbuff)
      END IF
      IF (ASSOCIATED(dlb_control%recv_info)) THEN
         DO i = 1, SIZE(dlb_control%recv_info)
            IF (ASSOCIATED(dlb_control%recv_info(i)%matrix_info)) &
               DEALLOCATE (dlb_control%recv_info(i)%matrix_info)
            IF (ASSOCIATED(dlb_control%recv_info(i)%target_list)) &
               DEALLOCATE (dlb_control%recv_info(i)%target_list)
         END DO
         DEALLOCATE (dlb_control%recv_info)
      END IF
      IF (ASSOCIATED(dlb_control%bo)) &
         DEALLOCATE (dlb_control%bo)
      IF (ASSOCIATED(dlb_control%expected_work)) &
         DEALLOCATE (dlb_control%expected_work)
      IF (ASSOCIATED(dlb_control%prediction_error)) &
         DEALLOCATE (dlb_control%prediction_error)
      IF (ASSOCIATED(dlb_control%target_list)) &
         DEALLOCATE (dlb_control%target_list)
      IF (ASSOCIATED(dlb_control%cavity)) &
         DEALLOCATE (dlb_control%cavity)
      IF (ASSOCIATED(dlb_control%weight)) &
         DEALLOCATE (dlb_control%weight)
      IF (ASSOCIATED(dlb_control%gradients)) &
         DEALLOCATE (dlb_control%gradients)
      DEALLOCATE (dlb_control)

   END SUBROUTINE mixed_cdft_dlb_release

! **************************************************************************************************
!> \brief releases the given buffers
!> \param buffer the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_buffers_release(buffer)
      TYPE(buffers)                                      :: buffer

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_buffers_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(buffer%cavity)) &
         DEALLOCATE (buffer%cavity)
      IF (ASSOCIATED(buffer%weight)) &
         DEALLOCATE (buffer%weight)
      IF (ASSOCIATED(buffer%gradients)) &
         DEALLOCATE (buffer%gradients)

   END SUBROUTINE mixed_cdft_buffers_release

! **************************************************************************************************
!> \brief releases the given pointer of buffers
!> \param p_buffer the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_p_buffers_release(p_buffer)
      TYPE(p_buffers)                                    :: p_buffer

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_p_buffers_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(p_buffer%buffs)) THEN
         DO i = 1, SIZE(p_buffer%buffs)
            CALL mixed_cdft_buffers_release(p_buffer%buffs(i))
         END DO
         DEALLOCATE (p_buffer%buffs)
      END IF

   END SUBROUTINE mixed_cdft_p_buffers_release

END MODULE mixed_cdft_types
