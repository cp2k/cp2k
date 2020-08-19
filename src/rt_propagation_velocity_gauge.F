!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to perform the RTP in the velocity gauge
! **************************************************************************************************

MODULE rt_propagation_velocity_gauge
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: one,&
                                              twopi,&
                                              zero
   USE physcon,                         ONLY: c_light_au
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_velocity_gauge'

   PUBLIC :: velocity_gauge_ks_matrix

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE velocity_gauge_ks_matrix(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'velocity_gauge_ks_matrix'

      INTEGER                                            :: handle, idir, ispin
      LOGICAL                                            :: do_hfx
      REAL(KIND=dp), DIMENSION(3)                        :: vec_pot
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_im, matrix_s, &
                                                            momentum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control)

      IF (dft_control%qs_control%gapw) THEN

         NULLIFY (cell, input, matrix_s, matrix_ks, matrix_ks_im, sab_orb)
         CALL get_qs_env(qs_env, &
                         cell=cell, &
                         input=input, &
                         matrix_s=matrix_s, &
                         matrix_ks=matrix_ks, &
                         matrix_ks_im=matrix_ks_im, &
                         sab_orb=sab_orb)

         !get vector potential
         CALL get_vector_potential(dft_control, cell, vec_pot)

         ! allocate and build matrices for linear momentum term
         NULLIFY (momentum)
         CALL dbcsr_allocate_matrix_set(momentum, 3)
         DO idir = 1, 3
            CALL dbcsr_init_p(momentum(idir)%matrix)
            CALL dbcsr_create(momentum(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(momentum(idir)%matrix, sab_orb)
            CALL dbcsr_set(momentum(idir)%matrix, zero)
         ENDDO
         CALL build_lin_mom_matrix(qs_env, momentum)

         ! set imaginary part of KS matrix to zero
         ! if not done for hfx in subroutine hfx_ks_matrix()
         hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
         CALL section_vals_get(hfx_sections, explicit=do_hfx)
         IF (.NOT. do_hfx) THEN
            DO ispin = 1, SIZE(matrix_ks_im)
               CALL dbcsr_set(matrix_ks_im(ispin)%matrix, zero)
            ENDDO
         ENDIF

         ! add linear term in vector potential to imaginary part of KS-matrix
         DO ispin = 1, SIZE(matrix_ks_im)
            DO idir = 1, 3
               CALL dbcsr_add(matrix_ks_im(ispin)%matrix, momentum(idir)%matrix, one, -2.0_dp*vec_pot(idir))
            ENDDO
         ENDDO

         CALL dbcsr_deallocate_matrix_set(momentum)

         ! add quadratic term to real part of KS matrix
         DO ispin = 1, SIZE(matrix_ks)
            DO idir = 1, 3
               CALL dbcsr_add(matrix_ks(ispin)%matrix, matrix_s(1)%matrix, one, vec_pot(idir)**2/c_light_au)
            ENDDO
         ENDDO

      ELSE
         CPABORT(" Not implemented for GPW yet!")
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE velocity_gauge_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param cell ...
!> \param vec_pot ...
! **************************************************************************************************
   PURE SUBROUTINE get_vector_potential(dft_control, cell, vec_pot)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: vec_pot

      CHARACTER(len=*), PARAMETER :: routineN = 'get_vector_potential'

      ! fieldstrength
      ! for a delta pulse: vec_pot = - c * kvec
      ! the speed of light cancels in the linear term in the Hamiltonian
      vec_pot(:) = cell%h_inv(1, :)*dft_control%rtp_control%delta_pulse_direction(1) + &
                   cell%h_inv(2, :)*dft_control%rtp_control%delta_pulse_direction(2) + &
                   cell%h_inv(3, :)*dft_control%rtp_control%delta_pulse_direction(3)
      vec_pot = vec_pot*twopi*dft_control%rtp_control%delta_pulse_scale

   END SUBROUTINE get_vector_potential

END MODULE rt_propagation_velocity_gauge
