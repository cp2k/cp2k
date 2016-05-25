!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utility subroutine for qs energy calculation
!> \par History
!>      none
!> \author MK (29.10.2002)
! *****************************************************************************
MODULE qs_energy_matrix_w
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                              cp_dbcsr_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform
   USE kpoint_types,                    ONLY: kpoint_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: calculate_w_matrix,&
                                              calculate_w_matrix_ot
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_matrix_w'

   PUBLIC :: qs_energies_compute_w

CONTAINS

! *****************************************************************************
!> \brief Refactoring of qs_energies_scf. Moves computation of matrix_w
!>        into separate subroutine
!> \param qs_env ...
!> \param calc_forces ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_compute_w(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_compute_w', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, is, ispin, nao, nspin
      LOGICAL                                            :: do_kpoints, has_unit_metric
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s, matrix_w, &
                                                            matrix_w_mp2, mo_derivs, rho_ao
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_s_kp, matrix_w_kp
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fmwork
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fmstruct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      ! if calculate forces, time to compute the w matrix
      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)

      IF (calc_forces .AND. .NOT. has_unit_metric) THEN
         CALL get_qs_env(qs_env, do_kpoints=do_kpoints)

         IF (do_kpoints) THEN

            CALL get_qs_env(qs_env, &
                            matrix_w_kp=matrix_w_kp, &
                            matrix_s_kp=matrix_s_kp, &
                            sab_orb=sab_nl, &
                            mos=mos, &
                            kpoints=kpoints)

            CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
            CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, nrow_global=nao, ncol_global=nao, &
                                     template_fmstruct=mo_coeff%matrix_struct)

            ALLOCATE (fmwork(2))
            DO is = 1, SIZE(fmwork)
               NULLIFY (fmwork(is)%matrix)
               CALL cp_fm_create(fmwork(is)%matrix, matrix_struct=ao_ao_fmstruct)
            END DO
            CALL cp_fm_struct_release(ao_ao_fmstruct)

            ! energy weighted density matrices in k-space
            CALL kpoint_density_matrices(kpoints, energy_weighted=.TRUE.)
            ! energy weighted density matrices in real space
            CALL kpoint_density_transform(kpoints, matrix_w_kp, .TRUE., &
                                          matrix_s_kp(1, 1)%matrix, sab_nl, fmwork)

            DO is = 1, SIZE(fmwork)
               CALL cp_fm_release(fmwork(is)%matrix)
            END DO
            DEALLOCATE (fmwork)

         ELSE

            NULLIFY (dft_control, rho_ao)
            CALL get_qs_env(qs_env, &
                            matrix_w=matrix_w, &
                            matrix_ks=matrix_ks, &
                            matrix_s=matrix_s, &
                            matrix_w_mp2=matrix_w_mp2, &
                            mo_derivs=mo_derivs, &
                            scf_control=scf_control, &
                            mos=mos, &
                            rho=rho, &
                            dft_control=dft_control)

            CALL qs_rho_get(rho, rho_ao=rho_ao)

            nspin = SIZE(mos)
            DO ispin = 1, nspin
               mo_set => mos(ispin)%mo_set
               IF (dft_control%roks) THEN
                  IF (scf_control%use_ot) THEN
                     IF (ispin > 1) THEN
                        ! not very elegant, indeed ...
                        CALL cp_dbcsr_set(matrix_w(ispin)%matrix, 0.0_dp)
                     ELSE
                        CALL calculate_w_matrix_ot(mo_set, mo_derivs(ispin)%matrix, &
                                                   matrix_w(ispin)%matrix, matrix_s(1)%matrix)
                     END IF
                  ELSE
                     CALL calculate_w_matrix(mo_set=mo_set, &
                                             matrix_ks=matrix_ks(ispin)%matrix, &
                                             matrix_p=rho_ao(ispin)%matrix, &
                                             matrix_w=matrix_w(ispin)%matrix)
                  END IF
               ELSE
                  IF (scf_control%use_ot) THEN
                     CALL calculate_w_matrix_ot(mo_set, mo_derivs(ispin)%matrix, &
                                                matrix_w(ispin)%matrix, matrix_s(1)%matrix)
                  ELSE
                     CALL calculate_w_matrix(mo_set, matrix_w(ispin)%matrix)
                  END IF
               END IF
            END DO

         END IF

      END IF
      CALL timestop(handle)

   END SUBROUTINE qs_energies_compute_w

END MODULE qs_energy_matrix_w
