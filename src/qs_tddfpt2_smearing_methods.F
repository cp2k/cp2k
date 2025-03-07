!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_smearing_methods
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              smeared_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE fermi_utils,                     ONLY: Fermi
   USE input_constants,                 ONLY: smear_fermi_dirac
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE scf_control_types,               ONLY: scf_control_type,&
                                              smear_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_smearing_methods'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   PUBLIC :: tddfpt_smeared_occupation, &
             add_smearing_aterm, compute_fermib, &
             orthogonalize_smeared_occupation, deallocate_fermi_params

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gs_mos ...
!> \param log_unit ...
! **************************************************************************************************
   SUBROUTINE tddfpt_smeared_occupation(qs_env, gs_mos, log_unit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(IN), POINTER                             :: gs_mos
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_smeared_occupation'

      INTEGER                                            :: handle, iocc, ispin, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nocc, nvirt
      REAL(kind=dp), DIMENSION(:), POINTER               :: mo_evals, occup
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(smear_type), POINTER                          :: smear

      CALL timeset(routineN, handle)

      nspins = SIZE(gs_mos)

      NULLIFY (mos, scf_control)
      CALL get_qs_env(qs_env, mos=mos, scf_control=scf_control)
      NULLIFY (smear)
      IF (ASSOCIATED(qs_env%scf_control%smear)) THEN
         smear => qs_env%scf_control%smear
      ELSE
         CPABORT("Smeared input section no longer associated.")
      END IF

      IF (debug_this_module) THEN
         NULLIFY (mo_evals, occup)
         ALLOCATE (nocc(nspins), nvirt(nspins))
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), eigenvalues=mo_evals, occupation_numbers=occup)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nocc(ispin))
            CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, ncol_global=nvirt(ispin))
            IF (log_unit > 0) THEN
               DO iocc = 1, nocc(ispin)
                  WRITE (log_unit, '(A,F14.5)') "Occupation numbers", occup(iocc)
               END DO
            END IF
         END DO
      END IF

      CALL allocate_fermi_params(qs_env, gs_mos)
      CALL compute_fermia(qs_env, gs_mos, log_unit)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_smeared_occupation
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gs_mos ...
!> \param log_unit ...
! **************************************************************************************************
   SUBROUTINE compute_fermia(qs_env, gs_mos, log_unit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(IN), POINTER                             :: gs_mos
      INTEGER, INTENT(IN)                                :: log_unit

      CHARACTER(len=*), PARAMETER                        :: routineN = 'compute_fermia'

      INTEGER                                            :: handle, iocc, ispin, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nocc
      REAL(kind=dp)                                      :: maxvalue, mu
      REAL(kind=dp), DIMENSION(:), POINTER               :: fermia, mo_evals, occup
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(smear_type), POINTER                          :: smear
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (mos, dft_control, tddfpt_control, scf_control)
      CALL get_qs_env(qs_env, dft_control=dft_control, scf_control=scf_control)
      tddfpt_control => dft_control%tddfpt2_control

      CALL get_qs_env(qs_env, mos=mos)
      nspins = SIZE(mos)

      NULLIFY (smear)
      IF (ASSOCIATED(qs_env%scf_control%smear)) THEN
         smear => qs_env%scf_control%smear
      ELSE
         CPABORT("Smeared input section no longer associated.")
      END IF

      IF (debug_this_module .AND. (log_unit > 0)) THEN
         WRITE (log_unit, '(A,F14.5)') "Smearing temperature", smear%electronic_temperature
      END IF

      NULLIFY (mo_evals, occup)
      ALLOCATE (nocc(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin), eigenvalues=mo_evals, occupation_numbers=occup, mu=mu)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nocc(ispin))
      END DO

      DO ispin = 1, nspins
         fermia => tddfpt_control%smeared_occup(ispin)%fermia
         DO iocc = 1, nocc(ispin)
            maxvalue = mu + 3.0_dp*smear%electronic_temperature - mo_evals(iocc)

            fermia(iocc) = MAX(0.0_dp, maxvalue)
            IF (debug_this_module) THEN
               IF (log_unit > 0) WRITE (log_unit, '(A,F14.5)') "Fermi smearing parameter alpha", fermia(iocc)
            END IF
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE compute_fermia
! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param Aop_evects ...
!> \param evects ...
!> \param S_evects ...
!> \param mos_occ ...
!> \param fermia ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE add_smearing_aterm(qs_env, Aop_evects, evects, S_evects, mos_occ, fermia, matrix_s)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(in)                       :: Aop_evects, evects, S_evects, mos_occ
      REAL(kind=dp), DIMENSION(:), POINTER               :: fermia
      TYPE(dbcsr_type), POINTER                          :: matrix_s

      CHARACTER(len=*), PARAMETER :: routineN = 'add_smearing_aterm'

      INTEGER                                            :: handle, iounit, nactive, nao
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct, matrix_struct_tmp
      TYPE(cp_fm_type)                                   :: CCSXvec, CSXvec, Cvec, SCCSXvec
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      NULLIFY (dft_control, tddfpt_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt2_control

      CALL cp_fm_get_info(matrix=evects, matrix_struct=matrix_struct, &
                          nrow_global=nao, ncol_global=nactive)
      CALL cp_fm_create(CCSXvec, matrix_struct)
      CALL cp_fm_create(SCCSXvec, matrix_struct)
      CALL cp_fm_create(Cvec, matrix_struct)

      NULLIFY (matrix_struct_tmp)
      CALL cp_fm_struct_create(matrix_struct_tmp, nrow_global=nactive, &
                               ncol_global=nactive, template_fmstruct=matrix_struct)
      CALL cp_fm_create(CSXvec, matrix_struct_tmp)
      CALL cp_fm_struct_release(fmstruct=matrix_struct_tmp)

      CALL parallel_gemm('T', 'N', nactive, nactive, nao, 1.0_dp, &
                         mos_occ, S_evects, 0.0_dp, CSXvec)
      CALL cp_fm_to_fm(mos_occ, Cvec)

      CALL cp_fm_column_scale(Cvec, fermia)
      CALL parallel_gemm('N', 'N', nao, nactive, nactive, 1.0_dp, &
                         Cvec, CSXvec, 0.0_dp, CCSXvec)

      ! alpha S C C^T S X
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, CCSXvec, &
                                   SCCSXvec, ncol=nactive, alpha=1.0_dp, beta=0.0_dp)
      CALL cp_fm_scale_and_add(1.0_dp, Aop_evects, 1.0_dp, SCCSXvec)

      CALL cp_fm_release(SCCSXvec)
      CALL cp_fm_release(CCSXvec)
      CALL cp_fm_release(CSXvec)
      CALL cp_fm_release(Cvec)

      CALL timestop(handle)
   END SUBROUTINE add_smearing_aterm
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gs_mos ...
!> \param evals ...
! **************************************************************************************************
   SUBROUTINE compute_fermib(qs_env, gs_mos, evals)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(IN)                                      :: gs_mos
      REAL(kind=dp), INTENT(in)                          :: evals

      CHARACTER(len=*), PARAMETER                        :: routineN = 'compute_fermib'

      INTEGER                                            :: handle, iocc, iounit, ispin, jocc, &
                                                            nactive, nspins
      REAL(KIND=dp)                                      :: dummykTS, nelec
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: interocc, occup_im
      REAL(kind=dp), DIMENSION(:), POINTER               :: fermia, mo_evals, occup, occuptmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fermib
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(smear_type), POINTER                          :: smear
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      NULLIFY (mos, scf_control)
      CALL get_qs_env(qs_env, mos=mos, scf_control=scf_control)
      NULLIFY (smear)
      IF (ASSOCIATED(qs_env%scf_control%smear)) THEN
         smear => qs_env%scf_control%smear
      ELSE
         CPABORT("Smeared input section no longer associated.")
      END IF

      NULLIFY (dft_control, tddfpt_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt2_control

      NULLIFY (fermib, fermia)
      nspins = SIZE(gs_mos)

      DO ispin = 1, nspins

         fermib => dft_control%tddfpt2_control%smeared_occup(ispin)%fermib
         fermia => dft_control%tddfpt2_control%smeared_occup(ispin)%fermia
         fermib = 0.0_dp

         ! get theta_Fi
         NULLIFY (mo_evals, occup)
         CALL get_mo_set(mos(ispin), eigenvalues=mo_evals, occupation_numbers=occup)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nactive)

         IF (smear%fixed_mag_mom == -1.0_dp) THEN
            nelec = REAL(mos(ispin)%nelectron, dp)
         ELSE
            nelec = mos(ispin)%n_el_f
         END IF

         ! compute theta_im
         NULLIFY (occuptmp)
         CALL get_mo_set(mos(ispin), occupation_numbers=occuptmp)
         ALLOCATE (occup_im(nactive, nactive), interocc(nactive, nactive))

         DO iocc = 1, nactive
            IF (smear%method .EQ. smear_fermi_dirac) THEN
               ! Different prefactor in comparison to Octopus !
               CALL Fermi(occuptmp, nelec, dummykTS, mos(ispin)%eigenvalues, mos(ispin)%eigenvalues(iocc), &
                          smear%electronic_temperature, mos(ispin)%maxocc)
            ELSE
               CPABORT("TDDFT with smearing only works with Fermi-Dirac distribution.")
            END IF
            DO jocc = 1, nactive
               occup_im(iocc, jocc) = occuptmp(jocc)
            END DO
         END DO

         ! compute fermib
         DO iocc = 1, nactive
            DO jocc = 1, nactive
               interocc(iocc, jocc) = (occup(iocc) - occup(jocc))/(mo_evals(iocc) - mo_evals(jocc) - evals)
               fermib(iocc, jocc) = fermib(iocc, jocc) + occup(iocc)*occup_im(iocc, jocc) &
                                    + occup(jocc)*occup_im(jocc, iocc) &
                                    + fermia(jocc)*interocc(iocc, jocc)*occup_im(jocc, iocc)
            END DO
         END DO

         IF (debug_this_module .AND. (iounit > 0)) THEN
            DO iocc = 1, nactive
               DO jocc = 1, nactive
                  WRITE (iounit, '(A,F14.5)') "Fermi smearing parameter beta,", fermib(iocc, jocc)
               END DO
            END DO
         END IF

         DEALLOCATE (occup_im)
         DEALLOCATE (interocc)

      END DO  ! ispin

      CALL timestop(handle)
   END SUBROUTINE compute_fermib
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gs_mos ...
! **************************************************************************************************
   SUBROUTINE allocate_fermi_params(qs_env, gs_mos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(IN), POINTER                             :: gs_mos

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_fermi_params'

      INTEGER                                            :: handle, ispin, nocc, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(smeared_type), DIMENSION(:), POINTER          :: smeared_occup
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, tddfpt_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt2_control

      NULLIFY (smeared_occup)
      smeared_occup => dft_control%tddfpt2_control%smeared_occup
      nspins = SIZE(gs_mos)
      ALLOCATE (smeared_occup(nspins))

      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nocc)
         ALLOCATE (smeared_occup(ispin)%fermia(nocc))
         ALLOCATE (smeared_occup(ispin)%fermib(nocc, nocc))
      END DO
      dft_control%tddfpt2_control%smeared_occup => smeared_occup

      CALL timestop(handle)
   END SUBROUTINE allocate_fermi_params
! **************************************************************************************************
!> \brief ...
!> \param smeared_occup ...
! **************************************************************************************************
   SUBROUTINE deallocate_fermi_params(smeared_occup)

      TYPE(smeared_type), DIMENSION(:), POINTER          :: smeared_occup

      INTEGER                                            :: ispin

      IF (ASSOCIATED(smeared_occup)) THEN
         DO ispin = 1, SIZE(smeared_occup)
            IF (ASSOCIATED(smeared_occup(ispin)%fermia)) THEN
               DEALLOCATE (smeared_occup(ispin)%fermia)
               DEALLOCATE (smeared_occup(ispin)%fermib)
               NULLIFY (smeared_occup(ispin)%fermia, smeared_occup(ispin)%fermib)
            END IF
         END DO
         DEALLOCATE (smeared_occup)
         NULLIFY (smeared_occup)
      END IF

   END SUBROUTINE deallocate_fermi_params
! **************************************************************************************************
!> \brief ...
!> \param evects ...
!> \param qs_env ...
!> \param mos_occ ...
!> \param S_C0 ...
! **************************************************************************************************
   SUBROUTINE orthogonalize_smeared_occupation(evects, qs_env, mos_occ, S_C0)

      TYPE(cp_fm_type), DIMENSION(:), INTENT(in)         :: evects
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(in)         :: mos_occ, S_C0

      CHARACTER(LEN=*), PARAMETER :: routineN = 'orthogonalize_smeared_occupation'

      INTEGER                                            :: handle, iocc, iounit, ispin, nactive, &
                                                            nao, nspins
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: bscale
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: evortholocal, subevects, subevectsresult
      REAL(kind=dp), DIMENSION(:), POINTER               :: occup
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fermib
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type)                                   :: betaSCC, Cvec, evortho, subevectsfm, &
                                                            subevectsresultfm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      NULLIFY (mos)
      CALL get_qs_env(qs_env, mos=mos)
      NULLIFY (dft_control, tddfpt_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt2_control

      CALL cp_fm_get_info(matrix=evects(1), matrix_struct=matrix_struct, &
                          nrow_global=nao, ncol_global=nactive)
      CALL cp_fm_create(evortho, matrix_struct)

      nspins = SIZE(evects)
      NULLIFY (para_env)
      CALL cp_fm_get_info(evects(1), para_env=para_env, context=blacs_env_global)

      NULLIFY (matrix_struct)
      CALL cp_fm_struct_create(matrix_struct, nrow_global=nao, ncol_global=nao, context=blacs_env_global)
      CALL cp_fm_create(betaSCC, matrix_struct)
      CALL cp_fm_struct_release(fmstruct=matrix_struct)

      ALLOCATE (evortholocal(nao, nactive))
      ALLOCATE (bscale(nactive))

      DO ispin = 1, nspins
         NULLIFY (matrix_struct)
         CALL cp_fm_get_info(matrix=mos_occ(ispin), matrix_struct=matrix_struct)
         CALL cp_fm_create(Cvec, matrix_struct)

         NULLIFY (occup)
         CALL get_mo_set(mos(ispin), occupation_numbers=occup, mo_coeff=mo_coeff)

         NULLIFY (fermib)
         IF (.NOT. ASSOCIATED(dft_control%tddfpt2_control%smeared_occup)) THEN
            CPABORT("Smeared occupation intermediates not associated.")
         END IF
         fermib => dft_control%tddfpt2_control%smeared_occup(ispin)%fermib

         DO iocc = 1, nactive
            CALL cp_fm_copy_general(mos_occ(ispin), Cvec, para_env)
            bscale(:) = fermib(iocc, :)
            CALL cp_fm_column_scale(Cvec, bscale)

            CALL parallel_gemm('N', 'T', nao, nao, nactive, 1.0_dp, Cvec, S_C0(ispin), 0.0_dp, betaSCC)

            ! get ith column of X
            NULLIFY (matrix_struct)
            CALL cp_fm_struct_create(matrix_struct, nrow_global=nao, ncol_global=1, context=blacs_env_global)
            CALL cp_fm_create(subevectsfm, matrix_struct)
            CALL cp_fm_create(subevectsresultfm, matrix_struct)
            CALL cp_fm_struct_release(fmstruct=matrix_struct)

            ALLOCATE (subevects(nao, 1))
            ALLOCATE (subevectsresult(nao, 1))
            CALL cp_fm_get_submatrix(fm=evects(1), target_m=subevects, &
                                     start_row=1, start_col=iocc, n_rows=nao, n_cols=1)
            CALL cp_fm_set_submatrix(subevectsfm, subevects, &
                                     start_row=1, start_col=1, n_rows=nao, n_cols=1)

            CALL parallel_gemm('N', 'N', nao, 1, nao, 1.0_dp, betaSCC, &
                               subevectsfm, 0.0_dp, subevectsresultfm)

            CALL cp_fm_get_submatrix(fm=subevectsresultfm, target_m=subevectsresult, &
                                     start_row=1, start_col=1, n_rows=nao, n_cols=1)
            CALL cp_fm_set_submatrix(evortho, subevectsresult, &
                                     start_row=1, start_col=iocc, n_rows=nao, n_cols=1)
            DEALLOCATE (subevects, subevectsresult)
            CALL cp_fm_release(subevectsfm)
            CALL cp_fm_release(subevectsresultfm)
         END DO ! iocc
      END DO ! nspins

      CALL cp_fm_column_scale(evects(1), occup)
      CALL cp_fm_scale_and_add(1.0_dp, evects(1), -1.0_dp, evortho)

      DEALLOCATE (bscale)
      DEALLOCATE (evortholocal)

      CALL cp_fm_release(betaSCC)
      CALL cp_fm_release(Cvec)

      CALL cp_fm_release(evortho)

      CALL timestop(handle)
   END SUBROUTINE orthogonalize_smeared_occupation

END MODULE qs_tddfpt2_smearing_methods
