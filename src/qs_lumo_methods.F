!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities to construct additional unoccupied GPW/GAPW orbitals.
! **************************************************************************************************
MODULE qs_lumo_methods
   USE admm_types,                      ONLY: admm_type
   USE admm_utils,                      ONLY: admm_correct_for_eigenvalues,&
                                              admm_uncorrect_for_eigenvalues
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_init_random,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_constants,                 ONLY: ot_precond_full_all
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE preconditioner_types,            ONLY: preconditioner_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_lumo_methods'

   PUBLIC :: make_lumo_gpw, make_lumo_gpw_spin

CONTAINS

! **************************************************************************************************
!> \brief Gets the LUMOs and their eigenvalues for all spin channels.
!> \param qs_env ...
!> \param scf_env ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param nlumo ...
!> \param nlumos ...
! **************************************************************************************************
   SUBROUTINE make_lumo_gpw(qs_env, scf_env, unoccupied_orbs, unoccupied_evals, nlumo, nlumos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: unoccupied_orbs
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: unoccupied_evals
      INTEGER, INTENT(IN)                                :: nlumo
      INTEGER, INTENT(OUT)                               :: nlumos

      CHARACTER(len=*), PARAMETER                        :: routineN = 'make_lumo_gpw'

      INTEGER                                            :: handle, ispin
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      NULLIFY (mos, dft_control)
      CALL get_qs_env(qs_env, mos=mos, dft_control=dft_control)

      DO ispin = 1, dft_control%nspins
         CALL make_lumo_gpw_spin(qs_env, scf_env, mos(ispin), unoccupied_orbs(ispin), &
                                 unoccupied_evals(ispin)%array, nlumo, nlumos, ispin)
      END DO

      CALL timestop(handle)

   END SUBROUTINE make_lumo_gpw

! **************************************************************************************************
!> \brief Gets the LUMOs and their eigenvalues for one spin channel.
!> \param qs_env ...
!> \param scf_env ...
!> \param mo_set ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param nlumo ...
!> \param nlumos ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE make_lumo_gpw_spin(qs_env, scf_env, mo_set, unoccupied_orbs, unoccupied_evals, &
                                 nlumo, nlumos, ispin)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(mo_set_type), INTENT(IN)                      :: mo_set
      TYPE(cp_fm_type), INTENT(OUT)                      :: unoccupied_orbs
      REAL(KIND=dp), DIMENSION(:), POINTER               :: unoccupied_evals
      INTEGER, INTENT(IN)                                :: nlumo
      INTEGER, INTENT(OUT)                               :: nlumos
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'make_lumo_gpw_spin'

      INTEGER                                            :: handle, homo, n, nao, nmo, output_unit
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(preconditioner_type), POINTER                 :: local_preconditioner
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (ks_rmpv, matrix_s, scf_control, dft_control, admm_env, para_env, blacs_env)
      CALL get_qs_env(qs_env, &
                      matrix_ks=ks_rmpv, &
                      matrix_s=matrix_s, &
                      scf_control=scf_control, &
                      dft_control=dft_control, &
                      admm_env=admm_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      NULLIFY (unoccupied_evals)
      IF (output_unit > 0) WRITE (output_unit, *) " "
      IF (output_unit > 0) WRITE (output_unit, *) " Lowest Eigenvalues of the unoccupied subspace spin ", ispin
      IF (output_unit > 0) WRITE (output_unit, FMT='(1X,A)') "-----------------------------------------------------"
      CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, homo=homo, nao=nao, nmo=nmo)
      CALL cp_fm_get_info(mo_coeff, nrow_global=n)
      nlumos = MAX(1, MIN(nlumo, nao - nmo))
      IF (nlumo == -1) nlumos = nao - nmo
      ALLOCATE (unoccupied_evals(nlumos))
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=n, ncol_global=nlumos)
      CALL cp_fm_create(unoccupied_orbs, fm_struct_tmp, name="lumos")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_init_random(unoccupied_orbs, nlumos)

      ! the full_all preconditioner makes not much sense for lumos search
      NULLIFY (local_preconditioner)
      IF (ASSOCIATED(scf_env)) THEN
         IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
            local_preconditioner => scf_env%ot_preconditioner(1)%preconditioner
            ! this one can for sure not be right (as it has to match a given C0)
            IF (local_preconditioner%in_use == ot_precond_full_all) THEN
               NULLIFY (local_preconditioner)
            END IF
         END IF
      END IF

      ! If we do ADMM, we add have to modify the Kohn-Sham matrix
      IF (dft_control%do_admm) THEN
         CALL admm_correct_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix)
      END IF

      CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix, matrix_s=matrix_s(1)%matrix, &
                          matrix_c_fm=unoccupied_orbs, &
                          matrix_orthogonal_space_fm=mo_coeff, &
                          eps_gradient=scf_control%eps_lumos, &
                          preconditioner=local_preconditioner, &
                          iter_max=scf_control%max_iter_lumos, &
                          size_ortho_space=nmo)

      CALL calculate_subspace_eigenvalues(unoccupied_orbs, ks_rmpv(ispin)%matrix, &
                                          unoccupied_evals, scr=output_unit, &
                                          ionode=output_unit > 0)

      ! If we do ADMM, we restore the original Kohn-Sham matrix
      IF (dft_control%do_admm) THEN
         CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE make_lumo_gpw_spin

END MODULE qs_lumo_methods
