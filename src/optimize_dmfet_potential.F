!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE optimize_dmfet_potential
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_trace,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE embed_types,                     ONLY: opt_dmfet_pot_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_dmfet_potential'

   PUBLIC :: build_full_dm, prepare_dmfet_opt, release_dmfet_opt, subsys_spin, check_dmfet

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
!> \param opt_dmfet_section ...
! **************************************************************************************************
   SUBROUTINE read_opt_dmfet_section(opt_dmfet, opt_dmfet_section)
      TYPE(opt_dmfet_pot_type)                           :: opt_dmfet
      TYPE(section_vals_type), POINTER                   :: opt_dmfet_section

      ! Read keywords

      CALL section_vals_val_get(opt_dmfet_section, "N_ITER", &
                                i_val=opt_dmfet%n_iter)

      CALL section_vals_val_get(opt_dmfet_section, "TRUST_RAD", &
                                r_val=opt_dmfet%trust_rad)

      CALL section_vals_val_get(opt_dmfet_section, "DM_CONV_MAX", &
                                r_val=opt_dmfet%conv_max)

      CALL section_vals_val_get(opt_dmfet_section, "DM_CONV_INT", &
                                r_val=opt_dmfet%conv_int)

      CALL section_vals_val_get(opt_dmfet_section, "BETA_DM_CONV_MAX", &
                                r_val=opt_dmfet%conv_max_beta)

      CALL section_vals_val_get(opt_dmfet_section, "BETA_DM_CONV_INT", &
                                r_val=opt_dmfet%conv_int_beta)

      CALL section_vals_val_get(opt_dmfet_section, "READ_DMFET_POT", &
                                l_val=opt_dmfet%read_dmfet_pot)

   END SUBROUTINE read_opt_dmfet_section

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION subsys_spin(qs_env) RESULT(subsys_open_shell)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: subsys_open_shell

      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      subsys_open_shell = .FALSE.
      IF (dft_control%nspins == 2) subsys_open_shell = .TRUE.

   END FUNCTION subsys_spin

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param opt_dmfet ...
!> \param opt_dmfet_section ...
! **************************************************************************************************
   SUBROUTINE prepare_dmfet_opt(qs_env, opt_dmfet, opt_dmfet_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(opt_dmfet_pot_type)                           :: opt_dmfet
      TYPE(section_vals_type), POINTER                   :: opt_dmfet_section

      INTEGER                                            :: diff_size, nao
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      ! Read the input
      CALL read_opt_dmfet_section(opt_dmfet, opt_dmfet_section)

      ! Get the orbital coefficients
      CALL get_qs_env(qs_env=qs_env, mos=mos, para_env=para_env)
      CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)

      ! Make cp_fm matrices
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      NULLIFY (opt_dmfet%dmfet_pot, opt_dmfet%dm_1, opt_dmfet%dm_2, opt_dmfet%dm_total, opt_dmfet%dm_diff)
      NULLIFY (fm_struct)

      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=nao, ncol_global=nao)
      CALL cp_fm_create(opt_dmfet%dmfet_pot, fm_struct, name="dmfet_pot")
      CALL cp_fm_create(opt_dmfet%dm_subsys, fm_struct, name="dm_subsys")
      CALL cp_fm_create(opt_dmfet%dm_total, fm_struct, name="dm_total")
      CALL cp_fm_create(opt_dmfet%dm_diff, fm_struct, name="dm_diff")

      CALL cp_fm_set_all(opt_dmfet%dmfet_pot, 0.0_dp)
      CALL cp_fm_set_all(opt_dmfet%dm_subsys, 0.0_dp)
      CALL cp_fm_set_all(opt_dmfet%dm_total, 0.0_dp)
      CALL cp_fm_set_all(opt_dmfet%dm_diff, 0.0_dp)

      ! Beta spin counterparts
      IF (opt_dmfet%open_shell_embed) THEN
         CALL cp_fm_create(opt_dmfet%dmfet_pot_beta, fm_struct, name="dmfet_pot_beta")
         CALL cp_fm_create(opt_dmfet%dm_subsys_beta, fm_struct, name="dm_subsys_beta")
         CALL cp_fm_create(opt_dmfet%dm_total_beta, fm_struct, name="dm_total_beta")
         CALL cp_fm_create(opt_dmfet%dm_diff_beta, fm_struct, name="dm_diff_beta")

         CALL cp_fm_set_all(opt_dmfet%dmfet_pot_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_dmfet%dm_subsys_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_dmfet%dm_total_beta, 0.0_dp)
         CALL cp_fm_set_all(opt_dmfet%dm_diff_beta, 0.0_dp)
      END IF

      ! Release structure and context
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_blacs_env_release(blacs_env)

      ! Array to store functional values
      ALLOCATE (opt_dmfet%w_func(opt_dmfet%n_iter))
      opt_dmfet%w_func = 0.0_dp

      ! Allocate max_diff and int_diff
      diff_size = 1
      IF (opt_dmfet%open_shell_embed) diff_size = 2
      ALLOCATE (opt_dmfet%max_diff(diff_size))
      ALLOCATE (opt_dmfet%int_diff(diff_size))

   END SUBROUTINE prepare_dmfet_opt

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
! **************************************************************************************************
   SUBROUTINE release_dmfet_opt(opt_dmfet)
      TYPE(opt_dmfet_pot_type)                           :: opt_dmfet

      CALL cp_fm_release(opt_dmfet%dmfet_pot)
      CALL cp_fm_release(opt_dmfet%dm_subsys)
      CALL cp_fm_release(opt_dmfet%dm_total)
      CALL cp_fm_release(opt_dmfet%dm_diff)

      IF (opt_dmfet%open_shell_embed) THEN
         CALL cp_fm_release(opt_dmfet%dmfet_pot_beta)
         CALL cp_fm_release(opt_dmfet%dm_subsys_beta)
         CALL cp_fm_release(opt_dmfet%dm_total_beta)
         CALL cp_fm_release(opt_dmfet%dm_diff_beta)
      END IF

      DEALLOCATE (opt_dmfet%w_func)
      DEALLOCATE (opt_dmfet%max_diff)
      DEALLOCATE (opt_dmfet%int_diff)
      DEALLOCATE (opt_dmfet%all_nspins)

   END SUBROUTINE release_dmfet_opt

! **************************************************************************************************
!> \brief Builds density matrices from MO coefficients in full matrix format
!> \param qs_env ...
!> \param dm ...
!> \param open_shell  Subsystem is open shell
!> \param open_shell_embed  Embedding is open shell
!> \param dm_beta ...
! **************************************************************************************************
   SUBROUTINE build_full_dm(qs_env, dm, open_shell, open_shell_embed, dm_beta)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: dm
      LOGICAL                                            :: open_shell, open_shell_embed
      TYPE(cp_fm_type), POINTER                          :: dm_beta

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_full_dm'

      INTEGER                                            :: handle, homo, nao
      REAL(KIND=dp)                                      :: coeff
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL timeset(routineN, handle)

      coeff = 2.0_dp
      IF (open_shell_embed) coeff = 1.0_dp

      ! Get the orbital coefficients
      CALL get_qs_env(qs_env=qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, homo=homo)

      ! Build the density matrix
      CALL cp_gemm(transa="N", transb="T", m=nao, n=nao, &
                   k=homo, alpha=coeff, &
                   matrix_a=mo_coeff, matrix_b=mo_coeff, &
                   beta=0.0_dp, matrix_c=dm)

      ! Open shell
      IF (open_shell) THEN
         CALL get_mo_set(mo_set=mos(2)%mo_set, mo_coeff=mo_coeff, nao=nao, homo=homo)

         ! Build the density matrix
         CALL cp_gemm(transa="N", transb="T", m=nao, n=nao, &
                      k=homo, alpha=coeff, &
                      matrix_a=mo_coeff, matrix_b=mo_coeff, &
                      beta=0.0_dp, matrix_c=dm_beta)
      END IF

      ! If embedding is open shell, but subsystem is not, copy alpha-spin DM intou beta-spin DM
      IF (open_shell_embed .AND. (.NOT. open_shell)) THEN
         CALL get_qs_env(qs_env=qs_env, para_env=para_env)
         CALL cp_fm_copy_general(dm, dm_beta, para_env)

      END IF

      CALL timestop(handle)

   END SUBROUTINE build_full_dm

! **************************************************************************************************
!> \brief ...
!> \param opt_dmfet ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE check_dmfet(opt_dmfet, qs_env)
      TYPE(opt_dmfet_pot_type)                           :: opt_dmfet
      TYPE(qs_environment_type), POINTER                 :: qs_env

      REAL(KIND=dp)                                      :: max_diff, max_diff_beta, trace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: diff_dbcsr, dm_s

      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      NULLIFY (diff_dbcsr)
      CALL dbcsr_init_p(diff_dbcsr)
      CALL dbcsr_create(matrix=diff_dbcsr, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (dm_s)
      CALL dbcsr_init_p(dm_s)
      CALL dbcsr_create(matrix=dm_s, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(opt_dmfet%dm_diff, diff_dbcsr, keep_sparsity=.FALSE.)

      CALL dbcsr_multiply("N", "N", 1.0_dp, diff_dbcsr, matrix_s(1)%matrix, &
                          0.0_dp, dm_s)

      CALL dbcsr_trace(dm_s, trace)

      IF (opt_dmfet%open_shell_embed) THEN
         CALL copy_fm_to_dbcsr(opt_dmfet%dm_diff_beta, diff_dbcsr, keep_sparsity=.FALSE.)

         CALL dbcsr_multiply("N", "N", 1.0_dp, diff_dbcsr, matrix_s(1)%matrix, &
                             0.0_dp, dm_s)

         CALL dbcsr_trace(dm_s, trace)

      END IF

      ! Release dbcsr objects
      CALL dbcsr_release(diff_dbcsr)
      CALL dbcsr_release(dm_s)

      ! Find the absolute maximum value of the DM difference
      CALL cp_fm_maxabsval(opt_dmfet%dm_diff, max_diff)
      IF (opt_dmfet%open_shell_embed) THEN
         CALL cp_fm_maxabsval(opt_dmfet%dm_diff_beta, max_diff_beta)
      END IF

   END SUBROUTINE check_dmfet

END MODULE optimize_dmfet_potential
