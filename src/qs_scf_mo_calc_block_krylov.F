!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_mo_calc_block_krylov
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_create
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_scf_abstract_mo_calc,         ONLY: AbstractMOCalc
   USE qs_scf_diagonalization,          ONLY: do_block_krylov_diag,&
                                              do_general_diag,&
                                              do_scf_diag_subspace
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! CASE (block_krylov_diag_method_nr, general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr)

   TYPE, PUBLIC, EXTENDS(AbstractMOCalc) :: MOCalcBlockKrylov
   PRIVATE

   TYPE(qs_environment_type), POINTER         :: qs_env => null()
   TYPE(qs_scf_env_type), POINTER             :: scf_env => null()
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_ks => null(), matrix_s => null()
   TYPE(section_vals_type), POINTER           :: scf_section => null()
   TYPE(scf_control_type), POINTER            :: scf_control => null()
   TYPE(qs_ks_env_type), POINTER              :: ks_env => null()
   TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos => null()
   TYPE(qs_rho_type), POINTER                 :: rho => null()

CONTAINS
   PROCEDURE, PUBLIC, PASS(self) :: init => block_krylov_init
   PROCEDURE, PUBLIC, PASS(self) :: pre => block_krylov_pre
   PROCEDURE, PUBLIC, PASS(self) :: run => block_krylov_run
END TYPE

CONTAINS
SUBROUTINE block_krylov_init(self, qs_env, scf_env, scf_section)
   CLASS(MOCalcBlockKrylov), INTENT(INOUT)             :: self
   TYPE(qs_environment_type), POINTER, INTENT(IN) :: qs_env
   TYPE(qs_scf_env_type), POINTER, INTENT(IN)     :: scf_env
   TYPE(section_vals_type), POINTER, INTENT(IN)   :: scf_section

   self%qs_env => qs_env
   self%scf_env => scf_env
   self%scf_section => scf_section

   CALL get_qs_env(qs_env, &
                   scf_control=self%scf_control, &
                   matrix_ks=self%matrix_ks, &
                   matrix_s=self%matrix_s, &
                   ks_env=self%ks_env, &
                   rho=self%rho, &
                   mos=self%mos)

   IF (.NOT. scf_env%skip_diis) THEN
      IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
         CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=self%scf_control%max_diis)
      END IF
   END IF
END SUBROUTINE block_krylov_init

SUBROUTINE block_krylov_pre(self)
   CLASS(MOCalcBlockKrylov), INTENT(INOUT) :: self

   CALL qs_diis_b_clear(self%scf_env%scf_diis_buffer)
END SUBROUTINE block_krylov_pre

SUBROUTINE block_krylov_run(self, diis_step, energy_only)
   !! the inner loop of scf, specific to iterative diagonalization using OT
   !! with S matrix; basically, in goes the ks matrix out goes a new p matrix
   CLASS(MOCalcBlockKrylov), INTENT(INOUT)   :: self
   LOGICAL, INTENT(INOUT)               :: diis_step
   LOGICAL, INTENT(INOUT)               :: energy_only
   LOGICAL :: skip_diag_sub

   MARK_USED(energy_only)

   ASSOCIATE (scf_env=>self%scf_env)
      IF ((scf_env%krylov_space%eps_std_diag > 0.0_dp) .AND. &
          (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%krylov_space%eps_std_diag)) THEN
         IF (scf_env%krylov_space%always_check_conv) THEN
            ! something is fishy with do_block_krylov_diag (missing dummy argument specs?!)
            CALL do_block_krylov_diag(self%scf_env, self%mos, self%matrix_ks, &
                                      self%scf_control, self%scf_section, check_moconv_only=.TRUE.)
         END IF
         CALL do_general_diag(self%scf_env, self%mos, self%matrix_ks, &
                              self%matrix_s, self%scf_control, self%scf_section, diis_step)
      ELSE
         CALL do_block_krylov_diag(self%scf_env, self%mos, self%matrix_ks, &
                                   self%scf_control, self%scf_section)
      END IF
      IF (self%scf_control%do_diag_sub) THEN
         skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                         (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%subspace_env%eps_diag_sub)
         IF (.NOT. skip_diag_sub) THEN
            CALL do_scf_diag_subspace(self%qs_env, self%scf_env, scf_env%subspace_env, self%mos, self%rho, &
                                      self%ks_env, self%scf_section, self%scf_control)
         END IF
      END IF
   END ASSOCIATE
END SUBROUTINE block_krylov_run

END MODULE qs_scf_mo_calc_block_krylov
