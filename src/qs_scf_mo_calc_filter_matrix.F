!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_mo_calc_filter_matrix
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_create
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fb_env_methods,               ONLY: fb_env_do_diag
   USE qs_scf_abstract_mo_calc,         ONLY: AbstractMOCalc
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! CASE (filter_matrix_diag_method_nr, general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr)

   TYPE, PUBLIC, EXTENDS(AbstractMOCalc) :: MOCalcFilterMatrix
   PRIVATE

   TYPE(qs_environment_type), POINTER         :: qs_env => null()
   TYPE(qs_scf_env_type), POINTER             :: scf_env => null()
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_ks => null(), matrix_s => null()
   TYPE(section_vals_type), POINTER           :: scf_section

CONTAINS
   PROCEDURE, PUBLIC, PASS(self) :: init => filter_matrix_init
   PROCEDURE, PUBLIC, PASS(self) :: pre => filter_matrix_pre
   PROCEDURE, PUBLIC, PASS(self) :: run => filter_matrix_run
END TYPE

CONTAINS
SUBROUTINE filter_matrix_init(self, qs_env, scf_env, scf_section)
   CLASS(MOCalcFilterMatrix), INTENT(INOUT)             :: self
   TYPE(qs_environment_type), POINTER, INTENT(IN) :: qs_env
   TYPE(qs_scf_env_type), POINTER, INTENT(IN)     :: scf_env
   TYPE(section_vals_type), POINTER, INTENT(IN)   :: scf_section
   TYPE(scf_control_type), POINTER            :: scf_control => null()

   self%qs_env => qs_env
   self%scf_env => scf_env
   self%scf_section => scf_section

   CALL get_qs_env(qs_env, &
                   scf_control=scf_control, &
                   matrix_ks=self%matrix_ks, &
                   matrix_s=self%matrix_s)

   IF (.NOT. scf_env%skip_diis) THEN
      IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
         CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=scf_control%max_diis)
      END IF
   END IF
END SUBROUTINE filter_matrix_init

SUBROUTINE filter_matrix_pre(self)
   CLASS(MOCalcFilterMatrix), INTENT(INOUT) :: self

   CALL qs_diis_b_clear(self%scf_env%scf_diis_buffer)
END SUBROUTINE filter_matrix_pre

SUBROUTINE filter_matrix_run(self, diis_step, energy_only)
   !! the inner loop of scf, specific to iterative diagonalization using OT
   !! with S matrix; basically, in goes the ks matrix out goes a new p matrix
   CLASS(MOCalcFilterMatrix), INTENT(INOUT)   :: self
   LOGICAL, INTENT(INOUT)               :: diis_step
   LOGICAL, INTENT(INOUT)               :: energy_only

   MARK_USED(energy_only)

   CALL fb_env_do_diag(self%scf_env%filter_matrix_env, self%qs_env, &
                       self%matrix_ks, self%matrix_s, self%scf_section, diis_step)
END SUBROUTINE filter_matrix_run

END MODULE qs_scf_mo_calc_filter_matrix
