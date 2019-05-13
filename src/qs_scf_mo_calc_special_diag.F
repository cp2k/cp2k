!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_mo_calc_special_diag
   USE cp_control_types,                ONLY: dft_control_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_create
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_scf_abstract_mo_calc,         ONLY: AbstractMOCalc
   USE qs_scf_diagonalization,          ONLY: do_roks_diag,&
                                              do_special_diag
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! CASE (special_diag_diag_method_nr, general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr)

   TYPE, PUBLIC, EXTENDS(AbstractMOCalc) :: MOCalcSpecialDiag
   PRIVATE

   TYPE(qs_environment_type), POINTER         :: qs_env => null()
   TYPE(qs_scf_env_type), POINTER             :: scf_env => null()
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_ks => null(), matrix_s => null()
   TYPE(section_vals_type), POINTER           :: scf_section => null()
   TYPE(scf_control_type), POINTER            :: scf_control => null()
   TYPE(dft_control_type), POINTER            :: dft_control => null()
   TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos => null()
   LOGICAL :: has_unit_metric

CONTAINS
   PROCEDURE, PUBLIC, PASS(self) :: init => special_diag_init
   PROCEDURE, PUBLIC, PASS(self) :: pre => special_diag_pre
   PROCEDURE, PUBLIC, PASS(self) :: run => special_diag_run
END TYPE

CONTAINS
SUBROUTINE special_diag_init(self, qs_env, scf_env, scf_section)
   CLASS(MOCalcSpecialDiag), INTENT(INOUT)             :: self
   TYPE(qs_environment_type), POINTER, INTENT(IN) :: qs_env
   TYPE(qs_scf_env_type), POINTER, INTENT(IN)     :: scf_env
   TYPE(section_vals_type), POINTER, INTENT(IN)   :: scf_section

   self%qs_env => qs_env
   self%scf_env => scf_env
   self%scf_section => scf_section

   CALL get_qs_env(qs_env, &
                   scf_control=self%scf_control, &
                   dft_control=self%dft_control, &
                   mos=self%mos, &
                   matrix_ks=self%matrix_ks, &
                   matrix_s=self%matrix_s, &
                   has_unit_metric=self%has_unit_metric)

   IF (.NOT. scf_env%skip_diis) THEN
      IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
         CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=self%scf_control%max_diis)
      END IF
   END IF
END SUBROUTINE special_diag_init

SUBROUTINE special_diag_pre(self)
   CLASS(MOCalcSpecialDiag), INTENT(INOUT) :: self

   CALL qs_diis_b_clear(self%scf_env%scf_diis_buffer)
END SUBROUTINE special_diag_pre

SUBROUTINE special_diag_run(self, diis_step, energy_only)
   !! the inner loop of scf, specific to iterative diagonalization using OT
   !! with S matrix; basically, in goes the ks matrix out goes a new p matrix
   CLASS(MOCalcSpecialDiag), INTENT(INOUT)   :: self
   LOGICAL, INTENT(INOUT)               :: diis_step
   LOGICAL, INTENT(INOUT)               :: energy_only

   MARK_USED(energy_only)

   IF (self%dft_control%roks) THEN
      CALL do_roks_diag(self%scf_env, self%mos, self%matrix_ks, self%matrix_s, &
                        self%scf_control, self%scf_section, diis_step, &
                        self%has_unit_metric)
   ELSE
      CALL do_special_diag(self%scf_env, self%mos, self%matrix_ks, &
                           self%scf_control, self%scf_section, &
                           diis_step)
   END IF
END SUBROUTINE special_diag_run

END MODULE qs_scf_mo_calc_special_diag
