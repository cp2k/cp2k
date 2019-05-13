!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_mo_calc_general_diag
   USE cp_control_types,                ONLY: dft_control_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_create
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_mom_methods,                  ONLY: do_mom_diag
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_scf_abstract_mo_calc,         ONLY: AbstractMOCalc
   USE qs_scf_diagonalization,          ONLY: do_general_diag,&
                                              do_roks_diag,&
                                              do_scf_diag_subspace
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC, EXTENDS(AbstractMOCalc) :: MOCalcGeneralDiag
   PRIVATE

   TYPE(qs_environment_type), POINTER         :: qs_env => null()
   TYPE(dft_control_type), POINTER            :: dft_control => null()
   TYPE(qs_scf_env_type), POINTER             :: scf_env => null()
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_ks => null(), matrix_s => null()
   TYPE(section_vals_type), POINTER           :: scf_section
   TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos => null()
   TYPE(qs_rho_type), POINTER                 :: rho => null()
   TYPE(qs_ks_env_type), POINTER              :: ks_env => null()
   TYPE(scf_control_type), POINTER            :: scf_control => null()
   LOGICAL :: has_unit_metric

CONTAINS
   PROCEDURE, PUBLIC, PASS(self) :: init => general_diag_init
   PROCEDURE, PUBLIC, PASS(self) :: pre => general_diag_pre
   PROCEDURE, PUBLIC, PASS(self) :: run => general_diag_run
END TYPE

CONTAINS
SUBROUTINE general_diag_init(self, qs_env, scf_env, scf_section)
   CLASS(MOCalcGeneralDiag), INTENT(INOUT)             :: self
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
                   ks_env=self%ks_env, &
                   rho=self%rho, &
                   has_unit_metric=self%has_unit_metric)

   IF (.NOT. scf_env%skip_diis) THEN
      IF (.NOT. ASSOCIATED(scf_env%scf_diis_buffer)) THEN
         CALL qs_diis_b_create(scf_env%scf_diis_buffer, nbuffer=self%scf_control%max_diis)
      END IF
   END IF
END SUBROUTINE general_diag_init

SUBROUTINE general_diag_pre(self)
   CLASS(MOCalcGeneralDiag), INTENT(INOUT) :: self

   CALL qs_diis_b_clear(self%scf_env%scf_diis_buffer)
END SUBROUTINE general_diag_pre

SUBROUTINE general_diag_run(self, diis_step, energy_only)
   CLASS(MOCalcGeneralDiag), INTENT(INOUT)   :: self
   LOGICAL, INTENT(INOUT)               :: diis_step
   LOGICAL, INTENT(INOUT)               :: energy_only
   LOGICAL :: skip_diag_sub

   MARK_USED(energy_only)

   ASSOCIATE (scf_env=>self%scf_env, &
              scf_control=>self%scf_control, &
              scf_section=>self%scf_section, &
              matrix_s=>self%matrix_s, &
              matrix_ks=>self%matrix_ks, &
              mos=>self%mos)

      IF (self%dft_control%roks) THEN
         CALL do_roks_diag(self%scf_env, self%mos, self%matrix_ks, self%matrix_s, &
                           self%scf_control, self%scf_section, diis_step, &
                           self%has_unit_metric)
      ELSE
         IF (self%scf_control%diagonalization%mom) THEN
            CALL do_mom_diag(self%scf_env, self%mos, self%matrix_ks, &
                             self%matrix_s, self%scf_control, self%scf_section, &
                             diis_step)
         ELSE
            CALL do_general_diag(self%scf_env, self%mos, self%matrix_ks, &
                                 self%matrix_s, self%scf_control, self%scf_section, &
                                 diis_step)
         END IF
         IF (self%scf_control%do_diag_sub) THEN
            skip_diag_sub = (scf_env%subspace_env%eps_diag_sub > 0.0_dp) .AND. &
                            (scf_env%iter_count == 1 .OR. scf_env%iter_delta > scf_env%subspace_env%eps_diag_sub)
            IF (.NOT. skip_diag_sub) THEN
               CALL do_scf_diag_subspace(self%qs_env, self%scf_env, scf_env%subspace_env, self%mos, self%rho, &
                                         self%ks_env, self%scf_section, self%scf_control)
            END IF
         END IF
      END IF
   END ASSOCIATE
END SUBROUTINE general_diag_run

END MODULE qs_scf_mo_calc_general_diag
