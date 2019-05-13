!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_mo_calc_block_davidson
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE qs_block_davidson_types,         ONLY: block_davidson_deallocate
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: mo_set_p_type
   USE qs_scf_abstract_mo_calc,         ONLY: AbstractMOCalc
   USE qs_scf_diagonalization,          ONLY: do_block_davidson_diag
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! CASE (block_davidson_diag_method_nr, general_diag_method_nr, special_diag_method_nr, block_krylov_diag_method_nr)

   TYPE, PUBLIC, EXTENDS(AbstractMOCalc) :: MOCalcBlockDavidson
   PRIVATE

   TYPE(qs_environment_type), POINTER         :: qs_env => null()
   TYPE(qs_scf_env_type), POINTER             :: scf_env => null()
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_ks => null(), matrix_s => null()
   TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos => null()
   TYPE(section_vals_type), POINTER           :: scf_section => null()
   TYPE(scf_control_type), POINTER            :: scf_control => null()

CONTAINS
   PROCEDURE, PUBLIC, PASS(self) :: init => block_davidson_init
   PROCEDURE, PUBLIC, PASS(self) :: pre => block_davidson_pre
   PROCEDURE, PUBLIC, PASS(self) :: run => block_davidson_run
   FINAL :: block_davidson_finalize
END TYPE

CONTAINS
SUBROUTINE block_davidson_init(self, qs_env, scf_env, scf_section)
   CLASS(MOCalcBlockDavidson), INTENT(INOUT)             :: self
   TYPE(qs_environment_type), POINTER, INTENT(IN) :: qs_env
   TYPE(qs_scf_env_type), POINTER, INTENT(IN)     :: scf_env
   TYPE(section_vals_type), POINTER, INTENT(IN)   :: scf_section

   self%qs_env => qs_env
   self%scf_env => scf_env
   self%scf_section => scf_section

   CALL get_qs_env(qs_env, &
                   scf_control=self%scf_control, &
                   mos=self%mos, &
                   matrix_ks=self%matrix_ks, &
                   matrix_s=self%matrix_s)
END SUBROUTINE block_davidson_init

SUBROUTINE block_davidson_pre(self)
   CLASS(MOCalcBlockDavidson), INTENT(INOUT) :: self

   MARK_USED(self)

END SUBROUTINE block_davidson_pre

SUBROUTINE block_davidson_run(self, diis_step, energy_only)
   !! the inner loop of scf, specific to iterative diagonalization using OT
   !! with S matrix; basically, in goes the ks matrix out goes a new p matrix
   CLASS(MOCalcBlockDavidson), INTENT(INOUT)   :: self
   LOGICAL, INTENT(INOUT)               :: diis_step
   LOGICAL, INTENT(INOUT)               :: energy_only

   MARK_USED(diis_step)
   MARK_USED(energy_only)

   CALL do_block_davidson_diag(self%qs_env, self%scf_env, self%mos, self%matrix_ks, self%matrix_s, self%scf_control, &
                               self%scf_section, .FALSE.)
END SUBROUTINE block_davidson_run

SUBROUTINE block_davidson_finalize(self)
      TYPE(MOCalcBlockDavidson), INTENT(INOUT)           :: self

   CALL block_davidson_deallocate(self%scf_env%block_davidson_env)
END SUBROUTINE

END MODULE qs_scf_mo_calc_block_davidson
