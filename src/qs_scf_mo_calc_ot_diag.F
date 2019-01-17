!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_mo_calc_ot_diag
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE preconditioner,                  ONLY: prepare_preconditioner,&
                                              restart_preconditioner
   USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                              qs_diis_b_create,&
                                              qs_diis_b_step
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                              calculate_subspace_eigenvalues
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_scf_abstract_mo_calc,         ONLY: AbstractMOCalc
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC, EXTENDS(AbstractMOCalc) :: MOCalcOTDiag
   PRIVATE

   TYPE(qs_environment_type), POINTER         :: qs_env => null()
   TYPE(dft_control_type), POINTER            :: dft_control => null()
   TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos => null()
   TYPE(scf_control_type), POINTER            :: scf_control => null()
   TYPE(qs_scf_env_type), POINTER             :: scf_env => null()
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_ks => null(), matrix_s => null()
   TYPE(section_vals_type), POINTER           :: scf_section

CONTAINS
   PROCEDURE, PUBLIC, PASS(self) :: init => ot_diag_init
   PROCEDURE, PUBLIC, PASS(self) :: pre => ot_diag_pre
   PROCEDURE, PUBLIC, PASS(self) :: run => ot_diag_run
END TYPE

CONTAINS
SUBROUTINE ot_diag_init(self, qs_env, scf_env, scf_section)
   CLASS(MOCalcOTDiag), INTENT(INOUT)             :: self
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
                   matrix_s=self%matrix_s)

   ! disable DFTB, SE, xTB for now
   IF (self%dft_control%qs_control%dftb .OR. &
       self%dft_control%qs_control%xtb .OR. &
       self%dft_control%qs_control%semi_empirical) THEN
      CPABORT("DFTB and SE not available with OT/DIAG")
   END IF
END SUBROUTINE ot_diag_init

SUBROUTINE ot_diag_pre(self)
   CLASS(MOCalcOTDiag), INTENT(INOUT) :: self

   IF (.NOT. self%scf_env%skip_diis) THEN
      IF (.NOT. ASSOCIATED(self%scf_env%scf_diis_buffer)) THEN
         CALL qs_diis_b_create(self%scf_env%scf_diis_buffer, nbuffer=self%scf_control%max_diis)
      END IF
      CALL qs_diis_b_clear(self%scf_env%scf_diis_buffer)
   END IF

   ASSOCIATE (settings=>self%scf_control%diagonalization%ot_settings)

      ! if an old preconditioner is still around (i.e. outer SCF is active),
      ! remove it if this could be worthwhile
      CALL restart_preconditioner(self%qs_env, self%scf_env%ot_preconditioner, &
                                  settings%preconditioner_type, &
                                  self%dft_control%nspins)

      CALL prepare_preconditioner(self%qs_env, self%mos, self%matrix_ks, self%matrix_s, &
                                  self%scf_env%ot_preconditioner, &
                                  settings%preconditioner_type, &
                                  settings%precond_solver_type, &
                                  settings%energy_gap, &
                                  self%dft_control%nspins)

   END ASSOCIATE
END SUBROUTINE ot_diag_pre

SUBROUTINE ot_diag_run(self, diis_step, energy_only)
   !! the inner loop of scf, specific to iterative diagonalization using OT
   !! with S matrix; basically, in goes the ks matrix out goes a new p matrix
   CLASS(MOCalcOTDiag), INTENT(INOUT)   :: self
   LOGICAL, INTENT(INOUT)               :: diis_step
   LOGICAL, INTENT(INOUT)               :: energy_only
   INTEGER                              :: homo, ispin, nmo, nspin
   REAL(KIND=dp)                        :: diis_error, eps_iter
   REAL(KIND=dp), DIMENSION(:), POINTER :: eigenvalues
   TYPE(cp_fm_type), POINTER            :: mo_coeff

   NULLIFY (eigenvalues)
   MARK_USED(energy_only)

   ! note: associate'd pointers get the TARGET attribute, making it impossible
   !       to pass them to our functions without changing the function definition.
   !       we can therefore only alias members we are not passing along directly.
   ASSOCIATE (scf_env=>self%scf_env, &
              scf_control=>self%scf_control, &
              scf_section=>self%scf_section, &
              matrix_ks=>self%matrix_ks)

      nspin = SIZE(matrix_ks)

      DO ispin = 1, nspin
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, &
                               scf_env%scf_work1(ispin)%matrix)
      END DO

      IF ((scf_env%iter_count > 1) .AND. (.NOT. self%scf_env%skip_diis)) THEN
         CALL qs_diis_b_step(scf_env%scf_diis_buffer, self%mos, scf_env%scf_work1, &
                             scf_env%scf_work2, scf_env%iter_delta, diis_error, diis_step, &
                             scf_control%eps_diis, scf_control%nmixing, &
                             s_matrix=self%matrix_s, &
                             scf_section=self%scf_section)
      ELSE
         diis_step = .FALSE.
      END IF

      eps_iter = scf_control%diagonalization%eps_iter

      IF (diis_step) THEN
         scf_env%iter_param = diis_error
         scf_env%iter_method = "DIIS/OTdiag"
         DO ispin = 1, nspin
            CALL copy_fm_to_dbcsr(scf_env%scf_work1(ispin)%matrix, &
                                  matrix_ks(ispin)%matrix, keep_sparsity=.TRUE.)
         END DO
         eps_iter = MAX(eps_iter, scf_control%diagonalization%eps_adapt*diis_error)
      ELSE
         IF (scf_env%mixing_method == 1) THEN
            scf_env%iter_param = scf_env%p_mix_alpha
            scf_env%iter_method = "P_Mix/OTdiag."
         ELSEIF (scf_env%mixing_method > 1) THEN
            scf_env%iter_param = scf_env%mixing_store%alpha
            scf_env%iter_method = TRIM(scf_env%mixing_store%iter_method)//"/OTdiag."
         END IF
      END IF

      scf_env%iter_delta = 0.0_dp

      DO ispin = 1, nspin
         CALL get_mo_set(self%mos(ispin)%mo_set, &
                         mo_coeff=mo_coeff, &
                         eigenvalues=eigenvalues, &
                         nmo=nmo, &
                         homo=homo)
         CALL ot_eigensolver(matrix_h=matrix_ks(ispin)%matrix, &
                             matrix_s=self%matrix_s(1)%matrix, &
                             matrix_c_fm=mo_coeff, &
                             preconditioner=scf_env%ot_preconditioner(1)%preconditioner, &
                             eps_gradient=eps_iter, &
                             iter_max=scf_control%diagonalization%max_iter, &
                             silent=.TRUE., &
                             ot_settings=scf_control%diagonalization%ot_settings)
         CALL calculate_subspace_eigenvalues(mo_coeff, matrix_ks(ispin)%matrix, &
                                             evals_arg=eigenvalues, &
                                             do_rotation=.TRUE.)
         !MK WRITE(*,*) routinen//' copy_dbcsr_to_fm'
         CALL copy_fm_to_dbcsr(self%mos(ispin)%mo_set%mo_coeff, &
                               self%mos(ispin)%mo_set%mo_coeff_b)
         !fm->dbcsr
      END DO

      CALL set_mo_occupation(mo_array=self%mos, smear=self%scf_control%smear)

      DO ispin = 1, nspin
         ! does not yet handle k-points
         CALL calculate_density_matrix(self%mos(ispin)%mo_set, &
                                       scf_env%p_mix_new(ispin, 1)%matrix)
      END DO

   END ASSOCIATE
END SUBROUTINE ot_diag_run

END MODULE qs_scf_mo_calc_ot_diag
