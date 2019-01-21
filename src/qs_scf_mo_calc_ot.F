!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) CP2K developers group                                                            !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_mo_calc_ot
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: ot_precond_full_all,&
                                              ot_precond_full_single,&
                                              ot_precond_full_single_inverse,&
                                              ot_precond_none,&
                                              ot_precond_s_inverse
   USE input_section_types,             ONLY: section_vals_type
   USE preconditioner,                  ONLY: prepare_preconditioner,&
                                              restart_preconditioner
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                              make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_ot,                           ONLY: qs_ot_new_preconditioner
   USE qs_ot_scf,                       ONLY: ot_scf_destroy,&
                                              ot_scf_init,&
                                              ot_scf_mini,&
                                              ot_scf_read_input
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_abstract_mo_calc,         ONLY: AbstractMOCalc
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type,&
                                              smear_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC, EXTENDS(AbstractMOCalc) :: MOCalcOT
   PRIVATE

   TYPE(qs_environment_type), POINTER         :: qs_env => null()
   TYPE(dft_control_type), POINTER            :: dft_control => null()
   TYPE(mo_set_p_type), DIMENSION(:), POINTER :: mos => null()
   TYPE(scf_control_type), POINTER            :: scf_control => null()
   TYPE(qs_scf_env_type), POINTER             :: scf_env => null()
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER  :: matrix_ks => null(), matrix_s => null()
   TYPE(section_vals_type), POINTER           :: scf_section

CONTAINS
   PROCEDURE, PUBLIC, PASS(self) :: init => ot_init
   PROCEDURE, PUBLIC, PASS(self) :: pre => ot_pre
   PROCEDURE, PUBLIC, PASS(self) :: run => ot_run
   FINAL :: ot_finalize
END TYPE

CONTAINS
SUBROUTINE ot_init(self, qs_env, scf_env, scf_section)
   CLASS(MOCalcOT), INTENT(INOUT)                 :: self
   TYPE(qs_environment_type), POINTER, INTENT(IN) :: qs_env
   TYPE(qs_scf_env_type), POINTER, INTENT(IN)     :: scf_env
   TYPE(section_vals_type), POINTER, INTENT(IN)   :: scf_section
   LOGICAL                                        :: has_unit_metric
   INTEGER                                        :: ispin, nmo, number_of_OT_envs
   TYPE(cp_fm_type), POINTER                      :: mo_coeff

   self%qs_env => qs_env
   self%scf_env => scf_env
   self%scf_section => scf_section

   CALL get_qs_env(qs_env, &
                   scf_control=self%scf_control, &
                   dft_control=self%dft_control, &
                   mos=self%mos, &
                   matrix_s=self%matrix_s, &
                   matrix_ks=self%matrix_ks, &
                   has_unit_metric=has_unit_metric)

   ASSOCIATE (scf_control=>self%scf_control, &
              dft_control=>self%dft_control)

      ! reortho the wavefunctions if we are having an outer scf and
      ! this is not the first iteration
      ! this is useful to avoid the build-up of numerical noise
      ! however, we can not play this trick if restricted (don't mix non-equivalent orbs)
      IF (scf_control%do_outer_scf_reortho) THEN
         IF (scf_control%outer_scf%have_scf .AND. .NOT. dft_control%restricted) THEN
            IF (scf_env%outer_scf%iter_count > 0) THEN
               DO ispin = 1, dft_control%nspins
                  CALL get_mo_set(mo_set=self%mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
                  IF (has_unit_metric) THEN
                     CALL make_basis_simple(mo_coeff, nmo)
                  ELSE
                     CALL make_basis_sm(mo_coeff, nmo, self%matrix_s(1)%matrix)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSE
         ! dont need any dirty trick for the numerically stable irac algorithm.
      ENDIF

      ! restricted calculations require just one set of OT orbitals
      number_of_OT_envs = dft_control%nspins
      IF (dft_control%restricted) number_of_OT_envs = 1

      ALLOCATE (scf_env%qs_ot_env(number_of_OT_envs))

      ! XXX Joost XXX should disentangle reading input from this part
      CALL ot_scf_read_input(scf_env%qs_ot_env, scf_section)

      ! keep a note that we are restricted
      IF (dft_control%restricted) THEN
         scf_env%qs_ot_env(1)%restricted = .TRUE.
         ! requires rotation
         IF (.NOT. scf_env%qs_ot_env(1)%settings%do_rotation) &
            CALL cp_abort(__LOCATION__, &
                          "Restricted calculation with OT requires orbital rotation. Please "// &
                          "activate the OT%ROTATION keyword!")
      ELSE
         scf_env%qs_ot_env(:)%restricted = .FALSE.
      ENDIF

      ! might need the KS matrix to init properly
      CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., &
                               calculate_forces=.FALSE.)

   END ASSOCIATE
END SUBROUTINE ot_init

SUBROUTINE ot_pre(self)
   CLASS(MOCalcOT), INTENT(INOUT) :: self

   TYPE(dbcsr_type), POINTER                      :: orthogonality_metric
   LOGICAL                                        :: do_rotation, has_unit_metric, is_full_all
   INTEGER                                        :: ispin

   CALL get_qs_env(self%qs_env, &
                   has_unit_metric=has_unit_metric)

   ASSOCIATE (scf_control=>self%scf_control, &
              dft_control=>self%dft_control, &
              scf_env=>self%scf_env)

      do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation

      ! another safety check
      IF (dft_control%low_spin_roks) THEN
         CPASSERT(do_rotation)
      ENDIF

      ! if an old preconditioner is still around (i.e. outer SCF is active),
      ! remove it if this could be worthwhile
      CALL restart_preconditioner(self%qs_env, self%scf_env%ot_preconditioner, &
                                  self%scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                  dft_control%nspins)

      !
      ! preconditioning still needs to be done correctly with has_unit_metric
      ! notice that a big part of the preconditioning (S^-1) is fine anyhow
      !
      IF (has_unit_metric) THEN
         NULLIFY (orthogonality_metric)
      ELSE
         orthogonality_metric => self%matrix_s(1)%matrix
      ENDIF

      CALL prepare_preconditioner(self%qs_env, self%mos, self%matrix_ks, self%matrix_s, scf_env%ot_preconditioner, &
                                  scf_env%qs_ot_env(1)%settings%preconditioner_type, &
                                  scf_env%qs_ot_env(1)%settings%precond_solver_type, &
                                  scf_env%qs_ot_env(1)%settings%energy_gap, dft_control%nspins, &
                                  has_unit_metric=has_unit_metric, &
                                  chol_type=scf_env%qs_ot_env(1)%settings%cholesky_type)

      CALL ot_scf_init(mo_array=self%mos, matrix_s=orthogonality_metric, &
                       broyden_adaptive_sigma=self%qs_env%broyden_adaptive_sigma, &
                       qs_ot_env=scf_env%qs_ot_env, matrix_ks=self%matrix_ks(1)%matrix)

      SELECT CASE (scf_env%qs_ot_env (1)%settings%preconditioner_type)
      CASE (ot_precond_none)
      CASE (ot_precond_full_all, ot_precond_full_single_inverse)
         ! this will rotate the MOs to be eigen states, which is not compatible with rotation
         ! e.g. mo_derivs here do not yet include potentially different occupations numbers
         ! only full all needs rotation
         is_full_all = scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all
         CPASSERT(.NOT. (do_rotation .AND. is_full_all))
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                          scf_env%ot_preconditioner(ispin)%preconditioner)
         ENDDO
      CASE (ot_precond_s_inverse, ot_precond_full_single)
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                          scf_env%ot_preconditioner(1)%preconditioner)
         ENDDO
      CASE DEFAULT
         DO ispin = 1, SIZE(scf_env%qs_ot_env)
            CALL qs_ot_new_preconditioner(scf_env%qs_ot_env(ispin), &
                                          scf_env%ot_preconditioner(1)%preconditioner)
         ENDDO
      END SELECT

      ! if we have non-uniform occupations we should be using rotation
      do_rotation = scf_env%qs_ot_env(1)%settings%do_rotation
      DO ispin = 1, SIZE(self%mos)
         IF (.NOT. self%mos(ispin)%mo_set%uniform_occupation) THEN
            CPASSERT(do_rotation)
         ENDIF
      ENDDO

   END ASSOCIATE
END SUBROUTINE ot_pre

SUBROUTINE ot_finalize(self)
      TYPE(MOCalcOT)                                     :: self

      INTEGER                                            :: ispin

   DO ispin = 1, SIZE(self%scf_env%qs_ot_env)
      CALL ot_scf_destroy(self%scf_env%qs_ot_env(ispin))
   ENDDO

   DEALLOCATE (self%scf_env%qs_ot_env)
END SUBROUTINE

SUBROUTINE ot_run(self, diis_step, energy_only)
   !! The inner loop of scf, specific to using to the orbital transformation method
   !! basically, in goes the ks matrix out goes a new p matrix
   CLASS(MOCalcOT), INTENT(INOUT) :: self
   LOGICAL, INTENT(INOUT)         :: diis_step
   LOGICAL, INTENT(INOUT)                             :: energy_only

   TYPE(smear_type), POINTER                          :: smear
   TYPE(qs_rho_type), POINTER                         :: rho
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
   LOGICAL                                            :: has_unit_metric
   INTEGER                                            :: ispin
   TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
   TYPE(dbcsr_type), POINTER                          :: orthogonality_metric
   TYPE(qs_energy_type), POINTER                      :: energy
   NULLIFY (rho_ao)

   MARK_USED(diis_step) ! should diis_step be removed?!

   CALL qs_rho_get(rho, rho_ao=rho_ao)
   CALL get_qs_env(self%qs_env, &
                   energy=energy, &
                   has_unit_metric=has_unit_metric) ! store this as a member variable?!

   IF (has_unit_metric) THEN
      NULLIFY (orthogonality_metric)
   ELSE
      orthogonality_metric => matrix_s(1)%matrix
   END IF

   ! in case of LSD the first spin qs_ot_env will drive the minimization
   ! in the case of a restricted calculation, it will make sure the spin orbitals are equal

   CALL ot_scf_mini(self%mos, mo_derivs, smear, orthogonality_metric, &
                    energy%total, energy_only, self%scf_env%iter_delta, &
                    self%scf_env%qs_ot_env, self%qs_env%input)

   DO ispin = 1, SIZE(self%mos)
      CALL set_mo_occupation(mo_set=self%mos(ispin)%mo_set, &
                             smear=smear)
   ENDDO

   DO ispin = 1, SIZE(self%mos)
      CALL calculate_density_matrix(self%mos(ispin)%mo_set, &
                                    rho_ao(ispin)%matrix, &
                                    use_dbcsr=.TRUE.)
   END DO

   self%scf_env%iter_method = self%scf_env%qs_ot_env(1)%OT_METHOD_FULL
   self%scf_env%iter_param = self%scf_env%qs_ot_env(1)%ds_min
   self%qs_env%broyden_adaptive_sigma = self%scf_env%qs_ot_env(1)%broyden_adaptive_sigma

END SUBROUTINE ot_run

END MODULE qs_scf_mo_calc_ot
