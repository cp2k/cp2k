!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! *****************************************************************************
!> \brief ...
!> \author ...
! *****************************************************************************
MODULE qs_basis_gradient

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_copy,&
                                              dbcsr_init,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
   USE qs_energy_matrix_w,              ONLY: qs_energies_compute_w
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              deallocate_qs_force,&
                                              qs_force_type,&
                                              replicate_qs_force,&
                                              zero_qs_force
   USE qs_gspace_mixing,                ONLY: mixing_allocate
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_methods,                   ONLY: qs_ks_allocate_basics,&
                                              qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                              make_basis_sm
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_subsys_types,                 ONLY: qs_subsys_set,&
                                              qs_subsys_type
   USE qs_update_s_mstruct,             ONLY: qs_env_update_s_mstruct
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_basis_gradient'

! *** Public subroutines ***

   PUBLIC :: qs_basis_center_gradient, qs_update_basis_center_pos, &
             return_basis_center_gradient_norm

CONTAINS

! *****************************************************************************
! for development of floating basis functions
! *****************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE qs_basis_center_gradient(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_basis_center_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iatom, ikind, img, ispin, &
                                                            natom, nimg, nkind, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, natom_of_kind
      LOGICAL                                            :: floating, has_unit_metric
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_w_kp
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: basis_force, force, qs_force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)
      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! get the gradient array
      CALL get_qs_env(qs_env, scf_env=scf_env, natom=natom)
      IF (ASSOCIATED(scf_env%floating_basis%gradient)) THEN
         gradient => scf_env%floating_basis%gradient
         CPASSERT(SIZE(gradient) == 3*natom)
      ELSE
         ALLOCATE (gradient(3, natom))
         scf_env%floating_basis%gradient => gradient
      END IF
      gradient = 0.0_dp

      ! init the force environment
      CALL get_qs_env(qs_env, force=force, subsys=subsys, atomic_kind_set=atomic_kind_set)
      IF (ASSOCIATED(force)) THEN
         qs_force => force
      ELSE
         NULLIFY (qs_force)
      END IF
      ! Allocate the force data structure
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
      CALL allocate_qs_force(basis_force, natom_of_kind)
      DEALLOCATE (natom_of_kind)
      CALL qs_subsys_set(subsys, force=basis_force)
      CALL zero_qs_force(basis_force)

      ! get atom mapping
      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)

      ! allocate energy weighted density matrices, if needed
      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)
      IF (.NOT. has_unit_metric) THEN
         CALL get_qs_env(qs_env, matrix_w_kp=matrix_w_kp)
         IF (.NOT. ASSOCIATED(matrix_w_kp)) THEN
            NULLIFY (matrix_w_kp)
            CALL get_qs_env(qs_env, ks_env=ks_env, matrix_s_kp=matrix_s, dft_control=dft_control)
            nspin = dft_control%nspins
            nimg = dft_control%nimages
            matrix => matrix_s(1, 1)%matrix
            CALL dbcsr_allocate_matrix_set(matrix_w_kp, nspin, nimg)
            DO ispin = 1, nspin
               DO img = 1, nimg
                  ALLOCATE (matrix_w_kp(ispin, img)%matrix)
                  CALL dbcsr_init(matrix_w_kp(ispin, img)%matrix)
                  CALL dbcsr_copy(matrix_w_kp(ispin, img)%matrix, matrix, name="W MATRIX")
                  CALL dbcsr_set(matrix_w_kp(ispin, img)%matrix, 0.0_dp)
               END DO
            END DO
            CALL set_ks_env(ks_env, matrix_w_kp=matrix_w_kp)
         ENDIF
      ENDIF
      ! time to compute the w matrix
      CALL qs_energies_compute_w(qs_env, .TRUE.)

      ! core hamiltonian forces
      CALL build_core_hamiltonian_matrix(qs_env=qs_env, calculate_forces=.TRUE.)
      ! Compute grid-based forces
      CALL qs_ks_update_qs_env(qs_env, calculate_forces=.TRUE.)

      ! replicate forces
      CALL replicate_qs_force(basis_force, para_env)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         i = atom_of_kind(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), floating=floating)
         IF (floating) gradient(1:3, iatom) = -basis_force(ikind)%total(1:3, i)
      END DO
      ! clean up force environment and reinitialize qs_force
      IF (ASSOCIATED(basis_force)) CALL deallocate_qs_force(basis_force)
      CALL qs_subsys_set(subsys, force=qs_force)
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL set_ks_env(ks_env, forces_up_to_date=.FALSE.)

      DEALLOCATE (atom_of_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE qs_basis_center_gradient

! *****************************************************************************
!> \brief ... returns the norm of the gradient vector, taking only floating
!>             components into account
!> \param qs_env ...
!> \retval norm ...
! **************************************************************************************************
   FUNCTION return_basis_center_gradient_norm(qs_env) RESULT(norm)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp)                                      :: norm

      INTEGER                                            :: iatom, ikind, natom, nfloat
      LOGICAL                                            :: floating
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      norm = 0.0_dp
      CALL get_qs_env(qs_env, scf_env=scf_env, particle_set=particle_set, qs_kind_set=qs_kind_set)
      gradient => scf_env%floating_basis%gradient
      natom = SIZE(particle_set)
      nfloat = 0
      DO iatom = 1, natom
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), floating=floating)
         IF (floating) THEN
            nfloat = nfloat+1
            norm = norm+SUM(ABS(gradient(1:3, iatom)))
         END IF
      END DO
      IF (nfloat > 0) THEN
         norm = norm/(3.0_dp*REAL(nfloat, KIND=dp))
      END IF

   END FUNCTION return_basis_center_gradient_norm

! *****************************************************************************
!> \brief move atoms with kind float according to gradient
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE qs_update_basis_center_pos(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_update_basis_center_pos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ikind, natom
      LOGICAL                                            :: floating
      REAL(KIND=dp)                                      :: alpha
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)
      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! update positions
      CALL get_qs_env(qs_env, scf_env=scf_env, particle_set=particle_set, qs_kind_set=qs_kind_set)
      gradient => scf_env%floating_basis%gradient
      natom = SIZE(particle_set)
      alpha = 0.50_dp
      DO iatom = 1, natom
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), floating=floating)
         IF (floating) THEN
            particle_set(iatom)%r(1:3) = particle_set(iatom)%r(1:3)+alpha*gradient(1:3, iatom)
         END IF
      END DO

      CALL qs_basis_reinit_energy(qs_env)

      CALL timestop(handle)

   END SUBROUTINE qs_update_basis_center_pos

! *****************************************************************************
!> \brief rebuilds the structures after a floating basis update
!> \param qs_env ...
!> \par History
!>      05.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_basis_reinit_energy(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_basis_reinit_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nmo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ! rebuild neighbor lists
      CALL get_qs_env(qs_env, input=input, para_env=para_env)
      CALL build_qs_neighbor_lists(qs_env, para_env, molecular=.FALSE., &
                                   force_env_section=input)
      ! update core hamiltonian
      CALL build_core_hamiltonian_matrix(qs_env=qs_env, calculate_forces=.FALSE.)
      ! update structures
      CALL qs_env_update_s_mstruct(qs_env)
      ! KS matrices
      CALL qs_ks_allocate_basics(qs_env)
      ! reinit SCF task matrices
      NULLIFY (scf_env)
      CALL get_qs_env(qs_env, scf_env=scf_env, dft_control=dft_control)
      IF (scf_env%mixing_method > 0) THEN
         CALL mixing_allocate(qs_env, scf_env%mixing_method, scf_env%p_mix_new, &
                              scf_env%p_delta, dft_control%nspins, &
                              scf_env%mixing_store)
      ELSE
         NULLIFY (scf_env%p_mix_new)
      END IF
      CALL get_qs_env(qs_env, mos=mos, rho=rho, matrix_s_kp=matrix_s_kp)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      DO ispin = 1, SIZE(mos)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
         ! reorthogonalize MOs
         CALL make_basis_sm(mo_coeff, nmo, matrix_s_kp(1, 1)%matrix)
         ! update density matrix
         CALL calculate_density_matrix(mos(ispin)%mo_set, rho_ao_kp(ispin, 1)%matrix)
      END DO
      CALL qs_rho_set(rho, rho_r_valid=.FALSE., drho_r_valid=.FALSE., rho_g_valid=.FALSE., &
                      drho_g_valid=.FALSE., tau_r_valid=.FALSE., tau_g_valid=.FALSE.)
      CALL qs_rho_update_rho(rho, qs_env)

      CALL timestop(handle)

   END SUBROUTINE qs_basis_reinit_energy

END MODULE qs_basis_gradient
