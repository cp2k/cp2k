!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate MAO's and analyze wavefunctions
!> \par History
!>      03.2016 created [JGH]
!>      12.2016 split into four modules [JGH]
!> \author JGH
! **************************************************************************************************
MODULE mao_basis
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE dbcsr_api,                       ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix_set,&
                                              dbcsr_distribution_type,&
                                              dbcsr_p_type,&
                                              dbcsr_reserve_diag_blocks,&
                                              dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE mao_methods,                     ONLY: mao_build_q
   USE mao_optimizer,                   ONLY: mao_optimize
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_basis'

   PUBLIC ::  mao_generate_basis

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mao_coef ...
!> \param ref_basis_set ...
!> \param pmat_external ...
!> \param smat_external ...
!> \param molecular ...
!> \param max_iter ...
!> \param eps_grad ...
!> \param nmao_external ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mao_generate_basis(qs_env, mao_coef, ref_basis_set, pmat_external, smat_external, &
                                 molecular, max_iter, eps_grad, nmao_external, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      CHARACTER(len=*), OPTIONAL                         :: ref_basis_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: pmat_external, smat_external
      LOGICAL, INTENT(IN), OPTIONAL                      :: molecular
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iter
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_grad
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: nmao_external
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_generate_basis', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=10)                                  :: mao_basis_set
      INTEGER                                            :: handle, iab, iatom, ikind, ispin, iw, &
                                                            mao_max_iter, natom, nbas, nimages, &
                                                            nkind, nmao, nspin
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: do_nmao_external, molecule
      REAL(KIND=dp)                                      :: electra(2), eps_filter, eps_fun, &
                                                            mao_eps_grad
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_q, matrix_smm, matrix_smo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: mao_basis_set_list, orb_basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: smm_list, smo_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      ! k-points?
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nimages = dft_control%nimages
      CPASSERT(nimages == 1)

      ! output
      IF (PRESENT(unit_nr)) THEN
         iw = unit_nr
      ELSE
         iw = -1
      END IF

      ! molecules
      IF (PRESENT(molecular)) THEN
         molecule = molecular
      ELSE
         molecule = .FALSE.
      END IF

      ! iterations
      IF (PRESENT(max_iter)) THEN
         mao_max_iter = max_iter
      ELSE
         mao_max_iter = 0
      END IF

      ! threshold
      IF (PRESENT(eps_grad)) THEN
         mao_eps_grad = eps_grad
      ELSE
         mao_eps_grad = 0.00001_dp
      END IF
      eps_fun = 10._dp*mao_eps_grad

      do_nmao_external = .FALSE.
      ! external number of MAOs per atom
      IF (PRESENT(nmao_external)) THEN
         do_nmao_external = .TRUE.
      END IF

      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
         WRITE (UNIT=iw, FMT="(T37,A)") "MAO BASIS"
         WRITE (iw, '(T2,A)') '!-----------------------------------------------------------------------------!'
      END IF

      ! Reference basis set
      IF (PRESENT(ref_basis_set)) THEN
         mao_basis_set = ref_basis_set
      ELSE
         mao_basis_set = "ORB"
      END IF
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (mao_basis_set_list(nkind), orb_basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         NULLIFY (mao_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (orb_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (basis_set_a, basis_set_b)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type="ORB")
         IF (ASSOCIATED(basis_set_a)) orb_basis_set_list(ikind)%gto_basis_set => basis_set_a
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_b, basis_type=mao_basis_set)
         IF (ASSOCIATED(basis_set_b)) mao_basis_set_list(ikind)%gto_basis_set => basis_set_b
      END DO
      IF (iw > 0) THEN
         DO ikind = 1, nkind
            IF (.NOT. ASSOCIATED(mao_basis_set_list(ikind)%gto_basis_set)) THEN
               WRITE (UNIT=iw, FMT="(T2,A,I4)") &
                  "WARNING: No MAO basis set associated with Kind ", ikind
            ELSE
               IF (do_nmao_external) THEN
                  nmao = nmao_external(ikind)
               ELSE
                  CALL get_qs_kind(qs_kind_set(ikind), mao=nmao)
               END IF
               nbas = mao_basis_set_list(ikind)%gto_basis_set%nsgf
               WRITE (UNIT=iw, FMT="(T2,A,I4,T30,A,I2,T56,A,I10)") &
                  "MAO basis set Kind ", ikind, " Number of MAO:", nmao, " Number of BSF:", nbas
            END IF
         END DO
      END IF

      ! neighbor lists
      NULLIFY (smm_list, smo_list)
      CALL setup_neighbor_list(smm_list, mao_basis_set_list, molecular=molecule, qs_env=qs_env)
      CALL setup_neighbor_list(smo_list, mao_basis_set_list, orb_basis_set_list, &
                               molecular=molecule, qs_env=qs_env)

      ! overlap matrices
      NULLIFY (matrix_smm, matrix_smo)
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL build_overlap_matrix_simple(ks_env, matrix_smm, &
                                       mao_basis_set_list, mao_basis_set_list, smm_list)
      CALL build_overlap_matrix_simple(ks_env, matrix_smo, &
                                       mao_basis_set_list, orb_basis_set_list, smo_list)

      ! get reference density matrix and overlap matrix
      IF (PRESENT(pmat_external)) THEN
         matrix_p => pmat_external
      ELSE
         CALL get_qs_env(qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      END IF
      IF (PRESENT(smat_external)) THEN
         matrix_s => smat_external
      ELSE
         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
      END IF

      nspin = SIZE(matrix_p, 1)
      eps_filter = 0.0_dp
      ! Q matrix
      CALL mao_build_q(matrix_q, matrix_p, matrix_s, matrix_smm, matrix_smo, smm_list, electra, eps_filter)

      ! MAO matrices
      CALL get_qs_env(qs_env=qs_env, natom=natom)
      CALL get_ks_env(ks_env=ks_env, particle_set=particle_set, dbcsr_dist=dbcsr_dist)
      NULLIFY (mao_coef)
      CALL dbcsr_allocate_matrix_set(mao_coef, nspin)
      ALLOCATE (row_blk_sizes(natom), col_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, &
                            basis=mao_basis_set_list)
      IF (do_nmao_external) THEN
         DO iatom = 1, natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 kind_number=ikind)
            col_blk_sizes(iatom) = nmao_external(ikind)
         END DO
      ELSE
         CALL get_particle_set(particle_set, qs_kind_set, nmao=col_blk_sizes)
      END IF

      ! check if MAOs have been specified
      DO iab = 1, natom
         IF (col_blk_sizes(iab) < 0) &
            CPABORT("Number of MAOs has to be specified in KIND section for all elements")
      END DO
      DO ispin = 1, nspin
         ! coefficients
         ALLOCATE (mao_coef(ispin)%matrix)
         CALL dbcsr_create(matrix=mao_coef(ispin)%matrix, &
                           name="MAO_COEF", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL dbcsr_reserve_diag_blocks(matrix=mao_coef(ispin)%matrix)
      END DO
      DEALLOCATE (row_blk_sizes, col_blk_sizes)

      ! optimize MAOs
      CALL mao_optimize(mao_coef, matrix_q, matrix_smm, electra, mao_max_iter, mao_eps_grad, iw)

      ! Deallocate the neighbor list structure
      IF (ASSOCIATED(smm_list)) THEN
         DO iab = 1, SIZE(smm_list)
            CALL deallocate_neighbor_list_set(smm_list(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (smm_list)
      END IF
      IF (ASSOCIATED(smo_list)) THEN
         DO iab = 1, SIZE(smo_list)
            CALL deallocate_neighbor_list_set(smo_list(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (smo_list)
      END IF

      DEALLOCATE (mao_basis_set_list, orb_basis_set_list)

      IF (ASSOCIATED(matrix_smm)) CALL dbcsr_deallocate_matrix_set(matrix_smm)
      IF (ASSOCIATED(matrix_smo)) CALL dbcsr_deallocate_matrix_set(matrix_smo)
      IF (ASSOCIATED(matrix_q)) CALL dbcsr_deallocate_matrix_set(matrix_q)

      CALL timestop(handle)

   END SUBROUTINE mao_generate_basis

END MODULE mao_basis
