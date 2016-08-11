!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods used by pao_main.F
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_methods
   USE ao_util,                         ONLY: exp_radius
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE bibliography,                    ONLY: Kolafa2004,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_add, cp_dbcsr_binary_read, cp_dbcsr_checksum, cp_dbcsr_col_block_sizes, &
        cp_dbcsr_complete_redistribute, cp_dbcsr_copy_into_existing, cp_dbcsr_create, &
        cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_distribution_new, &
        cp_dbcsr_get_block_p, cp_dbcsr_get_info, cp_dbcsr_init, cp_dbcsr_iterator, &
        cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
        cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, &
        cp_dbcsr_reserve_diag_blocks, cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_set, &
        cp_dbcsr_trace, cp_dbcsr_type, dbcsr_distribution_init, dbcsr_distribution_mp, &
        dbcsr_distribution_obj, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_obj
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dm_ls_scf_methods,               ONLY: density_matrix_trs4,&
                                              ls_scf_init_matrix_S
   USE dm_ls_scf_qs,                    ONLY: ls_scf_dm_to_ks,&
                                              ls_scf_qs_atomic_guess,&
                                              matrix_decluster,&
                                              matrix_qs_to_ls
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                              ls_scf_env_type
   USE iterate_matrix,                  ONLY: purify_mcweeny
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_walltime
   USE mathlib,                         ONLY: binomial,&
                                              diamat_all
   USE message_passing,                 ONLY: mp_max,&
                                              mp_sum
   USE pao_ml,                          ONLY: pao_ml_forces
   USE pao_param,                       ONLY: pao_calc_U,&
                                              pao_param_count,&
                                              pao_update_AB
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_initial_guess,                ONLY: calculate_atomic_fock_matrix
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              pao_descriptor_type,&
                                              pao_potential_type,&
                                              qs_kind_type,&
                                              set_qs_kind
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type

!$ USE OMP_LIB, ONLY: omp_get_level
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_methods'

   PUBLIC :: pao_print_atom_info, pao_init_kinds
   PUBLIC :: pao_build_orthogonalizer, pao_build_selector
   PUBLIC :: pao_build_diag_distribution
   PUBLIC :: pao_build_matrix_X, pao_build_core_hamiltonian
   PUBLIC :: pao_test_convergence
   PUBLIC :: pao_calc_energy, pao_check_trace_ps
   PUBLIC :: pao_store_P, pao_add_forces, pao_guess_initial_P
   PUBLIC :: pao_calc_outer_grad_lnv
   PUBLIC :: pao_check_grad

CONTAINS

! **************************************************************************************************
!> \brief Initialize qs kinds
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_init_kinds(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_init_kinds'

      INTEGER                                            :: handle, i, ikind, pao_basis_size
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(pao_descriptor_type), DIMENSION(:), POINTER   :: pao_descriptors
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)

      DO ikind = 1, SIZE(qs_kind_set)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          basis_set=basis_set, &
                          pao_basis_size=pao_basis_size, &
                          pao_potentials=pao_potentials, &
                          pao_descriptors=pao_descriptors)

         IF (pao_basis_size < 1) THEN
            ! pao disabled for ikind, set pao_basis_size to size of primary basis
            CALL set_qs_kind(qs_kind_set(ikind), pao_basis_size=basis_set%nsgf)
         ENDIF

         ! initialize radii of Gaussians to speedup screeing later on
         DO i = 1, SIZE(pao_potentials)
            pao_potentials(i)%beta_radius = exp_radius(0, pao_potentials(i)%beta, pao%eps_pgf, 1.0_dp)
         ENDDO
         DO i = 1, SIZE(pao_descriptors)
            pao_descriptors(i)%beta_radius = exp_radius(0, pao_descriptors(i)%beta, pao%eps_pgf, 1.0_dp)
            pao_descriptors(i)%screening_radius = exp_radius(0, pao_descriptors(i)%screening, pao%eps_pgf, 1.0_dp)
         ENDDO
      ENDDO
      CALL timestop(handle)
   END SUBROUTINE pao_init_kinds

! **************************************************************************************************
!> \brief Prints a one line summary for each atom.
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_print_atom_info(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_print_atom_info', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, natoms
      INTEGER, DIMENSION(:), POINTER                     :: pao_basis, param_cols, param_rows, &
                                                            pri_basis

      pri_basis => cp_dbcsr_row_block_sizes(pao%matrix_Y)
      pao_basis => cp_dbcsr_col_block_sizes(pao%matrix_Y)

      CPASSERT(SIZE(pao_basis) == SIZE(pri_basis))
      natoms = SIZE(pao_basis)

      param_rows => cp_dbcsr_row_block_sizes(pao%matrix_X)
      param_cols => cp_dbcsr_col_block_sizes(pao%matrix_X)
      CPASSERT(SIZE(param_rows) == natoms .AND. SIZE(param_cols) == natoms)

      IF (pao%iw_atoms > 0) THEN
         DO iatom = 1, natoms
            WRITE (pao%iw_atoms, "(A,I7,T20,A,I3,T45,A,I3,T65,A,I3)") &
               " PAO| atom: ", iatom, &
               " prim_basis: ", pri_basis(iatom), &
               " pao_basis: ", pao_basis(iatom), &
               " pao_params: ", (param_cols(iatom)*param_rows(iatom))
         ENDDO
      ENDIF
   END SUBROUTINE pao_print_atom_info

! **************************************************************************************************
!> \brief Constructs matrix_N and its inverse.
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_build_orthogonalizer(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_orthogonalizer'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, N
      LOGICAL                                            :: found
      REAL(dp)                                           :: v, w
      REAL(dp), DIMENSION(:), POINTER                    :: evals
      REAL(dp), DIMENSION(:, :), POINTER                 :: A, block_N, block_N_inv, block_S
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      CALL cp_dbcsr_init(pao%matrix_N)
      CALL cp_dbcsr_create(pao%matrix_N, template=matrix_s(1)%matrix, name="PAO matrix_N")
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_N)

      CALL cp_dbcsr_init(pao%matrix_N_inv)
      CALL cp_dbcsr_create(pao%matrix_N_inv, template=matrix_s(1)%matrix, name="PAO matrix_N_inv")
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_N_inv)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,matrix_s) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_N,block_N_inv,block_S,found,N,A,evals,k,i,j,w,v)
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_N)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_N)
         iatom = arow; CPASSERT(arow == acol)

         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_N_inv, row=iatom, col=iatom, block=block_N_inv, found=found)
         CPASSERT(ASSOCIATED(block_N_inv))

         CALL cp_dbcsr_get_block_p(matrix=matrix_s(1)%matrix, row=iatom, col=iatom, block=block_S, found=found)
         CPASSERT(ASSOCIATED(block_S))

         N = SIZE(block_S, 1); CPASSERT(SIZE(block_S, 1) == SIZE(block_S, 2)) ! primary basis size
         ALLOCATE (A(N, N), evals(N))

         ! take square root of atomic overlap matrix
         A = block_S
         CALL diamat_all(A, evals) !afterwards A contains the eigenvectors
         block_N = 0.0_dp
         block_N_inv = 0.0_dp
         DO k = 1, N
            ! NOTE: To maintain a consistent notation with the Berghold paper,
            ! the "_inv" is swapped: N^{-1}=sqrt(S); N=sqrt(S)^{-1}
            w = 1.0_dp/SQRT(evals(k))
            v = SQRT(evals(k))
            DO i = 1, N
               DO j = 1, N
                  block_N(i, j) = block_N(i, j)+w*A(i, k)*A(j, k)
                  block_N_inv(i, j) = block_N_inv(i, j)+v*A(i, k)*A(j, k)
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE (A, evals)
      END DO
      CALL cp_dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      ! store a copy that is distributed according to pao%diag_distribution
      CALL cp_dbcsr_init(pao%matrix_N_diag)
      CALL cp_dbcsr_create(pao%matrix_N_diag, &
                           name="PAO matrix_N_diag", &
                           dist=pao%diag_distribution, &
                           template=matrix_s(1)%matrix)
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_N_diag)
      CALL cp_dbcsr_complete_redistribute(pao%matrix_N, pao%matrix_N_diag)

      CALL timestop(handle)
   END SUBROUTINE pao_build_orthogonalizer

! **************************************************************************************************
!> \brief Build rectangular matrix to converert between primary and PAO basis.
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_build_selector(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_selector'

      INTEGER                                            :: acol, arow, handle, i, iatom, ikind, M, &
                                                            natoms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_aux, blk_sizes_pri
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_Y
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      natom=natoms, &
                      matrix_s=matrix_s, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      blk_sizes_pri => cp_dbcsr_col_block_sizes(matrix_s(1)%matrix)
      ALLOCATE (blk_sizes_aux(natoms))
      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), pao_basis_size=M)
         CPASSERT(M > 0)
         IF (blk_sizes_pri(iatom) < M) &
            CPABORT("PAO basis size exceeds primary basis size.")
         blk_sizes_aux(iatom) = M
      ENDDO

      CALL cp_dbcsr_init(pao%matrix_Y)
      CALL cp_dbcsr_create(pao%matrix_Y, &
                           template=matrix_s(1)%matrix, &
                           matrix_type="N", &
                           row_blk_size=blk_sizes_pri, &
                           col_blk_size=blk_sizes_aux, &
                           name="PAO matrix_Y")
      DEALLOCATE (blk_sizes_aux)

      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_Y)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao) &
!$OMP PRIVATE(iter,arow,acol,block_Y,i,M)
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_Y)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_Y)
         M = SIZE(block_Y, 2) ! size of pao basis
         block_Y = 0.0_dp
         DO i = 1, M
            block_Y(i, i) = 1.0_dp
         ENDDO
      END DO
      CALL cp_dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE pao_build_selector

! **************************************************************************************************
!> \brief Creates new DBCSR distribution which spreads diagonal blocks evenly across ranks
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_build_diag_distribution(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_diag_distribution'

      INTEGER                                            :: handle, iatom, natoms, pgrid_cols, &
                                                            pgrid_rows
      INTEGER, DIMENSION(:), POINTER                     :: diag_col_dist, diag_row_dist, no_clusters
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dbcsr_distribution_obj)                       :: main_dist
      TYPE(dbcsr_mp_obj)                                 :: mp_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, natom=natoms, matrix_s=matrix_s)

      ! get processor grid from matrix_s
      CALL cp_dbcsr_get_info(matrix=matrix_s(1)%matrix, distribution=main_dist)
      mp_env = dbcsr_distribution_mp(main_dist)
      pgrid_rows = dbcsr_mp_nprows(mp_env)
      pgrid_cols = dbcsr_mp_npcols(mp_env)

      ! create new mapping of matrix-grid to processor-grid
      ALLOCATE (diag_row_dist(natoms), diag_col_dist(natoms))
      DO iatom = 1, natoms
         diag_row_dist(iatom) = MOD(iatom-1, pgrid_rows)
         diag_col_dist(iatom) = MOD((iatom-1)/pgrid_rows, pgrid_cols)
      ENDDO
      NULLIFY (no_clusters)

      ! instanciate distribution object
      CALL dbcsr_distribution_init(pao%diag_distribution)
      CALL cp_dbcsr_distribution_new(pao%diag_distribution, &
                                     mp_env=mp_env, &
                                     row_dist_block=diag_row_dist, &
                                     col_dist_block=diag_col_dist, &
                                     row_dist_cluster=no_clusters, &
                                     col_dist_cluster=no_clusters)

      DEALLOCATE (diag_row_dist, diag_col_dist)

      CALL timestop(handle)
   END SUBROUTINE pao_build_diag_distribution

! **************************************************************************************************
!> \brief Creates the matrix_X
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_build_matrix_X(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_build_matrix_X'

      INTEGER                                            :: handle, iatom, ikind, natoms
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      natom=natoms, &
                      particle_set=particle_set)

      ! determine block-sizes of matrix_X
      ALLOCATE (row_blk_size(natoms), col_blk_size(natoms))
      col_blk_size = 1
      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL pao_param_count(pao, qs_env, ikind, nparams=row_blk_size(iatom))
      ENDDO

      ! build actual matrix_X
      CALL cp_dbcsr_init(pao%matrix_X)
      CALL cp_dbcsr_create(pao%matrix_X, &
                           name="PAO matrix_X", &
                           dist=pao%diag_distribution, &
                           matrix_type="N", &
                           row_blk_size=row_blk_size, &
                           col_blk_size=col_blk_size)
      DEALLOCATE (row_blk_size, col_blk_size)

      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_X)
      CALL cp_dbcsr_set(pao%matrix_X, 0.0_dp)

      CALL timestop(handle)
   END SUBROUTINE pao_build_matrix_X

! **************************************************************************************************
!> \brief Creates the matrix_H0 which contains the core hamiltonian
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_build_core_hamiltonian(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, &
                      matrix_s=matrix_s, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set)

      ! allocate matrix_H0
      CALL cp_dbcsr_init(pao%matrix_H0)
      CALL cp_dbcsr_create(pao%matrix_H0, &
                           name="PAO matrix_H0", &
                           dist=pao%diag_distribution, &
                           template=matrix_s(1)%matrix)
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_H0)

      ! calculate inital atomic fock matrix H0
      ! Can't use matrix_ks from ls_scf_qs_atomic_guess(), because it's not rotationally invariant.
      ! getting H0 directly from the atomic code
      CALL calculate_atomic_fock_matrix(pao%matrix_H0, &
                                        particle_set, &
                                        atomic_kind_set, &
                                        qs_kind_set, &
                                        output_unit=pao%iw)

   END SUBROUTINE pao_build_core_hamiltonian

! **************************************************************************************************
!> \brief Test whether the PAO optimization has reached convergence
!> \param pao ...
!> \param ls_scf_env ...
!> \param new_energy ...
!> \param is_converged ...
! **************************************************************************************************
   SUBROUTINE pao_test_convergence(pao, ls_scf_env, new_energy, is_converged)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(ls_scf_env_type)                              :: ls_scf_env
      REAL(KIND=dp), INTENT(IN)                          :: new_energy
      LOGICAL, INTENT(OUT)                               :: is_converged

      REAL(KIND=dp)                                      :: energy_diff, loop_eps, now, time_diff

      ! calculate progress
      energy_diff = new_energy-pao%energy_prev
      pao%energy_prev = new_energy
      now = m_walltime()
      time_diff = now-pao%step_start_time
      pao%step_start_time = now

      ! convergence criterion
      loop_eps = pao%norm_G/ls_scf_env%nelectron_total
      is_converged = loop_eps < pao%eps_pao

      IF (pao%istep > 1) THEN
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| energy improvement:", energy_diff
         ! IF(energy_diff>0.0_dp) CPWARN("PAO| energy increased")

         ! print one-liner
         IF (pao%iw > 0) WRITE (pao%iw, '(A,I6,11X,F20.9,1X,E10.3,1X,E10.3,1X,F9.3)') &
            " PAO| step ", &
            pao%istep, &
            new_energy, &
            loop_eps, &
            pao%linesearch%step_size, & !prev step, which let to the current energy
            time_diff
      ENDIF
   END SUBROUTINE pao_test_convergence

! **************************************************************************************************
!> \brief Calculate the pao energy
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE pao_calc_energy(pao, qs_env, ls_scf_env, energy)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      REAL(KIND=dp), INTENT(OUT)                         :: energy

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: penalty, trace_PH

      CALL timeset(routineN, handle)

      ! calculate matrix U, which determines the pao basis
      CALL pao_update_AB(pao, qs_env, ls_scf_env%ls_mstruct, penalty=penalty)

      ! calculat S, S_inv, S_sqrt, and S_sqrt_inv in the new pao basis
      CALL pao_rebuild_S(qs_env, ls_scf_env)

      ! calulate the density matrix P in the pao basis
      CALL pao_dm_trs4(qs_env, ls_scf_env)

      ! calculate the energy from the trace(PH) in the pao basis
      energy = 0.0_dp
      DO ispin = 1, ls_scf_env%nspins
         CALL cp_dbcsr_trace(ls_scf_env%matrix_p(ispin), ls_scf_env%matrix_ks(ispin), trace_PH)
         energy = energy+trace_PH
      ENDDO

      ! add penalty term
      energy = energy+penalty

      IF (pao%iw > 0) THEN
         WRITE (pao%iw, *) ""
         WRITE (pao%iw, *) "PAO| energy:", energy, "penalty:", penalty
      ENDIF
      CALL timestop(handle)
   END SUBROUTINE pao_calc_energy

! **************************************************************************************************
!> \brief Ensure that the number of electrons is correct.
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_check_trace_PS(ls_scf_env)
      TYPE(ls_scf_env_type)                              :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_check_trace_PS', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: tmp, trace_PS
      TYPE(cp_dbcsr_type)                                :: matrix_S_desym

      CALL timeset(routineN, handle)
      CALL cp_dbcsr_init(matrix_S_desym)
      CALL cp_dbcsr_create(matrix_S_desym, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL cp_dbcsr_desymmetrize(ls_scf_env%matrix_s, matrix_S_desym)

      trace_PS = 0.0_dp
      DO ispin = 1, ls_scf_env%nspins
         CALL cp_dbcsr_trace(ls_scf_env%matrix_p(ispin), matrix_S_desym, tmp)
         trace_PS = trace_PS+tmp
      ENDDO

      CALL cp_dbcsr_release(matrix_S_desym)

      IF (ABS(ls_scf_env%nelectron_total-trace_PS) > 0.5) &
         CPABORT("Number of electrons wrong. Trace(PS) ="//cp_to_string(trace_PS))

      CALL timestop(handle)
   END SUBROUTINE pao_check_trace_PS

! **************************************************************************************************
!> \brief Read primary density matrix from file.
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_read_preopt_dm(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_read_preopt_dm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: cs_pos
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s, rho_ao
      TYPE(cp_dbcsr_type)                                :: matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      rho=rho, &
                      energy=energy)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      IF (dft_control%nspins /= 1) CPABORT("open shell not yet implemented")

      DO ispin = 1, dft_control%nspins
         CALL cp_dbcsr_init(matrix_tmp) ! read into tmp-matrix to preserve sparsity of rho_ao
         CALL cp_dbcsr_binary_read(pao%preopt_dm_file, matrix_new=matrix_tmp, &
                                   distribution=cp_dbcsr_distribution(matrix_s(1)%matrix))
         cs_pos = cp_dbcsr_checksum(matrix_tmp, pos=.TRUE.)
         IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Read restart DM "// &
            TRIM(pao%preopt_dm_file)//" with checksum: ", cs_pos
         CALL cp_dbcsr_copy_into_existing(rho_ao(ispin)%matrix, matrix_tmp)
         CALL cp_dbcsr_release(matrix_tmp)
      ENDDO

      ! calculate corresponding ks matrix
      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
      CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., &
                               just_energy=.FALSE., print_active=.TRUE.)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Quickstep energy from restart density:", energy%total

      CALL timestop(handle)

   END SUBROUTINE pao_read_preopt_dm

! **************************************************************************************************
!> \brief Rebuilds S, S_inv, S_sqrt, and S_sqrt_inv in the pao basis
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_rebuild_S(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_rebuild_S', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s

      CALL timeset(routineN, handle)

      CALL cp_dbcsr_release(ls_scf_env%matrix_s_inv)
      CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt)
      CALL cp_dbcsr_release(ls_scf_env%matrix_s_sqrt_inv)

      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL ls_scf_init_matrix_s(matrix_s(1)%matrix, ls_scf_env)

      CALL timestop(handle)
   END SUBROUTINE pao_rebuild_S

! **************************************************************************************************
!> \brief Calculate density matrix using TRS4 purification
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_dm_trs4(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_dm_trs4', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: project_name
      INTEGER                                            :: handle, ispin, nelectron_spin_real, nspin
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: homo_spin, lumo_spin, mu_spin
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      project_name = logger%iter_info%project_name
      nspin = ls_scf_env%nspins

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
      DO ispin = 1, nspin
         CALL matrix_qs_to_ls(ls_scf_env%matrix_ks(ispin), matrix_ks(ispin)%matrix, &
                              ls_scf_env%ls_mstruct, covariant=.TRUE.)

         nelectron_spin_real = ls_scf_env%nelectron_spin(ispin)
         IF (ls_scf_env%nspins == 1) nelectron_spin_real = nelectron_spin_real/2
         CALL density_matrix_trs4(ls_scf_env%matrix_p(ispin), ls_scf_env%matrix_ks(ispin), &
                                  ls_scf_env%matrix_s_sqrt_inv, &
                                  nelectron_spin_real, ls_scf_env%eps_filter, homo_spin, lumo_spin, mu_spin, &
                                  dynamic_threshold=.FALSE., converged=converged, &
                                  max_iter_lanczos=ls_scf_env%max_iter_lanczos, &
                                  eps_lanczos=ls_scf_env%eps_lanczos)
         IF (.NOT. converged) CPABORT("TRS4 did not converge")
      ENDDO

      IF (nspin == 1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(1), 2.0_dp)

      CALL timestop(handle)
   END SUBROUTINE pao_dm_trs4

! **************************************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dU
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M_diag the derivate, matrix uses pao%diag_distribution
! **************************************************************************************************
   SUBROUTINE pao_calc_outer_grad_lnv(qs_env, ls_scf_env, matrix_M_diag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      TYPE(cp_dbcsr_type)                                :: matrix_M_diag

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_outer_grad_lnv', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nspin
      INTEGER, DIMENSION(:), POINTER                     :: pao_blk_sizes
      REAL(KIND=dp)                                      :: filter_eps
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks, matrix_s, rho_ao
      TYPE(cp_dbcsr_type) :: matrix_HB, matrix_HPS, matrix_M, matrix_M1, matrix_M1_dc, matrix_M2, &
         matrix_M2_dc, matrix_M3, matrix_M3_dc, matrix_NHB, matrix_NHBM2, matrix_NPA, &
         matrix_NPAM1, matrix_NSB, matrix_NSBM3, matrix_PA, matrix_PH, matrix_PHP, matrix_PSP, &
         matrix_SB, matrix_SP
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      ls_mstruct => ls_scf_env%ls_mstruct
      pao => ls_scf_env%pao_env

      CALL get_qs_env(qs_env, &
                      rho=rho, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      nspin = dft_control%nspins
      filter_eps = ls_scf_env%eps_filter
      pao_blk_sizes => cp_dbcsr_col_block_sizes(ls_mstruct%matrix_A)

      IF (nspin /= 1) CPABORT("open shell not yet implemented")
      !TODO: handle openshell case properly

      ! notation according to pao_math_lnv.pdf

      ! calculation uses distribution of matrix_s, after we redistribute using pao%diag_distribution
      CALL cp_dbcsr_init(matrix_M)
      CALL cp_dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")
      CALL cp_dbcsr_reserve_diag_blocks(matrix_M)

      !---------------------------------------------------------------------------
      ! calculate need products in pao basis
      CALL cp_dbcsr_init(matrix_PH)
      CALL cp_dbcsr_create(matrix_PH, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), ls_scf_env%matrix_ks(1), &
                             0.0_dp, matrix_PH, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_PHP)
      CALL cp_dbcsr_create(matrix_PHP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_PH, ls_scf_env%matrix_p(1), &
                             0.0_dp, matrix_PHP, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_SP)
      CALL cp_dbcsr_create(matrix_SP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s, ls_scf_env%matrix_p(1), &
                             0.0_dp, matrix_SP, filter_eps=filter_eps)

      IF (SIZE(ls_scf_env%matrix_p) == 1) CALL cp_dbcsr_scale(matrix_SP, 0.5_dp)

      CALL cp_dbcsr_init(matrix_HPS)
      CALL cp_dbcsr_create(matrix_HPS, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, ls_scf_env%matrix_ks(1), matrix_SP, &
                             0.0_dp, matrix_HPS, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_PSP)
      CALL cp_dbcsr_create(matrix_PSP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), matrix_SP, &
                             0.0_dp, matrix_PSP, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! M1 = dE_lnv / dP_pao
      CALL cp_dbcsr_init(matrix_M1)
      CALL cp_dbcsr_create(matrix_M1, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL cp_dbcsr_multiply("N", "T", 3.0_dp, ls_scf_env%matrix_ks(1), matrix_SP, &
                             1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "N", 3.0_dp, matrix_SP, ls_scf_env%matrix_ks(1), &
                             1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_HPS, matrix_SP, &
                             1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_SP, matrix_HPS, &
                             1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_SP, matrix_HPS, &
                             1.0_dp, matrix_M1, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL cp_dbcsr_init(matrix_M1_dc)
      CALL cp_dbcsr_create(matrix_M1_dc, &
                           template=matrix_s(1)%matrix, &
                           row_blk_size=pao_blk_sizes, &
                           col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M1_dc, matrix_M1, ls_mstruct)

      !---------------------------------------------------------------------------
      ! M2 = dE_lnv / dH
      CALL cp_dbcsr_init(matrix_M2)
      CALL cp_dbcsr_create(matrix_M2, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL cp_dbcsr_add(matrix_M2, matrix_PSP, 1.0_dp, 3.0_dp)

      CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PSP, matrix_SP, &
                             1.0_dp, matrix_M2, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL cp_dbcsr_init(matrix_M2_dc)
      CALL cp_dbcsr_create(matrix_M2_dc, &
                           template=matrix_s(1)%matrix, &
                           row_blk_size=pao_blk_sizes, &
                           col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M2_dc, matrix_M2, ls_mstruct)

      !---------------------------------------------------------------------------
      ! M3 = dE_lnv / dS
      CALL cp_dbcsr_init(matrix_M3)
      CALL cp_dbcsr_create(matrix_M3, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL cp_dbcsr_add(matrix_M3, matrix_PHP, 1.0_dp, 3.0_dp)

      CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PHP, matrix_SP, &
                             1.0_dp, matrix_M3, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "T", -2.0_dp, matrix_PSP, matrix_PH, &
                             1.0_dp, matrix_M3, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL cp_dbcsr_init(matrix_M3_dc)
      CALL cp_dbcsr_create(matrix_M3_dc, &
                           template=matrix_s(1)%matrix, &
                           row_blk_size=pao_blk_sizes, &
                           col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M3_dc, matrix_M3, ls_mstruct)

      !---------------------------------------------------------------------------
      ! combine M1 with matrices from primary basis
      CALL cp_dbcsr_init(matrix_PA)
      CALL cp_dbcsr_create(matrix_PA, template=ls_mstruct%matrix_A, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, rho_ao(1)%matrix, ls_mstruct%matrix_A, &
                             0.0_dp, matrix_PA, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_NPA)
      CALL cp_dbcsr_create(matrix_NPA, template=ls_mstruct%matrix_A, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N_inv, matrix_PA, &
                             0.0_dp, matrix_NPA, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_NPAM1)
      CALL cp_dbcsr_create(matrix_NPAM1, template=ls_mstruct%matrix_A, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NPA, matrix_M1_dc, &
                             0.0_dp, matrix_NPAM1, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NPAM1, pao%matrix_Y, &
                             1.0_dp, matrix_M, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! combine M2 with matrices from primary basis
      CALL cp_dbcsr_init(matrix_HB)
      CALL cp_dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B, &
                             0.0_dp, matrix_HB, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_NHB)
      CALL cp_dbcsr_create(matrix_NHB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_HB, &
                             0.0_dp, matrix_NHB, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_NHBM2)
      CALL cp_dbcsr_create(matrix_NHBM2, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NHB, matrix_M2_dc, &
                             0.0_dp, matrix_NHBM2, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NHBM2, pao%matrix_Y, &
                             1.0_dp, matrix_M, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! combine M3 with matrices from primary basis
      CALL cp_dbcsr_init(matrix_SB)
      CALL cp_dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B, &
                             0.0_dp, matrix_SB, filter_eps=filter_eps)

      IF (SIZE(ls_scf_env%matrix_p) == 1) CALL cp_dbcsr_scale(matrix_SB, 0.5_dp)

      CALL cp_dbcsr_init(matrix_NSB)
      CALL cp_dbcsr_create(matrix_NSB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_SB, &
                             0.0_dp, matrix_NSB, filter_eps=filter_eps)

      CALL cp_dbcsr_init(matrix_NSBM3)
      CALL cp_dbcsr_create(matrix_NSBM3, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_NSB, matrix_M3_dc, &
                             0.0_dp, matrix_NSBM3, filter_eps=filter_eps)

      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_NSBM3, pao%matrix_Y, &
                             1.0_dp, matrix_M, filter_eps=filter_eps)

      IF (SIZE(ls_scf_env%matrix_p) == 1) CALL cp_dbcsr_scale(matrix_M, 2.0_dp)

      !---------------------------------------------------------------------------
      ! redistribute using pao%diag_distribution
      CALL cp_dbcsr_init(matrix_M_diag)
      CALL cp_dbcsr_create(matrix_M_diag, &
                           name="PAO matrix_M", &
                           matrix_type="N", &
                           dist=pao%diag_distribution, &
                           template=matrix_s(1)%matrix)
      CALL cp_dbcsr_reserve_diag_blocks(matrix_M_diag)
      CALL cp_dbcsr_complete_redistribute(matrix_M, matrix_M_diag)

      !---------------------------------------------------------------------------
      ! cleanup: TODO release matrices as early as possible
      CALL cp_dbcsr_release(matrix_PH)
      CALL cp_dbcsr_release(matrix_PHP)
      CALL cp_dbcsr_release(matrix_SP)
      CALL cp_dbcsr_release(matrix_HPS)
      CALL cp_dbcsr_release(matrix_PSP)
      CALL cp_dbcsr_release(matrix_M)
      CALL cp_dbcsr_release(matrix_M1)
      CALL cp_dbcsr_release(matrix_M2)
      CALL cp_dbcsr_release(matrix_M3)
      CALL cp_dbcsr_release(matrix_M1_dc)
      CALL cp_dbcsr_release(matrix_M2_dc)
      CALL cp_dbcsr_release(matrix_M3_dc)
      CALL cp_dbcsr_release(matrix_PA)
      CALL cp_dbcsr_release(matrix_NPA)
      CALL cp_dbcsr_release(matrix_NPAM1)
      CALL cp_dbcsr_release(matrix_HB)
      CALL cp_dbcsr_release(matrix_NHB)
      CALL cp_dbcsr_release(matrix_NHBM2)
      CALL cp_dbcsr_release(matrix_SB)
      CALL cp_dbcsr_release(matrix_NSB)
      CALL cp_dbcsr_release(matrix_NSBM3)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_outer_grad_lnv

! **************************************************************************************************
!> \brief Debugging routine for checking the analytic gradient.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_check_grad(pao, qs_env, ls_scf_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_check_grad', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iatom, j, natoms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_col, blk_sizes_row
      LOGICAL                                            :: found
      REAL(dp)                                           :: delta, delta_max, eps, Gij_num
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G, block_X
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct

      IF (pao%check_grad_tol < 0.0_dp) RETURN ! no checking

      CALL timeset(routineN, handle)

      ls_mstruct => ls_scf_env%ls_mstruct

      CALL get_qs_env(qs_env, para_env=para_env, natom=natoms)

      eps = pao%num_grad_eps
      delta_max = 0.0_dp

      blk_sizes_col => cp_dbcsr_col_block_sizes(pao%matrix_X)
      blk_sizes_row => cp_dbcsr_row_block_sizes(pao%matrix_X)

      ! can not use an iterator here, because other DBCSR routines are called within loop.
      DO iatom = 1, natoms
         IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO| checking gradient of atom ', iatom
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=iatom, col=iatom, block=block_X, found=found)

         IF (ASSOCIATED(block_X)) THEN !only one node actually has the block
            CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
            CPASSERT(ASSOCIATED(block_G))
         ENDIF

         DO i = 1, blk_sizes_row(iatom)
            DO j = 1, blk_sizes_col(iatom)
               SELECT CASE (pao%num_grad_order)
               CASE (2) ! calculate derivative to 2th order
                  Gij_num = -eval_point(block_X, i, j, -eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num+eval_point(block_X, i, j, +eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num/(2.0_dp*eps)

               CASE (4) ! calculate derivative to 4th order
                  Gij_num = eval_point(block_X, i, j, -2_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num-8_dp*eval_point(block_X, i, j, -1_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num+8_dp*eval_point(block_X, i, j, +1_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num-eval_point(block_X, i, j, +2_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num/(12.0_dp*eps)

               CASE (6) ! calculate derivative to 6th order
                  Gij_num = -1_dp*eval_point(block_X, i, j, -3_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num+9_dp*eval_point(block_X, i, j, -2_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num-45_dp*eval_point(block_X, i, j, -1_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num+45_dp*eval_point(block_X, i, j, +1_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num-9_dp*eval_point(block_X, i, j, +2_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num+1_dp*eval_point(block_X, i, j, +3_dp*eps, pao, ls_scf_env, qs_env)
                  Gij_num = Gij_num/(60.0_dp*eps)

               CASE DEFAULT
                  CPABORT("Unsupported numerical derivative order: "//cp_to_string(pao%num_grad_order))
               END SELECT

               IF (ASSOCIATED(block_X)) THEN
                  delta = ABS(Gij_num-block_G(i, j))
                  delta_max = MAX(delta_max, delta)
                  !WRITE (*,*) "gradient check", iatom, i, j, Gij_num, block_G(i,j), delta
               ENDIF
            ENDDO
         ENDDO
      END DO

      CALL mp_max(delta_max, para_env%group)
      IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO| checked gradient, max delta:', delta_max
      IF (delta_max > pao%check_grad_tol) CALL cp_abort(__LOCATION__, &
                                                        "Analytic and numeric gradients differ too much:"//cp_to_string(delta_max))

      CALL timestop(handle)
   END SUBROUTINE pao_check_grad

! **************************************************************************************************
!> \brief Helper routine for pao_check_grad()
!> \param block_X ...
!> \param i ...
!> \param j ...
!> \param eps ...
!> \param pao ...
!> \param ls_scf_env ...
!> \param qs_env ...
!> \retval energy ...
! **************************************************************************************************
   FUNCTION eval_point(block_X, i, j, eps, pao, ls_scf_env, qs_env) RESULT(energy)
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X
      INTEGER, INTENT(IN)                                :: i, j
      REAL(dp), INTENT(IN)                               :: eps
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp)                                           :: energy

      CHARACTER(len=*), PARAMETER :: routineN = 'eval_point', routineP = moduleN//':'//routineN

      REAL(dp)                                           :: old_Xij

      IF (ASSOCIATED(block_X)) THEN
         old_Xij = block_X(i, j) ! backup old block_X
         block_X(i, j) = block_X(i, j)+eps ! add pertubation
      ENDIF

      ! calculate energy
      CALL pao_calc_energy(pao, qs_env, ls_scf_env, energy)

      ! restore old block_X
      IF (ASSOCIATED(block_X)) THEN
         block_X(i, j) = old_Xij
      ENDIF

   END FUNCTION eval_point

! **************************************************************************************************
!> \brief Stores density matrix as initial guess for next SCF optimization.
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_store_P(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_store_P'

      INTEGER                                            :: handle, ispin, istore
      INTEGER, DIMENSION(:), POINTER                     :: pao_blk_sizes
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_dbcsr_type)                                :: matrix_p_dc, matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao

      IF (ls_scf_env%scf_history%nstore == 0) RETURN
      CALL timeset(routineN, handle)
      ls_mstruct => ls_scf_env%ls_mstruct
      pao => ls_scf_env%pao_env
      CALL get_qs_env(qs_env, dft_control=dft_control, matrix_s=matrix_s)

      ls_scf_env%scf_history%istore = ls_scf_env%scf_history%istore+1
      istore = MOD(ls_scf_env%scf_history%istore-1, ls_scf_env%scf_history%nstore)+1
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Storing density matrix for ASPC guess in slot:", istore

      ! initialize storage
      IF (ls_scf_env%scf_history%istore <= ls_scf_env%scf_history%nstore) THEN
         DO ispin = 1, dft_control%nspins
            CALL cp_dbcsr_init(ls_scf_env%scf_history%matrix(ispin, istore))
            CALL cp_dbcsr_create(ls_scf_env%scf_history%matrix(ispin, istore), template=matrix_s(1)%matrix)
         ENDDO
      ENDIF

      CALL cp_dbcsr_init(matrix_p_dc)
      pao_blk_sizes => cp_dbcsr_col_block_sizes(ls_mstruct%matrix_A)
      CALL cp_dbcsr_create(matrix_p_dc, &
                           template=matrix_s(1)%matrix, &
                           row_blk_size=pao_blk_sizes, &
                           col_blk_size=pao_blk_sizes)

      CALL cp_dbcsr_init(matrix_tmp)
      CALL cp_dbcsr_create(matrix_tmp, template=ls_mstruct%matrix_B)

      ! We are storing the density matrix in the non-orthonormal primary basis.
      ! While the orthonormal basis would yield better extrapolations,
      ! we simply can not affort to calculat S_sqrt in the primary basis.
      ! We can't call matrix_ls_to_qs because it uses cp_dbcsr_copy_into_existing()
      DO ispin = 1, dft_control%nspins
         ! reverse possible molecular clustering
         CALL matrix_decluster(matrix_p_dc, ls_scf_env%matrix_p(ispin), ls_mstruct)

         ! transform into primary basis
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_mstruct%matrix_B, matrix_p_dc, &
                                0.0_dp, matrix_tmp)
         CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, ls_mstruct%matrix_B, &
                                0.0_dp, ls_scf_env%scf_history%matrix(ispin, istore))
      ENDDO
      CALL cp_dbcsr_release(matrix_tmp)
      CALL cp_dbcsr_release(matrix_p_dc)

      CALL timestop(handle)
   END SUBROUTINE pao_store_P

! **************************************************************************************************
!> \brief Provide an initial guess for the density matrix
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_guess_initial_P(pao, qs_env, ls_scf_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      IF (ls_scf_env%scf_history%istore > 0) THEN
         CALL pao_aspc_guess_P(pao, qs_env, ls_scf_env)
         pao%need_initial_scf = .TRUE.
      ELSE
         IF (LEN_TRIM(pao%preopt_dm_file) > 0) THEN
            CALL pao_read_preopt_dm(pao, qs_env)
            pao%need_initial_scf = .FALSE.
            pao%preopt_dm_file = "" ! load only for first MD step
         ELSE
            CALL ls_scf_qs_atomic_guess(qs_env, ls_scf_env%energy_init)
            IF (pao%iw > 0) WRITE (pao%iw, '(A,F20.9)') " PAO| Energy from initial atomic guess:", ls_scf_env%energy_init
            pao%need_initial_scf = .TRUE.
         ENDIF
      ENDIF

   END SUBROUTINE pao_guess_initial_P

! **************************************************************************************************
!> \brief Run the Always Stable Predictor-Corrector to guess an initial density matrix
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_aspc_guess_P(pao, qs_env, ls_scf_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_aspc_guess_P'

      INTEGER                                            :: handle, iaspc, ispin, istore, naspc
      REAL(dp)                                           :: alpha
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_dbcsr_type)                                :: matrix_P, matrix_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct

      CALL timeset(routineN, handle)
      ls_mstruct => ls_scf_env%ls_mstruct
      CPASSERT(ls_scf_env%scf_history%istore > 0)
      CALL cite_reference(Kolafa2004)
      CALL get_qs_env(qs_env, dft_control=dft_control, matrix_s=matrix_s)

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Calculating initial guess with ASPC"

      CALL cp_dbcsr_init(matrix_tmp)
      CALL cp_dbcsr_create(matrix_tmp, template=ls_mstruct%matrix_A)
      CALL cp_dbcsr_init(matrix_P)
      CALL cp_dbcsr_create(matrix_P, template=matrix_s(1)%matrix)

      naspc = MIN(ls_scf_env%scf_history%istore, ls_scf_env%scf_history%nstore)
      DO ispin = 1, dft_control%nspins
         ! actual extrapolation
         CALL cp_dbcsr_set(matrix_P, 0.0_dp)
         DO iaspc = 1, naspc
            alpha = (-1.0_dp)**(iaspc+1)*REAL(iaspc, KIND=dp)* &
                    binomial(2*naspc, naspc-iaspc)/binomial(2*naspc-2, naspc-1)
            istore = MOD(ls_scf_env%scf_history%istore-iaspc, ls_scf_env%scf_history%nstore)+1
            CALL cp_dbcsr_add(matrix_P, ls_scf_env%scf_history%matrix(ispin, istore), 1.0_dp, alpha)
         ENDDO

         ! transform back from primary basis into pao basis
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_P, ls_mstruct%matrix_A, 0.0_dp, matrix_tmp)
         CALL cp_dbcsr_multiply("T", "N", 1.0_dp, ls_mstruct%matrix_A, matrix_tmp, 0.0_dp, ls_scf_env%matrix_p(ispin))
      ENDDO

      CALL cp_dbcsr_release(matrix_tmp)
      CALL cp_dbcsr_release(matrix_P)

      ! linear combination of P's is not idempotent. A bit of McWeeny is needed to ensure it is again
      DO ispin = 1, dft_control%nspins
         IF (dft_control%nspins == 1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin), 0.5_dp)
         ! we could go to the orthonomal basis, but it seems not worth the trouble
         ! TODO : 10 iterations is a conservative upper bound, figure out when it fails
         CALL purify_mcweeny(ls_scf_env%matrix_p(ispin:ispin), ls_scf_env%matrix_s, ls_scf_env%eps_filter, 10)
         IF (dft_control%nspins == 1) CALL cp_dbcsr_scale(ls_scf_env%matrix_p(ispin), 2.0_dp)
      ENDDO

      ! compute corresponding energy and ks matrix
      CALL ls_scf_dm_to_ks(qs_env, ls_scf_env, ls_scf_env%energy_init, iscf=0)

      CALL timestop(handle)
   END SUBROUTINE pao_aspc_guess_P

! **************************************************************************************************
!> \brief Calculate the forces contributed by PAO
!> \param qs_env ...
!> \param ls_scf_env ...
! **************************************************************************************************
   SUBROUTINE pao_add_forces(qs_env, ls_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_add_forces'

      INTEGER                                            :: handle, iatom, natoms
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: forces
      TYPE(cp_dbcsr_type)                                :: matrix_M
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      pao => ls_scf_env%pao_env

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Adding forces."

      IF (pao%max_pao /= 0) THEN
         IF (pao%penalty_strength /= 0.0_dp) &
            CPABORT("PAO forces require PENALTY_STRENGTH or MAX_PAO set to zero")
         IF (pao%linpot_regu_strength /= 0.0_dp) &
            CPABORT("PAO forces require LINPOT_REGULARIZATION_STRENGTH or MAX_PAO set to zero")
         IF (pao%regularization /= 0.0_dp) &
            CPABORT("PAO forces require REGULARIZATION or MAX_PAO set to zero")
      ENDIF

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      natom=natoms)

      ALLOCATE (forces(natoms, 3))
      forces(:, :) = 0.0_dp
      CALL pao_calc_outer_grad_lnv(qs_env, ls_scf_env, matrix_M)
      CALL pao_calc_U(pao, qs_env, matrix_M, pao%matrix_G, forces=forces) ! without penalty terms
      CALL cp_dbcsr_release(matrix_M)

      IF (SIZE(pao%ml_training_set) > 0) &
         CALL pao_ml_forces(pao, qs_env, pao%matrix_G, forces)

      CALL mp_sum(forces, para_env%group)
      DO iatom = 1, natoms
         particle_set(iatom)%f = particle_set(iatom)%f+forces(iatom, :)
      ENDDO

      DEALLOCATE (forces)

      CALL timestop(handle)

   END SUBROUTINE pao_add_forces

END MODULE pao_methods
