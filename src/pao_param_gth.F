!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Parametrization based on GTH pseudo potentials
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_gth
   USE ai_overlap,                      ONLY: overlap_aab
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_arnoldi_extremal, cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_get_block_p, &
        cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
        cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, &
        cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_reserve_all_blocks, &
        cp_dbcsr_reserve_diag_blocks, cp_dbcsr_row_block_sizes, cp_dbcsr_set, cp_dbcsr_type, &
        dbcsr_distribution_mp, dbcsr_mp_group
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gamma1
   USE mathlib,                         ONLY: multinomial
   USE message_passing,                 ONLY: mp_min,&
                                              mp_sum
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              ncoset
   USE pao_param_fock,                  ONLY: pao_calc_U_block_fock,&
                                              pao_param_finalize_fock,&
                                              pao_param_init_fock
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              pao_potential_type,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_param_init_gth, pao_param_finalize_gth, pao_calc_U_gth
   PUBLIC :: pao_param_count_gth, pao_param_initguess_gth

CONTAINS

! **************************************************************************************************
!> \brief Initialize the linear potential parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init_gth(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init_gth'

      INTEGER                                            :: acol, arow, handle, iatom, idx, ikind, &
                                                            iterm, jatom, n, natoms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pri, col_blk_size, nterms, &
                                                            row_blk_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_V_term, vec_V_terms
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      CALL pao_param_init_fock(pao, qs_env)

      CALL get_qs_env(qs_env, &
                      matrix_s=matrix_s, &
                      particle_set=particle_set, &
                      natom=natoms)

      ALLOCATE (row_blk_size(natoms), col_blk_size(natoms), nterms(natoms))
      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL pao_param_count_gth(qs_env, ikind, nterms(iatom))
      ENDDO

      ! allocate matrix_V_terms
      blk_sizes_pri => cp_dbcsr_row_block_sizes(matrix_s(1)%matrix)
      col_blk_size = SUM(nterms)
      row_blk_size = blk_sizes_pri**2
      CALL cp_dbcsr_init(pao%matrix_V_terms)
      CALL cp_dbcsr_create(pao%matrix_V_terms, &
                           template=matrix_s(1)%matrix, &
                           matrix_type="N", &
                           row_blk_size=row_blk_size, &
                           col_blk_size=col_blk_size)
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_V_terms)
      CALL cp_dbcsr_set(pao%matrix_V_terms, 0.0_dp)

      ! calculate and store poential terms
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_V_terms)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, vec_V_terms)
         iatom = arow; CPASSERT(arow == acol)
         n = blk_sizes_pri(iatom)
         DO jatom = 1, natoms
            IF (jatom == iatom) CYCLE ! waste some storage to simplify things later
            DO iterm = 1, nterms(jatom)
               idx = SUM(nterms(1:jatom-1))+iterm
               block_V_term(1:n, 1:n) => vec_V_terms(:, idx) ! map column into matrix
               CALL gth_calc_term(qs_env, block_V_term, iatom, jatom, iterm)
            ENDDO
         ENDDO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      IF (pao%precondition) &
         CALL pao_param_gth_preconditioner(pao, qs_env, nterms)

      DEALLOCATE (row_blk_size, col_blk_size, nterms)
      CALL timestop(handle)
   END SUBROUTINE pao_param_init_gth

! **************************************************************************************************
!> \brief Finalize the GTH potential parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_gth(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CALL pao_param_finalize_fock(pao)
      CALL cp_dbcsr_release(pao%matrix_V_terms)
      IF (pao%precondition) THEN
         CALL cp_dbcsr_release(pao%matrix_precon)
         CALL cp_dbcsr_release(pao%matrix_precon_inv)
      ENDIF

   END SUBROUTINE pao_param_finalize_gth

! **************************************************************************************************
!> \brief Builds the preconditioner matrix_precon and matrix_precon_inv
!> \param pao ...
!> \param qs_env ...
!> \param nterms ...
! **************************************************************************************************
   SUBROUTINE pao_param_gth_preconditioner(pao, qs_env, nterms)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, DIMENSION(:), POINTER                     :: nterms

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_gth_preconditioner'

      INTEGER                                            :: acol, arow, group, handle, i, iatom, &
                                                            ioffset, j, jatom, joffset, m, n, &
                                                            natoms
      LOGICAL                                            :: arnoldi_converged, converged, found
      REAL(dp)                                           :: eval_max, eval_min
      REAL(dp), DIMENSION(:, :), POINTER                 :: block, block_overlap, block_V_term
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_type)                                :: matrix_gth_overlap
      TYPE(ls_scf_env_type), POINTER                     :: ls_scf_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, ls_scf_env=ls_scf_env)
      group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(pao%matrix_V_terms)))
      natoms = SIZE(nterms)

      CALL cp_dbcsr_init(matrix_gth_overlap)
      CALL cp_dbcsr_create(matrix_gth_overlap, &
                           template=pao%matrix_V_terms, &
                           matrix_type="N", &
                           row_blk_size=nterms, &
                           col_blk_size=nterms)
      CALL cp_dbcsr_reserve_all_blocks(matrix_gth_overlap)
      CALL cp_dbcsr_set(matrix_gth_overlap, 0.0_dp)

      DO iatom = 1, natoms
      DO jatom = 1, natoms
         ioffset = SUM(nterms(1:iatom-1))
         joffset = SUM(nterms(1:jatom-1))
         n = nterms(iatom)
         m = nterms(jatom)

         ALLOCATE (block(n, m))
         block = 0.0_dp

         CALL cp_dbcsr_iterator_start(iter, pao%matrix_V_terms)
         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
            CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_V_term)
            CPASSERT(arow == acol)
            DO i = 1, n
            DO j = 1, m
               block(i, j) = block(i, j)+SUM(block_V_term(:, ioffset+i)*block_V_term(:, joffset+j))
            ENDDO
            ENDDO
         ENDDO
         CALL cp_dbcsr_iterator_stop(iter)

         CALL mp_sum(block, group)

         CALL cp_dbcsr_get_block_p(matrix=matrix_gth_overlap, row=iatom, col=jatom, block=block_overlap, found=found)
         IF (ASSOCIATED(block_overlap)) &
            block_overlap = block

         DEALLOCATE (block)
      ENDDO
      ENDDO

      !TODO: good setting for arnoldi?
      CALL cp_dbcsr_arnoldi_extremal(matrix_gth_overlap, eval_max, eval_min, max_iter=100, &
                                     threshold=1e-2_dp, converged=arnoldi_converged)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| GTH-preconditioner converged, min, max, max/min:", &
         arnoldi_converged, eval_min, eval_max, eval_max/eval_min

      CALL cp_dbcsr_init(pao%matrix_precon)
      CALL cp_dbcsr_create(pao%matrix_precon, template=matrix_gth_overlap)
      CALL cp_dbcsr_init(pao%matrix_precon_inv)
      CALL cp_dbcsr_create(pao%matrix_precon_inv, template=matrix_gth_overlap)

      CALL matrix_sqrt_Newton_Schulz(pao%matrix_precon_inv, pao%matrix_precon, matrix_gth_overlap, &
                                     threshold=ls_scf_env%eps_filter, &
                                     order=ls_scf_env%sign_sqrt_order, &
                                     max_iter_lanczos=ls_scf_env%max_iter_lanczos, &
                                     eps_lanczos=ls_scf_env%eps_lanczos, &
                                     converged=converged)
      CALL cp_dbcsr_release(matrix_gth_overlap)

      IF (.NOT. converged) &
         CPABORT("PAO: Sqrt of GTH-preconditioner did not converge.")

      CALL timestop(handle)
   END SUBROUTINE pao_param_gth_preconditioner

! **************************************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param penalty ...
!> \param matrix_M1 ...
!> \param matrix_G ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U_gth(pao, penalty, matrix_M1, matrix_G)
      TYPE(pao_env_type), POINTER                        :: pao
      REAL(dp), INTENT(INOUT)                            :: penalty
      TYPE(cp_dbcsr_type), OPTIONAL                      :: matrix_M1, matrix_G

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_U_gth'

      INTEGER                                            :: acol, arow, group, handle, iatom, idx, &
                                                            iterm, n, natoms
      INTEGER, DIMENSION(:), POINTER                     :: nterms
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:), POINTER                    :: world_G, world_X
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G, block_M1, block_M2, block_U, &
                                                            block_V, block_V_term, block_X, &
                                                            vec_V_terms
      REAL(KIND=dp)                                      :: min_gap
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      min_gap = HUGE(dp)
      nterms => cp_dbcsr_row_block_sizes(pao%matrix_X)
      natoms = SIZE(nterms)
      group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(pao%matrix_X)))

      ! allocate arrays for world-view
      ALLOCATE (world_X(SUM(nterms)), world_G(SUM(nterms)))
      world_X = 0.0_dp; world_G = 0.0_dp

      ! collect world_X from atomic blocks
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         idx = SUM(nterms(1:iatom-1))
         world_X(idx+1:idx+nterms(iatom)) = block_X(:, 1)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)
      CALL mp_sum(world_X, group) ! sync world view across MPI ranks

      ! loop over atoms
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_U)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_U)
         iatom = arow; CPASSERT(arow == acol)
         n = SIZE(block_U, 1)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=vec_V_terms, found=found)
         CPASSERT(ASSOCIATED(vec_V_terms))

         ! calculate potential V of i'th atom
         ALLOCATE (block_V(n, n))
         block_V = 0.0_dp
         DO iterm = 1, SIZE(world_X)
            block_V_term(1:n, 1:n) => vec_V_terms(:, iterm) ! map column into matrix
            block_V = block_V+world_X(iterm)*block_V_term
         ENDDO

         ! calculate gradient block of i'th atom
         IF (.NOT. PRESENT(matrix_G)) THEN
            CALL pao_calc_U_block_fock(pao, iatom, penalty, min_gap, block_V, block_U)

         ELSE ! TURNING POINT (if calc grad) ------------------------------------
            CPASSERT(PRESENT(matrix_M1))
            CALL cp_dbcsr_get_block_p(matrix=matrix_M1, row=iatom, col=iatom, block=block_M1, found=found)
            ALLOCATE (block_M2(n, n))
            CALL pao_calc_U_block_fock(pao, iatom, penalty, min_gap, block_V, block_U, block_M1, block_M2)
            DO iterm = 1, SIZE(world_G)
               block_V_term(1:n, 1:n) => vec_V_terms(:, iterm) ! map column into matrix
               world_G(iterm) = world_G(iterm)+SUM(block_V_term*block_M2)
            ENDDO
            DEALLOCATE (block_M2)
         ENDIF
         DEALLOCATE (block_V)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      ! distribute world_G across atomic blocks
      IF (PRESENT(matrix_G)) THEN
         CALL mp_sum(world_G, group) ! sync world view across MPI ranks
         CALL cp_dbcsr_iterator_start(iter, matrix_G)
         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
            CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_G)
            iatom = arow; CPASSERT(arow == acol)
            idx = SUM(nterms(1:iatom-1))
            block_G(:, 1) = world_G(idx+1:idx+nterms(iatom))
         ENDDO
         CALL cp_dbcsr_iterator_stop(iter)
      ENDIF

      DEALLOCATE (world_X, world_G)

      ! sum penalty energies across ranks
      CALL mp_sum(penalty, group)

      ! print smallest homo-lumo gap encountered by fock-layer
      CALL mp_min(min_gap, group)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| min_gap:", min_gap

      CALL timestop(handle)

   END SUBROUTINE pao_calc_U_gth

! **************************************************************************************************
!> \brief Returns the number of parameters for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count_gth(qs_env, ikind, nparams)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      INTEGER                                            :: max_projector, maxl, ncombis
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), pao_potentials=pao_potentials)

      IF (SIZE(pao_potentials) /= 1) &
         CPABORT("GTH parametrization requires exactly one PAO_POTENTIAL section per KIND")

      max_projector = pao_potentials(1)%max_projector
      maxl = pao_potentials(1)%maxl

      IF (maxl < 0) &
         CPABORT("GTH parametrization requires non-negative PAO_POTENTIAL%MAXL")

      IF (max_projector < 0) &
         CPABORT("GTH parametrization requires non-negative PAO_POTENTIAL%MAX_PROJECTOR")

      IF (MOD(maxl, 2) /= 0) &
         CPABORT("GTH parametrization requires even-numbered PAO_POTENTIAL%MAXL")

      ncombis = (max_projector+1)*(max_projector+2)/2
      nparams = ncombis*(maxl/2+1)

   END SUBROUTINE pao_param_count_gth

! **************************************************************************************************
!> \brief Fills the given block_V with the requested potential term
!> \param qs_env ...
!> \param block_V ...
!> \param iatom ...
!> \param jatom ...
!> \param kterm ...
! **************************************************************************************************
   SUBROUTINE gth_calc_term(qs_env, block_V, iatom, jatom, kterm)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V
      INTEGER, INTENT(IN)                                :: iatom, jatom, kterm

      INTEGER                                            :: c, jkind, lpot, max_shell, min_shell, &
                                                            pot_max_projector, pot_maxl
      REAL(KIND=dp)                                      :: beta, norm2, pot_beta, R_vdw
      TYPE(pao_potential_type), DIMENSION(:), POINTER    :: pao_potentials
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set)
      CALL get_atomic_kind(particle_set(jatom)%atomic_kind, rvdw=R_vdw, kind_number=jkind)
      CALL get_qs_kind(qs_kind_set(jkind), pao_potentials=pao_potentials)

      IF (SIZE(pao_potentials) /= 1) &
         CPABORT("GTH parametrization requires exactly one PAO_POTENTIAL section per KIND")

      pot_max_projector = pao_potentials(1)%max_projector
      pot_maxl = pao_potentials(1)%maxl
      pot_beta = pao_potentials(1)%beta

      c = 0
      outer: &
         DO lpot = 0, pot_maxl, 2
      DO max_shell = 0, pot_max_projector
      DO min_shell = 0, max_shell
         c = c+1
         IF (c == kterm) EXIT outer
      ENDDO
      ENDDO
      ENDDO outer

      beta = 0.5_dp/(R_vdw*pot_beta)**2

      CALL calc_gaussian(qs_env, block_V, iatom, jatom, lpot=lpot, beta=beta, &
                         min_shell=min_shell, max_shell=max_shell)

      ! normalize
      norm2 = (2.0_dp*beta)**(-0.5_dp-lpot)*gamma1(lpot)
      block_V = block_V/SQRT(norm2)

   END SUBROUTINE gth_calc_term

! **************************************************************************************************
!> \brief Calculate intial guess for matrix_X
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_gth(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      INTEGER                                            :: acol, arow
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
         CPASSERT(arow == acol)
         CPASSERT(SIZE(block_X, 2) == 1)

         ! a simplistic guess, which at least makes the atom visible to others
         block_X = 0.0_dp
         block_X(1, 1) = 0.01_dp
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

   END SUBROUTINE pao_param_initguess_gth

! **************************************************************************************************
!> \brief Calculates potential term of the form r**lpot * Exp(-beta*r**2)
!> \param qs_env ...
!> \param block_V potential term that is returned
!> \param iatom central atom
!> \param jatom remote atom on which the Gaussian is centered
!> \param lpot polynomial prefactor, r**lpot
!> \param beta exponent of the Gaussian
!> \param min_shell ...
!> \param max_shell ...
! **************************************************************************************************
   SUBROUTINE calc_gaussian(qs_env, block_V, iatom, jatom, lpot, beta, min_shell, max_shell)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V
      INTEGER, INTENT(IN)                                :: iatom, jatom, lpot
      REAL(dp), INTENT(IN)                               :: beta
      INTEGER, INTENT(IN)                                :: min_shell, max_shell

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_gaussian'
      REAL(dp), PARAMETER                                :: inf_dist = 0.1*HUGE(1.0)

      INTEGER :: handle, ic, ikind, iset, ishell, jset, jshell, la1_max, la1_min, la2_max, &
         la2_min, lb_max, lb_min, na1, na2, natoms, nb, ncfga1, ncfga2, ncfgb, npgfa1, npgfa2, &
         npgfb
      REAL(dp), DIMENSION(:), POINTER                    :: rpgfa1, rpgfa2, zeta1, zeta2
      REAL(dp), DIMENSION(:, :), POINTER                 :: sab
      REAL(KIND=dp)                                      :: coeff
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rpgfb, zetb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: saab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      CALL init_orbital_pointers(lpot)

      CALL get_qs_env(qs_env, &
                      natom=natoms, &
                      cell=cell, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)

      ! setup description of potential
      lb_min = lpot
      lb_max = lpot
      ncfgb = ncoset(lb_max)-ncoset(lb_min-1)
      npgfb = 1 ! number of exponents
      nb = npgfb*ncfgb

      ! initialize exponents, using R_vdw from neighbour atom
      ALLOCATE (rpgfb(npgfb), zetb(npgfb))
      rpgfb(1) = inf_dist ! range of functions, we dont want screening
      zetb(1) = beta

      ! distance vector A-B
      ra = particle_set(iatom)%r
      rb = particle_set(jatom)%r
      rab = pbc(ra, rb, cell)

      ! loop over all set/shell combination and fill block_V
      block_V(:, :) = 0.0_dp
      DO iset = 1, basis_set%nset
      DO jset = 1, basis_set%nset
      DO ishell = 1, basis_set%nshell(iset)
      DO jshell = 1, basis_set%nshell(jset)
         IF (MIN(basis_set%l(ishell, iset), basis_set%l(jshell, jset)) /= min_shell) CYCLE
         IF (MAX(basis_set%l(ishell, iset), basis_set%l(jshell, jset)) /= max_shell) CYCLE

         ! setup iset
         la1_max = basis_set%l(ishell, iset)
         la1_min = basis_set%l(ishell, iset)
         npgfa1 = basis_set%npgf(iset)
         ncfga1 = ncoset(la1_max)-ncoset(la1_min-1)
         na1 = npgfa1*ncfga1
         zeta1 => basis_set%zet(:, iset)
         !rpgfa1  => basis_set%pgf_radius(:,iset) !TODO: reenable screening
         ALLOCATE (rpgfa1(npgfa1))
         rpgfa1(:) = inf_dist ! range of functions, we dont want screening

         ! setup jset
         la2_max = basis_set%l(jshell, jset)
         la2_min = basis_set%l(jshell, jset)
         npgfa2 = basis_set%npgf(jset)
         ncfga2 = ncoset(la2_max)-ncoset(la2_min-1)
         na2 = npgfa2*ncfga2
         zeta2 => basis_set%zet(:, jset)
         !rpgfa2  => basis_set%pgf_radius(:,jset) !TODO: reenable screening
         ALLOCATE (rpgfa2(npgfa2))
         rpgfa2(:) = inf_dist ! range of functions, we dont want screening

         ! calculate integrals
         ALLOCATE (saab(na1, na2, nb))
         saab = 0.0_dp
         CALL overlap_aab(la1_max, la1_min, npgfa1, rpgfa1, zeta1, &
                          la2_max, la2_min, npgfa2, rpgfa2, zeta2, &
                          lb_max, lb_min, npgfb, rpgfb, zetb, &
                          rab=rab, saab=saab)
         DEALLOCATE (rpgfa1, rpgfa2)

         ! sum potential terms: POW(x**2 + y**2 + z**2, lpot/2)
         ALLOCATE (sab(na1, na2))
         sab = 0.0_dp
         DO ic = 1, ncfgb
            coeff = multinomial(lpot/2, indco(:, ncoset(lpot-1)+ic)/2)
            sab = sab+coeff*saab(:, :, ic)
         ENDDO

         CALL my_contract(sab, block_V, basis_set, iset, ishell, jset, jshell)
         DEALLOCATE (sab, saab)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DEALLOCATE (rpgfb, zetb)

      ! symmetrize
      block_V = 0.5_dp*(block_V+TRANSPOSE(block_V))

      CALL timestop(handle)
   END SUBROUTINE calc_gaussian

! **************************************************************************************************
!> \brief Helper routine, contracts a basis block
!> \param sab ...
!> \param block_V ...
!> \param basis_set ...
!> \param iset ...
!> \param ishell ...
!> \param jset ...
!> \param jshell ...
! **************************************************************************************************
   SUBROUTINE my_contract(sab, block_V, basis_set, iset, ishell, jset, jshell)
      REAL(dp), DIMENSION(:, :), POINTER                 :: sab
      REAL(dp), DIMENSION(:, :), INTENT(OUT), TARGET     :: block_V
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      INTEGER, INTENT(IN)                                :: iset, ishell, jset, jshell

      INTEGER                                            :: a, b, c, d, ipgf, jpgf, l1, l2, n1, n2, &
                                                            nn1, nn2, sgfa1, sgfa2, sgla1, sgla2
      REAL(dp), DIMENSION(:, :), POINTER                 :: S, T1, T2, V

! first and last indices of given shell in block_V.
! This matrix is in the contracted spherical basis.

      sgfa1 = basis_set%first_sgf(ishell, iset)
      sgla1 = basis_set%last_sgf(ishell, iset)
      sgfa2 = basis_set%first_sgf(jshell, jset)
      sgla2 = basis_set%last_sgf(jshell, jset)

      ! prepare the result matrix
      V => block_V(sgfa1:sgla1, sgfa2:sgla2)

      ! Calculate strides of sphi matrix.
      ! This matrix is in the uncontraced cartesian basis.
      ! It contains all shells of the set.
      ! It's index runs over all primitive gaussians of the set
      ! and then for each gaussian over all configurations of *the entire set*. (0->lmax)
      nn1 = ncoset(basis_set%lmax(iset))
      nn2 = ncoset(basis_set%lmax(jset))

      ! Calculate strides of sab matrix
      ! This matrix is also in the uncontraced cartensian basis,
      ! however it contains only a single shell.
      ! It's index runs over all primitive gaussians of the set
      ! and then for each gaussian over all configrations of *the given shell*.
      l1 = basis_set%l(ishell, iset)
      l2 = basis_set%l(jshell, jset)
      n1 = ncoset(l1)-ncoset(l1-1)
      n2 = ncoset(l2)-ncoset(l2-1)

      DO ipgf = 1, basis_set%npgf(iset)
      DO jpgf = 1, basis_set%npgf(jset)
         ! prepare first trafo-matrix
         a = (ipgf-1)*nn1+ncoset(l1-1)+1
         T1 => basis_set%sphi(a:a+n1-1, sgfa1:sgla1)

         ! prepare second trafo-matrix
         b = (jpgf-1)*nn2+ncoset(l2-1)+1
         T2 => basis_set%sphi(b:b+n2-1, sgfa2:sgla2)

         ! prepare SAB matrix
         c = (ipgf-1)*n1+1
         d = (jpgf-1)*n2+1
         S => sab(c:c+n1-1, d:d+n2-1)

         ! do the transformation
         V = V+MATMUL(TRANSPOSE(T1), MATMUL(S, T2))
      ENDDO
      ENDDO

   END SUBROUTINE my_contract

END MODULE pao_param_gth
