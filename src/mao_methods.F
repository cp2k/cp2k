!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate MAO's and analyze wavefunctions
!> \par History
!>      03.2016 created [JGH]
!>      12.2016 split into four modules [JGH]
!> \author JGH
! **************************************************************************************************
MODULE mao_methods
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: create_primitive_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type,&
                                              write_gto_basis_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_allocate_matrix_set, cp_dbcsr_create, cp_dbcsr_desymmetrize, &
        cp_dbcsr_distribution, cp_dbcsr_get_block_p, cp_dbcsr_get_info, cp_dbcsr_init, &
        cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
        cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
        cp_dbcsr_release, cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type, dbcsr_distribution_mp, &
        dbcsr_mp_group, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_plus_fm_fm_t
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: mao_basis_ext,&
                                              mao_basis_orb,&
                                              mao_basis_prim
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lapack,                          ONLY: lapack_ssyev,&
                                              lapack_ssygv
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_methods'

   PUBLIC ::  mao_initialization, mao_function, mao_function_gradient, mao_orthogonalization, &
             mao_project_gradient, mao_scalar_product, mao_build_q, mao_basis_analysis, &
             mao_reference_basis, calculate_p_gamma

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param pmat ...
!> \param smat ...
! **************************************************************************************************
   SUBROUTINE mao_initialization(mao_coef, pmat, smat)
      TYPE(cp_dbcsr_type)                                :: mao_coef, pmat, smat

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_initialization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, info, jatom, lwork, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, bmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cblock, pblock, sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, mao_coef)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, cblock)
         CPASSERT(iatom == jatom)
         m = SIZE(cblock, 2)
         NULLIFY (pblock, sblock)
         CALL cp_dbcsr_get_block_p(matrix=pmat, row=iatom, col=jatom, block=pblock, found=found)
         CPASSERT(found)
         CALL cp_dbcsr_get_block_p(matrix=smat, row=iatom, col=jatom, block=sblock, found=found)
         CPASSERT(found)
         n = SIZE(sblock, 1)
         lwork = MAX(n*n, 100)
         ALLOCATE (amat(n, n), bmat(n, n), w(n), work(lwork))
         amat(1:n, 1:n) = pblock(1:n, 1:n)
         bmat(1:n, 1:n) = sblock(1:n, 1:n)
         info = 0
         CALL lapack_ssygv(1, "V", "U", n, amat, n, bmat, n, w, work, lwork, info)
         CPASSERT(info == 0)
         cblock(1:n, 1:m) = amat(1:n, n:n-m+1:-1)
         DEALLOCATE (amat, bmat, w, work)
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

      CALL mao_orthogonalization(mao_coef, smat)

   END SUBROUTINE mao_initialization

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param fval ...
!> \param qmat ...
!> \param smat ...
!> \param binv ...
!> \param reuse ...
! **************************************************************************************************
   SUBROUTINE mao_function(mao_coef, fval, qmat, smat, binv, reuse)
      TYPE(cp_dbcsr_type)                                :: mao_coef
      REAL(KIND=dp), INTENT(OUT)                         :: fval
      TYPE(cp_dbcsr_type)                                :: qmat, smat, binv
      LOGICAL, INTENT(IN)                                :: reuse

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_function', routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: convergence, threshold
      TYPE(cp_dbcsr_type)                                :: bmat, scmat, tmat

      threshold = 1.e-8_dp
      convergence = 1.e-6_dp
      ! temp matrices
      CALL cp_dbcsr_init(scmat)
      CALL cp_dbcsr_create(scmat, template=mao_coef)
      CALL cp_dbcsr_init(bmat)
      CALL cp_dbcsr_create(bmat, template=binv)
      CALL cp_dbcsr_init(tmat)
      CALL cp_dbcsr_create(tmat, template=qmat)
      ! calculate B=C(T)*S*C matrix, S=(MAO,MAO) overlap
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, smat, mao_coef, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mao_coef, scmat, 0.0_dp, bmat)
      ! calculate inverse of B
      CALL invert_Hotelling(binv, bmat, threshold, use_inv_as_guess=reuse, &
                            norm_convergence=convergence, silent=.TRUE.)
      ! calculate Binv*C and T=C(T)*Binv*C
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mao_coef, binv, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, scmat, mao_coef, 0.0_dp, tmat)
      ! function = Tr(Q*T)
      CALL cp_dbcsr_trace(qmat, tmat, fval)
      ! free temp matrices
      CALL cp_dbcsr_release(scmat)
      CALL cp_dbcsr_release(bmat)
      CALL cp_dbcsr_release(tmat)

   END SUBROUTINE mao_function

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param fval ...
!> \param mao_grad ...
!> \param qmat ...
!> \param smat ...
!> \param binv ...
!> \param reuse ...
! **************************************************************************************************
   SUBROUTINE mao_function_gradient(mao_coef, fval, mao_grad, qmat, smat, binv, reuse)
      TYPE(cp_dbcsr_type)                                :: mao_coef
      REAL(KIND=dp), INTENT(OUT)                         :: fval
      TYPE(cp_dbcsr_type)                                :: mao_grad, qmat, smat, binv
      LOGICAL, INTENT(IN)                                :: reuse

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_function_gradient', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: convergence, threshold
      TYPE(cp_dbcsr_type)                                :: bmat, scmat, t2mat, tmat

      threshold = 1.e-8_dp
      convergence = 1.e-6_dp
      ! temp matrices
      CALL cp_dbcsr_init(scmat)
      CALL cp_dbcsr_create(scmat, template=mao_coef)
      CALL cp_dbcsr_init(bmat)
      CALL cp_dbcsr_create(bmat, template=binv)
      CALL cp_dbcsr_init(tmat)
      CALL cp_dbcsr_create(tmat, template=qmat)
      CALL cp_dbcsr_init(t2mat)
      CALL cp_dbcsr_create(t2mat, template=scmat)
      ! calculate B=C(T)*S*C matrix, S=(MAO,MAO) overlap
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, smat, mao_coef, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mao_coef, scmat, 0.0_dp, bmat)
      ! calculate inverse of B
      CALL invert_Hotelling(binv, bmat, threshold, use_inv_as_guess=reuse, &
                            norm_convergence=convergence, silent=.TRUE.)
      ! calculate R=C*Binv and T=C*Binv*C(T)=R*C(T)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mao_coef, binv, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, scmat, mao_coef, 0.0_dp, tmat)
      ! function = Tr(Q*T)
      CALL cp_dbcsr_trace(qmat, tmat, fval)
      ! Gradient part 1: g = 2*Q*C*Binv = 2*Q*R
      CALL cp_dbcsr_multiply("N", "N", 2.0_dp, qmat, scmat, 0.0_dp, mao_grad, &
                             retain_sparsity=.TRUE.)
      ! Gradient part 2: g = -2*S*T*X; X = Q*R
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, qmat, scmat, 0.0_dp, t2mat)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmat, t2mat, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("N", "N", -2.0_dp, smat, scmat, 1.0_dp, mao_grad, &
                             retain_sparsity=.TRUE.)
      ! free temp matrices
      CALL cp_dbcsr_release(scmat)
      CALL cp_dbcsr_release(bmat)
      CALL cp_dbcsr_release(tmat)
      CALL cp_dbcsr_release(t2mat)

      CALL mao_project_gradient(mao_coef, mao_grad, smat)

   END SUBROUTINE mao_function_gradient

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param smat ...
! **************************************************************************************************
   SUBROUTINE mao_orthogonalization(mao_coef, smat)
      TYPE(cp_dbcsr_type)                                :: mao_coef, smat

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_orthogonalization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, info, jatom, lwork, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, bmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cblock, sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, mao_coef)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, cblock)
         CPASSERT(iatom == jatom)
         m = SIZE(cblock, 2)
         n = SIZE(cblock, 1)
         NULLIFY (sblock)
         CALL cp_dbcsr_get_block_p(matrix=smat, row=iatom, col=jatom, block=sblock, found=found)
         CPASSERT(found)
         lwork = MAX(n*n, 100)
         ALLOCATE (amat(n, m), bmat(m, m), w(m), work(lwork))
         amat(1:n, 1:m) = MATMUL(sblock(1:n, 1:n), cblock(1:n, 1:m))
         bmat(1:m, 1:m) = MATMUL(TRANSPOSE(cblock(1:n, 1:m)), amat(1:n, 1:m))
         info = 0
         CALL lapack_ssyev("V", "U", m, bmat, m, w, work, lwork, info)
         CPASSERT(info == 0)
         CPASSERT(ALL(w > 0.0_dp))
         w = 1.0_dp/SQRT(w)
         DO i = 1, m
            amat(1:m, i) = bmat(1:m, i)*w(i)
         END DO
         bmat(1:m, 1:m) = MATMUL(amat(1:m, 1:m), TRANSPOSE(bmat(1:m, 1:m)))
         cblock(1:n, 1:m) = MATMUL(cblock(1:n, 1:m), bmat(1:m, 1:m))
         DEALLOCATE (amat, bmat, w, work)
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE mao_orthogonalization

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param mao_grad ...
!> \param smat ...
! **************************************************************************************************
   SUBROUTINE mao_project_gradient(mao_coef, mao_grad, smat)
      TYPE(cp_dbcsr_type)                                :: mao_coef, mao_grad, smat

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_project_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, jatom, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cblock, gblock, sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, mao_coef)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, cblock)
         CPASSERT(iatom == jatom)
         m = SIZE(cblock, 2)
         n = SIZE(cblock, 1)
         NULLIFY (sblock)
         CALL cp_dbcsr_get_block_p(matrix=smat, row=iatom, col=jatom, block=sblock, found=found)
         CPASSERT(found)
         NULLIFY (gblock)
         CALL cp_dbcsr_get_block_p(matrix=mao_grad, row=iatom, col=jatom, block=gblock, found=found)
         CPASSERT(found)
         ALLOCATE (amat(m, m))
         amat(1:m, 1:m) = MATMUL(TRANSPOSE(cblock(1:n, 1:m)), MATMUL(sblock(1:n, 1:n), gblock(1:n, 1:m)))
         gblock(1:n, 1:m) = gblock(1:n, 1:m)-MATMUL(cblock(1:n, 1:m), amat(1:m, 1:m))
         DEALLOCATE (amat)
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE mao_project_gradient

! **************************************************************************************************
!> \brief ...
!> \param fmat1 ...
!> \param fmat2 ...
!> \retval spro ...
! **************************************************************************************************
   FUNCTION mao_scalar_product(fmat1, fmat2) RESULT(spro)
      TYPE(cp_dbcsr_type)                                :: fmat1, fmat2
      REAL(KIND=dp)                                      :: spro

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_scalar_product', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, iatom, jatom, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ablock, bblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      spro = 0.0_dp

      CALL cp_dbcsr_iterator_start(dbcsr_iter, fmat1)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, ablock)
         CPASSERT(iatom == jatom)
         m = SIZE(ablock, 2)
         n = SIZE(ablock, 1)
         CALL cp_dbcsr_get_block_p(matrix=fmat2, row=iatom, col=jatom, block=bblock, found=found)
         CPASSERT(found)
         spro = spro+SUM(ablock(1:n, 1:m)*bblock(1:n, 1:m))
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

      group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(fmat1)))
      CALL mp_sum(spro, group)

   END FUNCTION mao_scalar_product

! **************************************************************************************************
!> \brief Calculate the density matrix at the Gamma point
!> \param pmat ...
!> \param ksmat ...
!> \param smat ...
!> \param kpoints      Kpoint environment
!> \param nmos         Number of occupied orbitals
!> \param occ          Maximum occupation per orbital
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE calculate_p_gamma(pmat, ksmat, smat, kpoints, nmos, occ)

      TYPE(cp_dbcsr_type)                                :: pmat, ksmat, smat
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: nmos
      REAL(KIND=dp), INTENT(IN)                          :: occ

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_p_gamma', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: norb
      REAL(KIND=dp)                                      :: de
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_dbcsr_type)                                :: tempmat
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: fmksmat, fmsmat, fmvec, fmwork

      ! FM matrices

      CALL cp_dbcsr_get_info(smat, nfullrows_total=norb)
      CALL cp_fm_struct_create(fmstruct=matrix_struct, context=kpoints%blacs_env_all, &
                               nrow_global=norb, ncol_global=norb)
      CALL cp_fm_create(fmksmat, matrix_struct)
      CALL cp_fm_create(fmsmat, matrix_struct)
      CALL cp_fm_create(fmvec, matrix_struct)
      CALL cp_fm_create(fmwork, matrix_struct)
      ALLOCATE (eigenvalues(norb))

      ! DBCSR matrix
      CALL cp_dbcsr_init(tempmat)
      CALL cp_dbcsr_create(tempmat, template=smat, matrix_type=dbcsr_type_no_symmetry)

      ! transfer to FM
      CALL cp_dbcsr_desymmetrize(smat, tempmat)
      CALL copy_dbcsr_to_fm(tempmat, fmsmat)
      CALL cp_dbcsr_desymmetrize(ksmat, tempmat)
      CALL copy_dbcsr_to_fm(tempmat, fmksmat)

      ! diagonalize
      CALL cp_fm_geeig(fmksmat, fmsmat, fmvec, eigenvalues, fmwork)
      de = eigenvalues(nmos+1)-eigenvalues(nmos)
      IF (de < 0.001_dp) THEN
         CALL cp_warn(__LOCATION__, "MAO: No band gap at "// &
                      "Gamma point. MAO analysis not reliable.")
      END IF
      ! density matrix
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=pmat, matrix_v=fmvec, ncol=nmos, alpha=occ)

      DEALLOCATE (eigenvalues)
      CALL cp_dbcsr_release(tempmat)
      CALL cp_fm_release(fmksmat)
      CALL cp_fm_release(fmsmat)
      CALL cp_fm_release(fmvec)
      CALL cp_fm_release(fmwork)
      CALL cp_fm_struct_release(matrix_struct)

   END SUBROUTINE calculate_p_gamma

! **************************************************************************************************
!> \brief Define the MAO reference basis set
!> \param qs_env ...
!> \param mao_basis ...
!> \param mao_basis_set_list ...
!> \param orb_basis_set_list ...
!> \param iunit ...
!> \param print_basis ...
!> \par History
!>      07.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE mao_reference_basis(qs_env, mao_basis, mao_basis_set_list, orb_basis_set_list, &
                                  iunit, print_basis)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: mao_basis
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: mao_basis_set_list, orb_basis_set_list
      INTEGER, INTENT(IN), OPTIONAL                      :: iunit
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_basis

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_reference_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, nbas, nkind, unit_nr
      REAL(KIND=dp)                                      :: eps_pgf_orb
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, pbasis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      ! Reference basis set
      CPASSERT(.NOT. ASSOCIATED(mao_basis_set_list))
      CPASSERT(.NOT. ASSOCIATED(orb_basis_set_list))

      ! options
      IF (PRESENT(iunit)) THEN
         unit_nr = iunit
      ELSE
         unit_nr = -1
      END IF

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (mao_basis_set_list(nkind), orb_basis_set_list(nkind))
      DO ikind = 1, nkind
         NULLIFY (mao_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (orb_basis_set_list(ikind)%gto_basis_set)
      END DO
      !
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
         IF (ASSOCIATED(basis_set)) orb_basis_set_list(ikind)%gto_basis_set => basis_set
      END DO
      !
      SELECT CASE (mao_basis)
      CASE (mao_basis_orb)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
            IF (ASSOCIATED(basis_set)) mao_basis_set_list(ikind)%gto_basis_set => basis_set
         END DO
      CASE (mao_basis_prim)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
            NULLIFY (pbasis)
            IF (ASSOCIATED(basis_set)) THEN
               CALL create_primitive_basis_set(basis_set, pbasis)
               CALL get_qs_env(qs_env, dft_control=dft_control)
               eps_pgf_orb = dft_control%qs_control%eps_pgf_orb
               CALL init_interaction_radii_orb_basis(pbasis, eps_pgf_orb)
               pbasis%kind_radius = basis_set%kind_radius
               mao_basis_set_list(ikind)%gto_basis_set => pbasis
               CALL add_basis_set_to_container(qs_kind%basis_sets, pbasis, "MAO")
            END IF
         END DO
      CASE (mao_basis_ext)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="MAO")
            IF (ASSOCIATED(basis_set)) THEN
               basis_set%kind_radius = orb_basis_set_list(ikind)%gto_basis_set%kind_radius
               mao_basis_set_list(ikind)%gto_basis_set => basis_set
            END IF
         END DO
      CASE DEFAULT
         CPABORT("Unknown option for MAO basis")
      END SELECT
      IF (unit_nr > 0) THEN
         DO ikind = 1, nkind
            IF (.NOT. ASSOCIATED(mao_basis_set_list(ikind)%gto_basis_set)) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,I4)") &
                  "WARNING: No MAO basis set associated with Kind ", ikind
            ELSE
               nbas = mao_basis_set_list(ikind)%gto_basis_set%nsgf
               WRITE (UNIT=unit_nr, FMT="(T2,A,I4,T56,A,I10)") &
                  "MAO basis set Kind ", ikind, " Number of BSF:", nbas
            END IF
         END DO
      END IF

      IF (PRESENT(print_basis)) THEN
         IF (print_basis) THEN
            DO ikind = 1, nkind
               basis_set => mao_basis_set_list(ikind)%gto_basis_set
               CALL write_gto_basis_set(basis_set, unit_nr, "MAO REFERENCE BASIS")
            END DO
         END IF
      END IF

   END SUBROUTINE mao_reference_basis

! **************************************************************************************************
!> \brief Analyze the MAO basis, projection on angular functions
!> \param mao_coef ...
!> \param matrix_smm ...
!> \param mao_basis_set_list ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param unit_nr ...
!> \param para_env ...
!> \par History
!>      07.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE mao_basis_analysis(mao_coef, matrix_smm, mao_basis_set_list, particle_set, &
                                 qs_kind_set, unit_nr, para_env)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coef, matrix_smm
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: mao_basis_set_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_basis_analysis', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      INTEGER                                            :: ia, iab, iatom, ikind, iset, ishell, &
                                                            ispin, l, lmax, lshell, m, ma, na, &
                                                            natom, nspin
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cmask, vec1, vec2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weight
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block, cmao
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      ! Analyze the MAO basis
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(/,A)") " Analyze angular momentum character of MAOs "
         WRITE (unit_nr, "(T7,A,T15,A,T20,A,T40,A,T50,A,T60,A,T70,A,T80,A)") &
            "ATOM", "Spin", "MAO", "S", "P", "D", "F", "G"
      END IF
      lmax = 4 ! analyze up to g-functions
      natom = SIZE(particle_set)
      nspin = SIZE(mao_coef)
      DO iatom = 1, natom
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                              element_symbol=element_symbol, kind_number=ikind)
         basis_set => mao_basis_set_list(ikind)%gto_basis_set
         CALL get_qs_kind(qs_kind_set(ikind), mao=na)
         CALL get_gto_basis_set(basis_set, nsgf=ma)
         ALLOCATE (cmask(ma), vec1(ma), vec2(ma), weight(0:lmax, na))
         weight = 0.0_dp
         CALL cp_dbcsr_get_block_p(matrix=matrix_smm(1)%matrix, row=iatom, col=iatom, &
                                   block=block, found=found)
         DO ispin = 1, nspin
            CALL cp_dbcsr_get_block_p(matrix=mao_coef(ispin)%matrix, row=iatom, col=iatom, &
                                      block=cmao, found=found)
            IF (found) THEN
               DO l = 0, lmax
                  cmask = 0.0_dp
                  iab = 0
                  DO iset = 1, basis_set%nset
                     DO ishell = 1, basis_set%nshell(iset)
                        lshell = basis_set%l(ishell, iset)
                        DO m = -lshell, lshell
                           iab = iab+1
                           IF (l == lshell) cmask(iab) = 1.0_dp
                        END DO
                     END DO
                  END DO
                  DO ia = 1, na
                     vec1(1:ma) = cmask*cmao(1:ma, ia)
                     vec2(1:ma) = MATMUL(block, vec1)
                     weight(l, ia) = SUM(vec1(1:ma)*vec2(1:ma))
                  END DO
               END DO
            END IF
            CALL mp_sum(weight, para_env%group)
            IF (unit_nr > 0) THEN
               DO ia = 1, na
                  IF (ispin == 1 .AND. ia == 1) THEN
                     WRITE (unit_nr, "(i6,T9,A2,T17,i2,T20,i3,T31,5F10.4)") &
                        iatom, element_symbol, ispin, ia, weight(0:lmax, ia)
                  ELSE
                     WRITE (unit_nr, "(T17,i2,T20,i3,T31,5F10.4)") ispin, ia, weight(0:lmax, ia)
                  END IF
               END DO
            END IF
         END DO
         DEALLOCATE (cmask, weight, vec1, vec2)
      END DO
   END SUBROUTINE mao_basis_analysis

! **************************************************************************************************
!> \brief Calculte the Q=APA(T) matrix, A=(MAO,ORB) overlap
!> \param matrix_q ...
!> \param matrix_p ...
!> \param matrix_s ...
!> \param matrix_smm ...
!> \param matrix_smo ...
!> \param smm_list ...
!> \param electra ...
!> \param eps_filter ...
!> \param nimages ...
!> \param kpoints ...
!> \param matrix_ks ...
!> \param sab_orb ...
!> \par History
!>      08.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE mao_build_q(matrix_q, matrix_p, matrix_s, matrix_smm, matrix_smo, smm_list, &
                          electra, eps_filter, nimages, kpoints, matrix_ks, sab_orb)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_q
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_p, matrix_s
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_smm, matrix_smo
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: smm_list
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: electra
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN), OPTIONAL                      :: nimages
      TYPE(kpoint_type), OPTIONAL, POINTER               :: kpoints
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_ks
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_build_q', routineP = moduleN//':'//routineN

      INTEGER                                            :: im, ispin, nim, nocc, norb, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp)                                      :: elex, xkp(3)
      TYPE(cp_dbcsr_type)                                :: ksmat, pmat, smat, tmat

      nim = 1
      IF (PRESENT(nimages)) nim = nimages
      IF (nim > 1) THEN
         CPASSERT(PRESENT(kpoints))
         CPASSERT(PRESENT(matrix_ks))
         CPASSERT(PRESENT(sab_orb))
      END IF

      ! Reference
      nspin = SIZE(matrix_p, 1)
      DO ispin = 1, nspin
         electra(ispin) = 0.0_dp
         DO im = 1, nim
            CALL cp_dbcsr_trace(matrix_p(ispin, im)%matrix, matrix_s(1, im)%matrix, elex)
            electra(ispin) = electra(ispin)+elex
         END DO
      END DO

      ! Q matrix
      NULLIFY (matrix_q)
      CALL cp_dbcsr_allocate_matrix_set(matrix_q, nspin)
      DO ispin = 1, nspin
         ALLOCATE (matrix_q(ispin)%matrix)
         CALL cp_dbcsr_init(matrix_q(ispin)%matrix)
         CALL cp_dbcsr_create(matrix_q(ispin)%matrix, template=matrix_smm(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_q(ispin)%matrix, smm_list)
      END DO
      ! temp matrix
      CALL cp_dbcsr_init(tmat)
      CALL cp_dbcsr_create(tmat, template=matrix_smo(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      ! Q=APA(T)
      DO ispin = 1, nspin
         IF (nim == 1) THEN
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_smo(1)%matrix, matrix_p(ispin, 1)%matrix, &
                                   0.0_dp, tmat, filter_eps=eps_filter)
            CALL cp_dbcsr_multiply("N", "T", 1.0_dp, tmat, matrix_smo(1)%matrix, &
                                   0.0_dp, matrix_q(ispin)%matrix, filter_eps=eps_filter)
         ELSE
            ! k-points
            CALL cp_dbcsr_init(pmat)
            CALL cp_dbcsr_create(pmat, template=matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_init(smat)
            CALL cp_dbcsr_create(smat, template=matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_init(ksmat)
            CALL cp_dbcsr_create(ksmat, template=matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(pmat, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(smat, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(ksmat, sab_orb)
            NULLIFY (cell_to_index)
            CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
            ! calculate density matrix at gamma point
            xkp = 0.0_dp
            ! transform KS and S matrices to the gamma point
            CALL cp_dbcsr_set(ksmat, 0.0_dp)
            CALL rskp_transform(rmatrix=ksmat, rsmat=matrix_ks, ispin=ispin, &
                                xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_orb)
            CALL cp_dbcsr_set(smat, 0.0_dp)
            CALL rskp_transform(rmatrix=smat, rsmat=matrix_s, ispin=1, &
                                xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_orb)
            norb = NINT(electra(ispin))
            nocc = MOD(2, nspin)+1
            CALL calculate_p_gamma(pmat, ksmat, smat, kpoints, norb, REAL(nocc, KIND=dp))
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_smo(1)%matrix, pmat, &
                                   0.0_dp, tmat, filter_eps=eps_filter)
            CALL cp_dbcsr_multiply("N", "T", 1.0_dp, tmat, matrix_smo(1)%matrix, &
                                   0.0_dp, matrix_q(ispin)%matrix, filter_eps=eps_filter)
            CALL cp_dbcsr_release(pmat)
            CALL cp_dbcsr_release(smat)
            CALL cp_dbcsr_release(ksmat)
         END IF
      END DO
      ! free temp matrix
      CALL cp_dbcsr_release(tmat)

   END SUBROUTINE mao_build_q

END MODULE mao_methods
