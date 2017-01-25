!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate localized minimal basis
!> \par History
!>      12.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE minbas_methods
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_release,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_allocate_matrix_set, dbcsr_create, dbcsr_deallocate_matrix_set, &
        dbcsr_distribution_type, dbcsr_filter, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_reserve_diag_blocks, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE lapack,                          ONLY: lapack_ssyev
   USE mao_basis,                       ONLY: mao_generate_basis
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'minbas_methods'

   PUBLIC ::  minbas_calculation

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mos ...
!> \param quambo ...
!> \param mao ...
!> \param iounit ...
!> \param full_ortho ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE minbas_calculation(qs_env, mos, quambo, mao, iounit, full_ortho, eps_filter)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: quambo
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mao
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_ortho
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'minbas_calculation', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, i, iab, ispin, nao, natom, &
                                                            ndep, nmao, nmo, nmx, np, np1, nspin, &
                                                            nvirt, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: do_minbas, my_full_ortho
      REAL(KIND=dp)                                      :: my_eps_filter
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dval, dvalo, dvalv, eigval
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_a, fm_struct_b, fm_struct_c, &
                                                            fm_struct_d, fm_struct_e
      TYPE(cp_fm_type), POINTER                          :: fm1, fm2, fm3, fm4, fm5, fm6, fm_mos, &
                                                            fma, fmb, fmwork
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type)                                   :: smao, sortho
      TYPE(dbcsr_type), POINTER                          :: smat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      IF (PRESENT(iounit)) THEN
         unit_nr = iounit
      ELSE
         unit_nr = -1
      END IF

      IF (PRESENT(full_ortho)) THEN
         my_full_ortho = full_ortho
      ELSE
         my_full_ortho = .FALSE.
      END IF

      IF (PRESENT(eps_filter)) THEN
         my_eps_filter = eps_filter
      ELSE
         my_eps_filter = 1.0e-10_dp
      END IF

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspin = dft_control%nspins

      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, natom=natom)
      CALL get_ks_env(ks_env=ks_env, particle_set=particle_set, dbcsr_dist=dbcsr_dist)
      ALLOCATE (row_blk_sizes(natom), col_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes)
      CALL get_particle_set(particle_set, qs_kind_set, nmao=col_blk_sizes)
      nmao = SUM(col_blk_sizes)
      ! check if MAOs have been specified
      DO iab = 1, natom
         IF (col_blk_sizes(iab) < 0) &
            CPABORT("Number of MAOs has to be specified in KIND section for all elements")
      END DO
      CALL get_mo_set(mo_set=mos(1)%mo_set, nao=nao, nmo=nmo)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,T71,I10)') 'Total Number of Atomic Basis Set Functions   :', nao
         WRITE (unit_nr, '(T2,A,T71,I10)') 'Total Number of Minimal Basis Set Functions  :', nmao
         IF (nspin == 1) THEN
            WRITE (unit_nr, '(T2,A,T71,I10)') 'Total Number of Molecular Orbitals available :', nmo
         ELSE
            DO ispin = 1, nspin
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, nmo=nmx)
               WRITE (unit_nr, '(T2,A,i2,T71,I10)') &
                  'Total Number of Molecular Orbitals available for Spin ', ispin, nmx
            END DO
         END IF
      END IF
      CPASSERT(nmao <= nao)
      DO ispin = 1, nspin
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, nmo=nmx)
         IF (nmx /= nmo) EXIT
      END DO
      do_minbas = .TRUE.
      IF (nmao > nmo) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') 'Localized Minimal Basis Analysis not possible'
         END IF
         do_minbas = .FALSE.
      ELSEIF (nmo /= nmx) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') 'Different Number of Alpha and Beta MOs'
            WRITE (unit_nr, '(T2,A)') 'Localized Minimal Basis Analysis not possible'
         END IF
         do_minbas = .FALSE.
      ELSE
         IF (nao > nmo) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') 'WARNING: Only a subset of MOs is available: Analysis depends on MOs'
            END IF
         END IF
      END IF

      IF (do_minbas) THEN
         ! initialize QUAMBOs
         NULLIFY (quambo)
         CALL dbcsr_allocate_matrix_set(quambo, nspin)
         DO ispin = 1, nspin
            ! coeficients
            ALLOCATE (quambo(ispin)%matrix)
            CALL dbcsr_create(matrix=quambo(ispin)%matrix, &
                              name="QUAMBO", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         END DO

         ! initialize MAOs
         ! optimze MAOs (mao_coef is allocated in the routine)
         CALL mao_generate_basis(qs_env, mao_coef)

         ! sortho (nmao x nmao)
         CALL dbcsr_create(sortho, name="SORTHO", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL dbcsr_reserve_diag_blocks(matrix=sortho)

         DEALLOCATE (row_blk_sizes, col_blk_sizes)

         ! temporary FM matrices
         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         NULLIFY (fm_struct_a, fm_struct_b)
         CALL cp_fm_struct_create(fm_struct_a, nrow_global=nao, ncol_global=nmao, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_struct_create(fm_struct_b, nrow_global=nmo, ncol_global=nmao, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(fm1, fm_struct_a)
         CALL cp_fm_create(fm2, fm_struct_b)
         CALL cp_fm_create(fma, fm_struct_b)
         CALL cp_fm_create(fmb, fm_struct_b)

         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
         smat => matrix_s(1, 1)%matrix
         DO ispin = 1, nspin

            ! SMAO = Overlap*MAOs
            CALL dbcsr_create(smao, name="S*MAO", template=mao_coef(1)%matrix)
            CALL dbcsr_multiply("N", "N", 1.0_dp, smat, mao_coef(ispin)%matrix, 0.0_dp, smao)
            ! a(nj)* = C(vn)(T) * SMAO(vj)
            CALL copy_dbcsr_to_fm(smao, fm1)
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=fm_mos)
            CALL cp_gemm("T", "N", nmo, nmao, nao, 1.0_dp, fm_mos, fm1, 0.0_dp, fm2)
            CALL dbcsr_release(smao)
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,T51,A,i2,T71,I10)') 'MOs in Occupied Valence Set', 'Spin ', ispin, homo
            END IF
            nvirt = nmo-homo
            NULLIFY (fm_struct_c)
            CALL cp_fm_struct_create(fm_struct_c, nrow_global=nvirt, ncol_global=nvirt, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(fm3, fm_struct_c)
            CALL cp_fm_create(fm4, fm_struct_c)
            ! B(vw) = a(vj)* * a(wj)*
            CALL cp_gemm("N", "T", nvirt, nvirt, nmao, 1.0_dp, fm2, fm2, 0.0_dp, fm3, &
                         a_first_row=homo+1, b_first_row=homo+1)
            ALLOCATE (eigval(nvirt))
            CALL choose_eigv_solver(fm3, fm4, eigval)
            ! SVD(B) -> select p largest eigenvalues and vectors
            np = nmao-homo
            np1 = nvirt-np+1
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,T51,A,i2,T71,I10)') 'MOs in Virtual Valence Set', 'Spin ', ispin, np
            END IF
            ! R(vw) = SUM_p T(vp)*T(wp)
            CALL cp_gemm("N", "T", nvirt, nvirt, np, 1.0_dp, fm4, fm4, 0.0_dp, fm3, &
                         a_first_col=np1, b_first_col=np1)
            !
            ALLOCATE (dval(nmao), dvalo(nmao), dvalv(nmao))
            NULLIFY (fm_struct_d)
            CALL cp_fm_struct_create(fm_struct_d, nrow_global=nvirt, ncol_global=nmao, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(fm5, fm_struct_d)
            NULLIFY (fm_struct_e)
            CALL cp_fm_struct_create(fm_struct_e, nrow_global=nmao, ncol_global=nmao, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(fm6, fm_struct_e)
            ! D(j) = SUM_n (a(nj)*)^2 + SUM_vw R(vw) * a(vj)* * a(wj)*
            CALL cp_gemm("N", "N", nvirt, nmao, nvirt, 1.0_dp, fm3, fm2, 0.0_dp, fm5, &
                         b_first_row=homo+1)
            CALL cp_gemm("T", "N", nmao, nmao, nvirt, 1.0_dp, fm2, fm5, 0.0_dp, fm6, &
                         a_first_row=homo+1)
            CALL cp_fm_get_diag(fm6, dvalv(1:nmao))
            CALL cp_gemm("T", "N", nmao, nmao, homo, 1.0_dp, fm2, fm2, 0.0_dp, fm6)
            CALL cp_fm_get_diag(fm6, dvalo(1:nmao))
            DO i = 1, nmao
               dval(i) = 1.0_dp/SQRT(dvalo(i)+dvalv(i))
            END DO
            ! scale intermediate expansion
            CALL cp_fm_to_fm_submat(fm2, fma, homo, nmao, 1, 1, 1, 1)
            CALL cp_fm_to_fm_submat(fm5, fma, nvirt, nmao, 1, 1, homo+1, 1)
            CALL cp_fm_column_scale(fma, dval)
            ! Orthogonalization
            CALL cp_fm_create(fmwork, fm_struct_e)
            CALL cp_gemm("T", "N", nmao, nmao, nmo, 1.0_dp, fma, fma, 0.0_dp, fm6)
            IF (my_full_ortho) THEN
               ! full orthogonalization
               CALL cp_fm_power(fm6, fmwork, -0.5_dp, 1.0e-12_dp, ndep)
               IF (ndep > 0 .AND. unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A,T71,I10)') 'Warning: linear dependent basis   ', ndep
               END IF
               CALL cp_gemm("N", "N", nmo, nmao, nmao, 1.0_dp, fma, fm6, 0.0_dp, fmb)
            ELSE
               ! orthogonalize on-atom blocks
               CALL copy_fm_to_dbcsr(fm6, sortho, keep_sparsity=.TRUE.)
               CALL diag_sqrt_invert(sortho)
               CALL copy_dbcsr_to_fm(sortho, fm6)
               CALL cp_gemm("N", "N", nmo, nmao, nmao, 1.0_dp, fma, fm6, 0.0_dp, fmb)
            END IF
            ! store as QUAMBO
            CALL cp_gemm("N", "N", nao, nmao, nmo, 1.0_dp, fm_mos, fmb, 0.0_dp, fm1)
            CALL copy_fm_to_dbcsr(fm1, quambo(ispin)%matrix)
            CALL dbcsr_filter(quambo(ispin)%matrix, my_eps_filter)
            !
            DEALLOCATE (eigval, dval, dvalo, dvalv)
            CALL cp_fm_release(fm3)
            CALL cp_fm_release(fm4)
            CALL cp_fm_release(fm5)
            CALL cp_fm_release(fm6)
            CALL cp_fm_release(fmwork)
            CALL cp_fm_struct_release(fm_struct_c)
            CALL cp_fm_struct_release(fm_struct_d)
            CALL cp_fm_struct_release(fm_struct_e)

         END DO

         ! clean up
         CALL cp_fm_release(fm1)
         CALL cp_fm_release(fm2)
         CALL cp_fm_release(fma)
         CALL cp_fm_release(fmb)
         CALL cp_fm_struct_release(fm_struct_a)
         CALL cp_fm_struct_release(fm_struct_b)
         CALL dbcsr_release(sortho)

         ! return MAOs if requested
         IF (PRESENT(mao)) THEN
            mao => mao_coef
         ELSE
            CALL dbcsr_deallocate_matrix_set(mao_coef)
         END IF

      ELSE
         NULLIFY (quambo)
      END IF

      CALL timestop(handle)

   END SUBROUTINE minbas_calculation

! **************************************************************************************************
!> \brief ...
!> \param sortho ...
! **************************************************************************************************
   SUBROUTINE diag_sqrt_invert(sortho)
      TYPE(dbcsr_type)                                   :: sortho

      CHARACTER(len=*), PARAMETER :: routineN = 'diag_sqrt_invert', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, info, jatom, lwork, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, bmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sblock
      TYPE(dbcsr_iterator_type)                          :: dbcsr_iter

      CALL dbcsr_iterator_start(dbcsr_iter, sortho)
      DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, sblock)
         CPASSERT(iatom == jatom)
         n = SIZE(sblock, 1)
         lwork = MAX(n*n, 100)
         ALLOCATE (amat(n, n), bmat(n, n), w(n), work(lwork))
         amat(1:n, 1:n) = sblock(1:n, 1:n)
         info = 0
         CALL lapack_ssyev("V", "U", n, amat, n, w, work, lwork, info)
         CPASSERT(info == 0)
         w(1:n) = 1._dp/SQRT(w(1:n))
         DO i = 1, n
            bmat(1:n, i) = amat(1:n, i)*w(i)
         END DO
         sblock(1:n, 1:n) = MATMUL(amat, TRANSPOSE(bmat))
         DEALLOCATE (amat, bmat, w, work)
      END DO
      CALL dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE diag_sqrt_invert

END MODULE minbas_methods
