!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief orbital transformations
!> \par History
!>      Added Taylor expansion based computation of the matrix functions (01.2004)
!>      added additional rotation variables for non-equivalent occupied orbs (08.2004)
!> \author Joost VandeVondele (06.2002)
! *****************************************************************************
MODULE qs_ot
  USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                             cp_dbcsr_cholesky_invert,&
                                             cp_dbcsr_cholesky_restore
  USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_heevd,&
                                             cp_dbcsr_syevd
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, &
       cp_dbcsr_distribution, cp_dbcsr_filter, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_gershgorin_norm, cp_dbcsr_get_block_p, cp_dbcsr_get_info, &
       cp_dbcsr_get_occupation, cp_dbcsr_hadamard_product, cp_dbcsr_init, &
       cp_dbcsr_init_p, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_release, &
       cp_dbcsr_release_p, cp_dbcsr_scale, cp_dbcsr_scale_by_vector, &
       cp_dbcsr_set, cp_dbcsr_transposed, cp_dbcsr_type
  USE dbcsr_api,                       ONLY: dbcsr_distribution_mp,&
                                             dbcsr_mp_group
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE preconditioner,                  ONLY: apply_preconditioner
  USE preconditioner_types,            ONLY: preconditioner_type
  USE qs_ot_types,                     ONLY: qs_ot_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE 

  PUBLIC  :: qs_ot_get_p
  PUBLIC  :: qs_ot_get_orbitals
  PUBLIC  :: qs_ot_get_derivative
  PUBLIC  :: qs_ot_get_orbitals_ref
  PUBLIC  :: qs_ot_get_derivative_ref
  PUBLIC  :: qs_ot_new_preconditioner
  PUBLIC  :: qs_ot_get_scp_nddo_derivative
  PUBLIC  :: qs_ot_get_scp_nddo_coeffs
  PRIVATE :: qs_ot_p2m_diag
  PRIVATE :: qs_ot_sinc
  PRIVATE :: qs_ot_ref_poly
  PRIVATE :: qs_ot_ref_chol
  PRIVATE :: qs_ot_ref_lwdn
  PRIVATE :: qs_ot_ref_decide
  PRIVATE :: qs_ot_ref_update
  PRIVATE :: qs_ot_refine
  PRIVATE :: qs_ot_on_the_fly_localize

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot'

CONTAINS

  ! gets ready to use the preconditioner/ or renew the preconditioner
  ! only keeps a pointer to the preconditioner.
  ! If you change the preconditioner, you have to call this routine
  ! you remain responsible of proper deallocate of your preconditioner
  ! (or you can reuse it on the next step of the computation)
! *****************************************************************************
  SUBROUTINE qs_ot_new_preconditioner(qs_ot_env,preconditioner,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(preconditioner_type), POINTER       :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_new_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, ncoef
    LOGICAL                                  :: failure, mixed_precision

    failure = .FALSE.

    qs_ot_env%preconditioner => preconditioner
    qs_ot_env%os_valid = .FALSE.
    IF (.NOT. ASSOCIATED(qs_ot_env%matrix_psc0)) THEN
       CALL cp_dbcsr_init_p(qs_ot_env%matrix_psc0, error=error)
       CALL cp_dbcsr_copy(qs_ot_env%matrix_psc0,qs_ot_env%matrix_sc0,'matrix_psc0',error=error)
    ENDIF

    mixed_precision = qs_ot_env%settings%mixed_precision

    IF (.NOT. qs_ot_env%use_dx) THEN
       qs_ot_env%use_dx=.TRUE.
       CALL cp_dbcsr_init_p(qs_ot_env%matrix_dx, error=error)
       CALL cp_dbcsr_copy(qs_ot_env%matrix_dx,qs_ot_env%matrix_gx,'matrix_dx',error=error)
       IF (qs_ot_env%settings%do_rotation) THEN
          CALL cp_dbcsr_init_p(qs_ot_env%rot_mat_dx, error=error)
          CALL cp_dbcsr_copy(qs_ot_env%rot_mat_dx,qs_ot_env%rot_mat_gx,'rot_mat_dx',error=error)
       ENDIF
       IF (qs_ot_env%settings%do_ener) THEN
          ncoef = SIZE ( qs_ot_env % ener_gx)
          ALLOCATE ( qs_ot_env%ener_dx ( ncoef ), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          qs_ot_env%ener_dx = 0.0_dp
       ENDIF
       ! ***SCP
       IF ( qs_ot_env % settings % scp_nddo ) THEN
          ALLOCATE(qs_ot_env % dxmat)
          CALL cp_dbcsr_init(qs_ot_env % dxmat, error=error)
          CALL cp_dbcsr_copy(qs_ot_env % dxmat, qs_ot_env % gxmat, "SCP_DXMAT", error=error)
          CALL cp_dbcsr_set (qs_ot_env % dxmat, 0.0_dp, error=error )
       ENDIF
       ! ***SCP
    ENDIF

  END SUBROUTINE qs_ot_new_preconditioner

! *****************************************************************************
  SUBROUTINE qs_ot_on_the_fly_localize(qs_ot_env, C_NEW, SC, G_OLD, D, error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_dbcsr_type), POINTER             :: C_NEW, SC, G_OLD, D
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_on_the_fly_localize', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: taylor_order = 50
    REAL(KIND=dp), PARAMETER                 :: alpha = 0.1_dp, &
                                                f2_eps = 0.01_dp, &
                                                rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: blk, col, col_size, handle, &
                                                i, k, n, output_unit, p, row, &
                                                row_size
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    REAL(KIND=dp)                            :: expfactor, f2, norm_fro, &
                                                norm_gct, tmp
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: C, Gp1, Gp2, GU, U
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    !
    !
    CALL cp_dbcsr_get_info(C_NEW,nfullrows_total=n,nfullcols_total=k)
    !
    ! C = C*expm(-G)
    GU => qs_ot_env%buf1_k_k_nosym ! a buffer
    U  => qs_ot_env%buf2_k_k_nosym ! a buffer
    Gp1=> qs_ot_env%buf3_k_k_nosym ! a buffer
    Gp2=> qs_ot_env%buf4_k_k_nosym ! a buffer
    C  => qs_ot_env%buf1_n_k       ! a buffer
    !
    ! compute the derivative of the norm
    !-------------------------------------------------------------------
    ! (x^2+eps)^1/2
    f2 = 0.0_dp
    CALL cp_dbcsr_copy(C,C_NEW,error=error)
    CALL cp_dbcsr_iterator_start(iter, C)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, DATA, blk,&
            row_size=row_size, col_size=col_size)
       DO p=1,col_size! p
       DO i=1,row_size! i
          tmp = SQRT( DATA(i,p)**2 + f2_eps )
          f2 = f2 + tmp
          DATA(i,p) = DATA(i,p) / tmp
       ENDDO
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
    CALL mp_sum(f2,dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(C))))
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)
    IF(output_unit>0) WRITE(output_unit,*) routineN//' f2 =',f2
    !
    !
    CALL cp_dbcsr_multiply('T','N',1.0_dp,C,C_NEW,0.0_dp,GU,error=error)
    !
    ! antisymetrize
    CALL cp_dbcsr_transposed (U, GU, shallow_data_copy=.FALSE., &
         use_distribution=cp_dbcsr_distribution(GU), &
         transpose_distribution=.FALSE., &
         error=error)
    CALL cp_dbcsr_add(GU,U,alpha_scalar=-0.5_dp,beta_scalar=0.5_dp,error=error)
    !-------------------------------------------------------------------
    !
    norm_fro = cp_dbcsr_frobenius_norm(GU)
    norm_gct = cp_dbcsr_gershgorin_norm(GU)
    !write(*,*) 'qs_ot_localize: ||P-I||_f=',norm_fro,' ||P-I||_GCT=',norm_gct
    !
    !kscale = CEILING(LOG(MIN(norm_fro,norm_gct))/LOG(2.0_dp))
    !scale  = LOG(MIN(norm_fro,norm_gct))/LOG(2.0_dp)
    !write(*,*) 'qs_ot_localize: scale=',scale,' kscale=',kscale
    !
    ! rescale for steepest descent
    CALL cp_dbcsr_scale(GU, -alpha, error=error)
    !
    ! compute unitary transform
    ! zeroth and first order
    expfactor = 1.0_dp
    CALL cp_dbcsr_copy(U,GU,error=error)
    CALL cp_dbcsr_scale(U,expfactor,error=error)
    CALL cp_dbcsr_add_on_diag(U,1.0_dp,error=error)
    ! other orders
    CALL cp_dbcsr_copy(Gp1,GU,error=error)
    DO i = 2,taylor_order
       ! new power of G
       CALL cp_dbcsr_multiply('N','N',1.0_dp,GU,Gp1,0.0_dp,Gp2,error=error)
       CALL cp_dbcsr_copy(Gp1,Gp2,error=error)
       ! add to the taylor expansion so far
       expfactor = expfactor / REAL(i,KIND=dp)
       CALL cp_dbcsr_add(U,Gp1,alpha_scalar=1.0_dp,beta_scalar=expfactor,error=error)
       norm_fro = cp_dbcsr_frobenius_norm(Gp1)
       !write(*,*) 'Taylor expansion i=',i,' norm(X^i)/i!=',norm_fro*expfactor
       IF(norm_fro*expfactor.LT.1.0E-10_dp) EXIT
    ENDDO
    !
    ! rotate MOs
    CALL cp_dbcsr_multiply('N','N',1.0_dp,C_NEW,U,0.0_dp,C,error=error)
    CALL cp_dbcsr_copy(C_NEW,C,error=error)
    !
    ! rotate SC
    CALL cp_dbcsr_multiply('N','N',1.0_dp,SC,U,0.0_dp,C,error=error)
    CALL cp_dbcsr_copy(SC,C,error=error)
    !
    ! rotate D_i
    CALL cp_dbcsr_multiply('N','N',1.0_dp,D,U,0.0_dp,C,error=error)
    CALL cp_dbcsr_copy(D,C,error=error)
    !
    ! rotate G_i-1
    IF(ASSOCIATED(G_OLD)) THEN
       CALL cp_dbcsr_multiply('N','N',1.0_dp,G_OLD,U,0.0_dp,C,error=error)
       CALL cp_dbcsr_copy(G_OLD,C,error=error)
    ENDIF
    !
    CALL timestop(handle)
  END SUBROUTINE qs_ot_on_the_fly_localize

! *****************************************************************************
  SUBROUTINE qs_ot_ref_chol(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update, error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_dbcsr_type)                      :: C_OLD, C_TMP, C_NEW, P, SC
    LOGICAL, INTENT(IN)                      :: update
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_ref_chol', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, k, n

    CALL timeset(routineN,handle)
    !
    CALL cp_dbcsr_get_info(C_NEW,nfullrows_total=n,nfullcols_total=k)
    !
    ! P = U'*U
    CALL cp_dbcsr_cholesky_decompose(P,k,qs_ot_env%para_env,qs_ot_env%blacs_env,error=error)
    !
    ! C_NEW = C_OLD*inv(U)
    CALL cp_dbcsr_cholesky_restore(C_OLD,k,P,C_NEW,op="SOLVE",pos="RIGHT",&
         transa="N",para_env=qs_ot_env%para_env,blacs_env=qs_ot_env%blacs_env,&
         error=error)
    !
    ! Update SC if needed
    IF(update) THEN
       CALL cp_dbcsr_cholesky_restore(SC,k,P,C_TMP,op="SOLVE",pos="RIGHT",&
            transa="N",para_env=qs_ot_env%para_env,blacs_env=qs_ot_env%blacs_env,error=error)
       CALL cp_dbcsr_copy(SC,C_TMP,error=error)
    ENDIF
    !
    CALL timestop(handle)
  END SUBROUTINE qs_ot_ref_chol

! *****************************************************************************
  SUBROUTINE qs_ot_ref_lwdn(qs_ot_env, C_OLD, C_TMP, C_NEW, P, SC, update, error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_dbcsr_type)                      :: C_OLD, C_TMP, C_NEW, P, SC
    LOGICAL, INTENT(IN)                      :: update
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_ref_lwdn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, istat, k, n
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eig, fun
    TYPE(cp_dbcsr_type), POINTER             :: V, W

    failure = .FALSE.

    CALL timeset(routineN,handle)
    !
    CALL cp_dbcsr_get_info(C_NEW,nfullrows_total=n,nfullcols_total=k)
    !
    V   => qs_ot_env%buf1_k_k_nosym ! a buffer
    W   => qs_ot_env%buf2_k_k_nosym ! a buffer
    ALLOCATE(eig(k), fun(k), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL cp_dbcsr_syevd(P,V,eig,qs_ot_env%para_env,qs_ot_env%blacs_env,error=error)
    !
    ! compute the P^(-1/2)
    DO i = 1,k
       IF(eig(i).LE.0.0_dp) &
            & CALL stop_program(routineN,moduleN,__LINE__,"P not positive definite")
       IF(eig(i).LT.1.0E-8_dp) THEN
          fun(i)=0.0_dp
       ELSE
          fun(i)=1.0_dp/SQRT(eig(i))
       ENDIF
    ENDDO
    CALL cp_dbcsr_copy(W,V,error=error)
    CALL cp_dbcsr_scale_by_vector(V,alpha=fun,side='right',error=error)
    CALL cp_dbcsr_multiply('N','T',1.0_dp,W,V,0.0_dp,P,error=error)
    !
    ! Update C
    CALL cp_dbcsr_multiply('N','N',1.0_dp,C_OLD,P,0.0_dp,C_NEW,error=error)
    !
    ! Update SC if needed
    IF(update) THEN
       CALL cp_dbcsr_multiply('N','N',1.0_dp,SC,P,0.0_dp,C_TMP,error=error)
       CALL cp_dbcsr_copy(SC,C_TMP,error=error)
    ENDIF
    !
    DEALLOCATE(eig, fun, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    !
    CALL timestop(handle)
  END SUBROUTINE qs_ot_ref_lwdn

! *****************************************************************************
  SUBROUTINE qs_ot_ref_poly(qs_ot_env,C_OLD,C_TMP,C_NEW,P,SC,norm_in,update,output_unit,error)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_dbcsr_type), POINTER             :: C_OLD, C_TMP, C_NEW, P
    TYPE(cp_dbcsr_type)                      :: SC
    REAL(dp), INTENT(IN)                     :: norm_in
    LOGICAL, INTENT(IN)                      :: update
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_ref_poly', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, irefine, k, n
    LOGICAL                                  :: quick_exit
    REAL(dp)                                 :: norm, norm_fro, norm_gct, &
                                                occ_in, occ_out, rescale
    TYPE(cp_dbcsr_type), POINTER             :: BUF1, BUF2, BUF_NOSYM, FT, FY

    CALL timeset(routineN,handle)
    !
    CALL cp_dbcsr_get_info(C_NEW,nfullrows_total=n,nfullcols_total=k)
    !
    BUF_NOSYM => qs_ot_env%buf1_k_k_nosym! a buffer
    BUF1      => qs_ot_env%buf1_k_k_sym  ! a buffer
    BUF2      => qs_ot_env%buf2_k_k_sym  ! a buffer
    FY        => qs_ot_env%buf3_k_k_sym  ! a buffer
    FT        => qs_ot_env%buf4_k_k_sym  ! a buffer
    !
    ! initialize the norm (already computed in qs_ot_get_orbitals_ref)
    norm = norm_in
    !
    ! can we do a quick exit?
    quick_exit = .FALSE.
    IF(norm.LT.qs_ot_env%settings%eps_irac_quick_exit) quick_exit = .TRUE.
    !
    ! lets refine
    rescale = 1.0_dp
    DO irefine = 1,qs_ot_env%settings%max_irac
       !
       ! rescaling
       IF(norm.GT.1.0_dp) THEN
          IF(output_unit>0) WRITE(output_unit,'(A,I3,A)') &
               & routineN,irefine,': we rescale (C+a*D)'
          CALL cp_dbcsr_scale(P,1.0_dp/norm,error=error)
          rescale = rescale/SQRT(norm)
       ENDIF
       !
       ! get the refinement polynomial
       CALL qs_ot_refine(P, FY, BUF1, BUF2, qs_ot_env%settings%irac_degree, &
            qs_ot_env%settings%eps_irac_filter_matrix, output_unit, error)
       !
       ! collect the transformation
       IF(irefine.EQ.1) THEN
         CALL cp_dbcsr_copy(FT, FY, name='FT', error=error)
       ELSE
          CALL cp_dbcsr_multiply('N', 'N', 1.0_dp, FT, FY, 0.0_dp, BUF1, error=error)
          IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
             occ_in = cp_dbcsr_get_occupation(buf1)
             CALL cp_dbcsr_filter(buf1,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
             occ_out = cp_dbcsr_get_occupation(buf1)
             IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
                  WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(BUF1): occ_in',occ_in,' occ_out',occ_out
          ENDIF
          CALL cp_dbcsr_copy(FT, BUF1, name='FT', error=error)
       ENDIF
       !
       ! quick exit if possible
       IF(quick_exit) THEN
          IF(output_unit>0) WRITE(output_unit,'(A,I3,A)') &
               & routineN,irefine,': quick exit!'
          EXIT
       ENDIF
       !
       ! P = FY^T * P * FY
       CALL cp_dbcsr_multiply('N', 'N', 1.0_dp, P, FY, 0.0_dp, BUF_NOSYM, error=error)
       IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(buf_nosym)
          CALL cp_dbcsr_filter(buf_nosym,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(buf_nosym)
          IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(BUF_NOSYM): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       CALL cp_dbcsr_multiply('N', 'N', 1.0_dp, FY, BUF_NOSYM, 0.0_dp, P, error=error)
       IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(p)
          CALL cp_dbcsr_filter(p,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(p)
          IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(P): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       !
       ! check ||P-1||_gct
       CALL cp_dbcsr_add_on_diag(P, -1.0_dp, error=error)
       norm_fro = cp_dbcsr_frobenius_norm(P)
       norm_gct = cp_dbcsr_gershgorin_norm(P)
       CALL cp_dbcsr_add_on_diag(P, 1.0_dp, error=error)
       norm = MIN(norm_gct,norm_fro)
       !
       ! printing
       IF(output_unit>0) WRITE(output_unit,'(A,I3,A,E12.5)') &
            & routineN,irefine,': ||P-I||=',norm
       !
       ! blows up
       IF (norm > 1.0E10_dp) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Refinement blows up! "//&
                            "We need you to improve the code, please post your input on "//&
                            "the forum http://www.cp2k.org/")
       END IF
       !
       ! can we do a quick exit next step?
       IF(norm.LT.qs_ot_env%settings%eps_irac_quick_exit) quick_exit = .TRUE.
       !
       ! are we done?
       IF(norm.LT.qs_ot_env%settings%eps_irac) EXIT
       !
    ENDDO
    !
    ! C_NEW = C_NEW * FT * rescale
    CALL cp_dbcsr_multiply('N', 'N', rescale, C_OLD, FT, 0.0_dp, C_NEW, error=error)
    IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
       occ_in = cp_dbcsr_get_occupation(c_new)
       CALL cp_dbcsr_filter(c_new,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
       occ_out = cp_dbcsr_get_occupation(c_new)
       IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
            WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(C_NEW): occ_in',occ_in,' occ_out',occ_out
    ENDIF
    !
    ! update SC = SC * FY * rescale
    IF(update) THEN
       CALL cp_dbcsr_multiply('N', 'N', rescale, SC, FT, 0.0_dp, C_TMP, error=error)
       IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(c_tmp)
          CALL cp_dbcsr_filter(c_tmp,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(c_tmp)
          IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(C_TMP): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       CALL cp_dbcsr_copy(SC, C_TMP, error=error)
    ENDIF
    !
    CALL timestop(handle)
  END SUBROUTINE qs_ot_ref_poly

! *****************************************************************************
  FUNCTION qs_ot_ref_update(qs_ot_env1) RESULT(update)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env1
    LOGICAL                                  :: update

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_ref_update', &
      routineP = moduleN//':'//routineN

    update = .FALSE.
    SELECT CASE(qs_ot_env1%settings%ot_method)
    CASE("CG")
       SELECT CASE(qs_ot_env1%settings%line_search_method)
       CASE("2PNT")
          IF(qs_ot_env1%line_search_count.EQ.2) update = .TRUE.
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,"NYI")
       END SELECT
    CASE("DIIS")
       update = .TRUE.
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"NYI")
    END SELECT
  END FUNCTION qs_ot_ref_update

! *****************************************************************************
  SUBROUTINE qs_ot_ref_decide(qs_ot_env1, norm_in, ortho_irac)
    !
    TYPE(qs_ot_type)                         :: qs_ot_env1
    REAL(dp), INTENT(IN)                     :: norm_in
    CHARACTER(LEN=*), INTENT(INOUT)          :: ortho_irac

    ortho_irac = qs_ot_env1%settings%ortho_irac
    IF(norm_in.LT.qs_ot_env1%settings%eps_irac_switch) ortho_irac = "POLY"
  END SUBROUTINE qs_ot_ref_decide

! *****************************************************************************
  SUBROUTINE qs_ot_get_orbitals_ref(matrix_c, matrix_s, matrix_x, matrix_sx, &
       &                            matrix_gx_old, matrix_dx, qs_ot_env, &
       &                            qs_ot_env1, output_unit, error)
    !
    TYPE(cp_dbcsr_type), POINTER             :: matrix_c, matrix_s, matrix_x, &
                                                matrix_sx, matrix_gx_old, &
                                                matrix_dx
    TYPE(qs_ot_type)                         :: qs_ot_env, qs_ot_env1
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals_ref', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    CHARACTER(LEN=4)                         :: ortho_irac
    INTEGER                                  :: handle, k, n
    LOGICAL                                  :: on_the_fly_loc, update
    REAL(dp)                                 :: norm, norm_fro, norm_gct, &
                                                occ_in, occ_out
    TYPE(cp_dbcsr_type), POINTER             :: C_NEW, C_OLD, C_TMP, D, &
                                                G_OLD, P, S, SC

    CALL timeset(routineN,handle)

    IF(output_unit>0.AND.qs_ot_env1%settings%eps_irac_filter_matrix.GT.0.0_dp) &
         WRITE(output_unit,*) routinen//' eps_irac_filter_matrix=',&
         qs_ot_env1%settings%eps_irac_filter_matrix

    CALL cp_dbcsr_get_info(matrix_c,nfullrows_total=n,nfullcols_total=k)
    !
    C_NEW => matrix_c
    C_OLD => matrix_x ! need to be carefully updated for the gradient !
    SC    => matrix_sx! need to be carefully updated for the gradient !
    G_OLD => matrix_gx_old ! need to be carefully updated for localization !
    D     => matrix_dx     ! need to be carefully updated for localization !
    S     => matrix_s

    P     => qs_ot_env%p_k_k_sym ! a buffer
    C_TMP => qs_ot_env%buf1_n_k  ! a buffer
    !
    ! do we need to update C_OLD and SC?
    update = qs_ot_ref_update(qs_ot_env1)
    !
    ! do we want to on the fly localize?
    ! for the moment this is set from the input,
    ! later we might want to localize every n-step or
    ! when the sparsity increases...
    on_the_fly_loc = qs_ot_env1%settings%on_the_fly_loc
    !
    ! compute SC = S*C
    IF(ASSOCIATED(S)) THEN
       CALL cp_dbcsr_multiply('N','N',1.0_dp,S,C_OLD,0.0_dp,SC,error=error)
       IF (qs_ot_env1%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(sc)
          CALL cp_dbcsr_filter(sc,qs_ot_env1%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(sc)
          IF(output_unit>0.AND.qs_ot_env1%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(SC): occ_in',occ_in,' occ_out',occ_out
       ENDIF
    ELSE
       CALL cp_dbcsr_copy(SC,C_OLD,error=error)
    ENDIF
    !
    ! compute P = C'*SC
    CALL cp_dbcsr_multiply('T','N',1.0_dp,C_OLD,SC,0.0_dp,P,error=error)
    IF (qs_ot_env1%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
       occ_in = cp_dbcsr_get_occupation(p)
       CALL cp_dbcsr_filter(p,qs_ot_env1%settings%eps_irac_filter_matrix,error=error)
       occ_out = cp_dbcsr_get_occupation(p)
       IF(output_unit>0.AND.qs_ot_env1%settings%eps_irac_filter_matrix.GT.0.0_dp) &
            WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(P): occ_in',occ_in,' occ_out',occ_out
    ENDIF
    !
    ! check ||P-1||_f and ||P-1||_gct
    CALL cp_dbcsr_add_on_diag(P, -1.0_dp, error=error)
    norm_fro = cp_dbcsr_frobenius_norm(P)
    norm_gct = cp_dbcsr_gershgorin_norm(P)
    CALL cp_dbcsr_add_on_diag(P, 1.0_dp, error=error)
    norm = MIN(norm_gct,norm_fro)
    CALL qs_ot_ref_decide(qs_ot_env1,norm,ortho_irac)
    IF(output_unit>0) WRITE(output_unit,'(A,I3,A,E12.5,A)') &
         & routineN,0,': ||P-I||=',norm,&
         & ', ortho_irac = '//ortho_irac
    !
    ! select the orthogonality method
    SELECT CASE(ortho_irac)
    CASE("CHOL")
       CALL qs_ot_ref_chol(qs_ot_env,C_OLD,C_TMP,C_NEW,P,SC,update,error)
    CASE("LWDN")
       CALL qs_ot_ref_lwdn(qs_ot_env,C_OLD,C_TMP,C_NEW,P,SC,update,error)
    CASE("POLY")
       CALL qs_ot_ref_poly(qs_ot_env,C_OLD,C_TMP,C_NEW,P,SC,norm,update,output_unit,error)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Wrong argument")
    END SELECT
    !
    ! We update the C_i+1 and localization
    IF(update) THEN
       IF(on_the_fly_loc) THEN
          IF(output_unit>0) WRITE(output_unit,'(A)') &
               & routineN//' we localize C'
          CALL qs_ot_on_the_fly_localize(qs_ot_env,C_NEW,SC,G_OLD,D,error)
       ENDIF
       CALL cp_dbcsr_copy(C_OLD,C_NEW,error=error)
    ENDIF
    !
    CALL timestop(handle)
  END SUBROUTINE qs_ot_get_orbitals_ref

  SUBROUTINE qs_ot_refine(P,FY,P2,T,irac_degree,eps_irac_filter_matrix,output_unit,error)
    !----------------------------------------------------------------------
    ! refinement polynomial of degree 2,3 and 4 (PRB 70, 193102 (2004))
    !----------------------------------------------------------------------

    TYPE(cp_dbcsr_type), INTENT(inout)       :: P, FY, P2, T
    INTEGER, INTENT(in)                      :: irac_degree
    REAL(dp), INTENT(in)                     :: eps_irac_filter_matrix
    INTEGER, INTENT(in)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_refine', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, k
    REAL(dp)                                 :: occ_in, occ_out, r

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(P,nfullcols_total=k)
    SELECT CASE(irac_degree)
    CASE(2)
       ! C_out = C_in * ( 15/8 * I - 10/8 * P + 3/8 * P^2)
       r =   3.0_dp/8.0_dp
       CALL cp_dbcsr_multiply('N', 'N', r, P, P, 0.0_dp, FY, error=error)
       IF (eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(fy)
          CALL cp_dbcsr_filter(fy,eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(fy)
          IF(output_unit>0.AND.eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(FY): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       r = -10.0_dp/8.0_dp
       CALL cp_dbcsr_add(FY, P, alpha_scalar=1.0_dp, beta_scalar=r, error=error)
       r =  15.0_dp/8.0_dp
       CALL cp_dbcsr_add_on_diag(FY, alpha_scalar=r, error=error)
    CASE(3)
       ! C_out = C_in * ( 35/16 * I - 35/16 * P + 21/16 * P^2 - 5/16 P^3)
       CALL cp_dbcsr_multiply('N', 'N', 1.0_dp, P, P, 0.0_dp, P2, error=error)
       IF (eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(p2)
          CALL cp_dbcsr_filter(p2,eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(p2)
          IF(output_unit>0.AND.eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(P2): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       r =  -5.0_dp/16.0_dp
       CALL cp_dbcsr_multiply('N', 'N', r, P2, P, 0.0_dp, FY, error=error)
       IF (eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(fy)
          CALL cp_dbcsr_filter(fy,eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(fy)
          IF(output_unit>0.AND.eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(FY): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       r =  21.0_dp/16.0_dp
       CALL cp_dbcsr_add(FY, P2, alpha_scalar=1.0_dp, beta_scalar=r, error=error)
       r = -35.0_dp/16.0_dp
       CALL cp_dbcsr_add(FY, P, alpha_scalar=1.0_dp, beta_scalar=r, error=error)
       r =  35.0_dp/16.0_dp
       CALL cp_dbcsr_add_on_diag(FY, alpha_scalar=r, error=error)
    CASE(4)
       ! C_out = C_in * ( 315/128 * I - 420/128 * P + 378/128 * P^2 - 180/128 P^3 + 35/128 P^4 )
       !       = C_in * ( 315/128 * I - 420/128 * P + 378/128 * P^2 + ( - 180/128 * P + 35/128 * P^2 ) * P^2 )
       CALL cp_dbcsr_multiply('N', 'N', 1.0_dp, P, P, 0.0_dp, P2, error=error)   ! P^2
       IF (eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(p2)
          CALL cp_dbcsr_filter(p2,eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(p2)
          IF(output_unit>0.AND.eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(P2): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       r = -180.0_dp/128.0_dp
       CALL cp_dbcsr_add(T, P, alpha_scalar=0.0_dp, beta_scalar=r, error=error)  ! T=-180/128*P
       r =   35.0_dp/128.0_dp
       CALL cp_dbcsr_add(T, P2, alpha_scalar=1.0_dp, beta_scalar=r, error=error) ! T=T+35/128*P^2
       CALL cp_dbcsr_multiply('N', 'N', 1.0_dp, T, P2, 0.0_dp, FY, error=error)  ! Y=T*P^2
       IF (eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(fy)
          CALL cp_dbcsr_filter(fy,eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(fy)
          IF(output_unit>0.AND.eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(FY): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       r =  378.0_dp/128.0_dp
       CALL cp_dbcsr_add(FY, P2, alpha_scalar=1.0_dp, beta_scalar=r, error=error)! Y=Y+378/128*P^2
       r = -420.0_dp/128.0_dp
       CALL cp_dbcsr_add(FY, P, alpha_scalar=1.0_dp, beta_scalar=r, error=error) ! Y=Y-420/128*P
       r =  315.0_dp/128.0_dp
       CALL cp_dbcsr_add_on_diag(FY, alpha_scalar=r, error=error)                ! Y=Y+315/128*I
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"This irac_order NYI")
    END SELECT
    CALL timestop(handle)
  END SUBROUTINE qs_ot_refine


! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative_ref(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       &                              qs_ot_env,output_unit,error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_ref', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n
    LOGICAL                                  :: mixed_precision
    REAL(dp)                                 :: occ_in, occ_out
    TYPE(cp_dbcsr_type), POINTER             :: C, CHC, G, G_dp, HC, SC

    CALL timeset(routineN,handle)

    mixed_precision = qs_ot_env%settings%mixed_precision
    CALL cp_dbcsr_get_info(matrix_x,nfullrows_total=n,nfullcols_total=k)
    !
    C   => matrix_x              ! NBsf*NOcc
    SC  => matrix_sx             ! NBsf*NOcc need to be up2date
    HC  => matrix_hc             ! NBsf*NOcc
    G   => matrix_gx             ! NBsf*NOcc
    CHC  => qs_ot_env%buf1_k_k_sym ! buffer
    G_dp => qs_ot_env%buf1_n_k_dp  ! buffer dp

    IF(mixed_precision) THEN
       ! C'*(H*C)
       CALL cp_dbcsr_multiply('T','N',1.0_dp,C,HC,rzero,CHC,error=error)
       IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(chc)
          CALL cp_dbcsr_filter(chc,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(chc)
          IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(CHC): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       ! (S*C)*(C'*H*C)
       CALL cp_dbcsr_multiply('N','N',1.0_dp,SC,CHC,0.0_dp,G_dp,error=error)
       IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(g_dp)
          CALL cp_dbcsr_filter(g_dp,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(g_dp)
          IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(G_dp): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       ! G = 2*(1-S*C*C')*H*C
       CALL cp_dbcsr_add(G_dp,HC,alpha_scalar=-1.0_dp,beta_scalar=1.0_dp,error=error)
       CALL cp_dbcsr_copy(G,G_dp,error=error)
    ELSE
       ! C'*(H*C)
       CALL cp_dbcsr_multiply('T','N',1.0_dp,C,HC,0.0_dp,CHC,error=error)
       IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(chc)
          CALL cp_dbcsr_filter(chc,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(chc)
          IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(CHC): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       ! (S*C)*(C'*H*C)
       CALL cp_dbcsr_multiply('N','N',1.0_dp,SC,CHC,0.0_dp,G,error=error)
       IF (qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) THEN
          occ_in = cp_dbcsr_get_occupation(g)
          CALL cp_dbcsr_filter(g,qs_ot_env%settings%eps_irac_filter_matrix,error=error)
          occ_out = cp_dbcsr_get_occupation(g)
          IF(output_unit>0.AND.qs_ot_env%settings%eps_irac_filter_matrix.GT.0.0_dp) &
               WRITE(output_unit,'(2(A,F8.5))') routinen//' filter(G): occ_in',occ_in,' occ_out',occ_out
       ENDIF
       ! G = 2*(1-S*C*C')*H*C
       CALL cp_dbcsr_add(G,HC,alpha_scalar=-1.0_dp,beta_scalar=1.0_dp,error=error)
    ENDIF
    !
    CALL timestop(handle)
  END SUBROUTINE qs_ot_get_derivative_ref

  ! computes p=x*S*x and the matrix functionals related matrices
! *****************************************************************************
  SUBROUTINE qs_ot_get_p(matrix_x,matrix_sx,qs_ot_env,error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix_x, matrix_sx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_p', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(matrix_x,nfullrows_total=n,nfullcols_total=k)

    ! get the overlap
    CALL cp_dbcsr_multiply('T','N',rone,matrix_x,matrix_sx,rzero,&
         qs_ot_env%matrix_p,error=error)
    ! get an upper bound for the largest eigenvalue
    qs_ot_env % largest_eval_upper_bound = cp_dbcsr_gershgorin_norm(qs_ot_env%matrix_p)
    CALL decide_strategy(qs_ot_env)
    IF (qs_ot_env % do_taylor) THEN
       CALL qs_ot_p2m_taylor(qs_ot_env,error=error)
    ELSE
       CALL qs_ot_p2m_diag(qs_ot_env,error=error)
    ENDIF

    IF (qs_ot_env % settings % do_rotation) THEN
       CALL qs_ot_generate_rotation(qs_ot_env,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE qs_ot_get_p

! *****************************************************************************
!> \brief computes the rotation matrix rot_mat_u that is associated to a given
!>      rot_mat_x using rot_mat_u=exp(rot_mat_x)
!> \param a valid qs_ot_env
!> \par History
!>      08.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_ot_generate_rotation(qs_ot_env,error)

    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_generate_rotation', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: cone = (1.0_dp,0.0_dp), &
                                                czero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: evals_exp
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: data_z
    INTEGER                                  :: blk, col, handle, k, row
    LOGICAL                                  :: found
    REAL(KIND=dp), DIMENSION(:), POINTER     :: data_d
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: cmat_u, cmat_x

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(qs_ot_env%rot_mat_x,nfullrows_total=k)

    IF (k/=0) THEN
       CALL cp_dbcsr_init(cmat_x, error=error)
       CALL cp_dbcsr_init(cmat_u, error=error)
       CALL cp_dbcsr_copy(cmat_x,qs_ot_env%rot_mat_evec,name='cmat_x',error=error)
       CALL cp_dbcsr_copy(cmat_u,qs_ot_env%rot_mat_evec,name='cmat_u',error=error)
       ALLOCATE(evals_exp(k))

       ! rot_mat_u = exp(rot_mat_x)
       ! i rot_mat_x is hermitian, so go over the complex variables for diag
       !vwCALL cp_cfm_get_info(cmat_x,local_data=local_data_c,error=error)
       !vwCALL cp_fm_get_info(qs_ot_env%rot_mat_x,local_data=local_data_r,error=error)
       !vwlocal_data_c=CMPLX(0.0_dp,local_data_r,KIND=dp)
       CALL cp_dbcsr_iterator_start(iter, cmat_x)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, row, col, data_z, blk)
          CALL cp_dbcsr_get_block_p(qs_ot_env%rot_mat_x, row, col, data_d, found)
          IF(.NOT.found) THEN
             WRITE(*,*) routineN//' .NOT.found'
             !stop
          ELSE
             data_z=CMPLX(0.0_dp,data_d,KIND=dp)
          ENDIF
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)


       CALL cp_dbcsr_heevd(cmat_x,qs_ot_env%rot_mat_evec,qs_ot_env%rot_mat_evals,&
            qs_ot_env%para_env, qs_ot_env%blacs_env, error=error)
       evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * qs_ot_env%rot_mat_evals(:) )
       CALL cp_dbcsr_copy(cmat_x,qs_ot_env%rot_mat_evec,error=error)
       CALL cp_dbcsr_scale_by_vector(cmat_x,alpha=evals_exp,side='right',error=error)
       CALL cp_dbcsr_multiply('N','C',cone,cmat_x,qs_ot_env%rot_mat_evec,czero,cmat_u,error=error)
       CALL cp_dbcsr_copy(qs_ot_env%rot_mat_u, cmat_u, keep_imaginary=.FALSE., error=error)
       CALL cp_dbcsr_release(cmat_x, error=error)
       CALL cp_dbcsr_release(cmat_u, error=error)
       DEALLOCATE(evals_exp)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_ot_generate_rotation

! *****************************************************************************
!> \brief computes the derivative fields with respect to rot_mat_x
!> \param valid qs_ot_env. In particular qs_ot_generate_rotation has to be called before
!>                        and the rot_mat_dedu matrix has to be up to date
!> \par History
!>      08.2004 created [ Joost VandeVondele ]
! *****************************************************************************
  SUBROUTINE qs_ot_rot_mat_derivative(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_rot_mat_derivative', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), PARAMETER              :: cI = (0.0_dp,1.0_dp), &
                                                cone = (1.0_dp,0.0_dp), &
                                                czero = (0.0_dp,0.0_dp)

    INTEGER                                  :: handle, i, j, k
    REAL(KIND=dp)                            :: e1, e2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: data_d
    TYPE(cp_dbcsr_type)                          :: cmat_buf1, cmat_buf2
    TYPE(cp_dbcsr_iterator)                     :: iter
    COMPLEX(dp), DIMENSION(:,:), POINTER :: data_z
    INTEGER::row, col, blk,row_offset, col_offset,row_size, col_size
    LOGICAL::found

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(qs_ot_env%rot_mat_u,nfullrows_total=k)
    IF (k/=0) THEN
       CALL cp_dbcsr_copy(qs_ot_env%matrix_buf1,qs_ot_env%rot_mat_dedu,error=error)
       ! now we get to the derivative wrt the antisymmetric matrix rot_mat_x
       CALL cp_dbcsr_init(cmat_buf1, error=error)
       CALL cp_dbcsr_init(cmat_buf2, error=error)
       CALL cp_dbcsr_copy(cmat_buf1,qs_ot_env%rot_mat_evec,"cmat_buf1",error=error)
       CALL cp_dbcsr_copy(cmat_buf2,qs_ot_env%rot_mat_evec,"cmat_buf2",error=error)

       ! init cmat_buf1
       !CALL cp_fm_get_info(qs_ot_env%matrix_buf1,matrix_struct=fm_struct, local_data=local_data_r,error=error)
       !CALL cp_cfm_get_info(cmat_buf1, nrow_local=nrow_local,   ncol_local=ncol_local, &
       !     row_indices=row_indices, col_indices=col_indices, &
       !     local_data=local_data_c,error=error)
       !local_data_c=local_data_r

       CALL cp_dbcsr_iterator_start(iter, cmat_buf1)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, row, col, data_z, blk)
          CALL cp_dbcsr_get_block_p(qs_ot_env%matrix_buf1, row, col, data_d, found)
          data_z=data_d
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)

       CALL cp_dbcsr_multiply('T','N',cone,cmat_buf1,qs_ot_env%rot_mat_evec,&
            czero,cmat_buf2,error=error)
       CALL cp_dbcsr_multiply('C','N',cone,qs_ot_env%rot_mat_evec,cmat_buf2,&
            czero,cmat_buf1,error=error)

       CALL cp_dbcsr_iterator_start(iter, cmat_buf1)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row, col, data_z,blk,&
               row_size=row_size, col_size=col_size, &
               row_offset=row_offset, col_offset=col_offset)
          DO j=1,col_size
          DO i=1,row_size
             e1=qs_ot_env%rot_mat_evals(row_offset+i-1)
             e2=qs_ot_env%rot_mat_evals(col_offset+j-1)
             data_z(i,j)=data_z(i,j)*cint(e1,e2)
          ENDDO
          ENDDO
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)

       CALL cp_dbcsr_multiply('N','N',cone,qs_ot_env%rot_mat_evec,cmat_buf1,&
            czero,cmat_buf2,error=error)
       CALL cp_dbcsr_multiply('N','C',cone,cmat_buf2,qs_ot_env%rot_mat_evec,&
            czero,cmat_buf1,error=error)

       CALL cp_dbcsr_copy(qs_ot_env%matrix_buf1,cmat_buf1,error=error)

       CALL cp_dbcsr_transposed(qs_ot_env%matrix_buf2,qs_ot_env%matrix_buf1,&
            shallow_data_copy=.FALSE., use_distribution=cp_dbcsr_distribution(qs_ot_env%matrix_buf3), &
            transpose_distribution=.FALSE.,error=error)
       CALL cp_dbcsr_add(qs_ot_env%matrix_buf1,qs_ot_env%matrix_buf2,&
            alpha_scalar=-1.0_dp,beta_scalar=+1.0_dp,error=error)
       CALL cp_dbcsr_copy(qs_ot_env%rot_mat_gx,qs_ot_env%matrix_buf1,error=error)
       CALL cp_dbcsr_release(cmat_buf1, error=error)
       CALL cp_dbcsr_release(cmat_buf2, error=error)
    END IF
    CALL timestop(handle)
  CONTAINS
! *****************************************************************************
    FUNCTION cint(e1,e2)
    REAL(KIND=dp)                            :: e1, e2
    COMPLEX(KIND=dp)                         :: cint

    COMPLEX(KIND=dp)                         :: l1, l2, x
    INTEGER                                  :: I

      l1=(0.0_dp,-1.0_dp)*e1
      l2=(0.0_dp,-1.0_dp)*e2
      IF (ABS(l1-l2) .GT. 0.5_dp) THEN
         cint=(EXP(l1)-EXP(l2))/(l1-l2)
      ELSE
         x=1.0_dp
         cint=0.0_dp
         DO I=1,16
            cint=cint+x
            x=x*(l1-l2)/REAL(I+1,KIND=dp)
         ENDDO
         cint=cint*EXP(l2)
      ENDIF
    END FUNCTION cint
  END SUBROUTINE qs_ot_rot_mat_derivative

  !
  ! decide strategy
  ! tries to decide if the taylor expansion of cos(sqrt(xsx)) converges rapidly enough
  ! to make a taylor expansion of the functions cos(sqrt(xsx)) and sin(sqrt(xsx))/sqrt(xsx)
  ! and their derivatives faster than their computation based on diagonalization
  ! since xsx can be very small, especially during dynamics, only a few terms might indeed be needed
  ! we find the necessary order N to have largest_eval_upper_bound**(N+1)/(2(N+1))! < eps_taylor
  !
! *****************************************************************************
  SUBROUTINE decide_strategy(qs_ot_env)
    TYPE(qs_ot_type)                         :: qs_ot_env

    INTEGER                                  :: N
    REAL(KIND=dp)                            :: num_error

    qs_ot_env % do_taylor = .FALSE.
    N=0
    num_error=qs_ot_env % largest_eval_upper_bound / ( 2.0_dp )
    DO WHILE (num_error > qs_ot_env % settings % eps_taylor .AND. N <= qs_ot_env % settings % max_taylor)
       N=N+1
       num_error=num_error * qs_ot_env % largest_eval_upper_bound / REAL(( 2*N+1 )*(2*N+2),KIND=dp)
    END DO
    qs_ot_env % taylor_order = N
    IF ( qs_ot_env % taylor_order <= qs_ot_env % settings % max_taylor) THEN
       qs_ot_env % do_taylor = .TRUE.
    ENDIF

  END SUBROUTINE decide_strategy

  ! c=(c0*cos(p^0.5)+x*sin(p^0.5)*p^(-0.5)) x rot_mat_u
  ! this assumes that x is already ortho to S*C0, and that p is x*S*x
  ! rot_mat_u is an optional rotation matrix
! *****************************************************************************
  SUBROUTINE qs_ot_get_orbitals(matrix_c, matrix_x, qs_ot_env, error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix_c, matrix_x
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_orbitals', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n
    TYPE(cp_dbcsr_type), POINTER             :: matrix_kk

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(matrix_x,nfullrows_total=n,nfullcols_total=k)

    ! rotate the multiplying matrices cosp and sinp instead of the result,
    ! this should be cheaper for large basis sets
    IF (qs_ot_env%settings%do_rotation) THEN
       matrix_kk => qs_ot_env%matrix_buf1
       CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_cosp, &
            qs_ot_env%rot_mat_u,rzero,matrix_kk,error=error)
    ELSE
       matrix_kk => qs_ot_env%matrix_cosp
    ENDIF

    CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_c0,matrix_kk, &
         rzero,matrix_c,error=error)

    IF (qs_ot_env%settings%do_rotation) THEN
       matrix_kk => qs_ot_env%matrix_buf1
       CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_sinp, &
            qs_ot_env%rot_mat_u,rzero,matrix_kk,error=error)
    ELSE
       matrix_kk => qs_ot_env%matrix_sinp
    ENDIF
    CALL cp_dbcsr_multiply('N','N',rone,matrix_x,matrix_kk, &
         rone ,matrix_c,error=error)

    CALL timestop(handle)

  END SUBROUTINE qs_ot_get_orbitals

! *****************************************************************************
  SUBROUTINE qs_ot_get_scp_nddo_coeffs ( qs_ot_env, pscp, error )
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_dbcsr_type), POINTER             :: pscp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL cp_dbcsr_copy(pscp,qs_ot_env%xmat,error=error)

  END SUBROUTINE qs_ot_get_scp_nddo_coeffs
! *****************************************************************************
  SUBROUTINE qs_ot_get_scp_nddo_derivative ( qs_ot_env, pscp, fscp, error )
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_dbcsr_type), POINTER             :: pscp, fscp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL cp_dbcsr_copy(qs_ot_env%xmat,pscp,error=error)
    CALL cp_dbcsr_copy(qs_ot_env%gxmat,fscp,error=error)

  END SUBROUTINE qs_ot_get_scp_nddo_derivative
! *****************************************************************************
  ! this routines computes dE/dx=dx, with dx ortho to sc0
  ! needs dE/dC=hc,C0,X,SX,p
  ! if preconditioned it will not be the derivative, but the lagrangian multiplier
  ! is changed so that P*dE/dx is the right derivative (i.e. in the allowed subspace)
! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       qs_ot_env,error)
    TYPE(cp_dbcsr_type), POINTER             :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n, ortho_k
    TYPE(cp_dbcsr_type), POINTER             :: matrix_hc_local, matrix_target

    CALL timeset(routineN,handle)

    NULLIFY(matrix_hc_local)

    CALL cp_dbcsr_get_info(matrix_x,nfullrows_total=n,nfullcols_total=k)

    ! could in principle be taken inside qs_ot_get_derivative_* for increased efficiency
    ! create a local rotated version of matrix_hc leaving matrix_hc untouched (needed
    ! for lagrangian multipliers)
    IF (qs_ot_env % settings % do_rotation) THEN
       CALL cp_dbcsr_copy(matrix_gx,matrix_hc,error=error) ! use gx as temporary
       CALL cp_dbcsr_init_p(matrix_hc_local, error=error)
       CALL cp_dbcsr_copy(matrix_hc_local,matrix_hc,name='matrix_hc_local',error=error)
       CALL cp_dbcsr_set(matrix_hc_local,0.0_dp,error=error)
       CALL cp_dbcsr_multiply('N','T',rone,matrix_gx,qs_ot_env%rot_mat_u,rzero,matrix_hc_local,error=error)
    ELSE
       matrix_hc_local=>matrix_hc
    ENDIF

    IF (qs_ot_env % do_taylor) THEN
       CALL qs_ot_get_derivative_taylor(matrix_hc_local,matrix_x,matrix_sx,matrix_gx,qs_ot_env,error=error)
    ELSE
       CALL qs_ot_get_derivative_diag(matrix_hc_local,matrix_x,matrix_sx,matrix_gx,qs_ot_env,error=error)
    ENDIF

    ! and make it orthogonal
    CALL cp_dbcsr_get_info(qs_ot_env%matrix_sc0,nfullcols_total=ortho_k)

    IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
       matrix_target => qs_ot_env%matrix_psc0
    ELSE
       matrix_target => qs_ot_env%matrix_sc0
    ENDIF
    ! first make the matrix os if not yet valid
    IF (.NOT. qs_ot_env%os_valid) THEN
       ! this assumes that the preconditioner is a single matrix
       ! that maps sc0 onto psc0

       IF (ASSOCIATED(qs_ot_env%preconditioner)) THEN
          CALL apply_preconditioner(qs_ot_env%preconditioner, qs_ot_env%matrix_sc0, &
               qs_ot_env%matrix_psc0 ,error=error)
       ENDIF
       CALL cp_dbcsr_multiply('T','N',rone,&
            qs_ot_env%matrix_sc0,matrix_target, &
            rzero,qs_ot_env%matrix_os,&
            error=error)
       CALL cp_dbcsr_cholesky_decompose(qs_ot_env%matrix_os,&
            para_env=qs_ot_env%para_env,blacs_env=qs_ot_env%blacs_env,error=error)
       CALL cp_dbcsr_cholesky_invert(qs_ot_env%matrix_os,&
            para_env=qs_ot_env%para_env,blacs_env=qs_ot_env%blacs_env,&
            upper_to_full=.TRUE.,error=error)
       qs_ot_env%os_valid=.TRUE.
    ENDIF
    CALL cp_dbcsr_multiply('T','N',rone,matrix_target,matrix_gx, &
         rzero,qs_ot_env%matrix_buf1_ortho, error=error)
    CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_os,&
         qs_ot_env%matrix_buf1_ortho, rzero,qs_ot_env%matrix_buf2_ortho,error=error)
    CALL cp_dbcsr_multiply('N','N',-rone,qs_ot_env%matrix_sc0, &
         qs_ot_env%matrix_buf2_ortho, rone,matrix_gx,error=error)
    ! also treat the rot_mat gradient here
    IF (qs_ot_env%settings%do_rotation) THEN
       CALL qs_ot_rot_mat_derivative(qs_ot_env,error=error)
    ENDIF

    IF (qs_ot_env % settings % do_rotation) THEN
       CALL cp_dbcsr_release_p(matrix_hc_local, error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE qs_ot_get_derivative

! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative_diag(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       qs_ot_env,error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, k, n

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(matrix_x,nfullrows_total=n,nfullcols_total=k)

    ! go for the derivative now
    ! this de/dc*(dX/dx)*sinp
    CALL cp_dbcsr_multiply('N','N',rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_gx,&
         error=error)
    ! overlap hc*x
    CALL cp_dbcsr_multiply('T','N',rone,matrix_hc,matrix_x,rzero,qs_ot_env%matrix_buf2,&
         error=error)
    ! get it in the basis of the eigenvectors
    CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r,&
         rzero,qs_ot_env%matrix_buf1,error=error)
    CALL cp_dbcsr_multiply('T','N',rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_buf2,error=error)

    ! get the schur product of O_uv*B_uv
    CALL cp_dbcsr_hadamard_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_sinp_b, &
         qs_ot_env%matrix_buf3,error=error)

    ! overlap hc*c0
    CALL cp_dbcsr_multiply('T','N',rone,matrix_hc,qs_ot_env%matrix_c0,rzero, &
         qs_ot_env%matrix_buf2,error=error)
    ! get it in the basis of the eigenvectors
    CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_buf2,qs_ot_env%matrix_r, &
         rzero,qs_ot_env%matrix_buf1,error=error)
    CALL cp_dbcsr_multiply('T','N',rone,qs_ot_env%matrix_r, qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_buf2,error=error)

    CALL cp_dbcsr_hadamard_product(qs_ot_env%matrix_buf2,qs_ot_env%matrix_cosp_b, &
         qs_ot_env%matrix_buf4,error=error)

    ! add the two bs and compute b+b^T
    CALL cp_dbcsr_add(qs_ot_env%matrix_buf3,qs_ot_env%matrix_buf4,&
         alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)

    ! get the b in the eigenvector basis
    CALL cp_dbcsr_multiply('N','T',rone,qs_ot_env%matrix_buf3,qs_ot_env%matrix_r, &
         rzero,qs_ot_env%matrix_buf1,error=error)
    CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_buf3,error=error)
    CALL cp_dbcsr_transposed(qs_ot_env%matrix_buf1,qs_ot_env%matrix_buf3,&
         shallow_data_copy=.FALSE., use_distribution=cp_dbcsr_distribution(qs_ot_env%matrix_buf3), &
         transpose_distribution=.FALSE.,error=error)
    CALL cp_dbcsr_add(qs_ot_env%matrix_buf3,qs_ot_env%matrix_buf1,&
         alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)

    ! and add to the derivative
    CALL cp_dbcsr_multiply('N','N',rone,matrix_sx,qs_ot_env%matrix_buf3, &
         rone,matrix_gx,error=error)
    CALL timestop(handle)

  END SUBROUTINE qs_ot_get_derivative_diag

  ! compute the derivative of the taylor expansion below
! *****************************************************************************
  SUBROUTINE qs_ot_get_derivative_taylor(matrix_hc,matrix_x,matrix_sx,matrix_gx, &
       qs_ot_env, error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix_hc, matrix_x, &
                                                matrix_sx, matrix_gx
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_get_derivative_taylor', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp , rzero = 0.0_dp

    INTEGER                                  :: handle, i, k, n
    REAL(KIND=dp)                            :: cosfactor, sinfactor
    TYPE(cp_dbcsr_type), POINTER             :: matrix_left, matrix_right

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(matrix_x,nfullrows_total=n,nfullcols_total=k)

    ! go for the derivative now
    ! this de/dc*(dX/dx)*sinp i.e. zeroth order
    CALL cp_dbcsr_multiply('N','N',rone,matrix_hc,qs_ot_env%matrix_sinp,rzero,matrix_gx,&
         error=error)

    IF (qs_ot_env % taylor_order .LE. 0) THEN
       CALL timestop(handle)
       RETURN
    ENDIF

    ! we store the matrix that will multiply sx in matrix_r
    CALL cp_dbcsr_set(qs_ot_env%matrix_r,rzero,error=error)

    ! just better names for matrix_cosp_b and matrix_sinp_b (they are buffer space here)
    matrix_left  => qs_ot_env%matrix_cosp_b
    matrix_right => qs_ot_env%matrix_sinp_b

    ! overlap hc*x and add its transpose to matrix_left
    CALL cp_dbcsr_multiply('T','N',rone,matrix_hc,matrix_x,rzero,matrix_left,&
         error=error)
    CALL cp_dbcsr_transposed(qs_ot_env%matrix_buf1,matrix_left,&
         shallow_data_copy=.FALSE., use_distribution=cp_dbcsr_distribution(matrix_left), &
         transpose_distribution=.FALSE., error=error)
    CALL cp_dbcsr_add(matrix_left,qs_ot_env%matrix_buf1,&
         alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)
    CALL cp_dbcsr_copy(matrix_right,matrix_left,error=error)

    ! first order
    sinfactor=-1.0_dp/(2.0_dp*3.0_dp)
    CALL cp_dbcsr_add(qs_ot_env%matrix_r,matrix_left,&
    alpha_scalar=1.0_dp,beta_scalar=sinfactor,error=error)

    !      M
    !    OM+MO
    ! OOM+OMO+MOO
    !   ...
    DO i=2, qs_ot_env % taylor_order
       sinfactor=sinfactor * (-1.0_dp)/REAL(2*i * (2*i+1),KIND=dp)
       CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_p,matrix_left,rzero,qs_ot_env%matrix_buf1,error=error)
       CALL cp_dbcsr_multiply('N','N',rone,matrix_right,qs_ot_env%matrix_p,rzero,matrix_left,error=error)
       CALL cp_dbcsr_copy(matrix_right,matrix_left,error=error)
       CALL cp_dbcsr_add(matrix_left,qs_ot_env%matrix_buf1,&
            1.0_dp,1.0_dp,error=error)
       CALL cp_dbcsr_add(qs_ot_env%matrix_r,matrix_left,&
            alpha_scalar=1.0_dp,beta_scalar=sinfactor,error=error)
    ENDDO

    ! overlap hc*c0 and add its transpose to matrix_left
    CALL cp_dbcsr_multiply('T','N',rone,matrix_hc,qs_ot_env%matrix_c0,rzero,matrix_left,error=error)
    CALL cp_dbcsr_transposed(qs_ot_env%matrix_buf1,matrix_left,&
         shallow_data_copy=.FALSE., use_distribution=cp_dbcsr_distribution(matrix_left),&
         transpose_distribution=.FALSE., error=error)
    CALL cp_dbcsr_add(matrix_left,qs_ot_env%matrix_buf1,1.0_dp,1.0_dp,error=error)
    CALL cp_dbcsr_copy(matrix_right,matrix_left,error=error)

    ! first order
    cosfactor=-1.0_dp/(1.0_dp*2.0_dp)
    CALL cp_dbcsr_add(qs_ot_env%matrix_r,matrix_left,&
         alpha_scalar=1.0_dp,beta_scalar=cosfactor,error=error)

    !      M
    !    OM+MO
    ! OOM+OMO+MOO
    !   ...
    DO i=2, qs_ot_env % taylor_order
       cosfactor=cosfactor * (-1.0_dp)/REAL(2*i * (2*i-1),KIND=dp)
       CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_p,matrix_left,rzero,qs_ot_env%matrix_buf1,error=error)
       CALL cp_dbcsr_multiply('N','N',rone,matrix_right,qs_ot_env%matrix_p,rzero,matrix_left,error=error)
       CALL cp_dbcsr_copy(matrix_right,matrix_left,error=error)
       CALL cp_dbcsr_add(matrix_left,qs_ot_env%matrix_buf1,1.0_dp,1.0_dp,error=error)
       CALL cp_dbcsr_add(qs_ot_env%matrix_r,matrix_left,&
            alpha_scalar=1.0_dp,beta_scalar=cosfactor,error=error)
    ENDDO

    ! and add to the derivative
    CALL cp_dbcsr_multiply('N','N',rone,matrix_sx,qs_ot_env%matrix_r,rone,matrix_gx,error=error)

    CALL timestop(handle)

  END SUBROUTINE qs_ot_get_derivative_taylor

  ! computes a taylor expansion.
! *****************************************************************************
  SUBROUTINE qs_ot_p2m_taylor(qs_ot_env,error)
    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_p2m_taylor', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: handle, i, k
    REAL(KIND=dp)                            :: cosfactor, sinfactor

    CALL timeset(routineN,handle)

    ! zeroth order
    CALL cp_dbcsr_set(qs_ot_env%matrix_cosp,rzero,error=error)
    CALL cp_dbcsr_set(qs_ot_env%matrix_sinp,rzero,error=error)
    CALL cp_dbcsr_add_on_diag(qs_ot_env%matrix_cosp,rone,error=error)
    CALL cp_dbcsr_add_on_diag(qs_ot_env%matrix_sinp,rone,error=error)

    IF (qs_ot_env% taylor_order .LE. 0) THEN
       CALL timestop(handle)
       RETURN
    ENDIF

    ! first order
    cosfactor=-1.0_dp/(1.0_dp*2.0_dp)
    sinfactor=-1.0_dp/(2.0_dp*3.0_dp)
    CALL cp_dbcsr_add(qs_ot_env%matrix_cosp,qs_ot_env%matrix_p,alpha_scalar=1.0_dp,beta_scalar=cosfactor,error=error)
    CALL cp_dbcsr_add(qs_ot_env%matrix_sinp,qs_ot_env%matrix_p,alpha_scalar=1.0_dp,beta_scalar=sinfactor,error=error)
    IF (qs_ot_env% taylor_order .LE. 1) THEN
       CALL timestop(handle)
       RETURN
    ENDIF

    ! other orders
    CALL cp_dbcsr_get_info(qs_ot_env%matrix_p,nfullrows_total=k)
    CALL cp_dbcsr_copy(qs_ot_env%matrix_r,qs_ot_env%matrix_p,error=error)

    DO i=2, qs_ot_env%taylor_order
       ! new power of p
       CALL cp_dbcsr_multiply('N','N',rone,qs_ot_env%matrix_p,qs_ot_env%matrix_r,&
            rzero,qs_ot_env%matrix_buf1,error=error)
       CALL cp_dbcsr_copy(qs_ot_env%matrix_r,qs_ot_env%matrix_buf1,error=error)
       ! add to the taylor expansion so far
       cosfactor=cosfactor * (-1.0_dp)/REAL(2*i * (2*i-1),KIND=dp)
       sinfactor=sinfactor * (-1.0_dp)/REAL(2*i * (2*i+1),KIND=dp)
       CALL cp_dbcsr_add(qs_ot_env%matrix_cosp,qs_ot_env%matrix_r,&
            alpha_scalar=1.0_dp,beta_scalar=cosfactor,error=error)
       CALL cp_dbcsr_add(qs_ot_env%matrix_sinp,qs_ot_env%matrix_r,&
            alpha_scalar=1.0_dp,beta_scalar=sinfactor,error=error)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE qs_ot_p2m_taylor

  ! given p, computes  - eigenstuff (matrix_r,evals)
  !                    - cos(p^0.5),p^(-0.5)*sin(p^0.5)
  !                    - the real b matrices, needed for the derivatives of these guys
  !                    cosp_b_ij=(1/(2pii) * int(cos(z^1/2)/((z-eval(i))*(z-eval(j))))
  !                    sinp_b_ij=(1/(2pii) * int(z^(-1/2)*sin(z^1/2)/((z-eval(i))*(z-eval(j))))
! *****************************************************************************
  SUBROUTINE qs_ot_p2m_diag(qs_ot_env,error)

    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_p2m_diag', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rone = 1.0_dp, rzero = 0.0_dp

    INTEGER                                  :: blk, col, col_offset, &
                                                col_size, handle, i, j, k, &
                                                row, row_offset, row_size
    REAL(dp), DIMENSION(:, :), POINTER       :: DATA
    REAL(KIND=dp)                            :: a, b
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_get_info(qs_ot_env%matrix_p,nfullrows_total=k)
    CALL cp_dbcsr_copy(qs_ot_env%matrix_buf1,qs_ot_env%matrix_p,error=error)
    CALL cp_dbcsr_syevd(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,qs_ot_env%evals,&
         qs_ot_env%para_env,qs_ot_env%blacs_env,error=error)
    DO i=1,k
       qs_ot_env%evals(i)=MAX(0.0_dp,qs_ot_env%evals(i))
    ENDDO

    !$OMP PARALLEL DO
    DO i=1,k
       qs_ot_env%dum(i)=COS(SQRT(qs_ot_env%evals(i)))
    ENDDO
    CALL cp_dbcsr_copy(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,error=error)
    CALL cp_dbcsr_scale_by_vector(qs_ot_env%matrix_buf1,alpha=qs_ot_env%dum,side='right',error=error)
    CALL cp_dbcsr_multiply('N','T',rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_cosp,error=error)

    !$OMP PARALLEL DO
    DO i=1,k
       qs_ot_env%dum(i)=qs_ot_sinc(SQRT(qs_ot_env%evals(i)))
    ENDDO
    CALL cp_dbcsr_copy(qs_ot_env%matrix_buf1,qs_ot_env%matrix_r,error=error)
    CALL cp_dbcsr_scale_by_vector(qs_ot_env%matrix_buf1,alpha=qs_ot_env%dum,side='right',error=error)
    CALL cp_dbcsr_multiply('N','T',rone,qs_ot_env%matrix_r,qs_ot_env%matrix_buf1, &
         rzero,qs_ot_env%matrix_sinp,error=error)

!!$OMP PARALLEL DO PRIVATE(i,j,a,b)
!DO j=1,ncol_local
!   DO i=1,nrow_local
!      a=(SQRT(qs_ot_env%evals(row_indices(i))) &
!           -SQRT(qs_ot_env%evals(col_indices(j))))/2.0_dp
!      b=(SQRT(qs_ot_env%evals(row_indices(i))) &
!           +SQRT(qs_ot_env%evals(col_indices(j))))/2.0_dp
!      qs_ot_env%matrix_cosp_b%local_data(i,j) = -0.5_dp*qs_ot_sinc(a)*qs_ot_sinc(b)
!   ENDDO
!ENDDO

    CALL cp_dbcsr_copy(qs_ot_env%matrix_cosp_b,qs_ot_env%matrix_cosp,error=error)
    CALL cp_dbcsr_iterator_start(iter, qs_ot_env%matrix_cosp_b)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, DATA,&
            block_number=blk, row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO j=1,col_size
       DO i=1,row_size
          a=(SQRT(qs_ot_env%evals( row_offset + i - 1 )) &
            -SQRT(qs_ot_env%evals( col_offset + j - 1 )))/2.0_dp
          b=(SQRT(qs_ot_env%evals( row_offset + i - 1 )) &
            +SQRT(qs_ot_env%evals( col_offset + j - 1 )))/2.0_dp
          DATA(i,j) = -0.5_dp*qs_ot_sinc(a)*qs_ot_sinc(b)
       ENDDO
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

!!$OMP PARALLEL DO PRIVATE(i,j,a,b)
!DO j=1,ncol_local
!   DO i=1,nrow_local
!      a=SQRT(qs_ot_env%evals(row_indices(i)))
!      b=SQRT(qs_ot_env%evals(col_indices(j)))
!      qs_ot_env%matrix_sinp_b%local_data(i,j)=qs_ot_sincf(a,b)
!   ENDDO
!ENDDO

    CALL cp_dbcsr_copy(qs_ot_env%matrix_sinp_b,qs_ot_env%matrix_sinp,error=error)
    CALL cp_dbcsr_iterator_start(iter, qs_ot_env%matrix_sinp_b)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, row, col, DATA,&
            block_number=blk, row_size=row_size, col_size=col_size, &
            row_offset=row_offset, col_offset=col_offset)
       DO j=1,col_size
       DO i=1,row_size
          a=SQRT(qs_ot_env%evals( row_offset + i - 1 ))
          b=SQRT(qs_ot_env%evals( col_offset + j - 1 ))
          DATA(i,j)=qs_ot_sincf(a,b)
       ENDDO
       ENDDO
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)

  END SUBROUTINE qs_ot_p2m_diag

  ! computes sin(x)/x for all values of the argument
! *****************************************************************************
  FUNCTION qs_ot_sinc(x)

    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: qs_ot_sinc

    REAL(KIND=dp), PARAMETER :: q1 = 1.0_dp, q2 = -q1/(2.0_dp *3.0_dp), &
      q3 = -q2/(4.0_dp *5.0_dp), q4 = -q3/(6.0_dp *7.0_dp), &
      q5 = -q4/(8.0_dp *9.0_dp), q6 = -q5/(10.0_dp*11.0_dp), &
      q7 = -q6/(12.0_dp*13.0_dp), q8 = -q7/(14.0_dp*15.0_dp), &
      q9 = -q8/(16.0_dp*17.0_dp), q10 = -q9/(18.0_dp*19.0_dp)

    REAL(KIND=dp)                            :: y

    IF (ABS(x)>0.5_dp) THEN
       qs_ot_sinc=SIN(x)/x
    ELSE
       y=x*x
       qs_ot_sinc=q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*(q9+y*(q10)))))))))
    ENDIF
  END FUNCTION qs_ot_sinc
  ! computes (1/(x^2-y^2))*(sinc(x)-sinc(y)) for all positive values of the arguments
! *****************************************************************************
  FUNCTION qs_ot_sincf(xa,ya)

    REAL(KIND=dp), INTENT(IN)                :: xa, ya
    REAL(KIND=dp)                            :: qs_ot_sincf

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ot_sincf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: a, b, rs, sf, x, xs, y, ybx, &
                                                ybxs

! this is currently a limit of the routine, could be removed rather easily

    IF (xa.lt.0) CALL stop_program(routineN,moduleN,__LINE__,"x is negative")
    IF (ya.lt.0) CALL stop_program(routineN,moduleN,__LINE__,"y is negative")

    IF (xa.lt.ya) THEN
       x=ya
       y=xa
    ELSE
       x=xa
       y=ya
    ENDIF

    IF ( x .LT. 0.5_dp ) THEN ! use series, keeping in mind that x,y,x+y,x-y can all be zero

       qs_ot_sincf=0.0_dp
       IF (x .GT. 0.0_dp) THEN
          ybx=y/x
       ELSE ! should be irrelevant  !?
          ybx=0.0_dp
       ENDIF

       sf=-1.0_dp/((1.0_dp+ybx)*6.0_dp)
       rs=1.0_dp
       ybxs=ybx
       xs=1.0_dp

       DO i=1,10
          qs_ot_sincf=qs_ot_sincf+sf*rs*xs*(1.0_dp+ybxs)
          sf=-sf/(REAL((2*i+2),dp)*REAL((2*i+3),dp))
          rs=rs+ybxs
          ybxs=ybxs*ybx
          xs=xs*x*x
       ENDDO

    ELSE ! no series expansion
       IF ( x-y .GT. 0.1_dp ) THEN  ! safe to use the normal form
          qs_ot_sincf=(qs_ot_sinc(x)-qs_ot_sinc(y))/((x+y)*(x-y))
       ELSE
          a=(x+y)/2.0_dp
          b=(x-y)/2.0_dp ! might be close to zero
          ! y (=(a-b)) can not be close to zero since it is close to x>0.5
          qs_ot_sincf=(qs_ot_sinc(b)*COS(a)-qs_ot_sinc(a)*COS(b))/(2*x*y)
       ENDIF
    ENDIF

  END FUNCTION qs_ot_sincf

END MODULE qs_ot
