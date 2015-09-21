!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module that contains the algorithms to perform an itrative
!>         diagonalization by the block-Lanczos approach
!> \par History
!>      05.2009 created [MI]
!> \author fawzi
! *****************************************************************************
MODULE qs_scf_lanczos

  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_qr_factorization,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_transpose,&
                                             cp_fm_triangular_multiply
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_diag,                      ONLY: choose_eigv_solver
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_submatrix, cp_fm_p_type, cp_fm_release, &
       cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type, &
       cp_fm_vectorsnorm
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_scf_types,                    ONLY: krylov_space_type
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_lanczos'

  PUBLIC :: krylov_space_allocate, lanczos_refinement, lanczos_refinement_2v

CONTAINS

! *****************************************************************************

! *****************************************************************************
!> \brief  allocates matrices and vectros used in the construction of
!>        the krylov space and for the lanczos refinement
!> \param krylov_space ...
!> \param scf_control ...
!> \param mos ...
!> \param
!> \par History
!>      05.2009 created [MI]
! *****************************************************************************

  SUBROUTINE krylov_space_allocate(krylov_space, scf_control, mos)

    TYPE(krylov_space_type), POINTER         :: krylov_space
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CHARACTER(LEN=*), PARAMETER :: routineN = 'krylov_space_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ik, ispin, max_nmo, &
                                                nao, nblock, ndim, nk, nmo, &
                                                nspin
    LOGICAL                                  :: failure
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    failure = .FALSE.

    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(krylov_space))

    IF(.NOT. ASSOCIATED(krylov_space%mo_conv)) THEN
        NULLIFY(fm_struct_tmp, mo_coeff)

        krylov_space%nkrylov = scf_control%diagonalization%nkrylov
        krylov_space%nblock = scf_control%diagonalization%nblock_krylov

        nk =  krylov_space%nkrylov
        nblock =  krylov_space%nblock
        nspin = SIZE(mos,1)

        ALLOCATE(krylov_space%mo_conv(nspin))
        ALLOCATE(krylov_space%mo_refine(nspin))
        ALLOCATE(krylov_space%chc_mat(nspin))
        ALLOCATE(krylov_space%c_vec(nspin))
        max_nmo = 0
        DO ispin = 1,nspin
           CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, nao=nao, nmo=nmo)
           CALL cp_fm_create(krylov_space%mo_conv(ispin)%matrix,mo_coeff%matrix_struct)
           CALL cp_fm_create(krylov_space%mo_refine(ispin)%matrix,mo_coeff%matrix_struct)
           NULLIFY(fm_struct_tmp)
           CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmo, ncol_global=nmo,&
                                       para_env=mo_coeff%matrix_struct%para_env, &
                                       context=mo_coeff%matrix_struct%context)
           CALL cp_fm_create(krylov_space%chc_mat(ispin)%matrix,fm_struct_tmp,"chc")
           CALL cp_fm_create(krylov_space%c_vec(ispin)%matrix,fm_struct_tmp,"vec")
           CALL cp_fm_struct_release(fm_struct_tmp)
           max_nmo = MAX(max_nmo,nmo)
        END DO


        !the use of max_nmo might not be ok, in this case allocate nspin matrices
        ALLOCATE(krylov_space%c_eval(max_nmo))

        ALLOCATE(krylov_space%v_mat(nk))

        NULLIFY(fm_struct_tmp)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nblock,&
                                 para_env=mo_coeff%matrix_struct%para_env, &
                                 context=mo_coeff%matrix_struct%context)
        DO ik = 1,nk
          CALL cp_fm_create(krylov_space%v_mat(ik)%matrix,matrix_struct=fm_struct_tmp,&
                name="v_mat_"//TRIM(ADJUSTL(cp_to_string(ik))))
        END DO
        CALL cp_fm_create(krylov_space%tmp_mat,matrix_struct=fm_struct_tmp,&
                name="tmp_mat")
        CALL cp_fm_struct_release(fm_struct_tmp)

    ! NOTE: the following matrices are small and could be defined
!           as standard array rather than istributed fm
        NULLIFY(fm_struct_tmp)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nblock, ncol_global=nblock,&
                                 para_env=mo_coeff%matrix_struct%para_env, &
                                 context=mo_coeff%matrix_struct%context)
        CALL cp_fm_create(krylov_space%block1_mat,matrix_struct=fm_struct_tmp,&
                name="a_mat_"//TRIM(ADJUSTL(cp_to_string(ik))))
        CALL cp_fm_create(krylov_space%block2_mat,matrix_struct=fm_struct_tmp,&
                name="b_mat_"//TRIM(ADJUSTL(cp_to_string(ik))))
        CALL cp_fm_create(krylov_space%block3_mat,matrix_struct=fm_struct_tmp,&
                name="b2_mat_"//TRIM(ADJUSTL(cp_to_string(ik))))
        CALL cp_fm_struct_release(fm_struct_tmp)

        ndim = nblock*nk
        NULLIFY(fm_struct_tmp)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ndim, ncol_global=ndim,&
                                 para_env=mo_coeff%matrix_struct%para_env, &
                                 context=mo_coeff%matrix_struct%context)
        CALL cp_fm_create(krylov_space%block4_mat,matrix_struct=fm_struct_tmp,&
                name="t_mat")
        CALL cp_fm_create(krylov_space%block5_mat,matrix_struct=fm_struct_tmp,&
                name="t_vec")
        CALL cp_fm_struct_release(fm_struct_tmp)
        ALLOCATE(krylov_space%t_eval(ndim))

    ELSE
       !Nothing should be done
    END IF

    CALL timestop(handle)

  END SUBROUTINE krylov_space_allocate


! *****************************************************************************
!> \brief lanczos refinement by blocks of not-converged MOs
!> \param krylov_space ...
!> \param ks ...
!> \param c0 ...
!> \param c1 ...
!> \param eval ...
!> \param nao ...
!> \param eps_iter ...
!> \param ispin ...
!> \param check_moconv_only ...
!> \param
!> \par History
!>      05.2009 created [MI]
! *****************************************************************************

  SUBROUTINE lanczos_refinement(krylov_space, ks, c0, c1, eval,  nao,  &
              eps_iter, ispin, check_moconv_only)

    TYPE(krylov_space_type), POINTER         :: krylov_space
    TYPE(cp_fm_type), POINTER                :: ks, c0, c1
    REAL(dp), DIMENSION(:), POINTER          :: eval
    INTEGER, INTENT(IN)                      :: nao
    REAL(dp), INTENT(IN)                     :: eps_iter
    INTEGER, INTENT(IN)                      :: ispin
    LOGICAL, INTENT(IN), OPTIONAL            :: check_moconv_only

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_refinement', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rmone = -1.0_dp, &
                                                rone = 1.0_dp, rzero = 0.0_dp

    INTEGER :: hand1, hand2, hand3, hand4, hand5, handle, ib, ik, imo, &
      imo_low, imo_up, istat, it, jt, nblock, ndim, nmo, nmo_converged, &
      nmo_nc, nmob, num_blocks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: itaken
    LOGICAL                                  :: failure, my_check_moconv_only
    REAL(dp)                                 :: max_norm, min_norm, vmax
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: q_mat, tblock, tvblock
    REAL(dp), DIMENSION(:), POINTER          :: c_res, t_eval
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: v_mat
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: a_mat, b2_mat, b_mat, c2_tmp, &
                                                c3_tmp, c_tmp, chc, evec, hc, &
                                                t_mat, t_vec

    CALL timeset(routineN,handle)
    failure = .FALSE.

    NULLIFY(fm_struct_tmp)
    NULLIFY(chc, evec)
    NULLIFY(c_res, t_eval)
    NULLIFY( hc, t_mat, t_vec, c2_tmp)
    NULLIFY(a_mat, b_mat, b2_mat, v_mat)

    nmo = SIZE(eval,1)
    my_check_moconv_only = .FALSE.
    IF(PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only

    chc => krylov_space%chc_mat(ispin)%matrix
    evec => krylov_space%c_vec(ispin)%matrix
    c_res => krylov_space%c_eval
    t_eval => krylov_space%t_eval

    NULLIFY(fm_struct_tmp, c_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmo,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context)
    CALL cp_fm_create(c_tmp,matrix_struct=fm_struct_tmp,&
                name="c_tmp")
    CALL cp_fm_create(hc,matrix_struct=fm_struct_tmp,&
                name="hc")
    CALL cp_fm_struct_release(fm_struct_tmp)

    !Compute (C^t)HC
    CALL cp_gemm('N','N',nao,nmo,nao,rone,ks,c0,rzero,hc)
    CALL cp_gemm('T','N',nmo,nmo,nao,rone,c0,hc,rzero,chc)

    !Diagonalize  (C^t)HC
    CALL timeset(routineN//"diag_chc",hand1)
    CALL choose_eigv_solver(chc,evec,eval)
    CALL timestop(hand1)

    !Rotate the C vectors
    CALL cp_gemm('N','N',nao,nmo,nmo,rone,c0,evec,rzero,c1)

    !Check for converged states
    CALL cp_gemm('N','N',nao,nmo,nmo,rone,hc,evec,rzero,c_tmp)
    CALL cp_fm_to_fm(c1,c0,nmo,1,1)
    CALL cp_fm_column_scale(c1,eval)
    CALL cp_fm_scale_and_add(1.0_dp,c_tmp,rmone,c1)
    CALL cp_fm_vectorsnorm(c_tmp,c_res)

    nmo_converged = 0
    nmo_nc = 0
    max_norm = 0.0_dp
    min_norm = 1.e10_dp
    CALL cp_fm_set_all(c1,rzero)
    DO imo = 1,nmo
      max_norm = MAX(max_norm,c_res(imo))
      min_norm = MIN(min_norm,c_res(imo))
    END DO
    DO  imo = 1,nmo
      IF(c_res(imo) <= eps_iter ) THEN
          nmo_converged = nmo_converged + 1
      ELSE
          nmo_nc = nmo -nmo_converged
          EXIT
      END IF
    END DO

    nblock = krylov_space%nblock
    num_blocks = nmo_nc/nblock

    krylov_space%nmo_nc = nmo_nc
    krylov_space%nmo_conv = nmo_converged
    krylov_space%max_res_norm = max_norm
    krylov_space%min_res_norm = min_norm

    IF(my_check_moconv_only) THEN
      CALL cp_fm_release(c_tmp)
      CALL cp_fm_release(hc)
      CALL timestop(handle)
      RETURN
    ELSE IF(krylov_space%nmo_nc>0) THEN

      CALL cp_fm_to_fm(c0,c1,nmo_nc,nmo_converged+1,1)

      nblock = krylov_space%nblock
      IF(MODULO(nmo_nc,nblock)>0.0_dp) THEN
        num_blocks = nmo_nc/nblock +1
      ELSE
        num_blocks = nmo_nc/nblock
      END IF

      DO ib = 1, num_blocks

        imo_low = (ib-1)*nblock + 1
        imo_up =  MIN(ib*nblock,nmo_nc)
        nmob = imo_up - imo_low + 1
        ndim = krylov_space%nkrylov*nmob

        NULLIFY(fm_struct_tmp)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=ndim,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context)
        CALL cp_fm_create(c2_tmp,matrix_struct=fm_struct_tmp,&
                      name="c2_tmp")
        CALL cp_fm_struct_release(fm_struct_tmp)
        NULLIFY(fm_struct_tmp)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmob, ncol_global=ndim,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context)
        CALL cp_fm_create(c3_tmp,matrix_struct=fm_struct_tmp,&
                      name="c3_tmp")
        CALL cp_fm_struct_release(fm_struct_tmp)

        ! Create local matrix of right size
        IF(nmob/=nblock) THEN
          NULLIFY(a_mat, b_mat, b2_mat, t_mat, t_vec, v_mat )
          NULLIFY(fm_struct_tmp)
          CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmob, ncol_global=nmob,&
                             para_env=chc%matrix_struct%para_env, &
                             context=chc%matrix_struct%context)
          CALL cp_fm_create(a_mat,matrix_struct=fm_struct_tmp,&
                name="a_mat")
          CALL cp_fm_create(b_mat,matrix_struct=fm_struct_tmp,&
                name="b_mat")
          CALL cp_fm_create(b2_mat,matrix_struct=fm_struct_tmp,&
                name="b2_mat")
          CALL cp_fm_struct_release(fm_struct_tmp)
          NULLIFY(fm_struct_tmp)
          CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ndim, ncol_global=ndim,&
                             para_env=chc%matrix_struct%para_env, &
                             context=chc%matrix_struct%context)
          CALL cp_fm_create(t_mat,matrix_struct=fm_struct_tmp,&
                name="t_mat")
          CALL cp_fm_create(t_vec,matrix_struct=fm_struct_tmp,&
                name="t_vec")
          CALL cp_fm_struct_release(fm_struct_tmp)
          NULLIFY(fm_struct_tmp)
          CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmob,&
                             para_env=c0%matrix_struct%para_env, &
                             context=c0%matrix_struct%context)
          ALLOCATE(v_mat(krylov_space%nkrylov))
          DO ik=1,krylov_space%nkrylov
            CALL cp_fm_create(v_mat(ik)%matrix,matrix_struct=fm_struct_tmp,&
                name="v_mat")
          END DO
          CALL cp_fm_struct_release(fm_struct_tmp)
        ELSE
           a_mat => krylov_space%block1_mat
           b_mat => krylov_space%block2_mat
           b2_mat => krylov_space%block3_mat
           t_mat => krylov_space%block4_mat
           t_vec => krylov_space%block5_mat
           v_mat => krylov_space%v_mat
        END IF

        ALLOCATE(tblock(nmob,nmob),STAT=istat)
        ALLOCATE(tvblock(nmob,ndim),STAT=istat)

    CALL timeset(routineN//"_kry_loop",hand2)
        CALL cp_fm_set_all(b_mat,rzero)
        CALL cp_fm_set_all(t_mat,rzero)
        CALL cp_fm_to_fm(c1,v_mat(1)%matrix,nmob,imo_low,1)

      !Compute A =(V^t)HV
        CALL cp_gemm('N','N',nao,nmob,nao,rone,ks,v_mat(1)%matrix,rzero,hc)
        CALL cp_gemm('T','N',nmob,nmob,nao,rone,v_mat(1)%matrix,hc,&
                   rzero,a_mat)

      !Compute the residual matrix R for next
      !factorisation
        CALL cp_gemm('N','N',nao,nmob,nmob,rone,v_mat(1)%matrix,a_mat,&
                  rzero,c_tmp)
        CALL cp_fm_scale_and_add(rmone,c_tmp,rone,hc)

       ! Build the block tridiagonal matrix
        CALL cp_fm_get_submatrix(a_mat,tblock,1,1,nmob,nmob)
        CALL cp_fm_set_submatrix(t_mat,tblock,1,1,nmob,nmob)

        DO ik = 2,krylov_space%nkrylov

       ! Call lapack for QR factorization
           CALL cp_fm_set_all(b_mat,rzero)
           CALL cp_fm_to_fm(c_tmp,v_mat(ik)%matrix,nmob,1,1)
           CALL cp_fm_qr_factorization(c_tmp, b_mat, nao, nmob, 1,1)

           CALL cp_fm_triangular_multiply(b_mat,  v_mat(ik)%matrix, side="R", invert_tr=.TRUE.,&
              n_rows=nao,n_cols=nmob)

       !Compute A =(V^t)HV
           CALL cp_gemm('N','N',nao,nmob,nao,rone,ks,v_mat(ik)%matrix,rzero,hc)
           CALL cp_gemm('T','N',nmob,nmob,nao,rone,v_mat(ik)%matrix,hc,rzero, a_mat)

       !Compute the !residual matrix R !for next !factorisation
           CALL cp_gemm('N','N',nao,nmob,nmob,rone,v_mat(ik)%matrix,a_mat,&
                   rzero,c_tmp)
           CALL cp_fm_scale_and_add(rmone,c_tmp,rone,hc)
           CALL cp_fm_to_fm(v_mat(ik-1)%matrix,hc,nmob,1,1)
           CALL cp_fm_triangular_multiply(b_mat,hc,side='R',transpose_tr=.TRUE.,&
                   n_rows=nao,n_cols=nmob,alpha=rmone)
           CALL cp_fm_scale_and_add(rone,c_tmp,rone,hc)


        ! Build the block tridiagonal matrix
           it = (ik-2)*nmob + 1
           jt = (ik-1)*nmob + 1

           CALL cp_fm_get_submatrix(a_mat,tblock,1,1,nmob,nmob)
           CALL cp_fm_set_submatrix(t_mat,tblock,jt,jt,nmob,nmob)
           CALL cp_fm_transpose(b_mat,a_mat)
           CALL cp_fm_get_submatrix(a_mat,tblock,1,1,nmob,nmob)
           CALL cp_fm_set_submatrix(t_mat,tblock,it,jt,nmob,nmob)

        END DO ! ik
    CALL timestop(hand2)

        DEALLOCATE(tblock, STAT=istat)

    CALL timeset(routineN//"_diag_tri",hand3)

        CALL choose_eigv_solver(t_mat,t_vec,t_eval)
        ! Diagonalize the block-tridiagonal matrix
    CALL timestop(hand3)

    CALL timeset(routineN//"_build_cnew",hand4)
!        !Compute the refined vectors
        CALL cp_fm_set_all(c2_tmp,rzero)
        DO ik = 1,krylov_space%nkrylov
          jt = (ik-1)*nmob
          CALL cp_gemm('N','N',nao,ndim,nmob,rone,v_mat(ik)%matrix,t_vec,rone,c2_tmp,&
               b_first_row=(jt+1))
        END DO
        DEALLOCATE(tvblock,STAT=istat)

        CALL cp_fm_set_all(c3_tmp,rzero)
        CALL cp_gemm('T','N',nmob,ndim,nao,rone,v_mat(1)%matrix,c2_tmp,rzero,c3_tmp)

        !Try to avoid linear dependencies
        ALLOCATE(q_mat(nmob,ndim),STAT=istat)
        !get max
        CALL cp_fm_get_submatrix(c3_tmp,q_mat,1,1,nmob,ndim)

        ALLOCATE(itaken(ndim),STAT=istat)
        itaken = 0
        DO it = 1,nmob
           vmax = 0.0_dp
           !select index ik
           DO jt = 1,ndim
              IF(itaken(jt)==0 .AND. ABS(q_mat(it,jt))> vmax) THEN
                vmax = ABS(q_mat(it,jt))
                ik=jt
              END IF
           END DO
           itaken(ik) = 1

          CALL cp_fm_to_fm(c2_tmp,v_mat(1)%matrix,1,ik,it)
        END DO
        DEALLOCATE(itaken,STAT=istat)
        DEALLOCATE(q_mat,STAT=istat)

      !Copy in the converged set to enlarge the converged subspace
        CALL cp_fm_to_fm(v_mat(1)%matrix,c0,nmob,1,(nmo_converged+imo_low))
    CALL timestop(hand4)

        IF(nmob<nblock) THEN
          CALL cp_fm_release(a_mat)
          CALL cp_fm_release(b_mat)
          CALL cp_fm_release(b2_mat)
          CALL cp_fm_release(t_mat)
          CALL cp_fm_release(t_vec)
          DO ik=1,krylov_space%nkrylov
             CALL cp_fm_release(v_mat(ik)%matrix)
          END DO
          DEALLOCATE(v_mat,STAT=istat)

        END IF
        CALL cp_fm_release(c2_tmp)
        CALL cp_fm_release(c3_tmp)
      END DO ! ib

    CALL timeset(routineN//"_ortho",hand5)
      CALL cp_gemm('T','N',nmo,nmo,nao,rone,c0,c0,rzero,chc)

      CALL cp_fm_cholesky_decompose(chc,nmo)
      CALL cp_fm_triangular_multiply(chc,c0,'R',invert_tr=.TRUE.)
    CALL timestop(hand5)

      CALL cp_fm_release(c_tmp)
      CALL cp_fm_release(hc)
    ELSE
      CALL cp_fm_release(c_tmp)
      CALL cp_fm_release(hc)
      CALL timestop(handle)
      RETURN
    END IF

    CALL timestop(handle)

  END SUBROUTINE lanczos_refinement

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! *****************************************************************************
!> \brief ...
!> \param krylov_space ...
!> \param ks ...
!> \param c0 ...
!> \param c1 ...
!> \param eval ...
!> \param nao ...
!> \param eps_iter ...
!> \param ispin ...
!> \param check_moconv_only ...
! *****************************************************************************
 SUBROUTINE lanczos_refinement_2v(krylov_space, ks, c0, c1, eval,  nao,  &
               eps_iter, ispin, check_moconv_only)

    TYPE(krylov_space_type), POINTER         :: krylov_space
    TYPE(cp_fm_type), POINTER                :: ks, c0, c1
    REAL(dp), DIMENSION(:), POINTER          :: eval
    INTEGER, INTENT(IN)                      :: nao
    REAL(dp), INTENT(IN)                     :: eps_iter
    INTEGER, INTENT(IN)                      :: ispin
    LOGICAL, INTENT(IN), OPTIONAL            :: check_moconv_only

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lanczos_refinement_2v', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: rmone = -1.0_dp, &
                                                rone = 1.0_dp, rzero = 0.0_dp

    INTEGER :: hand1, hand2, hand3, hand4, hand5, hand6, handle, i, ia, ib, &
      ik, imo, imo_low, imo_up, info, it, j, jt, liwork, lwork, nblock, ndim, &
      nmo, nmo_converged, nmo_nc, nmob, num_blocks
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: itaken
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure, my_check_moconv_only
    REAL(dp)                                 :: max_norm, min_norm, vmax
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a_block, b_block, q_mat, t_mat
    REAL(dp), DIMENSION(:), POINTER          :: c_res, t_eval
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a_loc, b_loc
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: v_mat
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: b_mat, c2_tmp, c_tmp, chc, &
                                                evec, hc, v_tmp

    CALL timeset(routineN,handle)
    failure = .FALSE.

    NULLIFY(fm_struct_tmp)
    NULLIFY(chc, evec)
    NULLIFY(c_res, t_eval)
    NULLIFY(hc, c_tmp, c2_tmp)
    NULLIFY(v_tmp)
    NULLIFY(b_mat, v_mat)
    NULLIFY(b_loc, a_loc)


    nmo = SIZE(eval,1)
    my_check_moconv_only = .FALSE.
    IF(PRESENT(check_moconv_only)) my_check_moconv_only = check_moconv_only

    chc => krylov_space%chc_mat(ispin)%matrix
    evec => krylov_space%c_vec(ispin)%matrix
    c_res => krylov_space%c_eval
    t_eval => krylov_space%t_eval

    NULLIFY(fm_struct_tmp, c_tmp)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmo,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context)
    CALL cp_fm_create(c_tmp,matrix_struct=fm_struct_tmp,&
                name="c_tmp")
    CALL cp_fm_create(hc,matrix_struct=fm_struct_tmp,&
                name="hc")
    CALL cp_fm_struct_release(fm_struct_tmp)

    !Compute (C^t)HC
    CALL cp_gemm('N','N',nao,nmo,nao,rone,ks,c0,rzero,hc)
    CALL cp_gemm('T','N',nmo,nmo,nao,rone,c0,hc,rzero,chc)

    !Diagonalize  (C^t)HC
    CALL timeset(routineN//"diag_chc",hand1)
    CALL choose_eigv_solver(chc,evec,eval)

    CALL timestop(hand1)

    CALL timeset(routineN//"check_conv",hand6)
    !Rotate the C vectors
    CALL cp_gemm('N','N',nao,nmo,nmo,rone,c0,evec,rzero,c1)

    !Check for converged states
    CALL cp_gemm('N','N',nao,nmo,nmo,rone,hc,evec,rzero,c_tmp)
    CALL cp_fm_to_fm(c1,c0,nmo,1,1)
    CALL cp_fm_column_scale(c1,eval)
    CALL cp_fm_scale_and_add(1.0_dp,c_tmp,rmone,c1)
    CALL cp_fm_vectorsnorm(c_tmp,c_res)

    nmo_converged = 0
    nmo_nc = 0
    max_norm = 0.0_dp
    min_norm = 1.e10_dp
    CALL cp_fm_set_all(c1,rzero)
    DO imo = 1,nmo
      max_norm = MAX(max_norm,c_res(imo))
      min_norm = MIN(min_norm,c_res(imo))
    END DO
    DO  imo = 1,nmo
      IF(c_res(imo) <= eps_iter ) THEN
          nmo_converged = nmo_converged + 1
      ELSE
          nmo_nc = nmo -nmo_converged
          EXIT
      END IF
    END DO
    CALL timestop(hand6)

    CALL cp_fm_release(c_tmp)
    CALL cp_fm_release(hc)

    krylov_space%nmo_nc = nmo_nc
    krylov_space%nmo_conv = nmo_converged
    krylov_space%max_res_norm = max_norm
    krylov_space%min_res_norm = min_norm

    IF(my_check_moconv_only) THEN
      ! Do nothing
    ELSE IF(krylov_space%nmo_nc>0) THEN

      CALL cp_fm_to_fm(c0,c1,nmo_nc,nmo_converged+1,1)

      nblock = krylov_space%nblock
      IF(MODULO(nmo_nc,nblock)>0.0_dp) THEN
        num_blocks = nmo_nc/nblock +1
      ELSE
        num_blocks = nmo_nc/nblock
      END IF

      DO ib = 1, num_blocks

        imo_low = (ib-1)*nblock + 1
        imo_up =  MIN(ib*nblock,nmo_nc)
        nmob = imo_up - imo_low + 1
        ndim = krylov_space%nkrylov*nmob

       ! Allocation
        CALL timeset(routineN//"alloc",hand6)
        ALLOCATE(a_block(nmob,nmob))
        ALLOCATE(b_block(nmob,nmob))
        ALLOCATE(t_mat(ndim,ndim))

        NULLIFY(fm_struct_tmp)
      ! by forcing ncol_block=nmo, the needed part of the matrix is distributed on a subset of processes
      ! this is due to the use of two-dimensional grids of processes
      ! nrow_global is distributed over num_pe(1)
      ! a local_data array is anyway allocated for the processes non included
      ! this should have a minimum size
      ! with ncol_local=1.
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmob,&
                       ncol_block=nmob,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context,&
                       force_block=.TRUE.)
        CALL cp_fm_create(c_tmp,matrix_struct=fm_struct_tmp,&
                name="c_tmp")
        CALL cp_fm_set_all(c_tmp,rzero)
        CALL cp_fm_create(v_tmp,matrix_struct=fm_struct_tmp,&
                name="v_tmp")
        CALL cp_fm_set_all(v_tmp,rzero)
        CALL cp_fm_struct_release(fm_struct_tmp)
        NULLIFY(fm_struct_tmp)
        ALLOCATE(v_mat(krylov_space%nkrylov))
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=nmob,&
                       ncol_block=nmob,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context, &
                       force_block=.TRUE.)
        DO ik=1,krylov_space%nkrylov
            CALL cp_fm_create(v_mat(ik)%matrix,matrix_struct=fm_struct_tmp,&
                name="v_mat")
        END DO
        CALL cp_fm_create(hc,matrix_struct=fm_struct_tmp,&
                name="hc")
        CALL cp_fm_struct_release(fm_struct_tmp)
        NULLIFY(fm_struct_tmp)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nao, ncol_global=ndim,&
                       ncol_block=ndim,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context,&
                       force_block=.TRUE.)
        CALL cp_fm_create(c2_tmp,matrix_struct=fm_struct_tmp,&
                      name="c2_tmp")
        CALL cp_fm_struct_release(fm_struct_tmp)

        NULLIFY(fm_struct_tmp)
        CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nmob, ncol_global=nmob,&
                       para_env=c0%matrix_struct%para_env, &
                       context=c0%matrix_struct%context)
        CALL cp_fm_create(b_mat,matrix_struct=fm_struct_tmp,&
                name="b_mat")
        CALL cp_fm_struct_release(fm_struct_tmp)
        CALL timestop(hand6)
        !End allocation

        CALL cp_fm_set_all(b_mat,rzero)
        CALL cp_fm_to_fm(c1,v_mat(1)%matrix,nmob,imo_low,1)

       ! Here starts the construction of krylov space
        CALL timeset(routineN//"_kry_loop",hand2)
        !Compute A =(V^t)HV
        CALL cp_gemm('N','N',nao,nmob,nao,rone,ks,v_mat(1)%matrix,rzero,hc)

        a_block = 0.0_dp
        a_loc => v_mat(1)%matrix%local_data
        b_loc => hc%local_data

        IF(SIZE(hc%local_data,2)==nmob) THEN
          ! this is a work around to avoid problems due to the two dimensional grid of processes
          CALL dgemm('T','N',nmob,nmob,SIZE(hc%local_data,1),1.0_dp,a_loc(1,1),&
               SIZE(hc%local_data,1),b_loc(1,1),SIZE(hc%local_data,1),0.0_dp,a_block(1,1),nmob)
        END IF
        CALL mp_sum(a_block, hc%matrix_struct%para_env%group)

      !Compute the residual matrix R for next
      !factorisation
        c_tmp%local_data=0.0_dp
        IF(SIZE(c_tmp%local_data,2)==nmob) THEN
          a_loc => v_mat(1)%matrix%local_data
          b_loc => c_tmp%local_data
          CALL dgemm('N','N',SIZE(c_tmp%local_data,1),nmob,nmob,1.0_dp,a_loc(1,1),&
               SIZE(c_tmp%local_data,1),a_block(1,1),nmob,0.0_dp,&
               b_loc(1,1), SIZE(c_tmp%local_data,1))
        END IF
        CALL cp_fm_scale_and_add(rmone,c_tmp,rone,hc)

       ! Build the block tridiagonal matrix
        t_mat = 0.0_dp
        DO i=1,nmob
          t_mat(1:nmob,i) = a_block(1:nmob,i)
        END DO

        DO ik = 2,krylov_space%nkrylov
       ! Call lapack for QR factorization
           CALL cp_fm_set_all(b_mat,rzero)
           CALL cp_fm_to_fm(c_tmp,v_mat(ik)%matrix,nmob,1,1)
           CALL cp_fm_qr_factorization(c_tmp, b_mat, nao, nmob, 1,1)

           CALL cp_fm_triangular_multiply(b_mat,  v_mat(ik)%matrix, side="R", invert_tr=.TRUE.,&
              n_rows=nao,n_cols=nmob)

           CALL cp_fm_get_submatrix(b_mat,b_block,1,1,nmob,nmob)

       !Compute A =(V^t)HV
           CALL cp_gemm('N','N',nao,nmob,nao,rone,ks,v_mat(ik)%matrix,rzero,hc)

           a_block = 0.0_dp
           IF(SIZE(hc%local_data,2)==nmob) THEN
             a_loc => v_mat(ik)%matrix%local_data
             b_loc => hc%local_data
             CALL dgemm('T','N',nmob,nmob,SIZE(hc%local_data,1),1.0_dp,a_loc(1,1),&
                 SIZE(hc%local_data,1),b_loc(1,1),SIZE(hc%local_data,1),0.0_dp,a_block(1,1),nmob)
           END IF
           CALL mp_sum(a_block, hc%matrix_struct%para_env%group)

         !Compute the residual matrix R for next
         !factorisation
           c_tmp%local_data=0.0_dp
           IF(SIZE(c_tmp%local_data,2)==nmob) THEN
             a_loc => v_mat(ik)%matrix%local_data
             b_loc => c_tmp%local_data
             CALL dgemm('N','N',SIZE(c_tmp%local_data,1),nmob,nmob,1.0_dp,a_loc(1,1),&
                 SIZE(c_tmp%local_data,1),a_block(1,1),nmob,0.0_dp,&
                 b_loc(1,1), SIZE(c_tmp%local_data,1))
             END IF
             CALL cp_fm_scale_and_add(rmone,c_tmp,rone,hc)

           IF(SIZE(c_tmp%local_data,2)==nmob) THEN
             a_loc => v_mat(ik-1)%matrix%local_data
             DO j = 1,nmob
               DO i = 1,j
                 DO ia = 1,SIZE(c_tmp%local_data,1)
                    b_loc(ia,i) = b_loc(ia,i) - a_loc(ia,j)*b_block(i,j)
                 END DO
               END DO
             END DO
           END IF

        ! Build the block tridiagonal matrix
           it = (ik-2)*nmob
           jt = (ik-1)*nmob
           DO j=1,nmob
             t_mat(jt+1:jt+nmob,jt+j) = a_block(1:nmob,j)
             DO i = 1,nmob
               t_mat(it+i,jt+j) = b_block(j,i)
               t_mat(jt+j,it+i) = b_block(j,i)
             END DO
           END DO
        END DO ! ik
        CALL timestop(hand2)

        CALL timeset(routineN//"_diag_tri",hand3)
        lwork = 1 + 6*ndim + 2*ndim**2 + 5000
        liwork = 5*ndim + 3
        ALLOCATE(work(lwork))
        ALLOCATE(iwork(liwork))

        ! Diagonalize the block-tridiagonal matrix
        CALL dsyevd('V', 'U', ndim, t_mat(1,1),ndim,t_eval(1), &
             work(1), lwork, iwork(1), liwork, info)
        DEALLOCATE(work)
        DEALLOCATE(iwork)
        CALL timestop(hand3)

        CALL timeset(routineN//"_build_cnew",hand4)
!        !Compute the refined vectors

        c2_tmp%local_data=0.0_dp
        ALLOCATE(q_mat(nmob,ndim))
        q_mat = 0.0_dp
        b_loc => c2_tmp%local_data
        DO ik = 1,krylov_space%nkrylov
          CALL cp_fm_to_fm(v_mat(ik)%matrix,v_tmp,nmob,1,1)
          IF(SIZE(c2_tmp%local_data,2)==ndim) THEN
!            a_loc => v_mat(ik)%matrix%local_data
            a_loc => v_tmp%local_data
            it = (ik-1)*nmob
            CALL dgemm('N','N',SIZE(c2_tmp%local_data,1),ndim,nmob,1.0_dp,a_loc(1,1),&
                 SIZE(c2_tmp%local_data,1),t_mat(it+1,1),ndim,1.0_dp,&
                 b_loc(1,1),SIZE(c2_tmp%local_data,1))
          END IF
        END DO !ik

        !Try to avoid linear dependencies
          CALL cp_fm_to_fm(v_mat(1)%matrix,v_tmp,nmob,1,1)
        IF(SIZE(c2_tmp%local_data,2)==ndim) THEN
!          a_loc => v_mat(1)%matrix%local_data
          a_loc => v_tmp%local_data
          b_loc => c2_tmp%local_data
          CALL dgemm('T','N',nmob,ndim,SIZE(v_tmp%local_data,1),1.0_dp,a_loc(1,1),&
                 SIZE(v_tmp%local_data,1),b_loc(1,1),SIZE(v_tmp%local_data,1),&
                 0.0_dp,q_mat(1,1),nmob)
        END IF
        CALL mp_sum(q_mat, hc%matrix_struct%para_env%group)

        ALLOCATE(itaken(ndim))
        itaken = 0
        DO it = 1,nmob
           vmax = 0.0_dp
           !select index ik
           DO jt = 1,ndim
              IF(itaken(jt)==0 .AND. ABS(q_mat(it,jt))> vmax) THEN
                vmax = ABS(q_mat(it,jt))
                ik=jt
              END IF
           END DO
           itaken(ik) = 1

          CALL cp_fm_to_fm(c2_tmp,v_mat(1)%matrix,1,ik,it)
        END DO
        DEALLOCATE(itaken)
        DEALLOCATE(q_mat)

       !Copy in the converged set to enlarge the converged subspace
        CALL cp_fm_to_fm(v_mat(1)%matrix,c0,nmob,1,(nmo_converged+imo_low))
        CALL timestop(hand4)


        CALL cp_fm_release(c2_tmp)
        CALL cp_fm_release(c_tmp)
        CALL cp_fm_release(hc)
        CALL cp_fm_release(v_tmp)
        CALL cp_fm_release(b_mat)

        DEALLOCATE(t_mat)
        DEALLOCATE(a_block)
        DEALLOCATE(b_block)
        DO ik=1,krylov_space%nkrylov
           CALL cp_fm_release(v_mat(ik)%matrix)
        END DO
        DEALLOCATE(v_mat)

      END DO ! ib

      CALL timeset(routineN//"_ortho",hand5)
      CALL cp_gemm('T','N',nmo,nmo,nao,rone,c0,c0,rzero,chc)

      CALL cp_fm_cholesky_decompose(chc,nmo)
      CALL cp_fm_triangular_multiply(chc,c0,'R',invert_tr=.TRUE.)
      CALL timestop(hand5)
    ELSE
      ! Do nothing
    END IF

    CALL timestop(handle)
  END SUBROUTINE lanczos_refinement_2v

END MODULE qs_scf_lanczos
