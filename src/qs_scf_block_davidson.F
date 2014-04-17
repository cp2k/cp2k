!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief module that contains the algorithms to perform an itrative
!>         diagonalization by the block-Davidson approach
!>         P. Blaha, et al J. Comp. Physics, 229, (2010), 453-460
!>         Iterative diagonalization in augmented plane wave based
!>         methods in electronic structure calculations
!> \par History
!>      05.2011 created [MI]
!> \author MI
! *****************************************************************************
MODULE qs_scf_block_davidson

  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_release
  USE cp_dbcsr_interface,              ONLY: &
       cp_create_bl_distribution, cp_dbcsr_add, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_distribution, &
       cp_dbcsr_distribution_release, cp_dbcsr_get_diag, cp_dbcsr_get_info, &
       cp_dbcsr_init_p, cp_dbcsr_multiply, cp_dbcsr_norm, cp_dbcsr_release_p, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale_by_vector
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_transpose,&
                                             cp_fm_triangular_invert,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                             cp_fm_power
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_release, &
       cp_fm_set_all, cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type, &
       cp_fm_vectorsnorm
  USE dbcsr_api,                       ONLY: &
       dbcsr_distribution_mp, dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_npcols, dbcsr_mp_nprows, &
       dbcsr_norm_column, dbcsr_type_no_symmetry, dbcsr_type_real_default, &
       dbcsr_type_symmetric
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_sum
  USE preconditioner,                  ONLY: apply_preconditioner
  USE preconditioner_types,            ONLY: preconditioner_type
  USE qs_block_davidson_types,         ONLY: davidson_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_block_davidson'

  PUBLIC :: generate_extended_space, generate_extended_space_sparse

CONTAINS

! *****************************************************************************

  SUBROUTINE ritz_coefficients(bdav_env,mo_coeff,matrix_sc,matrix_hc,ritz_coeff,error)

    TYPE(davidson_type)                      :: bdav_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff, matrix_sc, matrix_hc
    REAL(dp), DIMENSION(:)                   :: ritz_coeff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ritz_coefficients', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, istat, nao, nmo
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: chc_diag, csc_diag
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: block_mat, matrix_tmp

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(block_mat,fm_struct_tmp,matrix_tmp)
    CALL cp_fm_get_info(mo_coeff,nrow_global=nao,ncol_global=nmo,error=error)

    ALLOCATE(csc_diag(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(chc_diag(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! storage matrix of size mos x mos, only the diagonal elements are used
    CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmo,ncol_global=nmo, &
                             context=mo_coeff%matrix_struct%context, &
                             para_env=mo_coeff%matrix_struct%para_env,error=error)
    CALL cp_fm_create(matrix_tmp,fm_struct_tmp,name="matrix_tmp",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    ! since we only use diagonal elements this is a bit of a waste
    ! compute CSC
!   CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,mo_coeff,matrix_sc,0.0_dp,matrix_tmp,error=error)
!    CALL cp_fm_get_diag(matrix_tmp,csc_diag,error=error)
    ! set the top left part of S[C,Z] block matrix  CSC
    block_mat => bdav_env%S_block_mat
    CALL cp_fm_set_all(block_mat,0.0_dp, 1.0_dp,error=error)
!    CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1,1,error=error)

    ! compute CHC
    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,mo_coeff,matrix_hc,0.0_dp,matrix_tmp,error=error)
    CALL cp_fm_get_diag(matrix_tmp,chc_diag,error=error)
    ! set the top left part of H[C,Z] block matrix CHC
    block_mat => bdav_env%H_block_mat
    CALL cp_fm_to_fm_submat(matrix_tmp,block_mat,nmo,nmo,1,1,1,1,error=error)

    DO i=1,nmo
!      IF(ABS(csc_diag(i))>EPSILON(0.0_dp)) THEN
        ritz_coeff(i) = chc_diag(i)!/csc_diag(i)
!      END IF
    END DO
    CALL cp_fm_release(matrix_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE ritz_coefficients
! +++++++++++++++++++++++++++++

  SUBROUTINE generate_extended_space(bdav_env,mo_set,matrix_h,matrix_s,output_unit,&
             preconditioner,error)

    TYPE(davidson_type)                      :: bdav_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(preconditioner_type), OPTIONAL, &
      POINTER                                :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_extended_space', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i_first, i_last, imo, istat, iter, j, jj, &
      max_iter, n, nao, nmat, nmat2, nmo, nmo_converged, nmo_not_converged, &
      nset, nset_conv, nset_not_conv
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iconv, inotconv
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: iconv_set, inotconv_set
    LOGICAL                                  :: converged, &
                                                do_apply_preconditioner, &
                                                failure
    REAL(dp)                                 :: lambda, max_norm, min_norm, &
                                                t1, t2
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: ritz_coeff, vnorm
    REAL(dp), DIMENSION(:), POINTER          :: eig_not_conv, eigenvalues, &
                                                evals
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_b
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER :: c_conv, c_notconv, c_out, c_pz, c_z, &
      h_block, h_fm, m_hc, m_sc, m_tmp, mo_coeff, mt_tmp, s_block, s_fm, &
      v_block, w_block

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(mo_coeff,mo_coeff_b,eigenvalues)

    do_apply_preconditioner = .FALSE.
    IF(PRESENT(preconditioner)) do_apply_preconditioner=.TRUE.
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff,mo_coeff_b=mo_coeff_b,eigenvalues=eigenvalues,&
                    nao=nao,nmo=nmo, homo=homo)
    IF(do_apply_preconditioner) THEN
         max_iter =  bdav_env%max_iter
    ELSE
       max_iter = 1
    END IF

    NULLIFY(c_conv, c_notconv, c_out, c_z, c_pz, m_hc, m_sc, m_tmp, mt_tmp)
    NULLIFY(h_block, h_fm, s_block, s_fm, v_block, w_block)
    NULLIFY(evals,eig_not_conv )
    t1 = m_walltime()
    IF (output_unit > 0) THEN
          WRITE(output_unit,"(T15,A,T23,A,T36,A,T49,A,T60,A,/,T8,A)")   &
                " Cycle ", " conv. MOS ", " B2MAX ", " B2MIN ", " Time",  REPEAT("-",60)
    END IF

    ALLOCATE(iconv(nmo), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(inotconv(nmo), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ritz_coeff(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vnorm(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    converged=.FALSE.
    DO iter = 1, max_iter

     ! compute Ritz values
      ritz_coeff=0.0_dp
      CALL cp_fm_create(m_hc,mo_coeff%matrix_struct,name="hc",error=error)
      CALL cp_dbcsr_sm_fm_multiply(matrix_h,mo_coeff,m_hc,nmo,error=error)
      CALL cp_fm_create(m_sc,mo_coeff%matrix_struct,name="sc",error=error)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s,mo_coeff,m_sc,nmo,error=error)

      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmo,ncol_global=nmo, &
                             context=mo_coeff%matrix_struct%context, &
                             para_env=mo_coeff%matrix_struct%para_env,error=error)
      CALL cp_fm_create(m_tmp,fm_struct_tmp,name="matrix_tmp",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)

      CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,mo_coeff,m_hc,0.0_dp,m_tmp,error=error)
      CALL cp_fm_get_diag(m_tmp,ritz_coeff,error=error)
      CALL cp_fm_release(m_tmp,error=error)

      ! Check for converged eigenvectors
!      CALL cp_fm_create(c_z,mo_coeff%matrix_struct,name="tmp",error=error)
      c_z => bdav_env%matrix_z
      c_pz => bdav_env%matrix_pz
      CALL cp_fm_to_fm(m_sc,c_z,error=error)
      CALL cp_fm_column_scale(c_z,ritz_coeff)
      CALL cp_fm_scale_and_add(-1.0_dp,c_z,1.0_dp,m_hc,error=error)
      CALL cp_fm_vectorsnorm(c_z,vnorm,error=error)

      nmo_converged = 0
      nmo_not_converged = 0
      max_norm = 0.0_dp
      min_norm = 1.e10_dp
      DO imo = 1,nmo
        max_norm = MAX(max_norm,vnorm(imo))
        min_norm = MIN(min_norm,vnorm(imo))
      END DO
      iconv = 0
      inotconv = 0
      DO  imo = 1,nmo
        IF(vnorm(imo) <= bdav_env%eps_iter ) THEN
           nmo_converged = nmo_converged + 1
           iconv(nmo_converged)=imo
        ELSE
           nmo_not_converged = nmo_not_converged + 1
           inotconv(nmo_not_converged)=imo
        END IF
      END DO

      IF(nmo_converged>0) THEN
        ALLOCATE(iconv_set(nmo_converged,2), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(inotconv_set(nmo_not_converged,2), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        i_last = iconv(1)
        nset = 0
        DO j = 1,nmo_converged
          imo=iconv(j)

          IF(imo==i_last+1) THEN
           i_last=imo
           iconv_set(nset,2)=imo
          ELSE
           i_last=imo
           nset=nset+1
           iconv_set(nset,1)=imo
           iconv_set(nset,2)=imo
          END IF
        END DO
        nset_conv = nset

        i_last = inotconv(1)
        nset = 0
        DO j = 1,nmo_not_converged
          imo=inotconv(j)

          IF(imo==i_last+1) THEN
           i_last=imo
           inotconv_set(nset,2)=imo
          ELSE
           i_last=imo
           nset=nset+1
           inotconv_set(nset,1)=imo
           inotconv_set(nset,2)=imo
          END IF
        END DO
        nset_not_conv = nset
        CALL cp_fm_release(m_sc,error=error)
        CALL cp_fm_release(m_hc,error=error)
        NULLIFY(c_z, c_pz)
      END IF


      IF(REAL(nmo_converged,dp)/REAL(nmo,dp)>bdav_env%conv_percent) THEN
         converged=.TRUE.
         DEALLOCATE(iconv_set,STAT=istat)
         DEALLOCATE (inotconv_set,STAT=istat)
         t2 = m_walltime()
         IF (output_unit > 0) THEN
           WRITE(output_unit,'(T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')&
                 iter, nmo_converged,  max_norm, min_norm, t2-t1

           WRITE(output_unit,*)  " Reached convergence in ", iter, &
             " Davidson iterations"
         END IF

         EXIT
      END IF

      IF(nmo_converged>0) THEN
        CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nao,ncol_global=nao, &
                                 context=mo_coeff%matrix_struct%context, &
                                 para_env=mo_coeff%matrix_struct%para_env,error=error)
        !allocate h_fm
        CALL cp_fm_create(h_fm,fm_struct_tmp,name="matrix_tmp",error=error)
        !allocate s_fm
        CALL cp_fm_create(s_fm,fm_struct_tmp,name="matrix_tmp",error=error)
        !copy matrix_h in h_fm
        CALL copy_dbcsr_to_fm(matrix_h,h_fm,error=error)
        CALL cp_fm_upper_to_full(h_fm,s_fm,error=error)

        !copy matrix_s in s_fm
!        CALL cp_fm_set_all(s_fm,0.0_dp,error=error)
        CALL copy_dbcsr_to_fm(matrix_s,s_fm,error=error)

        !allocate c_out
        CALL cp_fm_create(c_out,fm_struct_tmp,name="matrix_tmp",error=error)
        ! set c_out to zero
        CALL cp_fm_set_all(c_out,0.0_dp,error=error)
        CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       !allocate c_conv
        CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nao,ncol_global=nmo_converged, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)
        CALL cp_fm_create(c_conv,fm_struct_tmp,name="c_conv",error=error)
        CALL cp_fm_set_all(c_conv,0.0_dp,error=error)
        !allocate m_tmp
        CALL cp_fm_create(m_tmp,fm_struct_tmp,name="m_tmp_nxmc",error=error)
        CALL cp_fm_struct_release(fm_struct_tmp,error=error)
      END IF

      !allocate c_notconv
      CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nao,ncol_global=nmo_not_converged, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)
      CALL cp_fm_create(c_notconv,fm_struct_tmp,name="c_notconv",error=error)
      CALL cp_fm_set_all(c_notconv,0.0_dp,error=error)
      IF(nmo_converged>0) THEN
        CALL cp_fm_create(m_hc,fm_struct_tmp,name="m_hc",error=error)
        CALL cp_fm_create(m_sc,fm_struct_tmp,name="m_sc",error=error)
        !allocate c_z
        CALL cp_fm_create(c_z,fm_struct_tmp,name="c_z",error=error)
        CALL cp_fm_create(c_pz,fm_struct_tmp,name="c_pz",error=error)
        CALL cp_fm_set_all(c_z,0.0_dp,error=error)

      ! sum contributions to c_out
        jj=1
        DO j=1,nset_conv
          i_first=iconv_set(j,1)
          i_last=iconv_set(j,2)
          n=i_last-i_first+1
          CALL cp_fm_to_fm_submat(mo_coeff,c_conv,nao,n,1,i_first,1,jj,error=error)
          jj=jj+n
        END DO
        CALL cp_fm_symm('L','U',nao,nmo_converged,1.0_dp,s_fm,c_conv,0.0_dp,m_tmp,error=error)
        CALL cp_fm_gemm('N','T',nao,nao,nmo_converged,1.0_dp,m_tmp,m_tmp,0.0_dp,c_out,error=error)

        ! project c_out out of H
        lambda = 100.0_dp*ABS(eigenvalues(homo))
        CALL cp_fm_scale_and_add(lambda,c_out,1.0_dp,h_fm,error=error)
        CALL cp_fm_release(m_tmp,error=error)
        CALL cp_fm_release(h_fm,error=error)

      END IF

      !allocate m_tmp
      CALL cp_fm_create(m_tmp,fm_struct_tmp,name="m_tmp_nxm",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)
      IF(nmo_converged>0) THEN
        ALLOCATE(eig_not_conv(nmo_not_converged))
        jj=1
        DO j=1,nset_not_conv
          i_first=inotconv_set(j,1)
          i_last=inotconv_set(j,2)
          n=i_last-i_first+1
          CALL cp_fm_to_fm_submat(mo_coeff,c_notconv,nao,n,1,i_first,1,jj,error=error)
          eig_not_conv(jj:jj+n-1) = ritz_coeff(i_first:i_last)
          jj=jj+n
        END DO
        CALL cp_fm_gemm('N','N',nao,nmo_not_converged,nao,1.0_dp,c_out,c_notconv,0.0_dp,m_hc,error=error)
        CALL cp_fm_symm('L','U',nao,nmo_not_converged,1.0_dp,s_fm,c_notconv,0.0_dp,m_sc,error=error)
        ! extend suspace using only the not converged vectors
         CALL cp_fm_to_fm(m_sc,m_tmp,error=error)
         CALL cp_fm_column_scale(m_tmp,eig_not_conv)
         CALL cp_fm_scale_and_add(-1.0_dp,m_tmp,1.0_dp,m_hc,error=error)
         DEALLOCATE(eig_not_conv)
         CALL cp_fm_to_fm(m_tmp,c_z,error=error)
      ELSE
        CALL cp_fm_to_fm(mo_coeff,c_notconv,error=error)
      END IF

      !preconditioner
       IF(do_apply_preconditioner) THEN
         IF (preconditioner%in_use/=0) THEN
           CALL apply_preconditioner(preconditioner,c_z,c_pz,error=error)
         ELSE
           CALL cp_fm_to_fm(c_z,c_pz,error=error)
         ENDIF
       ELSE
         CALL cp_fm_to_fm(c_z,c_pz,error=error)
       END IF
       CALL cp_fm_release(m_tmp,error=error)

       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmo_not_converged,ncol_global=nmo_not_converged, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)

       CALL cp_fm_create(m_tmp,fm_struct_tmp,name="m_tmp_mxm",error=error)
       CALL cp_fm_create(mt_tmp,fm_struct_tmp,name="mt_tmp_mxm",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       nmat = nmo_not_converged
       nmat2 = 2* nmo_not_converged
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmat2,ncol_global=nmat2, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)

       CALL cp_fm_create(s_block,fm_struct_tmp,name="sb",error=error)
       CALL cp_fm_create(h_block,fm_struct_tmp,name="hb",error=error)
       CALL cp_fm_create(v_block,fm_struct_tmp,name="vb",error=error)
       CALL cp_fm_create(w_block,fm_struct_tmp,name="wb",error=error)
       ALLOCATE(evals(nmat2))

       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       ! compute CSC
       CALL cp_fm_set_all(s_block,0.0_dp, 1.0_dp,error=error)

       ! compute CHC
       CALL cp_fm_gemm('T','N',nmat,nmat,nao,1.0_dp,c_notconv,m_hc,0.0_dp,m_tmp,error=error)
       CALL cp_fm_to_fm_submat(m_tmp,h_block,nmat,nmat,1,1,1,1,error=error)

       ! compute ZSC
       CALL cp_fm_gemm('T','N',nmat,nmat,nao,1.0_dp,c_pz,m_sc,0.0_dp,m_tmp,error=error)
       CALL cp_fm_to_fm_submat(m_tmp,s_block,nmat,nmat,1,1,1+nmat,1,error=error)
       CALL cp_fm_transpose(m_tmp,mt_tmp,error=error)
       CALL cp_fm_to_fm_submat(mt_tmp,s_block,nmat,nmat,1,1,1,1+nmat,error=error)
       ! compute ZHC
       CALL cp_fm_gemm('T','N',nmat,nmat,nao,1.0_dp,c_pz,m_hc,0.0_dp,m_tmp,error=error)
       CALL cp_fm_to_fm_submat(m_tmp,h_block,nmat,nmat,1,1,1+nmat,1,error=error)
       CALL cp_fm_transpose(m_tmp,mt_tmp,error=error)
       CALL cp_fm_to_fm_submat(mt_tmp,h_block,nmat,nmat,1,1,1,1+nmat,error=error)

       CALL cp_fm_release(mt_tmp,error=error)

       ! reuse m_sc and m_hc to computr HZ and SZ
       IF(nmo_converged>0) THEN
         CALL cp_fm_gemm('N','N',nao,nmat,nao,1.0_dp,c_out,c_pz,0.0_dp,m_hc,error=error)
         CALL cp_fm_symm('L','U',nao,nmo_not_converged,1.0_dp,s_fm,c_pz,0.0_dp,m_sc,error=error)

         CALL cp_fm_release(c_out,error=error)
         CALL cp_fm_release(c_conv,error=error)
         CALL cp_fm_release(s_fm,error=error)
       ELSE
         CALL cp_dbcsr_sm_fm_multiply(matrix_h,c_pz,m_hc,nmo,error=error)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s,c_pz,m_sc,nmo,error=error)
       END IF

       ! compute ZSZ
       CALL cp_fm_gemm('T','N',nmat,nmat,nao,1.0_dp,c_pz,m_sc,0.0_dp,m_tmp,error=error)
       CALL cp_fm_to_fm_submat(m_tmp,s_block,nmat,nmat,1,1,1+nmat,1+nmat,error=error)
       ! compute ZHZ
       CALL cp_fm_gemm('T','N',nmat,nmat,nao,1.0_dp,c_pz,m_hc,0.0_dp,m_tmp,error=error)
       CALL cp_fm_to_fm_submat(m_tmp,h_block,nmat,nmat,1,1,1+nmat,1+nmat,error=error)

       CALL cp_fm_release(m_sc,error=error)

       ! solution of the reduced eigenvalues problem
       CALL reduce_extended_space(s_block,h_block,v_block,w_block, evals, nmat2, error=error)

      ! extract egenvectors
       CALL cp_fm_to_fm_submat(v_block,m_tmp, nmat,nmat,1,1,1,1,error=error)
       CALL cp_fm_gemm('N','N',nao,nmat,nmat,1.0_dp,c_notconv,m_tmp,0.0_dp,m_hc,error=error)
       CALL cp_fm_to_fm_submat(v_block,m_tmp, nmat,nmat,1+nmat,1,1,1,error=error)
       CALL cp_fm_gemm('N','N',nao,nmat,nmat,1.0_dp,c_pz,m_tmp,1.0_dp,m_hc,error=error)

       CALL cp_fm_release(m_tmp,error=error)

       CALL cp_fm_release(c_notconv,error=error)
       CALL cp_fm_release(s_block,error=error)
       CALL cp_fm_release(h_block,error=error)
       CALL cp_fm_release(w_block,error=error)
       CALL cp_fm_release(v_block,error=error)

       IF(nmo_converged>0) THEN
         CALL cp_fm_release(c_z,error=error)
         CALL cp_fm_release(c_pz,error=error)
         jj=1
         DO j=1,nset_not_conv
           i_first=inotconv_set(j,1)
           i_last=inotconv_set(j,2)
           n=i_last-i_first+1
           CALL cp_fm_to_fm_submat(m_hc,mo_coeff,nao,n,1,jj,1,i_first,error=error)
           eigenvalues(i_first:i_last) = evals(jj:jj+n-1)
           jj=jj+n
         END DO
         DEALLOCATE(iconv_set,STAT=istat)
         DEALLOCATE (inotconv_set,STAT=istat)
       ELSE
         CALL cp_fm_to_fm(m_hc,mo_coeff,error=error)
         eigenvalues(1:nmo) = evals(1:nmo)
       END IF
       DEALLOCATE(evals,STAT=istat)

       CALL cp_fm_release(m_hc,error=error)

       CALL copy_fm_to_dbcsr(mo_coeff,mo_coeff_b,error=error)!fm->dbcsr

       t2 = m_walltime()
       IF (output_unit > 0) THEN
         WRITE(output_unit,'(T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')&
                 iter, nmo_converged,  max_norm, min_norm, t2-t1
       END IF
       t1=m_walltime()

    END DO  ! iter

    DEALLOCATE(iconv,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(inotconv,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ritz_coeff,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(vnorm,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)
  END SUBROUTINE generate_extended_space

! *****************************************************************************
  SUBROUTINE generate_extended_space_sparse(bdav_env,mo_set,matrix_h,matrix_s,output_unit,&
             preconditioner,error)

    TYPE(davidson_type)                      :: bdav_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h, matrix_s
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(preconditioner_type), OPTIONAL, &
      POINTER                                :: preconditioner
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'generate_extended_space_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i_first, i_last, imo, istat, iter, j, jj, k, &
      max_iter, n, nao, nmat, nmat2, nmo, nmo_converged, nmo_not_converged, &
      nset, nset_conv, nset_not_conv
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iconv, inotconv
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: iconv_set, inotconv_set
    LOGICAL                                  :: converged, &
                                                do_apply_preconditioner, &
                                                failure
    REAL(dp)                                 :: lambda, max_norm, min_norm, &
                                                t1, t2
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eig_not_conv, evals, &
                                                ritz_coeff, vnorm
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist, &
                                                row_blk_size, row_dist
    TYPE(cp_dbcsr_type), POINTER :: c_out, matrix_hc, matrix_mm, matrix_pz, &
      matrix_sc, matrix_z, mo_coeff_b, mo_conv, mo_notconv, smo_conv
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER :: h_block, matrix_mm_fm, matrix_mmt_fm, &
      matrix_nm_fm, matrix_z_fm, mo_coeff, mo_conv_fm, mo_notconv_fm, &
      s_block, v_block, w_block
    TYPE(dbcsr_distribution_obj)             :: dist

    failure=.FALSE.
    CALL timeset(routineN,handle)

    do_apply_preconditioner = .FALSE.
    IF(PRESENT(preconditioner)) do_apply_preconditioner=.TRUE.

    NULLIFY(mo_coeff,mo_coeff_b, matrix_hc, matrix_sc, matrix_z, matrix_pz, matrix_mm)
    NULLIFY(mo_conv_fm, mo_notconv_fm, mo_conv, mo_notconv, smo_conv,c_out)
    NULLIFY(matrix_mm_fm, matrix_mmt_fm, mo_coeff, matrix_nm_fm, matrix_z_fm)
    NULLIFY(h_block,s_block,v_block,w_block)
    NULLIFY(fm_struct_tmp)
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b,&
         eigenvalues=eigenvalues,homo=homo,nao=nao,nmo=nmo)
    IF(do_apply_preconditioner) THEN
         max_iter =  bdav_env%max_iter
    ELSE
       max_iter = 1
    END IF

    t1 = m_walltime()
    IF (output_unit > 0) THEN
          WRITE(output_unit,"(T15,A,T23,A,T36,A,T49,A,T60,A,/,T8,A)")   &
                " Cycle ", " conv. MOS ", " B2MAX ", " B2MIN ", " Time",  REPEAT("-",60)
    END IF


    ! Allocate array for Ritz values
    ALLOCATE(ritz_coeff(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(iconv(nmo), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(inotconv(nmo), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vnorm(nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    converged=.FALSE.
    DO iter = 1, max_iter
      NULLIFY(c_out, mo_conv, mo_notconv_fm, mo_notconv)
      ! Prepare HC and SC, using mo_coeff_b (sparse), these are still sparse
      CALL cp_dbcsr_init_p(matrix_hc,error=error)
      CALL cp_dbcsr_create(matrix_hc,"matrix_hc",cp_dbcsr_distribution(mo_coeff_b),&
         dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(mo_coeff_b),&
         cp_dbcsr_col_block_sizes(mo_coeff_b),0,0,dbcsr_type_real_default,error=error)
      CALL cp_dbcsr_init_p(matrix_sc,error=error)
      CALL cp_dbcsr_create(matrix_sc,"matrix_sc",cp_dbcsr_distribution(mo_coeff_b),&
         dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(mo_coeff_b),&
         cp_dbcsr_col_block_sizes(mo_coeff_b),0,0,dbcsr_type_real_default,error=error)

      CALL cp_dbcsr_get_info(mo_coeff_b,nfullrows_total=n,nfullcols_total=k)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_h,mo_coeff_b,0.0_dp,matrix_hc,last_column=k,error=error)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_s,mo_coeff_b,0.0_dp,matrix_sc,last_column=k,error=error)

      ! compute Ritz values
      ritz_coeff=0.0_dp
      ! Allocate Sparse matrices: nmoxnmo
      ! matrix_mm
      CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo, &
            dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_b))))
      CALL cp_create_bl_distribution (row_dist, row_blk_size, nmo, &
            dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_b))))
      CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp(cp_dbcsr_distribution(mo_coeff_b)),&
            row_dist,col_dist)
      CALL cp_dbcsr_init_p(matrix_mm,error=error)
      CALL cp_dbcsr_create(matrix_mm,"matrix_mm",dist,dbcsr_type_no_symmetry,&
           row_blk_size,col_blk_size,0,0,dbcsr_type_real_default,error=error)
      CALL cp_dbcsr_distribution_release (dist)
      CALL array_release (col_blk_size)
      CALL array_release (col_dist)
      CALL array_release (row_blk_size)
      CALL array_release (row_dist)


      CALL cp_dbcsr_multiply('t','n',1.0_dp,mo_coeff_b,matrix_hc,0.0_dp,matrix_mm, last_column=k, error=error)
      CALL cp_dbcsr_get_diag(matrix_mm,ritz_coeff,error=error)
      CALL mp_sum(ritz_coeff,mo_coeff%matrix_struct%para_env%group)


      ! extended subspace P Z = P [H - theta S]C  this ia another matrix of type and size as mo_coeff_b
      CALL cp_dbcsr_init_p(matrix_z,error=error)
      CALL cp_dbcsr_create(matrix_z,"matrix_z",cp_dbcsr_distribution(mo_coeff_b),&
           dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(mo_coeff_b),&
           cp_dbcsr_col_block_sizes(mo_coeff_b),0,0,dbcsr_type_real_default,error=error)
      CALL cp_dbcsr_copy(matrix_z,matrix_sc,error=error)
      CALL cp_dbcsr_scale_by_vector(matrix_z,ritz_coeff,side='right',error=error)
      CALL cp_dbcsr_add(matrix_z,matrix_hc,-1.0_dp,1.0_dp,error=error)

      ! Check for converged eigenvectors
      vnorm =0.0_dp
      CALL cp_dbcsr_norm(matrix_z,which_norm=dbcsr_norm_column,norm_vector=vnorm,error=error)
      nmo_converged = 0
      nmo_not_converged = 0
      max_norm = 0.0_dp
      min_norm = 1.e10_dp
      DO imo = 1,nmo
        max_norm = MAX(max_norm,vnorm(imo))
        min_norm = MIN(min_norm,vnorm(imo))
      END DO
      iconv = 0
      inotconv = 0

      DO  imo = 1,nmo
        IF(vnorm(imo) <= bdav_env%eps_iter ) THEN
           nmo_converged = nmo_converged + 1
           iconv(nmo_converged)=imo
        ELSE
           nmo_not_converged = nmo_not_converged + 1
           inotconv(nmo_not_converged)=imo
        END IF
      END DO


      IF(nmo_converged>0) THEN
        ALLOCATE(iconv_set(nmo_converged,2), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(inotconv_set(nmo_not_converged,2), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        i_last = iconv(1)
        nset = 0
        DO j = 1,nmo_converged
          imo=iconv(j)

          IF(imo==i_last+1) THEN
           i_last=imo
           iconv_set(nset,2)=imo
          ELSE
           i_last=imo
           nset=nset+1
           iconv_set(nset,1)=imo
           iconv_set(nset,2)=imo
          END IF
        END DO
        nset_conv = nset

        i_last = inotconv(1)
        nset = 0
        DO j = 1,nmo_not_converged
          imo=inotconv(j)

          IF(imo==i_last+1) THEN
           i_last=imo
           inotconv_set(nset,2)=imo
          ELSE
           i_last=imo
           nset=nset+1
           inotconv_set(nset,1)=imo
           inotconv_set(nset,2)=imo
          END IF
        END DO
        nset_not_conv = nset

        CALL cp_dbcsr_release_p(matrix_hc, error=error)
        CALL cp_dbcsr_release_p(matrix_sc, error=error)
        CALL cp_dbcsr_release_p(matrix_z, error=error)
        CALL cp_dbcsr_release_p(matrix_mm, error=error)
      END IF

      IF(REAL(nmo_converged,dp)/REAL(nmo,dp)>bdav_env%conv_percent) THEN
        DEALLOCATE(iconv_set,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        DEALLOCATE (inotconv_set,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        converged=.TRUE.
        t2 = m_walltime()
        IF (output_unit > 0)  THEN
            WRITE(output_unit,'(T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')&
                 iter, nmo_converged,  max_norm, min_norm, t2-t1

            WRITE(output_unit,*)  " Reached convergence in ", iter, &
             " Davidson iterations"
        END IF

        EXIT
      END IF

      IF(nmo_converged>0) THEN

        !allocate mo_conv_fm
        CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nao,ncol_global=nmo_converged, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)
        CALL cp_fm_create(mo_conv_fm,fm_struct_tmp,name="mo_conv_fm",error=error)

        CALL cp_fm_struct_release(fm_struct_tmp,error=error)

        ! extract mo_conv from mo_coeff full matrix
        jj=1
        DO j=1,nset_conv
          i_first=iconv_set(j,1)
          i_last=iconv_set(j,2)
          n=i_last-i_first+1
          CALL cp_fm_to_fm_submat(mo_coeff,mo_conv_fm,nao,n,1,i_first,1,jj,error=error)
          jj=jj+n
        END DO

        ! allocate c_out sparse matrix, to project out the converged MOS
        CALL cp_dbcsr_init_p(c_out,error=error)
        CALL cp_dbcsr_create(c_out,"c_out",cp_dbcsr_distribution(matrix_s),&
             dbcsr_type_symmetric,cp_dbcsr_row_block_sizes(matrix_s),&
             cp_dbcsr_col_block_sizes(matrix_s),0,0,dbcsr_type_real_default,error=error)

        ! allocate mo_conv sparse
        CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo_converged, &
              dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
        CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s)),&
              dbcsr_distribution_row_dist(cp_dbcsr_distribution(matrix_s)),col_dist)
        CALL cp_dbcsr_init_p(mo_conv,error=error)
        CALL cp_dbcsr_create(mo_conv,"mo_conv",dist,dbcsr_type_no_symmetry,&
              cp_dbcsr_row_block_sizes(matrix_s),col_blk_size,0,0,dbcsr_type_real_default,error=error)
        CALL cp_dbcsr_init_p(smo_conv,error=error)
        CALL cp_dbcsr_create(smo_conv,"smo_conv",dist,dbcsr_type_no_symmetry,&
              cp_dbcsr_row_block_sizes(matrix_s),col_blk_size,0,0,dbcsr_type_real_default,error=error)

        CALL copy_fm_to_dbcsr(mo_conv_fm,mo_conv,error=error)!fm->dbcsr
        CALL cp_dbcsr_distribution_release (dist)
        CALL array_release (col_blk_size)
        CALL array_release (col_dist)

        CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_s,mo_conv,0.0_dp,smo_conv,last_column=nmo_converged,error=error)
        CALL cp_dbcsr_multiply('n','t',1.0_dp,smo_conv,smo_conv,0.0_dp,c_out,last_column=nao,error=error)
       ! project c_out out of H
        lambda = 100.0_dp*ABS(eigenvalues(homo))
        CALL cp_dbcsr_add(c_out,matrix_h,lambda,1.0_dp,error=error)

        CALL cp_dbcsr_release_p(mo_conv, error=error)
        CALL cp_dbcsr_release_p(smo_conv, error=error)
        CALL cp_fm_release(mo_conv_fm, error=error)


        !allocate c_notconv_fm
        CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nao,ncol_global=nmo_not_converged, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)
        CALL cp_fm_create(mo_notconv_fm,fm_struct_tmp,name="mo_notconv_fm",error=error)
        CALL cp_fm_struct_release(fm_struct_tmp,error=error)

        ! extract mo_notconv from mo_coeff full matrix
        ALLOCATE(eig_not_conv(nmo_not_converged),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        jj=1
        DO j=1,nset_not_conv
          i_first=inotconv_set(j,1)
          i_last=inotconv_set(j,2)
          n=i_last-i_first+1
          CALL cp_fm_to_fm_submat(mo_coeff,mo_notconv_fm,nao,n,1,i_first,1,jj,error=error)
          eig_not_conv(jj:jj+n-1) = ritz_coeff(i_first:i_last)
          jj=jj+n
        END DO

        ! allocate mo_conv sparse
        CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo_not_converged, &
              dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
        CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s)),&
              dbcsr_distribution_row_dist(cp_dbcsr_distribution(matrix_s)),col_dist)
        CALL cp_dbcsr_init_p(mo_notconv,error=error)
        CALL cp_dbcsr_create(mo_notconv,"mo_notconv",dist,dbcsr_type_no_symmetry,&
              cp_dbcsr_row_block_sizes(matrix_s),col_blk_size,0,0,dbcsr_type_real_default,error=error)

        CALL cp_dbcsr_init_p(matrix_hc,error=error)
        CALL cp_dbcsr_create(matrix_hc,"matrix_hc",dist,dbcsr_type_no_symmetry,&
              cp_dbcsr_row_block_sizes(matrix_s),col_blk_size,0,0,dbcsr_type_real_default,error=error)

        CALL cp_dbcsr_init_p(matrix_sc,error=error)
        CALL cp_dbcsr_create(matrix_sc,"matrix_sc",dist,dbcsr_type_no_symmetry,&
              cp_dbcsr_row_block_sizes(matrix_s),col_blk_size,0,0,dbcsr_type_real_default,error=error)

        CALL cp_dbcsr_init_p(matrix_z,error=error)
        CALL cp_dbcsr_create(matrix_z,"matrix_z",dist,dbcsr_type_no_symmetry,&
              cp_dbcsr_row_block_sizes(matrix_s),col_blk_size,0,0,dbcsr_type_real_default,error=error)

        CALL copy_fm_to_dbcsr(mo_notconv_fm,mo_notconv,error=error)!fm->dbcsr
        CALL cp_dbcsr_distribution_release (dist)
        CALL array_release (col_blk_size)
        CALL array_release (col_dist)


        CALL cp_dbcsr_multiply('n','n',1.0_dp,c_out,mo_notconv,0.0_dp,matrix_hc,&
             last_column=nmo_not_converged,error=error)
        CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_s,mo_notconv,0.0_dp,matrix_sc,&
             last_column=nmo_not_converged,error=error)

        CALL cp_dbcsr_copy(matrix_z,matrix_sc,error=error)
        CALL cp_dbcsr_scale_by_vector(matrix_z,eig_not_conv,side='right',error=error)
        CALL cp_dbcsr_add(matrix_z,matrix_hc,-1.0_dp,1.0_dp,error=error)

        DEALLOCATE(eig_not_conv,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)


        ! matrix_mm
        CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo_not_converged, &
              dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
        CALL cp_create_bl_distribution (row_dist, row_blk_size, nmo_not_converged, &
              dbcsr_mp_nprows(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s))))
        CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s)),&
              row_dist,col_dist)
        CALL cp_dbcsr_init_p(matrix_mm,error=error)
        CALL cp_dbcsr_create(matrix_mm,"matrix_mm",dist,dbcsr_type_no_symmetry,&
             row_blk_size,col_blk_size,0,0,dbcsr_type_real_default,error=error)
        CALL cp_dbcsr_distribution_release (dist)
        CALL array_release (col_blk_size)
        CALL array_release (col_dist)
        CALL array_release (row_blk_size)
        CALL array_release (row_dist)

        CALL cp_dbcsr_multiply('t','n',1.0_dp,mo_notconv,matrix_hc,0.0_dp,matrix_mm,&
             last_column=nmo_not_converged, error=error)

      ELSE
        mo_notconv=>mo_coeff_b
        mo_notconv_fm=>mo_coeff
        c_out => matrix_h
      END IF

      ! allocate matrix_pz using as template matrix_z
      CALL cp_dbcsr_init_p(matrix_pz,error=error)
      CALL cp_dbcsr_create(matrix_pz,"matrix_pz",cp_dbcsr_distribution(matrix_z),&
           dbcsr_type_no_symmetry,cp_dbcsr_row_block_sizes(matrix_z),&
           cp_dbcsr_col_block_sizes(matrix_z),0,0,dbcsr_type_real_default,error=error)

      IF(do_apply_preconditioner) THEN
        IF(preconditioner%in_use/=0) THEN
          CALL apply_preconditioner(preconditioner,matrix_z,matrix_pz,error=error)
        ELSE
          CALL cp_dbcsr_copy(matrix_pz,matrix_z,error=error)
        END IF
      ELSE
        CALL cp_dbcsr_copy(matrix_pz,matrix_z,error=error)
      END IF

      !allocate NMOxNMO  full matrices
       nmat = nmo_not_converged
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmat,ncol_global=nmat, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)
       CALL cp_fm_create(matrix_mm_fm,fm_struct_tmp,name="m_tmp_mxm",error=error)
       CALL cp_fm_create(matrix_mmt_fm,fm_struct_tmp,name="mt_tmp_mxm",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

      !allocate 2NMOx2NMO full matrices
       nmat2 = 2* nmo_not_converged
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nmat2,ncol_global=nmat2, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)

       CALL cp_fm_create(s_block,fm_struct_tmp,name="sb",error=error)
       CALL cp_fm_create(h_block,fm_struct_tmp,name="hb",error=error)
       CALL cp_fm_create(v_block,fm_struct_tmp,name="vb",error=error)
       CALL cp_fm_create(w_block,fm_struct_tmp,name="wb",error=error)
       ALLOCATE(evals(nmat2))
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

      ! compute CSC
       CALL cp_fm_set_all(s_block,0.0_dp, 1.0_dp,error=error)
      ! compute CHC
       CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
       CALL cp_fm_to_fm_submat(matrix_mm_fm,h_block,nmat,nmat,1,1,1,1,error=error)

      ! compute the bottom left  ZSC (top right is transpose)
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_sc,0.0_dp,matrix_mm,last_column=nmat,error=error)
      !  set the bottom left part of S[C,Z] block matrix  ZSC
      !copy sparse to full
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,s_block,nmat,nmat,1,1,1+nmat,1,error=error)
      CALL cp_fm_transpose(matrix_mm_fm,matrix_mmt_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mmt_fm,s_block,nmat,nmat,1,1,1,1+nmat,error=error)

      ! compute the bottom left  ZHC (top right is transpose)
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_hc,0.0_dp,matrix_mm,last_column=nmat,error=error)
      ! set the bottom left part of S[C,Z] block matrix  ZHC
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,h_block,nmat,nmat,1,1,1+nmat,1,error=error)
      CALL cp_fm_transpose(matrix_mm_fm,matrix_mmt_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mmt_fm,h_block,nmat,nmat,1,1,1,1+nmat,error=error)

      CALL cp_fm_release(matrix_mmt_fm,error=error)

      ! (reuse matrix_sc and matrix_hc to computr HZ and SZ)
      CALL cp_dbcsr_get_info(matrix_pz,nfullrows_total=n,nfullcols_total=k)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,c_out,matrix_pz,0.0_dp,matrix_hc,last_column=k,error=error)
      CALL cp_dbcsr_multiply('n','n',1.0_dp,matrix_s,matrix_pz,0.0_dp,matrix_sc,last_column=k,error=error)

      ! compute the bottom right  ZSZ
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_sc,0.0_dp,matrix_mm,last_column=k,error=error)
      ! set the bottom right part of S[C,Z] block matrix  ZSZ
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,s_block,nmat,nmat,1,1,1+nmat,1+nmat,error=error)

      ! compute the bottom right  ZHZ
      CALL cp_dbcsr_multiply('t','n',1.0_dp,matrix_pz,matrix_hc,0.0_dp,matrix_mm,last_column=k,error=error)
      ! set the bottom right part of H[C,Z] block matrix  ZHZ
      CALL copy_dbcsr_to_fm(matrix_mm,matrix_mm_fm,error=error)
      CALL cp_fm_to_fm_submat(matrix_mm_fm,h_block,nmat,nmat,1,1,1+nmat,1+nmat,error=error)

      CALL cp_dbcsr_release_p(matrix_mm, error=error)
      CALL cp_dbcsr_release_p(matrix_sc, error=error)
      CALL cp_dbcsr_release_p(matrix_hc, error=error)

      CALL reduce_extended_space(s_block,h_block,v_block,w_block, evals, nmat2, error=error)

       ! allocate two (nao x nmat) full matrix
       CALL cp_fm_struct_create(fm_struct_tmp,nrow_global=nao,ncol_global=nmat, &
                               context=mo_coeff%matrix_struct%context, &
                               para_env=mo_coeff%matrix_struct%para_env,error=error)
       CALL cp_fm_create(matrix_nm_fm,fm_struct_tmp,name="m_nxm",error=error)
       CALL cp_fm_create(matrix_z_fm,fm_struct_tmp,name="m_nxm",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       CALL copy_dbcsr_to_fm(matrix_pz,matrix_z_fm,error=error)
      ! extract egenvectors
       CALL cp_fm_to_fm_submat(v_block,matrix_mm_fm, nmat,nmat,1,1,1,1,error=error)
       CALL cp_fm_gemm('N','N',nao,nmat,nmat,1.0_dp,mo_notconv_fm,matrix_mm_fm,0.0_dp,matrix_nm_fm,error=error)
       CALL cp_fm_to_fm_submat(v_block,matrix_mm_fm, nmat,nmat,1+nmat,1,1,1,error=error)
       CALL cp_fm_gemm('N','N',nao,nmat,nmat,1.0_dp,matrix_z_fm,matrix_mm_fm,1.0_dp,matrix_nm_fm,error=error)

       CALL cp_dbcsr_release_p(matrix_z, error=error)
       CALL cp_dbcsr_release_p(matrix_pz, error=error)
       CALL cp_fm_release(matrix_z_fm,error=error)
       CALL cp_fm_release(s_block,error=error)
       CALL cp_fm_release(h_block,error=error)
       CALL cp_fm_release(w_block,error=error)
       CALL cp_fm_release(v_block,error=error)
       CALL cp_fm_release(matrix_mm_fm,error=error)


      ! in case some vector are already converged only a subset of vectors are copied in the MOS
      IF (nmo_converged>0) THEN
       jj=1
         DO j=1,nset_not_conv
           i_first=inotconv_set(j,1)
           i_last=inotconv_set(j,2)
           n=i_last-i_first+1
           CALL cp_fm_to_fm_submat(matrix_nm_fm,mo_coeff,nao,n,1,jj,1,i_first,error=error)
           eigenvalues(i_first:i_last) = evals(jj:jj+n-1)
           jj=jj+n
         END DO
         DEALLOCATE(iconv_set,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE (inotconv_set,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

         CALL cp_dbcsr_release_p(mo_notconv, error=error)
         CALL cp_dbcsr_release_p(c_out, error=error)
         CALL cp_fm_release(mo_notconv_fm,error=error)
      ELSE
         CALL cp_fm_to_fm(matrix_nm_fm,mo_coeff,error=error)
         eigenvalues(1:nmo) = evals(1:nmo)
      END IF
      DEALLOCATE(evals,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL cp_fm_release(matrix_nm_fm,error=error)
      CALL copy_fm_to_dbcsr(mo_coeff,mo_coeff_b,error=error)!fm->dbcsr

      t2 = m_walltime()
      IF (output_unit > 0) THEN
        WRITE(output_unit,'(T16,I5,T24,I6,T33,E12.4,2x,E12.4,T60,F8.3)')&
                iter, nmo_converged,  max_norm, min_norm, t2-t1
      END IF
      t1=m_walltime()


    END DO ! iter


    DEALLOCATE(ritz_coeff,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(iconv,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(inotconv,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(vnorm,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE generate_extended_space_sparse

! *****************************************************************************
! *****************************************************************************

  SUBROUTINE  reduce_extended_space(s_block,h_block,v_block,w_block, evals, ndim, error)

    TYPE(cp_fm_type), POINTER                :: s_block, h_block, v_block, &
                                                w_block
    REAL(dp), DIMENSION(:)                   :: evals
    INTEGER                                  :: ndim
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reduce_extended_space', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, info
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

      CALL cp_fm_to_fm(s_block,w_block,error=error)
      CALL cp_fm_cholesky_decompose(s_block,info_out=info,error=error)
      IF(info==0) THEN
         CALL cp_fm_triangular_invert(s_block,error=error)
         CALL cp_fm_cholesky_restore(H_block,ndim,S_block,w_block,"MULTIPLY",pos="RIGHT",error=error)
         CALL cp_fm_cholesky_restore(w_block,ndim,S_block,H_block,"MULTIPLY",pos="LEFT",transa="T",error=error)
         CALL choose_eigv_solver(H_block,w_block,evals,error=error)
         CALL cp_fm_cholesky_restore(w_block,ndim,S_block,v_block,"MULTIPLY",error=error)
      ELSE
! S^(-1/2)
         CALL cp_fm_power(w_block,s_block,-0.5_dp,1.0E-5_dp,info,error=error)
         CALL cp_fm_to_fm(w_block,s_block,error=error)
         CALL cp_fm_gemm('N','N',ndim,ndim,ndim,1.0_dp,H_block,s_block,0.0_dp,w_block,error=error)
         CALL cp_fm_gemm('N','N',ndim,ndim,ndim,1.0_dp,s_block,w_block,0.0_dp,H_block,error=error)
         CALL choose_eigv_solver(H_block,w_block,evals,error=error)
         CALL cp_fm_gemm('N','N',ndim,ndim,ndim,1.0_dp,s_block,w_block,0.0_dp,v_block,error=error)
      END IF

    CALL timestop(handle)

  END SUBROUTINE reduce_extended_space


END MODULE qs_scf_block_davidson
