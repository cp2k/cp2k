!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to apply a delta pulse for RTP and EMD
! *****************************************************************************

MODULE rt_delta_pulse
  USE cell_types,                      ONLY: cell_type
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                             cp_cfm_gemm
  USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_init
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_moments,                      ONLY: build_berry_moment_matrix
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_delta_pulse'

  PUBLIC :: apply_delta_pulse_periodic,&
            apply_delta_pulse


CONTAINS

! *****************************************************************************
!> \brief uses perturbation theory to get the proper initial conditions
!> \author Joost & Martin (2011)
! *****************************************************************************

  SUBROUTINE apply_delta_pulse_periodic(qs_env,mos_old,mos_new,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_old, mos_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse_periodic', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: eigenvalues_sqrt
    INTEGER                                  :: handle, icol, idir, irow, &
                                                ispin, nao, ncol_local, nmo, &
                                                nrow_global, nrow_local, nvirt
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: factor
    REAL(KIND=dp), DIMENSION(3)              :: kvec
    REAL(kind=dp), DIMENSION(:), POINTER     :: eigenvalues
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_type), POINTER               :: oo_c, oo_v, oo_vt
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: eigenvectors, mat_ks, &
                                                mat_tmp, momentum, oo_1, &
                                                oo_2, S_chol, S_inv, virtuals
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    ! we need the overlap and ks matrix for a full diagionalization
    CALL get_qs_env(qs_env, cell=cell, mos=mos, rtp=rtp, matrix_s=matrix_s, matrix_ks=matrix_ks, error=error)
    CALL get_rtp(rtp=rtp,S_inv=S_inv,error=error)
    CALL cp_fm_get_info(S_inv,nrow_global=nrow_global,error=error)

    CALL cp_fm_create(S_chol, matrix_struct=S_inv%matrix_struct, name="S_chol", error=error)
    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,S_chol,error=error)
    CALL cp_fm_cholesky_decompose(S_chol,error=error)
    NULLIFY(mat_ks,eigenvectors,mat_tmp)
    CALL cp_fm_create(mat_ks, matrix_struct=S_inv%matrix_struct, name="mat_ks", error=error)
    CALL cp_fm_create(eigenvectors, matrix_struct=S_inv%matrix_struct, name="eigenvectors", error=error)

    DO ispin=1,SIZE(matrix_ks)
       ALLOCATE(eigenvalues(nrow_global))
       CALL cp_fm_create(mat_tmp, matrix_struct=S_inv%matrix_struct, name="mat_tmp", error=error)

       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,mat_ks,error=error)
       CALL cp_fm_cholesky_reduce(mat_ks,S_chol,error=error)
       CALL cp_fm_syevd(mat_ks,mat_tmp,eigenvalues,error=error)
       CALL cp_fm_cholesky_restore(mat_tmp,nrow_global,S_chol,eigenvectors,"SOLVE",error=error)

       ! virtuals
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, nao=nao, nmo=nmo)
       nvirt=nao-nmo
       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nrow_global, ncol_global=nvirt, error=error)
       CALL cp_fm_create(virtuals, matrix_struct=fm_struct_tmp, name="virtuals", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       CALL cp_fm_to_fm(eigenvectors,virtuals,nvirt,nmo+1,1)

       ! occupied
       CALL cp_fm_to_fm(eigenvectors,mos_old(2*ispin-1)%matrix,nmo,1,1)

       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nvirt, ncol_global=nmo, error=error)
       CALL cp_fm_create(momentum, matrix_struct=fm_struct_tmp, name="momentum", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       ! the momentum operator (in a given direction)
       CALL cp_fm_set_all( mos_new(2*ispin-1)%matrix,0.0_dp,error=error)

       ! the prefactor (strength of the electric field)
       kvec(:) = cell%h_inv(1,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(1)+&
                 cell%h_inv(2,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(2)+&
                 cell%h_inv(3,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(3)
       kvec = - kvec * twopi * qs_env%dft_control%rtp_control%delta_pulse_scale

       DO idir=1,3
          factor= kvec(idir)
          IF (factor.NE.0.0_dp) THEN
             CALL cp_dbcsr_sm_fm_multiply(matrix_s(idir+1)%matrix,mos_old(2*ispin-1)%matrix,&
                                          mos_old(2*ispin)%matrix,ncol=nmo,error=error)
             CALL cp_fm_scale_and_add(1.0_dp, mos_new(2*ispin-1)%matrix,factor,mos_old(2*ispin)%matrix,error)
          ENDIF
       ENDDO

       CALL cp_fm_gemm('T','N',nvirt,nmo,nao,1.0_dp, virtuals, mos_new(2*ispin-1)%matrix,0.0_dp,momentum,error=error)


       ! the tricky bit ... rescale by the eigenvalue difference
       CALL cp_fm_get_info(momentum, nrow_local=nrow_local, ncol_local=ncol_local,&
                           row_indices=row_indices,col_indices=col_indices,local_data=local_data,error=error)
       DO icol=1,ncol_local
          DO irow=1,nrow_local
             factor= 1 / (eigenvalues(col_indices(icol))-eigenvalues(nmo+row_indices(irow)))
             local_data(irow,icol)=factor*local_data(irow,icol)
          ENDDO
       ENDDO
       CALL cp_fm_release(mat_tmp,error=error)
       DEALLOCATE(eigenvalues)


       ! now obtain the initial condition in mos_old
       CALL cp_fm_to_fm(eigenvectors,mos_old(2*ispin-1)%matrix,nmo,1,1)
       CALL cp_fm_gemm("N","N",nao,nmo,nvirt,1.0_dp,virtuals,momentum,0.0_dp,mos_old(2*ispin)%matrix,error=error)

       CALL cp_fm_release(virtuals, error=error)
       CALL cp_fm_release(momentum, error=error)

       ! orthonormalize afterwards
       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nmo, ncol_global=nmo, error=error)
       CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1", error=error)
       CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp, error=error)

       CALL cp_fm_create(mat_tmp, matrix_struct=mos_old(2*ispin-1)%matrix%matrix_struct, name="tmp_mat", error=error)
       ! get the complex overlap matrix
       ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
       CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,mos_old(2*ispin-1)%matrix,&
                                          mat_tmp,ncol=nmo,error=error)

       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*ispin-1)%matrix,mat_tmp,0.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,-1.0_dp,mos_old(2*ispin)%matrix,mat_tmp,0.0_dp,oo_2,error=error)

       CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,mos_old(2*ispin)%matrix,&
                                          mat_tmp,ncol=nmo,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*ispin)%matrix,mat_tmp,1.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*ispin-1)%matrix,mat_tmp,1.0_dp,oo_2,error=error)
       CALL cp_fm_release(mat_tmp, error=error)

       CALL cp_cfm_create ( oo_c, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_v, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_vt, oo_1 % matrix_struct,error=error)
       oo_c % local_data = CMPLX(oo_1%local_data,oo_2%local_data,KIND=dp)

       ! compute inv(sqrt(overlap))
       ALLOCATE(eigenvalues(nmo))
       ALLOCATE(eigenvalues_sqrt(nmo))
       CALL cp_cfm_heevd(oo_c,oo_v,eigenvalues,error)
       eigenvalues_sqrt=1/SQRT(eigenvalues)
       CALL cp_cfm_to_cfm(oo_v,oo_vt,error=error)
       CALL cp_cfm_column_scale(oo_v,eigenvalues_sqrt)
       DEALLOCATE(eigenvalues)
       DEALLOCATE(eigenvalues_sqrt)
       CALL cp_cfm_gemm('N','C',nmo,nmo,nmo,(1.0_dp,0.0_dp),&
            oo_v,oo_vt,(0.0_dp,0.0_dp),oo_c,error=error)
       oo_1%local_data=REAL(oo_c%local_data,KIND=dp)
       oo_2%local_data=AIMAG(oo_c%local_data)
       CALL cp_cfm_release(oo_c,error=error)
       CALL cp_cfm_release(oo_v,error=error)
       CALL cp_cfm_release(oo_vt,error=error)

       ! use this to compute the orthonormal vectors
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*ispin-1)%matrix,oo_1,0.0_dp,mos_new(2*ispin-1)%matrix,error=error)
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*ispin-1)%matrix,oo_2,0.0_dp,mos_new(2*ispin)%matrix,error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,-1.0_dp,mos_old(2*ispin)%matrix,oo_2,0.0_dp,mos_old(2*ispin-1)%matrix,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,mos_old(2*ispin-1)%matrix,1.0_dp,mos_new(2*ispin-1)%matrix,error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*ispin)%matrix,oo_1,1.0_dp,mos_new(2*ispin)%matrix,error=error)
       CALL cp_fm_to_fm(mos_new(2*ispin)%matrix,mos_old(2*ispin)%matrix,error)

       CALL cp_fm_release(oo_1,error=error)
       CALL cp_fm_release(oo_2,error=error)
    END DO

    CALL cp_fm_release(S_chol, error=error)
    CALL cp_fm_release(mat_ks, error=error)
    CALL cp_fm_release(eigenvectors, error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_delta_pulse_periodic

! *****************************************************************************
!> \brief applies exp(ikr) to the wavefunction.... stored in mos_old...
!> \author Joost & Martin (2011)
! *****************************************************************************

  SUBROUTINE apply_delta_pulse(qs_env,mos_old,mos_new,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_old, mos_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: eigenvalues_sqrt
    INTEGER                                  :: handle, i, nao, nmo
    REAL(KIND=dp), DIMENSION(3)              :: kvec
    REAL(kind=dp), DIMENSION(:), POINTER     :: eigenvalues
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_type), POINTER               :: oo_c, oo_v, oo_vt
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mat_S, oo_1, oo_2, S_inv, tmp
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    ! we need the inverse overlap

    CALL get_qs_env(qs_env, mos=mos, rtp=rtp, matrix_s=matrix_s, error=error)
    CALL get_rtp(rtp=rtp,S_inv=S_inv,error=error)

    CALL cp_fm_create(tmp, matrix_struct=S_inv%matrix_struct, name="tmp_mat", error=error)

    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,S_inv,error=error)
    CALL cp_fm_cholesky_decompose(S_inv,error=error)
    CALL cp_fm_cholesky_invert(S_inv,error=error)
    CALL cp_fm_upper_to_full(S_inv,tmp,error=error)

    CALL cp_fm_create(mat_S, matrix_struct=S_inv%matrix_struct, name="mat_S", error=error)
    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,mat_S,error=error)
    CALL cp_fm_upper_to_full(mat_S,tmp,error=error)

    CALL cp_fm_release(tmp,error=error)

    ! we need the berry matrix
    CALL get_qs_env(qs_env, cell=cell, error=error)

    ! direction ... unscaled, this will yield a exp(ikr) that is periodic with the cell
    kvec(:) = cell%h_inv(1,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(1)+&
              cell%h_inv(2,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(2)+&
              cell%h_inv(3,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(3)
    kvec = - kvec * twopi
    ! scaling will make the things not periodic with the cell, which would only be good for gas phase systems ?
    kvec(:) = qs_env%dft_control%rtp_control%delta_pulse_scale * kvec

    ALLOCATE(cosmat, sinmat)
    CALL cp_dbcsr_init(cosmat, error=error)
    CALL cp_dbcsr_init(sinmat, error=error)
    CALL cp_dbcsr_copy(cosmat,matrix_s(1)%matrix,'COS MOM',error=error)
    CALL cp_dbcsr_copy(sinmat,matrix_s(1)%matrix,'SIN MOM',error=error)
    CALL build_berry_moment_matrix(qs_env,cosmat,sinmat,kvec,error)

    ! apply inv(S)*operator to C
    DO i=1,SIZE(mos)
       CALL get_mo_set(mos(i)%mo_set, nao=nao,nmo=nmo)
       CALL cp_dbcsr_sm_fm_multiply(cosmat, mos(i)%mo_set%mo_coeff, mos_new(2*i-1)%matrix, ncol=nmo, error=error)
       CALL cp_dbcsr_sm_fm_multiply(sinmat, mos(i)%mo_set%mo_coeff, mos_new(2*i)%matrix, ncol=nmo, error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,S_inv,mos_new(2*i-1)%matrix,0.0_dp,mos_old(2*i-1)%matrix,error=error)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,S_inv,mos_new(2*i)%matrix,0.0_dp,mos_old(2*i)%matrix,error=error)

       ! in a finite basis, unfortunately, inv(S)*operator is not unitary, so orthonormalize afterwards
       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nmo, ncol_global=nmo, error=error)
       CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1", error=error)
       CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp, error=error)

       CALL cp_fm_create(tmp, matrix_struct=mos_old(2*i-1)%matrix%matrix_struct, name="tmp_mat", error=error)
       ! get the complex overlap matrix
       ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,mat_S,mos_old(2*i-1)%matrix,0.0_dp,tmp,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*i-1)%matrix,tmp,0.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,-1.0_dp,mos_old(2*i)%matrix,tmp,0.0_dp,oo_2,error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,mat_S,mos_old(2*i)%matrix,0.0_dp,tmp,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*i)%matrix,tmp,1.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*i-1)%matrix,tmp,1.0_dp,oo_2,error=error)
       CALL cp_fm_release(tmp, error=error)

       CALL cp_cfm_create ( oo_c, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_v, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_vt, oo_1 % matrix_struct,error=error)
       oo_c % local_data = CMPLX(oo_1%local_data,oo_2%local_data,KIND=dp)

       ! compute inv(sqrt(overlap))
       ALLOCATE(eigenvalues(nmo))
       ALLOCATE(eigenvalues_sqrt(nmo))
       CALL cp_cfm_heevd(oo_c,oo_v,eigenvalues,error)
       eigenvalues_sqrt=1/SQRT(eigenvalues)
       CALL cp_cfm_to_cfm(oo_v,oo_vt,error=error)
       CALL cp_cfm_column_scale(oo_v,eigenvalues_sqrt)
       DEALLOCATE(eigenvalues)
       DEALLOCATE(eigenvalues_sqrt)
       CALL cp_cfm_gemm('N','C',nmo,nmo,nmo,(1.0_dp,0.0_dp),&
            oo_v,oo_vt,(0.0_dp,0.0_dp),oo_c,error=error)
       oo_1%local_data=REAL(oo_c%local_data,KIND=dp)
       oo_2%local_data=AIMAG(oo_c%local_data)
       CALL cp_cfm_release(oo_c,error=error)
       CALL cp_cfm_release(oo_v,error=error)
       CALL cp_cfm_release(oo_vt,error=error)

       ! use this to compute the orthonormal vectors
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*i-1)%matrix,oo_1,0.0_dp,mos_new(2*i-1)%matrix,error=error)
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*i-1)%matrix,oo_2,0.0_dp,mos_new(2*i)%matrix,error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,-1.0_dp,mos_old(2*i)%matrix,oo_2,0.0_dp,mos_old(2*i-1)%matrix,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,mos_old(2*i-1)%matrix,1.0_dp,mos_new(2*i-1)%matrix,error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*i)%matrix,oo_1,1.0_dp,mos_new(2*i)%matrix,error=error)
       CALL cp_fm_to_fm(mos_new(2*i)%matrix,mos_old(2*i)%matrix,error)

       CALL cp_fm_release(oo_1,error=error)
       CALL cp_fm_release(oo_2,error=error)
    END DO

    CALL cp_fm_release(mat_S, error=error)

    CALL cp_dbcsr_deallocate_matrix(cosmat, error)
    CALL cp_dbcsr_deallocate_matrix(sinmat, error)


    CALL timestop(handle)

  END SUBROUTINE apply_delta_pulse

END MODULE rt_delta_pulse
