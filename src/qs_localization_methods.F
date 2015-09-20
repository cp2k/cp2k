!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Localization methods such as 2x2 Jacobi rotations
!>                                   Steepest Decents
!>                                   Conjugate Gradient
!> \par History
!>      Initial parallellization of jacobi (JVDV 07.2003)
!>      direct minimization using exponential parametrization (JVDV 09.2003)
!>      crazy rotations go fast (JVDV 10.2003)
!> \author CJM (04.2003)
! *****************************************************************************
MODULE qs_localization_methods
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                             cp_cfm_gemm,&
                                             cp_cfm_schur_product
  USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_get_element,&
                                             cp_cfm_get_info,&
                                             cp_cfm_p_type,&
                                             cp_cfm_release,&
                                             cp_cfm_set_all,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_external_control,             ONLY: external_control
  USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace,&
                                             cp_fm_transpose
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_get_submatrix, &
       cp_fm_maxabsrownorm, cp_fm_maxabsval, cp_fm_p_type, cp_fm_release, &
       cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE matrix_exp,                      ONLY: exp_pade_real,&
                                             get_nsquare_norder
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast,&
                                             mp_max,&
                                             mp_sendrecv,&
                                             mp_sum,&
                                             mp_sync
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PUBLIC :: initialize_weights, crazy_rotations,&
            direct_mini, rotate_orbitals, approx_l1_norm_sd, jacobi_rotations

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_localization_methods'

  PRIVATE

   TYPE set_c_1d_type
     COMPLEX(KIND=dp), POINTER, DIMENSION(:) :: c_array
   END TYPE
   TYPE set_c_2d_type
     COMPLEX(KIND=dp), POINTER, DIMENSION(:,:) :: c_array
   END TYPE

CONTAINS
! *****************************************************************************
!> \brief ...
!> \param C ...
!> \param iterations ...
!> \param eps ...
!> \param converged ...
!> \param sweeps ...
! *****************************************************************************
  SUBROUTINE approx_l1_norm_sd( C, iterations, eps, converged, sweeps)
    TYPE(cp_fm_type), POINTER                :: C
    INTEGER, INTENT(IN)                      :: iterations
    REAL(KIND=dp), INTENT(IN)                :: eps
    LOGICAL, INTENT(INOUT)                   :: converged
    INTEGER, INTENT(INOUT)                   :: sweeps

    CHARACTER(len=*), PARAMETER :: routineN = 'approx_l1_norm_sd', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: taylor_order = 100
    REAL(KIND=dp), PARAMETER                 :: alpha = 0.1_dp, &
                                                f2_eps = 0.01_dp

    INTEGER                                  :: handle, i, istep, k, n, &
                                                ncol_local, nrow_local, &
                                                output_unit, p
    REAL(KIND=dp)                            :: expfactor, f2, f2old, gnorm, &
                                                tnorm
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_k_k
    TYPE(cp_fm_type), POINTER                :: CTmp, G, Gp1, Gp2, U
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)

    NULLIFY(logger,CTmp,U,G,Gp1,Gp2,context,para_env,fm_struct_k_k)

    logger => cp_get_default_logger()
    output_unit= cp_logger_get_default_io_unit(logger)

    CALL cp_fm_struct_get(C%matrix_struct, nrow_global=n, ncol_global=k, &
                          nrow_local = nrow_local, ncol_local = ncol_local, &
                          para_env=para_env, context=context)
    CALL cp_fm_struct_create(fm_struct_k_k, para_env=para_env, context=context, &
                             nrow_global=k, ncol_global=k)
    CALL cp_fm_create(CTmp, C%matrix_struct)
    CALL cp_fm_create(U   , fm_struct_k_k)
    CALL cp_fm_create(G   , fm_struct_k_k)
    CALL cp_fm_create(Gp1 , fm_struct_k_k)
    CALL cp_fm_create(Gp2 , fm_struct_k_k)
    !
    ! printing
    IF(output_unit>0) THEN
       WRITE(output_unit,'(1X)')
       WRITE(output_unit,'(2X,A)') '-----------------------------------------------------------------------------'
       WRITE(output_unit,'(A,I5)'   ) '      Nbr iterations =',iterations
       WRITE(output_unit,'(A,E10.2)') '     eps convergence =',eps
       WRITE(output_unit,'(A,I5)'   ) '    Max Taylor order =',taylor_order
       WRITE(output_unit,'(A,E10.2)') '              f2 eps =',f2_eps
       WRITE(output_unit,'(A,E10.2)') '               alpha =',alpha
       WRITE(output_unit,'(A)') '     iteration    approx_l1_norm    g_norm   rel_err'
    ENDIF
    !
    f2old = 0.0_dp
    converged = .FALSE.
    !
    ! Start the steepest descent
    DO istep=1,iterations
       !
       !-------------------------------------------------------------------
       ! compute f_2
       ! f_2(x)=(x^2+eps)^1/2
       f2 = 0.0_dp
       DO p=1,ncol_local   ! p
          DO i=1,nrow_local! i
             f2 = f2 + SQRT( C%local_data(i,p)**2 + f2_eps )
          ENDDO
       ENDDO
       CALL mp_sum(f2,C%matrix_struct%para_env%group)
       !write(*,*) 'qs_localize: f_2=',f2
       !-------------------------------------------------------------------
       ! compute the derivative of f_2
       ! f_2(x)=(x^2+eps)^1/2
       DO p=1,ncol_local   ! p
          DO i=1,nrow_local! i
             CTmp%local_data(i,p) = C%local_data(i,p) / SQRT( C%local_data(i,p)**2 + f2_eps )
          ENDDO
       ENDDO
       CALL cp_gemm('T','N',k,k,n,1.0_dp,CTmp,C,0.0_dp,G)
       ! antisymmetrize
       CALL cp_fm_transpose(G,U)
       CALL cp_fm_scale_and_add(-0.5_dp,G,0.5_dp,U)
       !
       !-------------------------------------------------------------------
       !
       CALL cp_fm_frobenius_norm(G,gnorm)
       !write(*,*) 'qs_localize: norm(G)=',gnorm
       !
       ! rescale for steepest descent
       CALL cp_fm_scale(-alpha, G)
       !
       ! compute unitary transform
       ! zeroth order
       CALL cp_fm_set_all(U,0.0_dp,1.0_dp)
       ! first order
       expfactor = 1.0_dp
       CALL cp_fm_scale_and_add(1.0_dp,U,expfactor,G)
       CALL cp_fm_frobenius_norm(G,tnorm)
       !write(*,*) 'Taylor expansion i=',1,' norm(X^i)/i!=',tnorm
       IF(tnorm.GT.1.0E-10_dp) THEN
          ! other orders
          CALL cp_fm_to_fm(G,Gp1)
          DO i = 2,taylor_order
             ! new power of G
             CALL cp_gemm('N','N',k,k,k,1.0_dp,G,Gp1,0.0_dp,Gp2)
             CALL cp_fm_to_fm(Gp2,Gp1)
             ! add to the taylor expansion so far
             expfactor = expfactor / REAL(i,KIND=dp)
             CALL cp_fm_scale_and_add(1.0_dp,U,expfactor,Gp1)
             CALL cp_fm_frobenius_norm(Gp1,tnorm)
             !write(*,*) 'Taylor expansion i=',i,' norm(X^i)/i!=',tnorm*expfactor
             IF(tnorm*expfactor.LT.1.0E-10_dp) EXIT
          ENDDO
       ENDIF
       !
       ! incrementaly rotate the MOs
       CALL cp_gemm('N','N',n,k,k,1.0_dp,C,U,0.0_dp,CTmp)
       CALL cp_fm_to_fm(CTmp,C)
       !
       ! printing
       IF(output_unit.GT.0) THEN
          WRITE(output_unit,'(10X,I4,E18.10,2E10.2)') istep,f2,gnorm,ABS((f2-f2old)/f2)
       ENDIF
       !
       ! Are we done?
       sweeps = istep
       !IF(gnorm.LE.grad_thresh.AND.ABS((f2-f2old)/f2).LE.f2_thresh.AND.istep.GT.1) THEN
       IF(ABS((f2-f2old)/f2).LE.eps.AND.istep.GT.1) THEN
          converged = .TRUE.
          EXIT
       ENDIF
       f2old = f2
    ENDDO
    !
    ! here we should do one refine step to enforce C'*S*C=1 for any case
    !
    ! Print the final result
    IF(output_unit.GT.0) WRITE(output_unit,'(A,E16.10)')' sparseness function f2 = ',f2
    !
    ! sparsity
    !DO i=1,size(thresh,1)
    !   gnorm = 0.0_dp
    !   DO o=1,ncol_local
    !      DO p=1,nrow_local
    !         IF(ABS(C%local_data(p,o)).GT.thresh(i)) THEN
    !            gnorm = gnorm + 1.0_dp
    !         ENDIF
    !      ENDDO
    !   ENDDO
    !   CALL mp_sum(gnorm,C%matrix_struct%para_env%group)
    !   IF(output_unit.GT.0) THEN
    !      WRITE(output_unit,*) 'qs_localize: ratio2=',gnorm / ( REAL(k,KIND=dp)*REAL(n,KIND=dp) ),thresh(i)
    !   ENDIF
    !ENDDO
    !
    ! deallocate
    CALL cp_fm_struct_release(fm_struct_k_k)
    CALL cp_fm_release(CTmp)
    CALL cp_fm_release(U)
    CALL cp_fm_release(G)
    CALL cp_fm_release(Gp1)
    CALL cp_fm_release(Gp2)

    CALL timestop(handle)

  END SUBROUTINE approx_l1_norm_sd
! *****************************************************************************
!> \brief ...
!> \param cell ...
!> \param weights ...
! *****************************************************************************
  SUBROUTINE initialize_weights ( cell, weights )
    TYPE(cell_type), INTENT(IN)              :: cell
    REAL(KIND=dp), DIMENSION(:)              :: weights

    REAL(KIND=dp), DIMENSION(3, 3)           :: metric

    metric = 0.0_dp
    CALL dgemm('T','N',3,3,3,1._dp,cell%hmat,3,cell%hmat,3,0.0_dp,metric,3)

    weights(1) = METRIC(1,1)-METRIC(1,2)-METRIC(1,3)
    weights(2) = METRIC(2,2)-METRIC(1,2)-METRIC(2,3)
    weights(3) = METRIC(3,3)-METRIC(1,3)-METRIC(2,3)
    weights(4) = METRIC(1,2)
    weights(5) = METRIC(1,3)
    weights(6) = METRIC(2,3)

  END SUBROUTINE initialize_weights

! *****************************************************************************
!> \brief wrapper for the jacobi routines, should be removed if jacobi_rot_para
!>        can deal with serial para_envs.
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param para_env ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param sweeps ...
!> \param out_each ...
!> \param target_time ...
!> \param start_time ...
!> \par History
!> \author Joost VandeVondele (02.2010)
! *****************************************************************************
  SUBROUTINE jacobi_rotations(weights, zij, vectors, para_env, max_iter, eps_localization, &
              sweeps, out_each, target_time, start_time)

    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: ZIJ( :, : )
    TYPE(cp_fm_type), POINTER                :: vectors
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: max_iter
    REAL(KIND=dp), INTENT(IN)                :: eps_localization
    INTEGER                                  :: sweeps
    INTEGER, INTENT(IN)                      :: out_each
    REAL(dp)                                 :: target_time, start_time

    IF (para_env%num_pe==1) THEN
        CALL jacobi_rotations_serial ( weights, zij, vectors, max_iter, eps_localization, sweeps, out_each)
    ELSE
        CALL jacobi_rot_para(weights, zij, vectors, para_env, max_iter, eps_localization, &
              sweeps, out_each, target_time, start_time)
    ENDIF

  END SUBROUTINE jacobi_rotations

! *****************************************************************************
!> \brief this routine, private to the module is a serial backup, till we have jacobi_rot_para to work in serial
!>        while the routine below works in parallel, it is too slow to be useful
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param sweeps ...
!> \param out_each ...
! *****************************************************************************
  SUBROUTINE jacobi_rotations_serial ( weights, zij, vectors, max_iter, eps_localization, sweeps, out_each)
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: ZIJ( :, : )
    TYPE(cp_fm_type), POINTER                :: vectors
    INTEGER, INTENT(IN)                      :: max_iter
    REAL(KIND=dp), INTENT(IN)                :: eps_localization
    INTEGER                                  :: sweeps
    INTEGER, INTENT(IN)                      :: out_each

    CHARACTER(len=*), PARAMETER :: routineN = 'jacobi_rotations_serial', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), POINTER                :: mii( : ), mij( : ), mjj( : )
    INTEGER                                  :: dim2, handle, idim, istat, &
                                                istate, jstate, nstate, &
                                                unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ct, st, t1, t2, theta, &
                                                tolerance
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: c_zij
    TYPE(cp_cfm_type), POINTER               :: c_rmat
    TYPE(cp_fm_type), POINTER                :: rmat

    CALL timeset(routineN,handle)

    failure=.FALSE.

    dim2 = SIZE(zij,2)
    NULLIFY(rmat,c_rmat,c_zij)
    ALLOCATE(c_zij(dim2),STAT=istat)
    CPASSERT(istat==0)
    NULLIFY(mii,mij,mjj)
    ALLOCATE(mii(dim2),mij(dim2),mjj(dim2),STAT=istat)
    CPASSERT(istat==0)

    CALL cp_fm_create ( rmat, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_set_all ( rmat, 0._dp, 1._dp)

    CALL cp_cfm_create ( c_rmat, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_set_all ( c_rmat, (0._dp,0._dp) , (1._dp,0._dp))
    DO idim=1,dim2
       NULLIFY(c_zij(idim)%matrix)
       CALL cp_cfm_create ( c_zij(idim)%matrix, zij ( 1, 1 ) % matrix % matrix_struct)
       c_zij(idim)%matrix% local_data = CMPLX (zij(1,idim) % matrix % local_data, &
                    zij(2,idim) % matrix % local_data, dp )
    ENDDO

    CALL cp_fm_get_info ( rmat , nrow_global = nstate)
    tolerance = 1.0e10_dp
    sweeps = 0
    unit_nr = -1
    IF (rmat%matrix_struct%para_env%mepos .EQ. rmat%matrix_struct%para_env%source ) THEN
      unit_nr = cp_logger_get_default_unit_nr()
      WRITE(unit_nr,'(T10,A )') "Localization by iterative Jacobi rotation"
      WRITE(unit_nr,'(T20,A12,T40,A12,A8 )') "Iteration", "Tolerance", " Time "
    END IF
! do jacobi sweeps until converged
    DO WHILE ( tolerance >= eps_localization .AND. sweeps < max_iter )
      sweeps = sweeps + 1
      t1=m_walltime()
      DO istate = 1, nstate
        DO jstate = istate + 1, nstate
          DO idim = 1,dim2
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, istate, istate, mii(idim) )
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, istate, jstate, mij(idim) )
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, jstate, jstate, mjj(idim) )
          END DO
          CALL get_angle ( mii, mjj, mij, weights, theta )
          st = SIN ( theta )
          ct = COS ( theta )
          CALL rotate_zij ( istate, jstate, st, ct, c_zij )
          CALL rotate_rmat ( istate, jstate, st, ct, c_rmat )
        END DO
      END DO
      CALL check_tolerance ( c_zij, weights, tolerance )
      t2 = m_walltime()
      IF (unit_nr>0 .AND. MODULO(sweeps,out_each)==0) THEN
         WRITE(unit_nr,'(T20,I12,T40,E12.4,F8.3)') sweeps,tolerance, t2-t1
         CALL m_flush(unit_nr)
      ENDIF
    END DO

    DO idim=1,dim2
       zij(1,idim) % matrix % local_data =  REAL ( c_zij(idim)%matrix% local_data , dp )
       zij(2,idim) % matrix % local_data = AIMAG ( c_zij(idim)%matrix% local_data )
       CALL cp_cfm_release( c_zij(idim)%matrix)
    ENDDO
    DEALLOCATE(c_zij,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(mii,mij,mjj, STAT=istat)
    CPASSERT(istat==0)
    rmat % local_data = REAL ( c_rmat%local_data, dp )
    CALL cp_cfm_release( c_rmat)
    CALL rotate_orbitals ( rmat, vectors )
    CALL cp_fm_release ( rmat)

    CALL timestop(handle)

  END SUBROUTINE jacobi_rotations_serial
! *****************************************************************************
!> \brief ...
!> \param istate ...
!> \param jstate ...
!> \param st ...
!> \param ct ...
!> \param zij ...
! *****************************************************************************
  SUBROUTINE rotate_zij ( istate, jstate, st, ct, zij )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: istate, jstate
    TYPE ( cp_cfm_p_type ) :: zij ( : )
    REAL ( KIND = dp ), INTENT ( IN ) :: st, ct
! Locals
    INTEGER :: idim, nstate
    COMPLEX ( KIND = dp ) :: st_cmplx
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#else
    INTEGER :: stride
#endif

    st_cmplx = CMPLX ( st, 0.0_dp, dp )
    CALL cp_cfm_get_info ( zij ( 1 ) % matrix, nrow_global = nstate)
    DO  idim = 1, SIZE(zij,1)
#if defined(__SCALAPACK)
      desc(:) = zij(idim) % matrix %matrix_struct%descriptor(:)
      CALL pzrot( nstate , zij(idim) % matrix % local_data (1,1) , 1 , istate, desc, 1, &
                  zij(idim) % matrix % local_data (1,1) , 1, jstate,  desc, 1, ct, st_cmplx )
      CALL pzrot( nstate , zij(idim) % matrix % local_data (1,1) , istate , 1, desc, nstate, &
                  zij(idim) % matrix % local_data (1,1) , jstate, 1, desc, nstate, ct, st_cmplx )
#else
      CALL zrot( nstate, zij(idim) % matrix % local_data ( 1, istate ), 1, &
                 zij(idim) % matrix % local_data ( 1, jstate ), 1, ct, st_cmplx )
      stride = SIZE(zij(idim) % matrix % local_data,1)
      CALL zrot( nstate, zij(idim) %matrix % local_data ( istate, 1 ), stride, &
                 zij(idim) %matrix % local_data ( jstate, 1 ), stride, ct, st_cmplx )
#endif
    END DO
  END SUBROUTINE rotate_zij
! *****************************************************************************
!> \brief ...
!> \param istate ...
!> \param jstate ...
!> \param st ...
!> \param ct ...
!> \param rmat ...
! *****************************************************************************
  SUBROUTINE rotate_rmat ( istate, jstate, st, ct, rmat )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: istate, jstate
    TYPE ( cp_cfm_type ), POINTER :: rmat
    REAL ( KIND = dp ), INTENT ( IN ) :: ct, st
! Locals
    INTEGER :: nstate
    COMPLEX ( KIND = dp ) :: st_cmplx
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

    st_cmplx = CMPLX ( st, 0.0_dp, dp )
    CALL cp_cfm_get_info (  rmat, nrow_global = nstate)
#if defined(__SCALAPACK)
    desc ( : ) = rmat % matrix_struct % descriptor(:)
    CALL pzrot( nstate , rmat % local_data (1,1) , 1 , istate,  desc, 1, &
                  rmat % local_data (1,1) , 1, jstate, desc, 1, ct, st_cmplx )
#else
    CALL zrot( nstate, rmat % local_data(1,istate), 1, rmat % local_data(1,jstate), 1, ct, st_cmplx )
#endif
  END SUBROUTINE rotate_rmat
! *****************************************************************************
!> \brief ...
!> \param mii ...
!> \param mjj ...
!> \param mij ...
!> \param weights ...
!> \param theta ...
! *****************************************************************************
  SUBROUTINE get_angle ( mii, mjj, mij,  weights, theta )
    COMPLEX(KIND=dp), POINTER                :: mii( : ), mjj( : ), mij( : )
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    REAL(KIND=dp), INTENT(OUT)               :: theta

    COMPLEX(KIND=dp)                         :: z11, z12, z22
    INTEGER                                  :: dim_m, idim
    REAL(KIND=dp)                            :: a12, b12, d2, ratio

    a12 = 0.0_dp
    b12 = 0.0_dp
    dim_m = SIZE(mii)
    DO idim = 1, dim_m
      z11 = mii(idim)
      z22 = mjj(idim)
      z12 = mij(idim)
      a12 = a12 + weights ( idim ) * REAL ( CONJG ( z12 ) * ( z11 - z22 ), KIND=dp )
      b12 = b12 + weights ( idim ) * REAL ( (  z12 * CONJG ( z12 ) -  &
            0.25_dp * ( z11 - z22 ) * ( CONJG ( z11 ) - CONJG ( z22 ) ) ), KIND=dp)
    END DO
    IF ( ABS ( b12 ) > 1.e-10_dp ) THEN
      ratio = -a12/b12
      theta = 0.25_dp * ATAN ( ratio )
    ELSEIF ( ABS ( b12 ) < 1.e-10_dp ) THEN
      b12 = 0.0_dp
      theta = 0.0_dp
    ELSE
      theta = 0.25_dp * pi
    ENDIF
! Check second derivative info
    d2 = a12 * SIN ( 4._dp * theta ) - b12 * COS ( 4._dp * theta )
    IF ( d2 <=  0._dp ) THEN ! go to the maximum, not the minimum
       IF (theta > 0.0_dp) THEN ! make theta as small as possible
        theta = theta - 0.25_dp * pi
       ELSE
        theta = theta + 0.25_dp * pi
       ENDIF
    ENDIF
  END SUBROUTINE get_angle
! *****************************************************************************
!> \brief ...
!> \param zij ...
!> \param weights ...
!> \param tolerance ...
! *****************************************************************************
  SUBROUTINE check_tolerance ( zij, weights, tolerance )
    TYPE(cp_cfm_p_type)                      :: zij( : )
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    REAL(KIND=dp), INTENT(OUT)               :: tolerance

    CHARACTER(len=*), PARAMETER :: routineN = 'check_tolerance', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_fm_type), POINTER                :: force

    CALL timeset(routineN,handle)

! compute gradient at t=0

    NULLIFY ( force )
    CALL cp_fm_create ( force, zij( 1 ) % matrix % matrix_struct)
    CALL cp_fm_set_all ( force, 0._dp)
    CALL grad_at_0 ( zij, weights, force )
    CALL cp_fm_maxabsval ( force, tolerance)
    CALL cp_fm_release ( force)

    CALL timestop(handle)

  END SUBROUTINE check_tolerance
! *****************************************************************************
!> \brief ...
!> \param rmat ...
!> \param vectors ...
! *****************************************************************************
  SUBROUTINE rotate_orbitals ( rmat, vectors )
    TYPE(cp_fm_type), POINTER                :: rmat, vectors

    INTEGER                                  :: k, n
    TYPE(cp_fm_type), POINTER                :: wf

    NULLIFY ( wf )
    CALL cp_fm_create ( wf, vectors % matrix_struct)
    CALL cp_fm_get_info ( vectors, nrow_global = n, ncol_global=k)
    CALL cp_gemm("N", "N", n, k, k, 1.0_dp, vectors, rmat, 0.0_dp, wf)
    CALL cp_fm_to_fm ( wf, vectors)
    CALL cp_fm_release ( wf)
  END SUBROUTINE rotate_orbitals
! *****************************************************************************
!> \brief ...
!> \param diag ...
!> \param weights ...
!> \param matrix ...
!> \param ndim ...
! *****************************************************************************
  SUBROUTINE gradsq_at_0 ( diag, weights, matrix,ndim )
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: diag
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_type), POINTER                :: matrix
    INTEGER, INTENT(IN)                      :: ndim

    COMPLEX(KIND=dp)                         :: zii, zjj
    INTEGER                                  :: idim, istate, jstate, &
                                                ncol_local, nrow_global, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: gradsq_ij

    CALL cp_fm_get_info ( matrix, nrow_local = nrow_local,  &
                   ncol_local = ncol_local, nrow_global = nrow_global, &
                   row_indices = row_indices, col_indices = col_indices)

    DO istate = 1, nrow_local
      DO jstate = 1, ncol_local
! get real and imaginary parts
        gradsq_ij = 0.0_dp
        DO idim = 1, ndim
          zii = diag(row_indices(istate),idim)
          zjj = diag(col_indices(jstate),idim)
          gradsq_ij = gradsq_ij + weights ( idim ) * &
                     4.0_dp * REAL( ( CONJG ( zii ) * zii + CONJG ( zjj ) * zjj ), KIND=dp)
        END DO
        matrix % local_data ( istate, jstate ) = gradsq_ij
      END DO
    END DO
  END SUBROUTINE gradsq_at_0
! *****************************************************************************
!> \brief ...
!> \param matrix_p ...
!> \param weights ...
!> \param matrix ...
! *****************************************************************************
  SUBROUTINE grad_at_0 ( matrix_p, weights, matrix )
    TYPE(cp_cfm_p_type)                      :: matrix_p( : )
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_type), POINTER                :: matrix

    COMPLEX(KIND=dp)                         :: zii, zij, zjj
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: diag
    INTEGER                                  :: dim_m, idim, istate, jstate, &
                                                ncol_local, nrow_global, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: grad_ij

    NULLIFY(diag)
    CALL cp_fm_get_info ( matrix, nrow_local = nrow_local,  &
                   ncol_local = ncol_local, nrow_global = nrow_global, &
                   row_indices = row_indices, col_indices = col_indices)
    dim_m = SIZE(matrix_p,1)
    ALLOCATE(diag(nrow_global,dim_m))

    DO idim = 1,dim_m
      DO istate = 1, nrow_global
         CALL cp_cfm_get_element(matrix_p(idim)%matrix,istate,istate,diag(istate,idim))
      ENDDO
    ENDDO

    DO istate = 1, nrow_local
      DO jstate = 1, ncol_local
! get real and imaginary parts
        grad_ij = 0.0_dp
        DO idim = 1, dim_m
          zii = diag(row_indices(istate),idim)
          zjj = diag(col_indices(jstate),idim)
          zij = matrix_p(idim)%matrix %local_data(istate,jstate)
          grad_ij = grad_ij + weights ( idim ) * &
                    REAL ( 4.0_dp * CONJG ( zij ) * ( zjj - zii ), dp )
        END DO
        matrix % local_data ( istate, jstate ) = grad_ij
      END DO
    END DO
    DEALLOCATE(diag)
  END SUBROUTINE grad_at_0

! return energy and maximum gradient in the current point
! *****************************************************************************
!> \brief ...
!> \param weights ...
!> \param zij ...
!> \param tolerance ...
!> \param value ...
! *****************************************************************************
  SUBROUTINE check_tolerance_new( weights, zij, tolerance, value)
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: ZIJ( :, : )
    REAL(KIND=dp)                            :: tolerance, value

    COMPLEX(KIND=dp)                         :: kii, kij, kjj
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: diag
    INTEGER                                  :: idim, istate, jstate, &
                                                ncol_local, nrow_global, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: grad_ij, ra, rb

    NULLIFY(diag)
    CALL cp_fm_get_info ( zij(1,1)%matrix, nrow_local = nrow_local,  &
                   ncol_local = ncol_local, nrow_global = nrow_global, &
                   row_indices = row_indices, col_indices = col_indices)
    ALLOCATE(diag(nrow_global,SIZE(zij,2)))
    value=0.0_dp
    DO idim = 1,SIZE(zij,2)
      DO istate=1, nrow_global
         CALL cp_fm_get_element(zij(1,idim)%matrix,istate,istate,ra)
         CALL cp_fm_get_element(zij(2,idim)%matrix,istate,istate,rb)
         diag(istate,idim)=CMPLX(ra,rb,dp)
         value=value+weights(idim)-weights(idim)*ABS(diag(istate,idim))**2
      ENDDO
    ENDDO
    tolerance=0.0_dp
    DO istate=1,nrow_local
      DO jstate=1,ncol_local
         grad_ij = 0.0_dp
         DO idim = 1, SIZE(zij,2)
            kii = diag(row_indices(istate),idim)
            kjj = diag(col_indices(jstate),idim)
            ra = zij(1,idim) % matrix % local_data(istate,jstate)
            rb = zij(2,idim) % matrix % local_data(istate,jstate)
            kij = CMPLX(ra,rb,dp)
            grad_ij = grad_ij + weights ( idim ) * &
                      REAL ( 4.0_dp * CONJG ( kij ) * ( kjj - kii ), dp )
          END DO
          tolerance=MAX(ABS(grad_ij),tolerance)
       ENDDO
    ENDDO
    CALL mp_max(tolerance,zij(1,1)%matrix%matrix_struct%para_env%group)

    DEALLOCATE(diag)

  END SUBROUTINE check_tolerance_new


! *****************************************************************************
!> \brief yet another crazy try, computes the angles needed to rotate the orbitals first
!>        and rotates them all at the same time (hoping for the best of course)
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param max_iter ...
!> \param max_crazy_angle ...
!> \param crazy_scale ...
!> \param crazy_use_diag ...
!> \param eps_localization ...
!> \param iterations ...
!> \param converged ...
! *****************************************************************************
  SUBROUTINE crazy_rotations( weights, zij, vectors,max_iter, max_crazy_angle, crazy_scale, crazy_use_diag, &
                              eps_localization, iterations, converged)
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: ZIJ( :, : )
    TYPE(cp_fm_type), POINTER                :: vectors
    INTEGER, INTENT(IN)                      :: max_iter
    REAL(KIND=dp), INTENT(IN)                :: max_crazy_angle
    REAL(KIND=dp)                            :: crazy_scale
    LOGICAL                                  :: crazy_use_diag
    REAL(KIND=dp), INTENT(IN)                :: eps_localization
    INTEGER                                  :: iterations
    LOGICAL, INTENT(out), OPTIONAL           :: converged

    CHARACTER(len=*), PARAMETER :: routineN = 'crazy_rotations', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: cone = (1.0_dp,0.0_dp), &
                                                czero = (0.0_dp,0.0_dp)

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: evals_exp
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: diag_z
    COMPLEX(KIND=dp), POINTER                :: mii( : ), mij( : ), mjj( : )
    INTEGER :: dim2, handle, i, icol, idim, irow, istat, method, ncol_global, &
      ncol_local, norder, nrow_global, nrow_local, nsquare, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: do_emd, failure
    REAL(KIND=dp)                            :: eps_exp, limit_crazy_angle, &
                                                maxeval, norm, ra, rb, theta, &
                                                tolerance, value
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals
    TYPE(cp_cfm_type), POINTER               :: cmat_A, cmat_R, cmat_t1
    TYPE(cp_fm_type), POINTER                :: mat_R, mat_t, mat_theta, mat_U

    failure=.FALSE.

    CALL timeset(routineN,handle)
    NULLIFY(row_indices,col_indices)
    NULLIFY(mat_U,mat_t,mat_R)
    NULLIFY(cmat_A,cmat_R,cmat_t1)
    CALL cp_fm_get_info(zij(1,1)%matrix,nrow_global=nrow_global, &
                        ncol_global=ncol_global, &
                        row_indices=row_indices, col_indices=col_indices, &
                        nrow_local=nrow_local, ncol_local=ncol_local)

    limit_crazy_angle=max_crazy_angle

    NULLIFY(diag_z,evals,evals_exp,mii,mij,mjj)
    dim2 = SIZE(zij,2)
    ALLOCATE(diag_z(nrow_global,dim2),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(evals(nrow_global),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(evals_exp(nrow_global),STAT=istat)
    CPASSERT(istat==0)

    CALL cp_cfm_create ( cmat_A, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_R, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_t1, zij ( 1, 1 ) % matrix % matrix_struct)

    CALL cp_fm_create ( mat_U, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create ( mat_t, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create ( mat_R, zij ( 1, 1 ) % matrix % matrix_struct)

    NULLIFY(mat_theta)
    CALL cp_fm_create ( mat_theta, zij ( 1, 1 ) % matrix % matrix_struct)

    CALL cp_fm_set_all( mat_R,0.0_dp,1.0_dp)
    CALL cp_fm_set_all( mat_t, 0.0_dp)
    ALLOCATE(mii(dim2),mij(dim2),mjj(dim2),STAT=istat)
    CPASSERT(istat==0)
    DO idim=1,dim2
       CALL cp_fm_scale_and_add(1.0_dp,mat_t,weights(idim),zij(1,idim) % matrix)
       CALL cp_fm_scale_and_add(1.0_dp,mat_t,weights(idim),zij(2,idim) % matrix)
    ENDDO
    CALL cp_fm_syevd(mat_t,mat_U,evals)
    DO idim=1,dim2
       ! rotate z's
       CALL cp_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(1,idim)%matrix,mat_U,0.0_dp,mat_t)
       CALL cp_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(1,idim)%matrix)
       CALL cp_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(2,idim)%matrix,mat_U,0.0_dp,mat_t)
       CALL cp_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(2,idim)%matrix)
    ENDDO
    ! collect rotation matrix
    CALL cp_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_R,mat_U,0.0_dp,mat_t)
    CALL cp_fm_to_fm(mat_t,mat_R)

    unit_nr=-1
    IF (cmat_A%matrix_struct%para_env%mepos .EQ. cmat_A%matrix_struct%para_env%source ) THEN
        unit_nr=cp_logger_get_default_unit_nr()
        WRITE(unit_nr,'(T2,A7,A6,1X,A20,A12,A12,A12)') &
               "CRAZY| ","Iter","value    ","gradient","Max. eval","limit"
    ENDIF

    iterations=0
    tolerance=1.0_dp

    DO
      iterations=iterations+1
      DO idim=1,dim2
         DO i=1,nrow_global
            CALL cp_fm_get_element(zij(1,idim) % matrix, i, i, ra)
            CALL cp_fm_get_element(zij(2,idim) % matrix, i, i, rb)
            diag_z(i,idim)=CMPLX(ra,rb,dp)
         ENDDO
      ENDDO
      DO irow=1,nrow_local
         DO icol=1,ncol_local
            DO idim=1,dim2
               ra=zij(1,idim) % matrix % local_data(irow,icol)
               rb=zij(2,idim) % matrix % local_data(irow,icol)
               mij(idim)=CMPLX(ra,rb,dp)
               mii(idim)=diag_z(row_indices(irow),idim)
               mjj(idim)=diag_z(col_indices(icol),idim)
            ENDDO
            IF (row_indices(irow).NE.col_indices(icol)) THEN
              CALL get_angle ( mii, mjj, mij, weights, theta )
              theta=crazy_scale*theta
              IF (theta.gt.limit_crazy_angle)  theta=limit_crazy_angle
              IF (theta.lt.-limit_crazy_angle) theta=-limit_crazy_angle
              IF (crazy_use_diag) THEN
               cmat_A % local_data(irow,icol) = -CMPLX(0.0_dp, theta, dp )
              ELSE
               mat_theta % local_data(irow,icol) = - theta
              ENDIF
            ELSE
              IF (crazy_use_diag) THEN
               cmat_A % local_data(irow,icol) = czero
              ELSE
               mat_theta % local_data(irow,icol) = 0.0_dp
              ENDIF
            ENDIF
         ENDDO
      ENDDO

      ! construct rotation matrix U based on A using diagonalization
      ! alternatively, exp based on repeated squaring could be faster
      IF (crazy_use_diag) THEN
         CALL cp_cfm_heevd(cmat_A,cmat_R,evals)
         maxeval=MAXVAL(ABS(evals))
         evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * evals(:) )
         CALL cp_cfm_to_cfm(cmat_R,cmat_t1)
         CALL cp_cfm_column_scale(cmat_t1,evals_exp)
         CALL cp_cfm_gemm('N','C',nrow_global,nrow_global,nrow_global,cone,&
              cmat_t1,cmat_R,czero,cmat_A)
         mat_U%local_data=REAL(cmat_A%local_data,KIND=dp) ! U is a real matrix
      ELSE
         do_emd=.FALSE.
         method=2
         eps_exp=1.0_dp*EPSILON(eps_exp)
         CALL cp_fm_maxabsrownorm(mat_theta,norm)
         maxeval=norm ! an upper bound
         CALL get_nsquare_norder(norm,nsquare,norder,eps_exp,method,do_emd)
         CALL exp_pade_real(mat_U,mat_theta,nsquare,norder)
      ENDIF

      DO idim=1,dim2
         ! rotate z's
         CALL cp_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(1,idim)%matrix,mat_U,0.0_dp,mat_t)
         CALL cp_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(1,idim)%matrix)
         CALL cp_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(2,idim)%matrix,mat_U,0.0_dp,mat_t)
         CALL cp_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(2,idim)%matrix)
      ENDDO
      ! collect rotation matrix
      CALL cp_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_R,mat_U,0.0_dp,mat_t)
      CALL cp_fm_to_fm(mat_t,mat_R)

      CALL check_tolerance_new ( weights, zij, tolerance, value)

      IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T2,A7,I6,1X,G20.15,E12.4,E12.4,E12.4)') &
                 "CRAZY| ",iterations,value,tolerance,maxeval,limit_crazy_angle
          CALL m_flush(unit_nr)
      ENDIF
      IF (tolerance .LT. eps_localization .OR. iterations.ge.max_iter) EXIT
    ENDDO

    IF (PRESENT(converged)) converged=(tolerance .LT. eps_localization)

    CALL cp_cfm_release(cmat_A)
    CALL cp_cfm_release(cmat_R)
    CALL cp_cfm_release(cmat_T1)

    CALL cp_fm_release(mat_U)
    CALL cp_fm_release(mat_T)
    CALL cp_fm_release(mat_theta)

    CALL rotate_orbitals(mat_R,vectors)

    CALL cp_fm_release(mat_R)
    DEALLOCATE(evals_exp,evals,diag_z,STAT=istat)
    DEALLOCATE(mii,mij,mjj,STAT=istat)
    CPASSERT(istat==0)

    CALL timestop(handle)

  END SUBROUTINE crazy_rotations

! *****************************************************************************
!> \brief use the exponential parametrization as described in to perform a direct mini
!>        Gerd Berghold et al. PRB 61 (15), pag. 10040 (2000)
!> none of the input is modified for the time being, just finds the rotations
!> that minimizes, and throws it away afterwards :-)
!> apart from being expensive and not cleaned, this works fine
!> useful to try different spread functionals
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param iterations ...
! *****************************************************************************
  SUBROUTINE direct_mini( weights, zij, vectors, max_iter ,eps_localization, iterations)
    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: ZIJ( :, : )
    TYPE(cp_fm_type), POINTER                :: vectors
    INTEGER, INTENT(IN)                      :: max_iter
    REAL(KIND=dp), INTENT(IN)                :: eps_localization
    INTEGER                                  :: iterations

    CHARACTER(len=*), PARAMETER :: routineN = 'direct_mini', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: cone = (1.0_dp,0.0_dp), &
                                                czero = (0.0_dp,0.0_dp)
    REAL(KIND=dp), PARAMETER                 :: gold_sec = 0.3819_dp

    COMPLEX(KIND=dp)                         :: lk, ll, tmp
    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: evals_exp
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: diag_z
    INTEGER :: handle, i, icol, idim, irow, istat, line_search_count, &
      line_searches, lsl, lsm, lsr, n, ncol_local, ndim, nrow_local, &
      output_unit
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, new_direction
    REAL(KIND=dp) :: a, b, beta_pr, c, denom, ds, ds_min, fa, fb, fc, nom, &
      normg, normg_cross, normg_old, npos, omega, tol, val, x0, x1, xa, xb, xc
    REAL(KIND=dp), DIMENSION(150)            :: energy, grad, pos
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals, fval, fvald
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: c_zij
    TYPE(cp_cfm_type), POINTER               :: cmat_A, cmat_B, cmat_M, &
                                                cmat_R, cmat_t1, cmat_t2, &
                                                cmat_U
    TYPE(cp_fm_type), POINTER                :: matrix_A, matrix_G, &
                                                matrix_G_old, &
                                                matrix_G_search, matrix_H, &
                                                matrix_R, matrix_T
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(evals,evals_exp,diag_z,fval,fvald,c_zij)
    NULLIFY(matrix_A,matrix_G,matrix_T,matrix_G_search,matrix_G_old)
    NULLIFY(cmat_A,cmat_U,cmat_R,cmat_t1,cmat_t2,cmat_B,cmat_M)

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    output_unit= cp_logger_get_default_io_unit(logger)

    n = zij(1,1)%matrix%matrix_struct%nrow_global
    ndim = (SIZE(zij,2))

    failure=.FALSE.

    IF ( output_unit > 0) THEN
      WRITE(output_unit,'(T10,A )') "Localization by direct minimization of the functional; "
      WRITE(output_unit,'(T5,2A13,A20,A20,A10 )') " Line search ", " Iteration ", " Functional ", " Tolerance "," ds Min "
    END IF

    ALLOCATE(evals(n),evals_exp(n),diag_z(n,ndim),fval(n),fvald(n))
    ALLOCATE(c_zij(ndim), STAT=istat)
    CPASSERT(istat==0)

    ! create the three complex matrices Z
    DO idim=1,ndim
       NULLIFY(c_zij(idim)%matrix)
       CALL cp_cfm_create ( c_zij(idim)%matrix, zij ( 1, 1 ) % matrix % matrix_struct)
       c_zij(idim)%matrix% local_data = CMPLX (zij(1,idim) % matrix % local_data, &
                    zij(2,idim) % matrix % local_data, dp )
    ENDDO

    CALL cp_fm_create  ( matrix_A, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create  ( matrix_G, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create  ( matrix_T, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create  ( matrix_H, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create  ( matrix_G_search, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create  ( matrix_G_old, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create  ( matrix_R, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_set_all(matrix_R,0.0_dp,1.0_dp)

    CALL cp_fm_set_all(matrix_A,0.0_dp)
!    CALL cp_fm_init_random ( matrix_A )

    CALL cp_cfm_create ( cmat_A, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_U, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_R, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_t1, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_t2, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_B, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_M, zij ( 1, 1 ) % matrix % matrix_struct)

    CALL cp_cfm_get_info ( cmat_B, nrow_local=nrow_local, ncol_local=ncol_local, &
                                   row_indices=row_indices, col_indices=col_indices)

    CALL cp_fm_set_all(matrix_G_old,0.0_dp)
    CALL cp_fm_set_all(matrix_G_search,0.0_dp)
    normg_old=1.0E30_dp
    ds_min=1.0_dp
    new_direction=.TRUE.
    Iterations=0
    line_searches=0
    line_search_count=0
    DO
         iterations = iterations + 1
         ! compute U,R,evals given A
         cmat_A % local_data = CMPLX ( 0.0_dp , matrix_A % local_data, dp ) ! cmat_A is hermitian, evals are reals
         CALL cp_cfm_heevd(cmat_A,cmat_R,evals)
         evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * evals(:) )
         CALL cp_cfm_to_cfm(cmat_R,cmat_t1)
         CALL cp_cfm_column_scale(cmat_t1,evals_exp)
         CALL cp_cfm_gemm('N','C',n,n,n,cone,cmat_t1,cmat_R,czero,cmat_U)
         cmat_U%local_data=REAL(cmat_U%local_data,KIND=dp) ! enforce numerics, U is a real matrix

         IF ( new_direction .AND. MOD(line_searches,20).EQ.5 ) THEN ! reset with A .eq. 0
            DO idim=1,ndim
               CALL cp_cfm_gemm('N','N',n,n,n,cone,c_zij(idim)%matrix,cmat_U,czero,cmat_t1)
               CALL cp_cfm_gemm('C','N',n,n,n,cone,cmat_U,cmat_t1,czero,c_zij(idim)%matrix)
            ENDDO
            ! collect rotation matrix
            matrix_H%local_data=REAL(cmat_U%local_data,KIND=dp)
            CALL cp_gemm('N','N',n,n,n,1.0_dp,matrix_R,matrix_H,0.0_dp,matrix_T)
            CALL cp_fm_to_fm(matrix_T,matrix_R)

            CALL cp_cfm_set_all(cmat_U,czero,cone)
            CALL cp_cfm_set_all(cmat_R,czero,cone)
            CALL cp_cfm_set_all(cmat_A,czero)
            CALL cp_fm_set_all(matrix_A,0.0_dp)
            evals(:)=0.0_dp
            evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * evals(:) )
            CALL cp_fm_set_all(matrix_G_old,0.0_dp)
            CALL cp_fm_set_all(matrix_G_search,0.0_dp)
            normg_old=1.0E30_dp
         ENDIF

         ! compute Omega and M
         CALL cp_cfm_set_all(cmat_M,czero)
         omega = 0.0_dp
         DO idim=1,ndim
            CALL cp_cfm_gemm('N','N',n,n,n,cone,c_zij(idim)%matrix,cmat_U,czero,cmat_t1) ! t1=ZU
            CALL cp_cfm_gemm('C','N',n,n,n,cone,cmat_U,cmat_t1,czero,cmat_t2) ! t2=(U^T)ZU
            DO i=1,n
               CALL cp_cfm_get_element(cmat_t2,i,i,diag_z(i,idim))
               SELECT CASE (2) ! allows for selection of different spread functionals
               CASE (1)
                 fval(i) =-weights(idim)*LOG(ABS(diag_z(i,idim))**2)
                 fvald(i)=-weights(idim)/(ABS(diag_z(i,idim))**2)
               CASE (2) ! corresponds to the jacobi setup
                 fval(i) =weights(idim)-weights(idim)*ABS(diag_z(i,idim))**2
                 fvald(i)=-weights(idim)
               END SELECT
               omega=omega+fval(i)
            ENDDO
            DO icol=1,ncol_local
               DO irow=1,nrow_local
                  tmp = cmat_t1%local_data(irow,icol)*CONJG(diag_z(col_indices(icol),idim))
                  cmat_M%local_data(irow,icol)=cmat_M%local_data(irow,icol) &
                                              +4.0_dp*fvald(col_indices(icol))*REAL(tmp,KIND=dp)
               ENDDO
            ENDDO
         ENDDO

         ! compute Hessian diagonal approximation for the preconditioner
         IF (.TRUE.) THEN
            CALL gradsq_at_0(diag_z,weights,matrix_H,ndim)
         ELSE
            CALL cp_fm_set_all(matrix_H,1.0_dp)
         ENDIF

         ! compute B
         DO icol=1,ncol_local
          DO irow=1,nrow_local
             ll=(0.0_dp,-1.0_dp)*evals(row_indices(irow))
             lk=(0.0_dp,-1.0_dp)*evals(col_indices(icol))
             IF (ABS(ll-lk).lt.0.5_dp) THEN ! use a series expansion to avoid loss of precision
                tmp=1.0_dp
                cmat_B%local_data(irow,icol)=0.0_dp
                DO i=1,16
                   cmat_B%local_data(irow,icol)=cmat_B%local_data(irow,icol)+tmp
                   tmp=tmp*(ll-lk)/(i+1)
                ENDDO
                cmat_B%local_data(irow,icol)=cmat_B%local_data(irow,icol)*EXP(lk)
             ELSE
                cmat_B%local_data(irow,icol)=(EXP(lk)-EXP(ll))/(lk-ll)
             ENDIF
          ENDDO
         ENDDO
         ! compute gradient matrix_G

         CALL cp_cfm_gemm('C','N',n,n,n,cone,cmat_M,cmat_R,czero,cmat_t1) ! t1=(M^T)(R^T)
         CALL cp_cfm_gemm('C','N',n,n,n,cone,cmat_R,cmat_t1,czero,cmat_t2) ! t2=(R)t1
         CALL cp_cfm_schur_product(cmat_t2,cmat_B,cmat_t1)
         CALL cp_cfm_gemm('N','C',n,n,n,cone,cmat_t1,cmat_R,czero,cmat_t2)
         CALL cp_cfm_gemm('N','N',n,n,n,cone,cmat_R,cmat_t2,czero,cmat_t1)
         matrix_G%local_data=REAL(cmat_t1%local_data,KIND=dp)
         CALL cp_fm_transpose(matrix_G,matrix_T)
         CALL cp_fm_scale_and_add(  -1.0_dp,matrix_G, 1.0_dp,matrix_T)
         CALL cp_fm_maxabsval(matrix_G,tol)

         ! from here on, minimizing technology
         IF ( new_direction ) THEN
            ! energy converged up to machine precision ?
            line_searches=line_searches+1
            ! DO i=1,line_search_count
            !   write(15,*) pos(i),energy(i)
            ! ENDDO
            ! write(15,*) ""
            ! CALL m_flush(15)
            !write(16,*) evals(:)
            !write(17,*) matrix_A%local_data(:,:)
            !write(18,*) matrix_G%local_data(:,:)
            IF (output_unit > 0) THEN
               WRITE(output_unit,'(T5,I10,T18,I10,T31,2F20.6,F10.3)') line_searches,Iterations,Omega,tol,ds_min
              CALL m_flush(output_unit)
            ENDIF
            IF (tol<eps_localization .OR. iterations>max_iter)  EXIT

            IF (.TRUE.) THEN ! do conjugate gradient CG
                CALL cp_fm_trace(matrix_G,matrix_G_old,normg_cross)
                normg_cross=normg_cross*0.5_dp ! takes into account the fact that A is antisymmetric
                ! apply the preconditioner
                DO icol=1,ncol_local
                   DO irow=1,nrow_local
                      matrix_G_old%local_data(irow,icol)=matrix_G%local_data(irow,icol)/matrix_H%local_data(irow,icol)
                   ENDDO
                ENDDO
                CALL cp_fm_trace(matrix_G,matrix_G_old,normg)
                normg=normg*0.5_dp
                beta_pr=(normg-normg_cross)/normg_old
                normg_old=normg
                beta_pr=MAX(beta_pr,0.0_dp)
                CALL cp_fm_scale_and_add(beta_pr,matrix_G_search,-1.0_dp,matrix_G_old)
                CALL cp_fm_trace(matrix_G_search,matrix_G_old,normg_cross)
                IF (normg_cross .GE. 0) THEN ! back to SD
                   IF (matrix_A%matrix_struct%para_env%mepos .EQ. &
                       matrix_A%matrix_struct%para_env%source) THEN
                      WRITE(cp_logger_get_default_unit_nr(),*) "!"
                   ENDIF
                   beta_pr=0.0_dp
                   CALL cp_fm_scale_and_add(beta_pr,matrix_G_search,-1.0_dp,matrix_G_old)
                ENDIF
            ELSE ! SD
                CALL cp_fm_scale_and_add(0.0_dp,matrix_G_search,-1.0_dp,matrix_G)
            ENDIF
            ! ds_min=1.0E-4_dp
            line_search_count=0
         END IF
         line_search_count=line_search_count+1
         energy(line_search_count)=Omega

         ! line search section
         SELECT CASE (3)
         CASE(1) ! two point line search
           SELECT CASE (line_search_count)
           CASE (1)
              pos(1)=0.0_dp
              pos(2)=ds_min
              CALL cp_fm_trace(matrix_G,matrix_G_search,grad(1))
              grad(1)=grad(1)/2.0_dp
              new_direction=.FALSE.
           CASE (2)
              new_direction=.TRUE.
              x0=pos(1) ! 0.0_dp
              c=energy(1)
              b=grad(1)
              x1=pos(2)
              a=(energy(2)-b*x1-c)/(x1**2)
              IF (a.le.0.0_dp) a=1.0E-15_dp
              npos=-b/(2.0_dp*a)
              val=a*npos**2+b*npos+c
              IF (val.lt.energy(1) .AND. val.le.energy(2)) THEN
                 ! we go to a minimum, but ...
                 ! we take a guard against too large steps
                 pos(3)=MIN(npos,MAXVAL(pos(1:2))*4.0_dp)
              ELSE  ! just take an extended step
                 pos(3)=MAXVAL(pos(1:2))*2.0_dp
              ENDIF
           END SELECT
         CASE(2) ! 3 point line search
              SELECT CASE(line_search_count)
              CASE(1)
                new_direction=.FALSE.
                pos(1)=0.0_dp
                pos(2)=ds_min*0.8_dp
              CASE(2)
                new_direction=.FALSE.
                IF (energy(2).gt.energy(1)) THEN
                   pos(3)=ds_min*0.7_dp
                ELSE
                   pos(3)=ds_min*1.4_dp
                ENDIF
              CASE(3)
                new_direction=.TRUE.
                xa=pos(1)
                xb=pos(2)
                xc=pos(3)
                fa=energy(1)
                fb=energy(2)
                fc=energy(3)
                nom  =(xb-xa)**2*(fb-fc) -  (xb-xc)**2*(fb-fa)
                denom=(xb-xa)*(fb-fc) -  (xb-xc)*(fb-fa)
                IF (ABS(denom) .LE. 1.0E-18_dp*MAX(ABS(fb-fc),ABS(fb-fa))) THEN
                    npos = xb
                ELSE
                    npos = xb-0.5_dp*nom/denom ! position of the stationary point
                ENDIF
                val = (npos-xa)*(npos-xb)*fc/((xc-xa)*(xc-xb))+ &
                      (npos-xb)*(npos-xc)*fa/((xa-xb)*(xa-xc))+ &
                      (npos-xc)*(npos-xa)*fb/((xb-xc)*(xb-xa))
                IF (val.lt.fa .AND. val.le.fb .AND. val.le.fc) THEN ! OK, we go to a minimum
                    ! we take a guard against too large steps
                    pos(4)=MAX(MAXVAL(pos(1:3))*0.01_dp, &
                                                 MIN(npos,MAXVAL(pos(1:3))*4.0_dp))
                ELSE  ! just take an extended step
                    pos(4)=MAXVAL(pos(1:3))*2.0_dp
                ENDIF
              END SELECT
         CASE(3) ! golden section hunt
            new_direction=.FALSE.
            IF (line_search_count.eq.1) THEN
              lsl=1
              lsr=0
              lsm=1
              pos(1)=0.0_dp
              pos(2)=ds_min/gold_sec
            ELSE
              IF (line_search_count .EQ. 150) STOP "Too many"
              IF (lsr.eq.0) THEN
                 IF (energy(line_search_count-1).lt.energy(line_search_count)) THEN
                    lsr = line_search_count
                    pos(line_search_count+1)=pos(lsm)+(pos(lsr)-pos(lsm))*gold_sec
                 ELSE
                    lsl = lsm
                    lsm = line_search_count
                    pos(line_search_count+1)=pos(line_search_count) / gold_sec
                 ENDIF
              ELSE
                 IF (pos(line_search_count) .LT. pos(lsm)) THEN
                    IF (energy(line_search_count).LT. energy(lsm)) THEN
                       lsr = lsm
                       lsm = line_search_count
                    ELSE
                       lsl = line_search_count
                    ENDIF
                 ELSE
                    IF (energy(line_search_count).LT. energy(lsm)) THEN
                       lsl = lsm
                       lsm = line_search_count
                    ELSE
                       lsr = line_search_count
                    ENDIF
                 ENDIF
                 IF ( pos(lsr)-pos(lsm) .GT. pos(lsm)-pos(lsl)) THEN
                    pos(line_search_count+1) = pos(lsm)+gold_sec*(pos(lsr)-pos(lsm))
                 ELSE
                    pos(line_search_count+1) = pos(lsl)+gold_sec*(pos(lsm)-pos(lsl))
                 ENDIF
                 IF ((pos(lsr)-pos(lsl)) .LT. 1.0E-3_dp*pos(lsr)) THEN
                    new_direction=.TRUE.
                 ENDIF
              ENDIF ! lsr .eq. 0
            ENDIF ! first step
         END SELECT
         ! now go to the suggested point
         ds_min=pos(line_search_count+1)
         ds=pos(line_search_count+1)-pos(line_search_count)
         CALL cp_fm_scale_and_add(1.0_dp,matrix_A,ds,matrix_G_search)
    ENDDO

    ! collect rotation matrix
    matrix_H%local_data=REAL(cmat_U%local_data,KIND=dp)
    CALL cp_gemm('N','N',n,n,n,1.0_dp,matrix_R,matrix_H,0.0_dp,matrix_T)
    CALL cp_fm_to_fm(matrix_T,matrix_R)
    CALL rotate_orbitals(matrix_R,vectors)
    CALL cp_fm_release  ( matrix_R)

    CALL cp_fm_release  ( matrix_A)
    CALL cp_fm_release  ( matrix_G)
    CALL cp_fm_release  ( matrix_H)
    CALL cp_fm_release  ( matrix_T)
    CALL cp_fm_release  ( matrix_G_search)
    CALL cp_fm_release  ( matrix_G_old)
    CALL cp_cfm_release ( cmat_A)
    CALL cp_cfm_release ( cmat_U)
    CALL cp_cfm_release ( cmat_R)
    CALL cp_cfm_release ( cmat_t1)
    CALL cp_cfm_release ( cmat_t2)
    CALL cp_cfm_release ( cmat_B)
    CALL cp_cfm_release ( cmat_M)

    DEALLOCATE(evals,evals_exp,fval,fvald)

    DO idim=1,SIZE(c_zij)
       zij(1,idim) % matrix % local_data =  REAL ( c_zij(idim)%matrix% local_data , dp )
       zij(2,idim) % matrix % local_data = AIMAG ( c_zij(idim)%matrix% local_data )
       CALL cp_cfm_release( c_zij(idim)%matrix)
    ENDDO
    DEALLOCATE(c_zij)
    DEALLOCATE(diag_z)

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief Parallel algorithm for jacobi rotations
!> \param weights ...
!> \param zij ...
!> \param vectors ...
!> \param para_env ...
!> \param max_iter ...
!> \param eps_localization ...
!> \param sweeps ...
!> \param out_each ...
!> \param target_time ...
!> \param start_time ...
!> \par History
!>      use allgather for improved performance
!> \author MI (11.2009)
! *****************************************************************************
  SUBROUTINE jacobi_rot_para(weights, zij, vectors, para_env, max_iter, eps_localization, &
              sweeps, out_each, target_time, start_time)

    REAL(KIND=dp), INTENT(IN)                :: weights( : )
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: ZIJ( :, : )
    TYPE(cp_fm_type), POINTER                :: vectors
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: max_iter
    REAL(KIND=dp), INTENT(IN)                :: eps_localization
    INTEGER                                  :: sweeps
    INTEGER, INTENT(IN)                      :: out_each
    REAL(dp)                                 :: target_time, start_time

    CHARACTER(len=*), PARAMETER :: routineN = 'jacobi_rot_para', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: zi, zj
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: c_array_me, c_array_partner
    COMPLEX(KIND=dp), POINTER                :: mii( : ), mij( : ), mjj( : )
    INTEGER :: dim2, handle, i, idim, ii, ik, il1, il2, il_recv, &
      il_recv_partner, ilow1, ilow2, ip, ip_has_i, ip_partner, ip_recv_from, &
      ip_recv_partner, ipair, iperm, istat, istate, iu1, iu2, iup1, iup2, j, &
      jj, jstate, k, kk, n1, n2, nblock, nblock_max, npair, nperm, ns_me, &
      ns_partner, ns_recv_from, ns_recv_partner, nstate, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: list_pair, ns_bound
    LOGICAL                                  :: failure, should_stop
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: gmat, rmat_loc, rmat_recv, &
                                                rmat_send, rotmat, &
                                                z_ij_loc_im, z_ij_loc_re
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rmat_recv_all
    REAL(KIND=dp)                            :: ct, func, gmax, grad, ri, rj, &
                                                st, t1, t2, theta, tolerance, &
                                                xlow, xstate, xup, zc, zr
    TYPE(cp_fm_type), POINTER                :: rmat
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(set_c_1d_type), DIMENSION(:), &
      POINTER                                :: zdiag_all, zdiag_me
    TYPE(set_c_2d_type), DIMENSION(:), &
      POINTER                                :: cz_ij_loc, xyz_mix, xyz_mix_ns

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(logger)
    logger => cp_get_default_logger()
    output_unit= cp_logger_get_default_io_unit(logger)

    NULLIFY(rmat, cz_ij_loc, zdiag_all, zdiag_me)
    NULLIFY(xyz_mix,xyz_mix_ns)
    NULLIFY(mii,mij,mjj)

    dim2 = SIZE(zij,2)
    ALLOCATE(mii(dim2),mij(dim2),mjj(dim2),STAT=istat)
    CPASSERT(istat==0)

    CALL cp_fm_create ( rmat, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_set_all ( rmat, 0._dp, 1._dp)

    CALL cp_fm_get_info ( rmat , nrow_global = nstate)

    ALLOCATE(rcount(para_env%num_pe),STAT=istat)
    ALLOCATE(rdispl(para_env%num_pe),STAT=istat)

    tolerance = 1.0e10_dp
    sweeps = 0

    ! number of processor pairs and number of permutations
    npair = (para_env%num_pe+1)/2
    nperm = para_env%num_pe - MOD( para_env%num_pe +1,2)
    ALLOCATE(list_pair(2,npair),STAT=istat)
    CPASSERT(istat==0)

    ! Distribution of the states (XXXXX safe against more pe than states ??? XXXXX)
    xstate = REAL(nstate,dp)/REAL(para_env%num_pe,dp)
    nblock_max = 0
    ALLOCATE(ns_bound(0:para_env%num_pe-1,2),STAT=istat)
    CPASSERT(istat==0)
    Xlow=0.0D0
    Xup=0.0D0
    DO ip=1,para_env%num_pe
      xup = xlow + xstate
      ns_bound(ip-1,1)=NINT(xlow)+1
      ns_bound(ip-1,2)=NINT(xup)
      IF(NINT(xup).GT.nstate) THEN
         ns_bound(ip-1,2)=nstate
      ENDIF
      IF(NINT(xlow).GT.nstate) THEN
        ns_bound(ip-1,1)=nstate+1
      ENDIF
      xlow = xup
    ENDDO
    DO ip=0,para_env%num_pe-1
      nblock=ns_bound(ip,2)-ns_bound(ip,1)+1
      nblock_max=MAX(nblock_max,nblock)
    ENDDO

    ! otbtain local part of the matrix (could be made faster, but is likely irrelevant).
    ALLOCATE(z_ij_loc_re(nstate,nblock_max),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(z_ij_loc_im(nstate,nblock_max),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(cz_ij_loc(dim2),STAT=istat)
    CPASSERT(istat==0)
    DO idim = 1,dim2
      DO ip=0,para_env%num_pe-1
        nblock=ns_bound(ip,2)-ns_bound(ip,1)+1
        CALL cp_fm_get_submatrix(zij(1,idim)%matrix,z_ij_loc_re,1,ns_bound(ip,1),nstate,nblock)
        CALL cp_fm_get_submatrix(zij(2,idim)%matrix,z_ij_loc_im,1,ns_bound(ip,1),nstate,nblock)
        IF(para_env%mepos==ip) THEN
          ns_me = nblock
          ALLOCATE (cz_ij_loc(idim)%c_array(nstate,ns_me),STAT=istat)
          CPASSERT(istat==0)
          DO i=1,ns_me
            DO j = 1,nstate
              cz_ij_loc(idim)%c_array(j,i) = CMPLX(z_ij_loc_re(j,i),z_ij_loc_im(j,i),dp)
            END DO
          END DO
        END IF
      END DO ! ip
    END DO
    DEALLOCATE(z_ij_loc_re,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(z_ij_loc_im,STAT=istat)
    CPASSERT(istat==0)


    ! initialize rotation matrix
    ALLOCATE(rotmat(nstate,2*nblock_max),STAT=istat)
    CPASSERT(istat==0)
    rotmat=0.0_dp
    DO i =ns_bound(para_env%mepos,1),ns_bound(para_env%mepos,2)
       ii = i - ns_bound(para_env%mepos,1) + 1
       rotmat(i,ii) = 1.0_dp
    END DO

    ALLOCATE(xyz_mix(dim2),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(xyz_mix_ns(dim2),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(zdiag_me(dim2),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(zdiag_all(dim2),STAT=istat)
    CPASSERT(istat==0)

    ns_me = ns_bound(para_env%mepos,2) - ns_bound(para_env%mepos,1) +1
    IF(ns_me/=0) THEN
      ALLOCATE(c_array_me(nstate,ns_me,dim2),STAT=istat)
      CPASSERT(istat==0)
      DO idim = 1,dim2
         ALLOCATE(xyz_mix_ns(idim)%c_array(nstate,ns_me),STAT=istat)
         CPASSERT(istat==0)
      END DO
      ALLOCATE(gmat(nstate,ns_me),STAT=istat)
      CPASSERT(istat==0)
    END IF

    DO idim = 1,dim2
       ALLOCATE(zdiag_me(idim)%c_array(nblock_max),STAT=istat)
       CPASSERT(istat==0)
       zdiag_me(idim)%c_array = CMPLX(0.0_dp, 0.0_dp, dp)
       ALLOCATE(zdiag_all(idim)%c_array(para_env%num_pe*nblock_max),STAT=istat)
       CPASSERT(istat==0)
       zdiag_all(idim)%c_array = CMPLX(0.0_dp, 0.0_dp, dp)
    END DO
    ALLOCATE(rmat_recv(nblock_max*2,nblock_max),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(rmat_send(nblock_max*2,nblock_max),STAT=istat)
    CPASSERT(istat==0)

    ! buffer for message passing
    ALLOCATE(rmat_recv_all(nblock_max*2,nblock_max,0:para_env%num_pe-1),STAT=istat)
    CPASSERT(istat==0)

    IF(output_unit>0) THEN
      WRITE(output_unit,'(T10,A )') "Localization by iterative distributed Jacobi rotation"
      WRITE(output_unit,'(T20,A12,T32, A22,T60, A12,A8 )') "Iteration", "Functional", "Tolerance"," Time "
    END IF

    DO sweeps = 1,max_iter+1
      IF(sweeps==max_iter+1) THEN
        IF(output_unit>0) THEN
          WRITE(output_unit,*) ' LOCALIZATION! loop did not converge within the maximum number of iterations.'
          WRITE(output_unit,*) '               Present  Max. gradient = ', tolerance
        END IF
        EXIT
      ENDIF
      t1= m_walltime()

      DO iperm = 1,nperm

        ! fix partners for this permutation, and get the number of states
        CALL eberlein(iperm,para_env,list_pair)
        ip_partner = -1
        ns_partner = 0
        DO ipair = 1,npair
          IF(list_pair(1,ipair) == para_env%mepos) THEN
             ip_partner = list_pair(2,ipair)
             EXIT
          ELSE IF(list_pair(2,ipair) == para_env%mepos) THEN
             ip_partner = list_pair(1,ipair)
             EXIT
          END IF
        END DO
        IF(ip_partner >=0) THEN
           ns_partner =  ns_bound(ip_partner,2) - ns_bound(ip_partner,1) +1
        ELSE
           ns_partner = 0
        END IF

        ! if there is a non-zero block connecting two partners, jacobi-sweep it.
        IF(ns_partner*ns_me /= 0) THEN

          ALLOCATE(rmat_loc(ns_me+ns_partner,ns_me+ns_partner),STAT=istat)
          CPASSERT(istat==0)
          rmat_loc = 0.0_dp
          DO i = 1, ns_me+ns_partner
            rmat_loc(i,i) = 1.0_dp
          END DO

          ALLOCATE(c_array_partner(nstate,ns_partner,dim2),STAT=istat)
          CPASSERT(istat==0)

          DO idim = 1,dim2
            ALLOCATE(xyz_mix(idim)%c_array(ns_me+ns_partner,ns_me+ns_partner),STAT=istat)
            CPASSERT(istat==0)
            DO i = 1,ns_me
               c_array_me(1:nstate,i,idim) = cz_ij_loc(idim)%c_array(1:nstate,i)
            END DO
          END DO

          CALL mp_sendrecv(msgin=c_array_me,dest=ip_partner,&
                msgout=c_array_partner,source=ip_partner,comm=para_env%group)

          n1=ns_me
          n2=ns_partner
          ilow1=ns_bound(para_env%mepos,1)
          iup1 = ns_bound(para_env%mepos,1)+n1-1
          ilow2 = ns_bound(ip_partner,1)
          iup2 = ns_bound(ip_partner,1)+n2-1
          IF(ns_bound(para_env%mepos,1)<ns_bound(ip_partner,1)) THEN
           il1=1
           iu1=n1
           iu1=n1
           il2=1+n1
           iu2=n1+n2
          ELSE
           il1=1+n2
           iu1=n1+n2
           iu1=n1+n2
           il2=1
           iu2=n2
          END IF

          DO idim= 1,dim2
             DO i = 1,n1
               xyz_mix(idim)%c_array(il1:iu1,il1+i-1) = c_array_me(ilow1:iup1,i,idim)
               xyz_mix(idim)%c_array(il2:iu2,il1+i-1) = c_array_me(ilow2:iup2,i,idim)
             END DO
             DO i = 1,n2
               xyz_mix(idim)%c_array(il2:iu2,il2+i-1) = c_array_partner(ilow2:iup2,i,idim)
               xyz_mix(idim)%c_array(il1:iu1,il2+i-1) = c_array_partner(ilow1:iup1,i,idim)
             END DO
          END DO

          DO istate = 1,n1+n2
            DO jstate = istate+1,n1+n2
              DO idim = 1,dim2
                mii(idim) = xyz_mix(idim)%c_array(istate,istate)
                mij(idim) = xyz_mix(idim)%c_array(istate,jstate)
                mjj(idim) = xyz_mix(idim)%c_array(jstate,jstate)
              END DO
              CALL get_angle ( mii, mjj, mij, weights, theta )
              st = SIN ( theta )
              ct = COS ( theta )
              DO idim = 1,dim2
                DO i=1,n1+n2
                  zi = ct*xyz_mix(idim)%c_array(i,istate)+st*xyz_mix(idim)%c_array(i,jstate)
                  zj = -st*xyz_mix(idim)%c_array(i,istate)+ct*xyz_mix(idim)%c_array(i,jstate)
                  xyz_mix(idim)%c_array(i,istate) = zi
                  xyz_mix(idim)%c_array(i,jstate) = zj
                END DO
                DO i=1,n1+n2
                  zi = ct*xyz_mix(idim)%c_array(istate,i) + st*xyz_mix(idim)%c_array(jstate,i)
                  zj = -st*xyz_mix(idim)%c_array(istate,i) + ct*xyz_mix(idim)%c_array(jstate,i)
                  xyz_mix(idim)%c_array(istate,i) = zi
                  xyz_mix(idim)%c_array(jstate,i) = zj
                END DO
              END DO

              DO i = 1,n1+n2
                ri= ct*rmat_loc(i,istate)+st*rmat_loc(i,jstate)
                rj= ct*rmat_loc(i,jstate)-st*rmat_loc(i,istate)
                rmat_loc(i,istate)=ri
                rmat_loc(i,jstate)=rj
              END DO
            END DO
          END DO

          k=nblock_max+1
          CALL mp_sendrecv(rotmat(1:nstate,1:ns_me),ip_partner,&
               rotmat(1:nstate,k:k+n2-1),ip_partner,para_env%group)

          IF(ilow1<ilow2) THEN
            CALL dgemm("N","N",nstate,n1,n2,1.0_dp,rotmat(1,k),nstate,rmat_loc(1+n1,1),n1+n2,0.0_dp,gmat,nstate)
            CALL dgemm("N","N",nstate,n1,n1,1.0_dp,rotmat(1,1),nstate,rmat_loc(1,1),n1+n2,1.0_dp,gmat,nstate)
          ELSE
            CALL dgemm("N","N",nstate,n1,n2,1.0_dp,rotmat(1,k),nstate,rmat_loc(1,n2+1),n1+n2,0.0_dp,gmat,nstate)
            CALL dgemm("N","N",nstate,n1,n1,1.0_dp,rotmat(1,1),nstate,rmat_loc(n2+1,n2+1),n1+n2,1.0_dp,gmat,nstate)
          END IF

          CALL dcopy(nstate*n1,gmat(1,1),1,rotmat(1,1),1)

          DO idim = 1,dim2
             DO i = 1,n1
               xyz_mix_ns(idim)%c_array(1:nstate,i) = CMPLX(0.0_dp,0.0_dp,dp)
             END DO

             DO istate = 1,n1
               DO jstate = 1,nstate
                  DO i = 1,n2
                     xyz_mix_ns(idim)%c_array(jstate,istate) = xyz_mix_ns(idim)%c_array(jstate,istate) +&
                               c_array_partner(jstate,i,idim)*rmat_loc(il2+i-1,il1+istate-1)
                  END DO
               END DO
             END DO
             DO istate = 1,n1
               DO jstate = 1,nstate
                  DO i = 1,n1
                     xyz_mix_ns(idim)%c_array(jstate,istate) = xyz_mix_ns(idim)%c_array(jstate,istate) +&
                               c_array_me(jstate,i,idim)*rmat_loc(il1+i-1,il1+istate-1)
                  END DO
               END DO
             END DO
          END DO  ! idim

          DEALLOCATE(c_array_partner,STAT=istat)
          CPASSERT(istat==0)

        ELSE ! save my data
          DO idim = 1,dim2
             DO i = 1,ns_me
               xyz_mix_ns(idim)%c_array(1:nstate,i) = cz_ij_loc(idim)%c_array(1:nstate,i)
             END DO
          END DO
        END IF

        DO idim=1,dim2
           DO i = 1,ns_me
             cz_ij_loc(idim)%c_array(1:nstate,i) = CMPLX(0.0_dp,0.0_dp,dp)
           END DO
        END DO

        IF(ns_partner*ns_me /= 0) THEN
          ! transpose rotation matrix rmat_loc
          DO i =1,ns_me+ns_partner
            DO j = i+1,ns_me+ns_partner
              ri = rmat_loc(i,j)
              rmat_loc(i,j)=rmat_loc(j,i)
              rmat_loc(j,i)=ri
            END DO
          END DO

          ! prepare for distribution
          DO i = 1,n1
             rmat_send(1:n1,i) = rmat_loc(il1:iu1,il1+i-1)
          END DO
          ik = nblock_max
          DO i = 1,n2
             rmat_send(ik+1:ik+n1,i) = rmat_loc(il1:iu1,il2+i-1)
          END DO
        ELSE
          rmat_send = 0.0_dp
        END IF

        ! collect data from all tasks (this takes some significant time)
        CALL mp_allgather(rmat_send,rmat_recv_all,para_env%group)

        ! update blocks everywhere
        DO ip = 0,para_env%num_pe-1

          ip_recv_from=MOD(para_env%mepos-IP+para_env%num_pe,para_env%num_pe)
          rmat_recv(:,:)=rmat_recv_all(:,:,ip_recv_from)

          ns_recv_from = ns_bound(ip_recv_from,2)-ns_bound(ip_recv_from,1)+1

          IF(ns_me/=0) THEN
            IF(ns_recv_from/=0) THEN
              !look for the partner of ip_recv_from
              ip_recv_partner = -1
              ns_recv_partner = 0
              DO ipair = 1,npair
                IF(list_pair(1,ipair) == ip_recv_from) THEN
                  ip_recv_partner = list_pair(2,ipair)
                  EXIT
                ELSE IF(list_pair(2,ipair) == ip_recv_from) THEN
                  ip_recv_partner = list_pair(1,ipair)
                  EXIT
                END IF
              END DO

              IF(ip_recv_partner>=0) THEN
                ns_recv_partner =  ns_bound(ip_recv_partner,2)-ns_bound(ip_recv_partner,1)+1
              END IF
              IF(ns_recv_partner>0) THEN
                il1 =  ns_bound(para_env%mepos,1)
                il_recv = ns_bound(ip_recv_from,1)
                il_recv_partner = ns_bound(ip_recv_partner,1)
                ik = nblock_max

                DO idim = 1,dim2
                  DO i = 1,ns_recv_from
                    ii = il_recv+i-1
                    DO j = 1,ns_me
                      jj=j
                      DO k = 1,ns_recv_from
                        kk = il_recv+k-1
                        cz_ij_loc(idim)%c_array(ii,jj) = cz_ij_loc(idim)%c_array(ii,jj) + &
                                  rmat_recv(i,k) * xyz_mix_ns(idim)%c_array(kk,j)
                      END DO
                    END DO
                  END DO
                  DO i = 1,ns_recv_from
                    ii = il_recv+i-1
                    DO j = 1,ns_me
                      jj=j
                      DO k = 1,ns_recv_partner
                        kk = il_recv_partner+k-1
                        cz_ij_loc(idim)%c_array(ii,jj) = cz_ij_loc(idim)%c_array(ii,jj) + &
                                  rmat_recv(ik+i,k) * xyz_mix_ns(idim)%c_array(kk,j)
                      END DO
                    END DO
                  END DO
                END DO  ! idim
              ELSE
                il1 =  ns_bound(para_env%mepos,1)
                il_recv = ns_bound(ip_recv_from,1)
                DO idim = 1,dim2
                   DO j = 1,ns_me
                     jj=j
                     DO i = 1,ns_recv_from
                       ii = il_recv+i-1
                       cz_ij_loc(idim)%c_array(ii,jj) = xyz_mix_ns(idim)%c_array(ii,j)
                     END DO
                   END DO
                END DO ! idim
              END IF
            END IF
          END IF ! ns_me
        END DO ! ip

        IF(ns_partner*ns_me /= 0) THEN
          DEALLOCATE(rmat_loc, STAT=istat)
          CPASSERT(istat==0)
          DO idim = 1,dim2
            DEALLOCATE(xyz_mix(idim)%c_array,STAT=istat)
            CPASSERT(istat==0)
          END DO
        END IF

      END DO  ! iperm

      ! calculate the max gradient
      DO idim = 1,dim2
        DO i = ns_bound(para_env%mepos,1), ns_bound(para_env%mepos,2)
          ii = i-ns_bound(para_env%mepos,1)+1
          zdiag_me(idim)%c_array(ii) =  cz_ij_loc(idim)%c_array(i,ii)
          zdiag_me(idim)%c_array(ii) =  cz_ij_loc(idim)%c_array(i,ii)
        END DO
        rcount(:) = SIZE(zdiag_me(idim)%c_array)
        rdispl(1)=0
        DO ip=2,para_env%num_pe
          rdispl(ip)=rdispl(ip-1)+rcount(ip-1)
        ENDDO
        ! collect all the diagonal elements in a replicated 1d array
        CALL mp_allgather(zdiag_me(idim)%c_array,zdiag_all(idim)%c_array,rcount,rdispl,para_env%group)
      END DO

      gmax = 0.0_dp
      DO j = ns_bound(para_env%mepos,1),  ns_bound(para_env%mepos,2)
        k= j-ns_bound(para_env%mepos,1)+1
        DO i = 1,j-1
          ! find the location of the diagonal element (i,i)
          DO ip = 0,para_env%num_pe-1
             IF(i>=ns_bound(ip,1) .AND. i<=ns_bound(ip,2)) THEN
               ip_has_i = ip
               EXIT
             END IF
          END DO
          ii = nblock_max*ip_has_i + i - ns_bound(ip_has_i,1) + 1
          ! mepos has the diagonal element (j,j), as well as the off diagonal (i,j)
          jj = nblock_max*para_env%mepos + j - ns_bound(para_env%mepos,1) +1
          grad = 0.0_dp
          DO idim = 1,dim2
            zi = zdiag_all(idim)%c_array(ii)
            zj = zdiag_all(idim)%c_array(jj)
            grad = grad + weights(idim) * REAL(4.0_dp*CONJG( cz_ij_loc(idim)%c_array(i,k))*(zj-zi),dp )
          END DO
          gmax = MAX(gmax,ABS(grad))
        END DO
      END DO

      CALL mp_max(gmax,para_env%group)
      tolerance = gmax

      func = 0.0_dp
      DO i=ns_bound(para_env%mepos,1),  ns_bound(para_env%mepos,2)
        k = i-ns_bound(para_env%mepos,1) + 1
        DO idim = 1,dim2
          zr=REAL(cz_ij_loc(idim)%c_array(i,k),dp)
          zc=AIMAG(cz_ij_loc(idim)%c_array(i,k))
          func = func + weights(idim) *(1.0_dp - (zr*zr+zc*zc))/twopi/twopi
        END DO
      END DO
      CALL mp_sum(func,para_env%group)
      t2 =  m_walltime()

      IF(output_unit>0  .AND. MODULO(sweeps,out_each)==0) THEN
         WRITE(output_unit,'(T20,I12,T35,F20.10,T60,E12.4,F8.3)') sweeps,func, tolerance,t2-t1
         CALL m_flush(output_unit)
      END IF
      IF(tolerance < eps_localization) EXIT
      CALL external_control(should_stop,"LOC",target_time=target_time, start_time=start_time)
      IF(should_stop) EXIT


    END DO ! sweeps

    ! buffer for message passing
    DEALLOCATE(rmat_recv_all,STAT=istat)
    CPASSERT(istat==0)

    DEALLOCATE(rmat_recv, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(rmat_send, STAT=istat)
    CPASSERT(istat==0)
    IF (ns_me>0) THEN
       DEALLOCATE(c_array_me,STAT=istat)
       CPASSERT(istat==0)
    ENDIF
    DO idim = 1,dim2
      DEALLOCATE(zdiag_me(idim)%c_array, STAT=istat)
      CPASSERT(istat==0)
      DEALLOCATE(zdiag_all(idim)%c_array, STAT=istat)
      CPASSERT(istat==0)
    END DO
    DEALLOCATE(zdiag_me, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(zdiag_all, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(xyz_mix,STAT=istat)
    CPASSERT(istat==0)
    DO idim = 1,dim2
      IF(ns_me/=0) THEN
        DEALLOCATE(xyz_mix_ns(idim)%c_array,STAT=istat)
        CPASSERT(istat==0)
      ENDIF
    END DO
    DEALLOCATE(xyz_mix_ns,STAT=istat)
    CPASSERT(istat==0)
    IF(ns_me/=0) THEN
      DEALLOCATE(gmat,STAT=istat)
      CPASSERT(istat==0)
    ENDIF
    DEALLOCATE(mii,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(mij,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(mjj,STAT=istat)
    CPASSERT(istat==0)

    ilow1=ns_bound(para_env%mepos,1)
    ns_me=ns_bound(para_env%mepos,2)-ns_bound(para_env%mepos,1)+1
    ALLOCATE(z_ij_loc_re(nstate,nblock_max),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(z_ij_loc_im(nstate,nblock_max),STAT=istat)
    CPASSERT(istat==0)
    DO idim = 1,dim2
      DO ip = 0,para_env%num_pe-1
        z_ij_loc_re = 0.0_dp
        z_ij_loc_im = 0.0_dp
        nblock = ns_bound(ip,2)-ns_bound(ip,1) +1
        IF(ip==para_env%mepos) THEN
          ns_me = nblock
          DO i=1,ns_me
            ii = ilow1 + i -1
            DO j = 1,nstate
              z_ij_loc_re(j,i)=REAL(cz_ij_loc(idim)%c_array(j,i),dp)
              z_ij_loc_im(j,i)=AIMAG(cz_ij_loc(idim)%c_array(j,i))
            END DO
          END DO
        END IF
        CALL mp_bcast(z_ij_loc_re,ip,para_env%group)
        CALL mp_bcast(z_ij_loc_im,ip,para_env%group)
        CALL cp_fm_set_submatrix(zij(1,idim)%matrix,z_ij_loc_re,1,ns_bound(ip,1),nstate,nblock)
        CALL cp_fm_set_submatrix(zij(2,idim)%matrix,z_ij_loc_im,1,ns_bound(ip,1),nstate,nblock)
      END DO ! ip
    END DO

    DO ip = 0,para_env%num_pe-1
      z_ij_loc_re = 0.0_dp
      nblock = ns_bound(ip,2)-ns_bound(ip,1) +1
      IF(ip==para_env%mepos) THEN
         ns_me = nblock
         DO i=1,ns_me
           ii = ilow1 + i -1
           DO j = 1,nstate
             z_ij_loc_re(j,i)=rotmat(j,i)
           END DO
         END DO
      END IF
      CALL mp_bcast(z_ij_loc_re,ip,para_env%group)
      CALL cp_fm_set_submatrix(rmat,z_ij_loc_re,1,ns_bound(ip,1),nstate,nblock)
    END DO

    DEALLOCATE(z_ij_loc_re,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(z_ij_loc_im,STAT=istat)
    CPASSERT(istat==0)
    DO idim = 1,dim2
       DEALLOCATE(cz_ij_loc(idim)%c_array,STAT=istat)
       CPASSERT(istat==0)
    END DO
    DEALLOCATE(cz_ij_loc,STAT=istat)
    CPASSERT(istat==0)

    CALL mp_sync(para_env%group)
    CALL rotate_orbitals ( rmat, vectors )
    CALL cp_fm_release ( rmat)

    DEALLOCATE(rotmat,STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(ns_bound, list_pair)

    CALL timestop(handle)

  END SUBROUTINE jacobi_rot_para

! *****************************************************************************
!> \brief ...
!> \param iperm ...
!> \param para_env ...
!> \param list_pair ...
! *****************************************************************************
  SUBROUTINE eberlein(iperm,para_env,list_pair)
    INTEGER, INTENT(IN)                      :: iperm
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, DIMENSION(:, :)                 :: list_pair

    CHARACTER(len=*), PARAMETER :: routineN = 'eberlein', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, jj, npair

      npair=(para_env%num_pe+1)/2
      IF(iperm==1) THEN
!..set up initial ordering
        DO I=0,para_env%num_pe-1
          II=((i+1)+1)/2
          JJ=MOD((i+1)+1,2)+1
          list_pair(JJ,II)=i
        ENDDO
        IF(MOD(para_env%num_pe,2)==1) list_pair(2,npair)=-1
      ELSEIF(MOD(iperm,2)==0) THEN
!..a type shift
        jj=list_pair(1,npair)
        DO I=npair,3,-1
          list_pair(1,I)=list_pair(1,I-1)
        ENDDO
        list_pair(1,2)=list_pair(2,1)
        list_pair(2,1)=jj
      ELSE
!..b type shift
        jj=list_pair(2,1)
        DO I=1,npair-1
          list_pair(2,I)=list_pair(2,I+1)
        ENDDO
        list_pair(2,npair)=jj
      ENDIF

  END SUBROUTINE eberlein

END MODULE qs_localization_methods
