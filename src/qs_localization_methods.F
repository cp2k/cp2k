!   CP2K: A general program to perform molecular dynamics simulations         !
!!****** cp2k/qs_localization_methods [1.0] *
!!
!!   NAME
!!     qs_localization_methods
!!
!!   FUNCTION
!!     Localization methods such as 2x2 Jacobi rotations 
!!                                  Steepest Decents
!!                                  Conjugate Gradient
!! 
!!
!!   AUTHORS
!!     CJM (04.2003)
!!
!!   MODIFICATION HISTORY
!!     Initial parallellization of jacobi (JVDV 07.2003)
!!     direct minimization using exponential parametrization (JVDV 09.2003)
!!     crazy rotations go fast (JVDV 10.2003)
!!
!!   SOURCE
!******************************************************************************
MODULE qs_localization_methods
! *****************************************************************************
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
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace,&
                                             cp_fm_transpose
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_init_random,&
                                             cp_fm_maxval,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type,&
                                             cp_fm_write
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_max
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PUBLIC :: jacobi_rotations, initialize_weights, crazy_rotations, direct_mini
  PRIVATE

! *****************************************************************************
CONTAINS
! *****************************************************************************
  SUBROUTINE initialize_weights ( cell, weights )
    IMPLICIT NONE
    TYPE ( cell_type ), INTENT ( IN ) :: cell
    REAL (KIND = dp), DIMENSION ( : ) :: weights
! locals
    INTEGER :: idim
    DO idim = 1, 3
      weights ( idim ) = cell % hmat ( idim, idim ) * cell % hmat ( idim, idim )
    ENDDO
  END SUBROUTINE initialize_weights 

! *****************************************************************************
  SUBROUTINE jacobi_rotations ( weights, zij, vectors, sweeps )
    IMPLICIT NONE
    REAL (KIND = dp), INTENT ( IN ) :: weights ( : )
    TYPE ( cp_fm_p_type ), INTENT ( INOUT ) :: ZIJ ( :, : )
    TYPE ( cp_fm_type ), POINTER :: vectors
    INTEGER :: sweeps
! Locals
    TYPE ( cp_fm_type ), POINTER :: rmat
    TYPE ( cp_cfm_type ), POINTER :: c_rmat
    TYPE ( cp_cfm_p_type ), POINTER, DIMENSION(:) :: c_zij

    REAL ( KIND = dp ) :: tolerance, st, ct 
    COMPLEX (KIND = dp ) :: mii ( 3 ), mjj ( 3 ), mij ( 3 )
    REAL ( KIND = dp ) :: theta
    INTEGER :: istate, jstate, idim, nstate, isos, handle
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_fm_type ), POINTER :: matrix_r, matrix_i
    COMPLEX( KIND = dp ) :: dum

    CALL timeset("wannier_jacobi","I","",handle)

    NULLIFY(rmat,c_rmat,c_zij)
    ALLOCATE(c_zij(SIZE(zij,2)))

    CALL cp_fm_create ( rmat, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_set_all ( rmat, 0._dp, 1._dp, error )

    CALL cp_cfm_create ( c_rmat, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_set_all ( c_rmat, (0._dp,0._dp) , (1._dp,0._dp) )
    DO idim=1,SIZE(c_zij)
       NULLIFY(c_zij(idim)%matrix)
       CALL cp_cfm_create ( c_zij(idim)%matrix, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
       c_zij(idim)%matrix% local_data = CMPLX (zij(1,idim) % matrix % local_data, &
                    zij(2,idim) % matrix % local_data, dp )
    ENDDO

    CALL cp_fm_get_info ( rmat , nrow_global = nstate, error = error )
    tolerance = 1.0e10_dp
    sweeps = 0
! do jacobi sweeps until converged
    DO WHILE ( tolerance >= 1.0e-4_dp ) ! should be already machine precision minimization
      sweeps = sweeps + 1
      DO istate = 1, nstate
        DO jstate = istate + 1, nstate
          DO idim = 1, SIZE ( c_zij, 1 )
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, istate, istate, mii(idim) )
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, istate, jstate, mij(idim) )
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, jstate, jstate, mjj(idim) )
          END DO 
          CALL get_angle ( mii, mjj, mij, weights, theta )
          st = DSIN ( theta )
          ct = DCOS ( theta )
          CALL rotate_zij ( istate, jstate, st, ct, c_zij )
          CALL rotate_rmat ( istate, jstate, st, ct, c_rmat )
        END DO
      END DO
      CALL check_tolerance ( c_zij, weights, tolerance )
    END DO

    DO idim=1,SIZE(c_zij)
       zij(1,idim) % matrix % local_data =  REAL ( c_zij(idim)%matrix% local_data , dp )
       zij(2,idim) % matrix % local_data = AIMAG ( c_zij(idim)%matrix% local_data )
       CALL cp_cfm_release( c_zij(idim)%matrix )
    ENDDO
    DEALLOCATE(c_zij)
    rmat % local_data = REAL ( c_rmat%local_data, dp )
    CALL cp_cfm_release( c_rmat )
    CALL rotate_orbitals ( rmat, vectors ) 
    CALL cp_fm_release ( rmat )

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE jacobi_rotations
! *****************************************************************************
  SUBROUTINE rotate_zij ( istate, jstate, st, ct, zij )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: istate, jstate
    TYPE ( cp_cfm_p_type ) :: zij ( : )
    REAL ( KIND = dp ), INTENT ( IN ) :: st, ct
! Locals
    TYPE ( cp_error_type ) :: error
    INTEGER :: idim, nstate, stride
    COMPLEX ( KIND = dp ) :: st_cmplx
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

    st_cmplx = CMPLX ( st, 0.0_dp, dp )
    CALL cp_cfm_get_info ( zij ( 1 ) % matrix, nrow_global = nstate, error = error )
    DO  idim = 1, 3
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
  SUBROUTINE rotate_rmat ( istate, jstate, st, ct, rmat )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: istate, jstate
    TYPE ( cp_cfm_type ), POINTER :: rmat
    REAL ( KIND = dp ), INTENT ( IN ) :: ct, st
! Locals
    TYPE ( cp_error_type ) :: error
    INTEGER :: nstate
    TYPE ( cp_cfm_type ), POINTER :: rmat_cmplx
    COMPLEX ( KIND = dp ) :: st_cmplx
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

    st_cmplx = CMPLX ( st, 0.0_dp, dp )
    CALL cp_cfm_get_info (  rmat, nrow_global = nstate, error = error )
#if defined(__SCALAPACK)
    desc ( : ) = rmat % matrix_struct % descriptor(:)
    CALL pzrot( nstate , rmat % local_data (1,1) , 1 , istate,  desc, 1, &
                  rmat % local_data (1,1) , 1, jstate, desc, 1, ct, st_cmplx )
#else
    CALL zrot( nstate, rmat % local_data(1,istate), 1, rmat % local_data(1,jstate), 1, ct, st_cmplx )
#endif
  END SUBROUTINE rotate_rmat
! *****************************************************************************
  SUBROUTINE get_angle ( mii, mjj, mij,  weights, theta )
    IMPLICIT NONE
    COMPLEX ( KIND = dp ), INTENT ( IN ) :: mii ( : ), mjj ( : ), mij ( : )
    REAL ( KIND = dp ), INTENT ( IN ) :: weights ( : )
    REAL ( KIND = dp ), INTENT ( OUT ) :: theta
! Locals
    COMPLEX ( KIND = dp ) :: z11, z12, z22
    REAL ( KIND = dp ) :: a12, b12, ratio, d2
    INTEGER :: idim
    a12 = 0.0_dp
    b12 = 0.0_dp
    DO idim = 1, 3
      z11 = mii(idim) 
      z22 = mjj(idim) 
      z12 = mij(idim)
      a12 = a12 + weights ( idim ) * REAL ( CONJG ( z12 ) * ( z11 - z22 ), dp )
      b12 = b12 + weights ( idim ) * (  z12 * CONJG ( z12 ) -  &
            0.25_dp * ( z11 - z22 ) * ( CONJG ( z11 ) - CONJG ( z22 ) ) ) 
    END DO
    IF ( ABS ( b12 ) > 1.e-10_dp ) THEN
      ratio = -a12/b12
      theta = 0.25_dp * DATAN ( ratio )
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
  SUBROUTINE check_tolerance ( zij, weights, tolerance )
    IMPLICIT NONE
    TYPE ( cp_cfm_p_type ) :: zij ( : ) 
    REAL (KIND = dp), INTENT ( OUT ) :: tolerance
    REAL (KIND = dp), INTENT ( IN ) :: weights ( : )
! locals
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_fm_type ), POINTER :: force
! compute gradient at t=0
    NULLIFY ( force )
    CALL cp_fm_create ( force, zij( 1 ) % matrix % matrix_struct, error = error  )
    CALL cp_fm_set_all ( force, 0._dp, error=error )
    CALL grad_at_0 ( zij, weights, force )
    CALL cp_fm_maxval ( force, tolerance, error = error )
    CALL cp_fm_release ( force )
  END SUBROUTINE check_tolerance
! *****************************************************************************
  SUBROUTINE rotate_orbitals ( rmat, vectors ) 
    IMPLICIT NONE
    TYPE ( cp_fm_type ), POINTER :: rmat
    TYPE ( cp_fm_type ), POINTER :: vectors
!locals
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_fm_type ), POINTER :: wf
    INTEGER :: n, k
    NULLIFY ( wf )
    CALL cp_fm_create ( wf, vectors % matrix_struct, error = error  )
    CALL cp_fm_get_info ( vectors, nrow_global = n, ncol_global=k )
    CALL cp_fm_gemm("N", "N", n, k, k, 1.0_dp, vectors, rmat, 0.0_dp, wf )
    CALL cp_fm_to_fm ( wf, vectors )
    CALL cp_fm_release ( wf )
  END SUBROUTINE rotate_orbitals 
! *****************************************************************************
  SUBROUTINE gradsq_at_0 ( diag, weights, matrix )
    IMPLICIT NONE
    TYPE ( cp_fm_type ), POINTER :: matrix
    REAL ( KIND = dp ), INTENT ( IN ) :: weights ( : )
! Locals
    TYPE ( cp_error_type ) :: error
    COMPLEX ( KIND = dp ) :: zii, zjj
    COMPLEX ( KIND = dp ), DIMENSION(:,:), POINTER  :: diag
    REAL ( KIND = dp ) ::  gradsq_ij
    INTEGER :: nrow_local, ncol_local, istate, jstate, idim, nrow_global
    INTEGER, POINTER, DIMENSION(:) :: row_indices,col_indices

    CALL cp_fm_get_info ( matrix, nrow_local = nrow_local,  &
                   ncol_local = ncol_local, nrow_global = nrow_global, &
                   row_indices = row_indices, col_indices = col_indices, error = error )

    DO istate = 1, nrow_local
      DO jstate = 1, ncol_local
! get real and imaginary parts
        gradsq_ij = 0.0_dp
        DO idim = 1, 3
          zii = diag(row_indices(istate),idim)
          zjj = diag(col_indices(jstate),idim)
          gradsq_ij = gradsq_ij + weights ( idim ) * &
                     4.0_dp * ( CONJG ( zii ) * zii + CONJG ( zjj ) * zjj )
        END DO
        matrix % local_data ( istate, jstate ) = gradsq_ij
      END DO
    END DO
  END SUBROUTINE gradsq_at_0
! *****************************************************************************
  SUBROUTINE grad_at_0 ( matrix_p, weights, matrix )
    IMPLICIT NONE
    TYPE ( cp_cfm_p_type ) :: matrix_p ( : )
    TYPE ( cp_fm_type ), POINTER :: matrix
    REAL ( KIND = dp ), INTENT ( IN ) :: weights ( : )
! Locals
    TYPE ( cp_error_type ) :: error
    COMPLEX ( KIND = dp ) :: zii, zjj, zij 
    COMPLEX ( KIND = dp ), DIMENSION(:,:), POINTER  :: diag
    REAL ( KIND = dp ) :: grad_ij
    INTEGER :: nrow_local, ncol_local, istate, jstate, idim, nrow_global
    INTEGER, POINTER, DIMENSION(:) :: row_indices,col_indices

    NULLIFY(diag)
    CALL cp_fm_get_info ( matrix, nrow_local = nrow_local,  &
                   ncol_local = ncol_local, nrow_global = nrow_global, &
                   row_indices = row_indices, col_indices = col_indices, error = error )
    ALLOCATE(diag(nrow_global,3))

    DO idim = 1, 3
      DO istate = 1, nrow_global
         CALL cp_cfm_get_element(matrix_p(idim)%matrix,istate,istate,diag(istate,idim))
      ENDDO
    ENDDO

    DO istate = 1, nrow_local
      DO jstate = 1, ncol_local
! get real and imaginary parts
        grad_ij = 0.0_dp
        DO idim = 1, 3
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

! *****************************************************************************
! return energy and maximum gradient in the current point
! *****************************************************************************
  SUBROUTINE check_tolerance_new( weights, zij, tolerance, value)
    IMPLICIT NONE
    REAL (KIND = dp), INTENT ( IN ) :: weights ( : )
    TYPE ( cp_fm_p_type ), INTENT ( INOUT ) :: ZIJ ( :, : )
    REAL(KIND = dp) :: tolerance,value

    COMPLEX(KIND = dp), DIMENSION(:,:), POINTER :: diag
    COMPLEX(KIND = dp) :: kij,kii,kjj
    REAL(KIND = dp) :: ra,rb,grad_ij
    INTEGER :: nrow_local, ncol_local, istate, jstate, idim, nrow_global
    INTEGER, POINTER, DIMENSION(:) :: row_indices,col_indices

    NULLIFY(diag)
    CALL cp_fm_get_info ( zij(1,1)%matrix, nrow_local = nrow_local,  &
                   ncol_local = ncol_local, nrow_global = nrow_global, &
                   row_indices = row_indices, col_indices = col_indices)
    ALLOCATE(diag(nrow_global,3))
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
         DO idim = 1, 3
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

  END SUBROUTINE
! *****************************************************************************
! yet another crazy try, computes the angles needed to rotate the orbitals first
! and rotates them all at the same time (hoping for the best of course)
! *****************************************************************************
  SUBROUTINE crazy_rotations( weights, zij, vectors, iterations )
    IMPLICIT NONE
    REAL (KIND = dp), INTENT ( IN ) :: weights ( : )
    TYPE ( cp_fm_p_type ), INTENT ( INOUT ) :: ZIJ ( :, : )
    TYPE ( cp_fm_type ), POINTER :: vectors
    INTEGER :: iterations

    COMPLEX(KIND = dp), DIMENSION(:,:), POINTER :: diag_z
    INTEGER :: nrow_global,ncol_global,nrow_local,ncol_local
    INTEGER, DIMENSION(:), POINTER :: row_indices, col_indices
    INTEGER :: i,idim,irow,icol,handle
    COMPLEX(KIND = dp), DIMENSION(:), POINTER :: evals_exp
    REAL(KIND = dp), DIMENSION(:), POINTER :: evals
    REAL(KIND = dp)                        :: rance,theta,maxeval,value,tolerance
    REAL(KIND = dp)                        :: ra,rb

    TYPE( cp_cfm_type), POINTER :: cmat_A,cmat_R,cmat_t1
    TYPE( cp_fm_type), POINTER :: mat_U,mat_t,mat_R
    COMPLEX (KIND = dp ) :: mii ( 3 ), mjj ( 3 ), mij ( 3 )
    COMPLEX(KIND = dp), PARAMETER :: czero=(0.0_dp,0.0_dp), cone=(1.0_dp,0.0_dp)

    CALL timeset("crazy_rotations","I","",handle)
    NULLIFY(row_indices,col_indices)
    NULLIFY(mat_U,mat_t,mat_R)
    NULLIFY(cmat_A,cmat_R,cmat_t1)
    CALL cp_fm_get_info(zij(1,1)%matrix,nrow_global=nrow_global, ncol_global=ncol_global, &
                                        row_indices=row_indices, col_indices=col_indices, &
                                        nrow_local=nrow_local, ncol_local=ncol_local)
  
    ALLOCATE(diag_z(nrow_global,3))
    ALLOCATE(evals(nrow_global))
    ALLOCATE(evals_exp(nrow_global))

    CALL cp_cfm_create ( cmat_A, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_R, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_cfm_create ( cmat_t1, zij ( 1, 1 ) % matrix % matrix_struct)

    CALL cp_fm_create ( mat_U, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create ( mat_t, zij ( 1, 1 ) % matrix % matrix_struct)
    CALL cp_fm_create ( mat_R, zij ( 1, 1 ) % matrix % matrix_struct)

    CALL cp_fm_set_all( mat_R,0.0_dp,1.0_dp )
    CALL cp_fm_set_all( mat_t, 0.0_dp)
    DO idim=1,SIZE(zij,2)
       CALL cp_fm_scale_and_add(1.0_dp,mat_t,1.0_dp,zij(1,idim) % matrix)
       CALL cp_fm_scale_and_add(1.0_dp,mat_t,1.0_dp,zij(2,idim) % matrix)
    ENDDO
    CALL cp_fm_syevd(mat_t,mat_U,evals)
    DO idim=1,SIZE(zij,2)
       ! rotate z's
       CALL cp_fm_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(1,idim)%matrix,mat_U,0.0_dp,mat_t)
       CALL cp_fm_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(1,idim)%matrix)
       CALL cp_fm_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(2,idim)%matrix,mat_U,0.0_dp,mat_t)
       CALL cp_fm_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(2,idim)%matrix)
    ENDDO
    ! collect rotation matrix
    CALL cp_fm_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_R,mat_U,0.0_dp,mat_t)
    CALL cp_fm_to_fm(mat_t,mat_R)

    iterations=0
    tolerance=1.0_dp
    DO 
      iterations=iterations+1
      DO idim=1,SIZE(zij,2)
         DO i=1,nrow_global 
            CALL cp_fm_get_element(zij(1,idim) % matrix, i, i, ra)
            CALL cp_fm_get_element(zij(2,idim) % matrix, i, i, rb)
            diag_z(i,idim)=CMPLX(ra,rb,dp)
         ENDDO
      ENDDO
      DO irow=1,nrow_local
         DO icol=1,ncol_local
            DO idim=1,SIZE(zij,2)
               ra=zij(1,idim) % matrix % local_data(irow,icol)
               rb=zij(2,idim) % matrix % local_data(irow,icol)
               mij(idim)=CMPLX(ra,rb,dp)
               mii(idim)=diag_z(row_indices(irow),idim)
               mjj(idim)=diag_z(col_indices(icol),idim)
            ENDDO
            IF (row_indices(irow).NE.col_indices(icol)) THEN
              CALL get_angle ( mii, mjj, mij, weights, theta )
              IF (theta.gt.0.2_dp)  theta=0.2_dp
              IF (theta.lt.-0.2_dp) theta=-0.2_dp
              cmat_A % local_data(irow,icol) = -CMPLX(0.0_dp, theta, dp )
            ELSE
              cmat_A % local_data(irow,icol) = czero
            ENDIF
         ENDDO
      ENDDO

      ! construct rotation matrix U based on A using diagonalization
      ! alternatively, exp based on repeated squaring could be faster
      CALL cp_cfm_heevd(cmat_A,cmat_R,evals)
      maxeval=MAXVAL(ABS(evals))
      evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * evals(:) )
      CALL cp_cfm_to_cfm(cmat_R,cmat_t1)
      CALL cp_cfm_column_scale(cmat_t1,evals_exp)
      CALL cp_cfm_gemm('N','C',nrow_global,nrow_global,nrow_global,cone,cmat_t1,cmat_R,czero,cmat_A)
      mat_U%local_data=REAL(cmat_A%local_data) ! U is a real matrix

      DO idim=1,SIZE(zij,2)
         ! rotate z's
         CALL cp_fm_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(1,idim)%matrix,mat_U,0.0_dp,mat_t)
         CALL cp_fm_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(1,idim)%matrix)
         CALL cp_fm_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,zij(2,idim)%matrix,mat_U,0.0_dp,mat_t)
         CALL cp_fm_gemm('T','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_U,mat_t,0.0_dp,zij(2,idim)%matrix)
      ENDDO
      ! collect rotation matrix
      CALL cp_fm_gemm('N','N',nrow_global,nrow_global,nrow_global,1.0_dp,mat_R,mat_U,0.0_dp,mat_t)
      CALL cp_fm_to_fm(mat_t,mat_R)

      CALL check_tolerance_new ( weights, zij, tolerance, value )
      IF (cmat_A%matrix_struct%para_env%mepos .EQ. 0 ) THEN 
          WRITE(6,'(I4,1X,G20.15,E12.4,E12.4)') iterations,value,tolerance,maxeval
      ENDIF
      IF (tolerance .LT. 1.0E-4_dp .OR. iterations.gt.200) EXIT
    ENDDO

    CALL cp_cfm_release(cmat_A)
    CALL cp_cfm_release(cmat_R)
    CALL cp_cfm_release(cmat_T1)

    CALL cp_fm_release(mat_U)
    CALL cp_fm_release(mat_T)

    CALL rotate_orbitals(mat_R,vectors)

    CALL cp_fm_release(mat_R)
    DEALLOCATE(evals_exp,evals,diag_z)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE crazy_rotations
! *****************************************************************************
!
! use the exponential parametrization as described in to perform a direct mini
!
! Gerd Berghold et al. PRB 61 (15), pag. 10040 (2000)
!
! none of the input is modified for the time being, just finds the rotations 
! that minimizes, and throws it away afterwards :-)
! apart from being expensive and not cleaned, this works fine
! useful to try different spread functionals
!
! *****************************************************************************
  SUBROUTINE direct_mini( weights, zij, vectors, iterations )
    IMPLICIT NONE
    REAL (KIND = dp), INTENT ( IN ) :: weights ( : )
    TYPE ( cp_fm_p_type ), INTENT ( INOUT ) :: ZIJ ( :, : )
    TYPE ( cp_fm_type ), POINTER :: vectors
    INTEGER :: iterations,handle,line_searches

    TYPE ( cp_fm_type ), POINTER :: matrix_A,matrix_G,matrix_T,matrix_G_search,matrix_G_old,matrix_H,matrix_R
    TYPE ( cp_cfm_p_type ), POINTER, DIMENSION(:) :: c_zij
    TYPE ( cp_cfm_type ), POINTER :: cmat_A,cmat_R,cmat_t1,cmat_t2,cmat_U, cmat_B, cmat_M
    COMPLEX(KIND = dp), DIMENSION(:), POINTER :: evals_exp
    COMPLEX(KIND = dp), DIMENSION(:,:), POINTER :: diag_z
    REAL(KIND = dp), DIMENSION(:), POINTER :: evals,fval,fvald
    INTEGER :: n,idim
    COMPLEX(KIND = dp), PARAMETER :: czero=(0.0_dp,0.0_dp), cone=(1.0_dp,0.0_dp)
    TYPE ( cp_error_type ) :: error
    REAL(KIND = dp) :: omega,tol
    INTEGER :: irow, icol, nrow_local, ncol_local, i
    INTEGER, DIMENSION(:), POINTER :: row_indices, col_indices
    COMPLEX(KIND = dp) :: ll,lk,tmp

    REAL(KIND = dp) :: normg,normg_cross,normg_old,beta_pr,ds_min,ds
    REAL(KIND = dp) :: a,b,c,x0,x1,val,npos,xa,xb,xc,fa,fb,fc,nom,denom
    REAL(KIND = dp), DIMENSION(150) :: energy,pos,grad
    INTEGER :: line_search_count,lsl,lsr,lsm
    LOGICAL  :: new_direction
    REAL(KIND = dp), PARAMETER                       :: gold_sec=0.3819_dp

    NULLIFY(evals,evals_exp,diag_z,fval,fvald,c_zij) 
    NULLIFY(matrix_A,matrix_G,matrix_T,matrix_G_search,matrix_G_old)
    NULLIFY(cmat_A,cmat_U,cmat_R,cmat_t1,cmat_t2,cmat_B,cmat_M)

    CALL timeset("direct_mini","I","",handle)
    n = zij(1,1)%matrix%matrix_struct%nrow_global
    ALLOCATE(evals(n),evals_exp(n),diag_z(n,3),fval(n),fvald(n))
    ALLOCATE(c_zij(SIZE(zij,2)))


    ! create the three complex matrices Z
    DO idim=1,SIZE(c_zij)
       NULLIFY(c_zij(idim)%matrix)
       CALL cp_cfm_create ( c_zij(idim)%matrix, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
       c_zij(idim)%matrix% local_data = CMPLX (zij(1,idim) % matrix % local_data, &
                    zij(2,idim) % matrix % local_data, dp )
    ENDDO

    CALL cp_fm_create  ( matrix_A, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_create  ( matrix_G, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_create  ( matrix_T, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_create  ( matrix_H, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_create  ( matrix_G_search, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_create  ( matrix_G_old, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_create  ( matrix_R, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_fm_set_all(matrix_R,0.0_dp,1.0_dp)
    
    CALL cp_fm_set_all(matrix_A,0.0_dp)
!    CALL cp_fm_init_random ( matrix_A )

    CALL cp_cfm_create ( cmat_A, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_create ( cmat_U, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_create ( cmat_R, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_create ( cmat_t1, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_create ( cmat_t2, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_create ( cmat_B, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_create ( cmat_M, zij ( 1, 1 ) % matrix % matrix_struct, error=error )

    CALL cp_cfm_get_info ( cmat_B, nrow_local=nrow_local, ncol_local=ncol_local, &
                                   row_indices=row_indices, col_indices=col_indices )

    CALL cp_fm_set_all(matrix_G_old,0.0_dp)
    CALL cp_fm_set_all(matrix_G_search,0.0_dp)
    normg_old=1.0E99_dp
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
         cmat_U%local_data=REAL(cmat_U%local_data) ! enforce numerics, U is a real matrix
         
         IF ( new_direction .AND. MOD(line_searches,20).EQ.5 ) THEN ! reset with A .eq. 0
            DO idim=1,SIZE(c_zij)
               CALL cp_cfm_gemm('N','N',n,n,n,cone,c_zij(idim)%matrix,cmat_U,czero,cmat_t1) 
               CALL cp_cfm_gemm('C','N',n,n,n,cone,cmat_U,cmat_t1,czero,c_zij(idim)%matrix) 
            ENDDO
            ! collect rotation matrix
            matrix_H%local_data=REAL(cmat_U%local_data)
            CALL cp_fm_gemm('N','N',n,n,n,1.0_dp,matrix_R,matrix_H,0.0_dp,matrix_T)
            CALL cp_fm_to_fm(matrix_T,matrix_R)

            CALL cp_cfm_set_all(cmat_U,czero,cone)
            CALL cp_cfm_set_all(cmat_R,czero,cone)
            CALL cp_cfm_set_all(cmat_A,czero)
            CALL cp_fm_set_all(matrix_A,0.0_dp)
            evals(:)=0.0_dp
            evals_exp(:)=EXP( (0.0_dp,-1.0_dp) * evals(:) )
            CALL cp_fm_set_all(matrix_G_old,0.0_dp)
            CALL cp_fm_set_all(matrix_G_search,0.0_dp)
            normg_old=1.0E99_dp
         ENDIF

         ! compute Omega and M
         CALL cp_cfm_set_all(cmat_M,czero)
         omega = 0.0_dp
         DO idim=1,SIZE(c_zij)
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
                                              +4.0_dp*fvald(col_indices(icol))*REAL(tmp)
               ENDDO
            ENDDO
         ENDDO

         ! compute Hessian diagonal approximation
         IF (.TRUE.) THEN
            CALL gradsq_at_0(diag_z,weights,matrix_H)
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
                cmat_B%local_data(irow,icol)=0.0
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
         matrix_G%local_data=REAL(cmat_t1%local_data)
         CALL cp_fm_transpose(matrix_G,matrix_T)
         CALL cp_fm_scale_and_add(  -1.0_dp,matrix_G, 1.0_dp,matrix_T)
         CALL cp_fm_maxval(matrix_G,tol)


         ! from here on, minimizing technology
         IF ( new_direction ) THEN
            ! energy converged up to machine precision ?
            line_searches=line_searches+1
            ! DO i=1,line_search_count
            !   write(15,*) pos(i),energy(i)
            ! ENDDO
            ! write(15,*) ""
            CALL m_flush(15)
            !write(16,*) evals(:)
            !write(17,*) matrix_A%local_data(:,:)
            !write(18,*) matrix_G%local_data(:,:)
            IF (matrix_A%matrix_struct%para_env%mepos .EQ. 0 ) WRITE(6,*) line_searches,Iterations,Omega,tol,ds_min
            IF (tol.lt.1.0E-3)  EXIT

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
                   WRITE(6,*) "!" 
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
    matrix_H%local_data=REAL(cmat_U%local_data)
    CALL cp_fm_gemm('N','N',n,n,n,1.0_dp,matrix_R,matrix_H,0.0_dp,matrix_T)
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
       CALL cp_cfm_release( c_zij(idim)%matrix )
    ENDDO
    DEALLOCATE(c_zij)
    DEALLOCATE(diag_z)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE
END MODULE qs_localization_methods
