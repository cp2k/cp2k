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
!!
!!   SOURCE
!******************************************************************************
MODULE qs_localization_methods
! *****************************************************************************
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_complex_full_matrix_type, &
                                             cp_complex_full_matrix_p_type, &
                                             cp_cfm_get_info, &
                                             cp_cfm_set_all, &
                                             cp_cfm_get_element
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_maxval,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE kinds,                           ONLY: wp => dp
  USE mathconstants,                   ONLY: pi
  USE qs_blacs,                        ONLY: write_blacs_matrix
  USE qs_environment_types,            ONLY: get_qs_env, qs_environment_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PUBLIC :: jacobi_rotations, initialize_weights
  PRIVATE

! *****************************************************************************
CONTAINS
! *****************************************************************************
  SUBROUTINE initialize_weights ( cell, weights )
    IMPLICIT NONE
    TYPE ( cell_type ), INTENT ( IN ) :: cell
    REAL ( wp ), DIMENSION ( : ) :: weights
! locals
    INTEGER :: idim
    DO idim = 1, 3
      weights ( idim ) = cell % hmat ( idim, idim ) * cell % hmat ( idim, idim )
    ENDDO
  END SUBROUTINE initialize_weights 
! *****************************************************************************
  SUBROUTINE jacobi_rotations ( weights, zij, vectors, sweeps )
    IMPLICIT NONE
    REAL ( wp ), INTENT ( IN ) :: weights ( : )
    TYPE ( cp_full_matrix_p_type ), INTENT ( INOUT ) :: ZIJ ( :, : )
    TYPE ( cp_full_matrix_type ), POINTER :: vectors
    INTEGER :: sweeps
! Locals
    TYPE ( cp_full_matrix_type ), POINTER :: rmat
    TYPE ( cp_complex_full_matrix_type ), POINTER :: c_rmat
    TYPE ( cp_complex_full_matrix_p_type ), POINTER, DIMENSION(:) :: c_zij

    REAL ( KIND = wp ) :: tolerance, st, ct, mii ( 2, 3 ), mjj ( 2, 3 ), mij ( 2, 3 )
    REAL ( KIND = wp ) :: theta
    INTEGER :: istate, jstate, idim, nstate, isos,sweeps
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_full_matrix_type ), POINTER :: matrix_r, matrix_i
    COMPLEX( KIND = wp ) :: dum

    NULLIFY(rmat,c_rmat,c_zij)
    ALLOCATE(c_zij(SIZE(zij,2)))

!Create the rotation matrix
    CALL cp_fm_create ( rmat, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
! Initialize rotation to the unit matrix
    CALL cp_fm_set_all ( rmat, 0._wp, 1._wp, error )

!Create the complex rotation matrix and z matrices, and avoid using z
    CALL cp_cfm_create ( c_rmat, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    CALL cp_cfm_set_all ( c_rmat, (0._wp,0._wp) , (1._wp,0._wp) )
    DO idim=1,SIZE(c_zij)
       NULLIFY(c_zij(idim)%matrix)
       CALL cp_cfm_create ( c_zij(idim)%matrix, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
       c_zij(idim)%matrix% local_data = CMPLX (zij(1,idim) % matrix % local_data, &
                    zij(2,idim) % matrix % local_data, wp )
    ENDDO

! get nstate
    CALL cp_fm_get_info ( rmat , nrow_global = nstate, error = error )
    tolerance = 1e10_wp
    sweeps = 0
    DO WHILE ( tolerance >= 1.0e-6_wp )
      sweeps = sweeps + 1
! ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE 
      DO istate = 1, nstate
! get real and imaginary parts
        DO idim = 1, SIZE ( c_zij, 1 )
          CALL cp_cfm_get_element ( c_zij(idim) % matrix, istate, istate, dum )
          mii (1,idim) = REAL(dum,wp)
          mii (2,idim) = AIMAG(dum)
        END DO
        DO jstate = istate + 1, nstate
! get real and imaginary parts
          DO idim = 1, SIZE ( c_zij, 1 )
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, istate, jstate, dum )
            mij (1,idim) = REAL(dum,wp)
            mij (2,idim) = AIMAG(dum)
            CALL cp_cfm_get_element ( c_zij(idim) % matrix, jstate, jstate, dum )
            mjj (1,idim) = REAL(dum,wp)
            mjj (2,idim) = AIMAG(dum)
          END DO 
! get the optimal angle
          CALL get_angle ( mii, mjj, mij, zij, weights, theta )
          st = DSIN ( theta )
          ct = DCOS ( theta )
! rotate the ZIJ
          CALL rotate_zij ( istate, jstate, st, ct, c_zij )
! update the Unitary matrix, rmat
          CALL rotate_rmat ( istate, jstate, st, ct, c_rmat )
        END DO
      END DO
! Check tolerance
      CALL check_tolerance ( c_zij, weights, tolerance )
    END DO
! Rotate the orbitals XXXXXXXXXXXXXXXXXX
    DO idim=1,SIZE(c_zij)
       zij(1,idim) % matrix % local_data =  REAL ( c_zij(idim)%matrix% local_data , wp )
       zij(2,idim) % matrix % local_data = AIMAG ( c_zij(idim)%matrix% local_data )
       CALL cp_cfm_release( c_zij(idim)%matrix )
    ENDDO
    DEALLOCATE(c_zij)
    rmat % local_data = REAL ( c_rmat%local_data, wp )
    CALL cp_cfm_release( c_rmat )
    CALL rotate_orbitals ( rmat, vectors ) 
    CALL cp_fm_release ( rmat )
! DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE 
  END SUBROUTINE jacobi_rotations
! *****************************************************************************
  SUBROUTINE rotate_zij ( istate, jstate, st, ct, zij )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: istate, jstate
    TYPE ( cp_complex_full_matrix_p_type ) :: zij ( : )
    REAL ( KIND = wp ), INTENT ( IN ) :: st, ct
! Locals
    TYPE ( cp_error_type ) :: error
    INTEGER :: idim, nstate, stride
    COMPLEX ( KIND = wp ) :: st_cmplx
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

    st_cmplx = CMPLX ( st, 0._wp, wp )
    ! XXXXXXXXXX is nstate equal M_X ??
    CALL cp_cfm_get_info ( zij ( 1 ) % matrix, nrow_global = nstate, error = error )
    DO  idim = 1, 3
#if defined(__SCALAPACK)
      desc(:) = zij(idim) % matrix %matrix_struct%descriptor(:)
      CALL pzrot( nstate , zij(idim) % matrix % local_data (1,1) , istate, 1, desc, 1, &
                  zij(idim) % matrix % local_data (1,1) , jstate, 1, desc, 1, ct, st_cmplx )
#else
      CALL zrot( nstate, zij(idim) % matrix % local_data ( 1, istate ), 1, &
                 zij(idim) % matrix % local_data ( 1, jstate ), 1, ct, st_cmplx )
#endif
#if defined(__SCALAPACK)
      desc(:) = zij(idim) % matrix %matrix_struct%descriptor(:)
      CALL pzrot( nstate , zij(idim) % matrix % local_data (1,1) , 1, istate, desc, nstate, &
                  zij(idim) % matrix % local_data (1,1) , 1, jstate, desc, nstate, ct, st_cmplx )
#else
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
    TYPE ( cp_complex_full_matrix_type ), POINTER :: rmat
    REAL ( KIND = wp ), INTENT ( IN ) :: ct, st
! Locals
    TYPE ( cp_error_type ) :: error
    INTEGER :: nstate
    TYPE ( cp_complex_full_matrix_type ), POINTER :: rmat_cmplx
    COMPLEX ( KIND = wp ) :: st_cmplx
#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

    st_cmplx = CMPLX ( st, 0._wp, wp )
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
  SUBROUTINE get_angle ( mii, mjj, mij, zij, weights, theta )
    IMPLICIT NONE
    TYPE ( cp_full_matrix_p_type ) :: zij ( :, : )
    REAL ( KIND = wp ), INTENT ( IN ) :: mii ( :, : ), mjj ( :, : ), mij ( :, : )
    REAL ( KIND = wp ), INTENT ( IN ) :: weights ( : )
    REAL ( KIND = wp ), INTENT ( OUT ) :: theta
! Locals
    COMPLEX ( KIND = wp ) :: z11, z12, z22
    REAL ( KIND = wp ) :: a12, b12, ratio, d2
    INTEGER :: idim
    a12 = 0._wp
    b12 = 0._wp
    DO idim = 1, 3
      z11 = CMPLX ( mii ( 1, idim ), mii ( 2, idim ), wp )
      z22 = CMPLX ( mjj ( 1, idim ), mjj ( 2, idim ), wp )
      z12 = CMPLX ( mij ( 1, idim ), mij ( 2, idim ), wp )
      a12 = a12 + weights ( idim ) * REAL ( CONJG ( z12 ) * ( z11 - z22 ), wp )
      b12 = b12 + weights ( idim ) * (  z12 * CONJG ( z12 ) -  &
            0.25_wp * ( z11 - z22 ) * ( CONJG ( z11 ) - CONJG ( z22 ) ) ) 
    END DO
    IF ( DABS ( b12 ) > 1.e-10_wp ) THEN
      ratio = -a12/b12
      theta = 0.25_wp * DATAN ( ratio )
    ELSEIF ( DABS ( b12 ) < 1.e-10_wp ) THEN
      b12 = 0.0_wp
      theta = 0.0_wp
    ELSE
      theta = .25_wp * pi
    ENDIF
! Check second derivative info
    d2 = a12 * DSIN ( 4._wp * theta ) - b12 * DCOS ( 4._wp * theta )
    IF ( d2 <=  0._wp ) theta = theta + 0.25_wp * pi
  END SUBROUTINE get_angle 
! *****************************************************************************
  SUBROUTINE check_tolerance ( zij, weights, tolerance )
    IMPLICIT NONE
    TYPE ( cp_complex_full_matrix_p_type ) :: zij ( : ) 
    REAL ( wp ), INTENT ( OUT ) :: tolerance
    REAL ( wp ), INTENT ( IN ) :: weights ( : )
! locals
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_full_matrix_type ), POINTER :: force
! compute gradient at t=0
    NULLIFY ( force )
    CALL cp_fm_create ( force, zij( 1 ) % matrix % matrix_struct, error = error  )
    CALL cp_fm_set_all ( force, 0._wp, error=error )
    CALL grad_at_0 ( zij, weights, force )
    CALL cp_fm_maxval ( force, tolerance, error = error )
    CALL cp_fm_release ( force )
  END SUBROUTINE check_tolerance
! *****************************************************************************
  SUBROUTINE rotate_orbitals ( rmat, vectors ) 
    IMPLICIT NONE
    TYPE ( cp_full_matrix_type ), POINTER :: rmat
    TYPE ( cp_full_matrix_type ), POINTER :: vectors
!locals
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_full_matrix_type ), POINTER :: wf
    INTEGER :: n, k
    NULLIFY ( wf )
    CALL cp_fm_create ( wf, vectors % matrix_struct, error = error  )
    CALL cp_fm_get_info ( vectors, nrow_global = n, ncol_global=k )
    CALL cp_fm_gemm("N", "N", n, k, k, 1.0_wp, vectors, rmat, 0.0_wp, wf )
    CALL cp_fm_release ( wf )
  END SUBROUTINE rotate_orbitals 
! *****************************************************************************
  SUBROUTINE grad_at_0 ( matrix_p, weights, matrix )
    IMPLICIT NONE
    TYPE ( cp_complex_full_matrix_p_type ) :: matrix_p ( : )
    TYPE ( cp_full_matrix_type ), POINTER :: matrix
    REAL ( KIND = wp ), INTENT ( IN ) :: weights ( : )
! Locals
    TYPE ( cp_error_type ) :: error
    COMPLEX ( KIND = wp ) :: zii, zjj, zij 
    COMPLEX ( KIND = wp ), DIMENSION(:,:), POINTER  :: diag
    REAL ( KIND = wp ) :: mii ( 2 ), mjj ( 2 ), mij ( 2 ), grad_ij
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
        grad_ij = 0.0_wp
        DO idim = 1, 3
          zii = diag(row_indices(istate),idim)
          zjj = diag(col_indices(jstate),idim)
          ! CALL cp_cfm_get_element(matrix_p(idim)%matrix,istate,istate,zii)
          ! CALL cp_cfm_get_element(matrix_p(idim)%matrix,jstate,jstate,zjj)
          CALL cp_cfm_get_element(matrix_p(idim)%matrix,istate,jstate,zij)
          grad_ij = grad_ij + weights ( idim ) * &
                    REAL ( 4.0_wp * CONJG ( zij ) * ( zjj - zii ), wp )
        END DO
        matrix % local_data ( istate, jstate ) = grad_ij
!       matrix % local_data ( jstate, istate ) = -grad_ij
      END DO
    END DO
    DEALLOCATE(diag)
!    CALL write_blacs_matrix(matrix,6,formatted=.true.)
!    stop
  END SUBROUTINE grad_at_0
! *****************************************************************************
END MODULE qs_localization_methods
