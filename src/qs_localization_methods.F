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
                                             cp_complex_full_matrix_type
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
    REAL ( KIND = wp ) :: tolerance, st, ct, mii ( 2, 3 ), mjj ( 2, 3 ), mij ( 2, 3 )
    REAL ( KIND = wp ) :: theta
    INTEGER :: istate, jstate, idim, nstate, isos,sweeps
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_full_matrix_type ), POINTER :: matrix_r, matrix_i

!Create the rotation matrix
    CALL cp_fm_create ( rmat, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
! Initialize rotation to the unit matrix
    CALL cp_fm_set_all ( rmat, 0._wp, 1._wp, error )
! get nstate
    CALL cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstate, error = error )
    tolerance = 1e10_wp
    sweeps = 0
    DO WHILE ( tolerance >= 1.0e-6_wp )
      sweeps = sweeps + 1
! ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE ROTATE 
      DO istate = 1, nstate
! get real and imaginary parts
        DO idim = 1, SIZE ( zij, 2 )
          matrix_r => zij ( 1, idim ) % matrix
          matrix_i => zij ( 2, idim ) % matrix
          CALL cp_fm_get_element ( matrix_r, istate, istate, mii ( 1, idim ) )
          CALL cp_fm_get_element ( matrix_i, istate, istate, mii ( 2, idim ) )
        END DO
        DO jstate = istate + 1, nstate
! get real and imaginary parts
          DO idim = 1, SIZE ( zij, 2 )
            matrix_r => zij ( 1, idim ) % matrix
            matrix_i => zij ( 2, idim ) % matrix
            CALL cp_fm_get_element ( matrix_r, istate, jstate, mij ( 1, idim ) )
            CALL cp_fm_get_element ( matrix_i, istate, jstate, mij ( 2, idim ) )
            CALL cp_fm_get_element ( matrix_r, jstate, jstate, mjj ( 1, idim ) )
            CALL cp_fm_get_element ( matrix_i, jstate, jstate, mjj ( 2, idim ) )
          END DO 
! get the optimal angle
          CALL get_angle ( mii, mjj, mij, zij, weights, theta )
          st = DSIN ( theta )
          ct = DCOS ( theta )
! rotate the ZIJ
          CALL rotate_zij ( istate, jstate, st, ct, zij )
! update the Unitary matrix, rmat
          CALL rotate_rmat ( istate, jstate, st, ct, rmat )
        END DO
      END DO
! Check tolerance
      CALL check_tolerance ( zij, weights, tolerance )
    END DO
! Rotate the orbitals
    CALL rotate_orbitals ( rmat, vectors ) 
! DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE DONE ROTATE 
  END SUBROUTINE jacobi_rotations
! *****************************************************************************
  SUBROUTINE rotate_zij ( istate, jstate, st, ct, zij )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: istate, jstate
    TYPE ( cp_full_matrix_p_type ) :: zij ( :, : )
    REAL ( KIND = wp ), INTENT ( IN ) :: st, ct
! Locals
    TYPE ( cp_error_type ) :: error
    INTEGER :: idim, nstate
    TYPE ( cp_full_matrix_type ), POINTER :: matrix_r, matrix_i
    TYPE ( cp_complex_full_matrix_type ), POINTER :: zij_cmplx
    COMPLEX ( KIND = wp ) :: st_cmplx
    COMPLEX ( KIND = wp ), POINTER :: x ( : ), y ( : )

#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

    st_cmplx = CMPLX ( st, 0._wp, wp )
    CALL cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstate, error = error )
    CALL cp_cfm_create ( zij_cmplx, zij ( 1, 1 ) % matrix % matrix_struct, error=error )
    DO  idim = 1, 3
      matrix_r => zij ( 1, idim ) % matrix
      matrix_i => zij ( 2, idim ) % matrix

      zij_cmplx % local_data = CMPLX ( matrix_r % local_data, matrix_i % local_data, wp )
   
      x => zij_cmplx % local_data ( :, istate )
      y => zij_cmplx % local_data ( :, jstate )
#if defined(__SCALAPACK)
      desc(:) = matrix_r%matrix_struct%descriptor(:)

      CALL pzrot( nstate , x, 1, 1, desc, 1, y, 1, 1, desc, 1, ct, st_cmplx )
#else
      CALL zrot( nstate, x, 1, y, 1, ct, st_cmplx )
#endif
      x => zij_cmplx % local_data ( istate, : )
      y => zij_cmplx % local_data ( jstate, : )
#if defined(__SCALAPACK)
      CALL pzrot( nstate , x, 1, 1, desc, 1, y, 1, 1, desc, 1, ct, st_cmplx )
#else
      CALL zrot( nstate, x, 1, y, 1, ct, st_cmplx )
#endif
      matrix_r % local_data = REAL ( zij_cmplx % local_data, wp )
      matrix_i % local_data = AIMAG ( zij_cmplx % local_data )
    END DO 
!    CALL write_blacs_matrix(ZIJ(1,1)%matrix,6,formatted=.true.)
    CALL cp_cfm_release ( zij_cmplx )
  END SUBROUTINE rotate_zij 
! *****************************************************************************
  SUBROUTINE rotate_rmat ( istate, jstate, st, ct, rmat )
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: istate, jstate
    TYPE ( cp_full_matrix_type ), POINTER :: rmat
    REAL ( KIND = wp ), INTENT ( IN ) :: ct, st
! Locals
    TYPE ( cp_error_type ) :: error
    INTEGER :: nstate
    TYPE ( cp_complex_full_matrix_type ), POINTER :: rmat_cmplx
    COMPLEX ( KIND = wp ) :: st_cmplx
    COMPLEX ( KIND = wp ), POINTER :: x ( : ), y ( : )
    REAL ( wp ) :: zi,zj
    INTEGER :: kstate

!    CALL cp_fm_get_info (  rmat, nrow_global = nstate, error = error )
!    DO kstate = 1, nstate
!      ZI=ct*rmat%local_data(kstate,istate)+st*rmat%local_data(kstate,jstate)
!      ZJ=-st*rmat%local_data(kstate,istate)+ct*rmat%local_data(kstate,jstate)
!      RMAT%local_data(kstate,istate)=ZI
!      RMAT%local_data(kstate,jstate)=ZJ
!    ENDDO

#if defined(__SCALAPACK)
    INTEGER, DIMENSION(9) :: desc
#endif

    st_cmplx = CMPLX ( st, 0._wp, wp )
    CALL cp_fm_get_info (  rmat, nrow_global = nstate, error = error )
    CALL cp_cfm_create ( rmat_cmplx, rmat % matrix_struct, error=error )
    rmat_cmplx % local_data = CMPLX ( rmat % local_data, 0._wp, wp )
   
    x => rmat_cmplx % local_data ( :, istate )
    y => rmat_cmplx % local_data ( :, jstate )

#if defined(__SCALAPACK)
    desc ( : ) = rmat_cmplx % matrix_struct % descriptor(:)

    CALL pzrot( nstate , x, 1, 1, desc, 1, y, 1, 1, desc, 1, ct, st_cmplx )
#else
    CALL zrot( nstate, x, 1, y, 1, ct, st_cmplx )
#endif
    rmat % local_data = REAL ( rmat_cmplx % local_data, wp )
!    CALL write_blacs_matrix(rmat,6,formatted=.true.)
    CALL cp_cfm_release ( rmat_cmplx, error=error )
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
    TYPE ( cp_full_matrix_p_type ) :: zij ( :, : ) 
    REAL ( wp ), INTENT ( OUT ) :: tolerance
    REAL ( wp ), INTENT ( IN ) :: weights ( : )
! locals
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_full_matrix_type ), POINTER :: force
! compute gradient at t=0
    NULLIFY ( force )
    CALL cp_fm_create ( force, zij( 1, 1 ) % matrix % matrix_struct, error = error  )
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
    TYPE ( cp_full_matrix_p_type ) :: matrix_p ( :, : )
    TYPE ( cp_full_matrix_type ), POINTER :: matrix
    REAL ( KIND = wp ), INTENT ( IN ) :: weights ( : )
! Locals
    TYPE ( cp_error_type ) :: error
    TYPE ( cp_full_matrix_type ), POINTER :: matrix_r, matrix_i
    COMPLEX ( KIND = wp ) :: zii, zjj, zij 
    REAL ( KIND = wp ) :: mii ( 2 ), mjj ( 2 ), mij ( 2 ), grad_ij
    INTEGER :: nrow_local, ncol_local, istate, jstate, idim

    CALL cp_fm_get_info ( matrix, nrow_local = nrow_local,  &
                   ncol_local = ncol_local, error = error )
    DO istate = 1, nrow_local
      DO jstate = 1, ncol_local
! get real and imaginary parts
        grad_ij = 0.0_wp
        DO idim = 1, 3
          matrix_r => matrix_p ( 1, idim ) % matrix  
          matrix_i => matrix_p ( 2, idim ) % matrix  
          mii ( 1 ) =  matrix_r % local_data ( istate, istate )
          mii ( 2 ) =  matrix_i % local_data ( istate, istate )
          mjj ( 1 ) =  matrix_r % local_data ( jstate, jstate )
          mjj ( 2 ) =  matrix_i % local_data ( jstate, jstate )
          mij ( 1 ) =  matrix_r % local_data ( istate, jstate )
          mij ( 2 ) =  matrix_i % local_data ( istate, jstate )
          zii = CMPLX ( mii ( 1 ), mii ( 2 ), wp )
          zjj = CMPLX ( mjj ( 1 ), mjj ( 2 ), wp )
          zij = CMPLX ( mij ( 1 ), mij ( 2 ), wp )
          grad_ij = grad_ij + weights ( idim ) * &
                    REAL ( 4.0_wp * CONJG ( zij ) * ( zjj - zii ), wp )
        END DO
        matrix % local_data ( istate, jstate ) = grad_ij
!        matrix % local_data ( jstate, istate ) = -grad_ij
      END DO
    END DO
!    CALL write_blacs_matrix(matrix,6,formatted=.true.)
!    stop
  END SUBROUTINE grad_at_0
! *****************************************************************************
END MODULE qs_localization_methods
