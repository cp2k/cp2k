!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE dkh_main
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_schur_product,&
                                             cp_fm_syrk,&
                                             cp_fm_transpose,&
                                             cp_fm_triangular_multiply,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_reduce
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_get,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: a_fine
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dkh_main'


  PUBLIC                               :: dkh_atom_transformation

CONTAINS

! *****************************************************************************
!> \brief 2th order DKH calculations
!>
!> \param qs_env ...
!> \param matrix_s ...
!> \param matrix_v ...
!> \param matrix_t ...
!> \param matrix_pVp ...
!> \param n ...
!> \param dkh_order ...
!> \param error ...
!> \par Literature
!>  M. Reiher, A. Wolf, J. Chem. Phys. 121 (2004) 10944-10956
!>  A. Wolf, M. Reiher, B. A. Hess, J. Chem. Phys. 117 (2002) 9215-9226
!>
!>\par Note
!>      based on subroutines for DKH1 to DKH5 by
!>       A. Wolf, M. Reiher, B. A. Hess
!>
!>  INPUT:
!>    qs_env (:)        The quickstep environment
!>    n                 Number of primitive gaussians
!>    matrix_s    (:,:) overlap matrix
!>    matrix_pVp  (:,:) pVp matrix
!>
!>  IN_OUT:
!>    matrix_v    (:,:) input: nonrelativistic potential energy matrix
!>                      output: (ev1+ev2)
!>    matrix_t    (:,:) input: kinetic energy matrix
!>                      output: kinetic part of hamiltonian
!>                      in position space
!>
!>  INTERNAL
!>    sinv (:,:) inverted, orthogonalized overlap matrix
!>    ev0t (:)   DKH-even0 matrix in T-basis
!>    e    (:)   e=SQRT(p^2c^2+c^4)
!>    eig  (:,:) eigenvectors of sinv' h sinv
!>    tt   (:)   eigenvalues of sinv' h sinv
!>    revt (:,:) reverse transformation matrix T-basis -> position space
!>    aa   (:)   kinematical factors f. DKH SQRT((c^2+e(i))/(2.0*e(i)))
!>    rr   (:)   kinematical factors f. DKH c/(c^2+e(i))
!>    vt   (:,:) non relativistic potential matrix in T-basis
!>    pvpt (:,:) pvp integral matrix in T-basis
!>    ev1t (:,:) DKH-even1 matrix in T-basis
!>    evt2 (:,:) DKH-even2 matrix in T-basis
!>    ev1  (:,:) DKH-even1 matrix in position space
!>    ev2  (:,:) DKH-even2 matrix in position space
!>    ove (:,:) scratch
!>    aux (:,:) scratch
!>    velit  velocity of light 137 a.u.
!>    prea   prefactor, 1/137^2
!>    con2   prefactor, 2/137^2
!>    con    prefactor, 137^2
!> \author
!>     Jens Thar, Barbara Kirchner: Uni Bonn (09/2006)
!>     Markus Reiher: ETH Zurich (09/2006)
!>
! *****************************************************************************
  SUBROUTINE DKH_full_transformation (qs_env,matrix_s,matrix_v,matrix_t,matrix_pVp,n,dkh_order,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: matrix_s, matrix_v, matrix_t, &
                                                matrix_pVp
    INTEGER, INTENT(IN)                      :: n, dkh_order
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'DKH_full_transformation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    REAL(KIND=dp)                            :: velit
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: aa, e, ev0t, rr, tt
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: matrix_full
    TYPE(cp_fm_type), POINTER :: matrix_aux, matrix_aux2, matrix_eig, &
      matrix_ev1, matrix_ev2, matrix_ev3, matrix_ev4, matrix_pe1p, &
      matrix_rev, matrix_se, matrix_sinv

    CALL timeset(routineN,handle)
    NULLIFY(blacs_env)
    !-----------------------------------------------------------------------
    !     Define velocity of light
    !-----------------------------------------------------------------------

    !     velit = 137.0359895_dp   ! 1/a_fine  physcon
    velit = 1._dp/a_fine

    !-----------------------------------------------------------------------
    !     Construct the matrix structure
    !-----------------------------------------------------------------------
    CALL get_qs_env(qs_env, blacs_env=blacs_env, error=error)
    CALL cp_fm_struct_create( fmstruct = matrix_full,&
         context = blacs_env,&
         nrow_global = n,&
         ncol_global = n,&
         error = error)

    !-----------------------------------------------------------------------
    !     Allocate some matrices
    !-----------------------------------------------------------------------

    ALLOCATE(e(n))
    ALLOCATE(aa(n))
    ALLOCATE(rr(n))
    ALLOCATE(tt(n))
    ALLOCATE(ev0t(n))

    NULLIFY(matrix_eig)
    NULLIFY(matrix_aux)
    NULLIFY(matrix_aux2)
    NULLIFY(matrix_rev)
    NULLIFY(matrix_se)
    NULLIFY(matrix_ev1)
    NULLIFY(matrix_ev2)
    NULLIFY(matrix_ev3)
    NULLIFY(matrix_ev4)
    NULLIFY(matrix_sinv)
    NULLIFY(matrix_pe1p)

    CALL cp_fm_create(matrix_eig,  matrix_full, error=error)
    CALL cp_fm_create(matrix_aux,  matrix_full, error=error)
    CALL cp_fm_create(matrix_aux2,  matrix_full, error=error)
    CALL cp_fm_create(matrix_rev,  matrix_full, error=error)
    CALL cp_fm_create(matrix_se,   matrix_full, error=error)
    CALL cp_fm_create(matrix_ev1,  matrix_full, error=error)
    CALL cp_fm_create(matrix_ev2,  matrix_full, error=error)
    CALL cp_fm_create(matrix_sinv, matrix_full, error=error)
    CALL cp_fm_create(matrix_ev3,  matrix_full, error=error)
    CALL cp_fm_create(matrix_ev4,  matrix_full, error=error)
    CALL cp_fm_create(matrix_pe1p,  matrix_full, error=error)


    !-----------------------------------------------------------------------
    !     Now with Cholesky decomposition
    !-----------------------------------------------------------------------

    CALL cp_fm_to_fm(matrix_s,matrix_sinv,error)
    CALL cp_fm_cholesky_decompose(matrix_sinv,n,error=error)

    !-----------------------------------------------------------------------
    !     Calculate matrix representation from nonrelativistic T matrix
    !-----------------------------------------------------------------------

    CALL cp_fm_cholesky_reduce(matrix_t,matrix_sinv,error=error)
    CALL cp_fm_syevd(matrix_t,matrix_eig,tt,error=error)

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of Hamiltonian in T-basis
    !-----------------------------------------------------------------------

    CALL kintegral (n,ev0t,tt,e,velit)

    !-----------------------------------------------------------------------
    !     Calculate reverse transformation matrix revt
    !-----------------------------------------------------------------------

    CALL cp_fm_to_fm(matrix_eig,matrix_rev,error)
    CALL cp_fm_triangular_multiply(matrix_sinv,matrix_rev,transpose_tr=.TRUE.,error=error)

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of the Hamiltonian
    !-----------------------------------------------------------------------

    CALL cp_fm_to_fm(matrix_rev,matrix_aux,error)
    CALL cp_fm_column_scale(matrix_aux,ev0t)
    CALL cp_gemm("N","T",n,n,n,1.0_dp,matrix_rev,matrix_aux,0.0_dp,matrix_t,error)

    !-----------------------------------------------------------------------
    !     Calculate kinematical factors for DKH
    !     only vectors present - will be done on every CPU
    !-----------------------------------------------------------------------

    DO i=1,n
       aa(i)=SQRT((velit*velit+e(i)) / (2.0_dp*e(i)))
       rr(i)=SQRT(velit*velit)/(velit*velit+e(i))
    END DO

    !-----------------------------------------------------------------------
    !     Transform v integrals to T-basis (v -> v(t))
    !-----------------------------------------------------------------------

    CALL cp_fm_cholesky_reduce(matrix_v,matrix_sinv,error=error)
    CALL cp_fm_upper_to_full(matrix_v,matrix_aux,error)
    CALL cp_gemm("T","N",n,n,n,1.0_dp,matrix_eig,matrix_v,0.0_dp,matrix_aux,error)
    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_aux,matrix_eig,0.0_dp,matrix_v,error)

    !-----------------------------------------------------------------------
    !     Transform pVp integrals to T-basis (pVp -> pVp(t))
    !-----------------------------------------------------------------------

    CALL cp_fm_cholesky_reduce(matrix_pVp,matrix_sinv,error=error)
    CALL cp_fm_upper_to_full(matrix_pVp,matrix_aux,error)
    CALL cp_gemm("T","N",n,n,n,1.0_dp,matrix_eig,matrix_pVp,0.0_dp,matrix_aux,error)
    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_aux,matrix_eig,0.0_dp,matrix_pVp,error)

    !-----------------------------------------------------------------------
    !     Calculate even1 in T-basis
    !-----------------------------------------------------------------------

    CALL even1(matrix_ev1,matrix_v,matrix_pvp,aa,rr,matrix_aux,matrix_aux2,error)

    !-----------------------------------------------------------------------
    !     Calculate even2 in T-basis
    !-----------------------------------------------------------------------

    CALL even2c (n,matrix_ev2,matrix_v,matrix_pVp,aa,rr,tt,e,matrix_aux,error)

    !-----------------------------------------------------------------------
    !     Calculate even3 in T-basis, only if requested
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.3) THEN
       CALL peven1p(n,matrix_pe1p,matrix_v,matrix_pvp,matrix_aux,matrix_aux2,aa,rr,tt,error)
       CALL even3b(n,matrix_ev3,matrix_ev1,matrix_pe1p,matrix_v,matrix_pvp,aa,rr,tt,e,matrix_aux,error)

       !-----------------------------------------------------------------------
       !     Transform even3 back to position space
       !-----------------------------------------------------------------------

       CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev3,0.0_dp,matrix_aux,error)
       CALL cp_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev3,error)

       !-----------------------------------------------------------------------
       !     Calculate even4 in T-basis, only if requested
       !-----------------------------------------------------------------------

       IF (dkh_order.ge.4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="DKH order greater than 3 not yet available", &
               error=error, error_level=cp_failure_level)
          !          CALL even4a(n,matrix_ev4%local_data,matrix_ev2%local_data,matrix_pe1p%local_data,matrix_v%local_data,&
          !                      matrix_pvp%local_data,aa,rr,tt,e)

          !-----------------------------------------------------------------------
          !     Transform even4 back to position space
          !-----------------------------------------------------------------------

          !        CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev4,0.0_dp,matrix_aux,error)
          !        CALL cp_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev4,error)

       END IF
    END IF

    !----------------------------------------------------------------------
    !     Transform even1 back to position space
    !----------------------------------------------------------------------

    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev1,0.0_dp,matrix_aux,error)
    CALL cp_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev1,error)

    !-----------------------------------------------------------------------
    !     Transform even2 back to position space
    !-----------------------------------------------------------------------

    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_rev,matrix_ev2,0.0_dp,matrix_aux,error)
    CALL cp_gemm("N","T",n,n,n,1.0_dp,matrix_aux,matrix_rev,0.0_dp,matrix_ev2,error)


    !-----------------------------------------------------------------------
    !     Calculate v in position space
    !-----------------------------------------------------------------------
    !
    CALL cp_fm_scale_and_add(1.0_dp,matrix_ev1,1.0_dp,matrix_ev2,error)
    CALL cp_fm_upper_to_full(matrix_ev1,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_ev1,matrix_v,error)
    IF(dkh_order.ge.3) THEN
       CALL cp_fm_scale_and_add(1.0_dp,matrix_v,1.0_dp,matrix_ev3,error)
       IF(dkh_order.ge.4) THEN
          CALL cp_fm_scale_and_add(1.0_dp,matrix_v,1.0_dp,matrix_ev4,error)
       END IF
    END IF

    !-----------------------------------------------------------------------

    CALL cp_fm_release(matrix_eig,  error=error)
    CALL cp_fm_release(matrix_aux,  error=error)
    CALL cp_fm_release(matrix_aux2, error=error)
    CALL cp_fm_release(matrix_rev,  error=error)
    CALL cp_fm_release(matrix_se,   error=error)
    CALL cp_fm_release(matrix_ev1,  error=error)
    CALL cp_fm_release(matrix_ev2,  error=error)
    CALL cp_fm_release(matrix_sinv, error=error)
    CALL cp_fm_release(matrix_ev3,  error=error)
    CALL cp_fm_release(matrix_ev4,  error=error)
    CALL cp_fm_release(matrix_pe1p, error=error)

    CALL cp_fm_struct_release(matrix_full, error=error)

    DEALLOCATE(ev0t,e,aa,rr,tt)

    CALL timestop(handle)

  END SUBROUTINE DKH_full_transformation

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param ev0t ...
!> \param tt ...
!> \param e ...
!> \param velit ...
! *****************************************************************************
  SUBROUTINE kintegral (n,ev0t,tt,e,velit)
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: ev0t
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tt
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: e
    REAL(KIND=dp), INTENT(IN)                :: velit

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: con, con2, prea, ratio, tv1, &
                                                tv2, tv3, tv4

    prea=1/(velit*velit)
    con2=prea+prea
    con=1.0_dp/prea

    DO i=1,n
       IF (tt(i).LT.0.0_dp) THEN
          WRITE(*,*) ' dkh_main.F | tt(',i,') = ',tt(i)
       END IF

       !       If T is sufficiently small, use series expansion to avoid
       !       cancellation, otherwise calculate SQRT directly

       ev0t(i)=tt(i)
       ratio=tt(i)/velit
       IF (ratio.LE.0.02_dp) THEN
          tv1=tt(i)
          tv2=-tv1*tt(i)*prea*0.5_dp
          tv3=-tv2*tt(i)*prea
          tv4=-tv3*tt(i)*prea*1.25_dp
          ev0t(i)=tv1+tv2+tv3+tv4
       ELSE
          ev0t(i)=con*(SQRT(1.0_dp+con2*tt(i))-1.0_dp)
       END IF
       e(i)=ev0t(i)+con
    END DO

    RETURN
  END SUBROUTINE kintegral

! *****************************************************************************
!> \brief ...
!> \param matrix_ev1 ...
!> \param matrix_v ...
!> \param matrix_pvp ...
!> \param aa ...
!> \param rr ...
!> \param matrix_aux ...
!> \param matrix_aux2 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE even1(matrix_ev1,matrix_v,matrix_pvp,aa,rr,matrix_aux,matrix_aux2,error)
    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     ev1  out  even1 output matrix                                    -
    !     v    in   potential matrix v in T-space                          -
    !     pvp  in   pvp matrix in T-space                                  -
    !     aa   in   A-factors (diagonal)                                   -
    !     rr   in   R-factors (diagonal)                                   -
    !                                                                      -
    !-----------------------------------------------------------------------
    TYPE(cp_fm_type), POINTER                :: matrix_ev1, matrix_v, &
                                                matrix_pVp
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr
    TYPE(cp_fm_type), POINTER                :: matrix_aux, matrix_aux2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL cp_fm_to_fm(matrix_v,matrix_aux,error)
    CALL cp_fm_column_scale(matrix_aux,aa)
    CALL cp_fm_transpose(matrix_aux,matrix_ev1,error)
    CALL cp_fm_column_scale(matrix_ev1,aa)

    CALL cp_fm_to_fm(matrix_pVp,matrix_aux,error)
    CALL cp_fm_column_scale(matrix_aux,aa)
    CALL cp_fm_column_scale(matrix_aux,rr)
    CALL cp_fm_transpose(matrix_aux,matrix_aux2,error)
    CALL cp_fm_column_scale(matrix_aux2,aa)
    CALL cp_fm_column_scale(matrix_aux2,rr)

    CALL cp_fm_scale_and_add(1.0_dp,matrix_ev1,1.0_dp,matrix_aux2,error)

    RETURN
  END SUBROUTINE even1

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param matrix_pe1p ...
!> \param matrix_v ...
!> \param matrix_pvp ...
!> \param matrix_aux ...
!> \param matrix_aux2 ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE peven1p(n,matrix_pe1p,matrix_v,matrix_pvp,matrix_aux,matrix_aux2,aa,rr,tt,error)

    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     n      in   dimension of matrices                                -
    !     pev1tp out  peven1p output matrix                                -
    !     vt     in   potential matrix v in T-space                        -
    !     pvpt   in   pvp matrix in T-space                                -
    !     aa     in   A-factors (diagonal)                                 -
    !     rr     in   R-factors (diagonal)                                 -
    !                                                                      -
    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_fm_type), POINTER                :: matrix_pe1p, matrix_v, &
                                                matrix_pvp, matrix_aux, &
                                                matrix_aux2
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    REAL(KIND=dp), DIMENSION(n)              :: vec_ar, vec_arrt
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: vec_full
    TYPE(cp_fm_type), POINTER                :: vec_a

    DO i=1,n
       vec_ar(i)=aa(i)*rr(i)
       vec_arrt(i)=vec_ar(i)*rr(i)*tt(i)
    END DO

    CALL cp_fm_struct_get(matrix_v%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = vec_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = 1,&
         error = error)

    NULLIFY(vec_a)
    CALL cp_fm_create(vec_a,     vec_full, error=error)

    CALL cp_fm_get_info(matrix_v, nrow_local=nrow_local, &
         row_indices=row_indices,error=error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = vec_arrt(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_aux2,error)
    CALL cp_fm_schur_product(matrix_v,matrix_aux,matrix_pe1p,error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = vec_ar(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_aux2,error)
    CALL cp_fm_schur_product(matrix_pvp,matrix_aux,matrix_aux2,error)

    CALL cp_fm_scale_and_add(4.0_dp,matrix_pe1p,1.0_dp,matrix_aux2,error)

    CALL cp_fm_release(vec_a,error=error)
    CALL cp_fm_struct_release(vec_full,error=error)

    RETURN
  END SUBROUTINE peven1p

  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param matrix_ev2 ...
!> \param matrix_v ...
!> \param matrix_pVp ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
!> \param e ...
!> \param matrix_aux ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE even2c (n,matrix_ev2,matrix_v,matrix_pVp,aa,rr,tt,e,matrix_aux,error)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified: 20.02.2002 - DKH2                 *
    !                                                                      *
    !     2nd order DK-approximation ( original DK-transformation with     *
    !                                       U = SQRT(1+W^2) + W        )   *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (6.2.2002)                                         *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !     2008       Jens Thar: transfer to CP2K                           *
    !                                                                      *
    !     ev2 = 1/2 [W1,O1]                                                *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n           in   Dimension of matrices                           *
    !     ev2         out  even2 output matrix = final result              *
    !     vv          in   potential v                                     *
    !     gg          in   pvp                                             *
    !     aa          in   A-Factors (DIAGONAL)                            *
    !     rr          in   R-Factors (DIAGONAL)                            *
    !     tt          in   Nonrel. kinetic Energy (DIAGONAL)               *
    !     e           in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)   *
    !     matrix_v    symmetric (n x n)-matrix containing (A V A)          *
    !     matrix_pvp  symmetric (n x n)-matrix containing (A P V P A)      *
    !     vh          symmetric (n x n)-matrix containing (A V~ A)         *
    !     pvph        symmetric (n x n)-matrix containing (A P V~ P A)     *
    !     w1o1        W1*O1 (2-component form)                             *
    !     o1w1        O1*W1 (2-component form)                             *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_fm_type), POINTER                :: matrix_ev2, matrix_v, &
                                                matrix_pVp
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: matrix_full
    TYPE(cp_fm_type), POINTER                :: matrix_apVpa, matrix_apVVpa, &
                                                matrix_aux2, matrix_ava, &
                                                matrix_avva

!     result  intermediate result of even2-calculation
!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH2
!-----------------------------------------------------------------------

    CALL cp_fm_struct_get(matrix_v%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = matrix_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = n,&
         error = error)

    NULLIFY(matrix_aux2)
    NULLIFY(matrix_ava)
    NULLIFY(matrix_avva)
    NULLIFY(matrix_apVpa)
    NULLIFY(matrix_apVVpa)

    CALL cp_fm_create(matrix_aux2,    matrix_full, error=error)
    CALL cp_fm_create(matrix_ava,     matrix_full, error=error)
    CALL cp_fm_create(matrix_avva,    matrix_full, error=error)
    CALL cp_fm_create(matrix_apVpa,   matrix_full, error=error)
    CALL cp_fm_create(matrix_apVVpa,  matrix_full, error=error)

    CALL cp_fm_to_fm(matrix_v,matrix_ava,error)
    CALL cp_fm_to_fm(matrix_v,matrix_avva,error)
    CALL cp_fm_to_fm(matrix_pVp,matrix_apVpa,error)
    CALL cp_fm_to_fm(matrix_pVp,matrix_apVVpa,error)

    !  Calculate  v = A V A:

    CALL mat_axa(matrix_v,matrix_ava,n,aa,matrix_aux,error)

    !  Calculate  pvp = A P V P A:

    CALL mat_arxra(matrix_pVp,matrix_apVpa,n,aa,rr,matrix_aux,error)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h(matrix_v,matrix_avva,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_avva,matrix_aux2,error)
    CALL mat_axa(matrix_aux2,matrix_avva,n,aa,matrix_aux,error)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h(matrix_pVp,matrix_apVVpa,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_apVVpa,matrix_aux2,error)
    CALL mat_arxra(matrix_aux2,matrix_apVVpa,n,aa,rr,matrix_aux,error)

    !  Calculate w1o1:

    CALL cp_gemm("N","N",n,n,n,-1.0_dp,matrix_apVVpa,matrix_ava,0.0_dp,matrix_aux2,error)
    CALL mat_muld(matrix_aux2,matrix_apVVpa,matrix_apVpa,n,  1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL mat_mulm(matrix_aux2,matrix_avva,  matrix_ava,n,    1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL cp_gemm("N","N",n,n,n,-1.0_dp,matrix_avva,matrix_apVpa,1.0_dp,matrix_aux2,error)

    !  Calculate o1w1 (already stored in ev2):

    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_apVpa,matrix_avva,0.0_dp,matrix_ev2,error)
    CALL mat_muld(matrix_ev2,matrix_apVpa,matrix_apVVpa,n,  -1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL mat_mulm(matrix_ev2,matrix_ava,  matrix_avva,n,    -1.0_dp,1.0_dp,tt,rr,matrix_aux,error)
    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_ava,matrix_apVVpa,1.0_dp,matrix_ev2,error)

    !-----------------------------------------------------------------------
    !     2.   1/2 [W1,O1] = 1/2 W1O1 -  1/2 O1W1
    !-----------------------------------------------------------------------

    CALL cp_fm_scale_and_add(-0.5_dp,matrix_ev2,0.5_dp,matrix_aux2,error)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!
    !-----------------------------------------------------------------------

    CALL cp_fm_release(matrix_aux2,   error=error)
    CALL cp_fm_release(matrix_ava,    error=error)
    CALL cp_fm_release(matrix_avva,   error=error)
    CALL cp_fm_release(matrix_apVpa,  error=error)
    CALL cp_fm_release(matrix_apVVpa, error=error)

    CALL cp_fm_struct_release(matrix_full, error=error)

!    WRITE (*,*) "CAW:  DKH2 with even2c (Alex)"
!    WRITE (*,*) "JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even2c

  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param matrix_ev3 ...
!> \param matrix_ev1 ...
!> \param matrix_pe1p ...
!> \param matrix_v ...
!> \param matrix_pVp ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
!> \param e ...
!> \param matrix_aux ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE even3b (n,matrix_ev3,matrix_ev1,matrix_pe1p,matrix_v,matrix_pVp,aa,rr,tt,e,matrix_aux,error)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified:  20.2.2002 - DKH3                 *
    !                                                                      *
    !     3rd order DK-approximation (generalised DK-transformation)       *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (7.2.2002)                                         *
    !                                                                      *
    !     ev3 = 1/2 [W1,[W1,E1]]                                           *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev3     out  even3 output matrix = final result                  *
    !     e1      in   E1 = even1-operator                                 *
    !     pe1p    in   pE1p                                                *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     result  intermediate result of even2-calculation
    !     vh      symmetric (n x n)-matrix containing (A V~ A)
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)
    !     e1      E1
    !     pe1p    pE1p
    !     w1w1    (W1)^2
    !     w1e1w1  W1*E1*W1
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,2)
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_fm_type), POINTER                :: matrix_ev3, matrix_ev1, &
                                                matrix_pe1p, matrix_v, &
                                                matrix_pVp
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: matrix_full
    TYPE(cp_fm_type), POINTER                :: matrix_apVVpa, matrix_aux2, &
                                                matrix_avva, matrix_w1e1w1, &
                                                matrix_w1w1

!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH3
!-----------------------------------------------------------------------

    CALL cp_fm_struct_get(matrix_v%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = matrix_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = n,&
         error = error)

    NULLIFY(matrix_aux2)
    NULLIFY(matrix_w1w1)
    NULLIFY(matrix_w1e1w1)
    NULLIFY(matrix_avva)
    NULLIFY(matrix_apVVpa)

    CALL cp_fm_create(matrix_aux2,    matrix_full, error=error)
    CALL cp_fm_create(matrix_w1w1,    matrix_full, error=error)
    CALL cp_fm_create(matrix_w1e1w1,  matrix_full, error=error)
    CALL cp_fm_create(matrix_avva,    matrix_full, error=error)
    CALL cp_fm_create(matrix_apVVpa,  matrix_full, error=error)

    CALL cp_fm_to_fm(matrix_v,matrix_avva,error)
    CALL cp_fm_to_fm(matrix_pVp,matrix_apVVpa,error)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h(matrix_v,matrix_avva,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_avva,matrix_aux2,error)
    CALL mat_axa(matrix_aux2,matrix_avva,n,aa,matrix_aux,error)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h(matrix_pVp,matrix_apVVpa,n,e,matrix_aux,error)
    CALL cp_fm_to_fm(matrix_apVVpa,matrix_aux2,error)
    CALL mat_arxra(matrix_aux2,matrix_apVVpa,n,aa,rr,matrix_aux,error)

    !  Calculate w1w1:

    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_apVVpa,matrix_avva,0.0_dp,matrix_w1w1,error)
    CALL mat_muld(matrix_w1w1,matrix_apVVpa,matrix_apVVpa,n,  -1.0_dp,1.0_dp,tt,rr,matrix_aux2,error)
    CALL mat_mulm(matrix_w1w1,matrix_avva,  matrix_avva,n,    -1.0_dp,1.0_dp,tt,rr,matrix_aux2,error)
    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_avva,matrix_apVVpa,1.0_dp,matrix_w1w1,error)

    !  Calculate w1e1w1: (warning: ev3 is scratch array)

    CALL mat_muld(matrix_aux,matrix_apVVpa,matrix_pe1p,n,  1.0_dp,0.0_dp,tt,rr,matrix_aux2,error)
    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_avva,matrix_pe1p,0.0_dp,matrix_aux2,error)
    CALL cp_gemm("N","N",n,n,n,1.0_dp,matrix_aux,matrix_avva,0.0_dp,matrix_w1e1w1,error)
    CALL mat_muld(matrix_w1e1w1,matrix_aux,matrix_apVVpa,n,  -1.0_dp,1.0_dp,tt,rr,matrix_ev3,error)
    CALL cp_gemm("N","N",n,n,n,-1.0_dp,matrix_aux2,matrix_avva,1.0_dp,matrix_w1e1w1,error)
    CALL mat_muld(matrix_w1e1w1,matrix_aux2,matrix_apVVpa,n,  1.0_dp,1.0_dp,tt,rr,matrix_ev3,error)

    !-----------------------------------------------------------------------
    !     2.   ev3 = 1/2 (W1^2)E1 + 1/2 E1(W1^2) - W1E1W1
    !-----------------------------------------------------------------------

    CALL cp_gemm("N","N",n,n,n,0.5_dp,matrix_w1w1,matrix_ev1,0.0_dp,matrix_ev3,error)
    CALL cp_gemm("N","N",n,n,n,0.5_dp,matrix_ev1,matrix_w1w1,1.0_dp,matrix_ev3,error)
    CALL cp_fm_scale_and_add(1.0_dp,matrix_ev3,-1.0_dp,matrix_w1e1w1,error)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!
    !-----------------------------------------------------------------------

    CALL cp_fm_release(matrix_aux2,   error=error)
    CALL cp_fm_release(matrix_avva,   error=error)
    CALL cp_fm_release(matrix_apVVpa, error=error)
    CALL cp_fm_release(matrix_w1w1,   error=error)
    CALL cp_fm_release(matrix_w1e1w1, error=error)

    CALL cp_fm_struct_release(matrix_full, error=error)

!    WRITE (*,*) "CAW:  DKH3 with even3b (Alex)"
!    WRITE (*,*) "JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even3b

  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param ev4 ...
!> \param e1 ...
!> \param pe1p ...
!> \param vv ...
!> \param gg ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
!> \param e ...
! *****************************************************************************
  SUBROUTINE even4a (n,ev4,e1,pe1p,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf,   last modified: 25.02.2002   --   DKH4          *
    !                                                                      *
    !     4th order DK-approximation (scalar = spin-free)                  *
    !                                                                      *
    !     Version: 1.2  (25.2.2002) :  Elegant (short) way of calculation  *
    !                                  included                            *
    !              1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (8.2.2002)                                         *
    !                                                                      *
    !     ev4  =  1/2 [W2,[W1,E1]] + 1/8 [W1,[W1,[W1,O1]]]  =              *
    !                                                                      *
    !          =      sum_1        +         sum_2                         *
    !                                                                      *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                (not working yet)                                     *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev4     out  even4 output matrix = final result                  *
    !     e1     in   E1                                                   *
    !     pe1p   in   p(E1)p                                               *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     v       symmetric (n x n)-matrix containing (A V A)              *
    !     pvp     symmetric (n x n)-matrix containing (A P V P A)          *
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         *
    !     w1w1    (W1)^2                                                   *
    !     w1o1    W1*O1      (2-component formulation)                     *
    !     o1w1    O1*W1      (2-component formulation)                     *
    !     e1      symmetric (n x n)-matrix containing E1                   *
    !     pe1p    symmetric (n x n)-matrix containing p(E1)p               *
    !     sum_i   2 addends defined above  (i=1,2)                         *
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,..,4)            *
    !     scrh_i  temp. (n x n)-scr.-mat. with energy-denom. (i=1,..,4)    *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev4
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e1, pe1p, vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: o1w1, pvp, pvph, scr_1, &
      scr_2, scr_3, scr_4, scrh_1, scrh_2, scrh_3, scrh_4, sum_1, sum_2, v, &
      vh, w1o1, w1w1

!C-----------------------------------------------------------------------
!C     1.   General Structures and Patterns for DKH4
!C-----------------------------------------------------------------------

    ALLOCATE(v(n,n))
    ALLOCATE(pVp(n,n))
    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    v=0.0_dp
    pVp=0.0_dp
    vh=0.0_dp
    pVph=0.0_dp
    v(1:n,1:n)=vv(1:n,1:n)
    vh(1:n,1:n)=vv(1:n,1:n)
    pvp(1:n,1:n)=gg(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev4=0.0_dp
    !  Calculate  v = A V A:

    !     CALL mat_axa(v,n,aa)

    !  Calculate  pvp = A P V P A:

    !     CALL mat_arxra(pvp,n,aa,rr)

    !  Calculate  vh = A V~ A:

    !     CALL mat_1_over_h(vh,n,e)
    !     CALL mat_axa(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    !     CALL mat_1_over_h(pvph,n,e)
    !     CALL mat_arxra(pvph,n,aa,rr)


    !  Create/Initialize necessary matrices:
    ALLOCATE(w1w1(n,n))
    w1w1 = 0.0_dp
    ALLOCATE(w1o1(n,n))
    w1o1 = 0.0_dp
    ALLOCATE(o1w1(n,n))
    o1w1 = 0.0_dp
    ALLOCATE(sum_1(n,n))
    sum_1 = 0.0_dp
    ALLOCATE(sum_2(n,n))
    sum_2 = 0.0_dp
    ALLOCATE(scr_1(n,n))
    scr_1 = 0.0_dp
    ALLOCATE(scr_2(n,n))
    scr_2 = 0.0_dp
    ALLOCATE(scr_3(n,n))
    scr_3 = 0.0_dp
    ALLOCATE(scr_4(n,n))
    scr_4 = 0.0_dp
    ALLOCATE(scrh_1(n,n))
    scrh_1 = 0.0_dp
    ALLOCATE(scrh_2(n,n))
    scrh_2 = 0.0_dp
    ALLOCATE(scrh_3(n,n))
    scrh_3 = 0.0_dp
    ALLOCATE(scrh_4(n,n))
    scrh_4 = 0.0_dp

    !  Calculate w1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,vh,n,0.0_dp,w1w1,n)
    !      CALL mat_muld(w1w1,pvph,pvph,n, -1.0_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(w1w1,vh,  vh,n,   -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pvph,n,1.0_dp,w1w1,n)

    !  Calculate w1o1:
    CALL dgemm("N","N",n,n,n,-1.0_dp,pvph,n,v,n,0.0_dp,w1o1,n)
    !      CALL mat_muld(w1o1,pvph,pvp,n,  1.0_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(w1o1,vh,  v,n,    1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,vh,n,pvp,n,1.0_dp,w1o1,n)
    !  Calculate o1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvp,n,vh,n,0.0_dp,o1w1,n)
    !      CALL mat_muld(o1w1,pvp,pvph,n,  -1.0_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(o1w1,v,  vh,n,    -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,v,n,pvph,n,1.0_dp,o1w1,n)

    !-----------------------------------------------------------------------
    !   2. sum_1 = 1/2 [W2,[W1,E1]] = 1/2 (W2W1E1 - W2E1W1 - W1E1W2 + E1W1W2)
    !-----------------------------------------------------------------------

    !  scr_i & scrh_i  for steps 2a (W2W1E1)  and 2b (W2E1W1):

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scr_1,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    !      CALL mat_muld(scr_4, pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)

    !      CALL mat_muld(scrh_1,pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    !      CALL mat_1_over_h(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scrh_2,n)
    !      CALL mat_1_over_h(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scrh_3,n)
    !      CALL mat_1_over_h(scrh_3,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scrh_4,n)
    !      CALL mat_1_over_h(scrh_4,n,e)

    !  2a)  sum_1 = 1/2 W2W1E1               ( [1]-[8] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_1,n,scr_1,n,0.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scrh_1,scr_2,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_2,n,scr_1,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scrh_2,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_3,n,scr_1,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scrh_3,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(sum_1,scrh_4,scr_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_2,n,1.0_dp,sum_1,n)


    !  2b)  sum_1 = - 1/2 W2E1W1 (+ sum_1)   ( [9]-[16] )

    !      CALL mat_muld(sum_1,scrh_1,scr_3,n,-0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_1,scr_4,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_2,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_2,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_3,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scrh_3,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_4,n,scr_4,n,1.0_dp,sum_1,n)


    !  scr_i & scrh_i  for steps 2c (W1E1W2)  and 2d (E1W1W2):

    !      CALL mat_muld(scr_1, pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scr_3,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scr_4,n)

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scrh_1,n)
    !      CALL mat_1_over_h(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scrh_2,n)
    !      CALL mat_1_over_h(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    !      CALL mat_1_over_h(scrh_3,n,e)
    !      CALL mat_muld(scrh_4,pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)
    !      CALL mat_1_over_h(scrh_4,n,e)

    !  2c)  sum_1 = - 1/2 W1E1W2 (+ sum_1)   ( [17]-[24] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scr_1,n,scrh_1,n,0.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_1,scrh_2,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_2,n,scrh_1,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_2,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_1,scrh_3,n,-0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_1,scrh_4,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_2,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_2,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)

    !  2d)  sum_1 = 1/2 E1W1W2 (+ sum_1)     ( [25]-[32] )

    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_3,n,scrh_1,n,0.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_3,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_mulm(sum_1,scr_4,scrh_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_2,n,1.0_dp,sum_1,n)
    !      CALL mat_muld(sum_1,scr_3,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    !      CALL mat_muld(sum_1,scr_3,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scr_4,n,scrh_4,n,1.0_dp,sum_1,n)


    !-----------------------------------------------------------------------
    !   3.  sum_2 = 1/8 [W1,[W1,[W1,O1]]] =
    !
    !             = 1/8 ( (W1^3)O1 - 3(W1^2)O1W1 + 3 W1O1(W1^2) - O1(W1^3) )
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,0.125_dp,w1w1,n,w1o1,n,0.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.375_dp,w1w1,n,o1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,0.375_dp,w1o1,n,w1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.125_dp,o1w1,n,w1w1,n,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   4.  result = sum_1 + sum_2
    !-----------------------------------------------------------------------

    CALL mat_add(ev4,1.0_dp,sum_1,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   5. Finish up the stuff!!
    !-----------------------------------------------------------------------

    DEALLOCATE(v,pvp,vh,pvph,w1w1,w1o1,o1w1,sum_1,sum_2)
    DEALLOCATE(scr_1,scr_2,scr_3,scr_4,scrh_1,scrh_2,scrh_3,scrh_4)

!    WRITE (*,*) "CAW:  DKH4 with even4a (Alex)"
!    WRITE (*,*) "JT:   Now available in cp2k"

    RETURN
  END SUBROUTINE even4a

  !-----------------------------------------------------------------------
  !                                                                      -
  !     Matrix routines for DKH-procedure                                -
  !     Alexander Wolf                                                   -
  !     modifed: Jens Thar: Mem manager deleted                          -
  !     This file contains the                                           -
  !      following subroutines:                                          -
  !                                 1. mat_1_over_h                      -
  !                                 2. mat_axa                           -
  !                                 3. mat_arxra                         -
  !                                 4. mat_mulm                          -
  !                                 5. mat_muld                          -
  !                                 6. mat_add                           -
  !                                                                      -
  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param matrix_p ...
!> \param matrix_pp ...
!> \param n ...
!> \param e ...
!> \param matrix_aux ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mat_1_over_h (matrix_p,matrix_pp,n,e,matrix_aux,error)

    !***********************************************************************
    !                                                                      *
    !   2. SR mat_1_over_h: Transform matrix p into matrix p/(e(i)+e(j))   *
    !                                                                      *
    !   p    in  REAL(:,:) :   matrix p                                    *
    !   e    in  REAL(:)   :   rel. energy (diagonal)                      *
    !   n    in  INTEGER                                                   *
    !                                                                      *
    !***********************************************************************

    TYPE(cp_fm_type), POINTER                :: matrix_p, matrix_pp
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, j, ncol_local, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices

    CALL cp_fm_get_info(matrix_aux, nrow_local=nrow_local, ncol_local=ncol_local,&
         row_indices=row_indices,col_indices=col_indices,error=error)

    DO j=1,ncol_local
       DO i=1,nrow_local
          matrix_aux%local_data(i,j)=1/(e(row_indices(i))+e(col_indices(j)))
       ENDDO
    ENDDO


    CALL cp_fm_schur_product(matrix_p,matrix_aux,matrix_pp,error)

    RETURN

  END SUBROUTINE mat_1_over_h
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param matrix_x ...
!> \param matrix_axa ...
!> \param n ...
!> \param a ...
!> \param matrix_aux ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mat_axa (matrix_x,matrix_axa,n,a,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   3. SR mat_axa: Transform matrix p into matrix  a*p*a               *
    !C                                                                      *
    !C   p    in  REAL(:,:):   matrix p                                     *
    !C   a    in  REAL(:)  :   A-factors (diagonal)                         *
    !CJT n    in  INTEGER  :   dimension of matrix p                        *
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_axa
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: vec_full
    TYPE(cp_fm_type), POINTER                :: vec_a

    CALL cp_fm_struct_get(matrix_x%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = vec_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = 1,&
         error = error)

    NULLIFY(vec_a)
    CALL cp_fm_create(vec_a,     vec_full, error=error)

    CALL cp_fm_get_info(matrix_x, nrow_local=nrow_local, &
         row_indices=row_indices,error=error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = a(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_axa,error)
    CALL cp_fm_schur_product(matrix_x,matrix_aux,matrix_axa,error)

    !     DO i=1,n
    !       DO j=1,n
    !          p(i,j)=p(i,j)*a(i)*a(j)
    !       ENDDO
    !     ENDDO


    CALL cp_fm_release(vec_a,error)
    CALL cp_fm_struct_release(vec_full,error)

    RETURN
  END SUBROUTINE mat_axa
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param matrix_x ...
!> \param matrix_axa ...
!> \param n ...
!> \param a ...
!> \param r ...
!> \param matrix_aux ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mat_arxra (matrix_x,matrix_axa,n,a,r,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   4. SR mat_arxra: Transform matrix p into matrix  a*r*p*r*a         *
    !C                                                                      *
    !C   p    in  REAL(:,:) :   matrix p                                    *
    !C   a    in  REAL(:)   :   A-factors (diagonal)                        *
    !C   r    in  REAL(:)   :   R-factors (diagonal)                        *
    !C   n    in  INTEGER   :   dimension of matrix p                       *
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_x, matrix_axa
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, r
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, nrow_local
    INTEGER, DIMENSION(:), POINTER           :: row_indices
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: vec_full
    TYPE(cp_fm_type), POINTER                :: vec_a

    CALL cp_fm_struct_get(matrix_x%matrix_struct,context=context, error = error)
    CALL cp_fm_struct_create( fmstruct = vec_full,&
         context = context,&
         nrow_global = n,&
         ncol_global = 1,&
         error = error)

    CALL cp_fm_get_info(matrix_aux, nrow_local=nrow_local, &
         row_indices=row_indices,error=error)


    NULLIFY(vec_a)
    CALL cp_fm_create(vec_a,     vec_full, error=error)

    DO i=1,nrow_local
       vec_a%local_data(i,1) = a(row_indices(i))*r(row_indices(i))
    END DO

    CALL cp_fm_syrk('U','N',1,1.0_dp,vec_a,1,1,0.0_dp,matrix_aux,error)
    CALL cp_fm_upper_to_full(matrix_aux,matrix_axa,error)
    CALL cp_fm_schur_product(matrix_x,matrix_aux,matrix_axa,error)


    CALL cp_fm_release(vec_a,error)
    CALL cp_fm_struct_release(vec_full,error)

    RETURN
  END SUBROUTINE mat_arxra
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param matrix_p ...
!> \param matrix_q ...
!> \param matrix_r ...
!> \param n ...
!> \param alpha ...
!> \param beta ...
!> \param t ...
!> \param rr ...
!> \param matrix_aux ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mat_mulm (matrix_p,matrix_q,matrix_r,n,alpha,beta,t,rr,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   5. SR mat_mulm:  Multiply matrices according to:                   *
    !C                                                                      *
    !C                      p = alpha*q*(..P^2..)*r + beta*p                *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,.):   matrix r                                  *
    !C   n      in   INTEGER  :   dimension n of matrices                   *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_p, matrix_q, matrix_r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(n)              :: vec

    CALL cp_fm_to_fm(matrix_q,matrix_aux,error)

    DO i=1,n
       vec(i)=2.0_dp*t(i)*rr(i)*rr(i)
    END DO
    CALL cp_fm_column_scale(matrix_aux,vec)

    CALL cp_gemm("N","N",n,n,n,alpha,matrix_aux,matrix_r,beta,matrix_p,error)
    RETURN

  END SUBROUTINE mat_mulm
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param matrix_p ...
!> \param matrix_q ...
!> \param matrix_r ...
!> \param n ...
!> \param alpha ...
!> \param beta ...
!> \param t ...
!> \param rr ...
!> \param matrix_aux ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mat_muld (matrix_p,matrix_q,matrix_r,n,alpha,beta,t,rr,matrix_aux,error)

    !C***********************************************************************
    !C                                                                      *
    !C   16. SR mat_muld:  Multiply matrices according to:                  *
    !C                                                                      *
    !C                      p = alpha*q*(..1/P^2..)*r + beta*p              *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,:):   matrix r                                  *
    !C   n      in   INTEGER  :   Dimension of all matrices                 *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    TYPE(cp_fm_type), POINTER                :: matrix_p, matrix_q, matrix_r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr
    TYPE(cp_fm_type), POINTER                :: matrix_aux
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(n)              :: vec

    CALL cp_fm_to_fm(matrix_q,matrix_aux,error)

    DO i=1,n
       vec(i)=0.5_dp/(t(i)*rr(i)*rr(i))
    END DO

    CALL cp_fm_column_scale(matrix_aux,vec)

    CALL cp_gemm("N","N",n,n,n,alpha,matrix_aux,matrix_r,beta,matrix_p,error)

    RETURN

  END SUBROUTINE mat_muld
  !C-----------------------------------------------------------------------


! *****************************************************************************
!> \brief ...
!> \param s ...
!> \param v ...
!> \param h ...
!> \param pVp ...
!> \param n ...
!> \param dkh_order ...
! *****************************************************************************
  SUBROUTINE DKH_atom_transformation (s,v,h,pVp,n,dkh_order)

    !-----------------------------------------------------------------------
    !                                                                      *
    !  INPUT:                                                              *
    !    n          Number of primitive gaussians                          *
    !    s    (:,:) overlap matrix                                         *
    !    pVp  (:,:) pVp matrix                                             *
    !                                                                      *
    !  IN_OUT:                                                             *
    !    v    (:,:) input: nonrelativistic potential energy matrix         *
    !               output: (ev1+ev2)                                      *
    !    h    (:,:) input: kinetic energy matrix                           *
    !               output: kinetic part of hamiltonian in position space  *
    !                                                                      *
    !  INTERNAL                                                            *
    !    sinv (:,:) inverted, orthogonalized overlap matrix                *
    !    ev0t (:)   DKH-even0 matrix in T-basis                            *
    !    e    (:)   e=SQRT(p^2c^2+c^4)                                     *
    !    eig  (:,:) eigenvectors of sinv' h sinv                           *
    !    tt   (:)   eigenvalues of sinv' h sinv                            *
    !    revt (:,:) reverse transformation matrix T-basis -> position space*
    !    aa   (:)   kinematical factors f. DKH SQRT((c^2+e(i))/(2.0*e(i))) *
    !    rr   (:)   kinematical factors f. DKH c/(c^2+e(i))                *
    !    vt   (:,:) non relativistic potential matrix in T-basis           *
    !    pvpt (:,:) pvp integral matrix in T-basis                         *
    !    ev1t (:,:) DKH-even1 matrix in T-basis                            *
    !    evt2 (:,:) DKH-even2 matrix in T-basis                            *
    !    ev1  (:,:) DKH-even1 matrix in position space                     *
    !    ev2  (:,:) DKH-even2 matrix in position space                     *
    !    ove (:,:) scratch                                                 *
    !    aux (:,:) scratch                                                 *
    !    velit  velocity of light 137 a.u.                                 *
    !    prea   prefactor, 1/137^2                                         *
    !    con2   prefactor, 2/137^2                                         *
    !    con    prefactor, 137^2                                           *
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: s, v, h, pVp
    INTEGER, INTENT(IN)                      :: n, dkh_order

    INTEGER                                  :: i, j, k
    REAL(KIND=dp)                            :: velit
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: aa, e, ev0t, rr, tt
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: aux, eig, ev1, ev1t, ev2, &
                                                ev2t, ev3, ev3t, ev4, ev4t, &
                                                ove, pev1tp, pVpt, revt, &
                                                sinv, vt

    IF (dkh_order.lt.0) RETURN

    !CAW  pp: p^2-values (in momentum-space), stored as matrix!!
    !-----------------------------------------------------------------------
    !     Define velocity of light
    !-----------------------------------------------------------------------
    !      velit = 137.0359895_dp

    velit = 1._dp/a_fine

    !-----------------------------------------------------------------------
    !     Allocate some matrices
    !-----------------------------------------------------------------------

    ALLOCATE(eig(n,n))
    ALLOCATE(sinv(n,n))
    ALLOCATE(revt(n,n))
    ALLOCATE(aux(n,n))
    ALLOCATE(ove(n,n))
    ALLOCATE(ev0t(n))
    ALLOCATE(e(n))
    ALLOCATE(aa(n))
    ALLOCATE(rr(n))
    ALLOCATE(tt(n))
    ALLOCATE(ev1t(n,n))
    ALLOCATE(ev2t(n,n))
    ALLOCATE(ev3t(n,n))
    ALLOCATE(ev4t(n,n))
    ALLOCATE(vt(n,n))
    ALLOCATE(pVpt(n,n))
    ALLOCATE(pev1tp(n,n))
    ALLOCATE(ev1(n,n))
    ALLOCATE(ev2(n,n))
    ALLOCATE(ev3(n,n))
    ALLOCATE(ev4(n,n))

    !-----------------------------------------------------------------------
    !     Schmidt-orthogonalize overlap matrix
    !-----------------------------------------------------------------------

    CALL sog (n,s,sinv)

    !-----------------------------------------------------------------------
    !     Calculate matrix representation from nonrelativistic T matrix
    !-----------------------------------------------------------------------

    CALL dkh_diag ( h,n,eig,tt,sinv,aux,0 )

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of Hamiltonian in T-basis
    !-----------------------------------------------------------------------

    CALL kintegral_a (n,ev0t,tt,e,velit)

    !-----------------------------------------------------------------------
    !     Calculate reverse transformation matrix revt
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,1.0_dp,sinv,n,eig,n,0.0_dp,aux,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,s,n,aux,n,0.0_dp,revt,n)

    !-----------------------------------------------------------------------
    !     Calculate kinetic part of the Hamiltonian
    !-----------------------------------------------------------------------

    h = 0.0_dp
    DO i=1,n
       DO j=1,i
          DO k=1,n
             h(i,j)=h(i,j)+revt(i,k)*revt(j,k)*ev0t(k)
             h(j,i)=h(i,j)
          END DO
       END DO
    END DO

    !-----------------------------------------------------------------------
    !     Calculate kinematical factors for DKH
    !-----------------------------------------------------------------------

    DO i=1,n
       aa(i)=SQRT((velit*velit+e(i)) / (2.0_dp*e(i)))
       rr(i)=SQRT(velit*velit)/(velit*velit+e(i))
    END DO

    !-----------------------------------------------------------------------
    !     Transform v integrals to T-basis (v -> vt)
    !-----------------------------------------------------------------------

    CALL trsm(v,sinv,ove,n,aux)
    CALL trsm(ove,eig,vt,n,aux)

    !-----------------------------------------------------------------------
    !     Transform pVp integrals to T-basis (pVp -> pVpt)
    !-----------------------------------------------------------------------

    CALL trsm(pVp,sinv,ove,n,aux)
    CALL trsm(ove,eig,pVpt,n,aux)

    !-----------------------------------------------------------------------
    !     Calculate even1 in T-basis
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.1) THEN
       CALL even1_a(n,ev1t,vt,pvpt,aa,rr)

       !----------------------------------------------------------------------
       !     Transform even1 back to position space
       !----------------------------------------------------------------------

       CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev1t,n,0.0_dp,aux,n)
       CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev1,n)
    END IF

    !-----------------------------------------------------------------------
    !     Calculate even2 in T-basis
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.2) THEN
       CALL even2c_a (n,ev2t,vt,pvpt,aa,rr,tt,e)

       !-----------------------------------------------------------------------
       !     Transform even2 back to position space
       !-----------------------------------------------------------------------

       aux=0.0_dp
       CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev2t,n,0.0_dp,aux,n)
       CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev2,n)
    END IF

    !-----------------------------------------------------------------------
    !     Calculate even3 in T-basis, only if requested
    !-----------------------------------------------------------------------

    IF (dkh_order.ge.3) THEN
       CALL peven1p_a(n,pev1tp,vt,pvpt,aa,rr,tt)
       CALL even3b_a(n,ev3t,ev1t,pev1tp,vt,pvpt,aa,rr,tt,e)

       !-----------------------------------------------------------------------
       !     Transform even3 back to position space
       !-----------------------------------------------------------------------
       aux=0.0_dp
       CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev3t,n,0.0_dp,aux,n)
       CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev3,n)

       !-----------------------------------------------------------------------
       !     Calculate even4 in T-basis, only if requested
       !-----------------------------------------------------------------------

       IF (dkh_order.ge.4) THEN
          CALL even4a_a(n,ev4t,ev1t,pev1tp,vt,pvpt,aa,rr,tt,e)

          !-----------------------------------------------------------------------
          !     Transform even4 back to position space
          !-----------------------------------------------------------------------
          aux=0.0_dp
          CALL dgemm("N","N",n,n,n,1.0_dp,revt,n,ev4t,n,0.0_dp,aux,n)
          CALL dgemm("N","T",n,n,n,1.0_dp,aux,n,revt,n,0.0_dp,ev4,n)
       END IF
    END IF

    IF (dkh_order.ge.4) THEN
       STOP "DKH 4"
    END IF
    !-----------------------------------------------------------------------
    !     Calculate v in position space
    !-----------------------------------------------------------------------

    IF(dkh_order.ge.1) THEN
       CALL mat_add2(v,0.0_dp,1.0_dp,ev1,n)
    END IF
    IF(dkh_order.ge.2) THEN
       CALL mat_add2(v,1.0_dp,1.0_dp,ev2,n)
    END IF
    IF(dkh_order.ge.3) THEN
       CALL mat_add2(v,1.0_dp,1.0_dp,ev3,n)
    END IF
    IF(dkh_order.ge.4) THEN
       CALL mat_add2(v,1.0_dp,1.0_dp,ev4,n)
    END IF

    !-----------------------------------------------------------------------

    DEALLOCATE(eig,sinv,revt,ove,aux,vt,pVpt,ev1,ev2,ev3,ev4,ev1t,ev2t,ev3t,ev4t,pev1tp)
    DEALLOCATE(ev0t,e,aa,rr,tt)

    RETURN
  END SUBROUTINE dkh_atom_transformation

  !-----------------------------------------------------------------------


! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param ev0t ...
!> \param tt ...
!> \param e ...
!> \param velit ...
! *****************************************************************************
  SUBROUTINE kintegral_a (n,ev0t,tt,e,velit)

    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: ev0t
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: tt
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: e
    REAL(KIND=dp), INTENT(IN)                :: velit

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: con, con2, prea, ratio, tv1, &
                                                tv2, tv3, tv4

!-----------------------------------------------------------------------

    DO i=1,n
       IF (tt(i).LT.0.0_dp) THEN
          WRITE(*,*) ' dkh_main.F | tt(',i,') = ',tt(i)
       END IF

       !       Calculate some constants

       prea=1/(velit*velit)
       con2=prea+prea
       con=1.0_dp/prea

       !       If T is sufficiently small, use series expansion to avoid
       !       cancellation, otherwise calculate SQRT directly

       ev0t(i)=tt(i)
       ratio=tt(i)/velit
       IF (ratio.LE.0.02_dp) THEN
          tv1=tt(i)
          tv2=-tv1*tt(i)*prea/2.0_dp
          tv3=-tv2*tt(i)*prea
          tv4=-tv3*tt(i)*prea*1.25_dp
          ev0t(i)=tv1+tv2+tv3+tv4
       ELSE
          ev0t(i)=con*(SQRT(1.0_dp+con2*tt(i))-1.0_dp)
       END IF
       e(i)=ev0t(i)+con
    END DO

    RETURN
  END SUBROUTINE kintegral_a


  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param ev1t ...
!> \param vt ...
!> \param pvpt ...
!> \param aa ...
!> \param rr ...
! *****************************************************************************
  SUBROUTINE even1_a(n,ev1t,vt,pvpt,aa,rr)

    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     n    in   dimension of matrices                                  -
    !     ev1t out  even1 output matrix                                    -
    !     vt   in   potential matrix v in T-space                          -
    !     pvpt in   pvp matrix in T-space                                  -
    !     aa   in   A-factors (diagonal)                                   -
    !     rr   in   R-factors (diagonal)                                   -
    !                                                                      -
    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev1t
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: vt, pvpt
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr

    INTEGER                                  :: i, j

!-----------------------------------------------------------------------

    DO i=1,n
       DO j=1,i
          ev1t(i,j)=vt(i,j)*aa(i)*aa(j)+pVpt(i,j)*aa(i)*rr(i)*aa(j)*rr(j)
          ev1t(j,i)=ev1t(i,j)
       END DO
    END DO

    RETURN
  END SUBROUTINE even1_a

  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param pev1tp ...
!> \param vt ...
!> \param pvpt ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
! *****************************************************************************
  SUBROUTINE peven1p_a(n,pev1tp,vt,pvpt,aa,rr,tt)

    !-----------------------------------------------------------------------
    !                                                                      -
    !     1st order DKH-approximation                                      -
    !                                                                      -
    !     n      in   dimension of matrices                                -
    !     pev1tp out  peven1p output matrix                                -
    !     vt     in   potential matrix v in T-space                        -
    !     pvpt   in   pvp matrix in T-space                                -
    !     aa     in   A-factors (diagonal)                                 -
    !     rr     in   R-factors (diagonal)                                 -
    !                                                                      -
    !-----------------------------------------------------------------------


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pev1tp
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: vt, pvpt
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt

    INTEGER                                  :: i, j

!-----------------------------------------------------------------------

    DO i=1,n
       DO j=1,i
          pev1tp(i,j)=4.0_dp*vt(i,j)*aa(i)*aa(j)*rr(i)*rr(i)*rr(j)*rr(j)*tt(i)*tt(j)+&
               pVpt(i,j)*aa(i)*rr(i)*aa(j)*rr(j)
          pev1tp(j,i)=pev1tp(i,j)
       END DO
    END DO

    RETURN
  END SUBROUTINE peven1p_a

  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param ev2 ...
!> \param vv ...
!> \param gg ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
!> \param e ...
! *****************************************************************************
  SUBROUTINE even2c_a (n,ev2,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified: 20.02.2002 - DKH2                 *
    !                                                                      *
    !     2nd order DK-approximation ( original DK-transformation with     *
    !                                       U = SQRT(1+W^2) + W        )   *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (6.2.2002)                                         *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !     ev2 = 1/2 [W1,O1]                                                *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev2     out  even2 output matrix = final result                  *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     result  intermediate result of even2-calculation                 *
    !     v       symmetric (n x n)-matrix containing (A V A)              *
    !     pvp     symmetric (n x n)-matrix containing (A P V P A)          *
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         *
    !     w1o1    W1*O1 (2-component form)                                 *
    !     o1w1    O1*W1 (2-component form)                                 *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev2
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: o1w1, pvp, pvph, v, vh, w1o1

!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH2
!-----------------------------------------------------------------------

    ALLOCATE(v(n,n))
    ALLOCATE(pVp(n,n))
    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    v=0.0_dp
    pVp=0.0_dp
    vh=0.0_dp
    pVph=0.0_dp
    v(1:n,1:n)=vv(1:n,1:n)
    vh(1:n,1:n)=vv(1:n,1:n)
    pvp(1:n,1:n)=gg(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev2=0.0_dp
    !  Calculate  v = A V A:

    CALL mat_axa_a(v,n,aa)

    !  Calculate  pvp = A P V P A:

    CALL mat_arxra_a(pvp,n,aa,rr)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h_a(vh,n,e)
    CALL mat_axa_a(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h_a(pvph,n,e)
    CALL mat_arxra_a(pvph,n,aa,rr)

    !  Create/Initialize necessary matrices:
    ALLOCATE(w1o1(n,n))
    ALLOCATE(o1w1(n,n))
    w1o1=0.0_dp
    o1w1=0.0_dp

    !  Calculate w1o1:
    CALL dgemm("N","N",n,n,n,-1.0_dp,pvph,n,v,n,0.0_dp,w1o1,n)
    CALL mat_muld_a(w1o1,pvph,pvp,n,  1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1o1,vh,  v,n,    1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,vh,n,pvp,n,1.0_dp,w1o1,n)
    !  Calculate o1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvp,n,vh,n,0.0_dp,o1w1,n)
    CALL mat_muld_a(o1w1,pvp,pvph,n,  -1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(o1w1,v,  vh,n,    -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,v,n,pvph,n,1.0_dp,o1w1,n)
    !  Calculate in symmetric pakets

    !-----------------------------------------------------------------------
    !     2.   1/2 [W1,O1] = 1/2 W1O1 -  1/2 O1W1
    !-----------------------------------------------------------------------

    CALL mat_add (ev2,0.5_dp,w1o1,-0.5_dp,o1w1,n)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!
    !-----------------------------------------------------------------------

    DEALLOCATE(v,vh,pvp,pvph,w1o1,o1w1)

!    WRITE (*,*) "CAW:  DKH2 with even2c (Alex)"
!    WRITE (*,*) "!JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even2c_a

  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param ev3 ...
!> \param e1 ...
!> \param pe1p ...
!> \param vv ...
!> \param gg ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
!> \param e ...
! *****************************************************************************
  SUBROUTINE even3b_a (n,ev3,e1,pe1p,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf, last modified:  20.2.2002 - DKH3                 *
    !                                                                      *
    !     3rd order DK-approximation (generalised DK-transformation)       *
    !                                                                      *
    !     Version: 1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (7.2.2002)                                         *
    !                                                                      *
    !     ev3 = 1/2 [W1,[W1,E1]]                                           *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev3     out  even3 output matrix = final result                  *
    !     e1      in   E1 = even1-operator                                 *
    !     pe1p    in   pE1p                                                *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     result  intermediate result of even2-calculation                 *
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         *
    !     e1      E1                                                       *
    !     pe1p    pE1p                                                     *
    !     w1w1    (W1)^2                                                   *
    !     w1e1w1  W1*E1*W1                                                 *
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,2)               *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev3
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e1, pe1p, vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pvph, scr_1, scr_2, vh, &
                                                w1e1w1, w1w1

!-----------------------------------------------------------------------
!     1.   General Structures and Patterns for DKH3
!-----------------------------------------------------------------------

    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    vh=0.0_dp
    pVph=0.0_dp
    vh(1:n,1:n)=vv(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev3=0.0_dp

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h_a(vh,n,e)
    CALL mat_axa_a(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h_a(pvph,n,e)
    CALL mat_arxra_a(pvph,n,aa,rr)

    !  Create/Initialize necessary matrices:
    ALLOCATE(w1w1(n,n))
    ALLOCATE(w1e1w1(n,n))
    ALLOCATE(scr_1(n,n))
    ALLOCATE(scr_2(n,n))
    w1w1=0.0_dp
    w1e1w1=0.0_dp
    scr_1=0.0_dp
    scr_2=0.0_dp

    !  Calculate w1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,vh,n,0.0_dp,w1w1,n)
    CALL mat_muld_a(w1w1,pvph,pvph,n,-1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1w1,vh,  vh,n,  -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pvph,n,1.0_dp,w1w1,n)

    !  Calculate w1e1w1:
    CALL mat_muld_a(scr_1 ,pvph ,pe1p,n, 1.0_dp,0.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,scr_1,n,vh,n,0.0_dp,w1e1w1,n)
    CALL mat_muld_a(w1e1w1,scr_1,pvph,n, -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,scr_2,n,vh,n,1.0_dp,w1e1w1,n)
    CALL mat_muld_a(w1e1w1,scr_2,pvph,n, 1.0_dp,1.0_dp,tt,rr)


    !-----------------------------------------------------------------------
    !     2.   ev3 = 1/2 (W1^2)E1 + 1/2 E1(W1^2) - W1E1W1
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,0.5_dp,w1w1,n,e1,n,0.0_dp,ev3,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,e1,n,w1w1,n,1.0_dp,ev3,n)
    CALL mat_add2 (ev3,1.0_dp,-1.0_dp,w1e1w1,n)

    !-----------------------------------------------------------------------
    !     3.   Finish up the stuff!!
    !-----------------------------------------------------------------------

    DEALLOCATE(vh,pvph,w1w1,w1e1w1,scr_1,scr_2)

!    WRITE (*,*) "CAW:  DKH3 with even3b (Alex)"
!    WRITE (*,*) "JT:  Now available in cp2k"
    RETURN
  END SUBROUTINE even3b_a

  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param ev4 ...
!> \param e1 ...
!> \param pe1p ...
!> \param vv ...
!> \param gg ...
!> \param aa ...
!> \param rr ...
!> \param tt ...
!> \param e ...
! *****************************************************************************
  SUBROUTINE even4a_a (n,ev4,e1,pe1p,vv,gg,aa,rr,tt,e)

    !***********************************************************************
    !                                                                      *
    !     Alexander Wolf,   last modified: 25.02.2002   --   DKH4          *
    !                                                                      *
    !     4th order DK-approximation (scalar = spin-free)                  *
    !                                                                      *
    !     Version: 1.2  (25.2.2002) :  Elegant (short) way of calculation  *
    !                                  included                            *
    !              1.1  (20.2.2002) :  Usage of SR mat_add included        *
    !              1.0  (8.2.2002)                                         *
    !                                                                      *
    !     ev4  =  1/2 [W2,[W1,E1]] + 1/8 [W1,[W1,[W1,O1]]]  =              *
    !                                                                      *
    !          =      sum_1        +         sum_2                         *
    !                                                                      *
    !                                                                      *
    !     Modification history:                                            *
    !     30.09.2006 Jens Thar: deleted obsolete F77 memory manager        *
    !                                                                      *
    !         ----  Meaning of Parameters  ----                            *
    !                                                                      *
    !     n       in   Dimension of matrices                               *
    !     ev4     out  even4 output matrix = final result                  *
    !     e1     in   E1                                                   *
    !     pe1p   in   p(E1)p                                               *
    !     vv      in   potential v                                         *
    !     gg      in   pvp                                                 *
    !     aa      in   A-Factors (DIAGONAL)                                *
    !     rr      in   R-Factors (DIAGONAL)                                *
    !     tt      in   Nonrel. kinetic Energy (DIAGONAL)                   *
    !     e       in   Rel. Energy = SQRT(p^2*c^2 + c^4)  (DIAGONAL)       *
    !     v       symmetric (n x n)-matrix containing (A V A)              *
    !     pvp     symmetric (n x n)-matrix containing (A P V P A)          *
    !     vh      symmetric (n x n)-matrix containing (A V~ A)             *
    !     pvph    symmetric (n x n)-matrix containing (A P V~ P A)         *
    !     w1w1    (W1)^2                                                   *
    !     w1o1    W1*O1      (2-component formulation)                     *
    !     o1w1    O1*W1      (2-component formulation)                     *
    !     e1      symmetric (n x n)-matrix containing E1                   *
    !     pe1p    symmetric (n x n)-matrix containing p(E1)p               *
    !     sum_i   2 addends defined above  (i=1,2)                         *
    !     scr_i   temporary (n x n)-scratch-matrices (i=1,..,4)            *
    !     scrh_i  temp. (n x n)-scr.-mat. with energy-denom. (i=1,..,4)    *
    !                                                                      *
    !***********************************************************************


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: ev4
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: e1, pe1p, vv, gg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: aa, rr, tt, e

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: o1w1, pvp, pvph, scr_1, &
      scr_2, scr_3, scr_4, scrh_1, scrh_2, scrh_3, scrh_4, sum_1, sum_2, v, &
      vh, w1o1, w1w1

!C-----------------------------------------------------------------------
!C     1.   General Structures and Patterns for DKH4
!C-----------------------------------------------------------------------

    ALLOCATE(v(n,n))
    ALLOCATE(pVp(n,n))
    ALLOCATE(vh(n,n))
    ALLOCATE(pVph(n,n))
    v=0.0_dp
    pVp=0.0_dp
    vh=0.0_dp
    pVph=0.0_dp
    v(1:n,1:n)=vv(1:n,1:n)
    vh(1:n,1:n)=vv(1:n,1:n)
    pvp(1:n,1:n)=gg(1:n,1:n)
    pvph(1:n,1:n)=gg(1:n,1:n)

    ev4=0.0_dp
    !  Calculate  v = A V A:

    CALL mat_axa_a(v,n,aa)

    !  Calculate  pvp = A P V P A:

    CALL mat_arxra_a(pvp,n,aa,rr)

    !  Calculate  vh = A V~ A:

    CALL mat_1_over_h_a(vh,n,e)
    CALL mat_axa_a(vh,n,aa)

    !  Calculate  pvph = A P V~ P A:

    CALL mat_1_over_h_a(pvph,n,e)
    CALL mat_arxra_a(pvph,n,aa,rr)


    !  Create/Initialize necessary matrices:
    ALLOCATE(w1w1(n,n))
    w1w1 = 0.0_dp
    ALLOCATE(w1o1(n,n))
    w1o1 = 0.0_dp
    ALLOCATE(o1w1(n,n))
    o1w1 = 0.0_dp
    ALLOCATE(sum_1(n,n))
    sum_1 = 0.0_dp
    ALLOCATE(sum_2(n,n))
    sum_2 = 0.0_dp
    ALLOCATE(scr_1(n,n))
    scr_1 = 0.0_dp
    ALLOCATE(scr_2(n,n))
    scr_2 = 0.0_dp
    ALLOCATE(scr_3(n,n))
    scr_3 = 0.0_dp
    ALLOCATE(scr_4(n,n))
    scr_4 = 0.0_dp
    ALLOCATE(scrh_1(n,n))
    scrh_1 = 0.0_dp
    ALLOCATE(scrh_2(n,n))
    scrh_2 = 0.0_dp
    ALLOCATE(scrh_3(n,n))
    scrh_3 = 0.0_dp
    ALLOCATE(scrh_4(n,n))
    scrh_4 = 0.0_dp

    !  Calculate w1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,vh,n,0.0_dp,w1w1,n)
    CALL mat_muld_a(w1w1,pvph,pvph,n, -1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1w1,vh,  vh,n,   -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pvph,n,1.0_dp,w1w1,n)

    !  Calculate w1o1:
    CALL dgemm("N","N",n,n,n,-1.0_dp,pvph,n,v,n,0.0_dp,w1o1,n)
    CALL mat_muld_a(w1o1,pvph,pvp,n,  1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(w1o1,vh,  v,n,    1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-1.0_dp,vh,n,pvp,n,1.0_dp,w1o1,n)
    !  Calculate o1w1:
    CALL dgemm("N","N",n,n,n,1.0_dp,pvp,n,vh,n,0.0_dp,o1w1,n)
    CALL mat_muld_a(o1w1,pvp,pvph,n,  -1.0_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(o1w1,v,  vh,n,    -1.0_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,v,n,pvph,n,1.0_dp,o1w1,n)

    !-----------------------------------------------------------------------
    !   2. sum_1 = 1/2 [W2,[W1,E1]] = 1/2 (W2W1E1 - W2E1W1 - W1E1W2 + E1W1W2)
    !-----------------------------------------------------------------------

    !  scr_i & scrh_i  for steps 2a (W2W1E1)  and 2b (W2E1W1):

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scr_1,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    CALL mat_muld_a(scr_4, pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)

    CALL mat_muld_a(scrh_1,pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    CALL mat_1_over_h_a(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scrh_2,n)
    CALL mat_1_over_h_a(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scrh_3,n)
    CALL mat_1_over_h_a(scrh_3,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scrh_4,n)
    CALL mat_1_over_h_a(scrh_4,n,e)

    !  2a)  sum_1 = 1/2 W2W1E1               ( [1]-[8] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_1,n,scr_1,n,0.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scrh_1,scr_2,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_2,n,scr_1,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scrh_2,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_3,n,scr_1,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scrh_3,scr_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(sum_1,scrh_4,scr_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_2,n,1.0_dp,sum_1,n)


    !  2b)  sum_1 = - 1/2 W2E1W1 (+ sum_1)   ( [9]-[16] )

    CALL mat_muld_a(sum_1,scrh_1,scr_3,n,-0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_1,scr_4,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_2,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_2,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_3,scr_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scrh_3,scr_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scrh_4,n,scr_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scrh_4,n,scr_4,n,1.0_dp,sum_1,n)


    !  scr_i & scrh_i  for steps 2c (W1E1W2)  and 2d (E1W1W2):

    CALL mat_muld_a(scr_1, pvph,pe1p,n,1.0_dp,0.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,pe1p,n,0.0_dp,scr_2,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,pvph,n,0.0_dp,scr_3,n)
    CALL dgemm("N","N",n,n,n,1.0_dp,e1,n,vh,n,0.0_dp,scr_4,n)

    CALL dgemm("N","N",n,n,n,1.0_dp,vh,n,e1,n,0.0_dp,scrh_1,n)
    CALL mat_1_over_h_a(scrh_1,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pvph,n,e1,n,0.0_dp,scrh_2,n)
    CALL mat_1_over_h_a(scrh_2,n,e)
    CALL dgemm("N","N",n,n,n,1.0_dp,pe1p,n,vh,n,0.0_dp,scr_3,n)
    CALL mat_1_over_h_a(scrh_3,n,e)
    CALL mat_muld_a(scrh_4,pe1p,pvph,n,1.0_dp,0.0_dp,tt,rr)
    CALL mat_1_over_h_a(scrh_4,n,e)

    !  2c)  sum_1 = - 1/2 W1E1W2 (+ sum_1)   ( [17]-[24] )

    CALL dgemm("N","N",n,n,n,0.5_dp,scr_1,n,scrh_1,n,0.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_1,scrh_2,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_2,n,scrh_1,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_2,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_1,scrh_3,n,-0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_1,scrh_4,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_2,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_2,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)

    !  2d)  sum_1 = 1/2 E1W1W2 (+ sum_1)     ( [25]-[32] )

    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_3,n,scrh_1,n,0.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_3,scrh_2,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_mulm_a(sum_1,scr_4,scrh_1,n, 0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_2,n,1.0_dp,sum_1,n)
    CALL mat_muld_a(sum_1,scr_3,scrh_3,n, 0.5_dp,1.0_dp,tt,rr)
    CALL mat_muld_a(sum_1,scr_3,scrh_4,n,-0.5_dp,1.0_dp,tt,rr)
    CALL dgemm("N","N",n,n,n,-0.5_dp,scr_4,n,scrh_3,n,1.0_dp,sum_1,n)
    CALL dgemm("N","N",n,n,n,0.5_dp,scr_4,n,scrh_4,n,1.0_dp,sum_1,n)


    !-----------------------------------------------------------------------
    !   3.  sum_2 = 1/8 [W1,[W1,[W1,O1]]] =
    !
    !             = 1/8 ( (W1^3)O1 - 3(W1^2)O1W1 + 3 W1O1(W1^2) - O1(W1^3) )
    !-----------------------------------------------------------------------

    CALL dgemm("N","N",n,n,n,0.125_dp,w1w1,n,w1o1,n,0.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.375_dp,w1w1,n,o1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,0.375_dp,w1o1,n,w1w1,n,1.0_dp,sum_2,n)
    CALL dgemm("N","N",n,n,n,-0.125_dp,o1w1,n,w1w1,n,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   4.  result = sum_1 + sum_2
    !-----------------------------------------------------------------------

    CALL mat_add(ev4,1.0_dp,sum_1,1.0_dp,sum_2,n)

    !-----------------------------------------------------------------------
    !   5. Finish up the stuff!!
    !-----------------------------------------------------------------------

    DEALLOCATE(v,pvp,vh,pvph,w1w1,w1o1,o1w1,sum_1,sum_2)
    DEALLOCATE(scr_1,scr_2,scr_3,scr_4,scrh_1,scrh_2,scrh_3,scrh_4)

!    WRITE (*,*) "CAW:  DKH4 with even4a (Alex)"
!    WRITE (*,*) "JT:   Now available in cp2k"

    RETURN
  END SUBROUTINE even4a_a

  !-----------------------------------------------------------------------
  !                                                                      -
  !     Matrix routines for DKH-procedure                                -
  !     Alexander Wolf                                                   -
  !     modifed: Jens Thar: Mem manager deleted                          -
  !     This file contains the                                           -
  !      following subroutines:                                          -
  !                                 1. mat_1_over_h                      -
  !                                 2. mat_axa                           -
  !                                 3. mat_arxra                         -
  !                                 4. mat_mulm                          -
  !                                 5. mat_muld                          -
  !                                 6. mat_add                           -
  !                                                                      -
  !-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param p ...
!> \param n ...
!> \param e ...
! *****************************************************************************
  SUBROUTINE mat_1_over_h_a (p,n,e)

    !***********************************************************************
    !                                                                      *
    !   2. SR mat_1_over_h: Transform matrix p into matrix p/(e(i)+e(j))   *
    !                                                                      *
    !   p    in  REAL(:,:) :   matrix p                                    *
    !   e    in  REAL(:)   :   rel. energy (diagonal)                      *
    !   n    in  INTEGER                                                   *
    !                                                                      *
    !***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e

    INTEGER                                  :: i, j

    DO i=1,n
       DO j=1,n
          p(i,j)=p(i,j)/(e(i)+e(j))
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_1_over_h_a
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param p ...
!> \param n ...
!> \param a ...
! *****************************************************************************
  SUBROUTINE mat_axa_a (p,n,a)

    !C***********************************************************************
    !C                                                                      *
    !C   3. SR mat_axa: Transform matrix p into matrix  a*p*a               *
    !C                                                                      *
    !C   p    in  REAL(:,:):   matrix p                                     *
    !C   a    in  REAL(:)  :   A-factors (diagonal)                         *
    !CJT n    in  INTEGER  :   dimension of matrix p                        *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a

    INTEGER                                  :: i, j

    DO i=1,n
       DO j=1,n
          p(i,j)=p(i,j)*a(i)*a(j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_axa_a
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param p ...
!> \param n ...
!> \param a ...
!> \param r ...
! *****************************************************************************
  SUBROUTINE mat_arxra_a (p,n,a,r)

    !C***********************************************************************
    !C                                                                      *
    !C   4. SR mat_arxra: Transform matrix p into matrix  a*r*p*r*a         *
    !C                                                                      *
    !C   p    in  REAL(:,:) :   matrix p                                    *
    !C   a    in  REAL(:)   :   A-factors (diagonal)                        *
    !C   r    in  REAL(:)   :   R-factors (diagonal)                        *
    !C   n    in  INTEGER   :   dimension of matrix p                       *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, r

    INTEGER                                  :: i, j

    DO i=1,n
       DO j=1,n
          p(i,j)=p(i,j)*a(i)*a(j)*r(i)*r(j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_arxra_a
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param p ...
!> \param q ...
!> \param r ...
!> \param n ...
!> \param alpha ...
!> \param beta ...
!> \param t ...
!> \param rr ...
! *****************************************************************************
  SUBROUTINE mat_mulm_a (p,q,r,n,alpha,beta,t,rr)

    !C***********************************************************************
    !C                                                                      *
    !C   5. SR mat_mulm:  Multiply matrices according to:                   *
    !C                                                                      *
    !C                      p = alpha*q*(..P^2..)*r + beta*p                *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,.):   matrix r                                  *
    !C   n      in   INTEGER  :   dimension n of matrices                   *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: p, q, r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr

    INTEGER                                  :: i, j
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: qtemp

    ALLOCATE(qtemp(n,n))

    DO i=1,n
       DO j=1,n
          qtemp(i,j)=q(i,j)*2.0_dp*t(j)*rr(j)*rr(j)
       ENDDO
    ENDDO

    CALL dgemm("N","N",n,n,n,alpha,qtemp,n,r,n,beta,p,n)
    RETURN
    DEALLOCATE(qtemp)
  END SUBROUTINE mat_mulm_a
  !C-----------------------------------------------------------------------



! *****************************************************************************
!> \brief ...
!> \param p ...
!> \param q ...
!> \param r ...
!> \param n ...
!> \param alpha ...
!> \param beta ...
!> \param t ...
!> \param rr ...
! *****************************************************************************
  SUBROUTINE mat_muld_a (p,q,r,n,alpha,beta,t,rr)

    !C***********************************************************************
    !C                                                                      *
    !C   16. SR mat_muld:  Multiply matrices according to:                  *
    !C                                                                      *
    !C                      p = alpha*q*(..1/P^2..)*r + beta*p              *
    !C                                                                      *
    !C   p      out  REAL(:,:):   matrix p                                  *
    !C   q      in   REAL(:,:):   matrix q                                  *
    !C   r      in   REAL(:,:):   matrix r                                  *
    !C   n      in   INTEGER  :   Dimension of all matrices                 *
    !C   alpha  in   REAL(dp) :                                             *
    !C   beta   in   REAL(dp) :                                             *
    !C   t      in   REAL(:)  :   non-rel. kinetic energy  (diagonal)       *
    !C   rr     in   REAL(:)  :   R-factors  (diagonal)                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: p, q, r
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: t, rr

    INTEGER                                  :: i, j
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: qtemp

    ALLOCATE(qtemp(n,n))

    DO i=1,n
       DO j=1,n
          qtemp(i,j)=q(i,j)*0.5_dp/(t(j)*rr(j)*rr(j))
       ENDDO
    ENDDO

    CALL dgemm("N","N",n,n,n,alpha,qtemp,n,r,n,beta,p,n)

    RETURN
    DEALLOCATE(qtemp)
  END SUBROUTINE mat_muld_a
  !C-----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param p ...
!> \param alpha ...
!> \param beta ...
!> \param r ...
!> \param n ...
! *****************************************************************************
  SUBROUTINE mat_add2 (p,alpha,beta,r,n)

    !C***********************************************************************
    !C                                                                      *
    !C   19. SR mat_add:  Add two matrices of the same size according to:   *
    !C                                                                      *
    !C                            p = alpha*p + beta*r                      *
    !C                                                                      *
    !C                    and store them in the first                       *
    !C   p      out  REAL(:,:)  :   matrix p                                *
    !C   r      in   REAL(:,:)  :   matrix r                                *
    !C   alpha  in   REAL(dp)                                               *
    !C   beta   in   REAL(dp)                                               *
    !C                                                                      *
    !C   Matrix p must already exist before calling this SR!!               *
    !C                                                                      *
    !C  [written by: Alexander Wolf,  20.2.2002,  v1.0]                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: p
    REAL(KIND=dp), INTENT(IN)                :: alpha, beta
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, j

!C  Add matrices:

    DO i=1,n
       DO j=1,n
          p(i,j) = alpha*p(i,j) + beta*r(i,j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_add2

  !CAW---------------------------------------------------------------------


! *****************************************************************************
!> \brief ...
!> \param p ...
!> \param alpha ...
!> \param q ...
!> \param beta ...
!> \param r ...
!> \param n ...
! *****************************************************************************
  SUBROUTINE mat_add (p,alpha,q,beta,r,n)

    !C***********************************************************************
    !C                                                                      *
    !C   19. SR mat_add:  Add two matrices of the same size according to:   *
    !C                                                                      *
    !C                            p = alpha*q + beta*r                      *
    !C                                                                      *
    !C   p      out  REAL(:,:)  :   matrix p                                *
    !C   q      in   REAL(:,:)  :   matrix q                                *
    !C   r      in   REAL(:,:)  :   matrix r                                *
    !C   alpha  in   REAL(dp)                                               *
    !C   beta   in   REAL(dp)                                               *
    !C                                                                      *
    !C   Matrix p must already exist before calling this SR!!               *
    !C                                                                      *
    !C  [written by: Alexander Wolf,  20.2.2002,  v1.0]                     *
    !C                                                                      *
    !C***********************************************************************


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: p
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: q
    REAL(KIND=dp), INTENT(IN)                :: beta
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: r
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, j

!C  Add matrices:

    DO i=1,n
       DO j=1,n
          p(i,j) = alpha*q(i,j) + beta*r(i,j)
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE mat_add

  !CAW---------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param W ...
!> \param B ...
!> \param C ...
!> \param N ...
!> \param H ...
! *****************************************************************************
  SUBROUTINE TRSM ( W,B,C,N,H)


    REAL(KIND=dp), DIMENSION(:, :)           :: W, B, C
    INTEGER                                  :: N
    REAL(KIND=dp), DIMENSION(:, :)           :: H

    INTEGER                                  :: I, IJ, J, K, L

!C
!C     TRANSFORM SYMMETRIC matrix A by UNITARY TRANSFORMATION
!C     IN B. RESULT IS IN C
!C
!CAW      C = B^{dagger} * A * B

    IJ=0
    DO I=1,N
       DO J=1,I
          IJ=IJ+1
          C(I,J)=0.0_dp
          C(J,I)=0.0_dp
          H(I,J)=0.0_dp
          H(J,I)=0.0_dp
       END DO
    END DO
    DO I=1,N
       DO L=1,N
          DO K=1,N
             H(I,L)=B(K,I)*W(K,L)+H(I,L)
          END DO
       END DO
    END DO

    IJ=0
    DO I=1,N
       DO J=1,I
          IJ=IJ+1
          DO L=1,N
             C(I,J)=H(I,L)*B(L,J)+C(I,J)
             C(J,I)=C(I,J)
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE TRSM

  !JT----------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param matrix_t_pgf ...
!> \param n ...
!> \param eig ...
!> \param ew ...
!> \param matrix_sinv_pgf ...
!> \param aux ...
!> \param ic ...
! *****************************************************************************
  SUBROUTINE dkh_diag (matrix_t_pgf,n,eig,ew,matrix_sinv_pgf,aux,ic)


    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: matrix_t_pgf
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: eig
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: ew
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix_sinv_pgf
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: aux
    INTEGER                                  :: ic

    INTEGER                                  :: n2

    eig = 0.0_dp
    aux = 0.0_dp

    CALL dgemm("N","N",n,n,n,1.0_dp,matrix_t_pgf,n,matrix_sinv_pgf,n,0.0_dp,eig,n)

    aux = 0.0_dp

    CALL dgemm("T","N",n,n,n,1.0_dp,matrix_sinv_pgf,n,eig,n,0.0_dp,aux,n)

    n2=3*n-1

    CALL JACOB2 ( AUX,EIG,EW,N,IC )
    RETURN
  END SUBROUTINE dkh_diag

  !CAW---------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param sogt ...
!> \param eigv ...
!> \param eigw ...
!> \param n ...
!> \param ic ...
! *****************************************************************************
  SUBROUTINE JACOB2 ( sogt,eigv,eigw,n,ic )


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), DIMENSION(n), INTENT(OUT) :: eigw
    REAL(KIND=dp), DIMENSION(n, n), &
      INTENT(OUT)                            :: eigv
    REAL(KIND=dp), DIMENSION(n, n), &
      INTENT(INOUT)                          :: sogt
    INTEGER, INTENT(IN)                      :: ic

    INTEGER                                  :: i, il, im, ind, j, k, l, ll, &
                                                m, mm
    REAL(KIND=dp)                            :: cost, cost2, ext_norm, sincs, &
                                                sint, sint2, thr, thr_min, &
                                                tol, u1, x, xy, y

    tol=1.0E-15
    ext_norm=0.0_dp
    u1=REAL(n)
    DO i=1,n
       eigv(i,i)=1.0_dp
       eigw(i)=sogt(i,i)
       DO j=1,i
          IF(i.ne.j) THEN
             eigv(i,j)=0.0_dp
             eigv(j,i)=0.0_dp
             ext_norm=ext_norm+sogt(i,j)*sogt(i,j)
          END IF
       END DO
    END DO

    IF (ext_norm.gt.0.0_dp) THEN
       ext_norm=SQRT(2.0_dp*ext_norm)
       thr_min=ext_norm*tol/u1
       ind=0
       thr=ext_norm

       DO
          thr=thr/u1
          DO
             l=1
             DO
                m=l+1
                DO
                   IF ((ABS(sogt(m,l))-thr).ge.0.0_dp) THEN
                      ind=1
                      x=0.5_dp*(eigw(l)-eigw(m))
                      y=-sogt(m,l)/SQRT(sogt(m,l)*sogt(m,l)+x*x)
                      IF (x.lt.0.0_dp) y=-y

                      IF (y.gt.1.0_dp) y=1.0_dp
                      IF (y.lt.-1.0_dp) y=-1.0_dp
                      xy=1.0_dp-y*y
                      sint=y/SQRT(2.0_dp*(1.0_dp+SQRT(xy)))
                      sint2=sint*sint
                      cost2=1.0_dp-sint2
                      cost=SQRT(cost2)
                      sincs=sint*cost

                      DO i=1,n
                         IF((i-m).ne.0) THEN
                            IF ((i-m).lt.0) THEN
                               im=m
                               mm=i
                            ELSE
                               im=i
                               mm=m
                            END IF
                            IF ((i-l).ne.0) THEN
                               IF ((i-l).lt.0) THEN
                                  il=l
                                  ll=i
                               ELSE
                                  il=i
                                  ll=l
                               END IF
                               x=sogt(il,ll)*cost-sogt(im,mm)*sint
                               sogt(im,mm)=sogt(il,ll)*sint+sogt(im,mm)*cost
                               sogt(il,ll)=x
                            END IF
                         END IF

                         x=eigv(i,l)*cost-eigv(i,m)*sint
                         eigv(i,m)=eigv(i,l)*sint+eigv(i,m)*cost
                         eigv(i,l)=x
                      END DO

                      x=2.0_dp*sogt(m,l)*sincs
                      y=eigw(l)*cost2+eigw(m)*sint2-x
                      x=eigw(l)*sint2+eigw(m)*cost2+x
                      sogt(m,l)=(eigw(l)-eigw(m))*sincs+sogt(m,l)*(cost2-sint2)
                      eigw(l)=y
                      eigw(m)=x
                   END IF
                   IF ((m-n).eq.0) EXIT
                   m=m+1
                END DO
                IF ((l-m+1).eq.0) EXIT
                l=l+1
             END DO
             IF((ind-1).ne.0.0_dp) EXIT
             ind=0
          END DO
          IF ((thr-thr_min).le.0.0_dp) EXIT
       END DO
    END IF

    IF (ic.ne.0) THEN
       DO i=1,n
          DO j=1,n
             IF ((eigw(i)-eigw(j)).gt.0.0_dp) THEN
                x=eigw(i)
                eigw(i)=eigw(j)
                eigw(j)=x
                DO k=1,n
                   y=eigv(k,i)
                   eigv(k,i)=eigv(k,j)
                   eigv(k,j)=y
                END DO
             END IF
          END DO
       END DO

    END IF
    RETURN
  END SUBROUTINE JACOB2

  !JT---------------------------------------------------------------------

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \param matrix_s_pgf ...
!> \param matrix_sinv_pgf ...
! *****************************************************************************
  SUBROUTINE SOG (n,matrix_s_pgf,matrix_sinv_pgf)


    INTEGER                                  :: n
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: matrix_s_pgf
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: matrix_sinv_pgf

    INTEGER                                  :: i, j, jn, k
    REAL(KIND=dp)                            :: diag_s, row_sum, scalar
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: a
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: g

!
!     SUBROUTINE TO CALCULATE TRANSFORMATION TO SCHMIDT-
!     ORTHOGONALIZED BASIS.
!     sinv-1*matrix_s_pgf*sinv = "orthogonal matrix"
!     n              dimension of matrices
!     matrix_s_pgf   original overlap matrix
!     matrix_sinv_pgf new overlap matrix
!     g              scratch
!     a              scratch
!

    ALLOCATE(a(n))
    ALLOCATE(g(n,n))

    DO jn=1,n
       diag_s = matrix_s_pgf(jn,jn)
       g(jn,jn)=1.0_dp

       IF(jn.ne.1) THEN
          DO j=1,jn-1
             scalar=0.0_dp
             DO i=1,j
                scalar=scalar+matrix_s_pgf(i,jn)*g(i,j)
             END DO
             diag_s=diag_s-scalar*scalar
             a(j) = scalar
          END DO

          DO j=1,jn-1
             row_sum=0.0_dp
             DO k=j,jn-1
                row_sum=row_sum+a(k)*g(j,k)
             END DO
             g(j,jn)=-row_sum
          END DO
       END IF

       diag_s=1.0_dp/SQRT(diag_s)
       DO i=1,jn
          g(i,jn)=g(i,jn)*diag_s
       END DO
    END DO

    DO j=1,n
       DO i=1,j
          matrix_sinv_pgf(j,i)=0.0_dp
          matrix_sinv_pgf(i,j)=g(i,j)
       END DO
    END DO
    DEALLOCATE(a,g)

    RETURN

  END SUBROUTINE SOG

END MODULE dkh_main
