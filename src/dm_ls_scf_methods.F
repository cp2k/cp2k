!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief lower level routines for linear scaling SCF
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf_methods
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_filter, cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_get_occupation, cp_dbcsr_init, cp_dbcsr_multiply, &
       cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_trace
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_dbcsr_util,                   ONLY: lanczos_alg_serial
  USE dbcsr_api,                       ONLY: &
       dbcsr_data_clear_pointer, dbcsr_data_init, dbcsr_data_new, &
       dbcsr_data_obj, dbcsr_data_release, dbcsr_get_data_type, &
       dbcsr_iterator, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
       dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_put_block, &
       dbcsr_type_no_symmetry
  USE dm_ls_scf_types,                 ONLY: ls_cluster_atomic,&
                                             ls_mstruct_type
  USE f77_blas
  USE input_constants,                 ONLY: ls_s_preconditioner_atomic,&
                                             ls_s_preconditioner_molecular,&
                                             ls_s_preconditioner_none
  USE iterate_matrix,                  ONLY: matrix_sign_Newton_Schulz,&
                                             matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_methods'

  PUBLIC density_matrix_sign, density_matrix_sign_fixed_mu
  PUBLIC compute_matrix_preconditioner, apply_matrix_preconditioner
  PUBLIC density_matrix_trs4, compute_homo_lumo

CONTAINS

! *****************************************************************************
!> \brief compute for a block positive definite matrix s (bs)
!>        the sqrt(bs) and inv(sqrt(bs))
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE compute_matrix_preconditioner(matrix_s,preconditioner_type, ls_mstruct, &
             matrix_bs_sqrt,matrix_bs_sqrt_inv,threshold,order,eps_lanczos, max_iter_lanczos, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_s
    INTEGER                                  :: preconditioner_type
    TYPE(ls_mstruct_type)                    :: ls_mstruct
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_bs_sqrt, &
                                                matrix_bs_sqrt_inv
    REAL(KIND=dp)                            :: threshold
    INTEGER                                  :: order
    REAL(KIND=dp)                            :: eps_lanczos
    INTEGER                                  :: max_iter_lanczos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'compute_matrix_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, iblock_row
    LOGICAL                                  :: block_needed, tr
    TYPE(cp_dbcsr_type)                      :: matrix_bs
    TYPE(dbcsr_data_obj)                     :: block_any
    TYPE(dbcsr_iterator)                     :: iter

    CALL timeset(routineN,handle)

    ! first generate a block diagonal copy of s
    CALL cp_dbcsr_init(matrix_bs,error=error)
    CALL cp_dbcsr_create(matrix_bs,template=matrix_s,error=error)

    CALL dbcsr_data_init (block_any)
    CALL dbcsr_data_new(block_any, dbcsr_get_data_type(matrix_s%matrix))

    SELECT CASE(preconditioner_type)
    CASE(ls_s_preconditioner_none)
    CASE(ls_s_preconditioner_atomic, ls_s_preconditioner_molecular)
       CALL dbcsr_iterator_start(iter, matrix_s%matrix)
       DO WHILE (dbcsr_iterator_blocks_left(iter))
          !CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, block_s, blk)
          CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, block_any, tr)

          ! do we need the block ?
          ! this depends on the preconditioner, but also the matrix clustering method employed
          ! for a clustered matrix, right now, we assume that atomic and molecular preconditioners
          ! are actually the same, and only require that the diagonal blocks (clustered) are present

          block_needed=.FALSE.

          IF (iblock_row==iblock_col) THEN
              block_needed=.TRUE.
          ELSE
              IF (preconditioner_type==ls_s_preconditioner_molecular .AND. &
                  ls_mstruct%cluster_type==ls_cluster_atomic) THEN
                 IF (ls_mstruct%atom_to_molecule(iblock_row)==ls_mstruct%atom_to_molecule(iblock_col)) block_needed=.TRUE.
              ENDIF
          ENDIF

          ! add it
          IF (block_needed) THEN
             !CALL cp_dbcsr_put_block(matrix=matrix_bs,&
             !                   row=iblock_row, col=iblock_col, block=block_s)
             CALL dbcsr_put_block(matrix=matrix_bs%matrix,&
                                row=iblock_row, col=iblock_col, block=block_any)

          ENDIF

       ENDDO
       CALL dbcsr_iterator_stop(iter)
    END SELECT

    CALL dbcsr_data_clear_pointer (block_any)
    CALL dbcsr_data_release (block_any)

    CALL cp_dbcsr_finalize(matrix_bs,error=error)

    SELECT CASE(preconditioner_type)
    CASE(ls_s_preconditioner_none)
       ! for now make it a simple identity matrix
       CALL cp_dbcsr_copy(matrix_bs_sqrt,matrix_bs,error=error)
       CALL cp_dbcsr_set(matrix_bs_sqrt,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(matrix_bs_sqrt,1.0_dp,error=error)

       ! for now make it a simple identity matrix
       CALL cp_dbcsr_copy(matrix_bs_sqrt_inv,matrix_bs,error=error)
       CALL cp_dbcsr_set(matrix_bs_sqrt_inv,0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(matrix_bs_sqrt_inv,1.0_dp,error=error)
    CASE(ls_s_preconditioner_atomic, ls_s_preconditioner_molecular)
       CALL cp_dbcsr_copy(matrix_bs_sqrt,matrix_bs,error=error)
       CALL cp_dbcsr_copy(matrix_bs_sqrt_inv,matrix_bs,error=error)
       ! XXXXXXXXXXX
       ! XXXXXXXXXXX the threshold here could be done differently,
       ! XXXXXXXXXXX using eps_filter is reducing accuracy for no good reason, this is cheap
       ! XXXXXXXXXXX
       CALL matrix_sqrt_Newton_Schulz(matrix_bs_sqrt,matrix_bs_sqrt_inv,matrix_bs,&
                                      threshold=MIN(threshold,1.0E-10_dp),order=order, &
                                      eps_lanczos=eps_lanczos, max_iter_lanczos=max_iter_lanczos, error=error)
    END SELECT

    CALL cp_dbcsr_release(matrix_bs,error=error)

    CALL timestop(handle)

  END SUBROUTINE compute_matrix_preconditioner

! *****************************************************************************
!> \brief apply a preconditioner either
!>        forward (precondition)            inv(sqrt(bs)) * A * inv(sqrt(bs))
!>        backward (restore to old form)        sqrt(bs)  * A * sqrt(bs)
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE apply_matrix_preconditioner(matrix,direction,matrix_bs_sqrt,matrix_bs_sqrt_inv,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix
    CHARACTER(LEN=*)                         :: direction
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_bs_sqrt, &
                                                matrix_bs_sqrt_inv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_matrix_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix,matrix_type=dbcsr_type_no_symmetry,error=error)

    SELECT CASE(direction)
    CASE("forward")
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_bs_sqrt_inv,&
                                0.0_dp, matrix_tmp, error=error)
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_bs_sqrt_inv, matrix_tmp,&
                                0.0_dp, matrix, error=error)
    CASE("backward")
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_bs_sqrt,&
                                0.0_dp, matrix_tmp, error=error)
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_bs_sqrt, matrix_tmp,&
                                0.0_dp, matrix, error=error)
    CASE DEFAULT
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL cp_dbcsr_release(matrix_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_matrix_preconditioner

! *****************************************************************************
!> \brief compute the density matrix with a trace that is close to nelectron.
!>        take a mu as input, and improve by bisection as needed.
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE density_matrix_sign(matrix_p,mu,fixed_mu,matrix_ks,matrix_s,matrix_s_inv,nelectron,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p
    REAL(KIND=dp), INTENT(INOUT)             :: mu
    LOGICAL                                  :: fixed_mu
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_ks, matrix_s, &
                                                matrix_s_inv
    INTEGER, INTENT(IN)                      :: nelectron
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_sign', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: initial_increment = 0.01_dp

    INTEGER                                  :: handle, iter, unit_nr
    LOGICAL                                  :: has_mu_high, has_mu_low
    REAL(KIND=dp)                            :: increment, mu_high, mu_low, &
                                                trace
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    increment=initial_increment

    has_mu_low=.FALSE.
    has_mu_high=.FALSE.

    ! bisect if both bounds are known, otherwise find the bounds with a linear search
    DO iter=1,30
       IF (has_mu_low .AND. has_mu_high) THEN
          mu=(mu_low+mu_high)/2
          IF (ABS(mu_high-mu_low)<threshold) EXIT
       ENDIF

       CALL density_matrix_sign_fixed_mu(matrix_p,trace,mu,matrix_ks,matrix_s,matrix_s_inv,threshold,error)
       IF (unit_nr>0) WRITE(unit_nr,'(T2,A,I2,1X,F13.9,1X,F15.9)') &
             "Density matrix:  iter, mu, trace error: ", iter, mu, trace-nelectron

       ! OK, we can skip early if we are as close as possible to the exact result
       ! smaller differences should be considered 'noise'
       IF (ABS(trace-nelectron)<0.5_dp .OR. fixed_mu) EXIT

       IF (trace<nelectron) THEN
          mu_low=mu
          mu=mu+increment
          has_mu_low=.TRUE.
          increment=increment*2
       ELSE
          mu_high=mu
          mu=mu-increment
          has_mu_high=.TRUE.
          increment=increment*2
       ENDIF
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE density_matrix_sign

! *****************************************************************************
!> \brief for a fixed mu, compute the corresponding density matrix and its trace
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE density_matrix_sign_fixed_mu(matrix_p,trace,mu,matrix_ks,matrix_s,matrix_s_inv,threshold,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p
    REAL(KIND=dp), INTENT(OUT)               :: trace
    REAL(KIND=dp), INTENT(INOUT)             :: mu
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_ks, matrix_s, &
                                                matrix_s_inv
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_sign_fixed_mu', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    REAL(KIND=dp)                            :: frob_matrix
    TYPE(cp_dbcsr_type)                      :: matrix_p_ud, matrix_sign, &
                                                matrix_sinv_ks, matrix_tmp
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! get inv(S)*H-I*mu
    CALL cp_dbcsr_init(matrix_sinv_ks,error=error)
    CALL cp_dbcsr_create(matrix_sinv_ks,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s_inv, matrix_ks,&
                           0.0_dp, matrix_sinv_ks, filter_eps=threshold,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_sinv_ks,-mu,error=error)

    ! compute sign(inv(S)*H-I*mu)
    CALL cp_dbcsr_init(matrix_sign,error=error)
    CALL cp_dbcsr_create(matrix_sign,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL matrix_sign_Newton_Schulz(matrix_sign,matrix_sinv_ks,threshold,error)
    CALL cp_dbcsr_release(matrix_sinv_ks,error=error)

    ! now construct the density matrix PS=0.5*(I-sign(inv(S)H-I*mu))
    CALL cp_dbcsr_init(matrix_p_ud,error=error)
    CALL cp_dbcsr_create(matrix_p_ud,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_copy(matrix_p_ud,matrix_sign,error=error)
    CALL cp_dbcsr_scale(matrix_p_ud,-0.5_dp,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_p_ud,0.5_dp,error=error)
    CALL cp_dbcsr_release(matrix_sign,error=error)

    ! we now have PS, lets get its trace
    CALL cp_dbcsr_trace(matrix_p_ud,trace,error=error)

    ! we can also check it is idempotent PS*PS=PS
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix_s,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p_ud, matrix_p_ud,&
                           0.0_dp, matrix_tmp, filter_eps=threshold,error=error)
    CALL cp_dbcsr_add(matrix_tmp, matrix_p_ud, 1.0_dp, -1.0_dp, error=error)
    frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp)
    CALL cp_dbcsr_release(matrix_tmp,error=error)
    IF (unit_nr>0) WRITE(unit_nr,'(T2,A,F20.12)') "Deviation from idempotency: ", frob_matrix

    ! get P=PS*inv(S)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p_ud, matrix_s_inv,&
                           0.0_dp, matrix_p, filter_eps=threshold,error=error)
    CALL cp_dbcsr_release(matrix_p_ud,error=error)

    CALL timestop(handle)

  END SUBROUTINE density_matrix_sign_fixed_mu

! *****************************************************************************
!> \brief compute the density matrix using a trace-resetting algorithm
!> \par History
!>       2012.06 created [Florian Thoele]
!> \author Florian Thoele 
! *****************************************************************************
  SUBROUTINE density_matrix_trs4(matrix_p, matrix_ks, matrix_s_sqrt_inv, &
                                 nelectron, threshold, e_homo, e_lumo, e_mu, &
                                 dynamic_threshold, matrix_ks_deviation,&
                                 max_iter_lanczos, eps_lanczos, error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_p
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_ks, matrix_s_sqrt_inv
    INTEGER, INTENT(IN)                      :: nelectron
    REAL(KIND=dp), INTENT(IN)                :: threshold
    REAL(KIND=dp), INTENT(INOUT)             :: e_homo, e_lumo, e_mu
    LOGICAL, INTENT(IN), OPTIONAL            :: dynamic_threshold
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_ks_deviation
    INTEGER, INTENT(IN)                      :: max_iter_lanczos
    REAL(KIND=dp), INTENT(IN)                :: eps_lanczos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_trs4', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_iter = 100
    REAL(KIND=dp), PARAMETER                 :: gamma_max = 6.0_dp, &
                                                gamma_min = 0.0_dp

    INTEGER                                  :: branch, estimated_steps, &
                                                handle, i, j, unit_nr
    INTEGER(kind=int_8)                      :: flop1, flop2
    LOGICAL                                  :: converged, do_dyn_threshold
    REAL(KIND=dp) :: current_threshold, eps_max, eps_min, est_threshold, &
      frob_id, frob_x, gam, homo, lumo, max_eig, max_threshold, maxdev, &
      maxev, min_eig, minev, mmin, mu, mu_a, mu_b, mu_c, mu_fa, mu_fc, &
      occ_matrix, scaled_homo_bound, scaled_lumo_bound, t1, t2, trace_fx, &
      trace_gx, xi
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: gamma_values
    TYPE(cp_dbcsr_type)                      :: matrix_k0, matrix_x, &
                                                matrix_xidsq, matrix_xsq, &
                                                tmp_gx
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    do_dyn_threshold = .FALSE.
    IF(PRESENT(dynamic_threshold)) do_dyn_threshold = dynamic_threshold


    ! init X = (eps_n*I - H)/(eps_n - eps_0)  ... H* = S^-1/2*H*S^-1/2
    CALL cp_dbcsr_init(matrix_x,error=error)
    CALL cp_dbcsr_create(matrix_x, template=matrix_ks, matrix_type=dbcsr_type_no_symmetry, error=error)

    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s_sqrt_inv, matrix_ks,&
                           0.0_dp, matrix_x, filter_eps=threshold,error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_x, matrix_s_sqrt_inv, &
                           0.0_dp, matrix_x, filter_eps=threshold,error=error)

    CALL cp_dbcsr_init(matrix_k0,error=error)
    CALL cp_dbcsr_create(matrix_k0,template=matrix_ks,matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_copy(matrix_k0, matrix_x, error=error)

    ! compute the deviation in the mixed matrix, as seen in the ortho basis
    CALL cp_dbcsr_add(matrix_ks_deviation, matrix_x , -1.0_dp, 1.0_dp, error=error)
    CALL lanczos_alg_serial(matrix_ks_deviation, maxev, minev, max_iter=max_iter_lanczos, threshold=eps_lanczos, &
                            converged=converged, error=error)
    maxdev = MAX(ABS(maxev), ABS(minev))
    IF (unit_nr>0) THEN
       WRITE(unit_nr, '(T6,A,1X,L12)')   "Lanczos converged:      ", converged
       WRITE(unit_nr, '(T6,A,1X,F12.5)') "change in mixed matrix: ", maxdev
       WRITE(unit_nr, '(T6,A,1X,F12.5)') "HOMO upper bound:       ", e_homo+maxdev
       WRITE(unit_nr, '(T6,A,1X,F12.5)') "LUMO lower bound:       ", e_lumo-maxdev
       WRITE(unit_nr, '(T6,A,1X,L12)')   "Predicts a gap ?        ", ((e_lumo-maxdev)-(e_homo+maxdev))>0
    ENDIF
    ! save the old mixed matrix
    CALL cp_dbcsr_copy(matrix_ks_deviation, matrix_x , error=error)

    ! get largest/smallest eigenvalues for scaling
    CALL lanczos_alg_serial(matrix_x, max_eig, min_eig, max_iter=max_iter_lanczos, threshold=eps_lanczos,&
                            converged=converged, error=error)
    IF (unit_nr>0) WRITE(unit_nr,'(T6,A,1X,2F12.5,1X,A,1X,L1)') "Est. extremal eigenvalues", &
                                                         min_eig, max_eig," converged: ",converged
    eps_max = max_eig
    eps_min = min_eig
    
    ! scale KS matrix
    CALL cp_dbcsr_scale(matrix_x, -1.0_dp, error=error)
    CALL cp_dbcsr_add_on_diag(matrix_x, eps_max, error=error)
    CALL cp_dbcsr_scale(matrix_x, 1/(eps_max-eps_min), error=error)

    ! scale bounds for HOMO/LUMO
    scaled_homo_bound = (eps_max-(e_homo+maxdev))/(eps_max-eps_min)
    scaled_lumo_bound = (eps_max-(e_lumo-maxdev))/(eps_max-eps_min)
    current_threshold = threshold

    CALL cp_dbcsr_init(matrix_xsq,error=error)
    CALL cp_dbcsr_create(matrix_xsq,template=matrix_ks,matrix_type=dbcsr_type_no_symmetry,error=error)

    CALL cp_dbcsr_init(matrix_xidsq,error=error)
    CALL cp_dbcsr_create(matrix_xidsq,template=matrix_ks,matrix_type=dbcsr_type_no_symmetry,error=error)

    CALL cp_dbcsr_init(tmp_gx,error=error)
    CALL cp_dbcsr_create(tmp_gx,template=matrix_ks,matrix_type=dbcsr_type_no_symmetry,error=error)

    ALLOCATE(gamma_values(max_iter))

    DO i=1,max_iter
      t1 = m_walltime()
      flop1 = 0; flop2 = 0

      ! get X*X
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_x, matrix_x,&
                             0.0_dp, matrix_xsq, &
                             filter_eps=current_threshold,flop=flop1,error=error)

      ! intermediate use matrix_xidsq to compute = X*X-X
      CALL cp_dbcsr_copy(matrix_xidsq, matrix_x,error=error)
      CALL cp_dbcsr_add(matrix_xidsq, matrix_xsq, -1.0_dp, 1.0_dp, error=error)
      frob_id = cp_dbcsr_frobenius_norm(matrix_xidsq)
      frob_x = cp_dbcsr_frobenius_norm(matrix_x)

      ! xidsq = (1-X)*(1-X) 
      ! use (1-x)*(1-x) = 1 + x*x - 2*x 
      CALL cp_dbcsr_copy(matrix_xidsq, matrix_x,error=error)
      CALL cp_dbcsr_add(matrix_xidsq, matrix_xsq, -2.0_dp, 1.0_dp, error=error)
      CALL cp_dbcsr_add_on_diag(matrix_xidsq, 1.0_dp, error=error)
    
      ! tmp_gx = 4X-3X*X
      CALL cp_dbcsr_copy(tmp_gx, matrix_x,error=error)
      CALL cp_dbcsr_add(tmp_gx, matrix_xsq, 4.0_dp, -3.0_dp, error=error)

      ! get gamma 
      ! Tr(F) = Tr(XX*tmp_gx) Tr(G) is equivalent
      CALL cp_dbcsr_trace(matrix_xsq, matrix_xidsq, trace_gx,  error=error)
      CALL cp_dbcsr_trace(matrix_xsq, tmp_gx, trace_fx, error=error)
      gam = (nelectron - trace_fx)/trace_gx 
      gamma_values(i) = gam


      ! quantities used for dynamic thresholding, when the estimated gap is larger than zero
      xi = (scaled_homo_bound-scaled_lumo_bound)
      IF (do_dyn_threshold .AND. xi > 0.0_dp) THEN
        mmin = 0.5*(scaled_homo_bound+scaled_lumo_bound)
        max_threshold = ABS(1-2*mmin)*xi
  
        scaled_homo_bound = evaluate_trs4_polynomial(scaled_homo_bound, gamma_values(i:), 1)
        scaled_lumo_bound = evaluate_trs4_polynomial(scaled_lumo_bound, gamma_values(i:), 1)
        estimated_steps = estimate_steps(scaled_homo_bound, scaled_lumo_bound, threshold)
      
        est_threshold = (threshold/(estimated_steps+i+1))*xi/(1+threshold/(estimated_steps+i+1))
        est_threshold = MIN(max_threshold, est_threshold)
        IF (i > 1) est_threshold = MAX(est_threshold, 0.1_dp * current_threshold)
        current_threshold = est_threshold
      ELSE
        current_threshold = threshold
      ENDIF

      IF (gam > gamma_max) THEN
        ! Xn+1 = 2X-X*X
        CALL cp_dbcsr_add(matrix_x, matrix_xsq, 2.0_dp, -1.0_dp, error=error)
        CALL cp_dbcsr_filter(matrix_x, current_threshold, error=error)
        branch=1
      ELSE IF (gam < gamma_min) THEN
        ! Xn+1 = X*X
        CALL cp_dbcsr_copy(matrix_x, matrix_xsq, error=error)
        branch=2
      ELSE
        ! Xn+1 = F(X) + gam*G(X)
        CALL cp_dbcsr_add(tmp_gx, matrix_xidsq, 1.0_dp, gam, error=error)
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_xsq, tmp_gx,&
                               0.0_dp, matrix_x, &
                               flop=flop2, filter_eps=current_threshold,error=error)
        branch=3
      ENDIF
   
      occ_matrix = cp_dbcsr_get_occupation(matrix_x) 
      t2 = m_walltime()
      IF (unit_nr>0) THEN
         WRITE(unit_nr,&
                      '(T6,A,I3,1X,F10.8,E12.3,F12.3,F13.3,E12.3)') "TRS4 it ", &
                      i, occ_matrix, ABS(trace_gx), t2-t1,&
                      (flop1+flop2)/(1.0E6_dp*MAX(t2-t1,0.001_dp)), current_threshold
         CALL m_flush(unit_nr)
      ENDIF

      ! a branch of 1 or 2 appears to lead to a less accurate electron number count and premature exit
      ! if it turns out this does not exit because we get stuck in branch 1/2 for a reason we need to refine further
      IF (ABS(frob_id/frob_x) < SQRT(threshold) .AND. branch==3) EXIT

    END DO

    occ_matrix = cp_dbcsr_get_occupation(matrix_x) 
    IF (unit_nr>0) WRITE(unit_nr, '(T6,A,I3,1X,F10.8,E12.3)') 'Final TRS4 iteration  ', i, occ_matrix,ABS(trace_gx)

    ! free some memory
    CALL cp_dbcsr_release(tmp_gx, error=error)
    CALL cp_dbcsr_release(matrix_xsq, error=error)
    CALL cp_dbcsr_release(matrix_xidsq, error=error)

    ! output to matrix_p, P = inv(S)^0.5 X inv(S)^0.5
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_x, matrix_s_sqrt_inv,&
                           0.0_dp, matrix_p, filter_eps=threshold,error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s_sqrt_inv, matrix_p,&
                           0.0_dp, matrix_p, filter_eps=threshold,error=error)
    
    ! calculate the chemical potential by doing a bisection of fk(x0)-0.5, where fk is evaluated using the stored values for gamma
    ! E. Rubensson et al., Chem Phys Lett 432, 2006, 591-594
    mu_a = 0.0_dp; mu_b = 1.0_dp; 
    mu_fa = evaluate_trs4_polynomial(mu_a, gamma_values, i-1)-0.5_dp
    DO j=1,40
      mu_c = 0.5*(mu_a + mu_b)
      mu_fc = evaluate_trs4_polynomial(mu_c, gamma_values, i-1)-0.5_dp ! i-1 because in the last iteration, only convergence is checked
      IF (ABS(mu_fc) < 1.0E-6_dp .OR. (mu_b-mu_a)/2 < 1.0E-6_dp) EXIT !TODO: define threshold values
      
      IF (mu_fc*mu_fa > 0) THEN
        mu_a = mu_c
        mu_fa = mu_fc
      ELSE
        mu_b = mu_c
      ENDIF
    ENDDO
    mu = (eps_min - eps_max)*mu_c + eps_max
    DEALLOCATE(gamma_values)
    IF (unit_nr > 0) THEN
       WRITE(unit_nr, '(T6,A,1X,F12.5)') 'Chemical potential (mu): ', mu
    ENDIF
    e_mu = mu

    CALL compute_homo_lumo(matrix_k0,matrix_x,eps_min,eps_max,threshold,max_iter_lanczos,eps_lanczos,homo,lumo,unit_nr,error)

    e_homo = homo
    e_lumo = lumo

    CALL cp_dbcsr_release(matrix_x, error=error)
    CALL cp_dbcsr_release(matrix_k0, error=error)
    CALL timestop(handle)

  END SUBROUTINE density_matrix_trs4

! *****************************************************************************
!> \brief compute the homo and lumo given a KS matrix and a density matrix in the orthonormalized basis
!>        and the eps_min and eps_max, min and max eigenvalue of the ks matrix
!> \par History
!>       2012.06 created [Florian Thoele]
!> \author Florian Thoele 
! *****************************************************************************
  SUBROUTINE compute_homo_lumo(matrix_k,matrix_p,eps_min,eps_max,threshold,max_iter_lanczos,eps_lanczos,homo,lumo,unit_nr,error)
    TYPE(cp_dbcsr_type)                      :: matrix_k, matrix_p
    REAL(KIND=dp)                            :: eps_min, eps_max, threshold
    INTEGER, INTENT(IN)                      :: max_iter_lanczos
    REAL(KIND=dp), INTENT(IN)                :: eps_lanczos
    REAL(KIND=dp)                            :: homo, lumo
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_homo_lumo', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: converged
    REAL(KIND=dp)                            :: max_eig, min_eig, shift1, &
                                                shift2
    TYPE(cp_dbcsr_type)                      :: tmp1, tmp2, tmp3

! temporary matrices used for HOMO/LUMO calculation

    CALL cp_dbcsr_init(tmp1,error=error)
    CALL cp_dbcsr_create(tmp1,template=matrix_k,matrix_type=dbcsr_type_no_symmetry,error=error)

    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=matrix_k,matrix_type=dbcsr_type_no_symmetry,error=error)

    CALL cp_dbcsr_init(tmp3,error=error)
    CALL cp_dbcsr_create(tmp3,template=matrix_k,matrix_type=dbcsr_type_no_symmetry,error=error)
    
    shift1 = -eps_min
    shift2 = eps_max

    ! find largest ev of P*(K+shift*1), where shift is the neg. val. of the smallest ev of K
    CALL cp_dbcsr_copy(tmp2, matrix_k, error=error)
    CALL cp_dbcsr_add_on_diag(tmp2, shift1, error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p, tmp2,&
                           0.0_dp, tmp1, filter_eps=threshold,error=error)
    CALL lanczos_alg_serial(tmp1, max_eig, min_eig, converged=converged, &
              threshold=eps_lanczos, max_iter=max_iter_lanczos, error=error)
    homo = max_eig-shift1
    IF (unit_nr > 0) THEN
       WRITE(unit_nr, '(T6,A,1X,L12)')   "Lanczos converged:      ", converged
    ENDIF

    ! -(1-P)*(K-shift*1) = (1-P)*(shift*1 - K), where shift is the largest ev of K
    CALL cp_dbcsr_copy(tmp3, matrix_p, error=error)
    CALL cp_dbcsr_scale(tmp3, -1.0_dp, error=error)
    CALL cp_dbcsr_add_on_diag(tmp3, 1.0_dp, error=error)  !tmp3 = 1-P
    CALL cp_dbcsr_copy(tmp2, matrix_k, error=error)
    CALL cp_dbcsr_add_on_diag(tmp2, -shift2, error=error)
    CALL cp_dbcsr_multiply("N", "N", -1.0_dp, tmp3, tmp2,&
                           0.0_dp, tmp1, filter_eps=threshold,error=error)
    CALL lanczos_alg_serial(tmp1, max_eig, min_eig, converged=converged,&
                            threshold=eps_lanczos, max_iter=max_iter_lanczos, error=error)
    lumo = -max_eig+shift2

    IF (unit_nr > 0) THEN
       WRITE(unit_nr, '(T6,A,1X,L12)')   "Lanczos converged:      ", converged
       WRITE(unit_nr, '(T6,A,1X,3F12.5)') 'HOMO/LUMO/gap', homo, lumo, lumo-homo
    ENDIF
    CALL cp_dbcsr_release(tmp1, error=error)
    CALL cp_dbcsr_release(tmp2, error=error)
    CALL cp_dbcsr_release(tmp3, error=error)

  END SUBROUTINE compute_homo_lumo
 
  FUNCTION evaluate_trs4_polynomial(x, gamma_values, i) RESULT (xr)
    REAL(KIND=dp)                            :: x
    REAL(KIND=dp), DIMENSION(:)              :: gamma_values
    INTEGER                                  :: i
    REAL(KIND=dp)                            :: xr

    REAL(KIND=dp), PARAMETER                 :: gam_max = 6.0_dp, &
                                                gam_min = 0.0_dp

    INTEGER                                  :: k

    xr = x
    DO k=1,i
      IF (gamma_values(k) > gam_max) THEN
        xr = 2*xr-xr**2
      ELSE IF(gamma_values(k) < gam_min) THEN
       xr = xr**2
      ELSE
        xr = (xr*xr)*(4*xr-3*xr*xr) + gamma_values(k)*xr*xr*((1-xr)**2)
      ENDIF
    ENDDO
  END FUNCTION evaluate_trs4_polynomial

  FUNCTION estimate_steps(homo, lumo, threshold) RESULT (steps)
    REAL(KIND=dp)                            :: homo, lumo, threshold
    INTEGER                                  :: steps

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: h, l, m

    l = lumo
    h = homo

    DO i=1,200
      IF (ABS(l) < threshold .AND. ABS(1-h) < threshold) EXIT
      m = 0.5_dp*(h+l)
      IF (m > 0.5_dp) THEN
        h = h**2
        l = l**2
      ELSE
        h = 2*h-h**2
        l = 2*l-l**2
      ENDIF
    ENDDO   
    steps = i-1
  END FUNCTION estimate_steps



 

END MODULE dm_ls_scf_methods


