!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Adaptive Simpson's rule algorithm to integrate a complex-valued function in a complex plane
! **************************************************************************************************
MODULE negf_integr_simpson
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale,&
                                              cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_retain,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              z_one,&
                                              z_zero
   USE negf_integr_utils,               ONLY: contour_shape_arc,&
                                              contour_shape_linear,&
                                              equidistant_nodes_a_b,&
                                              rescale_normalised_nodes
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_integr_simpson'
   ! adaptive Simpson method requires 5 points per subinterval for the error estimate.
   ! So, in principle, at the end we can compute the value of the integral using
   ! Boole's rule and possibly improve the actual accuracy by up to one order of magnitude.
   LOGICAL, PARAMETER, PRIVATE          :: is_boole = .FALSE.

   INTEGER, PARAMETER, PUBLIC :: sr_shape_linear = contour_shape_linear, &
                                 sr_shape_arc = contour_shape_arc

   PUBLIC :: simpsonrule_type
   PUBLIC :: simpsonrule_init, simpsonrule_release, simpsonrule_get_next_nodes, simpsonrule_refine_integral

! **************************************************************************************************
!> \brief A structure to store data for non-converged sub-interval.
! **************************************************************************************************
   TYPE simpsonrule_subinterval_type
      !> unscaled lower and upper boundaries within the interval [-1 .. 1]
      REAL(kind=dp)                                      :: lb, ub
      !> target accuracy for this sub-interval
      REAL(kind=dp)                                      :: conv
      !> estimated error value on this sub-interval
      REAL(kind=dp)                                      :: error
      !> integrand values at equally spaced points [a, b, c, d, e] located on the curve shape([lb .. ub])
      TYPE(cp_cfm_type), POINTER                         :: fa, fb, fc, fd, fe
   END TYPE simpsonrule_subinterval_type

! **************************************************************************************************
!> \brief A structure to store data needed for adaptive Simpson's rule algorithm.
! **************************************************************************************************
   TYPE simpsonrule_type
      !> lower and upper boundaries of the curve on the complex plane
      COMPLEX(kind=dp)                                   :: a, b
      !> ID number which determines the shape of a curve along which the integral will be evaluated
      INTEGER                                            :: shape_id
      !> target accuracy
      REAL(kind=dp)                                      :: conv
      !> estimated error value on the entire integration interval,
      !> as well as on converged sub-intervals only
      REAL(kind=dp)                                      :: error, error_conv
      !> the estimated value of the integral on the entire interval
      TYPE(cp_cfm_type), POINTER                         :: integral
      !> work matrix to store the contribution to the integral on converged sub-intervals
      TYPE(cp_cfm_type), POINTER                         :: integral_conv
      !> work matrices which stores approximated integral computed by using a/b/c, c/d/e, and a/c/e points respectively
      TYPE(cp_cfm_type), POINTER                         :: integral_abc, integral_cde, integral_ace
      !> work matrix to temporarily store error estimate of the integral on a sub-interval for every matrix element
      TYPE(cp_fm_type), POINTER                          :: error_fm
      !> weights associated with matrix elements; the final error is computed as Trace(error_fm * weights)
      TYPE(cp_fm_type), POINTER                          :: weights
      ! non-converged sub-intervals
      TYPE(simpsonrule_subinterval_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subintervals
      !> complete list of nodes over the normalised interval [-1 .. 1] needed to restart
      !> Useful when a series of similar integrals need to be computed at an identical set
      !> of points, so intermediate quantities can be saved and reused.
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: tnodes
   END TYPE simpsonrule_type

   COMPLEX(kind=dp), PARAMETER, PRIVATE :: z_four = 4.0_dp*z_one

CONTAINS

! **************************************************************************************************
!> \brief Initialise a Simpson's rule environment variable.
!> \param sr_env   Simpson's rule environment (initialised on exit)
!> \param xnodes   points at which an integrand needs to be computed (initialised on exit)
!> \param nnodes   initial number of points to compute (initialised on exit)
!> \param a        integral lower boundary
!> \param b        integral upper boundary
!> \param shape_id shape of a curve along which the integral will be evaluated
!> \param conv     convergence threshold
!> \param weights  weights associated with matrix elements; used to compute cumulative error
!> \param tnodes_restart list of nodes over the interval [-1 .. 1] from a previous integral evaluation.
!>                       If present, the same set of 'xnodes' will be used to compute this integral.
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
!> \note When we integrate the retarded Green's function times the Fermi function over the energy
!>       domain and pass the overlap matrix (S) as the 'weights' matrix, the convergence threshold
!>       ('conv') becomes the maximum error in the total number of electrons multiplied by pi.
! **************************************************************************************************
   SUBROUTINE simpsonrule_init(sr_env, xnodes, nnodes, a, b, shape_id, conv, weights, tnodes_restart)
      TYPE(simpsonrule_type), INTENT(out)                :: sr_env
      INTEGER, INTENT(inout)                             :: nnodes
      COMPLEX(kind=dp), DIMENSION(nnodes), INTENT(out)   :: xnodes
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      INTEGER, INTENT(in)                                :: shape_id
      REAL(kind=dp), INTENT(in)                          :: conv
      TYPE(cp_fm_type), POINTER                          :: weights
      REAL(kind=dp), DIMENSION(nnodes), INTENT(in), &
         OPTIONAL                                        :: tnodes_restart

      CHARACTER(len=*), PARAMETER :: routineN = 'simpsonrule_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icol, irow, ncols, nrows
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: w_data, w_data_my
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      CPASSERT(nnodes > 4)
      CPASSERT(ASSOCIATED(weights))

      ! ensure that MOD(nnodes-1, 4) == 0
      nnodes = 4*((nnodes - 1)/4) + 1

      sr_env%shape_id = shape_id
      sr_env%a = a
      sr_env%b = b
      sr_env%conv = conv
      sr_env%error = HUGE(0.0_dp)
      sr_env%error_conv = 0.0_dp

      NULLIFY (sr_env%error_fm, sr_env%weights)
      CALL cp_fm_get_info(weights, local_data=w_data, nrow_local=nrows, ncol_local=ncols, matrix_struct=fm_struct)
      CALL cp_fm_create(sr_env%error_fm, fm_struct)
      CALL cp_fm_create(sr_env%weights, fm_struct)
      CALL cp_fm_get_info(sr_env%weights, local_data=w_data_my)

      ! use the explicit loop to avoid temporary arrays. The magic constant 15.0 is due to Simpson's rule error analysis.
      DO icol = 1, ncols
         DO irow = 1, nrows
            w_data_my(irow, icol) = ABS(w_data(irow, icol))/15.0_dp
         END DO
      END DO

      NULLIFY (sr_env%integral, sr_env%integral_conv)
      NULLIFY (sr_env%integral_abc, sr_env%integral_cde, sr_env%integral_ace)

      ALLOCATE (sr_env%tnodes(nnodes))

      IF (PRESENT(tnodes_restart)) THEN
         sr_env%tnodes(1:nnodes) = tnodes_restart(1:nnodes)
      ELSE
         CALL equidistant_nodes_a_b(-1.0_dp, 1.0_dp, nnodes, sr_env%tnodes)
      END IF
      CALL rescale_normalised_nodes(nnodes, sr_env%tnodes, a, b, shape_id, xnodes)

      CALL timestop(handle)
   END SUBROUTINE simpsonrule_init

! **************************************************************************************************
!> \brief Release a Simpson's rule environment variable.
!> \param sr_env   Simpson's rule environment (modified on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE simpsonrule_release(sr_env)
      TYPE(simpsonrule_type), INTENT(inout)              :: sr_env

      CHARACTER(len=*), PARAMETER :: routineN = 'simpsonrule_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, interval

      CALL timeset(routineN, handle)
      IF (ALLOCATED(sr_env%subintervals)) THEN
         DO interval = SIZE(sr_env%subintervals), 1, -1
            IF (ASSOCIATED(sr_env%subintervals(interval)%fa)) CALL cp_cfm_release(sr_env%subintervals(interval)%fa)
            IF (ASSOCIATED(sr_env%subintervals(interval)%fb)) CALL cp_cfm_release(sr_env%subintervals(interval)%fb)
            IF (ASSOCIATED(sr_env%subintervals(interval)%fc)) CALL cp_cfm_release(sr_env%subintervals(interval)%fc)
            IF (ASSOCIATED(sr_env%subintervals(interval)%fd)) CALL cp_cfm_release(sr_env%subintervals(interval)%fd)
            IF (ASSOCIATED(sr_env%subintervals(interval)%fe)) CALL cp_cfm_release(sr_env%subintervals(interval)%fe)
         END DO

         DEALLOCATE (sr_env%subintervals)
      END IF

      IF (ASSOCIATED(sr_env%integral)) CALL cp_cfm_release(sr_env%integral)
      IF (ASSOCIATED(sr_env%integral_conv)) CALL cp_cfm_release(sr_env%integral_conv)
      IF (ASSOCIATED(sr_env%integral_abc)) CALL cp_cfm_release(sr_env%integral_abc)
      IF (ASSOCIATED(sr_env%integral_cde)) CALL cp_cfm_release(sr_env%integral_cde)
      IF (ASSOCIATED(sr_env%integral_ace)) CALL cp_cfm_release(sr_env%integral_ace)

      IF (ASSOCIATED(sr_env%error_fm)) CALL cp_fm_release(sr_env%error_fm)
      IF (ASSOCIATED(sr_env%weights)) CALL cp_fm_release(sr_env%weights)

      IF (ALLOCATED(sr_env%tnodes)) DEALLOCATE (sr_env%tnodes)

      CALL timestop(handle)
   END SUBROUTINE simpsonrule_release

! **************************************************************************************************
!> \brief Get the next set of nodes where to compute integrand.
!> \param sr_env      Simpson's rule environment (modified on exit)
!> \param xnodes_next list of additional points (initialised on exit)
!> \param nnodes      actual number of points to compute (modified on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
!> \note The number of nodes returned is limited by the initial value of the nnodes variable;
!>       un exit nnodes == 0 means that the target accuracy has been achieved.
! **************************************************************************************************
   SUBROUTINE simpsonrule_get_next_nodes(sr_env, xnodes_next, nnodes)
      TYPE(simpsonrule_type), INTENT(inout)              :: sr_env
      INTEGER, INTENT(inout)                             :: nnodes
      COMPLEX(kind=dp), DIMENSION(nnodes), INTENT(out)   :: xnodes_next

      CHARACTER(len=*), PARAMETER :: routineN = 'simpsonrule_get_next_nodes', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nnodes_old
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: tnodes, tnodes_old

      CALL timeset(routineN, handle)
      ALLOCATE (tnodes(nnodes))

      CALL simpsonrule_get_next_nodes_real(sr_env, tnodes, nnodes)
      IF (nnodes > 0) THEN
         CALL MOVE_ALLOC(sr_env%tnodes, tnodes_old)
         nnodes_old = SIZE(tnodes_old)

         ALLOCATE (sr_env%tnodes(nnodes_old + nnodes))
         sr_env%tnodes(1:nnodes_old) = tnodes_old(1:nnodes_old)
         sr_env%tnodes(nnodes_old + 1:nnodes_old + nnodes) = tnodes(1:nnodes)
         DEALLOCATE (tnodes_old)

         CALL rescale_normalised_nodes(nnodes, tnodes, sr_env%a, sr_env%b, sr_env%shape_id, xnodes_next)
      END IF

      DEALLOCATE (tnodes)
      CALL timestop(handle)
   END SUBROUTINE simpsonrule_get_next_nodes

! **************************************************************************************************
!> \brief Low level routine that returns unscaled nodes on interval [-1 .. 1].
!> \param sr_env       Simpson's rule environment
!> \param xnodes_unity list of additional unscaled nodes (initialised on exit)
!> \param nnodes       actual number of points to compute (initialised on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE simpsonrule_get_next_nodes_real(sr_env, xnodes_unity, nnodes)
      TYPE(simpsonrule_type), INTENT(in)                 :: sr_env
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: xnodes_unity
      INTEGER, INTENT(out)                               :: nnodes

      CHARACTER(len=*), PARAMETER :: routineN = 'simpsonrule_get_next_nodes_real', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, interval, nintervals

      CALL timeset(routineN, handle)

      IF (ALLOCATED(sr_env%subintervals)) THEN
         nintervals = SIZE(sr_env%subintervals)
      ELSE
         nintervals = 0
      END IF

      IF (nintervals > 0) THEN
         IF (SIZE(xnodes_unity) < 4*nintervals) &
            nintervals = SIZE(xnodes_unity)/4

         DO interval = 1, nintervals
            xnodes_unity(4*interval - 3) = 0.125_dp* &
                                           (7.0_dp*sr_env%subintervals(interval)%lb + sr_env%subintervals(interval)%ub)
            xnodes_unity(4*interval - 2) = 0.125_dp* &
                                           (5.0_dp*sr_env%subintervals(interval)%lb + 3.0_dp*sr_env%subintervals(interval)%ub)
            xnodes_unity(4*interval - 1) = 0.125_dp* &
                                           (3.0_dp*sr_env%subintervals(interval)%lb + 5.0_dp*sr_env%subintervals(interval)%ub)
            xnodes_unity(4*interval) = 0.125_dp*(sr_env%subintervals(interval)%lb + 7.0_dp*sr_env%subintervals(interval)%ub)
         END DO
      END IF

      nnodes = 4*nintervals
      CALL timestop(handle)
   END SUBROUTINE simpsonrule_get_next_nodes_real

! **************************************************************************************************
!> \brief Compute integral using the simpson's rules.
!> \param sr_env     Simpson's rule environment
!> \param zdata_next precomputed integrand values at points xnodes_next (nullified on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE simpsonrule_refine_integral(sr_env, zdata_next)
      TYPE(simpsonrule_type), INTENT(inout)              :: sr_env
      TYPE(cp_cfm_p_type), DIMENSION(:), INTENT(inout)   :: zdata_next

      CHARACTER(len=*), PARAMETER :: routineN = 'simpsonrule_refine_integral', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: zscale
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: error_zdata
      INTEGER                                            :: handle, interval, ipoint, nintervals, &
                                                            nintervals_exist, npoints
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp)                                      :: rscale
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: errors
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: error_rdata
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(simpsonrule_subinterval_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subintervals

      CALL timeset(routineN, handle)

      npoints = SIZE(zdata_next)
      IF (ASSOCIATED(sr_env%integral)) THEN
         ! we need 4 new points per subinterval (p, q, r, s)
         !   p   q   r   s
         ! a . b . c . d . e
         CPASSERT(npoints > 0 .AND. MOD(npoints, 4) == 0)
      ELSE
         ! first call: need 4*n+1 points
         ! a1 b1 c1 d1 e1
         !             a2 b2 c2 d2 e2
         !                         a3 b3 c3 d3 e3
         CPASSERT(npoints > 1 .AND. MOD(npoints, 4) == 1)
      END IF

      ! compute weights of new points on a complex contour according to their values of the 't' parameter
      nintervals_exist = SIZE(sr_env%tnodes)
      CPASSERT(nintervals_exist >= npoints)
      ALLOCATE (zscale(npoints))

      CALL rescale_normalised_nodes(npoints, sr_env%tnodes(nintervals_exist - npoints + 1:nintervals_exist), &
                                    sr_env%a, sr_env%b, sr_env%shape_id, weights=zscale)

      ! rescale integrand values
      DO ipoint = 1, npoints
         CALL cp_cfm_scale(zscale(ipoint), zdata_next(ipoint)%matrix)
      END DO

      DEALLOCATE (zscale)

      ! insert new points
      nintervals = npoints/4
      IF (ASSOCIATED(sr_env%integral)) THEN
         ! subdivide existing intervals
         nintervals_exist = SIZE(sr_env%subintervals)
         CPASSERT(nintervals <= nintervals_exist)

         ALLOCATE (subintervals(nintervals_exist + nintervals))

         DO interval = 1, nintervals
            subintervals(2*interval - 1)%lb = sr_env%subintervals(interval)%lb
            subintervals(2*interval - 1)%ub = 0.5_dp*(sr_env%subintervals(interval)%lb + sr_env%subintervals(interval)%ub)
            subintervals(2*interval - 1)%conv = 0.5_dp*sr_env%subintervals(interval)%conv
            subintervals(2*interval - 1)%fa => sr_env%subintervals(interval)%fa
            subintervals(2*interval - 1)%fb => zdata_next(4*interval - 3)%matrix
            subintervals(2*interval - 1)%fc => sr_env%subintervals(interval)%fb
            subintervals(2*interval - 1)%fd => zdata_next(4*interval - 2)%matrix
            subintervals(2*interval - 1)%fe => sr_env%subintervals(interval)%fc
            CALL cp_cfm_retain(subintervals(2*interval - 1)%fe)

            subintervals(2*interval)%lb = subintervals(2*interval - 1)%ub
            subintervals(2*interval)%ub = sr_env%subintervals(interval)%ub
            subintervals(2*interval)%conv = subintervals(2*interval - 1)%conv
            subintervals(2*interval)%fa => sr_env%subintervals(interval)%fc
            subintervals(2*interval)%fb => zdata_next(4*interval - 1)%matrix
            subintervals(2*interval)%fc => sr_env%subintervals(interval)%fd
            subintervals(2*interval)%fd => zdata_next(4*interval)%matrix
            subintervals(2*interval)%fe => sr_env%subintervals(interval)%fe

            NULLIFY (zdata_next(4*interval - 3)%matrix)
            NULLIFY (zdata_next(4*interval - 2)%matrix)
            NULLIFY (zdata_next(4*interval - 1)%matrix)
            NULLIFY (zdata_next(4*interval)%matrix)
         END DO

         DO interval = nintervals + 1, nintervals_exist
            subintervals(interval + nintervals) = sr_env%subintervals(interval)
         END DO
         DEALLOCATE (sr_env%subintervals)
      ELSE
         ! first time -- allocate matrices and create a new set of intervals
         CALL cp_cfm_get_info(zdata_next(1)%matrix, matrix_struct=fm_struct)
         CALL cp_cfm_create(sr_env%integral, fm_struct)
         CALL cp_cfm_create(sr_env%integral_conv, fm_struct)
         CALL cp_cfm_create(sr_env%integral_abc, fm_struct)
         CALL cp_cfm_create(sr_env%integral_cde, fm_struct)
         CALL cp_cfm_create(sr_env%integral_ace, fm_struct)

         CALL cp_cfm_set_all(sr_env%integral_conv, z_zero)

         ALLOCATE (subintervals(nintervals))

         rscale = 1.0_dp/REAL(nintervals, kind=dp)

         DO interval = 1, nintervals
            ! lower bound: point with indices 1, 5, 9, ..., 4*nintervals+1
            subintervals(interval)%lb = sr_env%tnodes(4*interval - 3)
            subintervals(interval)%ub = sr_env%tnodes(4*interval + 1)
            subintervals(interval)%conv = rscale*sr_env%conv

            subintervals(interval)%fa => zdata_next(4*interval - 3)%matrix
            subintervals(interval)%fb => zdata_next(4*interval - 2)%matrix
            subintervals(interval)%fc => zdata_next(4*interval - 1)%matrix
            subintervals(interval)%fd => zdata_next(4*interval)%matrix
            subintervals(interval)%fe => zdata_next(4*interval + 1)%matrix
            CALL cp_cfm_retain(subintervals(interval)%fe)
         END DO

         CALL cp_cfm_release(zdata_next(4*nintervals + 1)%matrix)
      END IF

      ! we kept the originals matrices for internal use, so nullify the pointers
      ! to prevent  alteration of the matrices from the outside
      DO ipoint = 1, npoints
         NULLIFY (zdata_next(ipoint)%matrix)
      END DO

      CALL cp_fm_get_info(sr_env%error_fm, local_data=error_rdata)
      CALL cp_cfm_get_info(sr_env%integral_ace, local_data=error_zdata)

      ! do actual integration
      CALL cp_cfm_to_cfm(sr_env%integral_conv, sr_env%integral)
      sr_env%error = sr_env%error_conv
      nintervals_exist = SIZE(subintervals)

      DO interval = 1, nintervals_exist
         rscale = subintervals(interval)%ub - subintervals(interval)%lb
         CALL do_simpson_rule(sr_env%integral_ace, &
                              subintervals(interval)%fa, &
                              subintervals(interval)%fc, &
                              subintervals(interval)%fe, &
                              -0.5_dp*rscale)
         CALL do_simpson_rule(sr_env%integral_abc, &
                              subintervals(interval)%fa, &
                              subintervals(interval)%fb, &
                              subintervals(interval)%fc, &
                              0.25_dp*rscale)
         CALL do_simpson_rule(sr_env%integral_cde, &
                              subintervals(interval)%fc, &
                              subintervals(interval)%fd, &
                              subintervals(interval)%fe, &
                              0.25_dp*rscale)

         CALL cp_cfm_scale_and_add(z_one, sr_env%integral_abc, z_one, sr_env%integral_cde)
         CALL cp_cfm_scale_and_add(z_one, sr_env%integral_ace, z_one, sr_env%integral_abc)

         IF (is_boole) THEN
            CALL do_boole_rule(sr_env%integral_abc, &
                               subintervals(interval)%fa, &
                               subintervals(interval)%fb, &
                               subintervals(interval)%fc, &
                               subintervals(interval)%fd, &
                               subintervals(interval)%fe, &
                               0.5_dp*rscale, sr_env%integral_cde)
         END IF

         CALL cp_cfm_scale_and_add(z_one, sr_env%integral, z_one, sr_env%integral_abc)

         ! sr_env%error_fm = ABS(sr_env%integral_ace); no temporary arrays as pointers have different types
         error_rdata(:, :) = ABS(error_zdata(:, :))
         CALL cp_fm_trace(sr_env%error_fm, sr_env%weights, subintervals(interval)%error)

         sr_env%error = sr_env%error + subintervals(interval)%error

         ! add contributions from converged subintervals, so we could drop them afterward
         IF (subintervals(interval)%error <= subintervals(interval)%conv) THEN
            CALL cp_cfm_scale_and_add(z_one, sr_env%integral_conv, z_one, sr_env%integral_abc)
            sr_env%error_conv = sr_env%error_conv + subintervals(interval)%error
         END IF
      END DO

      IF (sr_env%error <= sr_env%conv) THEN
         ! We have already reached the target accuracy, so we can drop all subintervals
         ! (even those where local convergence has not been achieved). From now on environment
         ! components 'sr_env%error' and 'sr_env%integral_conv' hold incorrect values,
         ! but they should not been accessed from the outside anyway
         ! (uncomment the following two lines if they are actually need)

         ! sr_env%error_conv = sr_env%error
         ! CALL cp_cfm_to_cfm(sr_env%integral, sr_env%integral_conv)

         DO interval = nintervals_exist, 1, -1
            CALL cp_cfm_release(subintervals(interval)%fa)
            CALL cp_cfm_release(subintervals(interval)%fb)
            CALL cp_cfm_release(subintervals(interval)%fc)
            CALL cp_cfm_release(subintervals(interval)%fd)
            CALL cp_cfm_release(subintervals(interval)%fe)
         END DO
      ELSE
         ! sort subinterval according to their convergence, and drop convergent ones
         ALLOCATE (errors(nintervals_exist), inds(nintervals_exist))

         nintervals = 0
         DO interval = 1, nintervals_exist
            errors(interval) = subintervals(interval)%error

            IF (subintervals(interval)%error > subintervals(interval)%conv) &
               nintervals = nintervals + 1
         END DO

         CALL sort(errors, nintervals_exist, inds)

         IF (nintervals > 0) &
            ALLOCATE (sr_env%subintervals(nintervals))

         nintervals = 0
         DO ipoint = nintervals_exist, 1, -1
            interval = inds(ipoint)

            IF (subintervals(interval)%error > subintervals(interval)%conv) THEN
               nintervals = nintervals + 1

               sr_env%subintervals(nintervals) = subintervals(interval)
            ELSE
               CALL cp_cfm_release(subintervals(interval)%fa)
               CALL cp_cfm_release(subintervals(interval)%fb)
               CALL cp_cfm_release(subintervals(interval)%fc)
               CALL cp_cfm_release(subintervals(interval)%fd)
               CALL cp_cfm_release(subintervals(interval)%fe)
            END IF
         END DO

         DEALLOCATE (errors, inds)
      END IF

      DEALLOCATE (subintervals)

      CALL timestop(handle)
   END SUBROUTINE simpsonrule_refine_integral

! **************************************************************************************************
!> \brief Approximate value of the integral on subinterval [a .. c] using the Simpson's rule.
!> \param integral   approximated integral = length / 6 * (fa + 4*fb + fc) (initialised on exit)
!> \param fa         integrand value at point a
!> \param fb         integrand value at point b = (a + c) / 2
!> \param fc         integrand value at point c
!> \param length     distance between points a and c [ABS(c-a)]
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE do_simpson_rule(integral, fa, fb, fc, length)
      TYPE(cp_cfm_type), POINTER                         :: integral, fa, fb, fc
      REAL(kind=dp), INTENT(in)                          :: length

      CALL cp_cfm_to_cfm(fa, integral)
      CALL cp_cfm_scale_and_add(z_one, integral, z_four, fb)
      CALL cp_cfm_scale_and_add(z_one, integral, z_one, fc)
      CALL cp_cfm_scale(length/6.0_dp, integral)
   END SUBROUTINE do_simpson_rule

! **************************************************************************************************
!> \brief Approximate value of the integral on subinterval [a .. e] using the Boole's rule.
!> \param integral   approximated integral = length / 90 * (7*fa + 32*fb + 12*fc + 32*fd + 7*fe)
!>                   (initialised on exit)
!> \param fa         integrand value at point a
!> \param fb         integrand value at point b = a + (e-a)/4
!> \param fc         integrand value at point c = a + (e-a)/2
!> \param fd         integrand value at point d = a + 3*(e-a)/4
!> \param fe         integrand value at point e
!> \param length     distance between points a and e [ABS(e-a)]
!> \param work       work matrix
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE do_boole_rule(integral, fa, fb, fc, fd, fe, length, work)
      TYPE(cp_cfm_type), POINTER                         :: integral, fa, fb, fc, fd, fe
      REAL(kind=dp), INTENT(in)                          :: length
      TYPE(cp_cfm_type), POINTER                         :: work

      REAL(kind=dp)                                      :: rscale

      rscale = length/90.0_dp

      CALL cp_cfm_to_cfm(fc, integral)
      CALL cp_cfm_scale(12.0_dp*rscale, integral)

      CALL cp_cfm_to_cfm(fa, work)
      CALL cp_cfm_scale_and_add(z_one, work, z_one, fe)
      CALL cp_cfm_scale(7.0_dp*rscale, work)
      CALL cp_cfm_scale_and_add(z_one, integral, z_one, work)

      CALL cp_cfm_to_cfm(fb, work)
      CALL cp_cfm_scale_and_add(z_one, work, z_one, fd)
      CALL cp_cfm_scale(32.0_dp*rscale, work)
      CALL cp_cfm_scale_and_add(z_one, integral, z_one, work)
   END SUBROUTINE do_boole_rule
END MODULE negf_integr_simpson
