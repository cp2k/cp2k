!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Adaptive Clenshaw-Curtis quadrature algorithm to integrate a complex-valued function in
!>        a complex plane
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
MODULE negf_integr_cc
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale,&
                                              cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE fft_tools,                       ONLY: fft_fw1d
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: z_one
   USE negf_integr_utils,               ONLY: contour_shape_arc,&
                                              contour_shape_linear,&
                                              equidistant_nodes_a_b,&
                                              rescale_nodes_cos,&
                                              rescale_normalised_nodes
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_integr_cc'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   INTEGER, PARAMETER, PUBLIC :: cc_interval_full = 0, &
                                 cc_interval_half = 1

   INTEGER, PARAMETER, PUBLIC :: cc_shape_linear = contour_shape_linear, &
                                 cc_shape_arc = contour_shape_arc

   PUBLIC :: ccquad_type

   PUBLIC :: ccquad_init, &
             ccquad_release, &
             ccquad_double_number_of_points, &
             ccquad_reduce_and_append_zdata, &
             ccquad_refine_integral

! **************************************************************************************************
!> \brief Adaptive Clenshaw-Curtis environment.
! **************************************************************************************************
   TYPE ccquad_type
      !> integration lower and upper bounds
      COMPLEX(kind=dp)                                   :: a, b
      !> integration interval:
      !>   cc_interval_full -- [a .. b],
      !>       grid density: 'a' .. .  .   .   .  . .. 'b';
      !>   cc_interval_half -- [a .. 2b-a], assuming int_{b}^{2b-a} f(x) dx = 0,
      !>       grid density: 'a' .. .  .   . 'b'
      INTEGER                                            :: interval_id
      !> integration shape
      INTEGER                                            :: shape_id
      !> estimated error
      REAL(kind=dp)                                      :: error
      !> approximate integral value
      TYPE(cp_cfm_type), POINTER                         :: integral
      !> error estimate for every element of the 'integral' matrix
      TYPE(cp_fm_type), POINTER                          :: error_fm
      !> weights associated with matrix elements; the 'error' variable contains the value Trace(error_fm * weights)
      TYPE(cp_fm_type), POINTER                          :: weights
      !> integrand value at grid points. Due to symmetry of Clenshaw-Curtis quadratures,
      !> we only need to keep the left half-interval
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: zdata_cache
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: tnodes
   END TYPE ccquad_type

CONTAINS

! **************************************************************************************************
!> \brief Initialise a Clenshaw-Curtis quadrature environment variable.
!> \param cc_env      environment variable to initialise
!> \param xnodes      points at which an integrand needs to be computed (initialised on exit)
!> \param nnodes      initial number of points to compute (initialised on exit)
!> \param a           integral lower bound
!> \param b           integral upper bound
!> \param interval_id full [-1 .. 1] or half [-1 .. 0] interval
!> \param shape_id    shape of a curve along which the integral will be evaluated
!> \param weights     weights associated with matrix elements; used to compute cumulative error
!> \param tnodes_restart list of nodes over the interval [-1 .. 1] from a previous integral evaluation.
!>                       If present, the same set of 'xnodes' will be used to compute this integral.
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
!> \note Clenshaw-Curtis quadratures are defined on the interval [-1 .. 1] and have non-uniforms node
!>       distribution which is symmetric and much sparse about 0. When the half-interval [-1 .. 0]
!>       is requested, the integrand value on another subinterval (0 .. 1] is assumed to be zero.
!>       Half interval mode is typically useful for rapidly decaying integrands (e.g. multiplied by
!>       Fermi function), so we do not actually need a fine grid spacing on this tail.
! **************************************************************************************************
   SUBROUTINE ccquad_init(cc_env, xnodes, nnodes, a, b, interval_id, shape_id, weights, tnodes_restart)
      TYPE(ccquad_type), INTENT(out)                     :: cc_env
      INTEGER, INTENT(inout)                             :: nnodes
      COMPLEX(kind=dp), DIMENSION(nnodes), INTENT(out)   :: xnodes
      COMPLEX(kind=dp), INTENT(in)                       :: a, b
      INTEGER, INTENT(in)                                :: interval_id, shape_id
      TYPE(cp_fm_type), POINTER                          :: weights
      REAL(kind=dp), DIMENSION(nnodes), INTENT(in), &
         OPTIONAL                                        :: tnodes_restart

      CHARACTER(len=*), PARAMETER :: routineN = 'ccquad_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icol, ipoint, irow, ncols, &
                                                            nnodes_half, nrows
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: w_data, w_data_my
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL timeset(routineN, handle)

      CPASSERT(nnodes > 2)
      CPASSERT(ASSOCIATED(weights))

      ! ensure that MOD(nnodes-1, 2) == 0
      nnodes = 2*((nnodes - 1)/2) + 1

      cc_env%interval_id = interval_id
      cc_env%shape_id = shape_id
      cc_env%a = a
      cc_env%b = b
      cc_env%error = HUGE(0.0_dp)

      NULLIFY (cc_env%integral, cc_env%error_fm, cc_env%weights)
      CALL cp_fm_get_info(weights, local_data=w_data, nrow_local=nrows, ncol_local=ncols, matrix_struct=fm_struct)
      CALL cp_fm_create(cc_env%weights, fm_struct)
      CALL cp_fm_get_info(cc_env%weights, local_data=w_data_my)

      ! use the explicit loop to avoid temporary arrays
      DO icol = 1, ncols
         DO irow = 1, nrows
            w_data_my(irow, icol) = ABS(w_data(irow, icol))
         END DO
      END DO

      SELECT CASE (interval_id)
      CASE (cc_interval_full)
         nnodes_half = nnodes/2 + 1
      CASE (cc_interval_half)
         nnodes_half = nnodes
      CASE DEFAULT
         CPABORT("Unimplemented interval type")
      END SELECT

      ALLOCATE (cc_env%tnodes(nnodes))

      IF (PRESENT(tnodes_restart)) THEN
         cc_env%tnodes(1:nnodes) = tnodes_restart(1:nnodes)
      ELSE
         CALL equidistant_nodes_a_b(-1.0_dp, 0.0_dp, nnodes_half, cc_env%tnodes)

         ! rescale all but the end-points, as they are transformed into themselves (-1.0 -> -1.0; 0.0 -> 0.0).
         ! Moreover, by applying this rescaling transformation to the end-points we cannot guarantee the exact
         ! result due to rounding errors in evaluation of COS function.
         IF (nnodes_half > 2) &
            CALL rescale_nodes_cos(nnodes_half - 2, cc_env%tnodes(2:))

         SELECT CASE (interval_id)
         CASE (cc_interval_full)
            ! reflect symmetric nodes
            DO ipoint = nnodes_half - 1, 1, -1
               cc_env%tnodes(nnodes_half + ipoint) = -cc_env%tnodes(nnodes_half - ipoint)
            END DO
         CASE (cc_interval_half)
            ! rescale half-interval : [-1 .. 0] -> [-1 .. 1]
            cc_env%tnodes(1:nnodes_half) = 2.0_dp*cc_env%tnodes(1:nnodes_half) + 1.0_dp
         END SELECT
      END IF

      CALL rescale_normalised_nodes(nnodes, cc_env%tnodes, a, b, shape_id, xnodes)

      CALL timestop(handle)
   END SUBROUTINE ccquad_init

! **************************************************************************************************
!> \brief Release a Clenshaw-Curtis quadrature environment variable.
!> \param cc_env   environment variable to release (modified on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE ccquad_release(cc_env)
      TYPE(ccquad_type), INTENT(inout)                   :: cc_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ccquad_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipoint

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(cc_env%error_fm)) &
         CALL cp_fm_release(cc_env%error_fm)

      IF (ASSOCIATED(cc_env%weights)) &
         CALL cp_fm_release(cc_env%weights)

      IF (ASSOCIATED(cc_env%integral)) &
         CALL cp_cfm_release(cc_env%integral)

      IF (ALLOCATED(cc_env%zdata_cache)) THEN
         DO ipoint = SIZE(cc_env%zdata_cache), 1, -1
            IF (ASSOCIATED(cc_env%zdata_cache(ipoint)%matrix)) &
               CALL cp_cfm_release(cc_env%zdata_cache(ipoint)%matrix)
         END DO

         DEALLOCATE (cc_env%zdata_cache)
      END IF

      IF (ALLOCATED(cc_env%tnodes)) DEALLOCATE (cc_env%tnodes)

      CALL timestop(handle)
   END SUBROUTINE ccquad_release

! **************************************************************************************************
!> \brief Get the next set of points at which the integrand needs to be computed. These points are
!>        then can be used to refine the integral approximation.
!> \param cc_env       environment variable (modified on exit)
!> \param xnodes_next  set of additional points (allocated and initialised on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE ccquad_double_number_of_points(cc_env, xnodes_next)
      TYPE(ccquad_type), INTENT(inout)                   :: cc_env
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(inout)                                   :: xnodes_next

      CHARACTER(len=*), PARAMETER :: routineN = 'ccquad_double_number_of_points', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipoint, nnodes_exist, &
                                                            nnodes_half, nnodes_next
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: tnodes, tnodes_old

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ALLOCATED(xnodes_next))
      CPASSERT(ASSOCIATED(cc_env%integral))
      CPASSERT(ASSOCIATED(cc_env%error_fm))
      CPASSERT(ALLOCATED(cc_env%zdata_cache))

      ! due to symmetry of Clenshaw-Curtis quadratures, we only need to keep the left half-interval [-1 .. 0]
      nnodes_exist = SIZE(cc_env%zdata_cache)
      ! new nodes will be placed between the existed ones, so the number of nodes
      ! on the left half-interval [-1 .. 0] is equal to nnodes_exist - 1
      nnodes_half = nnodes_exist - 1

      SELECT CASE (cc_env%interval_id)
      CASE (cc_interval_full)
         ! double number of nodes as we have 2 half-intervals [-1 .. 0] and [0 .. 1]
         nnodes_next = 2*nnodes_half
      CASE (cc_interval_half)
         nnodes_next = nnodes_half
      CASE DEFAULT
         CPABORT("Unimplemented interval type")
      END SELECT

      ALLOCATE (xnodes_next(nnodes_next))
      ALLOCATE (tnodes(nnodes_next))

      CALL equidistant_nodes_a_b(0.5_dp/REAL(nnodes_half, kind=dp) - 1.0_dp, &
                                 -0.5_dp/REAL(nnodes_half, kind=dp), &
                                 nnodes_half, tnodes)

      CALL rescale_nodes_cos(nnodes_half, tnodes)

      SELECT CASE (cc_env%interval_id)
      CASE (cc_interval_full)
         ! reflect symmetric nodes
         DO ipoint = 1, nnodes_half
            tnodes(nnodes_half + ipoint) = -tnodes(nnodes_half - ipoint + 1)
         END DO
      CASE (cc_interval_half)
         ! rescale half-interval : [-1 .. 0] -> [-1 .. 1]
         tnodes(1:nnodes_half) = 2.0_dp*tnodes(1:nnodes_half) + 1.0_dp
      END SELECT

      ! append new tnodes to the cache
      CALL MOVE_ALLOC(cc_env%tnodes, tnodes_old)
      nnodes_exist = SIZE(tnodes_old)

      ALLOCATE (cc_env%tnodes(nnodes_exist + nnodes_next))
      cc_env%tnodes(1:nnodes_exist) = tnodes_old(1:nnodes_exist)
      cc_env%tnodes(nnodes_exist + 1:nnodes_exist + nnodes_next) = tnodes(1:nnodes_next)
      DEALLOCATE (tnodes_old)

      ! rescale nodes [-1 .. 1] -> [a .. b] according to the shape
      CALL rescale_normalised_nodes(nnodes_next, tnodes, cc_env%a, cc_env%b, cc_env%shape_id, xnodes_next)

      DEALLOCATE (tnodes)
      CALL timestop(handle)
   END SUBROUTINE ccquad_double_number_of_points

! **************************************************************************************************
!> \brief Prepare Clenshaw-Curtis environment for the subsequent refinement of the integral.
!> \param cc_env       environment variable (modified on exit)
!> \param zdata_next   additional integrand value at additional points (modified on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
!> \note Due to symmetry of Clenshaw-Curtis quadratures (weight(x) == weight(-x)), we do not need to
!>       keep all the matrices from 'zdata_next', only 'zdata_next(x) + zdata_next(-x)' is needed.
!>       In order to reduce the number of matrix allocations, we move some of the matrices from the
!>       end of the 'zdata_new' array to the 'cc_env%zdata_cache' array, and nullify the corresponding
!>       pointers at 'zdata_next' array. So the calling subroutine need to release the remained
!>       matrices or reuse them but taking into account the missed ones.
! **************************************************************************************************
   SUBROUTINE ccquad_reduce_and_append_zdata(cc_env, zdata_next)
      TYPE(ccquad_type), INTENT(inout)                   :: cc_env
      TYPE(cp_cfm_p_type), DIMENSION(:), INTENT(inout)   :: zdata_next

      CHARACTER(len=*), PARAMETER :: routineN = 'ccquad_reduce_and_append_zdata', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: zscale
      INTEGER                                            :: handle, ipoint, nnodes_exist, &
                                                            nnodes_half, nnodes_next
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: zdata_tmp

      CALL timeset(routineN, handle)

      nnodes_next = SIZE(zdata_next)
      CPASSERT(nnodes_next > 0)

      ! compute weights of new points on a complex contour according to their values of the 't' parameter
      nnodes_exist = SIZE(cc_env%tnodes)
      CPASSERT(nnodes_exist >= nnodes_next)

      ALLOCATE (zscale(nnodes_next))
      CALL rescale_normalised_nodes(nnodes_next, cc_env%tnodes(nnodes_exist - nnodes_next + 1:nnodes_exist), &
                                    cc_env%a, cc_env%b, cc_env%shape_id, weights=zscale)

      IF (cc_env%interval_id == cc_interval_half) zscale(:) = 2.0_dp*zscale(:)

      ! rescale integrand values
      DO ipoint = 1, nnodes_next
         CALL cp_cfm_scale(zscale(ipoint), zdata_next(ipoint)%matrix)
      END DO
      DEALLOCATE (zscale)

      ! squash points with the same clenshaw-curtis weights together
      IF (ALLOCATED(cc_env%zdata_cache)) THEN
         nnodes_exist = SIZE(cc_env%zdata_cache)
      ELSE
         nnodes_exist = 0
      END IF

      SELECT CASE (cc_env%interval_id)
      CASE (cc_interval_full)
         IF (ALLOCATED(cc_env%zdata_cache)) THEN
            CPASSERT(nnodes_exist == nnodes_next/2 + 1)
            nnodes_half = nnodes_exist - 1
         ELSE
            CPASSERT(MOD(nnodes_next, 2) == 1)
            nnodes_half = nnodes_next/2 + 1
         END IF
      CASE (cc_interval_half)
         IF (ALLOCATED(cc_env%zdata_cache)) THEN
            CPASSERT(nnodes_exist == nnodes_next + 1)
         END IF

         nnodes_half = nnodes_next
      END SELECT

      IF (cc_env%interval_id == cc_interval_full) THEN
         DO ipoint = nnodes_next/2, 1, -1
            CALL cp_cfm_scale_and_add(z_one, zdata_next(ipoint)%matrix, z_one, zdata_next(nnodes_next - ipoint + 1)%matrix)
         END DO
      END IF

      IF (ALLOCATED(cc_env%zdata_cache)) THEN
         ! note that nnodes_half+1 == nnodes_exist for both half- and full-intervals
         ALLOCATE (zdata_tmp(nnodes_half + nnodes_exist))

         DO ipoint = 1, nnodes_half
            zdata_tmp(2*ipoint - 1)%matrix => cc_env%zdata_cache(ipoint)%matrix
            zdata_tmp(2*ipoint)%matrix => zdata_next(ipoint)%matrix
            NULLIFY (zdata_next(ipoint)%matrix)
         END DO
         zdata_tmp(nnodes_half + nnodes_exist)%matrix => cc_env%zdata_cache(nnodes_exist)%matrix

         DEALLOCATE (cc_env%zdata_cache)
         CALL MOVE_ALLOC(zdata_tmp, cc_env%zdata_cache)
      ELSE
         CALL cp_cfm_scale(2.0_dp, zdata_next(nnodes_half)%matrix)

         ALLOCATE (cc_env%zdata_cache(nnodes_half))

         DO ipoint = 1, nnodes_half
            cc_env%zdata_cache(ipoint)%matrix => zdata_next(ipoint)%matrix
            NULLIFY (zdata_next(ipoint)%matrix)
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE ccquad_reduce_and_append_zdata

! **************************************************************************************************
!> \brief Refine approximated integral.
!> \param cc_env       environment variable (modified on exit)
!> \par History
!>   * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE ccquad_refine_integral(cc_env)
      TYPE(ccquad_type), INTENT(inout)                   :: cc_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ccquad_refine_integral', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: ztmp, ztmp_dct
      INTEGER :: handle, icol, ipoint, irow, ncols_local, nintervals, nintervals_half, &
         nintervals_half_plus_1, nintervals_half_plus_2, nintervals_plus_2, nrows_local, stat
      LOGICAL                                            :: equiv
      REAL(kind=dp)                                      :: rscale
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

!      TYPE(fft_plan_type)                                :: fft_plan
!      INTEGER(kind=int_8)                                :: plan

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(cc_env%zdata_cache))

      nintervals_half_plus_1 = SIZE(cc_env%zdata_cache)
      nintervals_half = nintervals_half_plus_1 - 1
      nintervals_half_plus_2 = nintervals_half_plus_1 + 1
      nintervals = 2*nintervals_half
      nintervals_plus_2 = nintervals + 2
      CPASSERT(nintervals_half > 1)

      IF (.NOT. ASSOCIATED(cc_env%integral)) THEN
         CALL cp_cfm_get_info(cc_env%zdata_cache(1)%matrix, matrix_struct=fm_struct)
         equiv = cp_fm_struct_equivalent(fm_struct, cc_env%weights%matrix_struct)
         CPASSERT(equiv)

         CALL cp_cfm_create(cc_env%integral, fm_struct)
         CALL cp_fm_create(cc_env%error_fm, fm_struct)
      END IF

      IF (debug_this_module) THEN
         DO ipoint = 1, nintervals_half_plus_1
            equiv = cp_fm_struct_equivalent(cc_env%zdata_cache(ipoint)%matrix%matrix_struct, cc_env%integral%matrix_struct)
            CPASSERT(equiv)
         END DO
      END IF

      CALL cp_cfm_get_info(cc_env%integral, nrow_local=nrows_local, ncol_local=ncols_local)

      ALLOCATE (weights(nintervals_half))

      ! omit the trivial weights(1) = 0.5
      DO ipoint = 2, nintervals_half
         rscale = REAL(2*(ipoint - 1), kind=dp)
         weights(ipoint) = 1.0_dp/(1.0_dp - rscale*rscale)
      END DO
      ! weights(1) <- weights(intervals_half + 1)
      rscale = REAL(nintervals, kind=dp)
      weights(1) = 1.0_dp/(1.0_dp - rscale*rscale)

      ! 1.0 / nintervals
      rscale = 1.0_dp/rscale

      ALLOCATE (ztmp(nintervals, nrows_local, ncols_local))
      ALLOCATE (ztmp_dct(nintervals, nrows_local, ncols_local))

!$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(icol, ipoint, irow), &
!$OMP             SHARED(cc_env, ncols_local, nintervals_half, nintervals_half_plus_1, nintervals_half_plus_2, nrows_local, ztmp)
      DO icol = 1, ncols_local
         DO irow = 1, nrows_local
            DO ipoint = 1, nintervals_half_plus_1
               ztmp(ipoint, irow, icol) = cc_env%zdata_cache(ipoint)%matrix%local_data(irow, icol)
            END DO

            DO ipoint = 2, nintervals_half
               ztmp(nintervals_half + ipoint, irow, icol) = ztmp(nintervals_half_plus_2 - ipoint, irow, icol)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL fft_fw1d(nintervals, nrows_local*ncols_local, .FALSE., ztmp, ztmp_dct, 1.0_dp, stat)
      IF (stat /= 0) THEN
         CALL cp_abort(__LOCATION__, &
                       "An FFT library is required for Clenshaw-Curtis quadrature. "// &
                       "You can use an alternative integration method instead.")
      END IF

!$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(icol, ipoint, irow), &
!$OMP             SHARED(cc_env, rscale, ncols_local, nintervals_half, nintervals_half_plus_1, nintervals_plus_2), &
!$OMP             SHARED(nrows_local, weights, ztmp_dct)
      DO icol = 1, ncols_local
         DO irow = 1, nrows_local
            ztmp_dct(1, irow, icol) = 0.5_dp*ztmp_dct(1, irow, icol)
            DO ipoint = 2, nintervals_half
               ztmp_dct(ipoint, irow, icol) = 0.5_dp*weights(ipoint)*(ztmp_dct(ipoint, irow, icol) + &
                                                                      ztmp_dct(nintervals_plus_2 - ipoint, irow, icol))
            END DO
            ztmp_dct(nintervals_half_plus_1, irow, icol) = weights(1)*ztmp_dct(nintervals_half_plus_1, irow, icol)

            cc_env%integral%local_data(irow, icol) = rscale*accurate_sum(ztmp_dct(1:nintervals_half_plus_1, irow, icol))
            cc_env%error_fm%local_data(irow, icol) = rscale*ABS(ztmp_dct(nintervals_half_plus_1, irow, icol))
         END DO
      END DO
!$OMP END PARALLEL DO

      DEALLOCATE (ztmp, ztmp_dct)

      CALL cp_fm_trace(cc_env%error_fm, cc_env%weights, cc_env%error)

      DEALLOCATE (weights)
      CALL timestop(handle)
   END SUBROUTINE ccquad_refine_integral

END MODULE negf_integr_cc
