!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Minimax-Ewald (MME) method for calculating 2-center and 3-center
!>        electron repulsion integrals (ERI) of periodic systems using a
!>        Hermite Gaussian basis.
!>        The method relies on analytical Fourier transforms of Cartesian and
!>        Hermite Gaussian functions and Poisson summation formula to represent
!>        ERIs as a discrete sum over direct lattice vectors or reciprocal
!>        lattice vectors. The reciprocal space potential 1/G^2 is approximated
!>        by a linear combination of Gaussians employing minimax approximation.
!>        Not yet implemented: 3c ERIs for nonorthogonal cells.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_integrate
   USE ao_util,                         ONLY: exp_radius
   USE eri_mme_gaussian,                ONLY: hermite_gauss_norm
   USE eri_mme_lattice_summation,       ONLY: &
        ellipsoid_bounds, eri_mme_2c_get_bounds, eri_mme_2c_get_rads, eri_mme_3c_get_bounds, &
        eri_mme_3c_get_rads, get_l, pgf_sum_2c_gspace_1d, pgf_sum_2c_gspace_3d, &
        pgf_sum_2c_rspace_1d, pgf_sum_2c_rspace_3d, pgf_sum_3c_1d, pgf_sum_3c_3d
   USE eri_mme_types,                   ONLY: eri_mme_param,&
                                              get_minimax_from_cutoff
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_integrate'

   PUBLIC :: eri_mme_2c_integrate, eri_mme_3c_integrate

CONTAINS

! **************************************************************************************************
!> \brief Low-level integration routine for 2-center ERIs.
!> \param param ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param zeta ...
!> \param zetb ...
!> \param rab ...
!> \param hab ...
!> \param o1 ...
!> \param o2 ...
!> \param G_count ...
!> \param R_count ...
!> \param normalize calculate integrals w.r.t. normalized Hermite-Gaussians
!> \param potential use exact potential instead of minimax approx. (for testing only)
!> \param pot_par potential parameter
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate(param, la_min, la_max, lb_min, lb_max, zeta, zetb, rab, &
                                   hab, o1, o2, G_count, R_count, normalize, potential, pot_par)
      TYPE(eri_mme_param), INTENT(IN)                    :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: hab
      INTEGER, INTENT(IN)                                :: o1, o2
      INTEGER, INTENT(INOUT), OPTIONAL                   :: G_count, R_count
      LOGICAL, INTENT(IN), OPTIONAL                      :: normalize
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_2c_integrate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, bx, by, bz, grid, i_aw, &
                                                            i_xyz, ico, jco, l_max, la, lb, n_aw
      INTEGER(KIND=int_8), DIMENSION(2)                  :: n_sum_3d
      INTEGER(KIND=int_8), DIMENSION(2, 3)               :: n_sum_1d
      INTEGER, DIMENSION(3)                              :: la_xyz, lb_xyz
      LOGICAL                                            :: do_g_sum, exact, is_ortho, norm
      REAL(KIND=dp)                                      :: alpha_G, alpha_R, cutoff, G_rad, G_res, &
                                                            Imm, inv_lgth, Ixyz, lgth, max_error, &
                                                            prefac, R_rad, R_res, vol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: S_G_1, S_G_2, S_G_no, S_G_no_H
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: S_G
      REAL(KIND=dp), DIMENSION(3)                        :: G_bounds, R_bounds
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aw

      CPASSERT(param%is_valid)

      grid = 0
      CALL eri_mme_2c_get_rads(la_max, lb_max, zeta, zetb, 0.0_dp, param%G_min, param%R_min, param%sum_precision, &
                               G_rad=G_rad)
      cutoff = G_rad**2/2
      CALL get_minimax_from_cutoff(param%minimax_grid, cutoff, n_aw, aw, grid)

      CPASSERT(grid .GT. 0)

      ! cell info
      h_inv = param%h_inv
      hmat = param%hmat
      vol = param%vol

      IF (PRESENT(normalize)) THEN
         norm = normalize
      ELSE
         norm = .FALSE.
      ENDIF

      l_max = la_max + lb_max

      IF (PRESENT(potential)) THEN
         exact = .TRUE.
      ELSE
         exact = .FALSE.
      ENDIF

      IF (exact) THEN
         is_ortho = .FALSE.
      ELSE
         is_ortho = param%is_ortho
      ENDIF

      IF (is_ortho) THEN
         ALLOCATE (S_G(0:l_max, 3, n_aw))
         S_G = 0.0_dp

         IF (param%debug) THEN
            ALLOCATE (S_G_1(0:l_max))
            ALLOCATE (S_G_2(0:l_max))
         ENDIF
      ELSE
         ALLOCATE (S_G_no(ncoset(l_max)))
         S_G_no(:) = 0.0_dp
         ALLOCATE (S_G_no_H(ncoset(l_max)))
      ENDIF

      IF (exact) THEN
         alpha_G = 0.25_dp/zeta + 0.25_dp/zetb
         ! resolution for Gaussian width
         G_res = 0.5_dp*param%G_min
         R_res = 0.5_dp*param%R_min

         G_rad = exp_radius(l_max, alpha_G, 0.01*param%sum_precision, 1.0_dp, epsabs=G_res)
         G_bounds(:) = ellipsoid_bounds(G_rad, TRANSPOSE(hmat)/(2.0_dp*pi))
         CALL pgf_sum_2c_gspace_3d(S_G_no, l_max, -rab, alpha_G, h_inv, G_bounds, G_rad, vol, potential, pot_par)
      ELSE

         DO i_aw = 1, n_aw

            CALL eri_mme_2c_get_bounds(hmat, h_inv, vol, is_ortho, param%G_min, param%R_min, la_max, lb_max, &
                                       zeta, zetb, aw(i_aw), param%sum_precision, n_sum_1d, n_sum_3d, &
                                       G_bounds, G_rad, R_bounds, R_rad)
            alpha_G = aw(i_aw) + 0.25_dp/zeta + 0.25_dp/zetb
            alpha_R = 0.25_dp/alpha_G
            IF (is_ortho) THEN ! orthorhombic cell

               ! 1) precompute Ewald-like sum

               DO i_xyz = 1, 3
                  lgth = ABS(hmat(i_xyz, i_xyz))
                  inv_lgth = ABS(h_inv(i_xyz, i_xyz))

                  ! perform sum in R or G space. Choose the space in which less summands are required for convergence
                  do_g_sum = n_sum_1d(1, i_xyz) < n_sum_1d(2, i_xyz) !G_bounds < R_bounds

                  IF (do_g_sum) THEN
                     CALL pgf_sum_2c_gspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_G, inv_lgth, G_bounds(i_xyz))
                     IF (PRESENT(G_count)) G_count = G_count + 1
                  ELSE
                     CALL pgf_sum_2c_rspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_R, lgth, R_bounds(i_xyz))
                     IF (PRESENT(R_count)) R_count = R_count + 1
                  ENDIF

                  IF (param%debug) THEN
                     ! check consistency of summation methods
                     CALL pgf_sum_2c_gspace_1d(S_G_1, -rab(i_xyz), alpha_G, inv_lgth, G_bounds(i_xyz))
                     CALL pgf_sum_2c_rspace_1d(S_G_2, -rab(i_xyz), alpha_R, lgth, R_bounds(i_xyz))
                     max_error = MAXVAL(ABS(S_G_1 - S_G_2)/(0.5_dp*(ABS(S_G_1) + ABS(S_G_2)) + 1.0_dp))

                     CPASSERT(max_error .LE. param%debug_delta)
                  ENDIF
               ENDDO

            ELSE ! general cell

               do_g_sum = n_sum_3d(1) < n_sum_3d(2) !PRODUCT(2*R_bounds+1) .GT. PRODUCT(2*G_bounds+1)

               IF (do_g_sum) THEN
                  CALL pgf_sum_2c_gspace_3d(S_G_no_H, l_max, -rab, alpha_G, h_inv, G_bounds, G_rad, vol)
                  IF (PRESENT(G_count)) G_count = G_count + 1
               ELSE
                  CALL pgf_sum_2c_rspace_3d(S_G_no_H, l_max, -rab, alpha_R, hmat, h_inv, R_bounds, R_rad)
                  IF (PRESENT(R_count)) R_count = R_count + 1
               ENDIF
               S_G_no(:) = S_G_no(:) + aw(n_aw + i_aw)*S_G_no_H
            ENDIF
         ENDDO
      ENDIF

      ! prefactor for integral values (unnormalized Hermite Gaussians)
      prefac = SQRT(1.0_dp/(zeta*zetb))

      ! 2) Assemble integral values from Ewald sums
      DO jco = ncoset(lb_min - 1) + 1, ncoset(lb_max)
         CALL get_l(jco, lb, bx, by, bz)
         lb_xyz = [bx, by, bz]
         DO ico = ncoset(la_min - 1) + 1, ncoset(la_max)
            CALL get_l(ico, la, ax, ay, az)
            la_xyz = [ax, ay, az]
            IF (is_ortho) THEN
               Imm = 0.0_dp
               DO i_aw = 1, n_aw
                  Ixyz = 1.0_dp
                  DO i_xyz = 1, 3
                     Ixyz = Ixyz*S_G(la_xyz(i_xyz) + lb_xyz(i_xyz), i_xyz, i_aw)*prefac
                  ENDDO
                  Imm = Imm + aw(n_aw + i_aw)*Ixyz
               ENDDO
            ELSE
               Imm = S_G_no(coset(ax + bx, ay + by, az + bz))*prefac**3
            ENDIF
            IF (la + lb .EQ. 0 .AND. .NOT. exact) THEN
               Imm = Imm - SUM(aw(n_aw + 1:2*n_aw))*prefac**3/vol ! subtracting G = 0 term
            ENDIF
            IF (.NOT. norm) THEN
               ! rescaling needed due to Hermite Gaussians (such that they can be contracted same way as Cartesian Gaussians)
               ! and factor of 4 pi**4 (-1)**lb
               hab(o1 + ico, o2 + jco) = Imm*4.0_dp*pi**4/((2.0_dp*zeta)**la*(-2.0_dp*zetb)**lb)
            ELSE
               ! same thing for normalized Hermite Gaussians
               hab(o1 + ico, o2 + jco) = Imm*4.0_dp*pi**4*(-1.0_dp)**lb*hermite_gauss_norm(zeta, la_xyz)* &
                                         hermite_gauss_norm(zetb, lb_xyz)
            ENDIF
         ENDDO ! la
      ENDDO ! lb

   END SUBROUTINE eri_mme_2c_integrate

! **************************************************************************************************
!> \brief Low-level integration routine for 3-center ERIs
!> \param param ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param habc ...
!> \param o1 ...
!> \param o2 ...
!> \param o3 ...
!> \param GG_count ...
!> \param GR_count ...
!> \param RR_count ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_integrate(param, la_min, la_max, lb_min, lb_max, lc_min, lc_max, zeta, zetb, zetc, RA, RB, RC, &
                                   habc, o1, o2, o3, GG_count, GR_count, RR_count)
      TYPE(eri_mme_param), INTENT(IN)                    :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: RA, RB, RC
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN)                                :: o1, o2, o3
      INTEGER, INTENT(INOUT), OPTIONAL                   :: GG_count, GR_count, RR_count

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_3c_integrate', &
         routineP = moduleN//':'//routineN

      CPASSERT(param%is_valid)
      IF (param%is_ortho) THEN
         CALL eri_mme_3c_integrate_ortho(param, la_min, la_max, lb_min, lb_max, lc_min, lc_max, zeta, zetb, zetc, RA, RB, RC, &
                                         habc, o1, o2, o3, RR_count)

      ELSE
         CALL eri_mme_3c_integrate_nonortho(param, la_min, la_max, lb_min, lb_max, lc_min, lc_max, zeta, zetb, zetc, RA, RB, RC, &
                                            habc, o1, o2, o3, GG_count, GR_count, RR_count)

      ENDIF
   END SUBROUTINE eri_mme_3c_integrate

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param habc ...
!> \param o1 ...
!> \param o2 ...
!> \param o3 ...
!> \param RR_count ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_integrate_ortho(param, la_min, la_max, lb_min, lb_max, lc_min, lc_max, zeta, zetb, zetc, RA, RB, RC, &
                                         habc, o1, o2, o3, RR_count)
      TYPE(eri_mme_param), INTENT(IN)                    :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: RA, RB, RC
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN)                                :: o1, o2, o3
      INTEGER, INTENT(INOUT), OPTIONAL                   :: RR_count

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_3c_integrate_ortho', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: grid, i_aw, lmax_0, n_aw
      INTEGER(KIND=int_8), DIMENSION(3)                  :: n_sum_3d
      INTEGER(KIND=int_8), DIMENSION(3, 3)               :: n_sum_1d
      REAL(KIND=dp)                                      :: alpha_R_0, cutoff, G_res, lgth, prefac, &
                                                            R_rad_0, R_res, vol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: S_G_0
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: S_G
      REAL(KIND=dp), DIMENSION(2)                        :: R_rads_3
      REAL(KIND=dp), DIMENSION(2, 3)                     :: R_bounds_3
      REAL(KIND=dp), DIMENSION(3)                        :: G_rads_1, R_rads_2
      REAL(KIND=dp), DIMENSION(3, 3)                     :: G_bounds_1, h_inv, hmat, R_bounds_2
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aw

      grid = 0

      CALL eri_mme_3c_get_rads(la_max, lb_max, lc_max, zeta, zetb, zetc, 0.0_dp, param%G_min, param%R_min, &
                               param%sum_precision, G_rads_1=G_rads_1)

      cutoff = (MIN(G_rads_1(1), G_rads_1(2) + G_rads_1(3)))**2/2

      CALL get_minimax_from_cutoff(param%minimax_grid, cutoff, n_aw, aw, grid)

      CPASSERT(grid .GT. 0)

      ! minimax coeffs
      n_aw = param%minimax_grid(grid)%n_minimax
      aw => param%minimax_grid(grid)%minimax_aw

      ! cell info
      h_inv = param%h_inv
      hmat = param%hmat
      vol = param%vol

      ! prefactor for integral values (unnormalized Hermite Gaussians)
      prefac = (zeta*zetb*zetc)**(-1.5_dp)*pi**(11.0_dp/2.0_dp)*4.0_dp

      ! Preparations for G=0 component
      G_res = 0.5_dp*param%G_min
      R_res = 0.5_dp*param%R_min

      ALLOCATE (S_G(n_aw, 3, 0:la_max, 0:lb_max, 0:lc_max))

      ! G= 0 components
      IF (lc_min == 0) THEN
         ALLOCATE (S_G_0(0:la_max + lb_max, 3))

         alpha_R_0 = zeta*zetb/(zeta + zetb)
         lmax_0 = la_max + lb_max
         R_rad_0 = exp_radius(lmax_0, alpha_R_0, param%sum_precision, 1.0_dp, epsabs=R_res)

         lgth = ABS(hmat(1, 1))
         CALL pgf_sum_2c_rspace_1d(S_G_0(:, 1), RB(1) - RA(1), alpha_R_0, lgth, R_rad_0/lgth)
         lgth = ABS(hmat(2, 2))
         CALL pgf_sum_2c_rspace_1d(S_G_0(:, 2), RB(2) - RA(2), alpha_R_0, lgth, R_rad_0/lgth)
         lgth = ABS(hmat(3, 3))
         CALL pgf_sum_2c_rspace_1d(S_G_0(:, 3), RB(3) - RA(3), alpha_R_0, lgth, R_rad_0/lgth)
      ENDIF

      DO i_aw = 1, n_aw
         CALL eri_mme_3c_get_bounds(hmat, h_inv, vol, .TRUE., param%G_min, param%R_min, la_max, lb_max, lc_max, &
                                    zeta, zetb, zetc, aw(i_aw), param%sum_precision, n_sum_1d, n_sum_3d, &
                                    G_bounds_1, G_rads_1, R_bounds_2, R_rads_2, R_bounds_3, R_rads_3)
         CALL pgf_sum_3c_1d(S_G(i_aw, 1, :, :, :), RA(1), RB(1), RC(1), &
                            zeta, zetb, zetc, aw(i_aw), ABS(hmat(1, 1)), &
                            R_bounds_3(:, 1))

         CALL pgf_sum_3c_1d(S_G(i_aw, 2, :, :, :), RA(2), RB(2), RC(2), &
                            zeta, zetb, zetc, aw(i_aw), ABS(hmat(2, 2)), &
                            R_bounds_3(:, 2))

         CALL pgf_sum_3c_1d(S_G(i_aw, 3, :, :, :), RA(3), RB(3), RC(3), &
                            zeta, zetb, zetc, aw(i_aw), ABS(hmat(3, 3)), &
                            R_bounds_3(:, 3))

         IF (PRESENT(RR_count)) RR_count = RR_count + 3
      ENDDO

      CALL eri_mme_3c_collect_ortho(vol, la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                                    zeta, zetb, zetc, n_aw, aw, S_G, S_G_0, prefac, habc, o1, o2, o3)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param RA ...
!> \param RB ...
!> \param RC ...
!> \param habc ...
!> \param o1 ...
!> \param o2 ...
!> \param o3 ...
!> \param GG_count ...
!> \param GR_count ...
!> \param RR_count ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_integrate_nonortho(param, la_min, la_max, lb_min, lb_max, lc_min, lc_max, zeta, zetb, zetc, RA, RB, RC, &
                                            habc, o1, o2, o3, GG_count, GR_count, RR_count)

      TYPE(eri_mme_param), INTENT(IN)                    :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: RA, RB, RC
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN)                                :: o1, o2, o3
      INTEGER, INTENT(INOUT), OPTIONAL                   :: GG_count, GR_count, RR_count

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_3c_integrate_nonortho', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, bx, by, bz, cx, cy, cz, &
                                                            grid, i_aw, ico, ir, jco, kco, la, lb, &
                                                            lc, lmax_0, method, n_aw, nresults, &
                                                            sum_method
      INTEGER(KIND=int_8), DIMENSION(3)                  :: n_sum_3d
      INTEGER(KIND=int_8), DIMENSION(3, 3)               :: n_sum_1d
      LOGICAL                                            :: db_sum1, db_sum2, db_sum3, do_g_sum_0
      REAL(KIND=dp)                                      :: alpha_G_0, alpha_R_0, cutoff, G_rad_0, &
                                                            G_res, max_error, max_result, &
                                                            min_result, prefac, R_rad_0, R_res, vol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: results_no, S_G_0_no
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: S_G_no, S_G_no_1, S_G_no_2, S_G_no_3, &
                                                            S_G_no_H
      REAL(KIND=dp), DIMENSION(2)                        :: R_rads_3
      REAL(KIND=dp), DIMENSION(2, 3)                     :: R_bounds_3
      REAL(KIND=dp), DIMENSION(3)                        :: G_bound_0, G_rads_1, R_0, R_bound_0, &
                                                            R_rads_2
      REAL(KIND=dp), DIMENSION(3, 3)                     :: G_bounds_1, h_inv, hmat, R_bounds_2
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aw

      CPASSERT(param%is_valid)

      grid = 0

      CALL eri_mme_3c_get_rads(la_max, lb_max, lc_max, zeta, zetb, zetc, 0.0_dp, param%G_min, param%R_min, &
                               param%sum_precision, G_rads_1=G_rads_1)

      cutoff = (MIN(G_rads_1(1), G_rads_1(2) + G_rads_1(3)))**2/2

      CALL get_minimax_from_cutoff(param%minimax_grid, cutoff, n_aw, aw, grid)

      CPASSERT(grid .GT. 0)

      ! minimax coeffs
      n_aw = param%minimax_grid(grid)%n_minimax
      aw => param%minimax_grid(grid)%minimax_aw

      ! cell info
      h_inv = param%h_inv
      hmat = param%hmat
      vol = param%vol

      ! prefactor for integral values (unnormalized Hermite Gaussians)
      prefac = (zeta*zetb*zetc)**(-1.5_dp)*pi**(11.0_dp/2.0_dp)*4.0_dp

      IF (param%debug) THEN
         ALLOCATE (S_G_no_1(ncoset(la_max), ncoset(lb_max), ncoset(lc_max)))
         ALLOCATE (S_G_no_2(ncoset(la_max), ncoset(lb_max), ncoset(lc_max)))
         ALLOCATE (S_G_no_3(ncoset(la_max), ncoset(lb_max), ncoset(lc_max)))
      ENDIF

      ! Preparations for G=0 component
      G_res = 0.5_dp*param%G_min
      R_res = 0.5_dp*param%R_min

      ALLOCATE (S_G_no(ncoset(la_max), ncoset(lb_max), ncoset(lc_max)))

      S_G_no(:, :, :) = 0.0_dp
      IF (param%debug) THEN
         S_G_no_1(:, :, :) = -1.0_dp
         S_G_no_2(:, :, :) = -1.0_dp
         S_G_no_3(:, :, :) = -1.0_dp
      ENDIF
      ALLOCATE (S_G_no_H(ncoset(la_max), ncoset(lb_max), ncoset(lc_max)))

      ! G= 0 components
      IF (lc_min == 0) THEN
         ALLOCATE (S_G_0_no(ncoset(la_max + lb_max)))
         alpha_G_0 = 0.25_dp/zetb + 0.25_dp/zeta
         alpha_R_0 = 0.25_dp/alpha_G_0
         lmax_0 = la_max + lb_max
         R_0 = RB - RA
         G_rad_0 = exp_radius(lmax_0, alpha_G_0, param%sum_precision, 1.0_dp, epsabs=G_res)
         R_rad_0 = exp_radius(lmax_0, alpha_R_0, param%sum_precision, 1.0_dp, epsabs=R_res)
         G_bound_0 = ellipsoid_bounds(G_rad_0, TRANSPOSE(hmat)/(2.0_dp*pi))
         R_bound_0 = ellipsoid_bounds(R_rad_0, h_inv)
         do_g_sum_0 = PRODUCT(2*R_bound_0 + 1) .GT. PRODUCT(2*G_bound_0 + 1)
         IF (do_g_sum_0) THEN
            CALL pgf_sum_2c_gspace_3d(S_G_0_no, lmax_0, R_0, alpha_G_0, h_inv, G_bound_0, G_rad_0, vol)
         ELSE
            CALL pgf_sum_2c_rspace_3d(S_G_0_no, lmax_0, R_0, alpha_R_0, hmat, h_inv, R_bound_0, R_rad_0)
         ENDIF
      ENDIF

      DO i_aw = 1, n_aw
         CALL eri_mme_3c_get_bounds(hmat, h_inv, vol, .FALSE., param%G_min, param%R_min, la_max, lb_max, lc_max, &
                                    zeta, zetb, zetc, aw(i_aw), param%sum_precision, n_sum_1d, n_sum_3d, &
                                    G_bounds_1, G_rads_1, R_bounds_2, R_rads_2, R_bounds_3, R_rads_3)
         sum_method = MINLOC(n_sum_3d, DIM=1)

         CALL pgf_sum_3c_3d(S_G_no_H, la_max, lb_max, lc_max, RA, RB, RC, &
                            zeta, zetb, zetc, aw(i_aw), hmat, h_inv, vol, &
                            G_bounds_1, R_bounds_2, R_bounds_3, &
                            G_rads_1, R_rads_2, R_rads_3, &
                            method=sum_method, method_out=method)
         S_G_no(:, :, :) = S_G_no(:, :, :) + aw(n_aw + i_aw)*S_G_no_H(:, :, :)

         SELECT CASE (method)
         CASE (1)
            IF (PRESENT(GG_count)) GG_count = GG_count + 1
         CASE (2)
            IF (PRESENT(GR_count)) GR_count = GR_count + 1
         CASE (3)
            IF (PRESENT(RR_count)) RR_count = RR_count + 1
         CASE DEFAULT
            CPABORT("")
         END SELECT

         IF (param%debug) THEN
            nresults = 0
            ! check consistency of summation methods

            db_sum1 = (n_sum_3d(1)) .LT. INT(param%debug_nsum, KIND=int_8)
            db_sum2 = (n_sum_3d(2)) .LT. INT(param%debug_nsum, KIND=int_8)
            db_sum3 = (n_sum_3d(3)) .LT. INT(param%debug_nsum, KIND=int_8)

            IF (param%unit_nr > 0) THEN
               WRITE (param%unit_nr, *) "ERI_MME DEBUG | number of summands (GG / GR / RR)", n_sum_3d
               WRITE (param%unit_nr, *) "ERI_MME DEBUG | sum methods to be compared (GG / GR / RR)", db_sum1, db_sum2, db_sum3
            ENDIF

            S_G_no_1(:, :, :) = 0.0_dp
            S_G_no_2(:, :, :) = 0.0_dp
            S_G_no_3(:, :, :) = 0.0_dp

            IF (db_sum1) THEN
               CALL pgf_sum_3c_3d(S_G_no_1, la_max, lb_max, lc_max, RA, RB, RC, &
                                  zeta, zetb, zetc, aw(i_aw), hmat, h_inv, vol, &
                                  G_bounds_1, R_bounds_2, R_bounds_3, &
                                  G_rads_1, R_rads_2, R_rads_3, &
                                  method=1)
               nresults = nresults + 1
            ENDIF

            IF (db_sum2) THEN
               CALL pgf_sum_3c_3d(S_G_no_2, la_max, lb_max, lc_max, RA, RB, RC, &
                                  zeta, zetb, zetc, aw(i_aw), hmat, h_inv, vol, &
                                  G_bounds_1, R_bounds_2, R_bounds_3, &
                                  G_rads_1, R_rads_2, R_rads_3, &
                                  method=2)
               nresults = nresults + 1
            ENDIF

            IF (db_sum3) THEN
               CALL pgf_sum_3c_3d(S_G_no_3, la_max, lb_max, lc_max, RA, RB, RC, &
                                  zeta, zetb, zetc, aw(i_aw), hmat, h_inv, vol, &
                                  G_bounds_1, R_bounds_2, R_bounds_3, &
                                  G_rads_1, R_rads_2, R_rads_3, &
                                  method=3)
               nresults = nresults + 1
            ENDIF

            max_error = 0.0_dp
            ALLOCATE (results_no(nresults))

            DO kco = ncoset(lc_min - 1) + 1, ncoset(lc_max)
               CALL get_l(kco, lc, cx, cy, cz)
               DO jco = ncoset(lb_min - 1) + 1, ncoset(lb_max)
                  CALL get_l(jco, lb, bx, by, bz)
                  DO ico = ncoset(la_min - 1) + 1, ncoset(la_max)
                     CALL get_l(ico, la, ax, ay, az)

                     max_error = 0.0_dp
                     ir = 0
                     IF (db_sum1) THEN
                        ir = ir + 1
                        results_no(ir) = S_G_no_1(ico, jco, kco)
                     ENDIF

                     IF (db_sum2) THEN
                        ir = ir + 1
                        results_no(ir) = S_G_no_2(ico, jco, kco)
                     ENDIF

                     IF (db_sum3) THEN
                        ir = ir + 1
                        results_no(ir) = S_G_no_3(ico, jco, kco)
                     ENDIF

                     max_result = MAXVAL(results_no)
                     min_result = MINVAL(results_no)
                     IF (nresults > 0) max_error = MAX(max_error, &
                                                    (max_result - min_result)/(0.5_dp*(ABS(max_result) + ABS(min_result)) + 1.0_dp))
                  ENDDO
               ENDDO
            ENDDO

            CPASSERT(max_error .LE. param%debug_delta)
            DEALLOCATE (results_no)
         ENDIF
      ENDDO

      ! assemble integral values

      CALL eri_mme_3c_collect_nonortho(vol, la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                                       zeta, zetb, zetc, n_aw, aw, S_G_no, S_G_0_no, prefac, habc, o1, o2, o3)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param vol ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param n_aw ...
!> \param aw ...
!> \param S_G ...
!> \param S_G_0 ...
!> \param prefac ...
!> \param habc ...
!> \param o1 ...
!> \param o2 ...
!> \param o3 ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_collect_ortho(vol, la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                                       zeta, zetb, zetc, n_aw, aw, S_G, S_G_0, prefac, habc, o1, o2, o3)
      REAL(KIND=dp), INTENT(IN)                          :: vol
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc
      INTEGER, INTENT(IN)                                :: n_aw
      REAL(KIND=dp), DIMENSION(2*n_aw), INTENT(IN)       :: aw
      REAL(KIND=dp), DIMENSION(:, :, 0:, 0:, 0:), &
         INTENT(IN)                                      :: S_G
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: S_G_0
      REAL(KIND=dp), INTENT(IN)                          :: prefac
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN)                                :: o1, o2, o3

      INTEGER                                            :: ax, ay, az, bx, by, bz, cx, cy, cz, ico, &
                                                            jco, kco, la, la_prev, lb, lb_prev, &
                                                            lc, lc_prev
      REAL(KIND=dp)                                      :: Imm, Ixyz_0, mone_lb, resc_a, &
                                                            resc_a_init, resc_b, resc_b_init, &
                                                            resc_c, resc_c_init

      ! Initialization of rescaling factors due to Hermite Gaussians
      resc_a_init = (2.0_dp*zeta)**la_min
      resc_b_init = (2.0_dp*zetb)**lb_min
      resc_c_init = (2.0_dp*zetc)**lc_min

      resc_c = resc_c_init
      lc_prev = lc_min
      DO kco = ncoset(lc_min - 1) + 1, ncoset(lc_max)
         CALL get_l(kco, lc, cx, cy, cz)
         IF (lc_prev < lc) resc_c = resc_c*(2.0_dp*zetc)

         resc_b = resc_b_init
         lb_prev = lb_min
         DO jco = ncoset(lb_min - 1) + 1, ncoset(lb_max)
            CALL get_l(jco, lb, bx, by, bz)
            mone_lb = (-1.0_dp)**lb
            IF (lb_prev < lb) resc_b = resc_b*(2.0_dp*zetb)

            resc_a = resc_a_init
            la_prev = la_min
            DO ico = ncoset(la_min - 1) + 1, ncoset(la_max)
               CALL get_l(ico, la, ax, ay, az)

               IF (la_prev < la) resc_a = resc_a*(2.0_dp*zeta)
               Ixyz_0 = 0.0_dp
               IF (lc == 0) THEN
                  Ixyz_0 = S_G_0(ax + bx, 1)* &
                           S_G_0(ay + by, 2)* &
                           S_G_0(az + bz, 3) &
                           /vol*mone_lb
               ENDIF
               Imm = SUM(aw(n_aw + 1:2*n_aw)*( &
                         S_G(1:n_aw, 1, ax, bx, cx)* &
                         S_G(1:n_aw, 2, ay, by, cy)* &
                         S_G(1:n_aw, 3, az, bz, cz)) - Ixyz_0)

               ! rescaling needed due to Hermite Gaussians
               habc(o1 + ico, o2 + jco, o3 + kco) = Imm*prefac/(resc_a*resc_b*resc_c)
               la_prev = la
            ENDDO ! la
            lb_prev = lb
         ENDDO ! lb
         lc_prev = lc
      ENDDO ! lc

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param vol ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param n_aw ...
!> \param aw ...
!> \param S_G ...
!> \param S_G_0 ...
!> \param prefac ...
!> \param habc ...
!> \param o1 ...
!> \param o2 ...
!> \param o3 ...
! **************************************************************************************************
   PURE SUBROUTINE eri_mme_3c_collect_nonortho(vol, la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                                               zeta, zetb, zetc, n_aw, aw, S_G, S_G_0, prefac, habc, o1, o2, o3)
      REAL(KIND=dp), INTENT(IN)                          :: vol
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc
      INTEGER, INTENT(IN)                                :: n_aw
      REAL(KIND=dp), DIMENSION(2*n_aw), INTENT(IN)       :: aw
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: S_G
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: S_G_0
      REAL(KIND=dp), INTENT(IN)                          :: prefac
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN)                                :: o1, o2, o3

      INTEGER                                            :: ax, ay, az, bx, by, bz, cx, cy, cz, ico, &
                                                            ijco, jco, kco, la, la_prev, lb, &
                                                            lb_prev, lc, lc_prev
      REAL(KIND=dp)                                      :: Imm, mone_lb, resc_a, resc_a_init, &
                                                            resc_b, resc_b_init, resc_c, &
                                                            resc_c_init

      ! Initialization of rescaling factors due to Hermite Gaussians
      resc_a_init = (2.0_dp*zeta)**la_min
      resc_b_init = (2.0_dp*zetb)**lb_min
      resc_c_init = (2.0_dp*zetc)**lc_min

      resc_c = resc_c_init
      lc_prev = lc_min
      DO kco = ncoset(lc_min - 1) + 1, ncoset(lc_max)
         CALL get_l(kco, lc, cx, cy, cz)
         IF (lc_prev < lc) resc_c = resc_c*(2.0_dp*zetc)

         resc_b = resc_b_init
         lb_prev = lb_min
         DO jco = ncoset(lb_min - 1) + 1, ncoset(lb_max)
            CALL get_l(jco, lb, bx, by, bz)
            mone_lb = (-1.0_dp)**lb
            IF (lb_prev < lb) resc_b = resc_b*(2.0_dp*zetb)

            resc_a = resc_a_init
            la_prev = la_min
            DO ico = ncoset(la_min - 1) + 1, ncoset(la_max)
               CALL get_l(ico, la, ax, ay, az)

               IF (la_prev < la) resc_a = resc_a*(2.0_dp*zeta)
               IF (lc .GT. 0) THEN
                  Imm = S_G(ico, jco, kco)
               ELSE
                  ijco = coset(ax + bx, ay + by, az + bz)
                  Imm = S_G(ico, jco, kco) - SUM(aw(n_aw + 1:2*n_aw))*S_G_0(ijco)/vol*mone_lb
               ENDIF

               ! rescaling needed due to Hermite Gaussians
               habc(o1 + ico, o2 + jco, o3 + kco) = Imm*prefac/(resc_a*resc_b*resc_c)
               la_prev = la
            ENDDO ! la
            lb_prev = lb
         ENDDO ! lb
         lc_prev = lc
      ENDDO ! lc

   END SUBROUTINE

END MODULE eri_mme_integrate
