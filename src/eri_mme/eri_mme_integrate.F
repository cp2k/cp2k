!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
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
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              abc_contract
   USE ao_util,                         ONLY: exp_radius
   USE eri_mme_gaussian,                ONLY: hermite_gauss_norm
   USE eri_mme_lattice_summation,       ONLY: lattice_sum_prec,&
                                              pgf_sum_2c_gspace_1d,&
                                              pgf_sum_2c_gspace_3d,&
                                              pgf_sum_2c_rspace_1d,&
                                              pgf_sum_2c_rspace_3d,&
                                              pgf_sum_3c_1d
   USE eri_mme_types,                   ONLY: eri_mme_param
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_integrate'

   PUBLIC :: eri_mme_2c_integrate, eri_mme_3c_integrate, &
             eri_mme_2c_integrate_set, eri_mme_3c_integrate_set

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
!> \param normalize     calculate integrals w.r.t. normalized Hermite-Gaussians
!> \param exact_method  create reference values by exact method (no minimax approx.
!>                      and sum is converged to much higher precision)
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate(param, la_min, la_max, lb_min, lb_max, zeta, zetb, rab, &
                                   hab, o1, o2, G_count, R_count, normalize, exact_method)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      INTEGER, INTENT(IN)                                :: o1, o2
      INTEGER, INTENT(INOUT), OPTIONAL                   :: G_count, R_count
      LOGICAL, INTENT(IN), OPTIONAL                      :: normalize, exact_method

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_2c_integrate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, bx, by, bz, G_bound, handle, &
                                                            i_aw, i_xyz, ico, jco, k, l_max, la, &
                                                            lb, n_aw, R_bound1, R_bound2
      INTEGER, DIMENSION(3)                              :: gbounds_curr, la_xyz, lb_xyz, &
                                                            rbounds_curr
      LOGICAL                                            :: do_g_sum, exact, is_ortho, norm
      REAL(KIND=dp)                                      :: alpha_G, alpha_R, G_rad, G_res, Imm, &
                                                            inv_lgth, Ixyz, lgth, max_error, &
                                                            prefac, R_rad, R_res, vol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: S_G_1, S_G_2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: S_G, S_G_no, S_R_no
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: S_R_no_H
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, hmat, ht
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aw

      CALL timeset(routineN, handle)

      CPASSERT(param%is_valid)

      ! minimax coeffs
      n_aw = param%n_minimax
      aw => param%minimax_aw

      ! cell info
      h_inv = param%h_inv
      hmat = param%hmat
      vol = param%vol

      IF (PRESENT(normalize)) THEN
         norm = normalize
      ELSE
         norm = .FALSE.
      ENDIF

      l_max = la_max+lb_max

      ! resolution for Gaussian width
      G_res = 0.5_dp*param%G_min
      R_res = 0.5_dp*param%R_min

      IF (PRESENT(exact_method)) THEN
         exact = exact_method
      ELSE
         exact = .FALSE.
      ENDIF

      IF (exact) THEN
         is_ortho = .FALSE.
      ELSE
         is_ortho = param%is_ortho
      ENDIF

      IF (is_ortho) THEN ! orthorhombic cell
         ALLOCATE (S_G(0:l_max, 3, n_aw))
         S_G = 0.0_dp

         IF (param%debug) THEN
            ALLOCATE (S_G_1(0:l_max))
            ALLOCATE (S_G_2(0:l_max))
         ENDIF

         ! 1) precompute Ewald-like sum
         DO i_aw = 1, n_aw

            ! Gaussian width in G and R space
            alpha_G = aw(i_aw)+0.25_dp/zeta+0.25_dp/zetb
            alpha_R = 0.25_dp/alpha_G

            ! Radius for screening (in 1d)
            G_rad = exp_radius(l_max, alpha_G, lattice_sum_prec, 1.0_dp, epsin=G_res)
            R_rad = exp_radius(l_max, alpha_R, lattice_sum_prec, 1.0_dp, epsin=R_res)

            DO i_xyz = 1, 3
               lgth = ABS(hmat(i_xyz, i_xyz))
               inv_lgth = ABS(h_inv(i_xyz, i_xyz))
               G_bound = FLOOR(G_rad/(inv_lgth*twopi))
               R_bound1 = CEILING((rab(i_xyz)-R_rad)/lgth)
               R_bound2 = FLOOR((rab(i_xyz)+R_rad)/lgth)

               ! perform sum in R or G space. Choose the space in which less summands are required for convergence
               do_g_sum = (2*G_bound+1) < R_bound2-R_bound1

               IF (do_g_sum) THEN
                  CALL pgf_sum_2c_gspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_G, inv_lgth, &
                                            -G_bound, G_bound)
                  IF (PRESENT(G_count)) G_count = G_count+1
               ELSE
                  CALL pgf_sum_2c_rspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_R, lgth, &
                                            R_bound1, R_bound2)
                  IF (PRESENT(R_count)) R_count = R_count+1
               ENDIF

               IF (param%debug) THEN
                  ! check consistency of summation methods
                  CALL pgf_sum_2c_gspace_1d(S_G_1, -rab(i_xyz), alpha_G, inv_lgth, -G_bound, G_bound)
                  CALL pgf_sum_2c_rspace_1d(S_G_2, -rab(i_xyz), alpha_R, lgth, R_bound1, R_bound2)
                  max_error = MAXVAL(ABS(S_G_1-S_G_2)/(0.5_dp*(ABS(S_G_1)+ABS(S_G_2))+1.0_dp))

                  CPASSERT(max_error .LE. param%debug_delta)
               ENDIF
            ENDDO
         ENDDO

      ELSE ! general cell
         ! minimax approximation only used for R space evaluation (minimax approx. less efficient in G space)

         ht = twopi*TRANSPOSE(h_inv)

         alpha_G = 0.25_dp/zeta+0.25_dp/zetb
         alpha_R = 0.25_dp/alpha_G

         IF (exact) THEN
            G_rad = exp_radius(la_max+lb_max, alpha_G, 0.01*lattice_sum_prec, 1.0_dp, epsin=G_res)
         ELSE
            G_rad = exp_radius(la_max+lb_max, alpha_G, lattice_sum_prec, 1.0_dp, epsin=G_res)
         ENDIF
         R_rad = exp_radius(la_max+lb_max, alpha_R, lattice_sum_prec, 1.0_dp, epsin=R_res)
         DO k = 1, 3
            gbounds_curr(k) = FLOOR(G_rad/SQRT((ht(1, k)**2+ht(2, k)**2+ht(3, k)**2)))
            rbounds_curr(k) = CEILING(R_rad/SQRT(hmat(1, k)**2+hmat(2, k)**2+hmat(3, k)**2))
         ENDDO

         IF (exact) THEN
            do_g_sum = .TRUE.
         ELSE
            do_g_sum = PRODUCT(2*rbounds_curr+1) .GT. PRODUCT(2*gbounds_curr+1)
         ENDIF

         IF (do_g_sum) THEN
            ALLOCATE (S_G_no(0:l_max, 0:l_max, 0:l_max))
            IF (.NOT. exact .AND. PRESENT(G_count)) G_count = G_count+1
            CALL pgf_sum_2c_gspace_3d(S_G_no, -rab, alpha_G, h_inv, -gbounds_curr, gbounds_curr, vol)
         ELSE
            ALLOCATE (S_R_no_H(0:l_max, 0:l_max, 0:l_max, n_aw))
            ALLOCATE (S_R_no(0:l_max, 0:l_max, 0:l_max))
            S_R_no(:, :, :) = 0.0_dp
            IF (PRESENT(R_count)) R_count = R_count+1
            DO i_aw = 1, n_aw
               alpha_R = 1.0_dp/(4.0_dp*aw(i_aw)+1.0_dp/zeta+1.0_dp/zetb)
               R_rad = exp_radius(l_max, alpha_R, lattice_sum_prec, 1.0_dp, epsin=R_res)
               DO k = 1, 3
                  rbounds_curr(k) = CEILING(R_rad/SQRT(hmat(1, k)**2+hmat(2, k)**2+hmat(3, k)**2))
               ENDDO
               CALL pgf_sum_2c_rspace_3d(S_R_no_H(:, :, :, i_aw), -rab, alpha_R, hmat, -rbounds_curr, rbounds_curr)
               S_R_no(:, :, :) = S_R_no(:, :, :)+aw(n_aw+i_aw)*S_R_no_H(:, :, :, i_aw)
            ENDDO
         ENDIF
      ENDIF

      ! prefactor for integral values (unnormalized Hermite Gaussians)
      prefac = SQRT(1.0_dp/(zeta*zetb))

      ! 2) Assemble integral values from Ewald sums
      DO lb = lb_min, lb_max
      DO bx = 0, lb
      DO by = 0, lb-bx
         bz = lb-bx-by
         jco = coset(bx, by, bz)
         DO la = la_min, la_max
         DO ax = 0, la
         DO ay = 0, la-ax
            az = la-ax-ay
            la_xyz = [ax, ay, az]
            lb_xyz = [bx, by, bz]
            ico = coset(ax, ay, az)
            IF (is_ortho) THEN
               Imm = 0.0_dp
               DO i_aw = 1, n_aw
                  Ixyz = 1.0_dp
                  DO i_xyz = 1, 3
                     Ixyz = Ixyz*S_G(la_xyz(i_xyz)+lb_xyz(i_xyz), i_xyz, i_aw)*prefac
                  ENDDO
                  Imm = Imm+aw(n_aw+i_aw)*Ixyz
               ENDDO
            ELSE
               IF (do_g_sum) THEN
                  Imm = S_G_no(ax+bx, ay+by, az+bz)*prefac**3
               ELSE
                  Imm = S_R_no(ax+bx, ay+by, az+bz)*prefac**3
               ENDIF
            ENDIF
            IF ((la+lb .EQ. 0) .AND. (is_ortho .OR. .NOT. do_g_sum)) THEN
               Imm = Imm-SUM(aw(n_aw+1:2*n_aw))*prefac**3/vol ! subtracting G = 0 term
            ENDIF
            IF (.NOT. norm) THEN
               ! rescaling needed due to Hermite Gaussians (such that they can be contracted same way as Cartesian Gaussians)
               ! and factor of 4 pi**4 (-1)**lb
               hab(o1+ico, o2+jco) = Imm*4.0_dp*pi**4/((2.0_dp*zeta)**la*(-2.0_dp*zetb)**lb)
            ELSE
               ! same thing for normalized Hermite Gaussians
               hab(o1+ico, o2+jco) = Imm*4.0_dp*pi**4*(-1.0_dp)**lb*hermite_gauss_norm(zeta, la_xyz)* &
                                     hermite_gauss_norm(zetb, lb_xyz)
            ENDIF
         ENDDO ! la
         ENDDO
         ENDDO
      ENDDO ! lb
      ENDDO
      ENDDO

      CALL timestop(handle)
   END SUBROUTINE eri_mme_2c_integrate

! **************************************************************************************************
!> \brief Integrate set pair and contract with sphi matrix.
!> \param param ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param npgfb ...
!> \param zetb ...
!> \param ra ...
!> \param rb ...
!> \param hab ...
!> \param i1 ...
!> \param i2 ...
!> \param sphi_a ...
!> \param sgfa ...
!> \param nsgfa ...
!> \param sphi_b ...
!> \param sgfb ...
!> \param nsgfb ...
!> \param G_count ...
!> \param R_count ...
!> \note arrays required to be contiguous
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate_set(param, la_min, la_max, lb_min, lb_max, npgfa, zeta, npgfb, &
                                       zetb, ra, rb, hab, i1, i2, sphi_a, sgfa, nsgfa, sphi_b, sgfb, nsgfb, G_count, R_count)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, npgfa
      REAL(KIND=dp), DIMENSION(npgfa), INTENT(IN)        :: zeta
      INTEGER, INTENT(IN)                                :: npgfb
      REAL(KIND=dp), DIMENSION(npgfb), INTENT(IN)        :: zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: hab
      INTEGER, INTENT(IN)                                :: i1, i2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sphi_a
      INTEGER, INTENT(IN)                                :: sgfa, nsgfa
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sphi_b
      INTEGER, INTENT(IN)                                :: sgfb, nsgfb
      INTEGER, INTENT(INOUT), OPTIONAL                   :: G_count, R_count

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_2c_integrate_set', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipgf, jpgf, na, nb, ncoa, ncob
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: hab_herm

      ! note: tested only for one exponent per pair (npgfa = npgfb = 1)

      CALL timeset(routineN, handle)
      ncoa = npgfa*ncoset(la_max)
      ncob = npgfb*ncoset(lb_max)

      ALLOCATE (hab_herm(ncoa, ncob))

      hab_herm(:, :) = 0.0_dp

      DO ipgf = 1, npgfa
         na = (ipgf-1)*ncoset(la_max)
         DO jpgf = 1, npgfb
            nb = (jpgf-1)*ncoset(lb_max)
            CALL eri_mme_2c_integrate(param, &
                                      la_min, la_max, lb_min, lb_max, &
                                      zeta(ipgf), zetb(jpgf), ra-rb, hab_herm, na, nb, G_count, R_count)
         ENDDO
      ENDDO

      CALL ab_contract(hab(i1:i1+nsgfa-1, i2:i2+nsgfb-1), hab_herm, sphi_a(:, sgfa:), sphi_b(:, sgfb:), ncoa, ncob, nsgfa, nsgfb)

      CALL timestop(handle)
   END SUBROUTINE eri_mme_2c_integrate_set

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
!> \param transp ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_integrate(param, la_min, la_max, lb_min, lb_max, lc_min, lc_max, zeta, zetb, zetc, RA, RB, RC, &
                                   habc, o1, o2, o3, GG_count, GR_count, RR_count, transp)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, zetc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: RA, RB, RC
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN)                                :: o1, o2, o3
      INTEGER, INTENT(INOUT), OPTIONAL                   :: GG_count, GR_count, RR_count
      LOGICAL, INTENT(IN), OPTIONAL                      :: transp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_3c_integrate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, bx, by, bz, cx, cy, cz, &
                                                            G_bound_0, handle, i_aw, i_xyz, ico, &
                                                            jco, kco, la, lb, lc, lmax_0, method, &
                                                            n_aw, R_bound1_0, R_bound2_0
      INTEGER, DIMENSION(3)                              :: la_xyz, lb_xyz, lc_xyz
      LOGICAL                                            :: do_g_sum_0
      REAL(KIND=dp) :: alpha_G_0, alpha_R_0, G_rad_0, G_res, Imm, Ixyz, Ixyz_0, lgth, max_error, &
         max_result, min_result, prefac, R_0, R_rad_0, R_res, resc_a, resc_a_init, resc_b, &
         resc_b_init, resc_c, resc_c_init, vol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: S_G_0
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: S_G_1_1, S_G_1_2, S_G_1_3, S_G_2_1, &
                                                            S_G_2_2, S_G_2_3, S_G_3
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: S_G
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, hmat
      REAL(KIND=dp), DIMENSION(7)                        :: results
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aw

      CALL timeset(routineN, handle)

      CPASSERT(param%is_valid)

      ! minimax coeffs
      n_aw = param%n_minimax
      aw => param%minimax_aw

      ! cell info
      h_inv = param%h_inv
      hmat = param%hmat
      vol = param%vol

      ! prefactor for integral values (unnormalized Hermite Gaussians)
      prefac = (zeta*zetb*zetc)**(-0.5_dp)*pi**(11.0_dp/6.0_dp)*4.0_dp**(1.0_dp/3.0_dp)

      ALLOCATE (S_G(0:la_max, 0:lb_max, 0:lc_max, 3, n_aw))

      IF (param%debug) THEN
         ALLOCATE (S_G_1_1(0:la_max, 0:lb_max, 0:lc_max))
         ALLOCATE (S_G_1_2(0:la_max, 0:lb_max, 0:lc_max))
         ALLOCATE (S_G_1_3(0:la_max, 0:lb_max, 0:lc_max))
         ALLOCATE (S_G_2_1(0:la_max, 0:lb_max, 0:lc_max))
         ALLOCATE (S_G_2_2(0:la_max, 0:lb_max, 0:lc_max))
         ALLOCATE (S_G_2_3(0:la_max, 0:lb_max, 0:lc_max))
         ALLOCATE (S_G_3(0:la_max, 0:lb_max, 0:lc_max))
      ENDIF

      ALLOCATE (S_G_0(0:la_max+lb_max, 3))

      ! Preparations for G=0 component
      G_res = 0.5_dp*param%G_min
      R_res = 0.5_dp*param%R_min

      IF (param%is_ortho) THEN
         DO i_aw = 1, n_aw
            DO i_xyz = 1, 3

               lgth = ABS(hmat(i_xyz, i_xyz))

               CALL pgf_sum_3c_1d(S_G(:, :, :, i_xyz, i_aw), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                  zeta, zetb, zetc, aw(i_aw), lgth, method_out=method)

               SELECT CASE (method)
               CASE (1)
                  IF (PRESENT(GG_count)) GG_count = GG_count+1
               CASE (2)
                  IF (PRESENT(GR_count)) GR_count = GR_count+1
               CASE (3)
                  IF (PRESENT(RR_count)) RR_count = RR_count+1
               CASE DEFAULT
                  CPABORT("")
               END SELECT

               IF (param%debug) THEN
                  ! check consistency of summation methods
                  CALL pgf_sum_3c_1d(S_G_1_1(:, :, :), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                     zeta, zetb, zetc, aw(i_aw), lgth, method=1, order=1)
                  CALL pgf_sum_3c_1d(S_G_1_2(:, :, :), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                     zeta, zetb, zetc, aw(i_aw), lgth, method=1, order=2)
                  CALL pgf_sum_3c_1d(S_G_1_3(:, :, :), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                     zeta, zetb, zetc, aw(i_aw), lgth, method=1, order=3)

                  CALL pgf_sum_3c_1d(S_G_2_1(:, :, :), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                     zeta, zetb, zetc, aw(i_aw), lgth, method=2, order=1)
                  CALL pgf_sum_3c_1d(S_G_2_2(:, :, :), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                     zeta, zetb, zetc, aw(i_aw), lgth, method=2, order=2)
                  CALL pgf_sum_3c_1d(S_G_2_3(:, :, :), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                     zeta, zetb, zetc, aw(i_aw), lgth, method=2, order=3)
                  CALL pgf_sum_3c_1d(S_G_3(:, :, :), param%G_min, param%R_min, RA(i_xyz), RB(i_xyz), RC(i_xyz), &
                                     zeta, zetb, zetc, aw(i_aw), lgth, method=3)
                  max_error = 0.0_dp
                  DO la = 0, la_max
                     DO lb = 0, lb_max
                        DO lc = 0, lc_max
                           results = [S_G_1_1(la, lb, lc), S_G_1_2(la, lb, lc), S_G_1_3(la, lb, lc), &
                                      S_G_2_1(la, lb, lc), S_G_2_2(la, lb, lc), S_G_2_3(la, lb, lc), &
                                      S_G_3(la, lb, lc)]
                           max_result = MAXVAL(results)
                           min_result = MINVAL(results)
                           max_error = MAX(max_error, &
                                           (max_result-min_result)/(0.5_dp*(ABS(max_result)+ABS(min_result))+1.0_dp))
                        ENDDO
                     ENDDO
                  ENDDO

                  CPASSERT(max_error .LE. param%debug_delta)
               ENDIF

               ! compute G=0 component
               IF (i_aw .EQ. 1) THEN
                  alpha_G_0 = 0.25_dp/zetb+0.25_dp/zeta
                  alpha_R_0 = 0.25_dp/alpha_G_0
                  lmax_0 = la_max+lb_max
                  R_0 = RB(i_xyz)-RA(i_xyz)
                  G_rad_0 = exp_radius(lmax_0, alpha_G_0, lattice_sum_prec, 1.0_dp, epsin=G_res)
                  R_rad_0 = exp_radius(lmax_0, alpha_R_0, lattice_sum_prec, 1.0_dp, epsin=R_res)
                  G_bound_0 = FLOOR(G_rad_0*lgth/twopi)
                  R_bound1_0 = CEILING((-R_0-R_rad_0)/lgth)
                  R_bound2_0 = FLOOR((-R_0+R_rad_0)/lgth)
                  do_g_sum_0 = (2*G_bound_0+1) < R_bound2_0-R_bound1_0
                  IF (do_g_sum_0) THEN
                     CALL pgf_sum_2c_gspace_1d(S_G_0(:, i_xyz), R_0, alpha_G_0, 1.0_dp/lgth, -G_bound_0, G_bound_0)
                  ELSE
                     CALL pgf_sum_2c_rspace_1d(S_G_0(:, i_xyz), R_0, alpha_R_0, lgth, R_bound1_0, R_bound2_0)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         ! Initialization of rescaling factors due to Hermite Gaussians
         resc_a_init = (2.0_dp*zeta)**la_min
         resc_b_init = (2.0_dp*zetb)**lb_min
         resc_c_init = (2.0_dp*zetc)**lc_min

         ! assemble integral values
         resc_c = resc_c_init
         DO lc = lc_min, lc_max
            DO cx = 0, lc
            DO cy = 0, lc-cx
               cz = lc-cx-cy
               lc_xyz = [cx, cy, cz]
               kco = coset(cx, cy, cz)
               resc_b = resc_b_init
               DO lb = lb_min, lb_max
                  DO bx = 0, lb
                  DO by = 0, lb-bx
                     bz = lb-bx-by
                     lb_xyz = [bx, by, bz]
                     jco = coset(bx, by, bz)
                     resc_a = resc_a_init
                     DO la = la_min, la_max
                        DO ax = 0, la
                        DO ay = 0, la-ax
                           az = la-ax-ay
                           la_xyz = [ax, ay, az]
                           ico = coset(ax, ay, az)
                           Imm = 0.0_dp
                           DO i_aw = 1, n_aw
                              ! need to subtract G = 0 term only if lc == 0
                              Ixyz = 1.0_dp
                              DO i_xyz = 1, 3
                                 Ixyz = Ixyz*S_G(la_xyz(i_xyz), lb_xyz(i_xyz), lc_xyz(i_xyz), i_xyz, i_aw)*prefac
                              ENDDO
                              IF (lc .GT. 0) THEN
                                 Imm = Imm+aw(n_aw+i_aw)*(Ixyz)
                              ELSE
                                 Ixyz_0 = 1.0_dp
                                 DO i_xyz = 1, 3
                                    Ixyz_0 = Ixyz_0*S_G_0(la_xyz(i_xyz)+lb_xyz(i_xyz), i_xyz)/lgth*(-1.0_dp)**lb_xyz(i_xyz)*prefac
                                 ENDDO
                                 Imm = Imm+aw(n_aw+i_aw)*(Ixyz-Ixyz_0)
                              ENDIF

                              ! rescaling needed due to Hermite Gaussians
                              habc(o1+ico, o2+jco, o3+kco) = Imm/(resc_a*resc_b*resc_c)

                              ! transpose
                              IF (PRESENT(transp)) THEN
                                 IF (transp) THEN
                                    habc(o1+ico, o2+jco, o3+kco) = habc(o1+ico, o2+jco, o3+kco)*(-1.0_dp)**(la+lb+lc)
                                 ENDIF
                              ENDIF
                           ENDDO ! i_aw
                        ENDDO
                        ENDDO
                        resc_a = resc_a*(2.0_dp*zeta)
                     ENDDO ! la
                  ENDDO
                  ENDDO
                  resc_b = resc_b*(2.0_dp*zetb)
               ENDDO ! lb
            ENDDO
            ENDDO
            resc_c = resc_c*(2.0_dp*zetc)
         ENDDO ! lc
      ELSE
         CPABORT("general (non-ortho) cell not yet implemented")
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE eri_mme_3c_integrate

! **************************************************************************************************
!> \brief Integrate set triple and contract with sphi matrix
!> \param param ...
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param lc_min ...
!> \param lc_max ...
!> \param npgfa ...
!> \param npgfb ...
!> \param npgfc ...
!> \param zeta ...
!> \param zetb ...
!> \param zetc ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param habc ...
!> \param i1 ...
!> \param i2 ...
!> \param i3 ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sgfa ...
!> \param sgfb ...
!> \param sgfc ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param GG_count ...
!> \param GR_count ...
!> \param RR_count ...
!> \param transp ...
!> \note
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_integrate_set(param, &
                                       la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                                       npgfa, npgfb, npgfc, &
                                       zeta, zetb, zetc, &
                                       ra, rb, rc, &
                                       habc, i1, i2, i3, &
                                       sphi_a, sphi_b, sphi_c, &
                                       sgfa, sgfb, sgfc, &
                                       nsgfa, nsgfb, nsgfc, &
                                       GG_count, GR_count, RR_count, &
                                       transp)

      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max, npgfa, npgfb, npgfc
      REAL(KIND=dp), DIMENSION(npgfa), INTENT(IN)        :: zeta
      REAL(KIND=dp), DIMENSION(npgfb), INTENT(IN)        :: zetb
      REAL(KIND=dp), DIMENSION(npgfc), INTENT(IN)        :: zetc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rb, rc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: habc
      INTEGER, INTENT(IN)                                :: i1, i2, i3
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sphi_a, sphi_b, sphi_c
      INTEGER, INTENT(IN)                                :: sgfa, sgfb, sgfc, nsgfa, nsgfb, nsgfc
      INTEGER, INTENT(INOUT), OPTIONAL                   :: GG_count, GR_count, RR_count
      LOGICAL, INTENT(IN), OPTIONAL                      :: transp

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_3c_integrate_set', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipgf, jpgf, kpgf, na, nb, nc, &
                                                            ncoa, ncob, ncoc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ha_hb_hc

      CALL timeset(routineN, handle)

      ncoa = npgfa*ncoset(la_max)
      ncob = npgfb*ncoset(lb_max)
      ncoc = npgfc*ncoset(lc_max)

      ALLOCATE (ha_hb_hc(ncoa, ncob, ncoc)); ha_hb_hc(:, :, :) = 0.0_dp

      DO ipgf = 1, npgfa
         na = (ipgf-1)*ncoset(la_max)
         DO jpgf = 1, npgfb
            nb = (jpgf-1)*ncoset(lb_max)
            DO kpgf = 1, npgfc
               nc = (kpgf-1)*ncoset(lc_max)
               CALL eri_mme_3c_integrate(param, &
                                         la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                                         zeta(ipgf), zetb(jpgf), zetc(kpgf), ra, rb, rc, ha_hb_hc, na, nb, nc, &
                                         GG_count, GR_count, RR_count, &
                                         transp=transp)
            ENDDO
         ENDDO
      ENDDO

      CALL abc_contract(habc(i1:i1+nsgfa-1, i2:i2+nsgfb-1, i3:i3+nsgfc-1), ha_hb_hc, &
                        sphi_a(:, sgfa:), sphi_b(:, sgfb:), sphi_c(:, sgfc:), &
                        ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc)

      CALL timestop(handle)

   END SUBROUTINE eri_mme_3c_integrate_set

END MODULE eri_mme_integrate
