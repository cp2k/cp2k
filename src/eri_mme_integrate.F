!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Minimax-Ewald (MME) method for calculating electron repulsion integrals
!>        (ERI) of periodic systems using a Hermite Gaussian basis.
!>        The method relies on analytical Fourier transforms of Cartesian and
!>        Hermite Gaussian functions and Poisson summation formula to represent
!>        ERIs as a discrete sum over direct lattice vectors or reciprocal
!>        lattice vectors. The reciprocal space potential 1/G^2 is approximated
!>        by a linear combination of Gaussians employing minimax approximation.
!>        Implementation of 2-center ERIs for general cells but currently only
!>        efficient for orthorhombic cells.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_integrate

   USE ao_util,                         ONLY: exp_radius
   USE eri_mme_gaussian,                ONLY: hermite_gauss_norm
   USE eri_mme_lattice_summation,       ONLY: pgf_sum_gspace_1d,&
                                              pgf_sum_gspace_3d,&
                                              pgf_sum_rspace_1d,&
                                              pgf_sum_rspace_3d
   USE eri_mme_types,                   ONLY: dynamic_lattice_space,&
                                              eri_mme_param,&
                                              lattice_sum_prec,&
                                              rec_lattice_space
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_integrate'

   PUBLIC :: eri_mme_2c_integrate_low, &
             eri_mme_2c_integrate_set

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
!> \param normalize     calculate integrals w.r.t. normalized Hermite-Gaussians
!> \param exact_method  create reference values by exact method (no minimax approx.
!>                      and sum is converged to much higher precision)
!> \param lattice_space Can be optionally set to rec_lattice_space or dir_lattice_space
!>                      in order to only use one space representation (testing purposes).
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate_low(param, la_min, la_max, lb_min, lb_max, zeta, zetb, rab, &
                                       hab, o1, o2, normalize, exact_method, lattice_space)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: hab
      INTEGER, INTENT(IN)                                :: o1, o2
      LOGICAL, INTENT(IN), OPTIONAL                      :: normalize, exact_method
      INTEGER, INTENT(IN), OPTIONAL                      :: lattice_space

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_2c_integrate_low', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, bx, by, bz, G_bound, &
                                                            gr_space, i_aw, i_xyz, ico, jco, k, &
                                                            l_max, la, lb, n_aw, R_bound1, R_bound2
      INTEGER, DIMENSION(3)                              :: gbounds_curr, la_xyz, lb_xyz, &
                                                            rbounds_curr
      LOGICAL                                            :: do_g_sum, exact, is_ortho, norm
      REAL(KIND=dp)                                      :: alpha_G, alpha_R, G_rad, G_res, Imm, &
                                                            inv_lgth, Ixyz, lgth, prefac, &
                                                            prefac_screening, R_rad, R_res, vol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: S_G, S_G_no, S_R_no
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: S_R_no_H
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, hmat, ht
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aw

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

      ! prefactor for screening (normalized Hermite Gaussians)
      prefac_screening = eri_mme_max_prefactor_2c(zeta, zetb, la_max, la_min, lb_max, lb_min)

      ! prefactor for integral values (unnormalized Hermite Gaussians)
      prefac = SQRT(1.0_dp/(zeta*zetb))

      l_max = la_max+lb_max

      ! resolution for screening
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

      IF (PRESENT(lattice_space)) THEN
         CPASSERT(lattice_space .GE. 0 .AND. lattice_space .LE. 2)
         gr_space = lattice_space
      ELSE
         gr_space = dynamic_lattice_space
      ENDIF

      IF (is_ortho) THEN ! orthorhombic cell
         ALLOCATE (S_G(0:l_max, 3, n_aw))
         S_G = 0.0_dp

         ! 1) precompute Ewald-like sum
         DO i_aw = 1, n_aw

            ! Gaussian width in G and R space
            alpha_G = aw(i_aw)+0.25_dp/zeta+0.25_dp/zetb
            alpha_R = 0.25_dp/alpha_G

            ! Radius for screening (in 1d)
            G_rad = exp_radius(l_max, alpha_G, lattice_sum_prec, prefac_screening, epsin=G_res)
            R_rad = exp_radius(l_max, alpha_R, lattice_sum_prec, prefac_screening, epsin=R_res)

            ! G radius should not be larger than cutoff
            IF (G_rad .GT. SQRT(2.0*param%cutoff)) THEN
               G_rad = SQRT(2.0*param%cutoff)
            ENDIF

            DO i_xyz = 1, 3
               lgth = ABS(hmat(i_xyz, i_xyz))
               inv_lgth = ABS(h_inv(i_xyz, i_xyz))
               G_bound = FLOOR(G_rad/(inv_lgth*twopi))
               R_bound1 = CEILING((rab(i_xyz)-R_rad)/lgth)
               R_bound2 = FLOOR((rab(i_xyz)+R_rad)/lgth)

               ! perform sum in R or G space. Either choose the space in which less summands are required for convergence
               ! or use autocalibrated exponents to switch between spaces.
               IF (gr_space .EQ. dynamic_lattice_space) THEN
                  do_g_sum = (2*G_bound+1) < R_bound2-R_bound1
               ELSE
                  do_g_sum = gr_space .EQ. rec_lattice_space
               ENDIF

               IF (do_g_sum) THEN
                  CALL pgf_sum_gspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_G, inv_lgth, &
                                         -G_bound, G_bound)
                  param%G_count = param%G_count+1
               ELSE
                  CALL pgf_sum_rspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_R, lgth, &
                                         R_bound1, R_bound2)
                  param%R_count = param%R_count+1
               ENDIF
            ENDDO
         ENDDO

      ELSE ! general cell
         ! minimax approximation only used for R space evaluation (minimax approx. less efficient in G space)

         ht = twopi*TRANSPOSE(h_inv)

         alpha_G = 0.25_dp/zeta+0.25_dp/zetb
         alpha_R = 0.25_dp/alpha_G

         IF (exact) THEN
            G_rad = exp_radius(la_max+lb_max, alpha_G, 0.01*lattice_sum_prec, prefac_screening, epsin=G_res)
         ELSE
            G_rad = exp_radius(la_max+lb_max, alpha_G, lattice_sum_prec, prefac_screening, epsin=G_res)
         ENDIF
         R_rad = exp_radius(la_max+lb_max, alpha_R, lattice_sum_prec, prefac_screening, epsin=R_res)
         DO k = 1, 3
            gbounds_curr(k) = FLOOR(G_rad/SQRT((ht(1, k)**2+ht(2, k)**2+ht(3, k)**2)))
            rbounds_curr(k) = CEILING(R_rad/SQRT(hmat(1, k)**2+hmat(2, k)**2+hmat(3, k)**2))
         ENDDO

         IF (exact) THEN
            do_g_sum = .TRUE.
         ELSE IF (gr_space .EQ. dynamic_lattice_space) THEN
            do_g_sum = PRODUCT(2*rbounds_curr+1) .GT. PRODUCT(2*gbounds_curr+1)
         ELSE
            do_g_sum = gr_space .EQ. rec_lattice_space
         ENDIF

         IF (do_g_sum) THEN
            ALLOCATE (S_G_no(0:l_max, 0:l_max, 0:l_max))
            IF (.NOT. exact) param%G_count = param%G_count+1
            CALL pgf_sum_gspace_3d(S_G_no, -rab, alpha_G, h_inv, -gbounds_curr, gbounds_curr, vol)
         ELSE
            ALLOCATE (S_R_no_H(0:l_max, 0:l_max, 0:l_max, n_aw))
            ALLOCATE (S_R_no(0:l_max, 0:l_max, 0:l_max))
            S_R_no(:, :, :) = 0.0_dp
            param%R_count = param%R_count+1
            DO i_aw = 1, n_aw
               alpha_R = 1.0_dp/(4.0_dp*aw(i_aw)+1.0_dp/zeta+1.0_dp/zetb)
               R_rad = exp_radius(l_max, alpha_R, lattice_sum_prec, prefac_screening, epsin=R_res)
               DO k = 1, 3
                  rbounds_curr(k) = CEILING(R_rad/SQRT(hmat(1, k)**2+hmat(2, k)**2+hmat(3, k)**2))
               ENDDO
               CALL pgf_sum_rspace_3d(S_R_no_H(:, :, :, i_aw), -rab, alpha_R, hmat, -rbounds_curr, rbounds_curr)
               S_R_no(:, :, :) = S_R_no(:, :, :)+aw(n_aw+i_aw)*S_R_no_H(:, :, :, i_aw)
            ENDDO
         ENDIF
      ENDIF

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
               ! rescale to have same transformation for Hermite-Gaussians to Spherical Gaussians as
               ! for Cartesian Gaussians
               hab(o1+ico, o2+jco) = Imm*4.0_dp*pi**4/((2.0_dp*zeta)**la*(-2.0_dp*zetb)**lb)
            ELSE
               hab(o1+ico, o2+jco) = Imm*4.0_dp*pi**4*(-1.0_dp)**lb*hermite_gauss_norm(zeta, la_xyz)* &
                                     hermite_gauss_norm(zetb, lb_xyz)
            ENDIF
         ENDDO ! la
         ENDDO
         ENDDO
      ENDDO ! lb
      ENDDO
      ENDDO

   END SUBROUTINE eri_mme_2c_integrate_low

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
!> \note arrays required to be contiguous
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_integrate_set(param, la_min, la_max, lb_min, lb_max, npgfa, zeta, npgfb, &
                                       zetb, ra, rb, hab, i1, i2, sphi_a, sgfa, nsgfa, sphi_b, sgfb, nsgfb)
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

      INTEGER                                            :: ipgf, jpgf, na, nb, ncoa, ncob
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: hab_half, hab_herm

      ! note: tested only for one exponent per pair (npgfa = npgfb = 1)

      ncoa = npgfa*ncoset(la_max)
      ncob = npgfb*ncoset(lb_max)

      ALLOCATE (hab_herm(ncoa, ncob))

      hab_herm(:, :) = 0.0_dp

      DO ipgf = 1, npgfa
         na = (ipgf-1)*ncoset(la_max)
         DO jpgf = 1, npgfb
            nb = (jpgf-1)*ncoset(lb_max)
            CALL eri_mme_2c_integrate_low(param, &
                                          la_min, la_max, lb_min, lb_max, &
                                          zeta(ipgf), zetb(jpgf), ra-rb, hab_herm, na, nb)
         ENDDO
      ENDDO

      ALLOCATE (hab_half(nsgfa, ncob))
      hab_half(:, :) = 0.0_dp

      CALL dgemm("T", "N", nsgfa, ncob, ncoa, &
                 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                 hab_herm(1, 1), SIZE(hab_herm, 1), &
                 0.0_dp, hab_half(1, 1), SIZE(hab_half, 1))
      CALL dgemm("N", "N", nsgfa, nsgfb, ncob, &
                 1.0_dp, hab_half(1, 1), SIZE(hab_half, 1), &
                 sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                 0.0_dp, hab(i1, i2), SIZE(hab, 1))

   END SUBROUTINE eri_mme_2c_integrate_set

! **************************************************************************************************
!> \brief Maximum prefactor for 2c integrals w.r.t. to normalized Hermite Gaussians
!>        in 1d.
!> \param zeta ...
!> \param zetb ...
!> \param la_max ...
!> \param la_min ...
!> \param lb_max ...
!> \param lb_min ...
!> \retval prefac ...
! **************************************************************************************************
   FUNCTION eri_mme_max_prefactor_2c(zeta, zetb, la_max, la_min, lb_max, lb_min) RESULT(prefac)
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb
      INTEGER, INTENT(IN)                                :: la_max, la_min, lb_max, lb_min
      REAL(KIND=dp)                                      :: prefac

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_max_prefactor_2c', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ax, ay, az, bx, by, bz, la, lb
      INTEGER, DIMENSION(3)                              :: la_xyz, lb_xyz
      REAL(KIND=dp)                                      :: zet_fac

      zet_fac = 1.0_dp/SQRT(zeta*zetb)
      prefac = 0.0_dp

      DO lb = lb_min, lb_max
      DO bx = 0, lb
      DO by = 0, lb-bx
         bz = lb-bx-by
         DO la = la_min, la_max
         DO ax = 0, la
         DO ay = 0, la-ax
            az = la-ax-ay
            la_xyz = [ax, ay, az]
            lb_xyz = [bx, by, bz]
            prefac = MAX(prefac, zet_fac*hermite_gauss_norm(zeta, la_xyz)*hermite_gauss_norm(zetb, lb_xyz))
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
   END FUNCTION eri_mme_max_prefactor_2c

END MODULE eri_mme_integrate
