!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Low-level routines for Minimax-Ewald (MME) method for electron repulsion
!>        integrals. Don't include this routines directly, use the high level
!>        routines in eri_mme_methods.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! *****************************************************************************

MODULE eri_mme_methods_low

  USE ao_util,                         ONLY: exp_radius
  USE cell_types,                      ONLY: cell_type
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE eri_mme_types,                   ONLY: dynamic_lattice_space,&
                                             eri_mme_param,&
                                             lattice_sum_prec,&
                                             lmax,&
                                             min_Rc,&
                                             rec_lattice_space
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: gamma1,&
                                             gaussi,&
                                             pi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE minimax,                         ONLY: check_range,&
                                             get_minimax_coeff
  USE orbital_pointers,                ONLY: coset
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_methods_low'

  PUBLIC :: eri_mme_2c_integrate_low, &
            eri_mme_set_params_low

CONTAINS

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE eri_mme_2c_integrate_low(param, la_min, la_max, lb_min, lb_max, zeta, zetb, rab, &
                                   hab, o1, o2, normalize, exact_method, lattice_space)
    TYPE(eri_mme_param), INTENT(INOUT)       :: param
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, lb_max
    REAL(KIND=dp), INTENT(IN)                :: zeta, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: hab
    INTEGER, INTENT(IN)                      :: o1, o2
    LOGICAL, INTENT(IN), OPTIONAL            :: normalize, exact_method
    INTEGER, INTENT(IN), OPTIONAL            :: lattice_space

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_2c_integrate_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, bx, by, bz, &
                                                G_bound, gr_space, i_aw, &
                                                i_xyz, ico, jco, k, l_max, &
                                                la, lb, n_aw, R_bound1, &
                                                R_bound2
    INTEGER, DIMENSION(3)                    :: gbounds_curr, la_xyz, lb_xyz, &
                                                rbounds_curr
    LOGICAL                                  :: do_g_sum, exact, is_ortho, &
                                                norm
    REAL(KIND=dp) :: alpha_G, alpha_R, G_rad, G_res, Imm, inv_lgth, Ixyz, &
      lgth, prefac, prefac_screening, R_rad, R_res, vol
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: S_G, S_G_no, S_R_no
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: S_R_no_H
    REAL(KIND=dp), DIMENSION(3, 3)           :: h_inv, hmat, ht
    REAL(KIND=dp), DIMENSION(:), POINTER     :: aw

    CPASSERT(param%is_valid)

    ! minimax coeffs
    n_aw = param%n_minimax
    aw => param%minimax_aw

    ! cell info
    h_inv = param%h_inv
    hmat = param%hmat
    vol = param%vol

    IF(PRESENT(normalize)) THEN
      norm = normalize
    ELSE
      norm = .FALSE.
    ENDIF

    ! prefactor for screening (normalized Hermite Gaussians)
    prefac_screening = eri_mme_max_prefactor(zeta, zetb, la_max, la_min, lb_max, lb_min)

    ! prefactor for integral values (unnormalized Hermite Gaussians)
    prefac = SQRT(1.0_dp/(zeta*zetb))

    l_max = la_max + lb_max

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
      ALLOCATE(S_G(0:l_max, 3, n_aw))
      S_G = 0.0_dp

      ! 1) precompute Ewald-like sum
      DO i_aw = 1 , n_aw

        ! Gaussian width in G and R space
        alpha_G = aw(i_aw) + 0.25_dp/zeta + 0.25_dp/zetb
        alpha_R = 0.25_dp/alpha_G

        ! Radius for screening (in 1d)
        G_rad = exp_radius(l_max, alpha_G, lattice_sum_prec, prefac_screening, epsin = G_res)
        R_rad = exp_radius(l_max, alpha_R, lattice_sum_prec, prefac_screening, epsin = R_res)

        ! G radius should not be larger than cutoff
        IF(G_rad .GT. SQRT(2.0*param%cutoff)) THEN
          G_rad = SQRT(2.0*param%cutoff)
        ENDIF

        DO i_xyz = 1, 3
          lgth = ABS(hmat(i_xyz, i_xyz))
          inv_lgth = ABS(h_inv(i_xyz, i_xyz))
          G_bound = FLOOR(G_rad/(inv_lgth*twopi))
          R_bound1 = CEILING((rab(i_xyz) - R_rad)/lgth)
          R_bound2 = FLOOR((rab(i_xyz) + R_rad)/lgth)

          ! perform sum in R or G space. Either choose the space in which less summands are required for convergence
          ! or use autocalibrated exponents to switch between spaces.
          IF (gr_space .EQ. dynamic_lattice_space) THEN
            IF(param%do_calib_alpha) THEN
              do_g_sum = alpha_G .GE. param%alpha_c(l_max, i_xyz)
            ELSE
              do_g_sum = (2*G_bound + 1) < R_bound2 - R_bound1
            ENDIF
          ELSE
            do_g_sum = gr_space .EQ. rec_lattice_space
          ENDIF

          IF(do_g_sum) THEN
           CALL pgf_sum_gspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_G, inv_lgth, &
                                  -G_bound, G_bound)
           param%G_count = param%G_count + 1
          ELSE
           CALL pgf_sum_rspace_1d(S_G(:, i_xyz, i_aw), -rab(i_xyz), alpha_R, lgth, &
                                  R_bound1, R_bound2)
           param%R_count = param%R_count + 1
          ENDIF
        ENDDO
      ENDDO

    ELSE ! general cell
       ! minimax approximation only used for R space evaluation (minimax approx. less efficient in G space)

      ht = twopi*TRANSPOSE(h_inv)

      alpha_G = 0.25_dp/zeta + 0.25_dp/zetb
      alpha_R = 0.25_dp/alpha_G

      IF (exact) THEN
        G_rad = exp_radius(la_max + lb_max, alpha_G, 0.01*lattice_sum_prec, prefac_screening, epsin = G_res)
      ELSE
        G_rad = exp_radius(la_max + lb_max, alpha_G, lattice_sum_prec, prefac_screening, epsin = G_res)
      ENDIF
      R_rad = exp_radius(la_max + lb_max, alpha_R, lattice_sum_prec, prefac_screening, epsin = R_res)
      DO k = 1, 3
         gbounds_curr(k) = FLOOR(G_rad/SQRT((ht(1, k)**2 + ht(2, k)**2 + ht(3, k)**2)))
         rbounds_curr(k) = CEILING(R_rad/SQRT(hmat(1, k)**2 + hmat(2, k)**2 + hmat(3, k)**2))
      ENDDO

      IF (exact) THEN
        do_g_sum = .TRUE.
      ELSE IF (gr_space .EQ. dynamic_lattice_space) THEN
        do_g_sum = PRODUCT(2*rbounds_curr + 1) .GT. PRODUCT(2*gbounds_curr + 1)
      ELSE
        do_g_sum = gr_space .EQ. rec_lattice_space
      ENDIF

      IF(do_g_sum) THEN
        ALLOCATE(S_G_no(0:l_max, 0:l_max, 0:l_max))
        IF(.NOT. exact) param%G_count = param%G_count + 1
        CALL pgf_sum_gspace_3d(S_G_no, -rab, alpha_G, h_inv, -gbounds_curr, gbounds_curr, vol)
      ELSE
        ALLOCATE(S_R_no_H(0:l_max, 0:l_max, 0:l_max, n_aw))
        ALLOCATE(S_R_no(0:l_max, 0:l_max, 0:l_max))
        S_R_no(:, :, :) = 0.0_dp
        param%R_count = param%R_count + 1
        DO i_aw = 1, n_aw
          alpha_R = 1.0_dp/(4.0_dp*aw(i_aw) + 1.0_dp/zeta + 1.0_dp/zetb)
          R_rad = exp_radius(l_max, alpha_R, lattice_sum_prec, prefac_screening, epsin = R_res)
          DO k = 1, 3
            rbounds_curr(k) = CEILING(R_rad/SQRT(hmat(1, k)**2 + hmat(2, k)**2 + hmat(3, k)**2))
          ENDDO
          CALL pgf_sum_rspace_3d(S_R_no_H(:, :, :, i_aw), -rab, alpha_R, hmat, -rbounds_curr, rbounds_curr)
          S_R_no(:, :, :) = S_R_no(:, :, :) + aw(n_aw + i_aw)*S_R_no_H(:, :, :, i_aw)
        ENDDO
      ENDIF
    ENDIF

    ! 2) Assemble integral values from Ewald sums
    DO lb = lb_min, lb_max
    DO bx = 0, lb
    DO by = 0, lb - bx
      bz = lb - bx - by
      jco = coset(bx, by, bz)
      DO la = la_min, la_max
      DO ax = 0, la
      DO ay = 0, la - ax
        az = la - ax - ay
        la_xyz = [ax, ay, az]
        lb_xyz = [bx, by, bz]
        ico = coset(ax, ay, az)
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
          IF(do_g_sum) THEN
            Imm = S_G_no(ax + bx, ay + by, az + bz)*prefac**3
          ELSE
            Imm = S_R_no(ax + bx, ay + by, az + bz)*prefac**3
          ENDIF
        ENDIF
        IF ((la + lb .EQ. 0) .AND. (is_ortho .OR. .NOT. do_g_sum)) THEN
          Imm = Imm - SUM(aw(n_aw + 1:2*n_aw))*prefac**3/vol ! subtracting G = 0 term
        ENDIF
        IF(.NOT. norm) THEN
          ! rescale to have same transformation for Hermite-Gaussians to Spherical Gaussians as 
          ! for Cartesian Gaussians
          hab(o1 + ico, o2 + jco) = Imm*4.0_dp*pi**4/((-2.0_dp*zeta)**la*(2.0_dp*zetb)**lb)
        ELSE
          hab(o1 + ico, o2 + jco) = Imm*4.0_dp*pi**4*hermite_gauss_norm(zeta, la_xyz)*hermite_gauss_norm(zetb, lb_xyz)
        ENDIF
      ENDDO ! la
      ENDDO
      ENDDO
    ENDDO ! lb
    ENDDO
    ENDDO

  END SUBROUTINE eri_mme_2c_integrate_low

! *****************************************************************************
!> \brief Maximum prefactor for integrals w.r.t. to normalized Hermite Gaussians 
!>        in 1d.
!> \param zeta ...
!> \param zetb ...
!> \param la_max ...
!> \param la_min ...
!> \param lb_max ...
!> \param lb_min ...
!> \retval prefac ...
! *****************************************************************************
  FUNCTION eri_mme_max_prefactor(zeta, zetb, la_max, la_min, lb_max, lb_min) RESULT(prefac)
    REAL(KIND=dp), INTENT(IN)                :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_max, la_min, lb_max, lb_min
    REAL(KIND=dp)                            :: prefac

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_max_prefactor', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, bx, by, bz, la, lb
    INTEGER, DIMENSION(3)                    :: la_xyz, lb_xyz
    REAL(KIND=dp)                            :: zet_fac

    zet_fac = 1.0_dp/SQRT(zeta*zetb)
    prefac = 0.0_dp

    DO lb = lb_min, lb_max
    DO bx = 0, lb
    DO by = 0, lb - bx
      bz = lb - bx - by
      DO la = la_min, la_max
      DO ax = 0, la
      DO ay = 0, la - ax
        az = la - ax - ay
        la_xyz = [ax, ay, az]
        lb_xyz = [bx, by, bz]
        prefac = MAX(prefac,zet_fac*hermite_gauss_norm(zeta, la_xyz)*hermite_gauss_norm(zetb, lb_xyz))
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    ENDDO
  END FUNCTION eri_mme_max_prefactor

! *****************************************************************************
!> \brief Compute Ewald-like sum in G space in 1 dimension
!>        S_G(l, alpha) = i^l*inv_lgth*sum_G( C(l, alpha, G) exp(iGR) ), with
!>                        C(l, alpha, r) = r^l exp(-alpha*r^2),
!>        dG = inv_lgth*twopi and G = -G_bound*dG, (-G_bound + 1)*dG, ..., G_bound*dG
!>             for all l < = l_max.
!> \param S_G ...
!> \param R ...
!> \param alpha ...
!> \param inv_lgth ...
!> \param G_min ...
!> \param G_c ...
!> \note  S_G is real.
! *****************************************************************************
  SUBROUTINE pgf_sum_gspace_1d(S_G, R, alpha, inv_lgth, G_min, G_c)
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(INOUT)                          :: S_G
    REAL(KIND=dp), INTENT(IN)                :: R, alpha, inv_lgth
    INTEGER, INTENT(IN)                      :: G_min, G_c

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pgf_sum_gspace_1d', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: exp_tot
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: i_pow_l, S_G_c
    INTEGER                                  :: k, l, l_max
    REAL(KIND=dp)                            :: dG, G, G_pow_l

    dG = inv_lgth*twopi
    l_max = UBOUND(S_G, 1)

    ALLOCATE(i_pow_l(0:l_max))
    i_pow_l(0) = 1.0_dp
    DO l = 1, l_max
      i_pow_l(l) = i_pow_l(l - 1)*gaussi
    ENDDO

    ALLOCATE(S_G_c(0:l_max))
    S_G_c(:) = 0.0_dp
    DO k = G_min, G_c
      G = k*dG
      exp_tot = EXP(-alpha*G**2)*EXP(gaussi*G*R)
      G_pow_l = 1.0_dp
      DO l = 0, l_max
        !G_pow_l = G**l
        S_G_c(l) = S_G_c(l) + G_pow_l*exp_tot
        G_pow_l = G_pow_l*G
      ENDDO
    ENDDO

    S_G(:) = REAL(S_G_c(:)*i_pow_l(:))*inv_lgth
  END SUBROUTINE pgf_sum_gspace_1d

! *****************************************************************************
!> \brief Compute Ewald-like sum in R space in 1 dimension
!>        S_R(l, alpha) = SQRT(alpha/pi) sum_R'( H(l, alpha, R-R') ),
!>        with H(l, alpha, R) = (d/dR)^l exp(-alpha*R^2),
!>        dR = lgth and R' = -R_min*dR, (-R_min + 1)*dR, ..., R_max*dR,
!>        for all l < = l_max.
!> \param S_R ...
!> \param R ...
!> \param alpha ...
!> \param lgth ...
!> \param R_min ...
!> \param R_max ...
!> \note  result is equivalent to pgf_sum_gspace_1d with
!>              S_R(l, alpha) = S_G(l, 1/(4*alpha))
! *****************************************************************************
  SUBROUTINE pgf_sum_rspace_1d(S_R, R, alpha, lgth, R_min, R_max)
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(INOUT)                          :: S_R
    REAL(KIND=dp), INTENT(IN)                :: R, alpha, lgth
    INTEGER, INTENT(IN)                      :: R_min, R_max

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pgf_sum_rspace_1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, l, l_max
    REAL(KIND=dp)                            :: dR, exp_tot, R_pow_l, Rp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: h_to_c

    dR = lgth
    l_max = UBOUND(S_R, 1)

    ! 1) compute sum over C(l, alpha, R - R') instead of H(l, alpha, R - R')

    S_R(:) = 0.0_dp
    Rp = R + R_min*dR
    DO k = R_min, R_max
      Rp = R + k*dR
      exp_tot = EXP(-alpha*Rp**2)
      R_pow_l = 1.0_dp
      DO l = 0, l_max
        S_R(l) = S_R(l) + R_pow_l*exp_tot
        R_pow_l = R_pow_l*Rp
      ENDDO
    ENDDO

    ! 2) transform C --> H

    ALLOCATE(h_to_c(0:l_max, 0:l_max))
    CALL create_hermite_to_cartesian(alpha, h_to_c)
    S_R = MATMUL(h_to_c, S_R)*SQRT(alpha/pi)
  END SUBROUTINE pgf_sum_rspace_1d

! *****************************************************************************
!> \brief As pgf_sum_gspace_1d but 3d sum required for non-orthorhombic cells
!> \param S_G ...
!> \param R ...
!> \param alpha ...
!> \param h_inv ...
!> \param G_min ...
!> \param G_c ...
!> \param vol ...
!> \note  MMME Method is not very efficient for non-orthorhombic cells
! *****************************************************************************
  SUBROUTINE pgf_sum_gspace_3d(S_G, R, alpha, h_inv, G_min, G_c, vol)
    REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
      INTENT(INOUT)                          :: S_G
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: R
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: h_inv
    INTEGER, DIMENSION(3), INTENT(IN)        :: G_min, G_c
    REAL(KIND=dp), INTENT(IN)                :: vol

    COMPLEX(KIND=dp)                         :: exp_tot
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Ig
    INTEGER                                  :: gx, gy, gz, k, l, l_max, lx, &
                                                ly, lz
    INTEGER, DIMENSION(3)                    :: g_xyz, l_xyz
    REAL(KIND=dp)                            :: G_sq, Ik
    REAL(KIND=dp), DIMENSION(3)              :: G
    REAL(KIND=dp), DIMENSION(3, 3)           :: ht

    l_max = UBOUND(S_G, 1)
    ht = twopi*TRANSPOSE(h_inv)
    ALLOCATE(Ig(0:l_max, 0:l_max, 0:l_max))
    Ig(:, :, :) = 0.0_dp
    DO gx = G_min(1), G_c(1)
    DO gy = G_min(2), G_c(2)
    DO gz = G_min(3), G_c(3)
      g_xyz = [gx, gy, gz]
      G = MATMUL(ht, g_xyz)
      G_sq = G(1)**2 + G(2)**2 + G(3)**2
      exp_tot = EXP(-alpha*G_sq)*EXP(-gaussi*DOT_PRODUCT(G, -R))
      IF(.NOT.(gx .EQ. 0 .AND. gy .EQ. 0 .AND. gz .EQ. 0)) THEN
        DO l = 0, l_max
        DO lx = 0, l
        DO ly = 0, l - lx
          lz = l - lx - ly
          l_xyz = [lx, ly, lz]
          Ik = 1.0_dp
          DO k = 1, 3
            Ik = Ik*G(k)**(l_xyz(k))
          ENDDO
          Ig(lx, ly, lz) = Ig(lx, ly, lz) + exp_tot/G_sq*Ik*gaussi**l/vol
        ENDDO
        ENDDO
        ENDDO
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    S_G(:, :, :) = REAL(Ig(:, :, :), KIND = dp)
  END SUBROUTINE pgf_sum_gspace_3d

! *****************************************************************************
!> \brief As pgf_sum_rspace_1d but 3d sum required for non-orthorhombic cells
!> \param S_R ...
!> \param R ...
!> \param alpha ...
!> \param hmat ...
!> \param R_min ...
!> \param R_max ...
!> \note  MMME Method is not very efficient for non-orthorhombic cells
! *****************************************************************************
  SUBROUTINE pgf_sum_rspace_3d(S_R, R, alpha, hmat, R_min, R_max)
    REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
      INTENT(INOUT)                          :: S_R
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: R
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    INTEGER, DIMENSION(3), INTENT(IN)        :: R_min, R_max

    INTEGER                                  :: k, l, l_max, ll, llx, lly, &
                                                llz, lx, ly, lz, sx, sy, sz
    INTEGER, DIMENSION(3)                    :: l_xyz, s_xyz
    REAL(KIND=dp)                            :: exp_tot, Ik, R_sq
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: h_to_c
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: Is, S_R_C
    REAL(KIND=dp), DIMENSION(3)              :: Rp

    l_max = UBOUND(S_R, 1)
    ALLOCATE(S_R_C(0:l_max, 0:l_max, 0:l_max))
    ALLOCATE(Is(0:l_max, 0:l_max, 0:l_max))

    S_R(:, :, :) = 0.0_dp
    S_R_C(:, :, :) = 0.0_dp
    Is(:, :, :) = 0.0_dp

    DO sx = R_min(1), R_max(1)
    DO sy = R_min(2), R_max(2)
    DO sz = R_min(3), R_max(3)
      s_xyz = [sx, sy, sz]
      Rp = MATMUL(hmat, s_xyz)
      R_sq = (Rp(1) + R(1))**2 + (Rp(2) + R(2))**2 + (Rp(3) + R(3))**2
      exp_tot = EXP(-alpha*R_sq)
      DO l = 0, l_max
      DO lx = 0, l
      DO ly = 0, l - lx
        lz = l - lx - ly
        l_xyz = [lx, ly, lz]
        Ik = 1.0_dp
        DO k = 1, 3
          Ik = Ik*(Rp(k) + R(k))**(l_xyz(k))
        ENDDO
        Is(lx, ly, lz) = Is(lx, ly, lz) + exp_tot*Ik
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    ENDDO
    S_R_C(:, :, :) = Is(:, :, :)

    ALLOCATE(h_to_c(0:l_max, 0:l_max))
    CALL create_hermite_to_cartesian(alpha, h_to_c)
    DO l = 0, l_max
    DO lx = 0, l
    DO ly = 0, l - lx
      lz = l - lx - ly
      DO ll = 0, l_max
      DO llx = 0, ll
      DO lly = 0, ll - llx
        llz = ll - llx - lly
        S_R(lx, ly, lz) = S_R(lx, ly, lz) + &
                        h_to_c(lx, llx)*h_to_c(ly, lly)*h_to_c(lz, llz)*&
                        S_R_C(llx, lly, llz)
      ENDDO
      ENDDO
      ENDDO
    ENDDO
    ENDDO
    ENDDO
    S_R(:, :, :) = S_R(:, :, :)*(alpha/pi)**1.5_dp

 END SUBROUTINE pgf_sum_rspace_3d

! *****************************************************************************
!> \brief Create matrix to transform between cartesian and hermite gaussian
!>        basis functions.
!> \param zet    exponent
!> \param h_to_c transformation matrix with dimensions (0:l_max, 0:l_max)
!> \note  is idempotent, so transformation is the same
!>        in both directions.
! *****************************************************************************
  SUBROUTINE create_hermite_to_cartesian(zet, h_to_c)
    REAL(KIND=dp), INTENT(IN)                :: zet
    REAL(KIND=dp), DIMENSION(0:, 0:), &
      INTENT(INOUT)                          :: h_to_c

    INTEGER                                  :: k, l, l_max

    l_max = UBOUND(h_to_c, 1)

    CPASSERT((UBOUND(h_to_c, 2) .EQ. l_max))

    h_to_c(:, :) = 0.0_dp
    h_to_c(0, 0) = 1.0_dp
    DO l = 0, l_max - 1
      DO k = 0, l + 1
        IF(k .EQ. 0) THEN
                       h_to_c(l + 1, k) = (k + 1)*h_to_c(l, k + 1)
        ELSE
          IF(k .LT. l) h_to_c(l + 1, k) = (k + 1)*h_to_c(l, k + 1) - 2.0_dp*zet*h_to_c(l, k - 1)
          IF(k .GE. l) h_to_c(l + 1, k) =                          - 2.0_dp*zet*h_to_c(l, k - 1)
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE create_hermite_to_cartesian

! *****************************************************************************
!> \brief Norm of 1d Hermite-Gauss functions
!> \param zet ...
!> \param l ...
!> \retval norm ...
! *****************************************************************************
  FUNCTION hermite_gauss_norm(zet, l) RESULT(norm)
    REAL(KIND=dp), INTENT(IN)                :: zet
    INTEGER, DIMENSION(3), INTENT(IN)        :: l
    REAL(KIND=dp)                            :: norm

    norm = 1.0_dp/SQRT((2.0_dp*zet)**(SUM(l) - 1.5_dp)*(gamma1(l(1))*gamma1(l(2))*gamma1(l(3))))

  END FUNCTION hermite_gauss_norm

! *****************************************************************************
!> \brief Set parameters for MME method with manual specification of basis parameters.
!>        Takes care of cutoff calibration if requested.
!> \param param ...
!> \param cell ...
!> \param zet_err_minimax Exponent used to estimate error of minimax approximation.
!> \param zet_err_cutoff  Exponents used to estimate error of finite cutoff.
!> \param l_err_cutoff    Total ang. mom. quantum numbers l to be combined with exponents in
!>                        zet_err_cutoff.
!> \param lmax            Maximum total angular momentum quantum number
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE eri_mme_set_params_low(param, cell, zet_err_minimax, zet_err_cutoff, l_err_cutoff, lmax, para_env)
    TYPE(eri_mme_param), INTENT(INOUT)       :: param
    TYPE(cell_type), INTENT(IN)              :: cell
    REAL(KIND=dp), INTENT(IN)                :: zet_err_minimax
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zet_err_cutoff
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: l_err_cutoff
    INTEGER, INTENT(IN)                      :: lmax
    TYPE(cp_para_env_type), INTENT(IN), &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_mme_set_params_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, l, l_mm, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: l_c
    LOGICAL                                  :: calib_failed, min_cutoff, &
                                                s_only
    REAL(KIND=dp)                            :: C_mm, cutoff, E_mm, err_c, &
                                                err_mm
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zet_c
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN, handle)

    logger => cp_get_default_logger()
    unit_nr = cp_print_key_unit_nr(logger, param%mme_section, "ERI_MME_INFO", &
                                   extension = ".eri_mme")

    s_only = lmax .EQ. 0
    param%lmax = lmax

    CPASSERT(ASSOCIATED(param%mme_section))
    CPASSERT(PRESENT(l_err_cutoff) .OR. s_only)

    IF (PRESENT(l_err_cutoff)) THEN
      CPASSERT(SIZE(zet_err_cutoff) .EQ. SIZE(l_err_cutoff))
    ENDIF

    ! l values for minimax error estimate (l_mm) and for cutoff error estimate (l_c)
    l_mm = MERGE(0, 1, s_only)
    ALLOCATE(l_c(SIZE(zet_err_cutoff)), zet_c(SIZE(zet_err_cutoff)))
    IF(s_only) THEN
      l_c(:) = 0
    ELSE
      l_c(:) = l_err_cutoff(:)
    ENDIF

    ! cell info
    param%hmat = cell%hmat
    param%h_inv = cell%h_inv
    param%vol = cell%deth
    param%is_ortho = cell%orthorhombic

    ! Minimum lattice vectors
    param%G_min = G_abs_min(param%h_inv)
    param%R_min = R_abs_min(param%hmat)

    ! cutoff calibration not yet implemented for general cell
    IF (.NOT.param%is_ortho) THEN
      param%do_calib_alpha = .FALSE.
      param%do_calib_cutoff = .FALSE.
    ENDIF

    ! Cutoff calibration and error estimate for orthorhombic cell
    IF (param%do_calib_cutoff) THEN
      CALL calibrate_cutoff(param, param%hmat, param%h_inv, param%G_min, param%vol, &
                            zet_err_minimax, l_mm, zet_err_cutoff, l_c, param%n_minimax, &
                            param%cutoff_min, param%cutoff_max, param%cutoff_eps, &
                            param%cutoff_delta, cutoff, err_mm, err_c, &
                            C_mm, para_env)

      param%cutoff = cutoff
      param%error_tot = err_mm + err_c
    ELSE IF (param%is_ortho) THEN
      CALL cutoff_minimax_error(param%cutoff, param%hmat, param%h_inv, param%vol, param%G_min, &
                                zet_err_minimax, l_mm, zet_err_cutoff, l_c, param%n_minimax, &
                                param%minimax_aw, err_mm, err_c, C_mm, min_cutoff, para_env)
      param%error_tot = err_mm + err_c
    ENDIF

    CALL get_minimax_coeff_v_gspace(param%n_minimax, param%cutoff, param%G_min, param%minimax_aw, &
                                    E_mm, min_cutoff)

    IF (unit_nr>0) WRITE(unit_nr, '(T2, A, T71, F10.1)') "ERI_MME| Cutoff for ERIs [a.u.]:", param%cutoff
    IF(param%is_ortho) THEN
      IF (unit_nr>0) THEN
        WRITE(unit_nr, '(T2, A)') "ERI_MME| Estimated absolute error for normalized Hermite-Gaussian basis"
        WRITE(unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Minimax error:", err_mm
        WRITE(unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Cutoff error:", err_c
        WRITE(unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME|   Total error (minimax + cutoff):", err_mm + err_c
        IF(param%print_calib) &
          WRITE(unit_nr, '(T2, A, T68, F13.10)') "ERI_MME|   Minimax scaling constant in AM-GM estimate:", C_mm
      ENDIF
    ENDIF

    ! Calibration for switching between G and R space representation
    calib_failed = .FALSE.
    IF(param%do_calib_alpha) THEN
      CALL calibrate_gr_switch(param%alpha_c, param%hmat, param%h_inv, &
                               lmax, calib_failed, para_env)
      param%do_calib_alpha = .NOT. calib_failed
    ENDIF
    IF(calib_failed) THEN
      CALL cp_warn(__LOCATION__,"GR-switch calibration FAILED. "//&
           "Default switching between G-space and R-space will be applied")
    ENDIF

    IF(param%do_calib_alpha .AND. param%print_calib) THEN
      IF(unit_nr>0) THEN
        WRITE(unit_nr, '(/T2, A)') "ERI_MME| Calibrated exponents for switching between R- and G-space evaluation" 
        WRITE(unit_nr, '(T11, A)') "l, alpha(x), alpha(y), alpha(z)"
        DO l = 0, lmax
          WRITE(unit_nr, '(T11, I2, 1X, ES9.2, 1X, ES9.2, 1X, ES9.2)') l, &
                        param%alpha_c(l, 1), param%alpha_c(l, 2), param%alpha_c(l, 3)
        ENDDO
      ENDIF
    ENDIF

    IF(unit_nr>0) WRITE(unit_nr, *)

    param%G_count = 0
    param%R_count = 0

    param%is_valid = .TRUE.

    CALL cp_print_key_finished_output(unit_nr, logger, param%mme_section, "ERI_MME_INFO")
    CALL timestop(handle)
  END SUBROUTINE eri_mme_set_params_low

! *****************************************************************************
!> \brief Find minimum length of R vectors, for a general (not necessarily
!>        orthorhombic) cell.
!> \param hmat ...
!> \retval R_m ...
! *****************************************************************************
  FUNCTION R_abs_min(hmat) RESULT(R_m)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    REAL(KIND=dp)                            :: R_m

    INTEGER                                  :: sx, sy, sz
    INTEGER, DIMENSION(3)                    :: sxyz
    REAL(KIND=dp)                            :: R_sq
    REAL(KIND=dp), DIMENSION(3)              :: R

    R_m = 0.0_dp

    DO sx = -1, 1
    DO sy = -1, 1
    DO sz = -1, 1
      IF(.NOT.(sx == 0 .AND. sy == 0 .AND. sz == 0)) THEN
        sxyz = [sx, sy, sz]
        R = MATMUL(hmat, sxyz)
        R_sq = R(1)**2 + R(2)**2 + R(3)**2
        IF(R_sq < R_m .OR. R_m < EPSILON(R_m)) R_m = R_sq
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    R_m = SQRT(R_m)

  END FUNCTION R_abs_min

! *****************************************************************************
!> \brief Find minimum length of G vectors, for a general (not necessarily
!>        orthorhombic) cell.
!> \param h_inv ...
!> \retval G_m ...
! *****************************************************************************
  FUNCTION G_abs_min(h_inv) RESULT(G_m)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: h_inv
    REAL(KIND=dp)                            :: G_m

    INTEGER                                  :: gx, gy, gz
    INTEGER, DIMENSION(3)                    :: gxyz
    REAL(KIND=dp)                            :: G_sq
    REAL(KIND=dp), DIMENSION(3)              :: G
    REAL(KIND=dp), DIMENSION(3, 3)           :: H

    H = twopi*TRANSPOSE(h_inv)
    G_m = 0.0_dp

    DO gx = -1, 1
    DO gy = -1, 1
    DO gz = -1, 1
      IF(.NOT.(gx == 0 .AND. gy == 0 .AND. gz == 0)) THEN
        gxyz = [gx, gy, gz]
        G = MATMUL(H, gxyz)
        G_sq = G(1)**2 + G(2)**2 + G(3)**2
        IF(G_sq < G_m .OR. G_m < EPSILON(G_m)) G_m = G_sq
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    G_m = SQRT(G_m)

  END FUNCTION G_abs_min

! *****************************************************************************
!> \brief Find optimal cutoff minimizing errors due to minimax approximation and
!>        due to finite cutoff using bisection on the difference of the errors
!> \param param ...
!> \param hmat ...
!> \param h_inv ...
!> \param G_min ...
!> \param vol ...
!> \param zet_mm    Exponent to estimate minimax error
!> \param l_mm      Total ang. mom. quantum number to estimate minimax error
!> \param zet_c     Max. exponents to estimate cutoff error
!> \param l_c       Max. total ang. mom. quantum numbers to estimate cutoff error
!> \param n_minimax Number of terms in minimax approximation
!> \param cutoff_l  Initial guess of lower bound for cutoff
!> \param cutoff_r  Initial guess of upper bound for cutoff
!> \param tol       Tolerance (cutoff precision)
!> \param delta     to modify initial guess interval
!> \param cutoff    Best cutoff
!> \param err_mm    Minimax error
!> \param err_c     Cutoff error
!> \param C_mm      Scaling constant to generalize AM–GM upper bound estimate to
!>                  minimax approx.
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE calibrate_cutoff(param, hmat, h_inv, G_min, vol, zet_mm, l_mm, zet_c, l_c, &
                              n_minimax, cutoff_l, cutoff_r, tol, delta, &
                              cutoff, err_mm, err_c, C_mm, para_env)
    TYPE(eri_mme_param), INTENT(INOUT)       :: param
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat, h_inv
    REAL(KIND=dp), INTENT(IN)                :: G_min
    REAL(KIND=dp)                            :: vol
    REAL(KIND=dp), INTENT(IN)                :: zet_mm
    INTEGER, INTENT(IN)                      :: l_mm
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zet_c
    INTEGER, DIMENSION(:), INTENT(IN)        :: l_c
    INTEGER, INTENT(IN)                      :: n_minimax
    REAL(KIND=dp), INTENT(IN)                :: cutoff_l, cutoff_r, tol, delta
    REAL(KIND=dp), INTENT(OUT)               :: cutoff, err_mm, err_c, C_mm
    TYPE(cp_para_env_type), INTENT(IN), &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calibrate_cutoff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iter1, iter2, max_iter, &
                                                unit_nr
    LOGICAL                                  :: do_print, min_cutoff, &
                                                use_min_cutoff, valid_initial
    REAL(KIND=dp)                            :: cutoff_mid, delta_c_mid, &
                                                delta_mm_mid
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: minimax_aw
    REAL(KIND=dp), DIMENSION(2)              :: cutoff_lr, delta_c, delta_mm
    TYPE(cp_logger_type), POINTER            :: logger

    logger => cp_get_default_logger()
    unit_nr = cp_print_key_unit_nr(logger, param%mme_section, "ERI_MME_INFO", &
                                   extension = ".eri_mme")
    CPASSERT(SIZE(l_c) .EQ. SIZE(zet_c))

    IF (param%print_calib) THEN
      IF(unit_nr>0) THEN
        WRITE(unit_nr, '(/T2, A)') "ERI_MME| Basis set parameters for estimating minimax error"
        WRITE(unit_nr, '(T2, A, T67, ES12.2, 1X, I1)') "ERI_MME|   exp, l:", zet_mm, l_mm
        WRITE(unit_nr, '(T2, A)') "ERI_MME| Basis set parameters for estimating cutoff error"
        DO i = 1, SIZE(zet_c)
          WRITE(unit_nr, '(T2, A, T67, ES12.2, 1X, I1)') "ERI_MME|   exp, l:", zet_c(i), l_c(i)
        ENDDO
      ENDIF
    ENDIF

    max_iter = 100

    IF((cutoff_r - cutoff_l)/(0.5_dp*(cutoff_r + cutoff_l)) .LE. tol)&
       CALL cp_abort(__LOCATION__, "difference of boundaries for cutoff "//&
            "(MAX - MIN) must be greater than cutoff precision.")

    IF((delta .GE. 1.0_dp) .OR. (delta.LE.0.0_dp))&
       CPABORT("relative delta to modify initial cutoff interval (DELTA) must be in (0, 1)")

    cutoff_lr(1) = cutoff_l
    cutoff_lr(2) = cutoff_r

    ALLOCATE(minimax_aw(2*n_minimax))

    do_print = (unit_nr>0) .AND. (param%print_calib)
    IF(do_print) THEN
      WRITE(unit_nr, '(/T2, A)') "ERI_MME| Calibrating cutoff by bisecting error(minimax) - error(cutoff)"
      WRITE(unit_nr, '(T2, A, T72, ES9.2)') "ERI_MME| Rel. cutoff precision", tol
      WRITE(unit_nr, '(T2, A, T77, F4.1)') "ERI_MME| Rel. cutoff delta to modify initial interval", delta
    ENDIF

    ! 1) find valid initial values for bisection
    DO iter1 = 1, max_iter + 1
      IF(iter1.GT.max_iter)&
         CALL cp_abort(__LOCATION__, "Maximum number of iterations in bisection to determine initial "//&
              "cutoff interval has been exceeded.")

      cutoff_lr(1) = MAX(cutoff_lr(1), 0.5_dp*G_min**2)
      use_min_cutoff = .FALSE. ! Flag to stop bisection in case minimum cutoff (limited by minimax
                               ! approx.) is hit

      DO i = 1, 2
        CALL cutoff_minimax_error(cutoff_lr(i), hmat, h_inv, vol, G_min, zet_mm, l_mm, zet_c, l_c, &
                                  n_minimax, minimax_aw, delta_mm(i), delta_c(i), C_mm, min_cutoff, para_env)
        IF (min_cutoff) use_min_cutoff = .TRUE.
      ENDDO

      valid_initial = .TRUE.
      IF((delta_mm(1) - delta_c(1)).GT.0) THEN
        cutoff_lr(1) = cutoff_lr(1)*(1.0_dp - ABS(delta))
        valid_initial = .FALSE.
      ELSE
        use_min_cutoff = .FALSE.
      ENDIF
      IF((delta_mm(2) - delta_c(2)).LT.0) THEN
        cutoff_lr(2) = cutoff_lr(2)*(1.0_dp + ABS(delta))
        valid_initial = .FALSE.
      ENDIF

      IF (valid_initial .OR. use_min_cutoff) EXIT
    ENDDO

    ! 2) bisection to find cutoff s.t. err_minimax(cutoff) - err_cutoff(cutoff) = 0
    IF (use_min_cutoff) THEN
      cutoff = cutoff_lr(1)
      CALL cutoff_minimax_error(cutoff, hmat, h_inv, vol, G_min, zet_mm, l_mm, zet_c, l_c, &
                                n_minimax, minimax_aw, err_mm, err_c, C_mm, min_cutoff, para_env)
      iter2 = 1
    ELSE
      IF(do_print) WRITE(unit_nr, '(/T2, A)') &
                "ERI_MME| Step, cutoff (min, max, mid), err(minimax), err(cutoff), err diff"

      DO iter2 = 1, max_iter + 1
        IF(iter2.GT.max_iter)&
           CPABORT("Maximum number of iterations in bisection to determine cutoff has been exceeded")

        cutoff_mid = 0.5_dp*(cutoff_lr(1) + cutoff_lr(2))
        CALL cutoff_minimax_error(cutoff_mid, hmat, h_inv, vol, G_min, zet_mm, l_mm, zet_c, l_c, &
                                  n_minimax, minimax_aw, delta_mm_mid, delta_c_mid, C_mm, min_cutoff, para_env)
        IF(do_print) WRITE(unit_nr, '(T11, I2, F11.1, F11.1, F11.1, 3X, ES9.2, 3X, ES9.2, 3X, ES9.2)')&
                     iter2, cutoff_lr(1), cutoff_lr(2), cutoff_mid, &
                     delta_mm_mid, delta_c_mid, delta_mm_mid - delta_c_mid

        IF ((cutoff_lr(2) - cutoff_lr(1))/cutoff_mid .LT. tol) EXIT
        IF(delta_mm_mid - delta_c_mid .GT. 0) THEN
          cutoff_lr(2) = cutoff_mid
          delta_mm(2) = delta_mm_mid
          delta_c(2) = delta_c_mid
        ELSE
          cutoff_lr(1) = cutoff_mid
          delta_mm(1) = delta_mm_mid
          delta_c(1) = delta_c_mid
        ENDIF
      ENDDO
      err_mm = delta_mm_mid
      err_c = delta_c_mid
      cutoff = cutoff_mid
    ENDIF

    IF(do_print) THEN
      WRITE(unit_nr, '(/T2, A)') "ERI_MME| Cutoff calibration number of steps:"
      WRITE(unit_nr, '(T2, A, T79, I2)') "ERI_MME|   Steps for initial interval", iter1 - 1
      WRITE(unit_nr, '(T2, A, T79, I2/)') "ERI_MME|   Bisection iteration steps", iter2 - 1
    ENDIF

    CALL cp_print_key_finished_output(unit_nr, logger, param%mme_section, "ERI_MME_INFO")

  END SUBROUTINE calibrate_cutoff

! *****************************************************************************
!> \brief Compute upper bounds for the errors of 2-center ERI's (P|P) due
!>        to minimax approximation and due to finite cutoff, where P is a
!>        normalized Hermite Gaussian.
!> \param cutoff ...
!> \param hmat ...
!> \param h_inv ...
!> \param vol ...
!> \param G_min ...
!> \param zet_mm     Exponent of P to estimate minimax error
!> \param l_mm       total ang. mom. quantum number of P to estimate minimax error
!> \param zet_ctff   Max. exponents of P to estimate cutoff error
!> \param l_ctff     Max. total ang. mom. quantum numbers of P to estimate cutoff error
!> \param n_minimax  Number of terms in minimax approximation
!> \param minimax_aw Minimax coefficients
!> \param err_mm     Minimax error
!> \param err_ctff   Cutoff error
!> \param C_mm       Scaling constant to generalize AM–GM upper bound estimate to
!>                   minimax approx.
!> \param min_cutoff Whether cutoff (in) is less than minimum cutoff (limited
!>                   by minimax approx.). In this case cutoff (out) is minimum
!>                   cutoff.
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE cutoff_minimax_error(cutoff, hmat, h_inv, vol, G_min, zet_mm, l_mm, zet_ctff, l_ctff, &
                                  n_minimax, minimax_aw, err_mm, err_ctff, C_mm, min_cutoff, para_env)
    REAL(KIND=dp), INTENT(INOUT)             :: cutoff
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat, h_inv
    REAL(KIND=dp), INTENT(IN)                :: vol, G_min, zet_mm
    INTEGER, INTENT(IN)                      :: l_mm
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: zet_ctff
    INTEGER, DIMENSION(:), INTENT(IN)        :: l_ctff
    INTEGER, INTENT(IN)                      :: n_minimax
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: minimax_aw
    REAL(KIND=dp), INTENT(OUT)               :: err_mm, err_ctff, C_mm
    LOGICAL, INTENT(OUT)                     :: min_cutoff
    TYPE(cp_para_env_type), INTENT(IN), &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cutoff_minimax_error', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_aw, i_xyz, iG, iter, &
                                                max_iter, n, nG
    REAL(KIND=dp) :: C, dG, E_mm, eps_zet, err0, err1, err_c, err_d, G, G_1, &
      G_c, gr, prod_mm_k, zet_a, zet_b, zet_c, zet_d, zet_div
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zet_ctff_tmp
    REAL(KIND=dp), DIMENSION(2)              :: err_ctff_curr, err_ctff_prev

! parameters for finding exponent maximizing cutoff error

    eps_zet = 1.0E-05_dp ! tolerance for exponent
    zet_div = 2.0_dp     ! sampling constant for finding initial values of exponents
    max_iter = 100       ! maximum number of iterations in golden section search

    CPASSERT(SIZE(l_ctff) .EQ. SIZE(zet_ctff))

    ALLOCATE(zet_ctff_tmp(SIZE(zet_ctff)))
    zet_ctff_tmp(:) = zet_ctff(:)

    CALL get_minimax_coeff_v_gspace(n_minimax, cutoff, G_min, minimax_aw(:), E_mm, min_cutoff)
    G_c = SQRT(2.0*cutoff)

    ! 1) Minimax error, simple analytical formula
    prod_mm_k = 1.0_dp
    DO i_xyz = 1, 3
      prod_mm_k = prod_mm_k*(ABS(hmat(i_xyz, i_xyz))/twopi + &
                  MERGE(SQRT(2.0_dp/(zet_mm*pi))*EXP(-1.0_dp), 0.0_dp, l_mm .GT. 0))
    ENDDO
    err_mm = 32*pi**4/vol*E_mm*prod_mm_k

    ! 2) Cutoff error, estimating G > G_c part of Ewald sum by using C/3 * 1/(Gx^2*Gy^2*Gz^2)^1/3 as an
    !                  upper bound for 1/G^2 (AM-GM inequality) and its minimax approximation (factor C).
    !                  Note: usually, minimax approx. falls off faster than 1/G**2, so C should be approximately 1.
    !                  The error is calculated for all l up to l_max and golden section search algorithm is
    !                  applied to find the exponent that maximizes cutoff error.
    G_1 = SQRT(1.0_dp/(3.0_dp*MINVAL(minimax_aw(1:n_minimax))))

    C_mm = 0.0_dp
    IF(G_1 .GT. G_c) THEN
      dG = 1.0E-3
      nG = CEILING((G_1 - G_c)/dG)
      G = G_c
      DO iG = 1, nG
        G = MIN(G, G_c)
        C = 0.0_dp
        DO i_aw = 1, n_minimax
          C = C + 3.0_dp*minimax_aw(n_minimax + i_aw)*EXP(-3.0_dp*minimax_aw(i_aw)*G**2)*G**2
        ENDDO
        C_mm = MAX(C, C_mm)
        G = G + dG
      ENDDO
    ELSE
      DO i_aw = 1, n_minimax
        C_mm = C_mm + 3.0_dp*minimax_aw(n_minimax + i_aw)*EXP(-3.0_dp*minimax_aw(i_aw)*G_c**2)*G_c**2
      ENDDO
    ENDIF
    C = MAX(1.0_dp, C_mm)

    err_ctff_prev(:) = 0.0_dp
    gr = 0.5_dp*(SQRT(5.0_dp) - 1.0_dp) ! golden ratio
    DO n = 1, SIZE(l_ctff)
      ! Find valid starting values for golden section search
      DO iter = 1, max_iter + 1
        err_ctff_curr = 0.0_dp
        IF(iter.GT.max_iter)&
           CALL cp_abort(__LOCATION__, "Maximum number of iterations for finding "//&
                "exponent maximizing cutoff error has been exceeded.")

        CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_ctff(n), zet_ctff_tmp(n), C, err_ctff_curr(n), para_env)
        IF (err_ctff_prev(n) .GE. err_ctff_curr(n)) THEN
          zet_a = zet_ctff_tmp(n)
          zet_b = MIN(zet_ctff_tmp(n)*zet_div**2, zet_ctff(n))
          EXIT
        ELSE
          err_ctff_prev(n) = err_ctff_curr(n)
        ENDIF
        zet_ctff_tmp(n) = zet_ctff_tmp(n)/zet_div
      ENDDO

      ! Golden section search
      zet_c = zet_b - gr*(zet_b - zet_a)
      zet_d = zet_a + gr*(zet_b - zet_a)
      DO iter = 1, max_iter + 1
        IF(ABS(zet_c - zet_d) .LT. eps_zet*(zet_a + zet_b)) THEN
          CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_ctff(n), zet_a, C, err0, para_env)
          CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_ctff(n), zet_b, C, err1, para_env)
          err_ctff_curr(n) = MAX(err0, err1)
          EXIT
        ENDIF
        CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_ctff(n), zet_c, C, err_c, para_env)
        CALL cutoff_error_fixed_exp(cutoff, h_inv, G_min, l_ctff(n), zet_d, C, err_d, para_env)
        IF (err_c .GT. err_d) THEN
          zet_b = zet_d
          zet_d = zet_c
          zet_c = zet_b - gr*(zet_b - zet_a)
        ELSE
          zet_a = zet_c
          zet_c = zet_d
          zet_d = zet_a + gr*(zet_b - zet_a)
        ENDIF
      ENDDO
    ENDDO
    err_ctff = MAXVAL(err_ctff_curr)

  END SUBROUTINE cutoff_minimax_error

! *****************************************************************************
!> \brief Calculate cutoff error estimate by using C_mm/3 * 1/(Gx^2*Gy^2*Gz^2)^1/3
!>        as an upper bound for 1/G^2 (and its minimax approximation) for |G| > G_c.
!>        Error is referring to a basis function P with fixed exponent zet_c and
!>        max. angular momentum l_c.
!> \param cutoff ...
!> \param h_inv ...
!> \param G_min ...
!> \param l_c ...
!> \param zet_c ...
!> \param C_mm ...
!> \param err_c ...
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE cutoff_error_fixed_exp(cutoff,  h_inv, G_min, l_c, zet_c, C_mm, err_c, para_env)
    REAL(KIND=dp), INTENT(INOUT)             :: cutoff
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: h_inv
    REAL(KIND=dp), INTENT(IN)                :: G_min
    INTEGER, INTENT(IN)                      :: l_c
    REAL(KIND=dp), INTENT(IN)                :: zet_c, C_mm
    REAL(KIND=dp), INTENT(OUT)               :: err_c
    TYPE(cp_para_env_type), INTENT(IN), &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cutoff_error_fixed_exp', &
      routineP = moduleN//':'//routineN

    INTEGER :: ax, ay, az, G_l, G_u, Gl_first, Gl_last, Gu_first, Gu_last, &
      i_xyz, l, my_p, n_Gl, n_Gl_left, n_Gl_p, n_Gu, n_Gu_left, n_Gu_p, n_p
    REAL(KIND=dp)                            :: alpha_G, eps_G, err_c_l, G_c, &
                                                G_rad, G_res, inv_lgth, &
                                                prefactor, sum_G_diff
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: S_G_l, S_G_u

    G_c = SQRT(2.0_dp*cutoff)
    eps_G = TINY(eps_G) ! sum up to machine precision
    G_res = 0.5_dp*G_min ! resolution for screening

    err_c = 0.0_dp
    alpha_G = 1.0_dp/(2.0_dp*zet_c)
    prefactor = 1.0_dp/zet_c

    ALLOCATE(S_G_l(0:2*l_c, 3))
    ALLOCATE(S_G_u(0:2*l_c, 3))

    G_rad = exp_radius(2*l_c, alpha_G, eps_G, prefactor, epsin = G_res)

    ! Parallelization of sum over G vectors
    my_p = para_env%mepos ! mpi rank
    n_p = para_env%num_pe ! total number of processes

    DO i_xyz = 1, 3
      inv_lgth = ABS(h_inv(i_xyz, i_xyz))

      G_l = FLOOR(G_c/(inv_lgth*twopi))
      G_u = FLOOR(G_rad/(inv_lgth*twopi))

      IF (G_u .LT. G_l) G_u = G_l

      ! Serial code:
      ! !Sum |G| <= G_c
      ! CALL pgf_sum_gspace_1d_deltal(S_G_l(:, i_xyz), alpha_G, inv_lgth, -G_l, G_l, &
      !                               2.0_dp/3.0_dp, prefactor)
      ! !Sum |G| > G_c
      ! CALL pgf_sum_gspace_1d_deltal(S_G_u(:, i_xyz), alpha_G, inv_lgth, G_l + 1, G_u, &
      !                               2.0_dp/3.0_dp, prefactor)

      ! Parallel code:
      n_Gu = MAX((G_u - G_l), 0)
      n_Gl = 2*G_l + 1
      n_Gu_p = n_Gu/n_p
      n_Gl_p = n_Gl/n_p
      n_Gu_left = MOD(n_Gu, n_p)
      n_Gl_left = MOD(n_Gl, n_p)

      IF (my_p .LT. n_Gu_left) THEN
        Gu_first = G_l + 1 + (n_Gu_p + 1)*my_p
        Gu_last  = G_l + 1 + (n_Gu_p + 1)*(my_p + 1) - 1
      ELSE
        Gu_first = G_l + 1 + n_Gu_left + n_Gu_p*my_p
        Gu_last =  G_l + 1 + n_Gu_left + n_Gu_p*(my_p + 1) - 1
      ENDIF

      IF (my_p .LT. n_Gl_left) THEN
        Gl_first = -G_l + (n_Gl_p + 1)*my_p
        Gl_last  = -G_l + (n_Gl_p + 1)*(my_p + 1) - 1
      ELSE
        Gl_first = -G_l + n_Gl_left + n_Gl_p*my_p
        Gl_last =  -G_l + n_Gl_left + n_Gl_p*(my_p + 1) - 1
      ENDIF

      ! Sum |G| <= G_c
      CALL pgf_sum_gspace_1d_deltal(S_G_l(:, i_xyz), alpha_G, inv_lgth, Gl_first, Gl_last, &
                                    2.0_dp/3.0_dp, prefactor)

      ! Sum |G| > G_c
      CALL pgf_sum_gspace_1d_deltal(S_G_u(:, i_xyz), alpha_G, inv_lgth, Gu_first, Gu_last, &
                                    2.0_dp/3.0_dp, prefactor)
    ENDDO

    CALL mp_sum(S_G_l, para_env%group)
    CALL mp_sum(S_G_u, para_env%group)

    S_G_u = S_G_u * 2.0_dp ! to include negative values of G

    DO l = 0, l_c
    DO ax = 0, l
    DO ay = 0, l - ax
      az = l - ax - ay

      ! Compute prod_k (S_G_l(l_k,k) + S_G_u(l_k,k)) - prod_k (S_G_l(l_k,k)) with k in {x, y, z}
      ! Note: term by term multiplication to avoid subtraction for numerical stability
      sum_G_diff = S_G_u(2*ax, 1)*S_G_u(2*ay, 2)*S_G_u(2*az, 3) + &
                   S_G_u(2*ax, 1)*S_G_u(2*ay, 2)*S_G_l(2*az, 3) + &
                   S_G_u(2*ax, 1)*S_G_l(2*ay, 2)*S_G_u(2*az, 3) + &
                   S_G_l(2*ax, 1)*S_G_u(2*ay, 2)*S_G_u(2*az, 3) + &
                   S_G_u(2*ax, 1)*S_G_l(2*ay, 2)*S_G_l(2*az, 3) + &
                   S_G_l(2*ax, 1)*S_G_u(2*ay, 2)*S_G_l(2*az, 3) + &
                   S_G_l(2*ax, 1)*S_G_l(2*ay, 2)*S_G_u(2*az, 3)

      err_c_l = 4.0_dp*pi**4*hermite_gauss_norm(zet_c, [ax, ay, az])**2 * &
                C_mm/3.0_dp*sum_G_diff

      err_c = MAX(err_c, err_c_l)
    ENDDO
    ENDDO
    ENDDO

    DEALLOCATE(S_G_u, S_G_l)

  END SUBROUTINE cutoff_error_fixed_exp

! *****************************************************************************
!> \brief Get minimax coefficient a_i and w_i for approximating
!>        1/G^2 by sum_i w_i exp(-a_i G^2)
!> \param n_minimax   Number of minimax terms
!> \param cutoff      Plane Wave cutoff
!> \param G_min       Minimum absolute value of G
!> \param minimax_aw  Minimax coefficients a_i, w_i
!> \param err_minimax Maximum error MAX (|1/G^2-\sum_i w_i exp(-a_i G^2)|)
!> \param min_cutoff  Whether cutoff (in) is less than minimum cutoff (limited
!>                    by minimax approx.). In this case cutoff (out) is minimum
!>                    cutoff.
! *****************************************************************************
  SUBROUTINE get_minimax_coeff_v_gspace(n_minimax, cutoff, G_min, minimax_aw, err_minimax, min_cutoff)
    INTEGER, INTENT(IN)                      :: n_minimax
    REAL(KIND=dp), INTENT(INOUT)             :: cutoff
    REAL(KIND=dp), INTENT(IN)                :: G_min
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: minimax_aw
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: err_minimax
    LOGICAL, INTENT(OUT)                     :: min_cutoff

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_minimax_coeff_v_gspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_aw, ierr, iG, nG
    REAL(KIND=dp)                            :: dG, G, G_max, minimax_Rc, V_mm

    dG = 1.0E-3 ! Resolution in G to determine error of minimax approximation
    min_cutoff = .FALSE.

    ! Note: G_c = SQRT(2*cutoff) cutoff in 1 cartesian direction
    ! G_max = SQRT(3*G_c**2) maximum absolute value of G vector
    ! Minimax approx. needs to be valid in range [G_min, G_max]

    ! 1) compute minimax coefficients
    G_max = SQRT(3.0_dp*2.0_dp*cutoff)
    minimax_Rc = (G_max/G_min)**2
    IF (minimax_Rc.LE.min_Rc) THEN
      min_cutoff = .TRUE.
      minimax_Rc = min_Rc
      cutoff = min_Rc*G_min**2/(3.0_dp*2.0_dp)
    ENDIF

    G_max = SQRT(3.0_dp*2.0_dp*cutoff)
    CPASSERT(G_max .GT. G_min)

    CALL check_range(n_minimax, minimax_Rc, ierr)
    CPASSERT(ierr == 0)
    CALL get_minimax_coeff(n_minimax, minimax_Rc, minimax_aw, ierr)
    CPASSERT(ierr == 0)
    minimax_aw = minimax_aw/G_min**2

    CPASSERT(ALL(minimax_aw(:) .GT. 0))

    ! 2) sample minimax error
    IF (PRESENT(err_minimax)) THEN
      err_minimax = 0.0_dp
      nG = CEILING((G_max - G_min)/dG)
      G = G_min
      DO iG = 1, nG
        G = MIN(G, G_max)
        V_mm = 0.0_dp
        DO i_aw = 1, n_minimax
          V_mm = V_mm + minimax_aw(n_minimax + i_aw)*EXP(-minimax_aw(i_aw)*G**2)
        ENDDO
        err_minimax = MAX(err_minimax, ABS(V_mm - 1.0_dp/G**2))
        G = G + dG
      ENDDO
    ENDIF

  END SUBROUTINE get_minimax_coeff_v_gspace

! *****************************************************************************
!> \brief Find compound exponent at which G and R space representation
!>        need the same execution time (for each cartesian component of l).
!>        This works only for orthorhombic cells.
!> \param alpha_c ...
!> \param hmat ...
!> \param h_inv ...
!> \param lmax ...
!> \param failed ...
!> \param para_env ...
!> \note  should not be needed anymore except for testing purposes since
!>        switching between R and G space representation based on number of
!>        summands works well enough.
! *****************************************************************************
  SUBROUTINE calibrate_gr_switch(alpha_c, hmat, h_inv, lmax, failed, para_env)
    REAL(KIND=dp), DIMENSION(0:, :), &
      INTENT(INOUT)                          :: alpha_c
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat, h_inv
    INTEGER, INTENT(IN)                      :: lmax
    LOGICAL, INTENT(OUT)                     :: failed
    TYPE(cp_para_env_type), INTENT(IN), &
      POINTER                                :: para_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calibrate_gr_switch', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: G_bound, i, i_xyz, iter, l, &
                                                max_iter, n_points, n_start, &
                                                R_bound
    REAL(KIND=dp)                            :: alpha, alpha_eps, dt_G, dt_R, &
                                                G_rad, G_res, inv_lgth, lgth, &
                                                R_rad, R_res, t_0, t_1, t_2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: alpha_start, S_G, S_R
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: R
    REAL(KIND=dp), DIMENSION(2)              :: alpha_lr

    failed = .FALSE.
    alpha_eps = 1.0E-2 ! tolerance
    n_points = 10000 ! number of values for R (atomic distances) per process
    max_iter = 100

    ! define set of exponents alpha from which to select starting values
    n_start = 13
    ALLOCATE(alpha_start(n_start))
    alpha_start(:) = [1.0E-6_dp, 1.0E-5_dp, 1.0E-4_dp, 1.0E-3_dp, 1.0E-2_dp, 1.0E-1_dp, 1.0E0_dp, &
                      1.0E1_dp, 1.0E2_dp, 1.0E3_dp, 1.0E4_dp, 1.0E5_dp, 1.0E6_dp]

    CPASSERT(SIZE(alpha_c, 2) .EQ. 3)
    alpha_c = 0.0_dp

    ALLOCATE(R(n_points, 3)) ! n_points values of R are tested for each x, y, z
    DO i_xyz = 1, 3
      R(:, i_xyz) = (/ (i*ABS(hmat(i_xyz, i_xyz))/n_points, i = 1, n_points) /)
    ENDDO

    ! Resolution for screening
    G_res = 0.5_dp*G_abs_min(h_inv)
    R_res = 0.5_dp*R_abs_min(hmat)

    outerloop: DO i_xyz = 1, 3
      lgth = ABS(hmat(i_xyz, i_xyz))
      inv_lgth = ABS(h_inv(i_xyz, i_xyz))
      DO l = 0, lmax
        ALLOCATE(S_G(0:l), S_R(0:l))

        ! 1) Sample to find initial interval for exponent
        DO iter = 1, n_start + 1
          IF (iter .GT. n_start) THEN
            failed = .TRUE.
            CPWARN("GR-switch calibration FAILED: No valid initial interval was found.")
            EXIT outerloop
          ENDIF
          alpha = alpha_start(iter)
          G_rad = exp_radius(l, alpha, lattice_sum_prec, 1.0_dp, epsin = G_res)
          R_rad = exp_radius(l, 0.25_dp/alpha, lattice_sum_prec, 1.0_dp, epsin = R_res)
          G_bound = FLOOR(G_rad/(inv_lgth*twopi))
          R_bound = FLOOR(R_rad/lgth) + 1
          IF((G_bound .GT. 100) .OR. (R_bound .GT. 100)) CYCLE
          CALL CPU_TIME(t_0)
          DO i = 1, n_points
            CALL pgf_sum_gspace_1d(S_G, R(i, i_xyz), alpha, inv_lgth, -G_bound, G_bound)
          ENDDO
          CALL CPU_TIME(t_1)
          DO i = 1, n_points
            CALL pgf_sum_rspace_1d(S_R, R(i, i_xyz), 0.25_dp/alpha, lgth, -R_bound, R_bound)
          ENDDO
          CALL CPU_TIME(t_2)
          dt_G = t_1 - t_0
          dt_R = t_2 - t_1
          CALL mp_sum(dt_G, para_env%group)
          CALL mp_sum(dt_R, para_env%group)
          IF(dt_R - dt_G .GT. 0.0_dp) THEN
            IF(iter .EQ. 1) THEN
              CPWARN("GR-switch calibration FAILED: No valid initial interval was found.")
              EXIT outerloop
            ENDIF
            alpha_lr(1) = alpha_start(iter - 1)
            alpha_lr(2) = alpha_start(iter)
            EXIT
          ENDIF
        ENDDO

        ! 2) bisection
        DO iter = 1, max_iter + 1
          IF (iter .GT. max_iter) THEN
            failed = .TRUE.
            CPWARN("GR-switch calibration FAILED: Maximum number of iterations has been exceeded")
            EXIT outerloop
          ENDIF
          alpha = 0.5_dp*(alpha_lr(1) + alpha_lr(2))
          IF((alpha_lr(2) - alpha_lr(1))/alpha .LT. alpha_eps) THEN
            alpha_c(l, i_xyz) = (alpha_lr(2) + alpha_lr(1))/2
            EXIT
          ENDIF
          G_rad = exp_radius(l, alpha, lattice_sum_prec, 1.0_dp, epsin = G_res)
          R_rad = exp_radius(l, 0.25_dp/alpha, lattice_sum_prec, 1.0_dp, epsin = R_res)
          G_bound = FLOOR(G_rad/(inv_lgth*twopi))
          R_bound = FLOOR(R_rad/lgth) + 1
          CALL CPU_TIME(t_0)
          DO i = 1, n_points
            CALL pgf_sum_gspace_1d(S_G, R(i, i_xyz), alpha, inv_lgth, -G_bound, G_bound)
          ENDDO
          CALL CPU_TIME(t_1)
          DO i = 1, n_points
            CALL pgf_sum_rspace_1d(S_R, R(i, i_xyz), 0.25_dp/alpha, lgth, -R_bound, R_bound)
          ENDDO
          CALL CPU_TIME(t_2)
          dt_G = t_1 - t_0
          dt_R = t_2 - t_1
          CALL mp_sum(dt_G, para_env%group)
          CALL mp_sum(dt_R, para_env%group)
          IF(dt_R - dt_G .GT. 0.0_dp) THEN
            alpha_lr(2) = alpha
          ELSE
            alpha_lr(1) = alpha
          ENDIF
        ENDDO
        DEALLOCATE(S_G, S_R)
      ENDDO
    ENDDO outerloop

  END SUBROUTINE calibrate_gr_switch

! *****************************************************************************
!> \brief Compute 1d sum
!>        S_G(l, alpha) = inv_lgth*sum_G( C(l, alpha, delta_l, G) ) with
!>          C(l, alpha, delta_l, G) = prefactor*|G|^(l-delta_l) exp(-alpha*G^2)
!>          if G not equal 0
!>          C(l = 0, alpha, delta_l, 0) = 1, C(l>0, alpha, delta_l, 0) = 0
!>        dG = inv_lgth*twopi and G = -G_bound*dG, (-G_bound + 1)*dG, ..., G_bound*dG
!>        for all l < = l_max.
!> \param S_G ...
!> \param alpha ...
!> \param inv_lgth ...
!> \param G_min ...
!> \param G_c ...
!> \param delta_l ...
!> \param prefactor ...
!> \note  needed for cutoff error estimate
! *****************************************************************************
  SUBROUTINE pgf_sum_gspace_1d_deltal(S_G, alpha, inv_lgth, G_min, G_c, delta_l, prefactor)
    REAL(KIND=dp), DIMENSION(0:), &
      INTENT(INOUT)                          :: S_G
    REAL(KIND=dp), INTENT(IN)                :: alpha, inv_lgth
    INTEGER, INTENT(IN)                      :: G_min, G_c
    REAL(KIND=dp), INTENT(IN)                :: delta_l, prefactor

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pgf_sum_gspace_1d_deltal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, l, l_max
    REAL(KIND=dp)                            :: dG, exp_tot, G, prefac

    prefac = prefactor*inv_lgth
    dG = inv_lgth*twopi
    l_max = UBOUND(S_G, 1)

    S_G(:) = 0.0_dp
    DO k = G_min, G_c
      IF (k.NE.0) THEN
        G = k*dG
        exp_tot = EXP(-alpha*G**2)*prefac
        DO l = 0, l_max
          S_G(l) = S_G(l) + ABS(G)**(l - delta_l)*exp_tot
        ENDDO
      ELSE
        S_G(0) = S_G(0) + prefac
      ENDIF
    ENDDO

  END SUBROUTINE pgf_sum_gspace_1d_deltal

END MODULE eri_mme_methods_low
