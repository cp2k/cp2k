!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE grid_collocate
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: compute_cube_center,&
                                              cube_info_type,&
                                              return_cube,&
                                              return_cube_nonortho
   USE d3_poly,                         ONLY: poly_cp2k2d3
   USE gauss_colloc,                    ONLY: collocGauss
   USE grid_modify_pab_block,           ONLY: prepare_adb_m_dab,&
                                              prepare_ardb_m_darb,&
                                              prepare_dab_p_adb,&
                                              prepare_dadb,&
                                              prepare_diadib,&
                                              prepare_diiadiib,&
                                              prepare_dijadijb
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathconstants,                   ONLY: fac
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE realspace_grid_types,            ONLY: realspace_grid_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_collocate'

   PUBLIC :: collocate_pgf_product

   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_AB = 401
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DADB = 402
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADBmDAB = 403
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB = 404
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DABpADB = 405
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DXDY = 601
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DYDZ = 602
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZDX = 603
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DXDX = 604
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DYDY = 605
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZDZ = 606
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DX = 501
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DY = 502
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZ = 503

CONTAINS

! **************************************************************************************************
!> \brief low level collocation of primitive gaussian functions
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rab2 ...
!> \param scale ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param rsgrid ...
!> \param cell ...
!> \param cube_info ...
!> \param eps_rho_rspace ...
!> \param ga_gb_function ...
!> \param map_consistent ...
!> \param collocate_rho0 ...
!> \param rpgf0_s ...
!> \param idir ...
!> \param ir ...
!> \param rsgauge ...
!> \param rsbuf ...
!> \param use_subpatch ...
!> \param subpatch_pattern ...
!> \param lmax_global Maximum possible value of lmax used to dimension arrays
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product(la_max, zeta, la_min, &
                                    lb_max, zetb, lb_min, &
                                    ra, rab, rab2, scale, pab, o1, o2, &
                                    rsgrid, cell, cube_info, &
                                    eps_rho_rspace, ga_gb_function, &
                                    map_consistent, &
                                    collocate_rho0, &
                                    rpgf0_s, idir, ir, rsgauge, rsbuf, &
                                    use_subpatch, subpatch_pattern, &
                                    lmax_global)

      INTEGER, INTENT(IN)                      :: la_max
      REAL(KIND=dp), INTENT(IN)                :: zeta
      INTEGER, INTENT(IN)                      :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                :: zetb
      INTEGER, INTENT(IN)                      :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
      REAL(KIND=dp), INTENT(IN)                :: rab2, scale
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
      INTEGER, INTENT(IN)                      :: o1, o2
      TYPE(realspace_grid_type)                :: rsgrid
      TYPE(cell_type), POINTER                 :: cell
      TYPE(cube_info_type), INTENT(IN)         :: cube_info
      REAL(KIND=dp), INTENT(IN)                :: eps_rho_rspace
      INTEGER, INTENT(IN)                      :: ga_gb_function
      LOGICAL, INTENT(IN), OPTIONAL            :: map_consistent, collocate_rho0
      REAL(dp), INTENT(IN), OPTIONAL           :: rpgf0_s
      INTEGER, INTENT(IN), OPTIONAL            :: idir, ir
      TYPE(realspace_grid_type), POINTER, OPTIONAL :: rsgauge, rsbuf
      LOGICAL, OPTIONAL                        :: use_subpatch
      INTEGER(KIND=int_8), OPTIONAL, INTENT(IN):: subpatch_pattern
      INTEGER, INTENT(IN)                      :: lmax_global

      CHARACTER(len=*), PARAMETER :: routineN = 'collocate_pgf_product', &
                                     routineP = moduleN//':'//routineN

      INTEGER :: cmax, gridbounds(2, 3), i, ico, icoef, ider1, ider2, ig, &
                 jco, k, l, la_max_local, la_min_local, lb_max_local, lb_min_local, &
                 length, lxa, lxb, lxy, lxyz, lya, lyb, &
                 lza, lzb, o1_local, o2_local, offset, start
      INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                  ub_cube
      INTEGER, DIMENSION(:), POINTER           :: sphere_bounds
      LOGICAL                                  :: my_collocate_rho0, &
                                                  my_map_consistent, subpatch_collocate
      REAL(KIND=dp) :: a, b, binomial_k_lxa, binomial_l_lxb, cutoff, f, pg, &
                       prefactor, radius, rpg, zetp
      REAL(KIND=dp), DIMENSION(3)              :: dr, rap, rb, rbp, roffset, rp
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab_local
      REAL(KIND=dp), DIMENSION(:, :, :), &
         POINTER                                :: grid
      INTEGER :: lxp, lyp, lzp, lp, lxpm, iaxis
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: map
      REAL(kind=dp) :: p_ele
      REAL(kind=dp), DIMENSION(0:lmax_global*2, 0:lmax_global, 0:lmax_global, 3) :: alpha
      REAL(kind=dp), DIMENSION(((lmax_global*2 + 1)*(lmax_global*2 + 2)*(lmax_global*2 + 3))/6) :: coef_xyz
      REAL(kind=dp), DIMENSION(((lmax_global*2 + 1)*(lmax_global*2 + 2))/2) :: coef_xyt
      REAL(kind=dp), DIMENSION(0:lmax_global*2) :: coef_xtt

      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :) :: pol_z
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :) :: pol_y
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :) :: pol_x
      REAL(KIND=dp) :: t_exp_1, t_exp_2, t_exp_min_1, t_exp_min_2, t_exp_plus_1, t_exp_plus_2

      REAL(kind=dp), POINTER, DIMENSION(:, :, :) :: grid_tmp, gauge

      subpatch_collocate = .FALSE.

      IF (PRESENT(use_subpatch)) THEN
         IF (use_subpatch) THEN
            subpatch_collocate = .TRUE.
            CPASSERT(PRESENT(subpatch_pattern))
         ENDIF
      ENDIF

      ! use identical radii for integrate and collocate ?
      IF (PRESENT(map_consistent)) THEN
         my_map_consistent = map_consistent
      ELSE
         my_map_consistent = .FALSE.
      ENDIF

      IF (PRESENT(collocate_rho0) .AND. PRESENT(rpgf0_s)) THEN
         my_collocate_rho0 = collocate_rho0
      ELSE
         my_collocate_rho0 = .FALSE.
      END IF

      zetp = zeta + zetb
      f = zetb/zetp
      rap(:) = f*rab(:)
      rbp(:) = rap(:) - rab(:)
      rp(:) = ra(:) + rap(:)
      rb(:) = ra(:) + rab(:)

      ! check to avoid overflows
      a = MAXVAL(ABS(rsgrid%desc%dh))
      a = 300._dp/(a*a)
      !   CPASSERT(zetp < a)

      IF (my_map_consistent) THEN
         cutoff = 1.0_dp
         prefactor = EXP(-zeta*f*rab2)
         radius = exp_radius_very_extended(la_min, la_max, lb_min, lb_max, ra=ra, rb=rb, rp=rp, &
                                           zetp=zetp, eps=eps_rho_rspace, &
                                           prefactor=prefactor, cutoff=cutoff)
         prefactor = scale*EXP(-zeta*f*rab2)
      ELSE IF (my_collocate_rho0) THEN
         cutoff = 0.0_dp
         prefactor = 1.0_dp
         radius = rpgf0_s
      ELSE
         cutoff = 0.0_dp
         prefactor = scale*EXP(-zeta*f*rab2)
         radius = exp_radius_very_extended(la_min, la_max, lb_min, lb_max, pab, o1, o2, ra, rb, rp, &
                                           zetp, eps_rho_rspace, prefactor, cutoff)
      ENDIF

      a = MAXVAL(ABS(rsgrid%desc%dh))
      IF (radius .LT. a/2.0_dp) THEN
         RETURN
      END IF

      ! it's a choice to compute lX_min/max, pab here,
      ! this way we get the same radius as we use for the corresponding density
      SELECT CASE (ga_gb_function)
      CASE (GRID_FUNC_DADB)
         la_max_local = la_max + 1
         la_min_local = MAX(la_min - 1, 0)
         lb_max_local = lb_max + 1
         lb_min_local = MAX(lb_min - 1, 0)
         ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
         ! is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
         ! (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})
         ! cleaner would possibly be to touch pzyx directly (avoiding the following allocate)
         ALLOCATE (pab_local(ncoset(la_max_local), ncoset(lb_max_local)))
         pab_local = 0.0_dp
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb

                  ! this element of pab results in 12 elements of pab_local
                  CALL prepare_dadb(pab_local, pab, lxa, lya, lza, lxb, lyb, lzb, o1, o2, zeta, zetb)

               ENDDO
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         o1_local = 0
         o2_local = 0
         pab_local = pab_local*0.5_dp
      CASE (GRID_FUNC_ADBmDAB)
         CPASSERT(PRESENT(idir))
         la_max_local = la_max + 1
         la_min_local = MAX(la_min - 1, 0)
         lb_max_local = lb_max + 1
         lb_min_local = MAX(lb_min - 1, 0)
         ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
         ! is equivalent to mapping pab with
         !    pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
         ! ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) =
         !          pgf_a *(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x}) -
         !                   (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b

         ALLOCATE (pab_local(ncoset(la_max_local), ncoset(lb_max_local)))
         pab_local = 0.0_dp
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  ! this element of pab results in 4 elements of pab_local
                  CALL prepare_adb_m_dab(pab_local, pab, idir, &
                                         lxa, lya, lza, lxb, lyb, lzb, o1, o2, zeta, zetb)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         o1_local = 0
         o2_local = 0
      CASE (GRID_FUNC_DABpADB)
         CPASSERT(PRESENT(idir))
         la_max_local = la_max + 1
         la_min_local = MAX(la_min - 1, 0)
         lb_max_local = lb_max + 1
         lb_min_local = MAX(lb_min - 1, 0)
         ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
         ! is equivalent to mapping pab with
         !    pgf_a (nabla_{idir} pgf_b) + (nabla_{idir} pgf_a) pgf_b
         ! ( pgf_a ) (ddx pgf_b) + (ddx pgf_a)( pgf_b ) =
         !          pgf_a *(lbx pgf_{b-1x} + 2*zetb*pgf_{b+1x}) +
         !                   (lax pgf_{a-1x} + 2*zeta*pgf_{a+1x}) pgf_b

         ALLOCATE (pab_local(ncoset(la_max_local), ncoset(lb_max_local)))
         pab_local = 0.0_dp
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  ! this element of pab results in 4 elements of pab_local
                  CALL prepare_dab_p_adb(pab_local, pab, idir, &
                                         lxa, lya, lza, lxb, lyb, lzb, o1, o2, zeta, zetb)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         o1_local = 0
         o2_local = 0
      CASE (GRID_FUNC_DX, GRID_FUNC_DY, GRID_FUNC_DZ)
         ider1 = ga_gb_function - 500
         la_max_local = la_max + 1
         la_min_local = MAX(la_min - 1, 0)
         lb_max_local = lb_max + 1
         lb_min_local = MAX(lb_min - 1, 0)
         ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
         ! is equivalent to mapping pab with
         !   d_{ider1} pgf_a d_{ider1} pgf_b
         ! dx pgf_a dx pgf_b =
         !        (lax pgf_{a-1x})*(lbx pgf_{b-1x}) - 2*zetb*lax*pgf_{a-1x}*pgf_{b+1x} -
         !         lbx pgf_{b-1x}*2*zeta*pgf_{a+1x}+ 4*zeta*zetab*pgf_{a+1x}pgf_{b+1x}

         ALLOCATE (pab_local(ncoset(la_max_local), ncoset(lb_max_local)))
         pab_local = 0.0_dp
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  ! this element of pab results in 4 elements of pab_local
                  CALL prepare_dIadIb(pab_local, pab, ider1, &
                                      lxa, lya, lza, lxb, lyb, lzb, o1, o2, zeta, zetb)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         o1_local = 0
         o2_local = 0
      CASE (GRID_FUNC_DXDY, GRID_FUNC_DYDZ, GRID_FUNC_DZDX)
         ider1 = ga_gb_function - 600
         ider2 = ga_gb_function - 600 + 1
         IF (ider2 > 3) ider2 = ider1 - 2
         la_max_local = la_max + 2
         la_min_local = MAX(la_min - 2, 0)
         lb_max_local = lb_max + 2
         lb_min_local = MAX(lb_min - 2, 0)
         ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
         ! is equivalent to mapping pab with
         !   d_{ider1} pgf_a d_{ider1} pgf_b
         ALLOCATE (pab_local(ncoset(la_max_local), ncoset(lb_max_local)))
         pab_local = 0.0_dp
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  ! this element of pab results in 16 elements of pab_local
                  CALL prepare_dijadijb(pab_local, pab, ider1, ider2, &
                                        lxa, lya, lza, lxb, lyb, lzb, o1, o2, zeta, zetb)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         o1_local = 0
         o2_local = 0
      CASE (GRID_FUNC_DXDX, GRID_FUNC_DYDY, GRID_FUNC_DZDZ)
         ider1 = ga_gb_function - 603
         la_max_local = la_max + 2
         la_min_local = MAX(la_min - 2, 0)
         lb_max_local = lb_max + 2
         lb_min_local = MAX(lb_min - 2, 0)
         ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
         ! is equivalent to mapping pab with
         !   dd_{ider1} pgf_a dd_{ider1} pgf_b

         ALLOCATE (pab_local(ncoset(la_max_local), ncoset(lb_max_local)))
         pab_local = 0.0_dp
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb
                  ! this element of pab results in 9 elements of pab_local
                  CALL prepare_diiadiib(pab_local, pab, ider1, &
                                        lxa, lya, lza, lxb, lyb, lzb, o1, o2, zeta, zetb)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         o1_local = 0
         o2_local = 0
      CASE (GRID_FUNC_ARDBmDARB)
         CPASSERT(PRESENT(idir))
         CPASSERT(PRESENT(ir))
         la_max_local = la_max + 1
         la_min_local = MAX(la_min - 1, 0)
         lb_max_local = lb_max + 2
         lb_min_local = MAX(lb_min - 1, 0)
         ! create a new pab_local so that mapping pab_local with pgf_a pgf_b
         ! is equivalent to mapping pab with
         ! pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}  pgf_b
         ! ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b ) =
         !                        pgf_a *(lbx pgf_{b-1x+1ir} - 2*zetb*pgf_{b+1x+1ir}) -
         !                       (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}

         ALLOCATE (pab_local(ncoset(la_max_local), ncoset(lb_max_local)))
         pab_local = 0.0_dp
         DO lxa = 0, la_max
         DO lxb = 0, lb_max
            DO lya = 0, la_max - lxa
            DO lyb = 0, lb_max - lxb
               DO lza = MAX(la_min - lxa - lya, 0), la_max - lxa - lya
               DO lzb = MAX(lb_min - lxb - lyb, 0), lb_max - lxb - lyb

                  ! this element of pab results in 4 elements of pab_local
                  CALL prepare_ardb_m_darb(pab_local, pab, idir, ir, &
                                           lxa, lya, lza, lxb, lyb, lzb, o1, o2, zeta, zetb)
               END DO
               END DO
            END DO
            END DO
         END DO
         END DO
         o1_local = 0
         o2_local = 0
      CASE (GRID_FUNC_AB)
         la_max_local = la_max
         la_min_local = la_min
         lb_max_local = lb_max
         lb_min_local = lb_min
         pab_local => pab
         o1_local = o1
         o2_local = o2
      CASE DEFAULT
         CPASSERT(.FALSE.)
      END SELECT

      ng(:) = rsgrid%desc%npts(:)
      grid => rsgrid%r(:, :, :)
      gridbounds(1, 1) = LBOUND(GRID, 1)
      gridbounds(2, 1) = UBOUND(GRID, 1)
      gridbounds(1, 2) = LBOUND(GRID, 2)
      gridbounds(2, 2) = UBOUND(GRID, 2)
      gridbounds(1, 3) = LBOUND(GRID, 3)
      gridbounds(2, 3) = UBOUND(GRID, 3)

      IF (PRESENT(ir) .AND. PRESENT(rsgauge) .AND. PRESENT(rsbuf)) THEN
         grid => rsbuf%r(:, :, :)
         grid_tmp => rsgrid%r(:, :, :)
         gauge => rsgauge%r(:, :, :)
      ENDIF

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} *
!           (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
      lp = la_max_local + lb_max_local
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
      alpha(:, :, :, :) = 0.0_dp
      DO iaxis = 1, 3
      DO lxa = 0, la_max_local
      DO lxb = 0, lb_max_local
         binomial_k_lxa = 1.0_dp
         a = 1.0_dp
         DO k = 0, lxa
            binomial_l_lxb = 1.0_dp
            b = 1.0_dp
            DO l = 0, lxb
               alpha(lxa - l + lxb - k, lxa, lxb, iaxis) = alpha(lxa - l + lxb - k, lxa, lxb, iaxis) + &
                                                           binomial_k_lxa*binomial_l_lxb*a*b
               binomial_l_lxb = binomial_l_lxb*REAL(lxb - l, dp)/REAL(l + 1, dp)
               b = b*(rp(iaxis) - (ra(iaxis) + rab(iaxis)))
            ENDDO
            binomial_k_lxa = binomial_k_lxa*REAL(lxa - k, dp)/REAL(k + 1, dp)
            a = a*(-ra(iaxis) + rp(iaxis))
         ENDDO
      ENDDO
      ENDDO
      ENDDO

!
!   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and alpha(ls,lxa,lxb,1)
!   use a three step procedure
!   we don't store zeros, so counting is done using lxyz,lxy in order to have contiguous memory access in collocate_fast.F
!
      lxyz = 0
      DO lzp = 0, lp
      DO lyp = 0, lp - lzp
      DO lxp = 0, lp - lzp - lyp
         lxyz = lxyz + 1
         coef_xyz(lxyz) = 0.0_dp
      ENDDO
      ENDDO
      ENDDO
      DO lzb = 0, lb_max_local
      DO lza = 0, la_max_local
         lxy = 0
         DO lyp = 0, lp - lza - lzb
            DO lxp = 0, lp - lza - lzb - lyp
               lxy = lxy + 1
               coef_xyt(lxy) = 0.0_dp
            ENDDO
            lxy = lxy + lza + lzb
         ENDDO
         DO lyb = 0, lb_max_local - lzb
         DO lya = 0, la_max_local - lza
            lxpm = (lb_max_local - lzb - lyb) + (la_max_local - lza - lya)
            coef_xtt(0:lxpm) = 0.0_dp
            DO lxb = MAX(lb_min_local - lzb - lyb, 0), lb_max_local - lzb - lyb
            DO lxa = MAX(la_min_local - lza - lya, 0), la_max_local - lza - lya
               ico = coset(lxa, lya, lza)
               jco = coset(lxb, lyb, lzb)
               p_ele = prefactor*pab_local(o1_local + ico, o2_local + jco)
               DO lxp = 0, lxa + lxb
                  coef_xtt(lxp) = coef_xtt(lxp) + p_ele*alpha(lxp, lxa, lxb, 1)
               ENDDO
            ENDDO
            ENDDO
            lxy = 0
            DO lyp = 0, lya + lyb
               DO lxp = 0, lp - lza - lzb - lya - lyb
                  lxy = lxy + 1
                  coef_xyt(lxy) = coef_xyt(lxy) + alpha(lyp, lya, lyb, 2)*coef_xtt(lxp)
               ENDDO
               lxy = lxy + lza + lzb + lya + lyb - lyp
            ENDDO
         ENDDO
         ENDDO
         lxyz = 0
         DO lzp = 0, lza + lzb
            lxy = 0
            DO lyp = 0, lp - lza - lzb
               DO lxp = 0, lp - lza - lzb - lyp
                  lxy = lxy + 1; lxyz = lxyz + 1
                  coef_xyz(lxyz) = coef_xyz(lxyz) + alpha(lzp, lza, lzb, 3)*coef_xyt(lxy)
               ENDDO
               lxy = lxy + lza + lzb; lxyz = lxyz + lza + lzb - lzp
            ENDDO
            DO lyp = lp - lza - lzb + 1, lp - lzp
               DO lxp = 0, lp - lyp - lzp
                  lxyz = lxyz + 1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO

      IF (subpatch_collocate) THEN
         CALL collocate_general_subpatch()
      ELSE
         IF (rsgrid%desc%orthorhombic) THEN
            CALL collocate_ortho()
            ! CALL collocate_general()
         ELSE
            CALL collocate_general_wings()
            !CALL collocate_general_opt()
         END IF
      END IF

      IF (ga_gb_function /= GRID_FUNC_AB) THEN
         DEALLOCATE (pab_local)
      ENDIF

   CONTAINS

      !
      ! this treats efficiently the orthogonal case
      !
! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_ortho()

!   *** properties of the grid ***

         ! notice we're in the ortho case
         dr(1) = rsgrid%desc%dh(1, 1)
         dr(2) = rsgrid%desc%dh(2, 2)
         dr(3) = rsgrid%desc%dh(3, 3)

!   *** get the sub grid properties for the given radius ***
         CALL return_cube(cube_info, radius, lb_cube, ub_cube, sphere_bounds)
         cmax = MAXVAL(ub_cube)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)
!

         ALLOCATE (map(-cmax:cmax, 3))
         CALL compute_cube_center(cubecenter, rsgrid%desc, zeta, zetb, ra, rab)
         roffset(:) = rp(:) - REAL(cubecenter(:), dp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point
         DO i = 1, 3
            IF (rsgrid%desc%perd(i) == 1) THEN
               start = lb_cube(i)
               DO
                  offset = MODULO(cubecenter(i) + start, ng(i)) + 1 - start
                  length = MIN(ub_cube(i), ng(i) - offset) - start
                  DO ig = start, start + length
                     map(ig, i) = ig + offset
                  END DO
                  IF (start + length .GE. ub_cube(i)) EXIT
                  start = start + length + 1
               END DO
            ELSE
               ! this takes partial grid + border regions into account
               offset = MODULO(cubecenter(i) + lb_cube(i) + rsgrid%desc%lb(i) - rsgrid%lb_local(i), ng(i)) + 1 - lb_cube(i)
               ! check for out of bounds
               IF (ub_cube(i) + offset > UBOUND(grid, i) .OR. lb_cube(i) + offset < LBOUND(grid, i)) THEN
                  CPASSERT(.FALSE.)
               ENDIF
               DO ig = lb_cube(i), ub_cube(i)
                  map(ig, i) = ig + offset
               END DO
            END IF
         ENDDO
         ALLOCATE (pol_z(1:2, 0:lp, -cmax:0))
         ALLOCATE (pol_y(1:2, 0:lp, -cmax:0))
         ALLOCATE (pol_x(0:lp, -cmax:cmax))

         IF (PRESENT(ir) .AND. PRESENT(rsgauge)) CALL collocate_ortho_set_to_0()

#include "prep.f90"
#include "call_collocate.f90"

         IF (PRESENT(ir) .AND. PRESENT(rsgauge)) CALL collocate_gauge_ortho()

         ! deallocation needed to pass around a pgi bug..
         DEALLOCATE (pol_z)
         DEALLOCATE (pol_y)
         DEALLOCATE (pol_x)
         DEALLOCATE (map)

      END SUBROUTINE collocate_ortho

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_gauge_ortho()
      INTEGER                                            :: i, igmax, igmin, j, j2, jg, jg2, jgmin, &
                                                            k, k2, kg, kg2, kgmin, sci
      REAL(KIND=dp)                                      :: point(3, 4), res(4), x, y, y2, z, z2

! notice we're in the ortho case

         dr(1) = rsgrid%desc%dh(1, 1)
         dr(2) = rsgrid%desc%dh(2, 2)
         dr(3) = rsgrid%desc%dh(3, 3)
         !
         sci = 1
         kgmin = sphere_bounds(sci)
         sci = sci + 1
         DO kg = kgmin, 0
            kg2 = 1 - kg
            k = map(kg, 3)
            k2 = map(kg2, 3)
            jgmin = sphere_bounds(sci)
            sci = sci + 1
            z = (REAL(kg, dp) + REAL(cubecenter(3), dp))*dr(3)
            z2 = (REAL(kg2, dp) + REAL(cubecenter(3), dp))*dr(3)
            DO jg = jgmin, 0
               jg2 = 1 - jg
               j = map(jg, 2)
               j2 = map(jg2, 2)
               igmin = sphere_bounds(sci)
               sci = sci + 1
               igmax = 1 - igmin
               y = (REAL(jg, dp) + REAL(cubecenter(2), dp))*dr(2)
               y2 = (REAL(jg2, dp) + REAL(cubecenter(2), dp))*dr(2)
               DO ig = igmin, igmax
                  i = map(ig, 1)
                  x = (REAL(ig, dp) + REAL(cubecenter(1), dp))*dr(1)
                  point(1, 1) = x; point(2, 1) = y; point(3, 1) = z
                  point(1, 2) = x; point(2, 2) = y2; point(3, 2) = z
                  point(1, 3) = x; point(2, 3) = y; point(3, 3) = z2
                  point(1, 4) = x; point(2, 4) = y2; point(3, 4) = z2
                  !
                  res(1) = (point(ir, 1) - rb(ir)) - gauge(i, j, k)
                  res(2) = (point(ir, 2) - rb(ir)) - gauge(i, j2, k)
                  res(3) = (point(ir, 3) - rb(ir)) - gauge(i, j, k2)
                  res(4) = (point(ir, 4) - rb(ir)) - gauge(i, j2, k2)
                  !
                  grid_tmp(i, j, k) = grid_tmp(i, j, k) + grid(i, j, k)*res(1)
                  grid_tmp(i, j2, k) = grid_tmp(i, j2, k) + grid(i, j2, k)*res(2)
                  grid_tmp(i, j, k2) = grid_tmp(i, j, k2) + grid(i, j, k2)*res(3)
                  grid_tmp(i, j2, k2) = grid_tmp(i, j2, k2) + grid(i, j2, k2)*res(4)
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE collocate_gauge_ortho

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_ortho_set_to_0()
      INTEGER                                            :: i, igmax, igmin, j, j2, jg, jg2, jgmin, &
                                                            k, k2, kg, kg2, kgmin, sci

!

         dr(1) = rsgrid%desc%dh(1, 1)
         dr(2) = rsgrid%desc%dh(2, 2)
         dr(3) = rsgrid%desc%dh(3, 3)
         !
         sci = 1
         kgmin = sphere_bounds(sci)
         sci = sci + 1
         DO kg = kgmin, 0
            kg2 = 1 - kg
            k = map(kg, 3)
            k2 = map(kg2, 3)
            jgmin = sphere_bounds(sci)
            sci = sci + 1
            DO jg = jgmin, 0
               jg2 = 1 - jg
               j = map(jg, 2)
               j2 = map(jg2, 2)
               igmin = sphere_bounds(sci)
               sci = sci + 1
               igmax = 1 - igmin
               DO ig = igmin, igmax
                  i = map(ig, 1)
                  grid(i, j, k) = 0.0_dp
                  grid(i, j2, k) = 0.0_dp
                  grid(i, j, k2) = 0.0_dp
                  grid(i, j2, k2) = 0.0_dp
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE collocate_ortho_set_to_0

!
!   this is a general 'optimized' routine to do the collocation
!
! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_general_opt()

      INTEGER :: i, i_index, il, ilx, ily, ilz, index_max(3), index_min(3), ismax, ismin, j, &
         j_index, jl, jlx, jly, jlz, k, k_index, kl, klx, kly, klz, lpx, lpy, lpz, lx, ly, lz, &
         offset(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: grid_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: coef_map
      REAL(KIND=dp)                                      :: a, b, c, d, di, dip, dj, djp, dk, dkp, &
                                                            exp0i, exp1i, exp2i, gp(3), &
                                                            hmatgrid(3, 3), pointj(3), pointk(3), &
                                                            res, v(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coef_ijk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: hmatgridp

!
! transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
! sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
! sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
!

         ALLOCATE (coef_ijk(((lp + 1)*(lp + 2)*(lp + 3))/6))

         ! aux mapping array to simplify life
         ALLOCATE (coef_map(0:lp, 0:lp, 0:lp))
         coef_map = HUGE(coef_map)
         lxyz = 0
         DO lzp = 0, lp
         DO lyp = 0, lp - lzp
         DO lxp = 0, lp - lzp - lyp
            lxyz = lxyz + 1
            coef_ijk(lxyz) = 0.0_dp
            coef_map(lxp, lyp, lzp) = lxyz
         ENDDO
         ENDDO
         ENDDO

         ! cell hmat in grid points
         hmatgrid = rsgrid%desc%dh

         ! center in grid coords
         gp = MATMUL(rsgrid%desc%dh_inv, rp)
         cubecenter(:) = FLOOR(gp)

         ! transform using multinomials
         ALLOCATE (hmatgridp(3, 3, 0:lp))
         hmatgridp(:, :, 0) = 1.0_dp
         DO k = 1, lp
            hmatgridp(:, :, k) = hmatgridp(:, :, k - 1)*hmatgrid(:, :)
         ENDDO

         lpx = lp
         DO klx = 0, lpx
         DO jlx = 0, lpx - klx
         DO ilx = 0, lpx - klx - jlx
            lx = ilx + jlx + klx
            lpy = lp - lx
            DO kly = 0, lpy
            DO jly = 0, lpy - kly
            DO ily = 0, lpy - kly - jly
               ly = ily + jly + kly
               lpz = lp - lx - ly
               DO klz = 0, lpz
               DO jlz = 0, lpz - klz
               DO ilz = 0, lpz - klz - jlz
                  lz = ilz + jlz + klz

                  il = ilx + ily + ilz
                  jl = jlx + jly + jlz
                  kl = klx + kly + klz
                  coef_ijk(coef_map(il, jl, kl)) = &
                     coef_ijk(coef_map(il, jl, kl)) + coef_xyz(coef_map(lx, ly, lz))* &
                     hmatgridp(1, 1, ilx)*hmatgridp(1, 2, jlx)*hmatgridp(1, 3, klx)* &
                     hmatgridp(2, 1, ily)*hmatgridp(2, 2, jly)*hmatgridp(2, 3, kly)* &
                     hmatgridp(3, 1, ilz)*hmatgridp(3, 2, jlz)*hmatgridp(3, 3, klz)* &
                     fac(lx)*fac(ly)*fac(lz)/ &
                     (fac(ilx)*fac(ily)*fac(ilz)*fac(jlx)*fac(jly)*fac(jlz)*fac(klx)*fac(kly)*fac(klz))
               ENDDO
               ENDDO
               ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO

         CALL return_cube_nonortho(cube_info, radius, index_min, index_max, rp)

         offset(:) = MODULO(index_min(:) + rsgrid%desc%lb(:) - rsgrid%lb_local(:), ng(:)) + 1

         ALLOCATE (grid_map(index_min(1):index_max(1)))
         DO i = index_min(1), index_max(1)
            grid_map(i) = MODULO(i, ng(1)) + 1
            IF (rsgrid%desc%perd(1) == 1) THEN
               grid_map(i) = MODULO(i, ng(1)) + 1
            ELSE
               grid_map(i) = i - index_min(1) + offset(1)
            ENDIF
         ENDDO

         ! go over the grid, but cycle if the point is not within the radius
         DO k = index_min(3), index_max(3)
            dk = k - gp(3)
            pointk = hmatgrid(:, 3)*dk

            IF (rsgrid%desc%perd(3) == 1) THEN
               k_index = MODULO(k, ng(3)) + 1
            ELSE
               k_index = k - index_min(3) + offset(3)
            ENDIF

            coef_xyt = 0.0_dp
            lxyz = 0
            dkp = 1.0_dp
            DO kl = 0, lp
               lxy = 0
               DO jl = 0, lp - kl
                  DO il = 0, lp - kl - jl
                     lxyz = lxyz + 1; lxy = lxy + 1
                     coef_xyt(lxy) = coef_xyt(lxy) + coef_ijk(lxyz)*dkp
                  ENDDO
                  lxy = lxy + kl
               ENDDO
               dkp = dkp*dk
            ENDDO

            DO j = index_min(2), index_max(2)
               dj = j - gp(2)
               pointj = pointk + hmatgrid(:, 2)*dj
               IF (rsgrid%desc%perd(2) == 1) THEN
                  j_index = MODULO(j, ng(2)) + 1
               ELSE
                  j_index = j - index_min(2) + offset(2)
               ENDIF

               coef_xtt = 0.0_dp
               lxy = 0
               djp = 1.0_dp
               DO jl = 0, lp
                  DO il = 0, lp - jl
                     lxy = lxy + 1
                     coef_xtt(il) = coef_xtt(il) + coef_xyt(lxy)*djp
                  ENDDO
                  djp = djp*dj
               ENDDO

               ! find bounds for the inner loop
               ! based on a quadratic equation in i
               ! a*i**2+b*i+c=radius**2
               v = pointj - gp(1)*hmatgrid(:, 1)
               a = DOT_PRODUCT(hmatgrid(:, 1), hmatgrid(:, 1))
               b = 2*DOT_PRODUCT(v, hmatgrid(:, 1))
               c = DOT_PRODUCT(v, v)
               d = b*b - 4*a*(c - radius**2)

               IF (d < 0) THEN
                  CYCLE
               ELSE
                  d = SQRT(d)
                  ismin = CEILING((-b - d)/(2*a))
                  ismax = FLOOR((-b + d)/(2*a))
               ENDIF
               ! prepare for computing -zetp*rsq
               a = -zetp*a
               b = -zetp*b
               c = -zetp*c
               i = ismin - 1

               ! the recursion relation might have to be done
               ! from the center of the gaussian (in both directions)
               ! instead as the current implementation from an edge
               exp2i = EXP((a*i + b)*i + c)
               exp1i = EXP(2*a*i + a + b)
               exp0i = EXP(2*a)

               DO i = ismin, ismax
                  di = i - gp(1)

                  ! polynomial terms
                  res = 0.0_dp
                  dip = 1.0_dp
                  DO il = 0, lp
                     res = res + coef_xtt(il)*dip
                     dip = dip*di
                  ENDDO

                  ! the exponential recursion
                  exp2i = exp2i*exp1i
                  exp1i = exp1i*exp0i
                  res = res*exp2i

                  i_index = grid_map(i)
                  grid(i_index, j_index, k_index) = grid(i_index, j_index, k_index) + res
               ENDDO
            ENDDO
         ENDDO
         !t2=nanotime_ia32()
         !write(*,*) t2-t1
         ! deallocation needed to pass around a pgi bug..
         DEALLOCATE (coef_ijk)
         DEALLOCATE (coef_map)
         DEALLOCATE (hmatgridp)
         DEALLOCATE (grid_map)

      END SUBROUTINE collocate_general_opt

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_general_subpatch()
      INTEGER, DIMENSION(2, 3)                           :: local_b
      INTEGER, DIMENSION(3)                              :: local_s, periodic
      REAL(dp), DIMENSION((lp+1)*(lp+2)*(lp+3)/6)        :: poly_d3

         periodic = 1 ! cell%perd
         CALL poly_cp2k2d3(coef_xyz, lp, poly_d3)
         local_b(1, :) = rsgrid%lb_real - rsgrid%desc%lb
         local_b(2, :) = rsgrid%ub_real - rsgrid%desc%lb
         local_s = rsgrid%lb_real - rsgrid%lb_local
         IF (BTEST(subpatch_pattern, 0)) local_b(1, 1) = local_b(1, 1) - rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 1)) local_b(2, 1) = local_b(2, 1) + rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 2)) local_b(1, 2) = local_b(1, 2) - rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 3)) local_b(2, 2) = local_b(2, 2) + rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 4)) local_b(1, 3) = local_b(1, 3) - rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 5)) local_b(2, 3) = local_b(2, 3) + rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 0)) local_s(1) = local_s(1) - rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 2)) local_s(2) = local_s(2) - rsgrid%desc%border
         IF (BTEST(subpatch_pattern, 4)) local_s(3) = local_s(3) - rsgrid%desc%border
         CALL collocGauss(h=cell%hmat, h_inv=cell%h_inv, &
                          grid=grid, poly=poly_d3, alphai=zetp, posi=rp, max_r2=radius*radius, &
                          periodic=periodic, gdim=ng, local_bounds=local_b, local_shift=local_s)
         ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,

      END SUBROUTINE

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_general_wings()
      INTEGER, DIMENSION(2, 3)                           :: local_b
      INTEGER, DIMENSION(3)                              :: periodic
      REAL(dp), DIMENSION((lp+1)*(lp+2)*(lp+3)/6)        :: poly_d3
      REAL(dp), DIMENSION(3)                             :: local_shift, rShifted

         periodic = 1 ! cell%perd
         CALL poly_cp2k2d3(coef_xyz, lp, poly_d3)
         local_b(1, :) = 0
         local_b(2, :) = MIN(rsgrid%desc%npts - 1, rsgrid%ub_local - rsgrid%lb_local)
         local_shift = REAL(rsgrid%desc%lb - rsgrid%lb_local, dp)/REAL(rsgrid%desc%npts, dp)
         rShifted(1) = rp(1) + cell%hmat(1, 1)*local_shift(1) &
                       + cell%hmat(1, 2)*local_shift(2) &
                       + cell%hmat(1, 3)*local_shift(3)
         rShifted(2) = rp(2) + cell%hmat(2, 1)*local_shift(1) &
                       + cell%hmat(2, 2)*local_shift(2) &
                       + cell%hmat(2, 3)*local_shift(3)
         rShifted(3) = rp(3) + cell%hmat(3, 1)*local_shift(1) &
                       + cell%hmat(3, 2)*local_shift(2) &
                       + cell%hmat(3, 3)*local_shift(3)
         CALL collocGauss(h=cell%hmat, h_inv=cell%h_inv, &
                          grid=grid, poly=poly_d3, alphai=zetp, posi=rShifted, max_r2=radius*radius, &
                          periodic=periodic, gdim=ng, local_bounds=local_b)
         ! defaults: local_shift=(/0,0,0/),poly_shift=(/0.0_dp,0.0_dp,0.0_dp/),scale=1.0_dp,

      END SUBROUTINE

!
!   this is a general 'reference' routine to do the collocation
!
! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_general()

      INTEGER                                            :: i, index_max(3), index_min(3), &
                                                            ipoint(3), j, k
      REAL(KIND=dp)                                      :: inv_ng(3), point(3), primpt

! still hard-wired (see MODULO)

         CPASSERT(ALL(rsgrid%desc%perd == 1))

         CALL return_cube_nonortho(cube_info, radius, index_min, index_max, rp)

         inv_ng = 1.0_dp/ng

         ! go over the grid, but cycle if the point is not within the radius
         DO k = index_min(3), index_max(3)
         DO j = index_min(2), index_max(2)
         DO i = index_min(1), index_max(1)
            ! point in real space
            point = MATMUL(cell%hmat, REAL((/i, j, k/), KIND=dp)*inv_ng)
            ! primitive_value of point
            primpt = primitive_value(point)
            ! skip if outside of the sphere
            IF (SUM((point - rp)**2) > radius**2) CYCLE
            ! point on the grid (including pbc)
            ipoint = MODULO((/i, j, k/), ng) + 1
            ! add to grid
            grid(ipoint(1), ipoint(2), ipoint(3)) = grid(ipoint(1), ipoint(2), ipoint(3)) + primpt
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE collocate_general

! **************************************************************************************************
!> \brief ...
!> \param point ...
!> \return ...
! **************************************************************************************************
      FUNCTION primitive_value(point) RESULT(res)
      REAL(KIND=dp)                                      :: point(3), res

      REAL(KIND=dp)                                      :: dra(3), drap(3), drb(3), drbp(3), myexp, &
                                                            pdrap

         res = 0.0_dp
         myexp = EXP(-zetp*SUM((point - rp)**2))*prefactor
         dra = point - ra
         drb = point - rb
         drap(1) = 1.0_dp
         DO lxa = 0, la_max_local
            drbp(1) = 1.0_dp
            DO lxb = 0, lb_max_local
               drap(2) = 1.0_dp
               DO lya = 0, la_max_local - lxa
                  drbp(2) = 1.0_dp
                  DO lyb = 0, lb_max_local - lxb
                     drap(3) = 1.0_dp
                     DO lza = 1, MAX(la_min_local - lxa - lya, 0)
                        drap(3) = drap(3)*dra(3)
                     ENDDO
                     DO lza = MAX(la_min_local - lxa - lya, 0), la_max_local - lxa - lya
                        drbp(3) = 1.0_dp
                        DO lzb = 1, MAX(lb_min_local - lxb - lyb, 0)
                           drbp(3) = drbp(3)*drb(3)
                        ENDDO
                        ico = coset(lxa, lya, lza)
                        pdrap = PRODUCT(drap)
                        DO lzb = MAX(lb_min_local - lxb - lyb, 0), lb_max_local - lxb - lyb
                           jco = coset(lxb, lyb, lzb)
                           res = res + pab_local(ico + o1_local, jco + o2_local)*myexp*pdrap*PRODUCT(drbp)
                           drbp(3) = drbp(3)*drb(3)
                        ENDDO
                        drap(3) = drap(3)*dra(3)
                     ENDDO
                     drbp(2) = drbp(2)*drb(2)
                  ENDDO
                  drap(2) = drap(2)*dra(2)
               ENDDO
               drbp(1) = drbp(1)*drb(1)
            ENDDO
            drap(1) = drap(1)*dra(1)
         ENDDO

      END FUNCTION primitive_value

   END SUBROUTINE collocate_pgf_product

END MODULE grid_collocate
