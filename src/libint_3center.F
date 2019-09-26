!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief 3-center electron repulsion integral routines based on libint2
!>        Currently available operators: Coulomb, Truncated Coulomb, Short Range (erfc)
!> \author A. Bussy (05.2019)
! **************************************************************************************************

MODULE libint_3center
   USE gamma,                           ONLY: fgamma => fgamma_0
   USE hfx_libint_wrapper,              ONLY: cp_libint_get_3eris,&
                                              cp_libint_set_params_eri,&
                                              cp_libint_t,&
                                              prim_data_f_size
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              t_c_g0_n
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_3center'

   PUBLIC :: eri_3center

   TYPE :: params_3c
      INTEGER                               :: m_max
      REAL(dp)                              :: ZetaInv, EtaInv, ZetapEtaInv, Rho
      REAL(dp), DIMENSION(3)                :: Q, W
      REAL(dp), DIMENSION(prim_data_f_size) :: Fm
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Computes the 3-center electron repulsion integrals (ab|c) for a given set of cartesian
!>        gaussian orbitals
!> \param int_abc the integrals as array of cartesian orbitals (allocated before hand)
!> \param la_min ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param ra ...
!> \param lb_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rb ...
!> \param lc_min ...
!> \param lc_max ...
!> \param npgfc ...
!> \param zetc ...
!> \param rc ...
!> \param lib the libint_t object for evaluation (assume that it is initialized outside)
!> \param op which operator to use for the integrals
!> \param omega in case of short range erfc operator
!> \param r_cutoff the cutoff radius in case of truncated Coulomb
!> \note Prior to calling this routine, the cp_libint_t type passed as argument must be initialized,
!>       the libint library must be static initialized, and in case of truncated Coulomb operator,
!>       the latter must be initialized too
! **************************************************************************************************
   SUBROUTINE eri_3center(int_abc, la_min, la_max, npgfa, zeta, ra, &
                          lb_min, lb_max, npgfb, zetb, rb, &
                          lc_min, lc_max, npgfc, zetc, rc, &
                          lib, op, omega, r_cutoff)

      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: int_abc
      INTEGER, INTENT(IN)                                :: la_min, la_max, npgfa
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra
      INTEGER, INTENT(IN)                                :: lb_min, lb_max, npgfb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rb
      INTEGER, INTENT(IN)                                :: lc_min, lc_max, npgfc
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetc
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rc
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      INTEGER, INTENT(IN)                                :: op
      REAL(dp), INTENT(IN), OPTIONAL                     :: omega, r_cutoff

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_3center', routineP = moduleN//':'//routineN

      INTEGER                                            :: a_mysize(1), a_offset, a_start, &
                                                            b_offset, b_start, c_offset, c_start, &
                                                            i, ipgf, j, jpgf, k, kpgf, li, lj, lk, &
                                                            ncoa, ncob, ncoc, p1, p2, p3
      REAL(dp)                                           :: my_omega, my_r_cutoff, zeti, zetj, zetk
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      TYPE(params_3c), POINTER                           :: params

      NULLIFY (params, p_work)
      ALLOCATE (params)

      !Note: we want to compute all possible integrals based on the 3-centers (ab|c) before
      !      having to switch to (ba|c) (or the other way around) due to angular momenta in libint
      !      For a triplet of centers (k|ji), we can only compute integrals for which lj >= li

      my_r_cutoff = 0.0_dp
      IF (PRESENT(r_cutoff)) my_r_cutoff = r_cutoff

      my_omega = 0.0_dp
      IF (PRESENT(omega)) my_omega = omega

      !Looping over the pgfs
      DO ipgf = 1, npgfa
         zeti = zeta(ipgf)
         a_start = (ipgf-1)*ncoset(la_max)

         DO jpgf = 1, npgfb
            zetj = zetb(jpgf)
            b_start = (jpgf-1)*ncoset(lb_max)

            DO kpgf = 1, npgfc
               zetk = zetc(kpgf)
               c_start = (kpgf-1)*ncoset(lc_max)

               !start with all the (c|ba) integrals (standard order) and keep to lb >= la
               CALL set_params(lib, ra, rb, rc, la_max, lb_max, lc_max, zeti, zetj, zetk, op, &
                               params_out=params, r_cutoff=my_r_cutoff, omega=my_omega)

               DO li = la_min, la_max
                  a_offset = a_start+ncoset(li-1)
                  ncoa = nco(li)
                  DO lj = MAX(li, lb_min), lb_max
                     b_offset = b_start+ncoset(lj-1)
                     ncob = nco(lj)
                     DO lk = lc_min, lc_max
                        c_offset = c_start+ncoset(lk-1)
                        ncoc = nco(lk)

                        a_mysize(1) = ncoa*ncob*ncoc
                        CALL cp_libint_get_3eris(li, lj, lk, lib, p_work, a_mysize)

                        DO k = 1, ncoc
                           p1 = (k-1)*ncob
                           DO j = 1, ncob
                              p2 = (p1+j-1)*ncoa
                              DO i = 1, ncoa
                                 p3 = p2+i
                                 int_abc(a_offset+i, b_offset+j, c_offset+k) = p_work(p3)
                              END DO
                           END DO
                        END DO

                     END DO !lk
                  END DO !lj
               END DO !li

               !swap centers 3 and 4 to compute (c|ab) with lb < la
               CALL set_params(lib, rb, ra, rc, params_in=params)

               DO lj = lb_min, lb_max
                  b_offset = b_start+ncoset(lj-1)
                  ncob = nco(lj)
                  DO li = MAX(lj+1, la_min), la_max
                     a_offset = a_start+ncoset(li-1)
                     ncoa = nco(li)
                     DO lk = lc_min, lc_max
                        c_offset = c_start+ncoset(lk-1)
                        ncoc = nco(lk)

                        a_mysize(1) = ncoa*ncob*ncoc
                        CALL cp_libint_get_3eris(lj, li, lk, lib, p_work, a_mysize)

                        DO k = 1, ncoc
                           p1 = (k-1)*ncoa
                           DO i = 1, ncoa
                              p2 = (p1+i-1)*ncob
                              DO j = 1, ncob
                                 p3 = p2+j
                                 int_abc(a_offset+i, b_offset+j, c_offset+k) = p_work(p3)
                              END DO
                           END DO
                        END DO

                     END DO !lk
                  END DO !li
               END DO !lj

            END DO !kpgf
         END DO !jpgf
      END DO !ipgf

      DEALLOCATE (params)

   END SUBROUTINE eri_3center

! **************************************************************************************************
!> \brief Sets the internals of the cp_libint_t object for integrals of type (k|ji)
!> \param lib ..
!> \param ri ...
!> \param rj ...
!> \param rk ...
!> \param li_max ...
!> \param lj_max ...
!> \param lk_max ...
!> \param zeti ...
!> \param zetj ...
!> \param zetk ...
!> \param op ...
!> \param r_cutoff ...
!> \param omega ...
!> \param params_in external parameters to use for libint
!> \param params_out returns the libint parameters computed based on the other arguments
!> \note The use of params_in and params_out comes from the fact that one might have to swap
!>       centers 3 and 4 because of angular momenta and pretty much all the parameters of libint
!>       remain the same upon such a change => might avoid recomputing things over and over again
! **************************************************************************************************
   SUBROUTINE set_params(lib, ri, rj, rk, li_max, lj_max, lk_max, zeti, zetj, zetk, &
                         op, r_cutoff, omega, params_in, params_out)

      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ri, rj, rk
      INTEGER, INTENT(IN), OPTIONAL                      :: li_max, lj_max, lk_max
      REAL(dp), INTENT(IN), OPTIONAL                     :: zeti, zetj, zetk
      INTEGER, INTENT(IN), OPTIONAL                      :: op
      REAL(dp), INTENT(IN), OPTIONAL                     :: r_cutoff, omega
      TYPE(params_3c), OPTIONAL, POINTER                 :: params_in, params_out

      CHARACTER(len=*), PARAMETER :: routineN = 'set_params', routineP = moduleN//':'//routineN

      INTEGER                                            :: l
      LOGICAL                                            :: use_gamma
      REAL(dp)                                           :: gammaq, omega2, omega_corr, omega_corr2, &
                                                            prefac, R, S1234, T, tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: Fm
      TYPE(params_3c), POINTER                           :: params

      !Assume that one of params_in or params_out is present, and that in the latter case, all
      !other optinal arguments are here

      !The internal structure of libint2 is based on 4-center integrals
      !For 3-center, one of those is a dummy center
      !The integral is assumed to be (k|ji) where the centers are ordered as:
      !k -> 1, j -> 3 and i -> 4 (the center #2 is the dummy center)

      !If external parameters are given, just use them
      IF (PRESENT(params_in)) THEN
         params => params_in

         !If no external parameters to use, compute them
      ELSE
         params => params_out

         !Note: some variable of 4-center integrals simplify with a dummy center:
         !      P -> rk, gammap -> zetk
         params%m_max = li_max+lj_max+lk_max
         gammaq = zeti+zetj
         params%ZetaInv = 1._dp/zetk; params%EtaInv = 1._dp/gammaq
         params%ZetapEtaInv = 1._dp/(zetk+gammaq)

         params%Q = (zeti*ri+zetj*rj)*params%EtaInv
         params%W = (zetk*rk+gammaq*params%Q)*params%ZetapEtaInv
         params%Rho = zetk*gammaq/(zetk+gammaq)

         SELECT CASE (op)
         CASE (do_potential_coulomb)
            T = params%Rho*SUM((params%Q-rk)**2)
            S1234 = EXP(-zeti*zetj*params%EtaInv*SUM((rj-ri)**2))
            prefac = 2._dp*pi/params%Rho*SQRT((pi*params%ZetapEtaInv)**3)*S1234

            params%Fm = 0.0_dp
            CALL fgamma(params%m_max, T, params%Fm)
            params%Fm = prefac*params%Fm
         CASE (do_potential_truncated)
            R = r_cutoff*SQRT(params%Rho)
            T = params%Rho*SUM((params%Q-rk)**2)
            S1234 = EXP(-zeti*zetj*params%EtaInv*SUM((rj-ri)**2))
            prefac = 2._dp*pi/params%Rho*SQRT((pi*params%ZetapEtaInv)**3)*S1234

            CPASSERT(get_lmax_init() .GE. params%m_max) !check if truncated coulomb init correctly
            CALL t_c_g0_n(params%Fm, use_gamma, R, T, params%m_max)
            IF (use_gamma) CALL fgamma(params%m_max, T, params%Fm)
            params%Fm = prefac*params%Fm
         CASE (do_potential_short)
            T = params%Rho*SUM((params%Q-rk)**2)
            S1234 = EXP(-zeti*zetj*params%EtaInv*SUM((rj-ri)**2))
            prefac = 2._dp*pi/params%Rho*SQRT((pi*params%ZetapEtaInv)**3)*S1234

            params%Fm = 0.0_dp
            CALL fgamma(params%m_max, T, params%Fm)

            omega2 = omega**2
            omega_corr2 = omega2/(omega2+params%Rho)
            omega_corr = SQRT(omega_corr2)
            T = T*omega_corr2
            ALLOCATE (Fm(prim_data_f_size))

            CALL fgamma(params%m_max, T, Fm)
            tmp = -omega_corr
            DO l = 1, params%m_max+1
               params%Fm(l) = params%Fm(l)+Fm(l)*tmp
               tmp = tmp*omega_corr2
            END DO
            params%Fm = prefac*params%Fm
         CASE DEFAULT
            CPABORT("Requested operator NYI")
         END SELECT

      END IF

      CALL cp_libint_set_params_eri(lib, rk, rk, rj, ri, params%ZetaInv, params%EtaInv, &
                                    params%ZetapEtaInv, params%Rho, rk, params%Q, params%W, &
                                    params%m_max, params%Fm)

   END SUBROUTINE set_params

END MODULE libint_3center

