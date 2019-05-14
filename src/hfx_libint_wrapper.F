!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper.
!> \par History
!>      11.2007 created [Manuel Guidon]
!>      10.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_libint_wrapper

#if(__LIBINT)
#include <libint2/libint2_params.h>
#include <libint2/config.h>
#endif

! maximum angular momentum to be supported in CP2K-LIBINT interface
#:set libint_max_am_supported = 8

   USE ISO_C_BINDING,                   ONLY: C_F_POINTER,&
                                              C_F_PROCPOINTER,&
                                              C_NULL_PTR
   USE kinds,                           ONLY: dp
#if(__LIBINT)
   USE libint_f,                        ONLY: &
        libint2_build, libint2_build_eri, libint2_build_eri1, libint2_cleanup_eri, &
        libint2_cleanup_eri1, libint2_init_eri, libint2_init_eri1, libint2_static_cleanup, &
        libint2_static_init, libint_t, libint2_max_am_eri, libint2_init_3eri, libint2_cleanup_3eri,&
        libint2_build_3eri
#endif
   USE orbital_pointers,                ONLY: nco
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: cp_libint_t, prim_data_f_size, build_eri_size, build_deriv1_eri_size, &
             libint_max_am, libderiv_max_am1, cp_libint_get_eris, cp_libint_get_derivs, &
             cp_libint_init_eri, cp_libint_init_eri1, cp_libint_cleanup_eri, &
             cp_libint_cleanup_eri1, cp_libint_static_init, cp_libint_static_cleanup, &
             get_ssss_f_val,cp_libint_set_contrdepth,cp_libint_set_params_eri_screen,&
             cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, &
             cp_libint_init_3eri, cp_libint_cleanup_3eri, cp_libint_get_3eris

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_libint_wrapper'

#if(__LIBINT)
   INTEGER, PARAMETER :: libint_max_am = libint2_max_am_eri
#else
   INTEGER, PARAMETER :: libint_max_am = 0
#endif

   INTEGER, PARAMETER :: libderiv_max_am1 = libint_max_am
   INTEGER, PARAMETER :: prim_data_f_size = 4*(libint_max_am)+1
   INTEGER, PARAMETER :: libint_vrr_classes_size = 2*(libint_max_am)+1
   INTEGER, PARAMETER :: libint_dvrr_classes_size = 2*(libderiv_max_am1)+1
   INTEGER, PARAMETER :: build_eri_size = libint_max_am
   INTEGER, PARAMETER :: build_deriv1_eri_size = libderiv_max_am1

   TYPE :: cp_libint_t
      PRIVATE
#if(__LIBINT)
      TYPE(libint_t), DIMENSION(1) :: prv
#else
      INTEGER :: unused = -1
#endif
   END TYPE

CONTAINS


   SUBROUTINE cp_libint_set_params_eri_screen(libint, A, B, C, D, P, Q, W, ZetaInv, EtaInv, ZetapEtaInv, Rho, m_max, F)
      TYPE(cp_libint_t)                       :: libint
      REAL(KIND=dp), INTENT(IN), DIMENSION(3) :: A, B, C, D, P, Q, W
      REAL(KIND=dp), INTENT(IN)               :: ZetaInv, EtaInv, ZetapEtaInv, Rho
      INTEGER, INTENT(IN)                     :: m_max
      REAL(KIND=dp), DIMENSION(:)             :: F

#if(__LIBINT)
      libint%prv(1)%AB_x(1) = A(1)-B(1)
      libint%prv(1)%AB_y(1) = A(2)-B(2)
      libint%prv(1)%AB_z(1) = A(3)-B(3)

      libint%prv(1)%CD_x(1) = C(1)-D(1)
      libint%prv(1)%CD_y(1) = C(2)-D(2)
      libint%prv(1)%CD_z(1) = C(3)-D(3)

      libint%prv(1)%PA_x(1) = P(1)-A(1)
      libint%prv(1)%PA_y(1) = P(2)-A(2)
      libint%prv(1)%PA_z(1) = P(3)-A(3)

      libint%prv(1)%QC_x(1) = Q(1)-C(1)
      libint%prv(1)%QC_y(1) = Q(2)-C(2)
      libint%prv(1)%QC_z(1) = Q(3)-C(3)

      libint%prv(1)%WP_x(1) = W(1)-P(1)
      libint%prv(1)%WP_y(1) = W(2)-P(2)
      libint%prv(1)%WP_z(1) = W(3)-P(3)

      libint%prv(1)%WQ_x(1) = W(1)-Q(1)
      libint%prv(1)%WQ_y(1) = W(2)-Q(2)
      libint%prv(1)%WQ_z(1) = W(3)-Q(3)

      libint%prv(1)%oo2z(1) = 0.5_dp*ZetaInv
      libint%prv(1)%oo2e(1) = 0.5_dp*EtaInv
      libint%prv(1)%oo2ze(1) = 0.5_dp*ZetapEtaInv
      libint%prv(1)%roz(1) = Rho*ZetaInv
      libint%prv(1)%roe(1) = Rho*EtaInv

#:for m_max in range(0, 4*libint_max_am_supported)
#if 4*LIBINT2_MAX_AM_eri > ${m_max}$ - 1
      IF(${m_max}$ .LE. m_max) &
         libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_${m_max}$(1) &
         = F(${m_max}$+1)

#endif
#:endfor

#else
      MARK_USED(libint)
      MARK_USED(A)
      MARK_USED(B)
      MARK_USED(C)
      MARK_USED(D)
      MARK_USED(P)
      MARK_USED(Q)
      MARK_USED(W)
      MARK_USED(ZetaInv)
      MARK_USED(EtaInv)
      MARK_USED(ZetapEtaInv)
      MARK_USED(Rho)
      MARK_USED(m_max)
      MARK_USED(F)

      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE

   SUBROUTINE cp_libint_set_params_eri_deriv(libint, A, B, C, D, P, Q, W, zeta_A, zeta_B, zeta_C, zeta_D, &
                                             ZetaInv, EtaInv, ZetapEtaInv, Rho, m_max, F)
      TYPE(cp_libint_t)                       :: libint
      REAL(KIND=dp), INTENT(IN), DIMENSION(3) :: A, B, C, D, P, Q, W
      REAL(KIND=dp), INTENT(IN)               :: zeta_A, zeta_B, zeta_C, zeta_D, ZetaInv, EtaInv, ZetapEtaInv, Rho

      INTEGER, INTENT(IN)                     :: m_max
      REAL(KIND=dp), DIMENSION(:)             :: F

#if(__LIBINT)
      REAL(KIND=dp)                           :: gammap, gammaq, gammapq, rhop, rhoq
      libint%prv(1)%AB_x(1) = A(1)-B(1)
      libint%prv(1)%AB_y(1) = A(2)-B(2)
      libint%prv(1)%AB_z(1) = A(3)-B(3)

      libint%prv(1)%CD_x(1) = C(1)-D(1)
      libint%prv(1)%CD_y(1) = C(2)-D(2)
      libint%prv(1)%CD_z(1) = C(3)-D(3)

      libint%prv(1)%PA_x(1) = P(1)-A(1)
      libint%prv(1)%PA_y(1) = P(2)-A(2)
      libint%prv(1)%PA_z(1) = P(3)-A(3)

      libint%prv(1)%PB_x(1) = P(1)-B(1)
      libint%prv(1)%PB_y(1) = P(2)-B(2)
      libint%prv(1)%PB_z(1) = P(3)-B(3)

      libint%prv(1)%QC_x(1) = Q(1)-C(1)
      libint%prv(1)%QC_y(1) = Q(2)-C(2)
      libint%prv(1)%QC_z(1) = Q(3)-C(3)

      libint%prv(1)%WP_x(1) = W(1)-P(1)
      libint%prv(1)%WP_y(1) = W(2)-P(2)
      libint%prv(1)%WP_z(1) = W(3)-P(3)

      libint%prv(1)%WQ_x(1) = W(1)-Q(1)
      libint%prv(1)%WQ_y(1) = W(2)-Q(2)
      libint%prv(1)%WQ_z(1) = W(3)-Q(3)

      libint%prv(1)%two_alpha0_bra(1) = 2.0_dp*Zeta_A
      libint%prv(1)%two_alpha0_ket(1) = 2.0_dp*Zeta_B
      libint%prv(1)%two_alpha1_ket(1) = 2.0_dp*Zeta_D

      gammap = Zeta_A + Zeta_B
      gammaq = Zeta_C + Zeta_D
      gammapq = gammap*gammaq/(gammap + gammaq)
      libint%prv(1)%alpha1_rho_over_zeta2(1) = Zeta_A*gammapq/(gammap*gammap)
      libint%prv(1)%alpha2_rho_over_zeta2(1) = Zeta_B*gammapq/(gammap*gammap)
      libint%prv(1)%alpha4_rho_over_eta2(1) = Zeta_D*gammapq/(gammaq*gammaq)
      libint%prv(1)%alpha1_over_zetapluseta(1) = Zeta_A/(gammap+gammaq)
      libint%prv(1)%alpha2_over_zetapluseta(1) = Zeta_B/(gammap+gammaq)
      libint%prv(1)%alpha4_over_zetapluseta(1) = Zeta_D/(gammap+gammaq)

      rhop = Zeta_A*Zeta_B/gammap
      rhoq = Zeta_C*Zeta_D/gammaq
      libint%prv(1)%rho12_over_alpha1(1) = rhop/Zeta_A

      libint%prv(1)%rho34_over_alpha3(1) = rhoq/Zeta_C


      libint%prv(1)%oo2z(1) = 0.5_dp*ZetaInv
      libint%prv(1)%oo2e(1) = 0.5_dp*EtaInv
      libint%prv(1)%oo2ze(1) = 0.5_dp*ZetapEtaInv
      libint%prv(1)%roz(1) = Rho*ZetaInv
      libint%prv(1)%roe(1) = Rho*EtaInv

#:for m_max in range(0, 4*libint_max_am_supported)
#if 4*LIBINT2_MAX_AM_eri > ${m_max}$ - 1
      IF(${m_max}$ .LE. m_max) &
         libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_${m_max}$(1) & ! ERROR: __LIBINT_MAX_AM is too large
         = F(${m_max}$+1)
#endif
#:endfor

#else
      MARK_USED(libint)
      MARK_USED(A)
      MARK_USED(B)
      MARK_USED(C)
      MARK_USED(D)
      MARK_USED(P)
      MARK_USED(Q)
      MARK_USED(W)
      MARK_USED(zeta_A)
      MARK_USED(zeta_B)
      MARK_USED(zeta_C)
      MARK_USED(zeta_D)
      MARK_USED(ZetaInv)
      MARK_USED(EtaInv)
      MARK_USED(ZetapEtaInv)
      MARK_USED(Rho)
      MARK_USED(m_max)
      MARK_USED(F)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE

   SUBROUTINE cp_libint_set_params_eri(libint, A, B, C, D, ZetaInv, EtaInv, ZetapEtaInv, Rho, P, Q, W, m_max, F)
      TYPE(cp_libint_t)                       :: libint
      REAL(KIND=dp), INTENT(IN), DIMENSION(3) :: A, B, C, D, P, Q, W
      REAL(KIND=dp), INTENT(IN)               :: ZetaInv, EtaInv, ZetapEtaInv, Rho
      REAL(KIND=dp), DIMENSION(:)             :: F

      INTEGER, INTENT(IN)                     :: m_max

#if(__LIBINT)
      libint%prv(1)%AB_x(1) = A(1)-B(1)
      libint%prv(1)%AB_y(1) = A(2)-B(2)
      libint%prv(1)%AB_z(1) = A(3)-B(3)

      libint%prv(1)%CD_x(1) = C(1)-D(1)
      libint%prv(1)%CD_y(1) = C(2)-D(2)
      libint%prv(1)%CD_z(1) = C(3)-D(3)

      libint%prv(1)%PA_x(1) = P(1)-A(1)
      libint%prv(1)%PA_y(1) = P(2)-A(2)
      libint%prv(1)%PA_z(1) = P(3)-A(3)

      libint%prv(1)%QC_x(1) = Q(1)-C(1)
      libint%prv(1)%QC_y(1) = Q(2)-C(2)
      libint%prv(1)%QC_z(1) = Q(3)-C(3)

      libint%prv(1)%WP_x(1) = W(1) - P(1)
      libint%prv(1)%WP_y(1) = W(2) - P(2)
      libint%prv(1)%WP_z(1) = W(3) - P(3)

      libint%prv(1)%WQ_x(1) = W(1) - Q(1)
      libint%prv(1)%WQ_y(1) = W(2) - Q(2)
      libint%prv(1)%WQ_z(1) = W(3) - Q(3)

      libint%prv(1)%oo2z(1) = 0.5_dp*ZetaInv
      libint%prv(1)%oo2e(1) = 0.5_dp*EtaInv
      libint%prv(1)%oo2ze(1) = 0.5_dp*ZetapEtaInv
      libint%prv(1)%roz(1) = Rho*ZetaInv
      libint%prv(1)%roe(1) = Rho*EtaInv

#:for m_max in range(0, 4*libint_max_am_supported)
#if 4*LIBINT2_MAX_AM_eri > ${m_max}$ - 1

      IF(${m_max}$ .LE. m_max) &
         libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_${m_max}$(1) & ! ERROR: __LIBINT_MAX_AM is too large
         = F(${m_max}$+1)
#endif
#:endfor

#else
      MARK_USED(libint)
      MARK_USED(A)
      MARK_USED(B)
      MARK_USED(C)
      MARK_USED(D)
      MARK_USED(P)
      MARK_USED(Q)
      MARK_USED(W)
      MARK_USED(ZetaInv)
      MARK_USED(EtaInv)
      MARK_USED(ZetapEtaInv)
      MARK_USED(Rho)
      MARK_USED(m_max)
      MARK_USED(F)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_eris(n_d, n_c, n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      PROCEDURE(libint2_build), POINTER                  :: pbuild

      CALL C_F_PROCPOINTER(libint2_build_eri(n_d, n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=a_mysize)
#else
      MARK_USED(n_d)
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_eris

! **************************************************************************************************
!> \brief ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_3eris(n_c, n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      PROCEDURE(libint2_build), POINTER                  :: pbuild

      CALL C_F_PROCPOINTER(libint2_build_3eri(n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=a_mysize)
#else
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_3eris

! **************************************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param work_forces ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_derivs(n_d, n_c, n_b, n_a, lib, work_forces, a_mysize)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*nco(n_c)*nco&
         (n_d), 12)                                      :: work_forces
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

      PROCEDURE(libint2_build), POINTER                  :: pbuild
      INTEGER                                            :: i, k
#endif


#if(__LIBINT)
      CALL C_F_PROCPOINTER(libint2_build_eri1(n_d, n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      DO k = 1, 12
         IF (k == 4 .OR. k == 5 .OR. k == 6) CYCLE
         CALL C_F_POINTER(lib%prv(1)%targets(k), p_work, SHAPE=a_mysize)
         DO i = 1, a_mysize(1)
            work_forces(i, k) = p_work(i)
         ENDDO
      ENDDO
#else
      MARK_USED(n_d)
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(work_forces)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_derivs

   FUNCTION get_ssss_f_val(lib)
      TYPE(cp_libint_t) :: lib
      REAL(KIND=dp) :: get_ssss_f_val

      get_ssss_f_val = 0

#if(__LIBINT)
      get_ssss_f_val = lib%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0(1)
#else
      MARK_USED(lib)
      MARK_USED(get_ssss_f_val)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END FUNCTION

   SUBROUTINE cp_libint_init_eri(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_eri(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

   SUBROUTINE cp_libint_set_contrdepth(lib, contrdepth)
      TYPE(cp_libint_t) :: lib
      INTEGER           :: contrdepth
#if(__LIBINT)
#if LIBINT_CONTRACTED_INTS
      lib%prv(1)%contrdepth = contrdepth
#else
      MARK_USED(lib)
      MARK_USED(contrdepth)
#endif
#else
      MARK_USED(lib)
      MARK_USED(contrdepth)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE

   SUBROUTINE cp_libint_init_eri1(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_eri1(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

   SUBROUTINE cp_libint_init_3eri(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_3eri(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

   SUBROUTINE cp_libint_cleanup_eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

   SUBROUTINE cp_libint_cleanup_eri1(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_eri1(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

   SUBROUTINE cp_libint_cleanup_3eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_3eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

   SUBROUTINE cp_libint_static_init()
#if(__LIBINT)
      CALL libint2_static_init()
#else
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

   SUBROUTINE cp_libint_static_cleanup()
#if(__LIBINT)
      CALL libint2_static_cleanup()
#else
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE

END MODULE hfx_libint_wrapper
