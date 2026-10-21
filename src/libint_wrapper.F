!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the Libint-Library or a c++ wrapper.
!> \par History
!>      11.2007 created [Manuel Guidon]
!>      10.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
MODULE libint_wrapper

#if(__LIBINT)
#include <libint2/libint2_params.h>
#include <libint2/config.h>
#endif

! maximum angular momentum to be supported in CP2K-LIBINT interface
   #:set libint_max_am_supported = 8

   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_F_POINTER, &
                            C_F_PROCPOINTER, &
                            C_NULL_PTR, &
                            C_FUNPTR
   USE cp_log_handling, ONLY: cp_to_string
   USE kinds, ONLY: dp
#if(__LIBINT)
   USE libint_f, ONLY: &
      libint2_build, libint2_build_eri, libint2_build_eri1, libint2_cleanup_eri, &
      libint2_cleanup_eri1, libint2_init_eri, libint2_init_eri1, libint2_static_cleanup, &
      libint2_static_init, libint_t, libint2_max_am_eri, libint2_init_3eri, libint2_cleanup_3eri, &
      libint2_init_2eri, libint2_cleanup_2eri, &
      libint2_build_2eri, libint2_build_3eri, libint2_build_3eri1, libint2_cleanup_3eri1, libint2_init_3eri1, &
      libint2_build_2eri1, libint2_cleanup_2eri1, libint2_init_2eri1
#endif
   USE orbital_pointers, ONLY: nco
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: cp_libint_t, prim_data_f_size, &
             libint_max_am, cp_libint_get_eris, cp_libint_get_derivs, &
             cp_libint_init_eri, cp_libint_init_eri1, cp_libint_cleanup_eri, &
             cp_libint_cleanup_eri1, cp_libint_static_init, cp_libint_static_cleanup, &
             cp_libint_set_contrdepth, &
             cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, &
             cp_libint_init_3eri, cp_libint_cleanup_3eri, cp_libint_get_3eris, &
             cp_libint_init_2eri, cp_libint_cleanup_2eri, cp_libint_get_2eris, &
             cp_libint_get_3eri_derivs, cp_libint_init_3eri1, cp_libint_cleanup_3eri1, &
             cp_libint_get_2eri_derivs, cp_libint_init_2eri1, cp_libint_cleanup_2eri1

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_wrapper'

#if(__LIBINT)
   INTEGER, PARAMETER :: libint_max_am = libint2_max_am_eri
#else
   INTEGER, PARAMETER :: libint_max_am = 0
#endif

   INTEGER, PARAMETER :: prim_data_f_size = 4*(libint_max_am) + 1

!> \brief Prepared primitive in original center order for the legacy Libint backend.
!> \note Zero exponents denote dummy centers. Fm contains orders zero through prepared_max_order.
   TYPE, PUBLIC :: eri_primitive_data_type
      INTEGER :: prepared_max_order = -1
      REAL(dp) :: centers(3, 4) = 0.0_dp, exponents(4) = 0.0_dp
      REAL(dp) :: P(3) = 0.0_dp, Q(3) = 0.0_dp, W(3) = 0.0_dp
      REAL(dp) :: ZetaInv = 0.0_dp, EtaInv = 0.0_dp, ZetapEtaInv = 0.0_dp, Rho = 0.0_dp
      REAL(dp) :: Fm(prim_data_f_size) = 0.0_dp
   END TYPE eri_primitive_data_type

   TYPE :: cp_libint_t
      PRIVATE
#if(__LIBINT)
      TYPE(libint_t), DIMENSION(1) :: prv
#else
      INTEGER :: unused = -1
#endif
   END TYPE cp_libint_t

CONTAINS

! **************************************************************************************************
!> \brief Fill an ERI backend from one prepared primitive and a pair-symmetry permutation.
!> \param libint backend workspace
!> \param params primitive data in original center order
!> \param permutation original center index for each backend center; preserves the ERI pair structure
!> \param m_max requested radial order, not exceeding the prepared order
! **************************************************************************************************
   SUBROUTINE cp_libint_set_params_eri(libint, params, permutation, m_max)
      TYPE(cp_libint_t), INTENT(INOUT) :: libint
      TYPE(eri_primitive_data_type), INTENT(IN) :: params
      INTEGER, INTENT(IN) :: permutation(4), m_max

#if(__LIBINT)
      INTEGER :: pair_order(2)
      REAL(dp) :: product_centers(3, 2), inverse_exponents(2)
#endif

      CPASSERT(m_max <= params%prepared_max_order)

#if(__LIBINT)
      pair_order = [1, 2]
      IF (permutation(1) > 2) pair_order = [2, 1]
      product_centers(:, 1) = params%P
      product_centers(:, 2) = params%Q
      inverse_exponents = [params%ZetaInv, params%EtaInv]

      ASSOCIATE (A => params%centers(:, permutation(1)), B => params%centers(:, permutation(2)), &
                 C => params%centers(:, permutation(3)), D => params%centers(:, permutation(4)), &
                 P => product_centers(:, pair_order(1)), Q => product_centers(:, pair_order(2)), &
                 W => params%W, ZetaInv => inverse_exponents(pair_order(1)), &
                 EtaInv => inverse_exponents(pair_order(2)), ZetapEtaInv => params%ZetapEtaInv, &
                 Rho => params%Rho, F => params%Fm)
         libint%prv(1)%AB_x(1) = A(1) - B(1)
         libint%prv(1)%AB_y(1) = A(2) - B(2)
         libint%prv(1)%AB_z(1) = A(3) - B(3)

         libint%prv(1)%CD_x(1) = C(1) - D(1)
         libint%prv(1)%CD_y(1) = C(2) - D(2)
         libint%prv(1)%CD_z(1) = C(3) - D(3)

         libint%prv(1)%PA_x(1) = P(1) - A(1)
         libint%prv(1)%PA_y(1) = P(2) - A(2)
         libint%prv(1)%PA_z(1) = P(3) - A(3)

         libint%prv(1)%QC_x(1) = Q(1) - C(1)
         libint%prv(1)%QC_y(1) = Q(2) - C(2)
         libint%prv(1)%QC_z(1) = Q(3) - C(3)

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

            IF (${m_max}$ <= m_max) THEN
               libint%prv(1)%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_${m_max}$ (1) &
                  = F(${m_max}$+1)
            END IF
#endif
         #:endfor

      END ASSOCIATE
#else
      MARK_USED(libint)
      MARK_USED(params)
      MARK_USED(permutation)
      MARK_USED(m_max)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_set_params_eri

! **************************************************************************************************
!> \brief Fill common ERI fields and the additional first-derivative recurrence fields.
!> \param libint backend workspace
!> \param params primitive data in original center order
!> \param permutation original center index for each backend center; preserves the ERI pair structure
!> \param m_max requested radial order, including the first-derivative increment
! **************************************************************************************************
   SUBROUTINE cp_libint_set_params_eri_deriv(libint, params, permutation, m_max)
      TYPE(cp_libint_t), INTENT(INOUT) :: libint
      TYPE(eri_primitive_data_type), INTENT(IN) :: params
      INTEGER, INTENT(IN) :: permutation(4), m_max

#if(__LIBINT)
      REAL(dp) :: P(3), exponents(4), gammap, gammaq, gammapq, rhop, rhoq
#endif

      CALL cp_libint_set_params_eri(libint, params, permutation, m_max)

#if(__LIBINT)
      P = params%P
      IF (permutation(1) > 2) P = params%Q
      exponents = params%exponents(permutation)
      ASSOCIATE (B => params%centers(:, permutation(2)), &
                 Zeta_A => exponents(1), Zeta_B => exponents(2), &
                 Zeta_C => exponents(3), Zeta_D => exponents(4))
         libint%prv(1)%PB_x(1) = P(1) - B(1)
         libint%prv(1)%PB_y(1) = P(2) - B(2)
         libint%prv(1)%PB_z(1) = P(3) - B(3)

         libint%prv(1)%two_alpha0_bra(1) = 2.0_dp*Zeta_A
         libint%prv(1)%two_alpha0_ket(1) = 2.0_dp*Zeta_B
         libint%prv(1)%two_alpha1_ket(1) = 2.0_dp*Zeta_D

         gammap = Zeta_A + Zeta_B
         gammaq = Zeta_C + Zeta_D
         gammapq = gammap*gammaq/(gammap + gammaq)
         libint%prv(1)%alpha1_rho_over_zeta2(1) = Zeta_A*gammapq/(gammap*gammap)
         libint%prv(1)%alpha2_rho_over_zeta2(1) = Zeta_B*gammapq/(gammap*gammap)
         libint%prv(1)%alpha4_rho_over_eta2(1) = Zeta_D*gammapq/(gammaq*gammaq)
         libint%prv(1)%alpha1_over_zetapluseta(1) = Zeta_A/(gammap + gammaq)
         libint%prv(1)%alpha2_over_zetapluseta(1) = Zeta_B/(gammap + gammaq)
         libint%prv(1)%alpha4_over_zetapluseta(1) = Zeta_D/(gammap + gammaq)

         rhop = Zeta_A*Zeta_B/gammap
         rhoq = Zeta_C*Zeta_D/gammaq
         libint%prv(1)%rho12_over_alpha1(1) = rhop/Zeta_A

         libint%prv(1)%rho34_over_alpha3(1) = rhoq/Zeta_C

      END ASSOCIATE
#endif
   END SUBROUTINE cp_libint_set_params_eri_deriv
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
      CHARACTER(LEN=160)                                 :: error_msg
      PROCEDURE(libint2_build), POINTER                  :: pbuild
      TYPE(C_FUNPTR)                                     :: build_funptr

      build_funptr = libint2_build_3eri(n_c, n_b, n_a)
      IF (.NOT. C_ASSOCIATED(build_funptr)) THEN
         error_msg = "LIBINT has no 3-center ERI kernel for (lc,lb,la)=("// &
                     TRIM(cp_to_string(n_c))//","//TRIM(cp_to_string(n_b))//","// &
                     TRIM(cp_to_string(n_a))//")"
         CALL cp_abort(__LOCATION__, TRIM(error_msg))
      END IF
      CALL C_F_PROCPOINTER(build_funptr, pbuild)
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
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_3eri_derivs(n_c, n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      INTEGER                                            :: a_mysize(1)
      REAL(dp), DIMENSION(:, :), POINTER                    :: p_work

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work_tmp
      PROCEDURE(libint2_build), POINTER                  :: pbuild
      INTEGER :: i

      CALL C_F_PROCPOINTER(libint2_build_3eri1(n_c, n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      ALLOCATE (p_work(a_mysize(1), 9))

      !Derivatives 1-3 can be obtained using translational invariance
      DO i = 4, 9
         NULLIFY (p_work_tmp)
         CALL C_F_POINTER(lib%prv(1)%targets(i), p_work_tmp, SHAPE=a_mysize)
         p_work(:, i) = p_work_tmp
      END DO
#else
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_3eri_derivs

! **************************************************************************************************
!> \brief ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_2eri_derivs(n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      INTEGER                                            :: a_mysize(1)
      REAL(dp), DIMENSION(:, :), POINTER                 :: p_work

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work_tmp
      PROCEDURE(libint2_build), POINTER                  :: pbuild
      INTEGER :: i

      CALL C_F_PROCPOINTER(libint2_build_2eri1(n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      ALLOCATE (p_work(a_mysize(1), 6))

      !Derivatives 1-3 can be obtained using translational invariance
      DO i = 4, 6
         NULLIFY (p_work_tmp)
         CALL C_F_POINTER(lib%prv(1)%targets(i), p_work_tmp, SHAPE=a_mysize)
         p_work(:, i) = p_work_tmp
      END DO
#else
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_2eri_derivs

! **************************************************************************************************
!> \brief ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
!> \param a_mysize ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_2eris(n_b, n_a, lib, p_work, a_mysize)
      INTEGER, INTENT(IN)                                :: n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      INTEGER                                            :: a_mysize(1)

#if(__LIBINT)
      PROCEDURE(libint2_build), POINTER                  :: pbuild

      CALL C_F_PROCPOINTER(libint2_build_2eri(n_b, n_a), pbuild)
      CALL pbuild(lib%prv)

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=a_mysize)
#else
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      MARK_USED(a_mysize)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_2eris

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
      REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*nco(n_c)*nco &
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
         END DO
      END DO
      ! Libint omits the canonical second-center derivatives; recover them by translation invariance.
      DO k = 4, 6
         DO i = 1, a_mysize(1)
            work_forces(i, k) = -1.0_dp*(work_forces(i, k - 3) + &
                                         work_forces(i, k + 3) + work_forces(i, k + 6))
         END DO
      END DO
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
   END SUBROUTINE cp_libint_init_eri

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

   END SUBROUTINE cp_libint_set_contrdepth

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
   END SUBROUTINE cp_libint_init_eri1

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
   END SUBROUTINE cp_libint_init_3eri

   SUBROUTINE cp_libint_init_3eri1(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_3eri1(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_3eri1

   SUBROUTINE cp_libint_init_2eri1(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_2eri1(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_2eri1

   SUBROUTINE cp_libint_init_2eri(lib, max_am)
      TYPE(cp_libint_t)               :: lib
      INTEGER                         :: max_am
#if(__LIBINT)
      CALL libint2_init_2eri(lib%prv, max_am, C_NULL_PTR)
#else
      MARK_USED(lib)
      MARK_USED(max_am)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_init_2eri

   SUBROUTINE cp_libint_cleanup_eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_eri

   SUBROUTINE cp_libint_cleanup_eri1(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_eri1(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_eri1

   SUBROUTINE cp_libint_cleanup_3eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_3eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_3eri

   SUBROUTINE cp_libint_cleanup_3eri1(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_3eri1(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_3eri1

   SUBROUTINE cp_libint_cleanup_2eri1(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_2eri1(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_2eri1

   SUBROUTINE cp_libint_cleanup_2eri(lib)
      TYPE(cp_libint_t)               :: lib
#if(__LIBINT)
      CALL libint2_cleanup_2eri(lib%prv)
#else
      MARK_USED(lib)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_cleanup_2eri

   SUBROUTINE cp_libint_static_init()
#if(__LIBINT)
      CALL libint2_static_init()
#else
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_static_init

   SUBROUTINE cp_libint_static_cleanup()
#if(__LIBINT)
      CALL libint2_static_cleanup()
#else
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_static_cleanup

END MODULE libint_wrapper
