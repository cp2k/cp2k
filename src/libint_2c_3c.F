!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Two- and three-center electron repulsion integrals with owning backend workspaces.
!>        Currently available operators: Coulomb, Truncated Coulomb, Short Range (erfc), Overlap
!> \author A. Bussy (05.2019)
! **************************************************************************************************

MODULE libint_2c_3c
   USE ai_coulomb,                      ONLY: coulomb2,&
                                              coulomb3
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: active_integral_library,&
                                              coulomb_operator_type,&
                                              cutoff_screen_factor,&
                                              library_libint,&
                                              library_native
   USE kinds,                           ONLY: dp
   USE libint_engine,                   ONLY: cp_libint_engine_context_type,&
                                              cp_libint_engine_evaluate,&
                                              cp_libint_engine_init,&
                                              cp_libint_engine_release,&
                                              libint_max_am
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_2c_3c'

   PUBLIC :: compute_coulomb_2c, compute_coulomb_3c, compute_coulomb_3c_derivs, compute_coulomb_2c_derivs, &
             coulomb_integral_init, coulomb_integral_cleanup, coulomb_integral_max_am

   ! Compiled angular-momentum capability, exposed through the common Coulomb interface.
   INTEGER, PARAMETER :: coulomb_integral_max_am = libint_max_am

   TYPE, PUBLIC :: coulomb_integral_context_type
      PRIVATE
      TYPE(cp_libint_engine_context_type) :: engine_2c = cp_libint_engine_context_type(), &
                                             engine_3c = cp_libint_engine_context_type()
      INTEGER :: library_2c = library_native, library_3c = library_native
      INTEGER :: derivative_order_2c = -1, derivative_order_3c = -1

   END TYPE coulomb_integral_context_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize the requested per-thread tasks and record their backend and derivative order.
!> \param context owning task context; pass by reference without intrinsic copying
!> \param max_am_2c maximum two-center angular momentum, if requested
!> \param max_am_3c maximum three-center angular momentum, if requested
!> \param derivative_order zero for values or one for center derivatives
!> \param library optional explicit backend for existing library-specific workflows
! **************************************************************************************************
   SUBROUTINE coulomb_integral_init(context, max_am_2c, max_am_3c, derivative_order, library)
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      INTEGER, INTENT(IN), OPTIONAL                      :: max_am_2c, max_am_3c, derivative_order, &
                                                            library

      INTEGER                                            :: order

      order = 0
      IF (PRESENT(derivative_order)) order = derivative_order
      CPASSERT(order == 0 .OR. order == 1)
      IF (PRESENT(max_am_2c)) THEN
         CPASSERT(context%derivative_order_2c == -1)
         context%library_2c = active_integral_library%coulomb2_library
         IF (PRESENT(library)) context%library_2c = library
         SELECT CASE (context%library_2c)
         CASE (library_libint)
            CALL cp_libint_engine_init(context%engine_2c, 2, max_am_2c, order)
         CASE (library_native)
            IF (order /= 0) THEN
               CPABORT("Native two-center Coulomb derivatives are not implemented by this interface")
            END IF
         CASE DEFAULT
            CPABORT("Requested two-center Coulomb library is not implemented")
         END SELECT
         context%derivative_order_2c = order
      END IF
      IF (PRESENT(max_am_3c)) THEN
         CPASSERT(context%derivative_order_3c == -1)
         context%library_3c = active_integral_library%coulomb3_library
         IF (PRESENT(library)) context%library_3c = library
         SELECT CASE (context%library_3c)
         CASE (library_libint)
            CALL cp_libint_engine_init(context%engine_3c, 3, max_am_3c, order)
         CASE (library_native)
            IF (order /= 0) THEN
               CPABORT("Native three-center Coulomb derivatives are not implemented by this interface")
            END IF
         CASE DEFAULT
            CPABORT("Requested three-center Coulomb library is not implemented")
         END SELECT
         context%derivative_order_3c = order
      END IF
   END SUBROUTINE coulomb_integral_init

! **************************************************************************************************
!> \brief Release the initialized tasks and retained derivative scratch of their original owner.
!> \param context owning task context
! **************************************************************************************************
   SUBROUTINE coulomb_integral_cleanup(context)
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context

      CALL cp_libint_engine_release(context%engine_2c)
      CALL cp_libint_engine_release(context%engine_3c)
      context%derivative_order_2c = -1
      context%derivative_order_3c = -1
   END SUBROUTINE coulomb_integral_cleanup

! **************************************************************************************************
!> \brief Computes the 3-center electron repulsion integrals (ab|c) for a given set of cartesian
!>        gaussian orbitals
!> \param context owning per-thread task and derivative scratch
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
!> \param rpgfa ...
!> \param rpgfb ...
!> \param rpgfc ...
!> \param ra ...
!> \param rb ...
!> \param rc ...
!> \param int_abc the integrals as array of cartesian orbitals (allocated before hand)
!> \param potential_parameter the info about the potential
!> \param dab optional center distance for primitive screening
!> \param dac optional center distance for primitive screening
!> \param dbc optional center distance for primitive screening
!> \param int_abc_ext the extremal value of int_abc, i.e., MAXVAL(ABS(int_abc))
!>        Available for the Libint backend only.
!> \note Each thread owns its Engine context and reusable result scratch.
! **************************************************************************************************
   SUBROUTINE compute_coulomb_3c(context, la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
                                 npgfa, npgfb, npgfc, zeta, zetb, zetc, rpgfa, rpgfb, rpgfc, &
                                 ra, rb, rc, int_abc, potential_parameter, dab, dac, dbc, int_abc_ext)

      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max, npgfa, npgfb, npgfc
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta, zetb, zetc, rpgfa, rpgfb, rpgfc
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb, rc
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: int_abc
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      REAL(dp), INTENT(IN), OPTIONAL                     :: dab, dac, dbc
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: int_abc_ext

      INTEGER                                            :: a_offset, a_start, b_offset, b_start, &
                                                            c_offset, c_start, i, ipgf, j, jpgf, &
                                                            k, kpgf, li, lj, lk, ncoa, ncob, ncoc, &
                                                            op, p1, p2, p3
      REAL(dp)                                           :: distance_ab, distance_ac, distance_bc, &
                                                            dr_ab, dr_ac, dr_bc, zeti, zetj, zetk
      REAL(dp), POINTER                                  :: p_work(:, :)

      CPASSERT(context%derivative_order_3c == 0)

      IF (context%library_3c == library_native) THEN
         BLOCK
            INTEGER                                            :: kpgf, nc_end, nc_start, ncoa, ncob, ncoc
            REAL(dp)                                           :: rab2, rac2, rbc2
            REAL(dp), ALLOCATABLE, DIMENSION(:)                :: f_work, gccc, rpgfa_work, rpgfb_work, &
                                                                  rpgfc_work
            REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: vabc
            REAL(dp), ALLOCATABLE, DIMENSION(:, :, :, :)       :: v_work

            ncoa = npgfa*ncoset(la_max)
            ncob = npgfb*ncoset(lb_max)
            ncoc = ncoset(lc_max)
            rab2 = (rb(1) - ra(1))**2 + (rb(2) - ra(2))**2 + (rb(3) - ra(3))**2
            rac2 = (rc(1) - ra(1))**2 + (rc(2) - ra(2))**2 + (rc(3) - ra(3))**2
            rbc2 = (rc(1) - rb(1))**2 + (rc(2) - rb(2))**2 + (rc(3) - rb(3))**2

            ! Native evaluation does not provide Libint screening extrema.
            CPASSERT(.NOT. PRESENT(int_abc_ext))
            ALLOCATE (f_work(0:la_max + lb_max + lc_max + 2))
            f_work(:) = 0.0_dp
            ALLOCATE (v_work(ncoa, ncob, ncoc, la_max + lb_max + lc_max + 1))
            v_work(:, :, :, :) = 0.0_dp

            ! no screening
            ALLOCATE (rpgfa_work(npgfa))
            ALLOCATE (rpgfb_work(npgfb))
            ALLOCATE (rpgfc_work(npgfc))
            rpgfa_work(:) = 1.0E10_dp
            rpgfb_work(:) = 1.0E10_dp
            rpgfc_work(:) = 1.0E10_dp
            ALLOCATE (gccc(ncoc))
            gccc(:) = 0.0_dp
            ALLOCATE (vabc(ncoa, ncob))
            vabc(:, :) = 0.0_dp

            ! coulomb3 handles a single primitive Gaussian on center C per call.
            DO kpgf = 1, npgfc
               ! Each primitive occupies an ncoc-wide block in the uncontracted array.
               nc_start = (kpgf - 1)*ncoc + ncoset(lc_min - 1) + 1
               nc_end = kpgf*ncoc

               CALL coulomb3(la_max, npgfa, zeta(:), rpgfa_work(:), la_min, &
                             lb_max, npgfb, zetb(:), rpgfb_work(:), lb_min, &
                             lc_max, zetc(kpgf), rpgfc_work(kpgf), lc_min, &
                             gccc, rb - ra, rab2, rc - ra, rac2, rbc2, &
                             vabc, int_abc(:, :, nc_start:nc_end), v_work, f_work)
            END DO
         END BLOCK
         RETURN
      END IF
      CPASSERT(context%library_3c == library_libint)

      IF (PRESENT(dab)) THEN
         distance_ab = dab
      ELSE
         distance_ab = SQRT((rb(1) - ra(1))**2 + (rb(2) - ra(2))**2 + (rb(3) - ra(3))**2)
      END IF
      IF (PRESENT(dac)) THEN
         distance_ac = dac
      ELSE
         distance_ac = SQRT((rc(1) - ra(1))**2 + (rc(2) - ra(2))**2 + (rc(3) - ra(3))**2)
      END IF
      IF (PRESENT(dbc)) THEN
         distance_bc = dbc
      ELSE
         distance_bc = SQRT((rc(1) - rb(1))**2 + (rc(2) - rb(2))**2 + (rc(3) - rb(3))**2)
      END IF

      NULLIFY (p_work)

      dr_ab = 0.0_dp
      dr_bc = 0.0_dp
      dr_ac = 0.0_dp

      op = potential_parameter%potential_type

      IF (op == do_potential_truncated .OR. op == do_potential_short &
          .OR. op == do_potential_mix_cl_trunc) THEN
         dr_bc = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ac = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSE IF (op == do_potential_coulomb) THEN
         dr_bc = 1000000.0_dp
         dr_ac = 1000000.0_dp
      END IF

      IF (PRESENT(int_abc_ext)) THEN
         int_abc_ext = 0.0_dp
      END IF

      !Looping over the pgfs
      DO ipgf = 1, npgfa
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max)

         DO jpgf = 1, npgfb

            ! screening
            IF (rpgfa(ipgf) + rpgfb(jpgf) + dr_ab < distance_ab) CYCLE

            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max)

            DO kpgf = 1, npgfc

               ! screening
               IF (rpgfb(jpgf) + rpgfc(kpgf) + dr_bc < distance_bc) CYCLE
               IF (rpgfa(ipgf) + rpgfc(kpgf) + dr_ac < distance_ac) CYCLE

               zetk = zetc(kpgf)
               c_start = (kpgf - 1)*ncoset(lc_max)

               DO li = la_min, la_max
                  a_offset = a_start + ncoset(li - 1)
                  ncoa = nco(li)
                  DO lj = lb_min, lb_max
                     b_offset = b_start + ncoset(lj - 1)
                     ncob = nco(lj)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        CALL cp_libint_engine_evaluate(context%engine_3c, [lk, lj, li], &
                                                       RESHAPE([rc, rb, ra], [3, 3]), [1, 1, 1], &
                                                       RESHAPE([zetk, zetj, zeti], [1, 3]), &
                                                       RESHAPE([1.0_dp, 1.0_dp, 1.0_dp], [1, 3]), &
                                                       potential_parameter, p_work)

                        IF (PRESENT(int_abc_ext)) int_abc_ext = MAX(int_abc_ext, MAXVAL(ABS(p_work)))
                        DO k = 1, ncoc
                           p1 = (k - 1)*ncob
                           DO j = 1, ncob
                              p2 = (p1 + j - 1)*ncoa
                              DO i = 1, ncoa
                                 p3 = p2 + i
                                 int_abc(a_offset + i, b_offset + j, c_offset + k) = p_work(p3, 1)
                              END DO
                           END DO
                        END DO

                     END DO !lk
                  END DO !lj
               END DO !li

            END DO !kpgf
         END DO !jpgf
      END DO !ipgf

   END SUBROUTINE compute_coulomb_3c

! **************************************************************************************************
!> \brief Computes the derivatives of the 3-center electron repulsion integrals (ab|c) for a given
!>        set of cartesian gaussian orbitals. Returns x,y,z derivatives for 1st and 2nd center
!> \param der_abc_1 the derivatives for the 1st center (allocated before hand)
!> \param der_abc_2 the derivatives for the 2nd center (allocated before hand)
!> \param la_min ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param ra ...
!> \param lb_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param rb ...
!> \param lc_min ...
!> \param lc_max ...
!> \param npgfc ...
!> \param zetc ...
!> \param rpgfc ...
!> \param rc ...
!> \param dab ...
!> \param dac ...
!> \param dbc ...
!> \param context owning per-thread task and derivative scratch
!> \param potential_parameter the info about the potential
!> \param der_abc_1_ext the extremal value of der_abc_1, i.e., MAXVAL(ABS(der_abc_1))
!> \param der_abc_2_ext ...
!> \note Each thread owns its Engine context and reusable result scratch.
! **************************************************************************************************
   SUBROUTINE compute_coulomb_3c_derivs(der_abc_1, der_abc_2, &
                                        la_min, la_max, npgfa, zeta, rpgfa, ra, &
                                        lb_min, lb_max, npgfb, zetb, rpgfb, rb, &
                                        lc_min, lc_max, npgfc, zetc, rpgfc, rc, &
                                        dab, dac, dbc, context, potential_parameter, &
                                        der_abc_1_ext, der_abc_2_ext)

      REAL(dp), DIMENSION(:, :, :, :), INTENT(INOUT)     :: der_abc_1, der_abc_2
      INTEGER, INTENT(IN)                                :: la_min, la_max, npgfa
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta, rpgfa
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra
      INTEGER, INTENT(IN)                                :: lb_min, lb_max, npgfb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetb, rpgfb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rb
      INTEGER, INTENT(IN)                                :: lc_min, lc_max, npgfc
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetc, rpgfc
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rc
      REAL(KIND=dp), INTENT(IN)                          :: dab, dac, dbc
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      REAL(dp), DIMENSION(3), INTENT(OUT), OPTIONAL      :: der_abc_1_ext, der_abc_2_ext

      INTEGER :: a_offset, a_start, b_offset, b_start, c_offset, c_start, i, i_deriv, ipgf, j, &
         jpgf, k, kpgf, li, lj, lk, ncart, ncoa, ncob, ncoc, op, p1, p2, p3
      LOGICAL                                            :: do_ext
      REAL(dp)                                           :: dr_ab, dr_ac, dr_bc, zeti, zetj, zetk
      REAL(dp), DIMENSION(3)                             :: der_abc_1_ext_prv, der_abc_2_ext_prv
      REAL(dp), POINTER                                  :: p_work(:, :)

      CPASSERT(context%derivative_order_3c == 1)
      CPASSERT(context%library_3c == library_libint)

      dr_ab = 0.0_dp
      dr_bc = 0.0_dp
      dr_ac = 0.0_dp

      op = potential_parameter%potential_type

      IF (op == do_potential_truncated .OR. op == do_potential_short &
          .OR. op == do_potential_mix_cl_trunc) THEN
         dr_bc = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ac = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSE IF (op == do_potential_coulomb) THEN
         dr_bc = 1000000.0_dp
         dr_ac = 1000000.0_dp
      END IF

      do_ext = .FALSE.
      IF (PRESENT(der_abc_1_ext) .OR. PRESENT(der_abc_2_ext)) do_ext = .TRUE.
      der_abc_1_ext_prv = 0.0_dp
      der_abc_2_ext_prv = 0.0_dp

      !Looping over the pgfs
      DO ipgf = 1, npgfa
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max)

         DO jpgf = 1, npgfb

            ! screening
            IF (rpgfa(ipgf) + rpgfb(jpgf) + dr_ab < dab) CYCLE

            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max)

            DO kpgf = 1, npgfc

               ! screening
               IF (rpgfb(jpgf) + rpgfc(kpgf) + dr_bc < dbc) CYCLE
               IF (rpgfa(ipgf) + rpgfc(kpgf) + dr_ac < dac) CYCLE

               zetk = zetc(kpgf)
               c_start = (kpgf - 1)*ncoset(lc_max)

               DO li = la_min, la_max
                  a_offset = a_start + ncoset(li - 1)
                  ncoa = nco(li)
                  DO lj = lb_min, lb_max
                     b_offset = b_start + ncoset(lj - 1)
                     ncob = nco(lj)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        ncart = ncoa*ncob*ncoc

                        CALL cp_libint_engine_evaluate(context%engine_3c, [lk, lj, li], &
                                                       RESHAPE([rc, rb, ra], [3, 3]), [1, 1, 1], &
                                                       RESHAPE([zetk, zetj, zeti], [1, 3]), &
                                                       RESHAPE([1.0_dp, 1.0_dp, 1.0_dp], [1, 3]), &
                                                       potential_parameter, p_work)

                        DO i_deriv = 1, 3
                           IF (do_ext) THEN
                              der_abc_1_ext_prv(i_deriv) = MAX(der_abc_1_ext_prv(i_deriv), &
                                                               MAXVAL(ABS(p_work(1:ncart, i_deriv + 6))))
                              der_abc_2_ext_prv(i_deriv) = MAX(der_abc_2_ext_prv(i_deriv), &
                                                               MAXVAL(ABS(p_work(1:ncart, i_deriv + 3))))
                           END IF
                           DO k = 1, ncoc
                              p1 = (k - 1)*ncob
                              DO j = 1, ncob
                                 p2 = (p1 + j - 1)*ncoa
                                 DO i = 1, ncoa
                                    p3 = p2 + i
                                    der_abc_1(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                       p_work(p3, i_deriv + 6)
                                    der_abc_2(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                       p_work(p3, i_deriv + 3)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO !lk
                  END DO !lj
               END DO !li

            END DO !kpgf
         END DO !jpgf
      END DO !ipgf

      IF (PRESENT(der_abc_1_ext)) der_abc_1_ext = der_abc_1_ext_prv
      IF (PRESENT(der_abc_2_ext)) der_abc_2_ext = der_abc_2_ext_prv

   END SUBROUTINE compute_coulomb_3c_derivs

! **************************************************************************************************
!> \brief Computes the 2-center electron repulsion integrals (a|b) for a given set of cartesian
!>        gaussian orbitals
!> \param context owning per-thread task and derivative scratch
!> \param la_min ...
!> \param la_max ...
!> \param lb_min ...
!> \param lb_max ...
!> \param npgfa ...
!> \param npgfb ...
!> \param zeta ...
!> \param zetb ...
!> \param rpgfa ...
!> \param rpgfb ...
!> \param ra ...
!> \param rb ...
!> \param int_ab the integrals as array of cartesian orbitals (allocated before hand)
!> \param potential_parameter the info about the potential
!> \param dab optional center distance for primitive screening
!> \note Each thread owns its Engine context and reusable result scratch.
! **************************************************************************************************
   SUBROUTINE compute_coulomb_2c(context, la_min, la_max, lb_min, lb_max, npgfa, npgfb, &
                                 zeta, zetb, rpgfa, rpgfb, ra, rb, int_ab, potential_parameter, dab)

      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, npgfa, &
                                                            npgfb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta, zetb, rpgfa, rpgfb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: int_ab
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      REAL(dp), INTENT(IN), OPTIONAL                     :: dab

      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            ipgf, j, jpgf, li, lj, ncoa, ncob, p1, &
                                                            p2
      REAL(dp)                                           :: distance_ab, dr_ab, zeti, zetj
      REAL(dp), POINTER                                  :: p_work(:, :)

      CPASSERT(context%derivative_order_2c == 0)

      IF (context%library_2c == library_native) THEN
         BLOCK
            INTEGER                                            :: ncoa, ncob
            REAL(dp)                                           :: rab2
            REAL(dp), ALLOCATABLE, DIMENSION(:)                :: f_work
            REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: v_work

            ncoa = npgfa*ncoset(la_max)
            ncob = npgfb*ncoset(lb_max)
            rab2 = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2

            ALLOCATE (f_work(0:la_max + lb_max + 2))
            ALLOCATE (v_work(ncoa, ncob, la_max + lb_max + 1))

            CALL coulomb2(la_max, npgfa, zeta, rpgfa, la_min, lb_max, npgfb, zetb, rpgfb, lb_min, &
                          rb - ra, rab2, int_ab, v_work, f_work, screening=.FALSE.)
         END BLOCK
         RETURN
      END IF
      CPASSERT(context%library_2c == library_libint)

      IF (PRESENT(dab)) THEN
         distance_ab = dab
      ELSE
         distance_ab = SQRT((ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2)
      END IF

      NULLIFY (p_work)

      dr_ab = 0.0_dp

      IF (potential_parameter%potential_type == do_potential_truncated .OR. &
          potential_parameter%potential_type == do_potential_short .OR. &
          potential_parameter%potential_type == do_potential_mix_cl_trunc) THEN
         dr_ab = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSE IF (potential_parameter%potential_type == do_potential_coulomb) THEN
         dr_ab = 1000000.0_dp
      END IF

      !Looping over the pgfs
      DO ipgf = 1, npgfa
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max)

         DO jpgf = 1, npgfb
            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max)

            !screening
            IF (rpgfa(ipgf) + rpgfb(jpgf) + dr_ab < distance_ab) CYCLE

            DO li = la_min, la_max
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)

                  CALL cp_libint_engine_evaluate(context%engine_2c, [lj, li], &
                                                 RESHAPE([rb, ra], [3, 2]), [1, 1], &
                                                 RESHAPE([zetj, zeti], [1, 2]), &
                                                 RESHAPE([1.0_dp, 1.0_dp], [1, 2]), potential_parameter, p_work)

                  DO j = 1, ncob
                     p1 = (j - 1)*ncoa
                     DO i = 1, ncoa
                        p2 = p1 + i
                        int_ab(a_offset + i, b_offset + j) = p_work(p2, 1)
                     END DO
                  END DO

               END DO
            END DO

         END DO
      END DO

   END SUBROUTINE compute_coulomb_2c

!> \brief Computes the 2-center derivatives of the electron repulsion integrals (a|b) for a given
!>        set of cartesian gaussian orbitals. Returns the derivatives wrt to the first center
!> \param der_ab the derivatives as array of cartesian orbitals (allocated before hand)
!> \param la_min ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param ra ...
!> \param lb_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param rb ...
!> \param dab ...
!> \param context owning per-thread task and derivative scratch
!> \param potential_parameter the info about the potential
!> \note Each thread owns its Engine context and reusable result scratch.
! **************************************************************************************************
   SUBROUTINE compute_coulomb_2c_derivs(der_ab, la_min, la_max, npgfa, zeta, rpgfa, ra, &
                                        lb_min, lb_max, npgfb, zetb, rpgfb, rb, &
                                        dab, context, potential_parameter)

      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: der_ab
      INTEGER, INTENT(IN)                                :: la_min, la_max, npgfa
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta, rpgfa
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra
      INTEGER, INTENT(IN)                                :: lb_min, lb_max, npgfb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetb, rpgfb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rb
      REAL(dp), INTENT(IN)                               :: dab
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter

      INTEGER                                            :: a_offset, a_start, b_offset, b_start, i, &
                                                            i_deriv, ipgf, j, jpgf, li, lj, ncoa, &
                                                            ncob, p1, p2
      REAL(dp)                                           :: dr_ab, zeti, zetj
      REAL(dp), POINTER                                  :: p_work(:, :)

      CPASSERT(context%derivative_order_2c == 1)
      CPASSERT(context%library_2c == library_libint)

      dr_ab = 0.0_dp

      IF (potential_parameter%potential_type == do_potential_truncated .OR. &
          potential_parameter%potential_type == do_potential_short .OR. &
          potential_parameter%potential_type == do_potential_mix_cl_trunc) THEN
         dr_ab = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSE IF (potential_parameter%potential_type == do_potential_coulomb) THEN
         dr_ab = 1000000.0_dp
      END IF

      !Looping over the pgfs
      DO ipgf = 1, npgfa
         zeti = zeta(ipgf)
         a_start = (ipgf - 1)*ncoset(la_max)

         DO jpgf = 1, npgfb
            zetj = zetb(jpgf)
            b_start = (jpgf - 1)*ncoset(lb_max)

            !screening
            IF (rpgfa(ipgf) + rpgfb(jpgf) + dr_ab < dab) CYCLE

            DO li = la_min, la_max
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)

                  CALL cp_libint_engine_evaluate(context%engine_2c, [lj, li], &
                                                 RESHAPE([rb, ra], [3, 2]), [1, 1], &
                                                 RESHAPE([zetj, zeti], [1, 2]), &
                                                 RESHAPE([1.0_dp, 1.0_dp], [1, 2]), potential_parameter, p_work)

                  DO i_deriv = 1, 3
                     DO j = 1, ncob
                        p1 = (j - 1)*ncoa
                        DO i = 1, ncoa
                           p2 = p1 + i
                           der_ab(a_offset + i, b_offset + j, i_deriv) = p_work(p2, i_deriv + 3)
                        END DO
                     END DO
                  END DO
               END DO
            END DO

         END DO
      END DO

   END SUBROUTINE compute_coulomb_2c_derivs

END MODULE libint_2c_3c
