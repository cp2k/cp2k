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
                                              do_potential_id,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: active_integral_library,&
                                              coulomb_operator_type,&
                                              cutoff_screen_factor,&
                                              library_libint,&
                                              library_native
   USE kinds,                           ONLY: dp
   USE libint_wrapper,                  ONLY: &
        cp_libint_cleanup_2eri, cp_libint_cleanup_2eri1, cp_libint_cleanup_3eri, &
        cp_libint_cleanup_3eri1, cp_libint_get_2eri_derivs, cp_libint_get_2eris, &
        cp_libint_get_3eri_derivs, cp_libint_get_3eris, cp_libint_init_2eri, cp_libint_init_2eri1, &
        cp_libint_init_3eri, cp_libint_init_3eri1, cp_libint_set_contrdepth, &
        cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, cp_libint_t, &
        eri_primitive_data_type, libint_max_am, prim_data_f_size
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
   USE r12_fundamentals,                ONLY: compute_r12_fundamentals
   USE t_c_g0,                          ONLY: get_lmax_init
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
      TYPE(cp_libint_t), ALLOCATABLE :: lib_2c, lib_3c
      INTEGER :: library_2c = library_native, library_3c = library_native
      INTEGER :: derivative_order_2c = -1, derivative_order_3c = -1
      REAL(dp), ALLOCATABLE :: derivative_work_2c(:, :), derivative_work_3c(:, :)

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
            ALLOCATE (context%lib_2c)
            IF (order == 0) THEN
               CALL cp_libint_init_2eri(context%lib_2c, max_am_2c)
            ELSE
               CALL cp_libint_init_2eri1(context%lib_2c, max_am_2c)
            END IF
            CALL cp_libint_set_contrdepth(context%lib_2c, 1)
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
            ALLOCATE (context%lib_3c)
            IF (order == 0) THEN
               CALL cp_libint_init_3eri(context%lib_3c, max_am_3c)
            ELSE
               CALL cp_libint_init_3eri1(context%lib_3c, max_am_3c)
            END IF
            CALL cp_libint_set_contrdepth(context%lib_3c, 1)
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

      IF (context%library_2c == library_libint) THEN
         IF (context%derivative_order_2c == 0) CALL cp_libint_cleanup_2eri(context%lib_2c)
         IF (context%derivative_order_2c == 1) CALL cp_libint_cleanup_2eri1(context%lib_2c)
      END IF
      IF (context%library_3c == library_libint) THEN
         IF (context%derivative_order_3c == 0) CALL cp_libint_cleanup_3eri(context%lib_3c)
         IF (context%derivative_order_3c == 1) CALL cp_libint_cleanup_3eri1(context%lib_3c)
      END IF
      IF (ALLOCATED(context%lib_2c)) DEALLOCATE (context%lib_2c)
      IF (ALLOCATED(context%lib_3c)) DEALLOCATE (context%lib_3c)
      IF (ALLOCATED(context%derivative_work_2c)) DEALLOCATE (context%derivative_work_2c)
      IF (ALLOCATED(context%derivative_work_3c)) DEALLOCATE (context%derivative_work_3c)
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
!> \note Prepare shared radial tables before entering a threaded evaluation region.
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
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      TYPE(eri_primitive_data_type)                      :: params

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

      !Note: we want to compute all possible integrals based on the 3-centers (ab|c) before
      !      having to switch to (ba|c) (or the other way around) due to angular momenta in libint
      !      For a triplet of centers (k|ji), we can only compute integrals for which lj >= li

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

               !start with all the (c|ba) integrals (standard order) and keep to lb >= la
               CALL prepare_3c_primitive(ra, rb, rc, zeti, zetj, zetk, &
                                         la_max + lb_max + lc_max, potential_parameter, params)
               IF (la_min <= lb_max) THEN
                  CALL cp_libint_set_params_eri(context%lib_3c, params, [1, 2, 3, 4], params%prepared_max_order)
               END IF

               DO li = la_min, la_max
                  a_offset = a_start + ncoset(li - 1)
                  ncoa = nco(li)
                  DO lj = MAX(li, lb_min), lb_max
                     b_offset = b_start + ncoset(lj - 1)
                     ncob = nco(lj)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        CALL cp_libint_get_3eris(li, lj, lk, context%lib_3c, p_work)

                        IF (PRESENT(int_abc_ext)) int_abc_ext = MAX(int_abc_ext, MAXVAL(ABS(p_work)))
                        DO k = 1, ncoc
                           p1 = (k - 1)*ncob
                           DO j = 1, ncob
                              p2 = (p1 + j - 1)*ncoa
                              DO i = 1, ncoa
                                 p3 = p2 + i
                                 int_abc(a_offset + i, b_offset + j, c_offset + k) = p_work(p3)
                              END DO
                           END DO
                        END DO

                     END DO !lk
                  END DO !lj
               END DO !li

               !swap centers 3 and 4 to compute (c|ab) with lb < la
               IF (lb_min < la_max) THEN
                  CALL cp_libint_set_params_eri(context%lib_3c, params, [1, 2, 4, 3], params%prepared_max_order)
               END IF

               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)
                  DO li = MAX(lj + 1, la_min), la_max
                     a_offset = a_start + ncoset(li - 1)
                     ncoa = nco(li)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        CALL cp_libint_get_3eris(lj, li, lk, context%lib_3c, p_work)

                        IF (PRESENT(int_abc_ext)) int_abc_ext = MAX(int_abc_ext, MAXVAL(ABS(p_work)))
                        DO k = 1, ncoc
                           p1 = (k - 1)*ncoa
                           DO i = 1, ncoa
                              p2 = (p1 + i - 1)*ncob
                              DO j = 1, ncob
                                 p3 = p2 + j
                                 int_abc(a_offset + i, b_offset + j, c_offset + k) = p_work(p3)
                              END DO
                           END DO
                        END DO

                     END DO !lk
                  END DO !li
               END DO !lj

            END DO !kpgf
         END DO !jpgf
      END DO !ipgf

   END SUBROUTINE compute_coulomb_3c

! **************************************************************************************************
!> \brief Prepare geometry and radial seeds for 3-center integrals at the requested order.
!> \param ri ...
!> \param rj ...
!> \param rk ...
!> \param zeti ...
!> \param zetj ...
!> \param zetk ...
!> \param m_max maximum radial order, including any derivative increment
!> \param potential_parameter radial operator parameters
!> \param params reusable primitive parameter cache
! **************************************************************************************************
   SUBROUTINE prepare_3c_primitive(ri, rj, rk, zeti, zetj, zetk, m_max, potential_parameter, params)
      REAL(dp), INTENT(IN)                               :: ri(3), rj(3), rk(3), zeti, zetj, zetk
      INTEGER, INTENT(IN)                                :: m_max
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      TYPE(eri_primitive_data_type), INTENT(INOUT)       :: params

      REAL(dp) :: fundamentals_work(prim_data_f_size), gammaq, normalization, prefac, S1234, T

      CPASSERT(m_max < prim_data_f_size)

      ! Libint layout: (k, dummy k | j, i).
      params%centers(:, 1) = rk
      params%centers(:, 2) = rk
      params%centers(:, 3) = rj
      params%centers(:, 4) = ri
      params%exponents = [zetk, 0.0_dp, zetj, zeti]
      params%P = rk

      params%prepared_max_order = m_max
      gammaq = zeti + zetj
      params%ZetaInv = 1._dp/zetk; params%EtaInv = 1._dp/gammaq
      params%ZetapEtaInv = 1._dp/(zetk + gammaq)

      params%Q = (zeti*ri + zetj*rj)*params%EtaInv
      params%W = (zetk*rk + gammaq*params%Q)*params%ZetapEtaInv
      params%Rho = zetk*gammaq/(zetk + gammaq)

      params%Fm = 0.0_dp
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_coulomb, do_potential_truncated, do_potential_short, do_potential_mix_cl_trunc)
         IF (potential_parameter%potential_type == do_potential_truncated .OR. &
             potential_parameter%potential_type == do_potential_mix_cl_trunc) THEN
            CPASSERT(get_lmax_init() >= params%prepared_max_order)
         END IF
         T = params%Rho*SUM((params%Q - rk)**2)
         S1234 = EXP(-zeti*zetj*params%EtaInv*SUM((rj - ri)**2))
         CALL compute_r12_fundamentals(params%Fm, normalization, params%prepared_max_order, potential_parameter%potential_type, &
                                       params%Rho, T, fundamentals_work, &
                                       omega=potential_parameter%omega, &
                                       cutoff_radius=potential_parameter%cutoff_radius, &
                                       scale_coulomb=potential_parameter%scale_coulomb, &
                                       scale_longrange=potential_parameter%scale_longrange, &
                                       scale_gaussian=0.0_dp)
         prefac = normalization*SQRT((pi*params%ZetapEtaInv)**3)*S1234
         params%Fm = prefac*params%Fm
      CASE (do_potential_id)
         S1234 = EXP(-zeti*zetj*params%EtaInv*SUM((rj - ri)**2) &
                     - gammaq*zetk*params%ZetapEtaInv*SUM((params%Q - rk)**2))
         prefac = SQRT((pi*params%ZetapEtaInv)**3)*S1234

         params%Fm(:) = prefac
      CASE DEFAULT
         CPABORT("Requested operator NYI")
      END SELECT
   END SUBROUTINE prepare_3c_primitive

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
!> \note Prepare shared radial tables before entering a threaded evaluation region.
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
      TYPE(eri_primitive_data_type)                      :: params

      CPASSERT(context%derivative_order_3c == 1)
      CPASSERT(context%library_3c == library_libint)

      IF (ALLOCATED(context%derivative_work_3c)) THEN
         IF (SIZE(context%derivative_work_3c, 1) < nco(la_max)*nco(lb_max)*nco(lc_max)) DEALLOCATE (context%derivative_work_3c)
      END IF
      IF (.NOT. ALLOCATED(context%derivative_work_3c)) ALLOCATE (context%derivative_work_3c(nco(la_max)*nco(lb_max)*nco(lc_max), 6))

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

      !Note: we want to compute all possible integrals based on the 3-centers (ab|c) before
      !      having to switch to (ba|c) (or the other way around) due to angular momenta in libint
      !      For a triplet of centers (k|ji), we can only compute integrals for which lj >= li

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

               !start with all the (c|ba) integrals (standard order) and keep to lb >= la
               CALL prepare_3c_primitive(ra, rb, rc, zeti, zetj, zetk, &
                                         la_max + lb_max + lc_max + 1, potential_parameter, params)
               IF (la_min <= lb_max) THEN
                  CALL cp_libint_set_params_eri_deriv(context%lib_3c, params, [1, 2, 3, 4], params%prepared_max_order)
               END IF

               DO li = la_min, la_max
                  a_offset = a_start + ncoset(li - 1)
                  ncoa = nco(li)
                  DO lj = MAX(li, lb_min), lb_max
                     b_offset = b_start + ncoset(lj - 1)
                     ncob = nco(lj)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        ncart = ncoa*ncob*ncoc

                        CALL cp_libint_get_3eri_derivs(li, lj, lk, context%lib_3c, context%derivative_work_3c)

                        DO i_deriv = 1, 3
                           IF (do_ext) THEN
                              der_abc_1_ext_prv(i_deriv) = MAX(der_abc_1_ext_prv(i_deriv), &
                                                               MAXVAL(ABS(context%derivative_work_3c(1:ncart, i_deriv + 3))))
                              der_abc_2_ext_prv(i_deriv) = MAX(der_abc_2_ext_prv(i_deriv), &
                                                               MAXVAL(ABS(context%derivative_work_3c(1:ncart, i_deriv))))
                           END IF
                           DO k = 1, ncoc
                              p1 = (k - 1)*ncob
                              DO j = 1, ncob
                                 p2 = (p1 + j - 1)*ncoa
                                 DO i = 1, ncoa
                                    p3 = p2 + i
                                    der_abc_1(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                       context%derivative_work_3c(p3, i_deriv + 3)
                                    der_abc_2(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                       context%derivative_work_3c(p3, i_deriv)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO !lk
                  END DO !lj
               END DO !li

               !swap centers 3 and 4 to compute (c|ab) with lb < la
               IF (lb_min < la_max) THEN
                  CALL cp_libint_set_params_eri_deriv(context%lib_3c, params, [1, 2, 4, 3], params%prepared_max_order)
               END IF

               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)
                  DO li = MAX(lj + 1, la_min), la_max
                     a_offset = a_start + ncoset(li - 1)
                     ncoa = nco(li)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        ncart = ncoa*ncob*ncoc
                        CALL cp_libint_get_3eri_derivs(lj, li, lk, context%lib_3c, context%derivative_work_3c)

                        DO i_deriv = 1, 3
                           IF (do_ext) THEN
                              der_abc_1_ext_prv(i_deriv) = MAX(der_abc_1_ext_prv(i_deriv), &
                                                               MAXVAL(ABS(context%derivative_work_3c(1:ncart, i_deriv))))
                              der_abc_2_ext_prv(i_deriv) = MAX(der_abc_2_ext_prv(i_deriv), &
                                                               MAXVAL(ABS(context%derivative_work_3c(1:ncart, i_deriv + 3))))
                           END IF
                           DO k = 1, ncoc
                              p1 = (k - 1)*ncoa
                              DO i = 1, ncoa
                                 p2 = (p1 + i - 1)*ncob
                                 DO j = 1, ncob
                                    p3 = p2 + j
                                    der_abc_1(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                       context%derivative_work_3c(p3, i_deriv)
                                    der_abc_2(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                       context%derivative_work_3c(p3, i_deriv + 3)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO !lk
                  END DO !li
               END DO !lj

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
!> \note Prepare shared radial tables before entering a threaded evaluation region.
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
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      TYPE(eri_primitive_data_type)                      :: params

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

            CALL prepare_2c_primitive(ra, rb, zeti, zetj, la_max + lb_max, &
                                      potential_parameter, params)
            CALL cp_libint_set_params_eri(context%lib_2c, params, [1, 2, 3, 4], params%prepared_max_order)

            DO li = la_min, la_max
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)

                  CALL cp_libint_get_2eris(li, lj, context%lib_2c, p_work)

                  DO j = 1, ncob
                     p1 = (j - 1)*ncoa
                     DO i = 1, ncoa
                        p2 = p1 + i
                        int_ab(a_offset + i, b_offset + j) = p_work(p2)
                     END DO
                  END DO

               END DO
            END DO

         END DO
      END DO

   END SUBROUTINE compute_coulomb_2c

! **************************************************************************************************
!> \brief Prepare geometry and radial seeds for 2-center integrals at the requested order.
!> \param rj ...
!> \param rk ...
!> \param zetj ...
!> \param zetk ...
!> \param m_max maximum radial order, including any derivative increment
!> \param potential_parameter radial operator parameters
!> \param params reusable primitive parameter cache
! **************************************************************************************************
   SUBROUTINE prepare_2c_primitive(rj, rk, zetj, zetk, m_max, potential_parameter, params)
      REAL(dp), INTENT(IN)                               :: rj(3), rk(3), zetj, zetk
      INTEGER, INTENT(IN)                                :: m_max
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      TYPE(eri_primitive_data_type), INTENT(INOUT)       :: params

      REAL(dp) :: fundamentals_work(prim_data_f_size), normalization, prefac, T

      CPASSERT(m_max < prim_data_f_size)

      ! Libint layout: (k, dummy k | j, dummy j).
      params%centers(:, 1) = rk
      params%centers(:, 2) = rk
      params%centers(:, 3) = rj
      params%centers(:, 4) = rj
      params%exponents = [zetk, 0.0_dp, zetj, 0.0_dp]
      params%P = rk
      params%Q = rj

      params%prepared_max_order = m_max
      params%ZetaInv = 1._dp/zetk; params%EtaInv = 1._dp/zetj
      params%ZetapEtaInv = 1._dp/(zetk + zetj)

      params%W = (zetk*rk + zetj*rj)*params%ZetapEtaInv
      params%Rho = zetk*zetj/(zetk + zetj)

      params%Fm = 0.0_dp
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_coulomb, do_potential_truncated, do_potential_short, do_potential_mix_cl_trunc)
         IF (potential_parameter%potential_type == do_potential_truncated .OR. &
             potential_parameter%potential_type == do_potential_mix_cl_trunc) THEN
            CPASSERT(get_lmax_init() >= params%prepared_max_order)
         END IF
         T = params%Rho*SUM((rj - rk)**2)
         CALL compute_r12_fundamentals(params%Fm, normalization, params%prepared_max_order, potential_parameter%potential_type, &
                                       params%Rho, T, fundamentals_work, &
                                       omega=potential_parameter%omega, &
                                       cutoff_radius=potential_parameter%cutoff_radius, &
                                       scale_coulomb=potential_parameter%scale_coulomb, &
                                       scale_longrange=potential_parameter%scale_longrange, &
                                       scale_gaussian=0.0_dp)
         prefac = normalization*SQRT((pi*params%ZetapEtaInv)**3)
         params%Fm = prefac*params%Fm
      CASE (do_potential_id)

         prefac = SQRT((pi*params%ZetapEtaInv)**3)*EXP(-zetj*zetk*params%ZetapEtaInv*SUM((rk - rj)**2))
         params%Fm(:) = prefac
      CASE DEFAULT
         CPABORT("Requested operator NYI")
      END SELECT
   END SUBROUTINE prepare_2c_primitive

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
!> \note Prepare shared radial tables before entering a threaded evaluation region.
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
      TYPE(eri_primitive_data_type)                      :: params

      CPASSERT(context%derivative_order_2c == 1)
      CPASSERT(context%library_2c == library_libint)

      IF (ALLOCATED(context%derivative_work_2c)) THEN
         IF (SIZE(context%derivative_work_2c, 1) < nco(la_max)*nco(lb_max)) DEALLOCATE (context%derivative_work_2c)
      END IF
      IF (.NOT. ALLOCATED(context%derivative_work_2c)) ALLOCATE (context%derivative_work_2c(nco(la_max)*nco(lb_max), 3))

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

            CALL prepare_2c_primitive(ra, rb, zeti, zetj, la_max + lb_max + 1, &
                                      potential_parameter, params)
            CALL cp_libint_set_params_eri_deriv(context%lib_2c, params, [1, 2, 3, 4], params%prepared_max_order)

            DO li = la_min, la_max
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)

                  CALL cp_libint_get_2eri_derivs(li, lj, context%lib_2c, context%derivative_work_2c)

                  DO i_deriv = 1, 3
                     DO j = 1, ncob
                        p1 = (j - 1)*ncoa
                        DO i = 1, ncoa
                           p2 = p1 + i
                           der_ab(a_offset + i, b_offset + j, i_deriv) = context%derivative_work_2c(p2, i_deriv)
                        END DO
                     END DO
                  END DO
               END DO
            END DO

         END DO
      END DO

   END SUBROUTINE compute_coulomb_2c_derivs

END MODULE libint_2c_3c
