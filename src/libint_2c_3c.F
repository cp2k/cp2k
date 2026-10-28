!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief 2- and 3-center electron repulsion integral routines based on libint2
!>        Currently available operators: Coulomb, Truncated Coulomb, Short Range (erfc), Overlap
!> \author A. Bussy (05.2019)
! **************************************************************************************************

MODULE libint_2c_3c
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type,&
                                              cutoff_screen_factor
   USE kinds,                           ONLY: dp
   USE libint_wrapper,                  ONLY: &
        cp_libint_get_2eri_derivs, cp_libint_get_2eris, cp_libint_get_3eri_derivs, &
        cp_libint_get_3eris, cp_libint_set_params_eri, cp_libint_set_params_eri_deriv, &
        cp_libint_t, eri_primitive_data_type, prim_data_f_size
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset
   USE r12_fundamentals,                ONLY: compute_r12_fundamentals
   USE t_c_g0,                          ONLY: get_lmax_init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_2c_3c'

   PUBLIC :: eri_2center, eri_3center, eri_3center_derivs, eri_2center_derivs

CONTAINS

! **************************************************************************************************
!> \brief Computes the 3-center electron repulsion integrals (ab|c) for a given set of cartesian
!>        gaussian orbitals
!> \param int_abc the integrals as array of cartesian orbitals (allocated before hand)
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
!> \param lib the libint_t object for evaluation (assume that it is initialized outside)
!> \param potential_parameter the info about the potential
!> \param int_abc_ext the extremal value of int_abc, i.e., MAXVAL(ABS(int_abc))
!> \note Prior to calling this routine, the cp_libint_t type passed as argument must be initialized,
!>       the libint library must be static initialized, and in case of truncated Coulomb operator,
!>       the latter must be initialized too
! **************************************************************************************************
   SUBROUTINE eri_3center(int_abc, la_min, la_max, npgfa, zeta, rpgfa, ra, &
                          lb_min, lb_max, npgfb, zetb, rpgfb, rb, &
                          lc_min, lc_max, npgfc, zetc, rpgfc, rc, &
                          dab, dac, dbc, lib, potential_parameter, &
                          int_abc_ext)

      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: int_abc
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
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: int_abc_ext

      INTEGER :: a_mysize(1), a_offset, a_start, b_offset, b_start, c_offset, c_start, i, ipgf, j, &
         jpgf, k, kpgf, li, lj, lk, ncoa, ncob, ncoc, op, p1, p2, p3
      REAL(dp)                                           :: dr_ab, dr_ac, dr_bc, zeti, zetj, zetk
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      TYPE(eri_primitive_data_type)                      :: params

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
                                         la_max + lb_max + lc_max, potential_parameter, params)
               CALL cp_libint_set_params_eri(lib, params, [1, 2, 3, 4], params%prepared_max_order)

               DO li = la_min, la_max
                  a_offset = a_start + ncoset(li - 1)
                  ncoa = nco(li)
                  DO lj = MAX(li, lb_min), lb_max
                     b_offset = b_start + ncoset(lj - 1)
                     ncob = nco(lj)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        a_mysize(1) = ncoa*ncob*ncoc
                        CALL cp_libint_get_3eris(li, lj, lk, lib, p_work, a_mysize)

                        IF (PRESENT(int_abc_ext)) THEN
                           DO k = 1, ncoc
                              p1 = (k - 1)*ncob
                              DO j = 1, ncob
                                 p2 = (p1 + j - 1)*ncoa
                                 DO i = 1, ncoa
                                    p3 = p2 + i
                                    int_abc(a_offset + i, b_offset + j, c_offset + k) = p_work(p3)
                                    int_abc_ext = MAX(int_abc_ext, ABS(p_work(p3)))
                                 END DO
                              END DO
                           END DO
                        ELSE
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
                        END IF

                     END DO !lk
                  END DO !lj
               END DO !li

               !swap centers 3 and 4 to compute (c|ab) with lb < la
               CALL cp_libint_set_params_eri(lib, params, [1, 2, 4, 3], params%prepared_max_order)

               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)
                  DO li = MAX(lj + 1, la_min), la_max
                     a_offset = a_start + ncoset(li - 1)
                     ncoa = nco(li)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        a_mysize(1) = ncoa*ncob*ncoc
                        CALL cp_libint_get_3eris(lj, li, lk, lib, p_work, a_mysize)

                        IF (PRESENT(int_abc_ext)) THEN
                           DO k = 1, ncoc
                              p1 = (k - 1)*ncoa
                              DO i = 1, ncoa
                                 p2 = (p1 + i - 1)*ncob
                                 DO j = 1, ncob
                                    p3 = p2 + j
                                    int_abc(a_offset + i, b_offset + j, c_offset + k) = p_work(p3)
                                    int_abc_ext = MAX(int_abc_ext, ABS(p_work(p3)))
                                 END DO
                              END DO
                           END DO
                        ELSE
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
                        END IF

                     END DO !lk
                  END DO !li
               END DO !lj

            END DO !kpgf
         END DO !jpgf
      END DO !ipgf

   END SUBROUTINE eri_3center

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
!> \param lib the libint_t object for evaluation (assume that it is initialized outside)
!> \param potential_parameter the info about the potential
!> \param der_abc_1_ext the extremal value of der_abc_1, i.e., MAXVAL(ABS(der_abc_1))
!> \param der_abc_2_ext ...
!> \note Prior to calling this routine, the cp_libint_t type passed as argument must be initialized,
!>       the libint library must be static initialized, and in case of truncated Coulomb operator,
!>       the latter must be initialized too. Note that the derivative wrt to the third center
!>       can be obtained via translational invariance
! **************************************************************************************************
   SUBROUTINE eri_3center_derivs(der_abc_1, der_abc_2, &
                                 la_min, la_max, npgfa, zeta, rpgfa, ra, &
                                 lb_min, lb_max, npgfb, zetb, rpgfb, rb, &
                                 lc_min, lc_max, npgfc, zetc, rpgfc, rc, &
                                 dab, dac, dbc, lib, potential_parameter, &
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
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      REAL(dp), DIMENSION(3), INTENT(OUT), OPTIONAL      :: der_abc_1_ext, der_abc_2_ext

      INTEGER :: a_mysize(1), a_offset, a_start, b_offset, b_start, c_offset, c_start, i, i_deriv, &
         ipgf, j, jpgf, k, kpgf, li, lj, lk, ncoa, ncob, ncoc, op, p1, p2, p3
      INTEGER, DIMENSION(3)                              :: permute_1, permute_2
      LOGICAL                                            :: do_ext
      REAL(dp)                                           :: dr_ab, dr_ac, dr_bc, zeti, zetj, zetk
      REAL(dp), DIMENSION(3)                             :: der_abc_1_ext_prv, der_abc_2_ext_prv
      REAL(dp), DIMENSION(:, :), POINTER                 :: p_deriv
      TYPE(eri_primitive_data_type)                      :: params

      NULLIFY (p_deriv)

      permute_1 = [4, 5, 6]
      permute_2 = [7, 8, 9]

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
               CALL cp_libint_set_params_eri_deriv(lib, params, [1, 2, 3, 4], params%prepared_max_order)

               DO li = la_min, la_max
                  a_offset = a_start + ncoset(li - 1)
                  ncoa = nco(li)
                  DO lj = MAX(li, lb_min), lb_max
                     b_offset = b_start + ncoset(lj - 1)
                     ncob = nco(lj)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        a_mysize(1) = ncoa*ncob*ncoc

                        CALL cp_libint_get_3eri_derivs(li, lj, lk, lib, p_deriv, a_mysize)

                        IF (do_ext) THEN
                           DO i_deriv = 1, 3
                              DO k = 1, ncoc
                                 p1 = (k - 1)*ncob
                                 DO j = 1, ncob
                                    p2 = (p1 + j - 1)*ncoa
                                    DO i = 1, ncoa
                                       p3 = p2 + i

                                       der_abc_1(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_2(i_deriv))
                                       der_abc_1_ext_prv(i_deriv) = MAX(der_abc_1_ext_prv(i_deriv), &
                                                                        ABS(p_deriv(p3, permute_2(i_deriv))))

                                       der_abc_2(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_1(i_deriv))
                                       der_abc_2_ext_prv(i_deriv) = MAX(der_abc_2_ext_prv(i_deriv), &
                                                                        ABS(p_deriv(p3, permute_1(i_deriv))))

                                    END DO
                                 END DO
                              END DO
                           END DO
                        ELSE
                           DO i_deriv = 1, 3
                              DO k = 1, ncoc
                                 p1 = (k - 1)*ncob
                                 DO j = 1, ncob
                                    p2 = (p1 + j - 1)*ncoa
                                    DO i = 1, ncoa
                                       p3 = p2 + i

                                       der_abc_1(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_2(i_deriv))

                                       der_abc_2(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_1(i_deriv))
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END IF

                        DEALLOCATE (p_deriv)
                     END DO !lk
                  END DO !lj
               END DO !li

               !swap centers 3 and 4 to compute (c|ab) with lb < la
               CALL cp_libint_set_params_eri_deriv(lib, params, [1, 2, 4, 3], params%prepared_max_order)

               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)
                  DO li = MAX(lj + 1, la_min), la_max
                     a_offset = a_start + ncoset(li - 1)
                     ncoa = nco(li)
                     DO lk = lc_min, lc_max
                        c_offset = c_start + ncoset(lk - 1)
                        ncoc = nco(lk)

                        a_mysize(1) = ncoa*ncob*ncoc
                        CALL cp_libint_get_3eri_derivs(lj, li, lk, lib, p_deriv, a_mysize)

                        IF (do_ext) THEN
                           DO i_deriv = 1, 3
                              DO k = 1, ncoc
                                 p1 = (k - 1)*ncoa
                                 DO i = 1, ncoa
                                    p2 = (p1 + i - 1)*ncob
                                    DO j = 1, ncob
                                       p3 = p2 + j

                                       der_abc_1(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_1(i_deriv))

                                       der_abc_1_ext_prv(i_deriv) = MAX(der_abc_1_ext_prv(i_deriv), &
                                                                        ABS(p_deriv(p3, permute_1(i_deriv))))

                                       der_abc_2(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_2(i_deriv))

                                       der_abc_2_ext_prv(i_deriv) = MAX(der_abc_2_ext_prv(i_deriv), &
                                                                        ABS(p_deriv(p3, permute_2(i_deriv))))
                                    END DO
                                 END DO
                              END DO
                           END DO
                        ELSE
                           DO i_deriv = 1, 3
                              DO k = 1, ncoc
                                 p1 = (k - 1)*ncoa
                                 DO i = 1, ncoa
                                    p2 = (p1 + i - 1)*ncob
                                    DO j = 1, ncob
                                       p3 = p2 + j

                                       der_abc_1(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_1(i_deriv))

                                       der_abc_2(a_offset + i, b_offset + j, c_offset + k, i_deriv) = &
                                          p_deriv(p3, permute_2(i_deriv))
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END IF

                        DEALLOCATE (p_deriv)
                     END DO !lk
                  END DO !li
               END DO !lj

            END DO !kpgf
         END DO !jpgf
      END DO !ipgf

      IF (PRESENT(der_abc_1_ext)) der_abc_1_ext = der_abc_1_ext_prv
      IF (PRESENT(der_abc_2_ext)) der_abc_2_ext = der_abc_2_ext_prv

   END SUBROUTINE eri_3center_derivs

! **************************************************************************************************
!> \brief Computes the 2-center electron repulsion integrals (a|b) for a given set of cartesian
!>        gaussian orbitals
!> \param int_ab the integrals as array of cartesian orbitals (allocated before hand)
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
!> \param lib the libint_t object for evaluation (assume that it is initialized outside)
!> \param potential_parameter the info about the potential
!> \note Prior to calling this routine, the cp_libint_t type passed as argument must be initialized,
!>       the libint library must be static initialized, and in case of truncated Coulomb operator,
!>       the latter must be initialized too
! **************************************************************************************************
   SUBROUTINE eri_2center(int_ab, la_min, la_max, npgfa, zeta, rpgfa, ra, &
                          lb_min, lb_max, npgfb, zetb, rpgfb, rb, &
                          dab, lib, potential_parameter)

      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: int_ab
      INTEGER, INTENT(IN)                                :: la_min, la_max, npgfa
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta, rpgfa
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra
      INTEGER, INTENT(IN)                                :: lb_min, lb_max, npgfb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetb, rpgfb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rb
      REAL(dp), INTENT(IN)                               :: dab
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter

      INTEGER                                            :: a_mysize(1), a_offset, a_start, &
                                                            b_offset, b_start, i, ipgf, j, jpgf, &
                                                            li, lj, ncoa, ncob, p1, p2
      REAL(dp)                                           :: dr_ab, zeti, zetj
      REAL(dp), DIMENSION(:), POINTER                    :: p_work
      TYPE(eri_primitive_data_type)                      :: params

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
            IF (rpgfa(ipgf) + rpgfb(jpgf) + dr_ab < dab) CYCLE

            CALL prepare_2c_primitive(ra, rb, zeti, zetj, la_max + lb_max, &
                                      potential_parameter, params)
            CALL cp_libint_set_params_eri(lib, params, [1, 2, 3, 4], params%prepared_max_order)

            DO li = la_min, la_max
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)

                  a_mysize(1) = ncoa*ncob
                  CALL cp_libint_get_2eris(li, lj, lib, p_work, a_mysize)

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

   END SUBROUTINE eri_2center

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
!> \param lib the libint_t object for evaluation (assume that it is initialized outside)
!> \param potential_parameter the info about the potential
!> \note Prior to calling this routine, the cp_libint_t type passed as argument must be initialized,
!>       the libint library must be static initialized, and in case of truncated Coulomb operator,
!>       the latter must be initialized too
! **************************************************************************************************
   SUBROUTINE eri_2center_derivs(der_ab, la_min, la_max, npgfa, zeta, rpgfa, ra, &
                                 lb_min, lb_max, npgfb, zetb, rpgfb, rb, &
                                 dab, lib, potential_parameter)

      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: der_ab
      INTEGER, INTENT(IN)                                :: la_min, la_max, npgfa
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zeta, rpgfa
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra
      INTEGER, INTENT(IN)                                :: lb_min, lb_max, npgfb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: zetb, rpgfb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rb
      REAL(dp), INTENT(IN)                               :: dab
      TYPE(cp_libint_t), INTENT(INOUT)                   :: lib
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter

      INTEGER                                            :: a_mysize(1), a_offset, a_start, &
                                                            b_offset, b_start, i, i_deriv, ipgf, &
                                                            j, jpgf, li, lj, ncoa, ncob, p1, p2
      INTEGER, DIMENSION(3)                              :: permute
      REAL(dp)                                           :: dr_ab, zeti, zetj
      REAL(dp), DIMENSION(:, :), POINTER                 :: p_deriv
      TYPE(eri_primitive_data_type)                      :: params

      NULLIFY (p_deriv)

      permute = [4, 5, 6]

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
            CALL cp_libint_set_params_eri_deriv(lib, params, [1, 2, 3, 4], params%prepared_max_order)

            DO li = la_min, la_max
               a_offset = a_start + ncoset(li - 1)
               ncoa = nco(li)
               DO lj = lb_min, lb_max
                  b_offset = b_start + ncoset(lj - 1)
                  ncob = nco(lj)

                  a_mysize(1) = ncoa*ncob
                  CALL cp_libint_get_2eri_derivs(li, lj, lib, p_deriv, a_mysize)

                  DO i_deriv = 1, 3
                     DO j = 1, ncob
                        p1 = (j - 1)*ncoa
                        DO i = 1, ncoa
                           p2 = p1 + i
                           der_ab(a_offset + i, b_offset + j, i_deriv) = p_deriv(p2, permute(i_deriv))
                        END DO
                     END DO
                  END DO

                  DEALLOCATE (p_deriv)
               END DO
            END DO

         END DO
      END DO

   END SUBROUTINE eri_2center_derivs

END MODULE libint_2c_3c
