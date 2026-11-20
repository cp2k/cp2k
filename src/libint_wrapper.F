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
! Borrow the scalar kernel scratch for each group without per-record profiling or accumulation state.
#if LIBINT_CONTRACTED_INTS && LIBINT2_MAX_VECLEN == 1 && !LIBINT2_FLOP_COUNT && !defined(LIBINT_FLOP_COUNT)
#if !defined(LIBINT_PROFILE) && !defined(LIBINT2_PROFILE) && !defined(LIBINT_ACCUM_INTS) && !defined(LIBINT2_ACCUM_INTS)
#define CP2K_LIBINT_IMAGE_CONTRACTION 1
#endif
#endif
#endif

! maximum angular momentum to be supported in CP2K-LIBINT interface
   #:set libint_max_am_supported = 8

   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_F_POINTER, &
                            C_NULL_PTR, &
                            C_NULL_FUNPTR, &
                            C_FUNPTR
   USE kinds, ONLY: dp
#if(__LIBINT)
   USE libint_f, ONLY: &
      libint2_build_eri, libint2_build_eri1, libint2_cleanup_eri, &
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
   PUBLIC :: get_eri_permutation, reorder_eri_block
   PUBLIC :: libint_can_contract_images, cp_libint_image_workspace_type, &
             cp_libint_get_image_eris, cp_libint_get_image_derivs

#if CP2K_LIBINT_IMAGE_CONTRACTION
   LOGICAL, PARAMETER :: libint_can_contract_images = .TRUE.
#else
   LOGICAL, PARAMETER :: libint_can_contract_images = .FALSE.
#endif
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

!> \brief Storage for one group of native records. Its leading record borrows the kernel scratch.
   TYPE :: cp_libint_image_workspace_type
      PRIVATE
#if(__LIBINT)
      TYPE(libint_t), ALLOCATABLE :: records(:)
#else
      INTEGER :: unused = -1
#endif
   END TYPE cp_libint_image_workspace_type

#if(__LIBINT)
   INTERFACE
      ! Keep the indirect kernel call on the C side of the interoperability boundary.
      SUBROUTINE cp_libint_build(build, libint) BIND(C, NAME="cp2k_libint_build")
         IMPORT :: C_FUNPTR, libint_t
         TYPE(C_FUNPTR), VALUE :: build
         TYPE(libint_t), DIMENSION(*) :: libint
      END SUBROUTINE cp_libint_build
   END INTERFACE

   INTERFACE build_eri_kernel
      MODULE PROCEDURE build_eri_kernel_2c, build_eri_kernel_3c, build_eri_kernel_4c
   END INTERFACE build_eri_kernel
#endif

CONTAINS

! **************************************************************************************************
!> \brief Select the canonical center order required by generated Libint kernels.
!> \param la original first-center angular momentum
!> \param lb original second-center angular momentum
!> \param lc original third-center angular momentum
!> \param ld original fourth-center angular momentum
!> \param permutation original center index for each canonical center
! **************************************************************************************************
   PURE SUBROUTINE get_eri_permutation(la, lb, lc, ld, permutation)
      INTEGER, INTENT(IN)                                :: la, lb, lc, ld
      INTEGER, INTENT(OUT)                               :: permutation(4)

      permutation = [1, 2, 3, 4]
      IF (la < lb) permutation(1:2) = [2, 1]
      IF (lc < ld) permutation(3:4) = [4, 3]
      IF (la + lb > lc + ld) permutation = permutation([3, 4, 1, 2])
   END SUBROUTINE get_eri_permutation

! **************************************************************************************************
!> \brief Restore original Cartesian axes and derivative center labels.
!> \param source canonical block with the fourth Cartesian index varying fastest
!> \param destination original block with the first Cartesian index varying fastest
!> \param ncart Cartesian shell sizes in original center order
!> \param permutation original center index for each canonical center
!> \param ncomponents one integral component or twelve center derivatives
!> \param nplanes number of independent planes
! **************************************************************************************************
   PURE SUBROUTINE reorder_eri_block(source, destination, ncart, permutation, ncomponents, nplanes)
      INTEGER, INTENT(IN)                                :: ncart(4), permutation(4), ncomponents, &
                                                            nplanes
      REAL(dp), INTENT(OUT) :: destination(PRODUCT(ncart), ncomponents, nplanes)
      REAL(dp), INTENT(IN) :: source(PRODUCT(ncart), ncomponents, nplanes)

      INTEGER                                            :: component, extent(4), i, isource, j, k, &
                                                            l, output_component, p1, p2, p3, p4, &
                                                            plane, stride(4)

      stride = [1, ncart(1), ncart(1)*ncart(2), ncart(1)*ncart(2)*ncart(3)]
      stride = stride(permutation)
      extent = ncart(permutation)
      ! Read the Libint block contiguously and restore the original Fortran layout with strides.
      DO plane = 1, nplanes
         DO component = 1, ncomponents
            output_component = component
            IF (ncomponents == 12) THEN
               output_component = 3*(permutation((component - 1)/3 + 1) - 1) + MOD(component - 1, 3) + 1
            END IF
            isource = 0
            p1 = 1
            DO i = 1, extent(1)
               p2 = p1
               DO j = 1, extent(2)
                  p3 = p2
                  DO k = 1, extent(3)
                     p4 = p3
                     DO l = 1, extent(4)
                        isource = isource + 1
                        destination(p4, output_component, plane) = source(isource, component, plane)
                        p4 = p4 + stride(4)
                     END DO
                     p3 = p3 + stride(3)
                  END DO
                  p2 = p2 + stride(2)
               END DO
               p1 = p1 + stride(1)
            END DO
         END DO
      END DO
   END SUBROUTINE reorder_eri_block

#if(__LIBINT)
   #:for rank in range(2, 5)
! **************************************************************************************************
!> \brief Validate and invoke a generated Libint kernel without indexing outside its task table.
!> \param records native primitive records, starting with the one holding the scratch
!> \param kernels zero-based kernel table for this integral task and derivative order
!> \param angular_momenta indices in the table's own center order
!> \param kernel_name integral task used in diagnostics
! **************************************************************************************************
      SUBROUTINE build_eri_kernel_${rank}$c(records, kernels, angular_momenta, kernel_name)
         TYPE(libint_t), INTENT(INOUT) :: records(*)
         TYPE(C_FUNPTR), INTENT(IN) :: kernels(${", ".join(["0:"]*rank)}$)
         INTEGER, INTENT(IN) :: angular_momenta(${rank}$)
         CHARACTER(LEN=*), INTENT(IN) :: kernel_name

         CHARACTER(LEN=160) :: error_msg
         TYPE(C_FUNPTR) :: build_funptr

         build_funptr = C_NULL_FUNPTR
         IF (ALL(angular_momenta >= LBOUND(kernels)) .AND. ALL(angular_momenta <= UBOUND(kernels))) THEN
            build_funptr = kernels(${", ".join("angular_momenta(" + str(i) + ")" for i in range(1, rank + 1))}$)
         END IF
         IF (.NOT. C_ASSOCIATED(build_funptr)) THEN
            WRITE (error_msg, '(A,${rank}$(I0,1X))') &
               "LIBINT has no "//TRIM(kernel_name)//" kernel for angular momenta: ", angular_momenta
            CPABORT(TRIM(error_msg))
         END IF
         CALL cp_libint_build(build_funptr, records)
      END SUBROUTINE build_eri_kernel_${rank}$c

   #:endfor
#endif

   #! Emit the same primitive packing for the scalar setter and image groups.
   #:def pack_eri_fields(record, grouped)
      INTEGER :: pair_order(2)
      REAL(dp) :: product_centers(3, 2), inverse_exponents(2)

      CPASSERT(m_max <= params%prepared_max_order)

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
         #:if grouped
            ${record}$%AB_x(1) = displacements(1, 1)
            ${record}$%AB_y(1) = displacements(2, 1)
            ${record}$%AB_z(1) = displacements(3, 1)
            ${record}$%CD_x(1) = displacements(1, 2)
            ${record}$%CD_y(1) = displacements(2, 2)
            ${record}$%CD_z(1) = displacements(3, 2)
         #:else
            ${record}$%AB_x(1) = A(1) - B(1)
            ${record}$%AB_y(1) = A(2) - B(2)
            ${record}$%AB_z(1) = A(3) - B(3)
            ${record}$%CD_x(1) = C(1) - D(1)
            ${record}$%CD_y(1) = C(2) - D(2)
            ${record}$%CD_z(1) = C(3) - D(3)
         #:endif

         ${record}$%PA_x(1) = P(1) - A(1)
         ${record}$%PA_y(1) = P(2) - A(2)
         ${record}$%PA_z(1) = P(3) - A(3)

         ${record}$%QC_x(1) = Q(1) - C(1)
         ${record}$%QC_y(1) = Q(2) - C(2)
         ${record}$%QC_z(1) = Q(3) - C(3)

         ${record}$%WP_x(1) = W(1) - P(1)
         ${record}$%WP_y(1) = W(2) - P(2)
         ${record}$%WP_z(1) = W(3) - P(3)

         ${record}$%WQ_x(1) = W(1) - Q(1)
         ${record}$%WQ_y(1) = W(2) - Q(2)
         ${record}$%WQ_z(1) = W(3) - Q(3)

         ${record}$%oo2z(1) = 0.5_dp*ZetaInv
         ${record}$%oo2e(1) = 0.5_dp*EtaInv
         ${record}$%oo2ze(1) = 0.5_dp*ZetapEtaInv
         ${record}$%roz(1) = Rho*ZetaInv
         ${record}$%roe(1) = Rho*EtaInv

         #:for m_max in range(0, 4*libint_max_am_supported)
#if 4*LIBINT2_MAX_AM_eri > ${m_max}$ - 1

            IF (${m_max}$ <= m_max) THEN
               ${record}$%f_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_${m_max}$ (1) &
                  = F(${m_max}$+1)
            END IF
#endif
         #:endfor

      END ASSOCIATE
   #:enddef

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
      $:pack_eri_fields("libint%prv(1)", False)
#else
      MARK_USED(libint)
      MARK_USED(params)
      MARK_USED(permutation)
      MARK_USED(m_max)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif
   END SUBROUTINE cp_libint_set_params_eri

#if CP2K_LIBINT_IMAGE_CONTRACTION
! **************************************************************************************************
!> \brief Pack one primitive into a native Libint record.
!> \param record native primitive storage
!> \param params primitive data in original center order
!> \param permutation original center index for each backend center
!> \param m_max highest required radial order
!> \param displacements shared pair displacements in backend order
! **************************************************************************************************
   SUBROUTINE set_eri_primitive(record, params, permutation, m_max, displacements)
      TYPE(libint_t), INTENT(INOUT) :: record
      TYPE(eri_primitive_data_type), INTENT(IN) :: params
      INTEGER, INTENT(IN) :: permutation(4), m_max
      REAL(dp), INTENT(IN) :: displacements(3, 2)

      $:pack_eri_fields("record", True)
   END SUBROUTINE set_eri_primitive
#endif

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

      CALL cp_libint_set_params_eri(libint, params, permutation, m_max)
#if(__LIBINT)
      CALL set_eri_derivative(libint%prv(1), params, permutation)
#endif
   END SUBROUTINE cp_libint_set_params_eri_deriv

#if(__LIBINT)
! **************************************************************************************************
!> \brief Pack the remaining first-derivative recurrence fields.
!> \param record native primitive storage
!> \param params primitive data in original center order
!> \param permutation original center index for each backend center
! **************************************************************************************************
   SUBROUTINE set_eri_derivative(record, params, permutation)
      TYPE(libint_t), INTENT(INOUT) :: record
      TYPE(eri_primitive_data_type), INTENT(IN) :: params
      INTEGER, INTENT(IN) :: permutation(4)
      REAL(dp) :: P(3), exponents(4), gammap, gammaq, gammapq, rhop, rhoq

      P = params%P
      IF (permutation(1) > 2) P = params%Q
      exponents = params%exponents(permutation)
      ASSOCIATE (B => params%centers(:, permutation(2)), &
                 Zeta_A => exponents(1), Zeta_B => exponents(2), &
                 Zeta_C => exponents(3), Zeta_D => exponents(4))
         record%PB_x(1) = P(1) - B(1)
         record%PB_y(1) = P(2) - B(2)
         record%PB_z(1) = P(3) - B(3)

         record%two_alpha0_bra(1) = 2.0_dp*Zeta_A
         record%two_alpha0_ket(1) = 2.0_dp*Zeta_B
         record%two_alpha1_ket(1) = 2.0_dp*Zeta_D

         gammap = Zeta_A + Zeta_B
         gammaq = Zeta_C + Zeta_D
         gammapq = gammap*gammaq/(gammap + gammaq)
         record%alpha1_rho_over_zeta2(1) = Zeta_A*gammapq/(gammap*gammap)
         record%alpha2_rho_over_zeta2(1) = Zeta_B*gammapq/(gammap*gammap)
         record%alpha4_rho_over_eta2(1) = Zeta_D*gammapq/(gammaq*gammaq)
         record%alpha1_over_zetapluseta(1) = Zeta_A/(gammap + gammaq)
         record%alpha2_over_zetapluseta(1) = Zeta_B/(gammap + gammaq)
         record%alpha4_over_zetapluseta(1) = Zeta_D/(gammap + gammaq)

         rhop = Zeta_A*Zeta_B/gammap
         rhoq = Zeta_C*Zeta_D/gammaq
         record%rho12_over_alpha1(1) = rhop/Zeta_A

         record%rho34_over_alpha3(1) = rhoq/Zeta_C

      END ASSOCIATE
   END SUBROUTINE set_eri_derivative
#endif
! **************************************************************************************************
!> \brief ...
!> \param n_d ...
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work ...
! **************************************************************************************************
   SUBROUTINE cp_libint_get_eris(n_d, n_c, n_b, n_a, lib, p_work)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

#if(__LIBINT)
      CALL build_eri_kernel(lib%prv, libint2_build_eri, [n_d, n_c, n_b, n_a], "4-center ERI")

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=[nco(n_a)*nco(n_b)*nco(n_c)*nco(n_d)])
#else
      MARK_USED(n_d)
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
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
! **************************************************************************************************
   SUBROUTINE cp_libint_get_3eris(n_c, n_b, n_a, lib, p_work)
      INTEGER, INTENT(IN)                                :: n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

#if(__LIBINT)
      CALL build_eri_kernel(lib%prv, libint2_build_3eri, [n_c, n_b, n_a], "3-center ERI")

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=[nco(n_a)*nco(n_b)*nco(n_c)])
#else
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_3eris

! **************************************************************************************************
!> \brief Copy independent 3-center derivatives into caller-owned compact scratch.
!> \param n_c ...
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work leading Cartesian block rows and 6 columns receive Cartesian derivative triplets
! **************************************************************************************************
   SUBROUTINE cp_libint_get_3eri_derivs(n_c, n_b, n_a, lib, p_work)
      INTEGER, INTENT(IN)                                :: n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:, :), INTENT(OUT)              :: p_work

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work_tmp
      INTEGER :: i, ncart

      ncart = nco(n_a)*nco(n_b)*nco(n_c)
      CPASSERT(SIZE(p_work, 1) >= ncart)
      CPASSERT(SIZE(p_work, 2) >= 6)

      CALL build_eri_kernel(lib%prv, libint2_build_3eri1, [n_c, n_b, n_a], "3-center ERI derivative")

      !Derivatives 1-3 can be obtained using translational invariance
      DO i = 4, 9
         NULLIFY (p_work_tmp)
         CALL C_F_POINTER(lib%prv(1)%targets(i), p_work_tmp, SHAPE=[ncart])
         p_work(1:ncart, i - 3) = p_work_tmp
      END DO
#else
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_3eri_derivs

! **************************************************************************************************
!> \brief Copy independent 2-center derivatives into caller-owned compact scratch.
!> \param n_b ...
!> \param n_a ...
!> \param lib ...
!> \param p_work leading Cartesian block rows and 3 columns receive Cartesian derivative triplets
! **************************************************************************************************
   SUBROUTINE cp_libint_get_2eri_derivs(n_b, n_a, lib, p_work)
      INTEGER, INTENT(IN)                                :: n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:, :), INTENT(OUT)              :: p_work

#if(__LIBINT)
      REAL(dp), DIMENSION(:), POINTER                    :: p_work_tmp
      INTEGER :: i, ncart

      ncart = nco(n_a)*nco(n_b)
      CPASSERT(SIZE(p_work, 1) >= ncart)
      CPASSERT(SIZE(p_work, 2) >= 3)

      CALL build_eri_kernel(lib%prv, libint2_build_2eri1, [n_b, n_a], "2-center ERI derivative")

      !Derivatives 1-3 can be obtained using translational invariance
      DO i = 4, 6
         NULLIFY (p_work_tmp)
         CALL C_F_POINTER(lib%prv(1)%targets(i), p_work_tmp, SHAPE=[ncart])
         p_work(1:ncart, i - 3) = p_work_tmp
      END DO
#else
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
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
! **************************************************************************************************
   SUBROUTINE cp_libint_get_2eris(n_b, n_a, lib, p_work)
      INTEGER, INTENT(IN)                                :: n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(:), POINTER                    :: p_work

#if(__LIBINT)
      CALL build_eri_kernel(lib%prv, libint2_build_2eri, [n_b, n_a], "2-center ERI")

      CALL C_F_POINTER(lib%prv(1)%targets(1), p_work, SHAPE=[nco(n_a)*nco(n_b)])
#else
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(p_work)
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
! **************************************************************************************************
   SUBROUTINE cp_libint_get_derivs(n_d, n_c, n_b, n_a, lib, work_forces)
      INTEGER, INTENT(IN)                                :: n_d, n_c, n_b, n_a
      TYPE(cp_libint_t)                                  :: lib
      REAL(dp), DIMENSION(nco(n_a)*nco(n_b)*nco(n_c)*nco &
                          (n_d), 12)                                      :: work_forces

#if(__LIBINT)
      CALL build_eri_kernel(lib%prv, libint2_build_eri1, [n_d, n_c, n_b, n_a], "4-center ERI derivative")
      CALL copy_eri_derivatives(lib%prv(1), SIZE(work_forces, 1), work_forces)
#else
      MARK_USED(n_d)
      MARK_USED(n_c)
      MARK_USED(n_b)
      MARK_USED(n_a)
      MARK_USED(lib)
      MARK_USED(work_forces)
      CPABORT("This CP2K executable has not been linked against the required library libint.")
#endif

   END SUBROUTINE cp_libint_get_derivs

#if(__LIBINT)
! **************************************************************************************************
!> \brief Copy generated first derivatives and restore the omitted center.
!> \param record first native record of the evaluated batch
!> \param ncart number of Cartesian integrals for each component
!> \param work_forces derivative blocks in canonical center order
! **************************************************************************************************
   SUBROUTINE copy_eri_derivatives(record, ncart, work_forces)
      TYPE(libint_t), INTENT(IN) :: record
      INTEGER, INTENT(IN) :: ncart
      REAL(dp), INTENT(OUT) :: work_forces(ncart, 12)
      REAL(dp), POINTER :: p_work(:)
      INTEGER :: i, k

      DO k = 1, 12
         IF (k == 4 .OR. k == 5 .OR. k == 6) CYCLE
         CALL C_F_POINTER(record%targets(k), p_work, SHAPE=[ncart])
         DO i = 1, ncart
            work_forces(i, k) = p_work(i)
         END DO
      END DO
      ! Libint omits the canonical second-center derivatives; recover them by translation invariance.
      DO k = 4, 6
         DO i = 1, ncart
            work_forces(i, k) = -1.0_dp*(work_forces(i, k - 3) + &
                                         work_forces(i, k + 3) + work_forces(i, k + 6))
         END DO
      END DO
   END SUBROUTINE copy_eri_derivatives
#endif

! **************************************************************************************************
!> \brief Sum a group of image quartets with common pair displacements.
!> \param lib workspace owning the generated kernel scratch
!> \param images reusable storage for the native records
!> \param angular_momenta shell angular momenta in backend order
!> \param params image-resolved primitives in original center order
!> \param permutation original center index for each backend center
!> \param displacements shared pair displacements in original center order
!> \param p_work grouped Cartesian integrals, valid until the workspace is reused
! **************************************************************************************************
   SUBROUTINE cp_libint_get_image_eris(lib, images, angular_momenta, params, permutation, displacements, p_work)
      TYPE(cp_libint_t), INTENT(IN) :: lib
      TYPE(cp_libint_image_workspace_type), INTENT(INOUT) :: images
      INTEGER, INTENT(IN) :: angular_momenta(4), permutation(4)
      TYPE(eri_primitive_data_type), INTENT(IN) :: params(:)
      REAL(dp), INTENT(IN) :: displacements(3, 2)
      REAL(dp), POINTER :: p_work(:)

#if CP2K_LIBINT_IMAGE_CONTRACTION
      CALL pack_image_group(lib, images, params, permutation, SUM(angular_momenta), displacements, .FALSE.)
      CALL build_eri_kernel(images%records, libint2_build_eri, angular_momenta(4:1:-1), "4-center ERI")
      CALL C_F_POINTER(images%records(1)%targets(1), p_work, SHAPE=[PRODUCT(nco(angular_momenta))])
#else
      MARK_USED(lib)
      MARK_USED(images)
      MARK_USED(angular_momenta)
      MARK_USED(params)
      MARK_USED(permutation)
      MARK_USED(displacements)
      MARK_USED(p_work)
      CPABORT("This Libint build does not support image contraction.")
#endif
   END SUBROUTINE cp_libint_get_image_eris

! **************************************************************************************************
!> \brief Sum unweighted image derivatives with common pair displacements.
!> \param lib workspace owning the generated kernel scratch
!> \param images reusable storage for the native records
!> \param angular_momenta shell angular momenta in backend order
!> \param params image-resolved primitives in original center order
!> \param permutation original center index for each backend center
!> \param displacements shared pair displacements in original center order
!> \param work_forces grouped derivative blocks in canonical center order
! **************************************************************************************************
   SUBROUTINE cp_libint_get_image_derivs(lib, images, angular_momenta, params, permutation, displacements, work_forces)
      TYPE(cp_libint_t), INTENT(IN) :: lib
      TYPE(cp_libint_image_workspace_type), INTENT(INOUT) :: images
      INTEGER, INTENT(IN) :: angular_momenta(4), permutation(4)
      TYPE(eri_primitive_data_type), INTENT(IN) :: params(:)
      REAL(dp), INTENT(IN) :: displacements(3, 2)
      REAL(dp), INTENT(OUT) :: work_forces(PRODUCT(nco(angular_momenta)), 12)

#if CP2K_LIBINT_IMAGE_CONTRACTION
      CALL pack_image_group(lib, images, params, permutation, SUM(angular_momenta) + 1, displacements, .TRUE.)
      CALL build_eri_kernel(images%records, libint2_build_eri1, angular_momenta(4:1:-1), "4-center ERI derivative")
      CALL copy_eri_derivatives(images%records(1), SIZE(work_forces, 1), work_forces)
#else
      MARK_USED(lib)
      MARK_USED(images)
      MARK_USED(angular_momenta)
      MARK_USED(params)
      MARK_USED(permutation)
      MARK_USED(displacements)
      MARK_USED(work_forces)
      CPABORT("This Libint build does not support image contraction.")
#endif
   END SUBROUTINE cp_libint_get_image_derivs

#if CP2K_LIBINT_IMAGE_CONTRACTION
! **************************************************************************************************
!> \brief Prepare the native records of a group without allocating another kernel stack.
!> \param lib existing Libint workspace
!> \param images reusable storage for the native records
!> \param params image-resolved primitives in original center order
!> \param permutation original center index for each backend center
!> \param m_max highest required radial order
!> \param displacements shared pair displacements in original center order
!> \param derivative whether to populate first-derivative fields
! **************************************************************************************************
   SUBROUTINE pack_image_group(lib, images, params, permutation, m_max, displacements, derivative)
      TYPE(cp_libint_t), INTENT(IN) :: lib
      TYPE(cp_libint_image_workspace_type), INTENT(INOUT) :: images
      TYPE(eri_primitive_data_type), INTENT(IN) :: params(:)
      INTEGER, INTENT(IN) :: permutation(4), m_max
      REAL(dp), INTENT(IN) :: displacements(3, 2)
      LOGICAL, INTENT(IN) :: derivative

      INTEGER :: i, n
      REAL(dp) :: backend_displacements(3, 2)

      n = SIZE(params)
      IF (ALLOCATED(images%records)) THEN
         IF (SIZE(images%records) < n) DEALLOCATE (images%records)
      END IF
      IF (.NOT. ALLOCATED(images%records)) ALLOCATE (images%records(n))
      backend_displacements(:, 1) = displacements(:, (permutation(1) + 1)/2)
      backend_displacements(:, 2) = displacements(:, (permutation(3) + 1)/2)
      IF (permutation(1) > permutation(2)) backend_displacements(:, 1) = -backend_displacements(:, 1)
      IF (permutation(3) > permutation(4)) backend_displacements(:, 2) = -backend_displacements(:, 2)
      DO i = 1, n
         CALL set_eri_primitive(images%records(i), params(i), permutation, m_max, backend_displacements)
         IF (derivative) CALL set_eri_derivative(images%records(i), params(i), permutation)
      END DO
      ! Only the leading record uses the stack and owns the contraction depth.
      images%records(1)%stack = lib%prv(1)%stack
      images%records(1)%vstack = lib%prv(1)%vstack
      images%records(1)%contrdepth = n
   END SUBROUTINE pack_image_group
#endif

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
