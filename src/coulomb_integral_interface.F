!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common interface for two- and three-center Coulomb integrals.
! **************************************************************************************************
MODULE coulomb_integral_interface
   USE ai_coulomb,                      ONLY: coulomb2,&
                                              coulomb3
   USE integral_library_types,          ONLY: active_integral_library,&
                                              coulomb_operator_type,&
                                              library_libint,&
                                              library_native
   USE kinds,                           ONLY: dp
   USE libint_2c_3c,                    ONLY: eri_2center,&
                                              eri_3center
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_2eri,&
                                              cp_libint_cleanup_3eri,&
                                              cp_libint_init_2eri,&
                                              cp_libint_init_3eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE orbital_pointers,                ONLY: ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: coulomb_integral_cleanup, &
             coulomb_integral_init, &
             coulomb_operator_type, &
             compute_coulomb_2c, &
             compute_coulomb_3c

   TYPE, PUBLIC :: coulomb_integral_context_type
      PRIVATE
      TYPE(cp_libint_t) :: lib_2c
      TYPE(cp_libint_t) :: lib_3c
      LOGICAL           :: initialized_2c = .FALSE.
      LOGICAL           :: initialized_3c = .FALSE.
   END TYPE coulomb_integral_context_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize the engines needed by the requested Coulomb integral operations.
!> \param context operation context containing the selected-library engines
!> \param max_am_2c maximum angular momentum for the two-center operation
!> \param max_am_3c maximum angular momentum for the three-center operation
! **************************************************************************************************
   SUBROUTINE coulomb_integral_init(context, max_am_2c, max_am_3c)
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      INTEGER, INTENT(IN), OPTIONAL                      :: max_am_2c, max_am_3c

      IF (PRESENT(max_am_2c)) THEN
         SELECT CASE (active_integral_library%coulomb2_library)
         CASE (library_libint)
            CALL cp_libint_init_2eri(context%lib_2c, max_am_2c)
            CALL cp_libint_set_contrdepth(context%lib_2c, 1)
         CASE (library_native)
            ! Nothing to do.
         CASE DEFAULT
            CPABORT("Requested two-center Coulomb library is not implemented")
         END SELECT
         context%initialized_2c = .TRUE.
      END IF

      IF (PRESENT(max_am_3c)) THEN
         SELECT CASE (active_integral_library%coulomb3_library)
         CASE (library_libint)
            CALL cp_libint_init_3eri(context%lib_3c, max_am_3c)
            CALL cp_libint_set_contrdepth(context%lib_3c, 1)
         CASE (library_native)
            ! Nothing to do.
         CASE DEFAULT
            CPABORT("Requested three-center Coulomb library is not implemented")
         END SELECT
         context%initialized_3c = .TRUE.
      END IF
   END SUBROUTINE coulomb_integral_init

! **************************************************************************************************
!> \brief Release the engines initialized in a Coulomb integral context.
!> \param context operation context containing the selected-library engines
! **************************************************************************************************
   SUBROUTINE coulomb_integral_cleanup(context)
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context

      IF (context%initialized_2c) THEN
         SELECT CASE (active_integral_library%coulomb2_library)
         CASE (library_libint)
            CALL cp_libint_cleanup_2eri(context%lib_2c)
         CASE (library_native)
            ! Nothing to do.
         CASE DEFAULT
            CPABORT("Requested two-center Coulomb library is not implemented")
         END SELECT
         context%initialized_2c = .FALSE.
      END IF

      IF (context%initialized_3c) THEN
         SELECT CASE (active_integral_library%coulomb3_library)
         CASE (library_libint)
            CALL cp_libint_cleanup_3eri(context%lib_3c)
         CASE (library_native)
            ! Nothing to do.
         CASE DEFAULT
            CPABORT("Requested three-center Coulomb library is not implemented")
         END SELECT
         context%initialized_3c = .FALSE.
      END IF
   END SUBROUTINE coulomb_integral_cleanup

! **************************************************************************************************
!> \brief Evaluate one uncontracted two-center Coulomb integral block.
!> \param context operation context containing the selected-library engine
!> \param la_min minimum angular momentum on center A
!> \param la_max maximum angular momentum on center A
!> \param lb_min minimum angular momentum on center B
!> \param lb_max maximum angular momentum on center B
!> \param npgfa number of primitive functions on center A
!> \param npgfb number of primitive functions on center B
!> \param zeta exponents on center A
!> \param zetb exponents on center B
!> \param rpgfa primitive radii on center A
!> \param rpgfb primitive radii on center B
!> \param ra position of center A
!> \param rb position of center B
!> \param hab uncontracted integral block
!> \param potential_parameter Coulomb operator parameters
! **************************************************************************************************
   SUBROUTINE compute_coulomb_2c(context, la_min, la_max, lb_min, lb_max, npgfa, npgfb, zeta, zetb, rpgfa, rpgfb, &
                                 ra, rb, hab, potential_parameter)
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, npgfa, &
                                                            npgfb
      REAL(dp), DIMENSION(npgfa), INTENT(IN)             :: zeta
      REAL(dp), DIMENSION(npgfb), INTENT(IN)             :: zetb
      REAL(dp), DIMENSION(npgfa), INTENT(IN)             :: rpgfa
      REAL(dp), DIMENSION(npgfb), INTENT(IN)             :: rpgfb
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: hab
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter

      INTEGER                                            :: ncoa, ncob
      REAL(dp)                                           :: rab2
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: f_work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: v_work

      ncoa = npgfa*ncoset(la_max)
      ncob = npgfb*ncoset(lb_max)
      rab2 = (ra(1) - rb(1))**2 + (ra(2) - rb(2))**2 + (ra(3) - rb(3))**2

      SELECT CASE (active_integral_library%coulomb2_library)
      CASE (library_libint)
         CPASSERT(context%initialized_2c)
         CALL eri_2center(hab, la_min, la_max, npgfa, zeta, rpgfa, ra, &
                          lb_min, lb_max, npgfb, zetb, rpgfb, rb, SQRT(rab2), context%lib_2c, potential_parameter)
      CASE (library_native)
         ALLOCATE (f_work(0:la_max + lb_max + 2))
         ALLOCATE (v_work(ncoa, ncob, la_max + lb_max + 1))

         CALL coulomb2(la_max, npgfa, zeta, rpgfa, la_min, lb_max, npgfb, zetb, rpgfb, lb_min, &
                       rb - ra, rab2, hab, v_work, f_work, screening=.FALSE.)

         DEALLOCATE (v_work, f_work)
      CASE DEFAULT
         CPABORT("Requested two-center Coulomb library is not implemented")
      END SELECT
   END SUBROUTINE compute_coulomb_2c

! **************************************************************************************************
!> \brief Evaluate one uncontracted three-center Coulomb integral block.
!> \param context operation context containing the selected-library engine
!> \param la_min minimum angular momentum on center A
!> \param la_max maximum angular momentum on center A
!> \param lb_min minimum angular momentum on center B
!> \param lb_max maximum angular momentum on center B
!> \param lc_min minimum angular momentum on center C
!> \param lc_max maximum angular momentum on center C
!> \param npgfa number of primitive functions on center A
!> \param npgfb number of primitive functions on center B
!> \param npgfc number of primitive functions on center C
!> \param zeta exponents on center A
!> \param zetb exponents on center B
!> \param zetc exponents on center C
!> \param rpgfa primitive radii on center A
!> \param rpgfb primitive radii on center B
!> \param rpgfc primitive radii on center C
!> \param ra position of center A
!> \param rb position of center B
!> \param rc position of center C
!> \param habc uncontracted integral block
!> \param potential_parameter Coulomb operator parameters
! **************************************************************************************************
   SUBROUTINE compute_coulomb_3c(context, la_min, la_max, lb_min, lb_max, lc_min, lc_max, npgfa, npgfb, npgfc, &
                                 zeta, zetb, zetc, rpgfa, rpgfb, rpgfc, ra, rb, rc, habc, potential_parameter)
      TYPE(coulomb_integral_context_type), INTENT(INOUT) :: context
      INTEGER, INTENT(IN)                                :: la_min, la_max, lb_min, lb_max, lc_min, &
                                                            lc_max, npgfa, npgfb, npgfc
      REAL(dp), DIMENSION(npgfa), INTENT(IN)             :: zeta
      REAL(dp), DIMENSION(npgfb), INTENT(IN)             :: zetb
      REAL(dp), DIMENSION(npgfc), INTENT(IN)             :: zetc
      REAL(dp), DIMENSION(npgfa), INTENT(IN)             :: rpgfa
      REAL(dp), DIMENSION(npgfb), INTENT(IN)             :: rpgfb
      REAL(dp), DIMENSION(npgfc), INTENT(IN)             :: rpgfc
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rb, rc
      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: habc
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter

      INTEGER                                            :: kpgf, nc_end, nc_start, ncoa, ncob, ncoc
      REAL(dp)                                           :: rab2, rac2, rbc2
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: f_work, gccc, rpgfa_work, rpgfb_work, &
                                                            rpgfc_work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: vabc
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :, :)       :: v_work

      ncoa = npgfa*ncoset(la_max)
      ncob = npgfb*ncoset(lb_max)
      ncoc = npgfc*ncoset(lc_max)
      rab2 = (rb(1) - ra(1))**2 + (rb(2) - ra(2))**2 + (rb(3) - ra(3))**2
      rac2 = (rc(1) - ra(1))**2 + (rc(2) - ra(2))**2 + (rc(3) - ra(3))**2
      rbc2 = (rc(1) - rb(1))**2 + (rc(2) - rb(2))**2 + (rc(3) - rb(3))**2

      SELECT CASE (active_integral_library%coulomb3_library)
      CASE (library_libint)
         CPASSERT(context%initialized_3c)
         CALL eri_3center(habc, la_min, la_max, npgfa, zeta, rpgfa, ra, &
                          lb_min, lb_max, npgfb, zetb, rpgfb, rb, &
                          lc_min, lc_max, npgfc, zetc, rpgfc, rc, &
                          SQRT(rab2), SQRT(rac2), SQRT(rbc2), context%lib_3c, potential_parameter)
      CASE (library_native)
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

         ! in the RI basis, there is only a single primitive Gaussian
         kpgf = 1

         ! ncoset is indexed by -1,..., ,5
         ! e.g. for pure s: lc_min = lc_max = 0, nc_start = 1
         nc_start = ncoset(lc_min - 1) + 1
         nc_end = ncoset(lc_max)

         CALL coulomb3(la_max, npgfa, zeta(:), rpgfa_work(:), la_min, &
                       lb_max, npgfb, zetb(:), rpgfb_work(:), lb_min, &
                       lc_max, zetc(kpgf), rpgfc_work(kpgf), lc_min, &
                       gccc, rb - ra, rab2, rc - ra, rac2, rbc2, &
                       vabc, habc(:, :, nc_start:nc_end), v_work, f_work)

         DEALLOCATE (v_work, f_work, rpgfa_work, rpgfb_work, rpgfc_work, gccc, vabc)
      CASE DEFAULT
         CPABORT("Requested three-center Coulomb library is not implemented")
      END SELECT
   END SUBROUTINE compute_coulomb_3c

END MODULE coulomb_integral_interface
