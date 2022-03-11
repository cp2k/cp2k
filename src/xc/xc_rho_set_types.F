!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief contains the structure
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE xc_rho_set_types
   USE cp_array_utils,                  ONLY: cp_3d_r_p_type
   USE kinds,                           ONLY: dp
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_cr3d,&
                                              pw_pool_create_pw,&
                                              pw_pool_give_back_cr3d,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_spline_utils,                 ONLY: pw_spline_scale_deriv
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE xc_input_constants,              ONLY: xc_deriv_pw,&
                                              xc_deriv_spline2,&
                                              xc_deriv_spline3,&
                                              xc_rho_no_smooth
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_equal,&
                                              xc_rho_cflags_setall,&
                                              xc_rho_cflags_type
   USE xc_util,                         ONLY: xc_pw_derive,&
                                              xc_pw_laplace,&
                                              xc_pw_smooth
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_rho_set_types'

   PUBLIC :: xc_rho_set_type
   PUBLIC :: xc_rho_set_create, xc_rho_set_release, &
             xc_rho_set_update, xc_rho_set_get

! **************************************************************************************************
!> \brief represent a density, with all the representation and data needed
!>      to perform a functional evaluation
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param id_nr identification number (unique)
!> \param local_bounds the part of the 3d array on which the functional is
!>        computed
!> \param owns which components are owned by this structure (and should
!>        be deallocated
!> \param has which components are present and up to date
!> \param rho the density
!> \param drho the gradient of the density (x,y and z direction)
!> \param norm_drho the norm of the gradient of the density
!> \param rhoa , rhob: spin alpha and beta parts of the density in the LSD case
!> \param drhoa , drhob: gradient of the spin alpha and beta parts of the density
!>        in the LSD case (x,y and z direction)
!> \param norm_drhoa , norm_drhob: norm of the gradient of rhoa and rhob
!> \param rho_ 1_3: rho^(1.0_dp/3.0_dp)
!> \param rhoa_ 1_3, rhob_1_3: rhoa^(1.0_dp/3.0_dp), rhob^(1.0_dp/3.0_dp)
!> \param tau the kinetic (KohnSham) part of rho
!> \param tau_a the kinetic (KohnSham) part of rhoa
!> \param tau_b the kinetic (KohnSham) part of rhob
!> \note
!>      the use of 3d arrays is the result of trying to use only basic
!>      types (to be generic and independent from the method), and
!>      avoiding copies using the actual structure.
!>      only the part defined by local bounds is guaranteed to be present,
!>      and it is the only meaningful part.
!> \par History
!>      11.2003 created [fawzi & thomas]
!>      12.2008 added laplace parts [mguidon]
!> \author fawzi & thomas
! **************************************************************************************************
   TYPE xc_rho_set_type
      INTEGER, DIMENSION(2, 3) :: local_bounds
      REAL(kind=dp) :: rho_cutoff, drho_cutoff, tau_cutoff
      TYPE(xc_rho_cflags_type) :: owns, has
      ! for spin restricted systems
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: rho
      TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drho
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: norm_drho
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: rho_1_3
      REAL(kind=dp), DIMENSION(:, :, :), POINTER :: tau
      ! for UNrestricted systems
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: rhoa, rhob
      TYPE(cp_3d_r_p_type), DIMENSION(3)         :: drhoa, drhob
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: norm_drhoa, norm_drhob
      REAL(kind=dp), DIMENSION(:, :, :), POINTER :: rhoa_1_3, rhob_1_3
      REAL(kind=dp), DIMENSION(:, :, :), POINTER :: tau_a, tau_b
      REAL(kind=dp), DIMENSION(:, :, :), POINTER :: laplace_rho, laplace_rhoa, laplace_rhob
   END TYPE xc_rho_set_type

CONTAINS

! **************************************************************************************************
!> \brief allocates and does (minimal) initialization of a rho_set
!> \param rho_set the structure to allocate
!> \param local_bounds ...
!> \param rho_cutoff ...
!> \param drho_cutoff ...
!> \param tau_cutoff ...
! **************************************************************************************************
   SUBROUTINE xc_rho_set_create(rho_set, local_bounds, rho_cutoff, drho_cutoff, &
                                tau_cutoff)
      TYPE(xc_rho_set_type)                              :: rho_set
      INTEGER, DIMENSION(2, 3), INTENT(in)               :: local_bounds
      REAL(kind=dp), INTENT(in), OPTIONAL                :: rho_cutoff, drho_cutoff, tau_cutoff

      INTEGER                                            :: i

      rho_set%rho_cutoff = EPSILON(0.0_dp)
      IF (PRESENT(rho_cutoff)) rho_set%rho_cutoff = rho_cutoff
      rho_set%drho_cutoff = EPSILON(0.0_dp)
      IF (PRESENT(drho_cutoff)) rho_set%drho_cutoff = drho_cutoff
      rho_set%tau_cutoff = EPSILON(0.0_dp)
      IF (PRESENT(tau_cutoff)) rho_set%tau_cutoff = tau_cutoff
      rho_set%local_bounds = local_bounds
      CALL xc_rho_cflags_setall(rho_set%owns, .TRUE.)
      CALL xc_rho_cflags_setall(rho_set%has, .FALSE.)
      NULLIFY (rho_set%rho)
      DO i = 1, 3
         NULLIFY (rho_set%drho(i)%array)
      END DO
      NULLIFY (rho_set%rho_1_3)
      NULLIFY (rho_set%norm_drho, rho_set%rhoa, rho_set%rhob)
      DO i = 1, 3
         NULLIFY (rho_set%drhoa(i)%array, rho_set%drhob(i)%array)
      END DO
      NULLIFY (rho_set%norm_drhoa, rho_set%norm_drhob, &
               rho_set%rhoa_1_3, rho_set%rhob_1_3, &
               rho_set%tau, rho_set%tau_a, rho_set%tau_b, rho_set%laplace_rho, rho_set%laplace_rhoa, &
               rho_set%laplace_rhob)
   END SUBROUTINE xc_rho_set_create

! **************************************************************************************************
!> \brief releases the given rho_set
!> \param rho_set the structure to release
!> \param pw_pool the plae where to give back the arrays
! **************************************************************************************************
   SUBROUTINE xc_rho_set_release(rho_set, pw_pool)
      TYPE(xc_rho_set_type)                              :: rho_set
      TYPE(pw_pool_type), OPTIONAL, POINTER              :: pw_pool

      INTEGER                                            :: i

      IF (PRESENT(pw_pool)) THEN
         IF (ASSOCIATED(pw_pool)) THEN
            CALL xc_rho_set_clean(rho_set, pw_pool)
         ELSE
            CPABORT("pw_pool must be associated")
         END IF
      END IF

      rho_set%local_bounds(1, :) = -HUGE(0) ! we want to crash...
      rho_set%local_bounds(1, :) = HUGE(0)
      IF (rho_set%owns%rho .AND. ASSOCIATED(rho_set%rho)) THEN
         DEALLOCATE (rho_set%rho)
      END IF
      IF (rho_set%owns%rho_spin) THEN
         IF (ASSOCIATED(rho_set%rhoa)) THEN
            DEALLOCATE (rho_set%rhoa)
         END IF
         IF (ASSOCIATED(rho_set%rhob)) THEN
            DEALLOCATE (rho_set%rhob)
         END IF
      END IF
      IF (rho_set%owns%rho_1_3 .AND. ASSOCIATED(rho_set%rho_1_3)) THEN
         DEALLOCATE (rho_set%rho_1_3)
      END IF
      IF (rho_set%owns%rho_spin) THEN
         IF (ASSOCIATED(rho_set%rhoa_1_3)) THEN
            DEALLOCATE (rho_set%rhoa_1_3)
         END IF
         IF (ASSOCIATED(rho_set%rhob_1_3)) THEN
            DEALLOCATE (rho_set%rhob_1_3)
         END IF
      END IF
      IF (rho_set%owns%drho) THEN
         DO i = 1, 3
            IF (ASSOCIATED(rho_set%drho(i)%array)) THEN
               DEALLOCATE (rho_set%drho(i)%array)
            END IF
         END DO
      END IF
      IF (rho_set%owns%drho_spin) THEN
         DO i = 1, 3
            IF (ASSOCIATED(rho_set%drhoa(i)%array)) THEN
               DEALLOCATE (rho_set%drhoa(i)%array)
            END IF
            IF (ASSOCIATED(rho_set%drhob(i)%array)) THEN
               DEALLOCATE (rho_set%drhob(i)%array)
            END IF
         END DO
      END IF
      IF (rho_set%owns%laplace_rho .AND. ASSOCIATED(rho_set%laplace_rho)) THEN
         DEALLOCATE (rho_set%laplace_rho)
      END IF

      IF (rho_set%owns%norm_drho .AND. ASSOCIATED(rho_set%norm_drho)) THEN
         DEALLOCATE (rho_set%norm_drho)
      END IF
      IF (rho_set%owns%laplace_rho_spin) THEN
         IF (ASSOCIATED(rho_set%laplace_rhoa)) THEN
            DEALLOCATE (rho_set%laplace_rhoa)
         END IF
         IF (ASSOCIATED(rho_set%laplace_rhob)) THEN
            DEALLOCATE (rho_set%laplace_rhob)
         END IF
      END IF

      IF (rho_set%owns%norm_drho_spin) THEN
         IF (ASSOCIATED(rho_set%norm_drhoa)) THEN
            DEALLOCATE (rho_set%norm_drhoa)
         END IF
         IF (ASSOCIATED(rho_set%norm_drhob)) THEN
            DEALLOCATE (rho_set%norm_drhob)
         END IF
      END IF
      IF (rho_set%owns%tau .AND. ASSOCIATED(rho_set%tau)) THEN
         DEALLOCATE (rho_set%tau)
      END IF
      IF (rho_set%owns%tau_spin) THEN
         IF (ASSOCIATED(rho_set%tau_a)) THEN
            DEALLOCATE (rho_set%tau_a)
         END IF
         IF (ASSOCIATED(rho_set%tau_b)) THEN
            DEALLOCATE (rho_set%tau_b)
         END IF
      END IF
   END SUBROUTINE xc_rho_set_release

! **************************************************************************************************
!> \brief returns the various attributes of rho_set
!> \param rho_set the object you want info about
!> \param can_return_null if true the object returned can be null,
!>        if false (the default) it stops with an error if a requested
!>        component is not associated
!> \param rho ...
!> \param drho ...
!> \param norm_drho ...
!> \param rhoa ...
!> \param rhob ...
!> \param norm_drhoa ...
!> \param norm_drhob ...
!> \param rho_1_3 ...
!> \param rhoa_1_3 ...
!> \param rhob_1_3 ...
!> \param laplace_rho ...
!> \param laplace_rhoa ...
!> \param laplace_rhob ...
!> \param drhoa ...
!> \param drhob ...
!> \param rho_cutoff ...
!> \param drho_cutoff ...
!> \param tau_cutoff ...
!> \param tau ...
!> \param tau_a ...
!> \param tau_b ...
!> \param local_bounds ...
! **************************************************************************************************
   SUBROUTINE xc_rho_set_get(rho_set, can_return_null, rho, drho, norm_drho, &
                             rhoa, rhob, norm_drhoa, norm_drhob, rho_1_3, rhoa_1_3, &
                             rhob_1_3, laplace_rho, laplace_rhoa, laplace_rhob, drhoa, drhob, rho_cutoff, &
                             drho_cutoff, tau_cutoff, tau, tau_a, tau_b, local_bounds)
      TYPE(xc_rho_set_type)                              :: rho_set
      LOGICAL, INTENT(in), OPTIONAL                      :: can_return_null
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: rho
      TYPE(cp_3d_r_p_type), DIMENSION(3), OPTIONAL       :: drho
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: norm_drho, rhoa, rhob, norm_drhoa, &
                                                            norm_drhob, rho_1_3, rhoa_1_3, &
                                                            rhob_1_3, laplace_rho, laplace_rhoa, &
                                                            laplace_rhob
      TYPE(cp_3d_r_p_type), DIMENSION(3), OPTIONAL       :: drhoa, drhob
      REAL(kind=dp), INTENT(out), OPTIONAL               :: rho_cutoff, drho_cutoff, tau_cutoff
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: tau, tau_a, tau_b
      INTEGER, DIMENSION(2, 3), INTENT(OUT), OPTIONAL    :: local_bounds

      INTEGER                                            :: i
      LOGICAL                                            :: my_can_return_null

      my_can_return_null = .FALSE.
      IF (PRESENT(can_return_null)) my_can_return_null = can_return_null

      IF (PRESENT(rho)) THEN
         rho => rho_set%rho
         CPASSERT(my_can_return_null .OR. ASSOCIATED(rho))
      END IF
      IF (PRESENT(drho)) THEN
         DO i = 1, 3
            drho(i)%array => rho_set%drho(i)%array
            CPASSERT(my_can_return_null .OR. ASSOCIATED(rho_set%drho(i)%array))
         END DO
      END IF
      IF (PRESENT(norm_drho)) THEN
         norm_drho => rho_set%norm_drho
         CPASSERT(my_can_return_null .OR. ASSOCIATED(norm_drho))
      END IF
      IF (PRESENT(laplace_rho)) THEN
         laplace_rho => rho_set%laplace_rho
         CPASSERT(my_can_return_null .OR. ASSOCIATED(laplace_rho))
      END IF
      IF (PRESENT(rhoa)) THEN
         rhoa => rho_set%rhoa
         CPASSERT(my_can_return_null .OR. ASSOCIATED(rhoa))
      END IF
      IF (PRESENT(rhob)) THEN
         rhob => rho_set%rhob
         CPASSERT(my_can_return_null .OR. ASSOCIATED(rhob))
      END IF
      IF (PRESENT(drhoa)) THEN
         DO i = 1, 3
            drhoa(i)%array => rho_set%drhoa(i)%array
            CPASSERT(my_can_return_null .OR. ASSOCIATED(rho_set%drhoa(i)%array))
         END DO
      END IF
      IF (PRESENT(drhob)) THEN
         DO i = 1, 3
            drhob(i)%array => rho_set%drhob(i)%array
            CPASSERT(my_can_return_null .OR. ASSOCIATED(rho_set%drhob(i)%array))
         END DO
      END IF
      IF (PRESENT(laplace_rhoa)) THEN
         laplace_rhoa => rho_set%laplace_rhoa
         CPASSERT(my_can_return_null .OR. ASSOCIATED(laplace_rhoa))
      END IF
      IF (PRESENT(laplace_rhob)) THEN
         laplace_rhob => rho_set%laplace_rhob
         CPASSERT(my_can_return_null .OR. ASSOCIATED(laplace_rhob))
      END IF
      IF (PRESENT(norm_drhoa)) THEN
         norm_drhoa => rho_set%norm_drhoa
         CPASSERT(my_can_return_null .OR. ASSOCIATED(norm_drhoa))
      END IF
      IF (PRESENT(norm_drhob)) THEN
         norm_drhob => rho_set%norm_drhob
         CPASSERT(my_can_return_null .OR. ASSOCIATED(norm_drhob))
      END IF
      IF (PRESENT(rho_1_3)) THEN
         rho_1_3 => rho_set%rho_1_3
         CPASSERT(my_can_return_null .OR. ASSOCIATED(rho_1_3))
      END IF
      IF (PRESENT(rhoa_1_3)) THEN
         rhoa_1_3 => rho_set%rhoa_1_3
         CPASSERT(my_can_return_null .OR. ASSOCIATED(rhoa_1_3))
      END IF
      IF (PRESENT(rhob_1_3)) THEN
         rhob_1_3 => rho_set%rhob_1_3
         CPASSERT(my_can_return_null .OR. ASSOCIATED(rhob_1_3))
      END IF
      IF (PRESENT(tau)) THEN
         tau => rho_set%tau
         CPASSERT(my_can_return_null .OR. ASSOCIATED(tau))
      END IF
      IF (PRESENT(tau_a)) THEN
         tau_a => rho_set%tau_a
         CPASSERT(my_can_return_null .OR. ASSOCIATED(tau_a))
      END IF
      IF (PRESENT(tau_b)) THEN
         tau_b => rho_set%tau_b
         CPASSERT(my_can_return_null .OR. ASSOCIATED(tau_b))
      END IF
      IF (PRESENT(rho_cutoff)) rho_cutoff = rho_set%rho_cutoff
      IF (PRESENT(drho_cutoff)) drho_cutoff = rho_set%drho_cutoff
      IF (PRESENT(tau_cutoff)) tau_cutoff = rho_set%tau_cutoff
      IF (PRESENT(local_bounds)) local_bounds = rho_set%local_bounds
   END SUBROUTINE xc_rho_set_get

! **************************************************************************************************
!> \brief cleans (releases) most of the data stored in the rho_set giving back
!>      what it can to the pw_pool
!> \param rho_set the rho_set to be cleaned
!> \param pw_pool place to give back 3d arrays,...
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE xc_rho_set_clean(rho_set, pw_pool)
      TYPE(xc_rho_set_type), INTENT(INOUT)               :: rho_set
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      INTEGER                                            :: idir

      IF (rho_set%owns%rho) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%rho, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%rho)
      END IF
      IF (rho_set%owns%rho_1_3) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%rho_1_3, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%rho_1_3)
      END IF
      IF (rho_set%owns%drho) THEN
         DO idir = 1, 3
            CALL pw_pool_give_back_cr3d(pw_pool, rho_set%drho(idir)%array, &
                                        accept_non_compatible=.TRUE.)
         END DO
      ELSE
         DO idir = 1, 3
            NULLIFY (rho_set%drho(idir)%array)
         END DO
      END IF
      IF (rho_set%owns%norm_drho) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%norm_drho, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%norm_drho)
      END IF
      IF (rho_set%owns%laplace_rho) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%laplace_rho, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%laplace_rho)
      END IF
      IF (rho_set%owns%tau) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%tau, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%tau)
      END IF
      IF (rho_set%owns%rho_spin) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%rhoa, &
                                     accept_non_compatible=.TRUE.)
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%rhob, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%rhoa, rho_set%rhob)
      END IF
      IF (rho_set%owns%rho_spin_1_3) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%rhoa_1_3, &
                                     accept_non_compatible=.TRUE.)
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%rhob_1_3, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%rhoa_1_3, rho_set%rhob_1_3)
      END IF
      IF (rho_set%owns%drho_spin) THEN
         DO idir = 1, 3
            CALL pw_pool_give_back_cr3d(pw_pool, rho_set%drhoa(idir)%array, &
                                        accept_non_compatible=.TRUE.)
            CALL pw_pool_give_back_cr3d(pw_pool, rho_set%drhob(idir)%array, &
                                        accept_non_compatible=.TRUE.)
         END DO
      ELSE
         DO idir = 1, 3
            NULLIFY (rho_set%drhoa(idir)%array, rho_set%drhob(idir)%array)
         END DO
      END IF
      IF (rho_set%owns%laplace_rho_spin) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%laplace_rhoa, &
                                     accept_non_compatible=.TRUE.)
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%laplace_rhob, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%laplace_rhoa, rho_set%laplace_rhob)
      END IF
      IF (rho_set%owns%norm_drho_spin) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%norm_drhoa, &
                                     accept_non_compatible=.TRUE.)
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%norm_drhob, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%norm_drhoa, rho_set%norm_drhob)
      END IF
      IF (rho_set%owns%tau_spin) THEN
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%tau_a, &
                                     accept_non_compatible=.TRUE.)
         CALL pw_pool_give_back_cr3d(pw_pool, rho_set%tau_b, &
                                     accept_non_compatible=.TRUE.)
      ELSE
         NULLIFY (rho_set%tau_a, rho_set%tau_b)
      END IF

      CALL xc_rho_cflags_setall(rho_set%has, .FALSE.)
      CALL xc_rho_cflags_setall(rho_set%owns, .FALSE.)

   END SUBROUTINE xc_rho_set_clean

! **************************************************************************************************
!> \brief updates the given rho set with the density given by
!>      rho_r (and rho_g). The rho set will contain the components specified
!>      in needs
!> \param rho_set the rho_set to update
!> \param rho_r the new density (in r space)
!> \param rho_g the new density (in g space, needed for some
!>        derivatives)
!> \param tau ...
!> \param needs the components of rho that are needed
!> \param xc_deriv_method_id ...
!> \param xc_rho_smooth_id ...
!> \param pw_pool pool for the allocation of pw and cr3d
! **************************************************************************************************
   SUBROUTINE xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                                xc_deriv_method_id, xc_rho_smooth_id, pw_pool)
      TYPE(xc_rho_set_type), INTENT(INOUT)               :: rho_set
      TYPE(pw_p_type), DIMENSION(:), INTENT(IN)          :: rho_r
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g
      TYPE(pw_p_type), DIMENSION(:), INTENT(IN)          :: tau
      TYPE(xc_rho_cflags_type), INTENT(in)               :: needs
      INTEGER, INTENT(IN)                                :: xc_deriv_method_id, xc_rho_smooth_id
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      REAL(KIND=dp), PARAMETER                           :: f13 = (1.0_dp/3.0_dp)

      INTEGER                                            :: i, idir, ispin, j, k, nspins
      LOGICAL                                            :: gradient_f, my_rho_g_local, &
                                                            needs_laplace, needs_rho_g
      REAL(kind=dp)                                      :: rho_cutoff
      TYPE(pw_p_type), DIMENSION(2)                      :: laplace_rho_r, my_rho_r
      TYPE(pw_p_type), DIMENSION(3, 2)                   :: drho_r
      TYPE(pw_type), POINTER                             :: my_rho_g, tmp_g

      DO ispin = 1, 2
         NULLIFY (my_rho_r(ispin)%pw)
         DO idir = 1, 3
            NULLIFY (drho_r(idir, ispin)%pw)
         END DO
      END DO
      NULLIFY (tmp_g, my_rho_g)

      IF (ANY(rho_set%local_bounds /= pw_pool%pw_grid%bounds_local)) &
         CPABORT("pw_pool cr3d have different size than expected")
      nspins = SIZE(rho_r)
      rho_set%local_bounds = rho_r(1)%pw%pw_grid%bounds_local
      rho_cutoff = 0.5*rho_set%rho_cutoff

      my_rho_g_local = .FALSE.
      ! some checks
      SELECT CASE (nspins)
      CASE (1)
         CPASSERT(SIZE(rho_r) == 1)
         CPASSERT(ASSOCIATED(rho_r(1)%pw))
         CPASSERT(rho_r(1)%pw%in_use == REALDATA3D)
         CPASSERT(.NOT. needs%rho_spin)
         CPASSERT(.NOT. needs%drho_spin)
         CPASSERT(.NOT. needs%norm_drho_spin)
         CPASSERT(.NOT. needs%rho_spin_1_3)
         CPASSERT(.NOT. needs%tau_spin)
         CPASSERT(.NOT. needs%laplace_rho_spin)
      CASE (2)
         CPASSERT(SIZE(rho_r) == 2)
         CPASSERT(ASSOCIATED(rho_r(1)%pw))
         CPASSERT(ASSOCIATED(rho_r(2)%pw))
         CPASSERT(rho_r(1)%pw%in_use == REALDATA3D)
         CPASSERT(rho_r(2)%pw%in_use == REALDATA3D)
         CPASSERT(.NOT. needs%rho)
         CPASSERT(.NOT. needs%drho)
         CPASSERT(.NOT. needs%rho_1_3)
         CPASSERT(.NOT. needs%tau)
         CPASSERT(.NOT. needs%laplace_rho)
      CASE default
         CPABORT("Unknown number of spin states")
      END SELECT

      CALL xc_rho_set_clean(rho_set, pw_pool=pw_pool)

      needs_laplace = (needs%laplace_rho .OR. needs%laplace_rho_spin)
      gradient_f = (needs%drho_spin .OR. needs%norm_drho_spin .OR. &
                    needs%drho .OR. needs%norm_drho .OR. &
                    needs_laplace)
      needs_rho_g = ((xc_deriv_method_id == xc_deriv_spline3 .OR. &
                      xc_deriv_method_id == xc_deriv_spline2 .OR. &
                      xc_deriv_method_id == xc_deriv_pw)) .AND. (gradient_f .OR. needs_laplace)
      IF ((gradient_f .AND. needs_laplace) .AND. &
          (xc_deriv_method_id /= xc_deriv_pw)) THEN
         CALL cp_abort(__LOCATION__, &
                       "MGGA functionals that require the Laplacian are "// &
                       "only compatible with 'XC_DERIV PW' and 'XC_SMOOTH_RHO NONE'")
      END IF

      IF (needs_rho_g) THEN
         CALL pw_pool_create_pw(pw_pool, tmp_g, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END IF
      DO ispin = 1, nspins
         CALL pw_pool_create_pw(pw_pool, my_rho_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         ! introduce a smoothing kernel on the density
         IF (xc_rho_smooth_id == xc_rho_no_smooth) THEN
            IF (needs_rho_g) THEN
               IF (ASSOCIATED(rho_g)) THEN
                  my_rho_g_local = .FALSE.
                  my_rho_g => rho_g(ispin)%pw
               END IF
            END IF

            CALL pw_copy(rho_r(ispin)%pw, my_rho_r(ispin)%pw)
         ELSE
            CALL xc_pw_smooth(rho_r(ispin)%pw, my_rho_r(ispin)%pw, xc_rho_smooth_id)
         END IF

         IF (gradient_f) THEN ! calculate the grad of rho
            ! normally when you need the gradient you need the whole gradient
            ! (for the partial integration)
            ! deriv rho
            DO idir = 1, 3
               NULLIFY (drho_r(idir, ispin)%pw)
               CALL pw_pool_create_pw(pw_pool, drho_r(idir, ispin)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
            END DO
            IF (needs_rho_g) THEN
               IF (.NOT. ASSOCIATED(my_rho_g)) THEN
                  my_rho_g_local = .TRUE.
                  CALL pw_pool_create_pw(pw_pool, my_rho_g, &
                                         use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  CALL pw_transfer(my_rho_r(ispin)%pw, my_rho_g)
               END IF
               IF (.NOT. my_rho_g_local .AND. (xc_deriv_method_id == xc_deriv_spline2 .OR. &
                                               xc_deriv_method_id == xc_deriv_spline3)) THEN
                  CALL pw_pool_create_pw(pw_pool, my_rho_g, &
                                         use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  my_rho_g_local = .TRUE.
                  CALL pw_copy(rho_g(ispin)%pw, my_rho_g)
               END IF
            END IF
            IF (needs%laplace_rho .OR. needs%laplace_rho_spin) THEN
               NULLIFY (laplace_rho_r(ispin)%pw)
               CALL pw_pool_create_pw(pw_pool, laplace_rho_r(ispin)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL xc_pw_laplace(my_rho_g, pw_pool, xc_deriv_method_id, laplace_rho_r(ispin)%pw, tmp_g=tmp_g)
            END IF
            DO idir = 1, 3
               CALL pw_zero(drho_r(idir, ispin)%pw)
               CALL xc_pw_derive(my_rho_r(ispin)%pw, tmp_g, drho_r(idir, ispin)%pw, idir, xc_deriv_method_id, pw_g=my_rho_g)
            END DO

            IF (needs_rho_g) THEN
               IF (my_rho_g_local) THEN
                  my_rho_g_local = .FALSE.
                  CALL pw_pool_give_back_pw(pw_pool, my_rho_g)
               END IF
            END IF

            IF (xc_deriv_method_id /= xc_deriv_pw) THEN
               CALL pw_spline_scale_deriv(drho_r(:, ispin))
            END IF

         END IF

      END DO

      IF (ASSOCIATED(tmp_g)) CALL pw_pool_give_back_pw(pw_pool, tmp_g)

      SELECT CASE (nspins)
      CASE (1)
         IF (needs%rho_1_3) THEN
            CALL pw_pool_create_cr3d(pw_pool, rho_set%rho_1_3)
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rho_set,my_rho_r)
            DO k = rho_set%local_bounds(1, 3), rho_set%local_bounds(2, 3)
               DO j = rho_set%local_bounds(1, 2), rho_set%local_bounds(2, 2)
                  DO i = rho_set%local_bounds(1, 1), rho_set%local_bounds(2, 1)
                     rho_set%rho_1_3(i, j, k) = MAX(my_rho_r(1)%pw%cr3d(i, j, k), 0.0_dp)**f13
                  END DO
               END DO
            END DO
            rho_set%owns%rho_1_3 = .TRUE.
            rho_set%has%rho_1_3 = .TRUE.
         END IF
         IF (needs%rho) THEN
            rho_set%rho => my_rho_r(1)%pw%cr3d
            NULLIFY (my_rho_r(1)%pw%cr3d)
            rho_set%owns%rho = .TRUE.
            rho_set%has%rho = .TRUE.
         END IF
         IF (needs%norm_drho) THEN
            CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drho)
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rho_set,drho_r)
            DO k = rho_set%local_bounds(1, 3), rho_set%local_bounds(2, 3)
               DO j = rho_set%local_bounds(1, 2), rho_set%local_bounds(2, 2)
                  DO i = rho_set%local_bounds(1, 1), rho_set%local_bounds(2, 1)
                     rho_set%norm_drho(i, j, k) = SQRT( &
                                                  drho_r(1, 1)%pw%cr3d(i, j, k)**2 + &
                                                  drho_r(2, 1)%pw%cr3d(i, j, k)**2 + &
                                                  drho_r(3, 1)%pw%cr3d(i, j, k)**2)
                  END DO
               END DO
            END DO
            rho_set%owns%norm_drho = .TRUE.
            rho_set%has%norm_drho = .TRUE.
         END IF
         IF (needs%laplace_rho) THEN
            rho_set%laplace_rho => laplace_rho_r(1)%pw%cr3d
            NULLIFY (laplace_rho_r(1)%pw%cr3d)
            rho_set%owns%laplace_rho = .TRUE.
            rho_set%has%laplace_rho = .TRUE.
         END IF

         IF (needs%drho) THEN
            DO idir = 1, 3
               rho_set%drho(idir)%array => drho_r(idir, 1)%pw%cr3d
               NULLIFY (drho_r(idir, 1)%pw%cr3d)
            END DO
            rho_set%owns%drho = .TRUE.
            rho_set%has%drho = .TRUE.
         END IF
      CASE (2)
         IF (needs%rho_spin_1_3) THEN
            CALL pw_pool_create_cr3d(pw_pool, rho_set%rhoa_1_3)
            !assume that the bounds are the same?
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rho_set,my_rho_r)
            DO k = rho_set%local_bounds(1, 3), rho_set%local_bounds(2, 3)
               DO j = rho_set%local_bounds(1, 2), rho_set%local_bounds(2, 2)
                  DO i = rho_set%local_bounds(1, 1), rho_set%local_bounds(2, 1)
                     rho_set%rhoa_1_3(i, j, k) = MAX(my_rho_r(1)%pw%cr3d(i, j, k), 0.0_dp)**f13
                  END DO
               END DO
            END DO
            CALL pw_pool_create_cr3d(pw_pool, rho_set%rhob_1_3)
            !assume that the bounds are the same?
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rho_set,my_rho_r)
            DO k = rho_set%local_bounds(1, 3), rho_set%local_bounds(2, 3)
               DO j = rho_set%local_bounds(1, 2), rho_set%local_bounds(2, 2)
                  DO i = rho_set%local_bounds(1, 1), rho_set%local_bounds(2, 1)
                     rho_set%rhob_1_3(i, j, k) = MAX(my_rho_r(2)%pw%cr3d(i, j, k), 0.0_dp)**f13
                  END DO
               END DO
            END DO
            rho_set%owns%rho_spin_1_3 = .TRUE.
            rho_set%has%rho_spin_1_3 = .TRUE.
         END IF
         IF (needs%norm_drho) THEN

            CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drho)
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rho_set,drho_r)
            DO k = rho_set%local_bounds(1, 3), rho_set%local_bounds(2, 3)
               DO j = rho_set%local_bounds(1, 2), rho_set%local_bounds(2, 2)
                  DO i = rho_set%local_bounds(1, 1), rho_set%local_bounds(2, 1)
                     rho_set%norm_drho(i, j, k) = SQRT( &
                                                  (drho_r(1, 1)%pw%cr3d(i, j, k) + drho_r(1, 2)%pw%cr3d(i, j, k))**2 + &
                                                  (drho_r(2, 1)%pw%cr3d(i, j, k) + drho_r(2, 2)%pw%cr3d(i, j, k))**2 + &
                                                  (drho_r(3, 1)%pw%cr3d(i, j, k) + drho_r(3, 2)%pw%cr3d(i, j, k))**2)
                  END DO
               END DO
            END DO

            rho_set%owns%norm_drho = .TRUE.
            rho_set%has%norm_drho = .TRUE.
         END IF
         IF (needs%rho_spin) THEN

            rho_set%rhoa => my_rho_r(1)%pw%cr3d
            NULLIFY (my_rho_r(1)%pw%cr3d)

            rho_set%rhob => my_rho_r(2)%pw%cr3d
            NULLIFY (my_rho_r(2)%pw%cr3d)

            rho_set%owns%rho_spin = .TRUE.
            rho_set%has%rho_spin = .TRUE.
         END IF
         IF (needs%norm_drho_spin) THEN

            CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drhoa)
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rho_set,drho_r)
            DO k = rho_set%local_bounds(1, 3), rho_set%local_bounds(2, 3)
               DO j = rho_set%local_bounds(1, 2), rho_set%local_bounds(2, 2)
                  DO i = rho_set%local_bounds(1, 1), rho_set%local_bounds(2, 1)
                     rho_set%norm_drhoa(i, j, k) = SQRT( &
                                                   drho_r(1, 1)%pw%cr3d(i, j, k)**2 + &
                                                   drho_r(2, 1)%pw%cr3d(i, j, k)**2 + &
                                                   drho_r(3, 1)%pw%cr3d(i, j, k)**2)
                  END DO
               END DO
            END DO

            CALL pw_pool_create_cr3d(pw_pool, rho_set%norm_drhob)
            rho_set%owns%norm_drho_spin = .TRUE.
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(rho_set,drho_r)
            DO k = rho_set%local_bounds(1, 3), rho_set%local_bounds(2, 3)
               DO j = rho_set%local_bounds(1, 2), rho_set%local_bounds(2, 2)
                  DO i = rho_set%local_bounds(1, 1), rho_set%local_bounds(2, 1)
                     rho_set%norm_drhob(i, j, k) = SQRT( &
                                                   drho_r(1, 2)%pw%cr3d(i, j, k)**2 + &
                                                   drho_r(2, 2)%pw%cr3d(i, j, k)**2 + &
                                                   drho_r(3, 2)%pw%cr3d(i, j, k)**2)
                  END DO
               END DO
            END DO

            rho_set%owns%norm_drho_spin = .TRUE.
            rho_set%has%norm_drho_spin = .TRUE.
         END IF
         IF (needs%laplace_rho_spin) THEN
            rho_set%laplace_rhoa => laplace_rho_r(1)%pw%cr3d
            NULLIFY (laplace_rho_r(1)%pw%cr3d)

            rho_set%laplace_rhob => laplace_rho_r(2)%pw%cr3d
            NULLIFY (laplace_rho_r(2)%pw%cr3d)

            rho_set%owns%laplace_rho_spin = .TRUE.
            rho_set%has%laplace_rho_spin = .TRUE.
         END IF
         IF (needs%drho_spin) THEN
            DO idir = 1, 3
               rho_set%drhoa(idir)%array => drho_r(idir, 1)%pw%cr3d
               NULLIFY (drho_r(idir, 1)%pw%cr3d)
               rho_set%drhob(idir)%array => drho_r(idir, 2)%pw%cr3d
               NULLIFY (drho_r(idir, 2)%pw%cr3d)
            END DO
            rho_set%owns%drho_spin = .TRUE.
            rho_set%has%drho_spin = .TRUE.
         END IF
      END SELECT
      ! post cleanup
      DO ispin = 1, nspins
         IF (needs%laplace_rho .OR. needs%laplace_rho_spin) THEN
            CALL pw_pool_give_back_pw(pw_pool, laplace_rho_r(ispin)%pw, &
                                      accept_non_compatible=.TRUE.)
         END IF
         DO idir = 1, 3
            CALL pw_pool_give_back_pw(pw_pool, drho_r(idir, ispin)%pw, &
                                      accept_non_compatible=.TRUE.)
         END DO
      END DO
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(pw_pool, my_rho_r(ispin)%pw, &
                                   accept_non_compatible=.TRUE.)
      END DO

      ! tau part
      IF (needs%tau .OR. needs%tau_spin) THEN
         DO ispin = 1, nspins
            CPASSERT(ASSOCIATED(tau(ispin)%pw))
            CPASSERT(ASSOCIATED(tau(ispin)%pw%cr3d))
         END DO
      END IF
      IF (needs%tau) THEN
         rho_set%tau => tau(1)%pw%cr3d
         rho_set%owns%tau = .FALSE.
         rho_set%has%tau = .TRUE.
      END IF
      IF (needs%tau_spin) THEN
         rho_set%tau_a => tau(1)%pw%cr3d
         rho_set%tau_b => tau(2)%pw%cr3d
         rho_set%owns%tau_spin = .FALSE.
         rho_set%has%tau_spin = .TRUE.
      END IF

      CPASSERT(xc_rho_cflags_equal(rho_set%has, needs))

   END SUBROUTINE xc_rho_set_update

END MODULE xc_rho_set_types
