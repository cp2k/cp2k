!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief superstucture that hold various representations of the density and
!>      keeps track of which ones are valid
!> \par History
!>      08.2002 created [fawzi]
!>      08.2014 kpoints [JGH]
!>      11.2014 make qs_rho_type PRIVATE [Ole Schuett]
!>      11.2014 unified k-point and gamma-point code [Ole Schuett]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qs_rho_types
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: get_1d_pointer,&
                                              get_2d_pointer,&
                                              kpoint_transitional_release,&
                                              kpoint_transitional_type,&
                                              set_1d_pointer,&
                                              set_2d_pointer
   USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_release,&
                                              pw_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_types'

   PUBLIC :: qs_rho_p_type, qs_rho_type
   PUBLIC :: qs_rho_release, qs_rho_clear_pwpool, &
             qs_rho_get, qs_rho_set, qs_rho_clear, qs_rho_create, qs_rho_unset_rho_ao

! **************************************************************************************************
!> \brief keeps the density in various representations, keeping track of
!>      which ones are valid.
!> \param most attributes are array with either lda or lsd_alpha,lsd_beta.
!> \param rho_ao the filtered rho in the localized atom basis (to have rho(r)
!>        the filtered matrix is enough, but rho(r,r') is lost).
!> \param rho_ao_kp the filtered rho in the localized atom basis (to have rho(r)
!>        the filtered matrix is enough, but rho(r,r') is lost).
!>        for kpoints, in real space index form
!> \param rho_r grids with rho in the real space
!> \param tau_r grids with the kinetic energy density in real space
!> \param rho_g grids with rho in the g space
!> \param tau_g grids with the kinetic energy density in g space
!> \param rho_g_valid , rho_r_valid, tau_r_valid, tau_g_valid: if the
!>        corresponding component is valid
!> \param tot_rho_r the total charge in r space (valid only if rho_r is)
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE qs_rho_type
      PRIVATE
      TYPE(kpoint_transitional_type)                 :: rho_ao
      TYPE(kpoint_transitional_type)                 :: rho_ao_im
      TYPE(pw_type), DIMENSION(:), POINTER           :: rho_r => Null(), &
                                                        rho_g => Null(), &
                                                        tau_r => Null(), &
                                                        tau_g => Null()
      TYPE(pw_type), DIMENSION(:, :), POINTER :: drho_r => NULL(), drho_g => NULL()
      ! Final rho_iter of last SCCS cycle (r-space)
      TYPE(pw_type), POINTER                         :: rho_r_sccs => Null()
      LOGICAL                                        :: rho_g_valid = .FALSE., &
                                                        rho_r_valid = .FALSE., &
                                                        drho_r_valid = .FALSE., &
                                                        drho_g_valid = .FALSE., &
                                                        tau_r_valid = .FALSE., &
                                                        tau_g_valid = .FALSE., &
                                                        soft_valid = .FALSE., &
                                                        complex_rho_ao = .FALSE.
      REAL(KIND=dp), DIMENSION(:), POINTER           :: tot_rho_r => Null(), &
                                                        tot_rho_g => Null()
   END TYPE qs_rho_type

! **************************************************************************************************
   TYPE qs_rho_p_type
      TYPE(qs_rho_type), POINTER                     :: rho
   END TYPE qs_rho_p_type

CONTAINS

! **************************************************************************************************
!> \brief Allocates a new instance of rho.
!> \param rho ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE qs_rho_create(rho)
      TYPE(qs_rho_type), INTENT(OUT)                     :: rho

   END SUBROUTINE qs_rho_create

! **************************************************************************************************
!> \brief releases a rho_struct by decreasing the reference count by one
!>      and deallocating if it reaches 0 (to be called when you don't want
!>      anymore a shared copy)
!> \param rho_struct the structure to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_rho_release(rho_struct)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      CALL qs_rho_clear(rho_struct)

   END SUBROUTINE qs_rho_release

! **************************************************************************************************
!> \brief Deallocates all components, without deallocating rho_struct itself.
!> \param rho_struct ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE qs_rho_clear(rho_struct)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      INTEGER                                            :: i, j

      IF (ASSOCIATED(rho_struct%rho_r)) THEN
         DO i = 1, SIZE(rho_struct%rho_r)
            CALL pw_release(rho_struct%rho_r(i))
         END DO
         DEALLOCATE (rho_struct%rho_r)
      END IF
      IF (ASSOCIATED(rho_struct%drho_r)) THEN
         DO j = 1, SIZE(rho_struct%drho_r, 2)
            DO i = 1, SIZE(rho_struct%drho_r, 1)
               CALL pw_release(rho_struct%drho_r(i, j))
            END DO
         END DO
         DEALLOCATE (rho_struct%drho_r)
      END IF
      IF (ASSOCIATED(rho_struct%drho_g)) THEN
         DO i = 1, SIZE(rho_struct%drho_g, 2)
            DO j = 1, SIZE(rho_struct%drho_g, 1)
               CALL pw_release(rho_struct%drho_g(i, j))
            END DO
         END DO
         DEALLOCATE (rho_struct%drho_g)
      END IF
      IF (ASSOCIATED(rho_struct%tau_r)) THEN
         DO i = 1, SIZE(rho_struct%tau_r)
            CALL pw_release(rho_struct%tau_r(i))
         END DO
         DEALLOCATE (rho_struct%tau_r)
      END IF
      IF (ASSOCIATED(rho_struct%rho_g)) THEN
         DO i = 1, SIZE(rho_struct%rho_g)
            CALL pw_release(rho_struct%rho_g(i))
         END DO
         DEALLOCATE (rho_struct%rho_g)
      END IF
      IF (ASSOCIATED(rho_struct%tau_g)) THEN
         DO i = 1, SIZE(rho_struct%tau_g)
            CALL pw_release(rho_struct%tau_g(i))
         END DO
         DEALLOCATE (rho_struct%tau_g)
      END IF
      IF (ASSOCIATED(rho_struct%rho_r_sccs)) THEN
         CALL pw_release(rho_struct%rho_r_sccs)
         DEALLOCATE (rho_struct%rho_r_sccs)
      END IF

      CALL kpoint_transitional_release(rho_struct%rho_ao)

      CALL kpoint_transitional_release(rho_struct%rho_ao_im)

      IF (ASSOCIATED(rho_struct%tot_rho_r)) &
         DEALLOCATE (rho_struct%tot_rho_r)
      IF (ASSOCIATED(rho_struct%tot_rho_g)) &
         DEALLOCATE (rho_struct%tot_rho_g)

   END SUBROUTINE qs_rho_clear

! **************************************************************************************************
!> \brief Unsets the rho_ao / rho_ao_kp field without calling kpoint_transitional_release().
!> \param rho_struct ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE qs_rho_unset_rho_ao(rho_struct)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      rho_struct%rho_ao = kpoint_transitional_type()
   END SUBROUTINE qs_rho_unset_rho_ao

! **************************************************************************************************
!> \brief returns info about the density described by this object.
!>      If some representation is not available an error is issued
!> \param rho_struct ...
!> \param rho_ao ...
!> \param rho_ao_im ...
!> \param rho_ao_kp ...
!> \param rho_ao_im_kp ...
!> \param rho_r ...
!> \param drho_r ...
!> \param rho_g ...
!> \param drho_g ...
!> \param tau_r ...
!> \param tau_g ...
!> \param rho_r_valid ...
!> \param drho_r_valid ...
!> \param rho_g_valid ...
!> \param drho_g_valid ...
!> \param tau_r_valid ...
!> \param tau_g_valid ...
!> \param tot_rho_r ...
!> \param tot_rho_g ...
!> \param rho_r_sccs ...
!> \param soft_valid ...
!> \param complex_rho_ao ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_rho_get(rho_struct, rho_ao, rho_ao_im, rho_ao_kp, rho_ao_im_kp, rho_r, drho_r, &
                         rho_g, drho_g, tau_r, tau_g, rho_r_valid, drho_r_valid, rho_g_valid, &
                         drho_g_valid, tau_r_valid, tau_g_valid, tot_rho_r, tot_rho_g, &
                         rho_r_sccs, soft_valid, complex_rho_ao)
      TYPE(qs_rho_type), INTENT(IN)                      :: rho_struct
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho_ao, rho_ao_im
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: rho_ao_kp, rho_ao_im_kp
      TYPE(pw_type), DIMENSION(:), OPTIONAL, POINTER     :: rho_r
      TYPE(pw_type), DIMENSION(:, :), OPTIONAL, POINTER  :: drho_r
      TYPE(pw_type), DIMENSION(:), OPTIONAL, POINTER     :: rho_g
      TYPE(pw_type), DIMENSION(:, :), OPTIONAL, POINTER  :: drho_g
      TYPE(pw_type), DIMENSION(:), OPTIONAL, POINTER     :: tau_r, tau_g
      LOGICAL, INTENT(out), OPTIONAL                     :: rho_r_valid, drho_r_valid, rho_g_valid, &
                                                            drho_g_valid, tau_r_valid, tau_g_valid
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: tot_rho_r, tot_rho_g
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho_r_sccs
      LOGICAL, INTENT(out), OPTIONAL                     :: soft_valid, complex_rho_ao

      IF (PRESENT(rho_ao)) rho_ao => get_1d_pointer(rho_struct%rho_ao)
      IF (PRESENT(rho_ao_kp)) rho_ao_kp => get_2d_pointer(rho_struct%rho_ao)

      IF (PRESENT(rho_ao_im)) rho_ao_im => get_1d_pointer(rho_struct%rho_ao_im)
      IF (PRESENT(rho_ao_im_kp)) rho_ao_im_kp => get_2d_pointer(rho_struct%rho_ao_im)

      IF (PRESENT(rho_r)) rho_r => rho_struct%rho_r
      IF (PRESENT(drho_r)) drho_r => rho_struct%drho_r
      IF (PRESENT(rho_g)) rho_g => rho_struct%rho_g
      IF (PRESENT(drho_g)) drho_g => rho_struct%drho_g
      IF (PRESENT(tau_r)) tau_r => rho_struct%tau_r
      IF (PRESENT(tau_g)) tau_g => rho_struct%tau_g
      IF (PRESENT(rho_r_valid)) rho_r_valid = rho_struct%rho_r_valid
      IF (PRESENT(rho_g_valid)) rho_g_valid = rho_struct%rho_g_valid
      IF (PRESENT(drho_r_valid)) drho_r_valid = rho_struct%drho_r_valid
      IF (PRESENT(drho_g_valid)) drho_g_valid = rho_struct%drho_g_valid
      IF (PRESENT(tau_r_valid)) tau_r_valid = rho_struct%tau_r_valid
      IF (PRESENT(tau_g_valid)) tau_g_valid = rho_struct%tau_g_valid
      IF (PRESENT(soft_valid)) soft_valid = rho_struct%soft_valid
      IF (PRESENT(tot_rho_r)) tot_rho_r => rho_struct%tot_rho_r
      IF (PRESENT(tot_rho_g)) tot_rho_g => rho_struct%tot_rho_g
      IF (PRESENT(rho_r_sccs)) rho_r_sccs => rho_struct%rho_r_sccs
      IF (PRESENT(complex_rho_ao)) complex_rho_ao = rho_struct%complex_rho_ao

   END SUBROUTINE qs_rho_get

! **************************************************************************************************
!> \brief ...
!> \param rho_struct ...
!> \param rho_ao ...
!> \param rho_ao_im ...
!> \param rho_ao_kp ...
!> \param rho_ao_im_kp ...
!> \param rho_r ...
!> \param drho_r ...
!> \param rho_g ...
!> \param drho_g ...
!> \param tau_r ...
!> \param tau_g ...
!> \param rho_r_valid ...
!> \param drho_r_valid ...
!> \param rho_g_valid ...
!> \param drho_g_valid ...
!> \param tau_r_valid ...
!> \param tau_g_valid ...
!> \param tot_rho_r ...
!> \param tot_rho_g ...
!> \param rho_r_sccs ...
!> \param soft_valid ...
!> \param complex_rho_ao ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE qs_rho_set(rho_struct, rho_ao, rho_ao_im, rho_ao_kp, rho_ao_im_kp, rho_r, drho_r, &
                         rho_g, drho_g, tau_r, tau_g, rho_r_valid, drho_r_valid, rho_g_valid, &
                         drho_g_valid, tau_r_valid, tau_g_valid, tot_rho_r, tot_rho_g, &
                         rho_r_sccs, soft_valid, complex_rho_ao)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho_ao, rho_ao_im
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: rho_ao_kp, rho_ao_im_kp
      TYPE(pw_type), DIMENSION(:), OPTIONAL, POINTER     :: rho_r
      TYPE(pw_type), DIMENSION(:, :), OPTIONAL, POINTER  :: drho_r
      TYPE(pw_type), DIMENSION(:), OPTIONAL, POINTER     :: rho_g
      TYPE(pw_type), DIMENSION(:, :), OPTIONAL, POINTER  :: drho_g
      TYPE(pw_type), DIMENSION(:), OPTIONAL, POINTER     :: tau_r, tau_g
      LOGICAL, INTENT(in), OPTIONAL                      :: rho_r_valid, drho_r_valid, rho_g_valid, &
                                                            drho_g_valid, tau_r_valid, tau_g_valid
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: tot_rho_r, tot_rho_g
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho_r_sccs
      LOGICAL, INTENT(in), OPTIONAL                      :: soft_valid, complex_rho_ao

      IF (PRESENT(rho_ao)) CALL set_1d_pointer(rho_struct%rho_ao, rho_ao)
      IF (PRESENT(rho_ao_kp)) CALL set_2d_pointer(rho_struct%rho_ao, rho_ao_kp)

      IF (PRESENT(rho_ao_im)) CALL set_1d_pointer(rho_struct%rho_ao_im, rho_ao_im)
      IF (PRESENT(rho_ao_im_kp)) CALL set_2d_pointer(rho_struct%rho_ao_im, rho_ao_im_kp)

      IF (PRESENT(rho_r)) rho_struct%rho_r => rho_r
      IF (PRESENT(rho_g)) rho_struct%rho_g => rho_g
      IF (PRESENT(drho_r)) rho_struct%drho_r => drho_r
      IF (PRESENT(drho_g)) rho_struct%drho_g => drho_g
      IF (PRESENT(tau_r)) rho_struct%tau_r => tau_r
      IF (PRESENT(tau_g)) rho_struct%tau_g => tau_g
      IF (PRESENT(rho_r_valid)) rho_struct%rho_r_valid = rho_r_valid
      IF (PRESENT(rho_g_valid)) rho_struct%rho_g_valid = rho_g_valid
      IF (PRESENT(drho_r_valid)) rho_struct%drho_r_valid = drho_r_valid
      IF (PRESENT(drho_g_valid)) rho_struct%drho_g_valid = drho_g_valid
      IF (PRESENT(tau_r_valid)) rho_struct%tau_r_valid = tau_r_valid
      IF (PRESENT(tau_g_valid)) rho_struct%tau_g_valid = tau_g_valid
      IF (PRESENT(soft_valid)) rho_struct%soft_valid = soft_valid
      IF (PRESENT(tot_rho_r)) rho_struct%tot_rho_r => tot_rho_r
      IF (PRESENT(tot_rho_g)) rho_struct%tot_rho_g => tot_rho_g
      IF (PRESENT(rho_r_sccs)) rho_struct%rho_r_sccs => rho_r_sccs
      IF (PRESENT(complex_rho_ao)) rho_struct%complex_rho_ao = complex_rho_ao

   END SUBROUTINE qs_rho_set
! **************************************************************************************************
!> \brief ...
!> \param rho_struct ...
!> \param auxbas_pw_pool ...
! **************************************************************************************************
   SUBROUTINE qs_rho_clear_pwpool(rho_struct, auxbas_pw_pool)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: auxbas_pw_pool

      INTEGER                                            :: i

      IF (ASSOCIATED(rho_struct%rho_r)) THEN
         DO i = 1, SIZE(rho_struct%rho_r)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_struct%rho_r(i))
         END DO
         DEALLOCATE (rho_struct%rho_r)
         NULLIFY (rho_struct%rho_r)
      END IF
      IF (ASSOCIATED(rho_struct%rho_g)) THEN
         DO i = 1, SIZE(rho_struct%rho_g)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_struct%rho_g(i))
         END DO
         DEALLOCATE (rho_struct%rho_g)
         NULLIFY (rho_struct%rho_g)
      END IF

   END SUBROUTINE qs_rho_clear_pwpool

END MODULE qs_rho_types
