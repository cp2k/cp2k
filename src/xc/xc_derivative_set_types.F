!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief represent a group ofunctional derivatives
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi & thomas
! **************************************************************************************************
MODULE xc_derivative_set_types
   USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_dealloc,&
                                              cp_sll_xc_deriv_insert_el,&
                                              cp_sll_xc_deriv_next,&
                                              cp_sll_xc_deriv_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_self
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_create,&
                                              pw_grid_release
   USE pw_methods,                      ONLY: pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create,&
                                              pw_pool_create_cr3d,&
                                              pw_pool_create_pw,&
                                              pw_pool_release,&
                                              pw_pool_retain,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_create,&
                                              pw_type
   USE xc_derivative_desc,              ONLY: standardize_desc
   USE xc_derivative_types,             ONLY: xc_derivative_create,&
                                              xc_derivative_release,&
                                              xc_derivative_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_derivative_set_types'

   PUBLIC :: xc_derivative_set_type
   PUBLIC :: xc_dset_create, xc_dset_release, &
             xc_dset_get_derivative, xc_dset_zero_all, xc_dset_recover_pw

! **************************************************************************************************
!> \brief A derivative set contains the different derivatives of a xc-functional
!>      in form of a linked list
! **************************************************************************************************
   TYPE xc_derivative_set_type
      TYPE(pw_pool_type), POINTER, PRIVATE :: pw_pool
      TYPE(cp_sll_xc_deriv_type), POINTER :: derivs
   END TYPE xc_derivative_set_type

CONTAINS

! **************************************************************************************************
!> \brief returns the requested xc_derivative
!> \param derivative_set the set where to search for the derivative
!> \param description the description of the derivative you want to have
!> \param allocate_deriv if the derivative should be allocated when not present
!>                        Defaults to false.
!> \return ...
! **************************************************************************************************
   FUNCTION xc_dset_get_derivative(derivative_set, description, allocate_deriv) &
      RESULT(res)

      TYPE(xc_derivative_set_type), INTENT(IN)           :: derivative_set
      INTEGER, DIMENSION(:), INTENT(in)                  :: description
      LOGICAL, INTENT(in), OPTIONAL                      :: allocate_deriv
      TYPE(xc_derivative_type), POINTER                  :: res

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: std_deriv_desc
      LOGICAL                                            :: my_allocate_deriv
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: cr3d_ptr
      TYPE(cp_sll_xc_deriv_type), POINTER                :: pos
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      NULLIFY (pos, deriv_att, cr3d_ptr)

      my_allocate_deriv = .FALSE.
      IF (PRESENT(allocate_deriv)) my_allocate_deriv = allocate_deriv
      NULLIFY (res)
      CALL standardize_desc(description, std_deriv_desc)
      pos => derivative_set%derivs
      DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
         IF (SIZE(deriv_att%split_desc) == SIZE(std_deriv_desc)) THEN
         IF (ALL(deriv_att%split_desc == std_deriv_desc)) THEN
            res => deriv_att
            EXIT
         END IF
         END IF
      END DO
      IF (.NOT. ASSOCIATED(res) .AND. my_allocate_deriv) THEN
         CALL pw_pool_create_cr3d(derivative_set%pw_pool, cr3d_ptr)
         cr3d_ptr = 0.0_dp
         ALLOCATE (res)
         CALL xc_derivative_create(res, std_deriv_desc, &
                                   cr3d_ptr=cr3d_ptr)
         CALL cp_sll_xc_deriv_insert_el(derivative_set%derivs, res)
      END IF
   END FUNCTION xc_dset_get_derivative

! **************************************************************************************************
!> \brief creates a derivative set object
!> \param derivative_set the set where to search for the derivative
!> \param pw_pool pool where to get the cr3d arrays needed to store the
!>        derivatives
!> \param local_bounds ...
! **************************************************************************************************
   SUBROUTINE xc_dset_create(derivative_set, pw_pool, local_bounds)

      TYPE(xc_derivative_set_type), INTENT(OUT)          :: derivative_set
      TYPE(pw_pool_type), OPTIONAL, POINTER              :: pw_pool
      INTEGER, DIMENSION(2, 3), INTENT(IN), OPTIONAL     :: local_bounds

      TYPE(pw_grid_type), POINTER                        :: pw_grid

      NULLIFY (pw_grid)

      NULLIFY (derivative_set%derivs)
      IF (PRESENT(pw_pool)) THEN
         derivative_set%pw_pool => pw_pool
         CALL pw_pool_retain(pw_pool)
         IF (PRESENT(local_bounds)) THEN
            IF (ANY(pw_pool%pw_grid%bounds_local /= local_bounds)) &
               CPABORT("incompatible local_bounds and pw_pool")
         END IF
      ELSE
         !FM ugly hack, should be replaced by a pool only for 3d arrays
         CPASSERT(PRESENT(local_bounds))
         CALL pw_grid_create(pw_grid, mp_comm_self)
         pw_grid%bounds_local = local_bounds
         NULLIFY (derivative_set%pw_pool)
         CALL pw_pool_create(derivative_set%pw_pool, pw_grid)
         CALL pw_grid_release(pw_grid)
      END IF

   END SUBROUTINE xc_dset_create

! **************************************************************************************************
!> \brief releases a derivative set
!> \param derivative_set the set to release
! **************************************************************************************************
   SUBROUTINE xc_dset_release(derivative_set)

      TYPE(xc_derivative_set_type)                       :: derivative_set

      TYPE(cp_sll_xc_deriv_type), POINTER                :: pos
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      NULLIFY (deriv_att, pos)

      pos => derivative_set%derivs
      DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
         CALL xc_derivative_release(deriv_att, pw_pool=derivative_set%pw_pool)
         DEALLOCATE (deriv_att)
      END DO
      CALL cp_sll_xc_deriv_dealloc(derivative_set%derivs)
      IF (ASSOCIATED(derivative_set%pw_pool)) CALL pw_pool_release(derivative_set%pw_pool)

   END SUBROUTINE xc_dset_release

! **************************************************************************************************
!> \brief ...
!> \param deriv_set ...
! **************************************************************************************************
   SUBROUTINE xc_dset_zero_all(deriv_set)

      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set

      TYPE(cp_sll_xc_deriv_type), POINTER                :: pos
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      NULLIFY (pos, deriv_att)

      IF (ASSOCIATED(deriv_set%derivs)) THEN
         pos => deriv_set%derivs
         DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
            deriv_att%deriv_data = 0.0_dp
         END DO
      END IF

   END SUBROUTINE xc_dset_zero_all

! **************************************************************************************************
!> \brief Recovers a derivative on a pw_type, the caller is responsible to release the grid later
!>        If the derivative is not found, either creates a blank pw_type from pw_pool or leaves it unassociated
!> \param deriv_set ...
!> \param description ...
!> \param pw ...
!> \param pw_grid ...
!> \param pw_pool create pw from this pool if derivative not found
! **************************************************************************************************
   SUBROUTINE xc_dset_recover_pw(deriv_set, description, pw, pw_grid, pw_pool)
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: description
      TYPE(pw_type), POINTER                             :: pw
      TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
      TYPE(pw_pool_type), INTENT(IN), OPTIONAL, POINTER  :: pw_pool

      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      deriv_att => xc_dset_get_derivative(deriv_set, description)
      IF (ASSOCIATED(deriv_att)) THEN
         CALL pw_create(pw, pw_grid=pw_grid, cr3d_ptr=deriv_att%deriv_data, &
                        use_data=REALDATA3D, in_space=REALSPACE)
         NULLIFY (deriv_att%deriv_data)
      ELSE IF (PRESENT(pw_pool)) THEN
         CALL pw_pool_create_pw(pw_pool, pw, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(pw)
      ELSE
         NULLIFY (pw)
      END IF

   END SUBROUTINE xc_dset_recover_pw

END MODULE xc_derivative_set_types
