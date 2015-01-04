!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Datatype to translate between k-points (2d) and gamma-point (1d) code.
!> \note  In principle storing just the 2d pointer would be sufficient.
!>        However due to a bug in ifort with the deallocation of
!>        bounds-remapped pointers, we also have to store the original
!>        1d pointer used for allocation.
!>
!> \par History
!>      11.2014 created [Ole Schuett]
!> \author Ole Schuett
! *****************************************************************************
MODULE kpoint_transitional
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

 PUBLIC :: kpoint_transitional_type, kpoint_transitional_release
 PUBLIC :: get_1d_pointer, get_2d_pointer, set_1d_pointer, set_2d_pointer

 TYPE kpoint_transitional_type
     PRIVATE
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER          :: ptr_1d => Null()
     TYPE(cp_dbcsr_p_type), DIMENSION(:,:), POINTER        :: ptr_2d => Null()
     LOGICAL                                               :: set_as_1d = .FALSE.
 END TYPE kpoint_transitional_type

CONTAINS

! *****************************************************************************
!> \brief Smart getter, raises an error when called during a k-point calculation
!> \param this ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION get_1d_pointer(this) RESULT(res)
    TYPE(kpoint_transitional_type)           :: this
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: res

     IF(ASSOCIATED(this%ptr_1d)) THEN
        IF(SIZE(this%ptr_2d, 2) /= 1) &
           CALL cp_unimplemented_error("kpoint_transitional","Method not implemented for k-points")
     ENDIF

     res => this%ptr_1d
  END FUNCTION get_1d_pointer

! *****************************************************************************
!> \brief Simple getter, neede because of PRIVATE
!> \param this ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
 FUNCTION get_2d_pointer(this) RESULT(res)
    TYPE(kpoint_transitional_type)           :: this
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: res

     res => this%ptr_2d
  END FUNCTION get_2d_pointer

! *****************************************************************************
!> \brief Assigns a 1D pointer
!> \param this ...
!> \param ptr_1d ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE set_1d_pointer(this, ptr_1d)
    TYPE(kpoint_transitional_type)           :: this
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ptr_1d

    INTEGER                                  :: n

     IF(ASSOCIATED(ptr_1d)) THEN
        n = SIZE(ptr_1d)
        this%ptr_1d           => ptr_1d
        this%ptr_2d(1:n, 1:1) => ptr_1d
        this%set_as_1d        = .TRUE.
     ELSE
         this%ptr_1d => Null()
         this%ptr_2d => Null()
     ENDIF
  END SUBROUTINE set_1d_pointer

! *****************************************************************************
!> \brief Assigns a 2D pointer
!> \param this ...
!> \param ptr_2d ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE set_2d_pointer(this, ptr_2d)
    TYPE(kpoint_transitional_type)           :: this
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: ptr_2d

     IF(ASSOCIATED(ptr_2d)) THEN
        this%ptr_1d    => ptr_2d(:,1)
        this%ptr_2d    => ptr_2d
        this%set_as_1d = .FALSE.
     ELSE
         this%ptr_1d => Null()
         this%ptr_2d => Null()
     ENDIF
  END SUBROUTINE set_2d_pointer

! *****************************************************************************
!> \brief Release the matrix set, using the right pointer
!> \param this ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE kpoint_transitional_release(this, error)
    TYPE(kpoint_transitional_type)           :: this
    TYPE(cp_error_type), INTENT(inout)       :: error

    IF (ASSOCIATED(this%ptr_1d)) THEN
      IF (this%set_as_1d) THEN
        CALL cp_dbcsr_deallocate_matrix_set(this%ptr_1d, error)
      ELSE
        CALL cp_dbcsr_deallocate_matrix_set(this%ptr_2d, error)
      ENDIF
    ENDIF
    NULLIFY(this%ptr_1d, this%ptr_2d)
  END SUBROUTINE kpoint_transitional_release

END MODULE kpoint_transitional
