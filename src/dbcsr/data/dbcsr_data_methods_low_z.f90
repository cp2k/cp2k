!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Sets a data pointer.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
! **************************************************************************************************
  SUBROUTINE set_data_p_z (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_z', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------
    IF(area%d%data_type /= dbcsr_type_complex_8) &
       CPABORT("set_data_p_z: data-area has wrong type")

    area%d%c_dp => p
  END SUBROUTINE set_data_p_z

! **************************************************************************************************
!> \brief Sets a data pointer.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
! **************************************************************************************************
  SUBROUTINE set_data_p_2d_z (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_2d_z', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF(area%d%data_type /= dbcsr_type_complex_8_2d) &
       CPABORT("set_data_p_2d_z: data-area has wrong type")

    area%d%c2_dp => p
  END SUBROUTINE set_data_p_2d_z


! **************************************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[in] select_data_type   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval data          pointer to data
!> \par Calling
!>      This routine is hidden behind the dbcsr_get_data interface, hence the
!>      need for the select_data_type argument.
!>      see dbcsr_get_data_p_z
! **************************************************************************************************
  FUNCTION dbcsr_get_data_c_z (area, select_data_type, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_8), INTENT(IN)            :: select_data_type
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_c_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u

!   ---------------------------------------------------------------------------

    ! The select_data_type argument is needed to make this function unique
    ! enough to use in the interface.
    IF (KIND(select_data_type) .NE. KIND(DATA))&
       CPABORT("compiler borken")

    IF (ASSOCIATED (area%d)) THEN
       IF(area%d%data_type /= dbcsr_type_complex_8) &
          CPABORT("dbcsr_get_data_c_z: data-area has wrong type")
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             IF(l .LT. LBOUND (area%d%c_dp, 1))&
                CPABORT("Out of bounds")
             IF(u .GT. UBOUND (area%d%c_dp, 1))&
                CPABORT("Out of bounds")
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_z

! **************************************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \brief dbcsr_get_data_c_z
!> \param[in] area       data area
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval DATA pointer to data
!> \par Calling
!>      This routine can be called explicitly.
! **************************************************************************************************
  FUNCTION dbcsr_get_data_p_z (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF(area%d%data_type /= dbcsr_type_complex_8) &
          CPABORT("dbcsr_get_data_p_z: data-area has wrong type")
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             IF(l .LT. LBOUND (area%d%c_dp, 1))&
                CPABORT("Out of bounds")
             IF(u .GT. UBOUND (area%d%c_dp, 1))&
                CPABORT("Out of bounds")
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_z

! **************************************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \brief dbcsr_get_data_c_z
!> \param[in] area       data area
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval DATA pointer to data
!> \par Calling
!>      This routine can be called explicitly.
! **************************************************************************************************
  FUNCTION dbcsr_get_data_p_2d_z (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_2d_z', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF(area%d%data_type /= dbcsr_type_complex_8_2d) &
          CPABORT("dbcsr_get_data_p_2d_z: data-area has wrong type")
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c2_dp)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             IF(l(1) .LT. LBOUND (area%d%c2_dp, 1))&
                CPABORT("Out of bounds")
             IF(l(2) .LT. LBOUND (area%d%c2_dp, 2))&
                CPABORT("Out of bounds")
             IF(u(1) .GT. UBOUND (area%d%c2_dp, 1))&
                CPABORT("Out of bounds")
             IF(u(2) .GT. UBOUND (area%d%c2_dp, 2))&
                CPABORT("Out of bounds")
          ENDIF
          DATA => area%d%c2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%c2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_2d_z



! **************************************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] DATA pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! **************************************************************************************************
  SUBROUTINE get_data_z (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)  :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER    :: DATA
    INTEGER, INTENT(IN), OPTIONAL     :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF(area%d%data_type /= dbcsr_type_complex_8) &
          CPABORT("get_data_z: data-area has wrong type")
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             IF(l < LBOUND (area%d%c_dp, 1)) &
                CPABORT("Out of bounds")
             IF(u > UBOUND (area%d%c_dp, 1)) &
                CPABORT("Out of bounds")
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_z


! **************************************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] DATA pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! **************************************************************************************************
  SUBROUTINE get_data_2d_z (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_2d_z', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF(area%d%data_type /= dbcsr_type_complex_8_2d) &
          CPABORT("get_data_2d_z: data-area has wrong type")
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c2_dp)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             IF(l(1) < LBOUND (area%d%c2_dp, 1)) &
                CPABORT("Out of bounds")
             IF(l(2) < LBOUND (area%d%c2_dp, 2)) &
                CPABORT("Out of bounds")
             IF(u(1) > UBOUND (area%d%c2_dp, 1)) &
                CPABORT("Out of bounds")
             IF(u(2) > UBOUND (area%d%c2_dp, 2)) &
                CPABORT("Out of bounds")
          ENDIF
          DATA => area%d%c2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%c2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_2d_z

! **************************************************************************************************
!> \brief Sets a scalar in an encapsulated data structure
!> \param[in] scalar                    scalar to encapsulate
!> \retval encapsulated_scalar          encapsulated scalar 
! **************************************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_z (scalar) RESULT (encapsulated_scalar)
    COMPLEX(kind=real_8), INTENT(IN)       :: scalar
    TYPE(dbcsr_scalar_type)   :: encapsulated_scalar

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scalar_z', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    encapsulated_scalar = dbcsr_scalar_zero (dbcsr_type_complex_8)
    encapsulated_scalar%c_dp = scalar
  END FUNCTION dbcsr_scalar_z

! **************************************************************************************************
!> \brief Sets a scalar in an encapsulated data structure
!> \param[in] encapsulated_scalar          encapsulated scalar
!> \param[out] value                       value of the scalar
! **************************************************************************************************
  ELEMENTAL SUBROUTINE dbcsr_scalar_get_value_z (encapsulated_scalar, value)
    TYPE(dbcsr_scalar_type), INTENT(IN) :: encapsulated_scalar
    COMPLEX(kind=real_8), INTENT(OUT)                :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scalar_get_value_z', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    value = encapsulated_scalar%c_dp
  END SUBROUTINE dbcsr_scalar_get_value_z
