!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Sets a data pointer.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
! *****************************************************************************
  SUBROUTINE set_data_p_z (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_complex(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%c_dp => p
  END SUBROUTINE set_data_p_z

! *****************************************************************************
!> \brief Sets a data pointer.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
! *****************************************************************************
  SUBROUTINE set_data_p_2d_z (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER         :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_2d_complex(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%c2_dp => p
  END SUBROUTINE set_data_p_2d_z


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval data          pointer to data
!> \par Calling
!>      This routine is hidden behind the dbcsr_get_data interface, hence the
!>      need for the coersion argument.
!>      see dbcsr_get_data_p_z
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_z (area, coersion, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_8), INTENT(IN)            :: coersion
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_c_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    TYPE(dbcsr_error_type)         :: error

!   ---------------------------------------------------------------------------

    ! The coersion argument is needed to make this function unique
    ! enough to use in the interface.
    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l .GE. LBOUND (area%d%c_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u .LE. UBOUND (area%d%c_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_z

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \brief dbcsr_get_data_c_z
!> \param[in] area       data area
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval DATA pointer to data
!> \par Calling
!>      This routine can be called explicitly.
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_z (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l .GE. LBOUND (area%d%c_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u .LE. UBOUND (area%d%c_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_z

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \brief dbcsr_get_data_c_z
!> \param[in] area       data area
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval DATA pointer to data
!> \par Calling
!>      This routine can be called explicitly.
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_2d_z (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_2d_z', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c2_dp)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l(1) .GE. LBOUND (area%d%c2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (l(2) .GE. LBOUND (area%d%c2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(1) .LE. UBOUND (area%d%c2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(2) .LE. UBOUND (area%d%c2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%c2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%c2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_2d_z



! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] DATA pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_z (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)  :: area
    COMPLEX(kind=real_8), DIMENSION(:), POINTER    :: DATA
    INTEGER, INTENT(IN), OPTIONAL     :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_z', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l, "GE", LBOUND (area%d%c_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u, "LE", UBOUND (area%d%c_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%c_dp(l:u)
       ELSE
          DATA => area%d%c_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_z


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] DATA pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_2d_z (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_2d_z', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%c2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%c2_dp)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l(1), "GE", LBOUND (area%d%c2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (l(2), "GE", LBOUND (area%d%c2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(1), "LE", UBOUND (area%d%c2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(2), "LE", UBOUND (area%d%c2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%c2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%c2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_2d_z


! *****************************************************************************
!> \brief Returns the entire data for a matrix.
!> \param matrix ...
!> \param[out] DATA pointer to data
!> \par Warning
!>      This routine should only be used within DBCSR code.
! *****************************************************************************
  SUBROUTINE get_data_m_z (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: DATA

    CALL get_data_z (matrix%m%data_area, DATA)
  END SUBROUTINE get_data_m_z



! *****************************************************************************
!> \brief Sets a scalar in an encapsulated data structure
!> \param[in] scalar                    scalar to encapsulate
!> \retval encapsulated_scalar          encapsulated scalar 
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_z (scalar) RESULT (encapsulated_scalar)
    COMPLEX(kind=real_8), INTENT(IN)       :: scalar
    TYPE(dbcsr_scalar_type)   :: encapsulated_scalar

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scalar_complex(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    encapsulated_scalar = dbcsr_scalar_zero (dbcsr_type_complex_8)
    encapsulated_scalar%c_dp = scalar
  END FUNCTION dbcsr_scalar_z

! *****************************************************************************
!> \brief Sets a scalar in an encapsulated data structure
!> \param encapsulated_scalar ...
!> \param value ...
!> \params[in] encapsulated_scalar          encapsulated scalar
!> \params[out] value                       value of the scalar
! *****************************************************************************
  ELEMENTAL SUBROUTINE dbcsr_scalar_get_value_z (encapsulated_scalar, value)
    TYPE(dbcsr_scalar_type), INTENT(IN) :: encapsulated_scalar
    COMPLEX(kind=real_8), INTENT(OUT)                :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scalar_get_value_complex(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    value = encapsulated_scalar%c_dp
  END SUBROUTINE dbcsr_scalar_get_value_z


! *****************************************************************************
!> \brief Used to determine appropriate type for data.
!> \param[in] data                      data to query
!> \retval data_type                    appropriate data_type 
! *****************************************************************************
  PURE FUNCTION query_type_z_1d (DATA) RESULT (data_type)
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN) :: DATA
    INTEGER                           :: data_type
    data_type = dbcsr_type_complex_8
  END FUNCTION query_type_z_1d
! *****************************************************************************
!> \brief ...
!> \param data ...
!> \retval data_type ...
! *****************************************************************************
  PURE FUNCTION query_type_z_2d (DATA) RESULT (data_type)
    COMPLEX(kind=real_8), DIMENSION(:,:), INTENT(IN) :: DATA
    INTEGER                             :: data_type
    data_type = dbcsr_type_complex_8_2d
  END FUNCTION query_type_z_2d
