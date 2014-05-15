!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Sets a data pointer.
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_d (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    REAL(kind=real_8), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_real(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%r_dp => p
  END SUBROUTINE set_data_p_d

! *****************************************************************************
!> \brief Sets a data pointer.
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_2d_d (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    REAL(kind=real_8), DIMENSION(:,:), POINTER         :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_2d_real(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%r2_dp => p
  END SUBROUTINE set_data_p_2d_d


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine is hidden behind the dbcsr_get_data interface, hence the
!>      need for the coersion argument.
!> \sa dbcsr_get_data_p_d
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval data          pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_d (area, coersion, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(kind=real_8), INTENT(IN)            :: coersion
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub
    REAL(kind=real_8), DIMENSION(:), POINTER :: DATA

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_c_d', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    TYPE(dbcsr_error_type)         :: error

!   ---------------------------------------------------------------------------

    ! The coersion argument is needed to make this function unique
    ! enough to use in the interface.
    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l .GE. LBOUND (area%d%r_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u .LE. UBOUND (area%d%r_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r_dp(l:u)
       ELSE
          DATA => area%d%r_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_d

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine can be called explicitly.
!> \brief dbcsr_get_data_c_d
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_d (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(kind=real_8), DIMENSION(:), POINTER :: DATA
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_d', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l .GE. LBOUND (area%d%r_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u .LE. UBOUND (area%d%r_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r_dp(l:u)
       ELSE
          DATA => area%d%r_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_d

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine can be called explicitly.
!> \brief dbcsr_get_data_c_d
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_2d_d (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    REAL(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_2d_d', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r2_dp)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l(1) .GE. LBOUND (area%d%r2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (l(2) .GE. LBOUND (area%d%r2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(1) .LE. UBOUND (area%d%r2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(2) .LE. UBOUND (area%d%r2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%r2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_2d_d



! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_d (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)  :: area
    REAL(kind=real_8), DIMENSION(:), POINTER    :: DATA
    INTEGER, INTENT(IN), OPTIONAL     :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_d', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r_dp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r_dp, 1)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l, "GE", LBOUND (area%d%r_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u, "LE", UBOUND (area%d%r_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r_dp(l:u)
       ELSE
          DATA => area%d%r_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_d


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_2d_d (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    REAL(kind=real_8), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_2d_d', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r2_dp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r2_dp)
          IF (PRESENT (ub)) u = ub
          IF (debug_mod) THEN
             CALL dbcsr_assert (l(1), "GE", LBOUND (area%d%r2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (l(2), "GE", LBOUND (area%d%r2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(1), "LE", UBOUND (area%d%r2_dp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(2), "LE", UBOUND (area%d%r2_dp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r2_dp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%r2_dp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_2d_d


! *****************************************************************************
!> \brief Returns the entire data for a matrix.
!> \par Warning
!>      This routine should only be used within DBCSR code.
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_d (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(kind=real_8), DIMENSION(:), POINTER :: DATA

    CALL get_data_d (matrix%m%data_area, DATA)
  END SUBROUTINE get_data_m_d



! *****************************************************************************
!> \brief Sets a scalar in an encapsulated data structure
!> \param[in] scalar                    scalar to encapsulate
!> \result encapsulated_scalar          encapsulated scalar
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_d (scalar) RESULT (encapsulated_scalar)
    REAL(kind=real_8), INTENT(IN)       :: scalar
    TYPE(dbcsr_scalar_type)   :: encapsulated_scalar

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scalar_real(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    encapsulated_scalar = dbcsr_scalar_zero (dbcsr_type_real_8)
    encapsulated_scalar%r_dp = scalar
  END FUNCTION dbcsr_scalar_d

! *****************************************************************************
!> \brief Sets a scalar in an encapsulated data structure
!> \params[in] encapsulated_scalar          encapsulated scalar
!> \params[out] value                       value of the scalar
! *****************************************************************************
  ELEMENTAL SUBROUTINE dbcsr_scalar_get_value_d (encapsulated_scalar, value)
    TYPE(dbcsr_scalar_type), INTENT(IN) :: encapsulated_scalar
    REAL(kind=real_8), INTENT(OUT)                :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scalar_get_value_real(kind=real_8)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    value = encapsulated_scalar%r_dp
  END SUBROUTINE dbcsr_scalar_get_value_d


! *****************************************************************************
!> \brief Used to determine appropriate type for data.
!> \param[in] data                      data to query
!> \result data_type                    appropriate data_type
! *****************************************************************************
  PURE FUNCTION query_type_d_1d (DATA) RESULT (data_type)
    REAL(kind=real_8), DIMENSION(:), INTENT(IN) :: DATA
    INTEGER                           :: data_type
    data_type = dbcsr_type_real_8
  END FUNCTION query_type_d_1d
  PURE FUNCTION query_type_d_2d (DATA) RESULT (data_type)
    REAL(kind=real_8), DIMENSION(:,:), INTENT(IN) :: DATA
    INTEGER                             :: data_type
    data_type = dbcsr_type_real_8_2d
  END FUNCTION query_type_d_2d
