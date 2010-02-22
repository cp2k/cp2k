! *****************************************************************************
!> \brief Sets a data pointer.
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_s (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    REAL(kind=real_4), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_real(kind=real_4)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%r_sp => p
  END SUBROUTINE set_data_p_s

! *****************************************************************************
!> \brief Sets a data pointer.
!> \par Assumptions
!>      Assumes that no memory will be lost when repointing the
!>      pointer in the data area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_2d_s (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    REAL(kind=real_4), DIMENSION(:,:), POINTER         :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_2d_real(kind=real_4)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    area%d%r2_sp => p
  END SUBROUTINE set_data_p_2d_s


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine is hidden behind the dbcsr_get_data interface, hence the
!>      need for the coersion argument.
!> \sa dbcsr_get_data_p_s
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \retval data          pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_s (area, coersion, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(kind=real_4), INTENT(IN)            :: coersion
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub
    REAL(kind=real_4), DIMENSION(:), POINTER :: DATA

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_c_s', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    ! The coersion argument is needed to make this function unique
    ! enough to use in the interface. The following line is added to
    ! remove compiler warnings about unused arguments.
    l = coersion
    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r_sp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r_sp, 1)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL dbcsr_assert (l .GE. LBOUND (area%d%r_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u .LE. UBOUND (area%d%r_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r_sp(l:u)
       ELSE
          DATA => area%d%r_sp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_s

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine can be called explicitly.
!> \brief dbcsr_get_data_c_s
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_s (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(kind=real_4), DIMENSION(:), POINTER :: DATA
    INTEGER, INTENT(IN), OPTIONAL  :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_s', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r_sp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r_sp, 1)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL dbcsr_assert (l .GE. LBOUND (area%d%r_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u .LE. UBOUND (area%d%r_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r_sp(l:u)
       ELSE
          DATA => area%d%r_sp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_s

! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \par Calling
!>      This routine can be called explicitly.
!> \brief dbcsr_get_data_c_s
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_2d_s (area, lb, ub) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    REAL(kind=real_4), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_p_2d_s', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r2_sp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r2_sp)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL dbcsr_assert (l(1) .GE. LBOUND (area%d%r2_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (l(2) .GE. LBOUND (area%d%r2_sp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(1) .LE. UBOUND (area%d%r2_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(2) .LE. UBOUND (area%d%r2_sp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r2_sp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%r2_sp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_2d_s



! *****************************************************************************
!> \brief Returns the entire data for a matrix.
!> \par Warning
!>      This routine should only be used within DBCSR code.
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_s (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(kind=real_4), DIMENSION(:), POINTER :: DATA

    CALL get_data_s (matrix%m%data_area, DATA)
  END SUBROUTINE get_data_m_s


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_s (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)  :: area
    REAL(kind=real_4), DIMENSION(:), POINTER    :: DATA
    INTEGER, INTENT(IN), OPTIONAL     :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_s', &
      routineP = moduleN//':'//routineN

    INTEGER                        :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r_sp, 1)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r_sp, 1)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL dbcsr_assert (l .GE. LBOUND (area%d%r_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u .LE. UBOUND (area%d%r_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r_sp(l:u)
       ELSE
          DATA => area%d%r_sp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_s


! *****************************************************************************
!> \brief Returns the single/double precision real/complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
!> \param[in] lb         (optional) lower bound for pointer
!> \param[in] ub         (optional) upper bound for pointer
! *****************************************************************************
  SUBROUTINE get_data_2d_s (area, DATA, lb, ub)
    TYPE(dbcsr_data_obj), INTENT(IN)            :: area
    REAL(kind=real_4), DIMENSION(:,:), POINTER            :: DATA
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: lb, ub

    CHARACTER(len=*), PARAMETER :: routineN = 'get_data_2d_s', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(2)          :: l, u
    LOGICAL, PARAMETER             :: dbg = .TRUE.
    TYPE(dbcsr_error_type)         :: error
!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (area%d)) THEN
       IF (PRESENT (lb) .OR. PRESENT (ub)) THEN
          l = LBOUND (area%d%r2_sp)
          IF (PRESENT (lb)) l = lb
          u = UBOUND (area%d%r2_sp)
          IF (PRESENT (ub)) u = ub
          IF (dbg) THEN
             CALL dbcsr_assert (l(1) .GE. LBOUND (area%d%r2_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (l(2) .GE. LBOUND (area%d%r2_sp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(1) .LE. UBOUND (area%d%r2_sp, 1),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
             CALL dbcsr_assert (u(2) .LE. UBOUND (area%d%r2_sp, 2),&
                  dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
                  "Out of bounds",__LINE__,error)
          ENDIF
          DATA => area%d%r2_sp(l(1):u(1), l(2):u(2))
       ELSE
          DATA => area%d%r2_sp
       ENDIF
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_2d_s

! *****************************************************************************
!> \brief Sets a scalar in an encapsulated data structure
!> \param[in] scalar                    scalar to encapsulate
!> \result encapsulated_scalar          encapsulated scalar
! *****************************************************************************
  ELEMENTAL FUNCTION dbcsr_scalar_s (scalar) RESULT (encapsulated_scalar)
    REAL(kind=real_4), INTENT(IN)       :: scalar
    TYPE(dbcsr_scalar_type)   :: encapsulated_scalar

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_scalar_real(kind=real_4)', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    encapsulated_scalar%data_type = dbcsr_type_real_4
    encapsulated_scalar%r_sp = scalar
  END FUNCTION dbcsr_scalar_s


#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
#if defined(__PTR_RANK_REMAP)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using Fortran 2003 pointer
!>        rank remapping.
!>
!> Merge this with dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_s_rank_remap2 (r2p, d1, d2, r1p)
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: r2p
    REAL(kind=real_4), DIMENSION(d1*d2), &
      TARGET                                 :: r1p

    r2p(1:d1,1:d2) => r1p(1:d1*d2)
  END SUBROUTINE pointer_s_rank_remap2
#else
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
!>
!> Merge this with the same routine in dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_s_rank_remap2 (r2p, d1, d2, r1p)
    REAL(kind=real_4), DIMENSION(:, :), &
      POINTER                                :: r2p
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(kind=real_4), DIMENSION(d1, *), &
      TARGET                                 :: r1p

    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_s_rank_remap2
#endif
#endif

