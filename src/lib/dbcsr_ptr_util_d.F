! *****************************************************************************
!> \brief Returns a pointer with different bounds.
!> \param[in] original   original data pointer
!> \param[in] lb, ub     lower and upper bound for the new pointer view
!> \param[out] view      new pointer
! *****************************************************************************
  FUNCTION pointer_view_d (original, lb, ub) RESULT (view)
    REAL(kind=real_8), DIMENSION(:), POINTER :: original, view
    INTEGER, INTENT(IN)                  :: lb, ub
    view => original(lb:ub)
  END FUNCTION pointer_view_d


! *****************************************************************************
!> \brief Ensures that an array is appropriately large.
!> \param[in,out] array       array to verify and possibly resize
!> \param[in] lb    (optional) desired array lower bound
!> \param[in] ub    desired array upper bound
!> \param[in] factor          (optional) factor by which to exagerrate
!>                            enlargements
!> \param[in] nocopy          (optional) copy array on enlargement; default
!>                            is to copy
!> \param[in] special         (optional) use MPI-allocated memory; default
!>                            is no
!> \param[in] zero_pad        (optional) zero new allocations; default is to
!>                            write nothing
! *****************************************************************************
  SUBROUTINE ensure_array_size_d(array, lb, ub, factor,&
       nocopy, special, zero_pad, error)
    REAL(kind=real_8), DIMENSION(:), POINTER           :: array
    INTEGER, INTENT(IN), OPTIONAL            :: lb
    INTEGER, INTENT(IN)                      :: ub
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, special, zero_pad
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ensure_array_size_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: lb_new, lb_orig, stat, &
                                                ub_new, ub_orig, old_size,&
                                                size_increase, error_handler
    LOGICAL                                  :: dbg, docopy = .TRUE., &
                                                spec,&
                                                pad
    REAL(kind=real_8), DIMENSION(:), POINTER           :: newarray

!   ---------------------------------------------------------------------------
    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.

    docopy = .TRUE.
    IF (PRESENT (nocopy)) docopy = .NOT. nocopy
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    lb_new = 1
    IF (PRESENT (lb)) lb_new = lb
    pad = .FALSE.
    IF (PRESENT (zero_pad)) pad = zero_pad
    !> Creates a new array if it doesn't yet exist.
    IF (.NOT.ASSOCIATED(array)) THEN
       IF (spec) THEN
          CALL dbcsr_assert (lb.EQ.1, dbcsr_warning_level, dbcsr_internal_error,&
               routineN, "Special memory can only start at 1.",__LINE__,error)
          CALL mp_allocate(array, ub-lb_new+1, stat=stat)
       ELSE
          ALLOCATE(array(lb_new:ub), stat=stat)
       ENDIF
       CALL dbcsr_assert (stat == 0, dbcsr_warning_level, dbcsr_internal_error,&
            routineN, "array",__LINE__,error)
       IF (pad) array(:) = 0.0_real_8
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    lb_orig = LBOUND(array,1)
    ub_orig = UBOUND(array,1)
    old_size = ub_orig - lb_orig + 1
    ! The existing array is big enough.
    IF (lb_orig.LE.lb_new .AND. ub_orig.GE.ub) THEN
       CALL dbcsr_error_stop(error_handler, error)
       RETURN
    ENDIF
    ! A reallocation must be performed
    IF(dbg) WRITE(*,*)routineP//' Current bounds are',lb_orig,':',ub_orig,&
         '; special?',spec
    !CALL timeset(routineN,timing_handle)
    IF (lb_orig.GT.lb_new) THEN
       IF (PRESENT(factor)) THEN
          size_increase = lb_orig - lb_new
          size_increase = MAX (NINT(REAL(size_increase)*factor),&
                               NINT(REAL(old_size)*factor))
          lb_new = MIN (lb_orig, lb_new - size_increase)
       ELSE
          lb_new = lb_orig
       ENDIF
    ENDIF
    IF (ub_orig.LT.ub) THEN
       IF (PRESENT(factor)) THEN
          size_increase = ub - ub_orig
          size_increase = MAX (NINT(REAL(size_increase)*factor),&
                               NINT(REAL(old_size)*factor))
          ub_new = MAX (ub_orig, ub + size_increase)
       ELSE
          ub_new = ub
       ENDIF
    ELSE
       ub_new = ub
    ENDIF
    IF(dbg) WRITE(*,*)routineP//' Resizing to bounds',lb_new,':',ub_new
    IF(.NOT.docopy) THEN
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
    ENDIF
    IF (spec) THEN
       CALL dbcsr_assert (lb_new.EQ.1, dbcsr_warning_level, dbcsr_internal_error,&
            routineN, "Special memory can only start at 1.",__LINE__,error)
       CALL mp_allocate(newarray, ub_new-lb_new+1, stat=stat)
    ELSE
       ALLOCATE(newarray(lb_new:ub_new), stat=stat)
    ENDIF
    CALL dbcsr_assert (stat == 0, dbcsr_warning_level, dbcsr_internal_error,&
         routineN, "newarray",__LINE__,error)
    IF(docopy) THEN
       IF(dbg) CALL dbcsr_assert(lb_new.LE.lb_orig .AND. ub_new.GE.ub_orig,&
            dbcsr_failure_level, dbcsr_internal_error, routineP,&
            "Old extent exceeds the new one.",__LINE__,error)
       newarray(lb_orig:ub_orig) = array(lb_orig:ub_orig)
       IF (spec) THEN
          CALL mp_deallocate(array)
       ELSE
          DEALLOCATE(array)
       ENDIF
       IF (pad) THEN
          array(ub_orig+1:ub_new) = 0.0_real_8
          array(lb_new:lb_orig-1) = 0.0_real_8
       ENDIF
    ELSEIF (pad) THEN
       array(:) = 0.0_real_8
    ENDIF
    array => newarray
    IF (dbg) WRITE(*,*)routineP//' New array size', SIZE(array)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE ensure_array_size_d


! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to
!> \param[in] from  array to copy from
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_ad(to, from, offset, factor, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: to
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_ad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                ub_t, error_handler
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------
    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    spec = to%d%special
    s_f = SIZE(from, 1)
    IF (dbg) THEN
       WRITE(*,*)routineP//' from size',SIZE(from),';', SIZE(from,1)
       WRITE(*,*)routineP//' from lbound',LBOUND(from),';', LBOUND(from,1)
       WRITE(*,*)routineP//' from ubound',UBOUND(from),';', UBOUND(from,1)
    ENDIF
    IF (.NOT. ASSOCIATED (to%d%r_dp)) THEN
       IF (PRESENT (offset)) THEN
          start = offset
          last = offset + SIZE(from)-1
          CALL ensure_array_size(to%d%r_dp, lb=start,&
               ub=last,&
               factor=1.0, nocopy=.TRUE., special=spec, error=error)
       ELSE
          start = 1
          last = SIZE(from)
          CALL ensure_array_size(to%d%r_dp, ub=last,&
               factor=1.0, nocopy=.TRUE., special=spec, error=error)
       ENDIF
    ELSE
       ub_t = UBOUND(to%d%r_dp, 1)
       lb_t = LBOUND(to%d%r_dp, 1)
       s_t = SIZE(to%d%r_dp)
       start = lb_t
       IF (PRESENT (offset)) start = offset
       ! assert offset > lb_t
       last = start+s_f-1
       IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
            'starting at',start,'last',last,'source size', s_f
       IF (ub_t .LT. last) THEN
          IF (PRESENT (factor)) THEN
             f = factor
          ELSE
             f = 1.4142136
          ENDIF
          !ub_t = MAX(lb_t + INT(s_t * f), lb_t+s_f-1)
          ub_t = MAX(start + INT(s_t * f), last)
          IF(dbg) WRITE(*,*)routineP//' expanding array from :',&
               UBOUND(to%d%r_dp),&
               'to',ub_t
          !CALL reallocate(to, lb_t, ub_t)
          CALL ensure_array_size(to%d%r_dp, lb_t, ub_t,&
               factor=1.0, nocopy=.FALSE.,&
               special=spec, error=error)
       ENDIF
    ENDIF
    to%d%r_dp(start:last) = from(:)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE addto_array_ad

! *****************************************************************************
!> \brief Copies an array into another array, ensuring that the target array
!>        is large enough.
!> \par Makes use of FORTRAN's array bounds for array sections.
!> \param to        array to copy to
!> \param[in] from  array to copy from
!> \param[in] offset          starting position in the target array
!> \param[in] factor          factor for size increase
!> \param error     cp2k error
! *****************************************************************************
  SUBROUTINE addto_array_d(to, from, offset, factor, special, error)
    REAL(kind=real_8), DIMENSION(:), POINTER           :: to
    REAL(kind=real_8), DIMENSION(:), INTENT(IN)        :: from
    INTEGER, INTENT(IN), OPTIONAL            :: offset
    REAL, INTENT(IN), OPTIONAL               :: factor
    LOGICAL, INTENT(IN), OPTIONAL            :: special
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'addto_array_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: last, lb_t, s_f, s_t, start, &
                                                ub_t, error_handler
    LOGICAL                                  :: dbg, spec
    REAL                                     :: f

!   ---------------------------------------------------------------------------
    CALL dbcsr_error_set(routineN, error_handler, error)
    dbg = .FALSE.
    spec = .FALSE.
    IF (PRESENT (special)) spec = special
    s_f = SIZE(from, 1)
    IF (dbg) THEN
       WRITE(*,*)routineP//' from size',SIZE(from),';', SIZE(from,1)
       WRITE(*,*)routineP//' from lbound',LBOUND(from),';', LBOUND(from,1)
       WRITE(*,*)routineP//' from ubound',UBOUND(from),';', UBOUND(from,1)
    ENDIF
    ub_t = UBOUND(to, 1)
    lb_t = LBOUND(to, 1)
    s_t = SIZE(to)
    start = lb_t
    IF (PRESENT (offset)) start = offset
    ! assert offset > lb_t
    last = start+s_f-1
    IF(dbg) WRITE(*,*)routineP//' target array bounds',lb_t,':',ub_t,&
         'starting at',start,'last',last,'source size', s_f
    IF (ub_t .LT. last) THEN
       IF (PRESENT (factor)) THEN
          f = factor
       ELSE
          f = 1.4142136
       ENDIF
       !ub_t = MAX(lb_t + INT(s_t * f), lb_t+s_f-1)
       ub_t = MAX(start + INT(s_t * f), last)
       IF(dbg) WRITE(*,*)routineP//' expanding array from :',UBOUND(to),&
            'to',ub_t
       !CALL reallocate(to, lb_t, ub_t)
       CALL ensure_array_size(to, lb_t, ub_t, factor=1.0, nocopy=.FALSE.,&
            special=spec, error=error)
    ENDIF
    to(start:last) = from(:)
    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE addto_array_d



#if defined(__PTR_RANK_REMAP)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using Fortran 2003 pointer
!>        rank remapping.
! *****************************************************************************
  SUBROUTINE pointer_d_rank_remap2 (r2p, d1, d2, r1p)
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: r2p
    REAL(kind=real_8), DIMENSION(:), &
      POINTER                                :: r1p

    r2p(1:d1,1:d2) => r1p(1:d1*d2)
  END SUBROUTINE pointer_d_rank_remap2
#else
#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
! *****************************************************************************
  SUBROUTINE pointer_d_rank_remap2 (r2p, d1, d2, r1p)
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: r2p
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(kind=real_8), DIMENSION(d1, *), &
      TARGET                                 :: r1p

    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_d_rank_remap2
#else
! *****************************************************************************
!> \brief Not supported
! *****************************************************************************
  SUBROUTINE pointer_d_rank_remap2 (r2p, d1, d2, r1p)
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(kind=real_8), DIMENSION(:, :), &
      POINTER                                :: r2p
    REAL(kind=real_8), DIMENSION(d1*d2), &
      TARGET                                 :: r1p

!    r2p(1:d1,1:d2) => r1p(1:d1*d2)
    NULLIFY (r2p)
  END SUBROUTINE pointer_d_rank_remap2
#endif
#endif

