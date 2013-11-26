!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_history
  USE kinds,                           ONLY: dp

 IMPLICIT NONE
 PRIVATE

 
 TYPE history_fingerprint_type
    PRIVATE
    REAL(KIND=dp)                            :: Epot = 0.0
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: goedecker
  END TYPE history_fingerprint_type


  TYPE history_entry_type
     TYPE(history_fingerprint_type), POINTER :: p => Null()
     INTEGER                                 :: id = -1
  END TYPE history_entry_type


  TYPE history_type
    PRIVATE
    TYPE(history_entry_type),DIMENSION(:), POINTER :: entries => Null()
    INTEGER                              :: length = 0
    INTEGER                              :: iw = -1
    !.. TODO make config params
    REAL(KIND=dp)                        :: E_precision =  1.0e-5
    REAL(KIND=dp)                        :: FP_precision =  1.0e-2
  END TYPE history_type



  PUBLIC :: history_type, history_fingerprint_type
  PUBLIC :: history_init, history_finalize
  PUBLIC :: history_add, history_lookup
  PUBLIC :: history_fingerprint
  PUBLIC :: history_fingerprint_match

  LOGICAL, PARAMETER                     :: debug = .FALSE.
  INTEGER, PARAMETER                     :: history_grow_unit = 1000
  CONTAINS


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE history_init(history, iw)
    TYPE(history_type), INTENT(INOUT)        :: history
    INTEGER                                  :: iw

    ALLOCATE(history%entries(history_grow_unit))
    history%iw = iw

 END SUBROUTINE history_init


! *****************************************************************************
! *****************************************************************************
 FUNCTION history_fingerprint(history, Epot, pos) RESULT(fp)
    TYPE(history_type), INTENT(IN)           :: history
    REAL(KIND=dp), INTENT(IN)                :: Epot
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pos
    TYPE(history_fingerprint_type)           :: fp

    fp%Epot = Epot
    fp%goedecker = goedecker_fingerprint(pos)

 END FUNCTION history_fingerprint

! *****************************************************************************
! *****************************************************************************
 FUNCTION history_fingerprint_match(history, fp1, fp2) RESULT(res)
    TYPE(history_type), INTENT(IN)           :: history
    TYPE(history_fingerprint_type), &
      INTENT(IN)                             :: fp1, fp2
    LOGICAL                                  :: res

    res = (ABS(fp1%Epot - fp2%Epot) < history%E_precision) .AND. &
          (fingerprint_distance(fp1, fp2) < history%fp_precision)

 END FUNCTION history_fingerprint_match

! *****************************************************************************
! *****************************************************************************
 SUBROUTINE history_finalize(history)
    TYPE(history_type)                       :: history

    DEALLOCATE(history%entries)

 END SUBROUTINE history_finalize


! *****************************************************************************
! *****************************************************************************
 FUNCTION interpolation_search(history, Efind) RESULT(res)
    TYPE(history_type), INTENT(IN)           :: history
    REAL(KIND=dp), INTENT(IN)                :: Efind
    INTEGER                                  :: res

    INTEGER                                  :: high, low, mid
    REAL(KIND=dp)                            :: slope

    low = 1
    high = history%length

    DO WHILE(low <= high)
       !linear interpolation
       slope = (high - low) / (history%entries(high)%p%Epot - history%entries(low)%p%Epot)
       mid = low + slope * (Efind - history%entries(low)%p%Epot) 
       mid = MIN(MAX(mid, low), high)

       IF (history%entries(mid)%p%Epot < Efind) THEN
        low = mid +1
       ELSE
        high = mid -1
       END IF
    END DO
    res = low
 END FUNCTION interpolation_search


! *****************************************************************************
! *****************************************************************************
 SUBROUTINE history_add(history, fingerprint, id)
    TYPE(history_type), INTENT(INOUT)        :: history
    TYPE(history_fingerprint_type), &
      INTENT(IN)                             :: fingerprint
    INTEGER, INTENT(IN), OPTIONAL            :: id

    INTEGER                                  :: k, n
    TYPE(history_entry_type), DIMENSION(:), &
      POINTER                                :: tmp

    n = SIZE(history%entries)
    IF(n == history%length) THEN
        ! grow history%entries array
        tmp => history%entries
        ALLOCATE(history%entries(n+history_grow_unit))
        history%entries(1:n) = tmp(:)
        DEALLOCATE(tmp)
        n = n + history_grow_unit
    ENDIF

    k = interpolation_search(history, fingerprint%Epot)

    history%entries(k+1:) = history%entries(k:n-1)
    ALLOCATE(history%entries(k)%p)
    history%entries(k)%p = fingerprint
    IF(PRESENT(id)) &
      history%entries(k)%id = id
    history%length = history%length + 1

    IF(debug) THEN
       ! check history for correct order
       DO k=1, history%length
          !WRITE(*,*) "history: ", k, "Epot",history%entries(k)%p%Epot
          IF(k>1)THEN
             IF(history%entries(k-1)%p%Epot > history%entries(k)%p%Epot) &
               STOP "history in wrong order"
           END IF
       END DO
    ENDIF

 END SUBROUTINE history_add


! *****************************************************************************
! *****************************************************************************
  SUBROUTINE history_lookup(history, fingerprint, found, id)
    TYPE(history_type), INTENT(IN)           :: history
    TYPE(history_fingerprint_type)           :: fingerprint
    LOGICAL, INTENT(OUT)                     :: found
    INTEGER, INTENT(OUT), OPTIONAL           :: id

    INTEGER                                  :: i, k, k_max, k_min
    REAL(KIND=dp)                            :: best_match, dist, Epot

    found = .FALSE.
    IF(PRESENT(id)) id = -1
    best_match = HUGE(1.0_dp)

    IF(history%length==0) RETURN

    Epot = fingerprint%Epot
    k = interpolation_search(history, fingerprint%Epot)

    DO k_min=k-1, 1, -1
      IF(history%entries(k_min)%p%Epot < Epot - history%E_precision) EXIT
    ENDDO

    DO k_max=k+1, history%length
      IF(history%entries(k_max)%p%Epot > Epot + history%E_precision) EXIT
    ENDDO

    k_min = MAX(k_min, 1)
    k_max = MIN(k_max, history%length)

    DO i=k_min, k_max
       dist = fingerprint_distance(fingerprint, history%entries(i)%p)
       !WRITE(*,*) "entry ", i, " dist: ",dist
       IF(dist < history%fp_precision .AND. dist < best_match) THEN
          best_match = dist
          found = .TRUE.
          IF(PRESENT(id)) id = history%entries(i)%id
       END IF
    ENDDO


 END SUBROUTINE history_lookup

! *****************************************************************************
! Credit: Stefan Goedecker
! ***************************************************************************** 
 FUNCTION goedecker_fingerprint(pos) RESULT(res)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: pos
    REAL(KIND=dp), DIMENSION(SIZE(pos)/3)    :: res

    INTEGER                                  :: i, info, j, N
    REAL(KIND=dp)                            :: d2, t
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: matrix, work
    REAL(KIND=dp), DIMENSION(3)              :: d

    N = SIZE(pos) / 3 ! number of atoms
    ALLOCATE(matrix(N,N), work(N,N))
       DO i=1,N
         matrix(i,i) = 1.0
         DO j=i+1,N
         d = pos(3*i-2:3*i) - pos(3*j-2:3*j)
         d2 = SUM(d**2)
         t=EXP(-0.5 * d2)
         matrix(i,j) = t
         matrix(j, i) = t
         ENDDO
       ENDDO
       !TODO: call dsyv through cp2k wrappers
       !TODO: do we need to store lower triangle of matrix?

       CALL DSYEV('N','U',N,matrix,N, res, work, N**2, info)
       IF (info.ne.0) STOP 'info'

 END FUNCTION goedecker_fingerprint


! *****************************************************************************
! Credit: Stefan Goedecker
! ***************************************************************************** 
 PURE FUNCTION fingerprint_distance(fp1, fp2) RESULT(res)
    TYPE(history_fingerprint_type), &
      INTENT(IN)                             :: fp1, fp2
    REAL(KIND=dp)                            :: res

     res = SQRT(SUM( (fp1%goedecker - fp2%goedecker)**2 )/SIZE(fp1%goedecker))
 END FUNCTION fingerprint_distance

END MODULE glbopt_history
