    isize   = iend-istart+1
    ! Initialize the INDEX array only for the first row..
    IF (j==1) THEN
       DO i = 1,isize
          INDEX(i) = i
       ENDDO
    END IF

    ! Allocate scratch arrays
    ALLOCATE(work(isize), work2(isize), tmp_index(isize), bck_index(isize))
    ind = 0
    DO i = istart, iend
      ind = ind + 1
      work(ind)      = matrix(j,i)
      bck_index(ind) = INDEX(i)
    END DO

    ! Ordering row (j) interval istart..iend
    CALL sort(work, isize, tmp_index)

    ! Copy into global INDEX array with a proper mapping
    ind = 0
    DO i = istart, iend
       ind = ind + 1
       INDEX(i)=bck_index(tmp_index(ind))
       matrix(j,i)=work(ind)
    END DO

    ! Reorder the rest of the array according the present reordering
    DO k = j+1, jsize
       ind = 0
       DO i = istart, iend
          ind = ind + 1
          work2(ind)      = matrix(k,i)
       END DO
       ind = 0
       DO i = istart, iend
          ind = ind + 1
          matrix(k,i)=work2(tmp_index(ind))
       END DO
    END DO

    ! There are more rows to order..
    IF (j<jsize) THEN
       kstart = istart
       item   = work(1)
       ind    = 0
       DO i = istart, iend
          ind = ind + 1
          IF (item/=work(ind)) THEN
             kend = i-1
             IF (kstart/=kend) THEN
                CALL sort(matrix, kstart, kend, j+1, jsize, INDEX)
             END IF
             item   = work(ind)
             kstart = i
          END IF
       END DO
       kend = i-1
       IF (kstart/=kend) THEN
          CALL sort(matrix, kstart, kend, j+1, jsize, INDEX)
       END IF
    END IF
    DEALLOCATE(work, work2, tmp_index, bck_index)
