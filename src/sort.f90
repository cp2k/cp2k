    IF (n==0) RETURN
    DO i = 1, n
       INDEX(i) = i
    END DO
    !
    IF (ALL(arr==arr(1))) RETURN ! Nothing to order..
    jstack = 0
    l = 1
    ir = n
1   IF (ir-l<m) THEN
       DO j = l + 1, ir
          a = arr(j)
          ib = INDEX(j)
          DO i = j - 1, 1, -1
             IF (arr(i)<=a) GO TO 2
             arr(i+1) = arr(i)
             INDEX(i+1) = INDEX(i)
          END DO
          i = 0
2         arr(i+1) = a
          INDEX(i+1) = ib
       END DO
       IF (jstack==0) RETURN
       ir = istack(jstack)
       l = istack(jstack-1)
       jstack = jstack - 2
    ELSE
       k = (l+ir)/2
       temp = arr(k)
       arr(k) = arr(l+1)
       arr(l+1) = temp
       itemp = INDEX(k)
       INDEX(k) = INDEX(l+1)
       INDEX(l+1) = itemp
       IF (arr(l+1)>arr(ir)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l)>arr(ir)) THEN
          temp = arr(l)
          arr(l) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l)
          INDEX(l) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l+1)>arr(l)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(l)
          arr(l) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(l)
          INDEX(l) = itemp
       END IF
       i = l + 1
       j = ir
       a = arr(l)
       ib = INDEX(l)
3      CONTINUE
       i = i + 1
       IF (arr(i)<a) GO TO 3
4      CONTINUE
       j = j - 1
       IF (arr(j)>a) GO TO 4
       IF (j<i) GO TO 5
       temp = arr(i)
       arr(i) = arr(j)
       arr(j) = temp
       itemp = INDEX(i)
       INDEX(i) = INDEX(j)
       INDEX(j) = itemp
       GO TO 3
5      arr(l) = arr(j)
       arr(j) = a
       INDEX(l) = INDEX(j)
       INDEX(j) = ib
       jstack = jstack + 2
       IF (jstack>nstack) STOP ' Nstack too small in sortr'
       IF (ir-i+1>=j-l) THEN
          istack(jstack) = ir
          istack(jstack-1) = i
          ir = j - 1
       ELSE
          istack(jstack) = j - 1
          istack(jstack-1) = l
          l = i
       END IF
    END IF
    GO TO 1
