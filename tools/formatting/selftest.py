
content = """
MODULE prettify_selftest
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: dp, test_routine, &
             test_function, test_type, str_function
   INTEGER, PARAMETER ::dp = SELECTED_REAL_KIND(15, 307)
   TYPE test_type
      REAL(KIND=dp) :: r = 1.0D-3
      INTEGER :: i
   END TYPE test_type

CONTAINS

   SUBROUTINE test_routine( &
      r, i, j, k, l)
      INTEGER, INTENT(IN)                                :: r, i, j, k
      INTEGER, INTENT(OUT)                               :: l

      l = test_function(r, i, j, k)
   END &
      SUBROUTINE

   PURE FUNCTION test_function(r, i, j, &
                               k) &
      RESULT(l)
      INTEGER, INTENT(IN)                                :: r, i, j, k
      INTEGER                                            :: l

      l = r+i+j+k
   END FUNCTION
   FUNCTION &
      str_function(a) RESULT(l)
      CHARACTER(LEN=*)                                   :: a
      INTEGER                                            :: l

      IF (LEN(a) < 5) THEN
         l = 0
      ELSE
         l = 1
      ENDIF
   END FUNCTION

END MODULE

PROGRAM example_prog
 USE example, ONLY: dp, test_routine, test_function, test_type,str_function

   IMPLICIT NONE
   INTEGER :: r, i, j, k, l, my_integer, m
   INTEGER, DIMENSION(5) :: arr
   INTEGER, DIMENSION(20) :: big_arr
   INTEGER :: ENDIF
   TYPE(test_type) :: t
   REAL(KIND=dp) :: r1, r2, r3, r4, r5, r6
   INTEGER, POINTER :: point

   point => NULL()

! 1) white space formatting !
!***************************!
! example 1.1
   r = 1; i = -2; j = 3; k = 4; l = 5
   r2 = 0.0_dp; r3 = 1.0_dp; r4 = 2.0_dp; r5 = 3.0_dp; r6 = 4.0_dp
   r1 = -(r2**i*(r3+r5*(-r4)-r6))-2.e+2
   IF (r .EQ. 2 .AND. r <= 5) i = 3
   WRITE (*, *) (MERGE(3, 1, i <= 2))
   WRITE (*, *) test_function(r, i, j, k)
   t%r = 4.0_dp
   t%i = str_function("t  % i   =  ")

! example 1.2
   my_integer = 2
   i = 3
   j = 5

   big_arr = [1, 2, 3, 4, 5, &
              6, 7, 8, 9, 10, &
              11, 12, 13, 14, 15, &
              16, 17, 18, 19, 20]

! example 1.3: disabling auto-formatter:
   my_integer = 2 !&
   i          = 3 !&
   j          = 5 !&

!&<
   my_integer = 2
   i          = 3
   j          = 5
!&>

   big_arr = [ 1,  2,  3,  4,  5, & !&
               6,  7,  8,  9, 10, & !&
              11, 12, 13, 14, 15, & !&
              16, 17, 18, 19, 20] !&

! example 1.4:

   big_arr = [1, 2, 3, 4, 5,&
           &  6, 7, 8, 9, 10, &
           & 11, 12, 13, 14, 15,&
            &16, 17, 18, 19, 20]

! 2) auto indentation for loops !
!*******************************!

! example 2.1
   l = 0
   DO r = 1, 10
      SELECT CASE (r)
      CASE (1)
         do_label: DO i = 1, 100
            IF (i <= 2) THEN
               m = 0
               DO WHILE (m < 4)
                  m = m+1
                  DO k = 1, 3
                     IF (k == 1) l = l+1
                  END DO
               ENDDO
            ENDIF
         ENDDO do_label
      CASE (2)
         l = i+j+k
      END SELECT
   ENDDO

! example 2.2
   DO m = 1, 2
      DO r = 1, 3
         WRITE (*, *) r
         DO k = 1, 4
         DO l = 1, 3
         DO i = 4, 5
            DO my_integer = 1, 1
            DO j = 1, 2
               WRITE (*, *) test_function(m, r, k, l)+i
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDDO
   ENDDO

! 3) auto alignment for linebreaks   !
!************************************!

! example 3.1
   l = test_function(1, 2, test_function(1, 2, 3, 4), 4)+3*(2+1)

   l = test_function(1, 2, test_function(1, 2, 3, 4), 4)+ &
       3*(2+1)

   l = test_function(1, 2, &
                     test_function(1, 2, 3, 4), 4)+ &
       3*(2+1)

   l = test_function(1, 2, &
                     test_function(1, 2, 3, &
                                   4), 4)+ &
       3*(2+1)

! example 3.2
   arr = [1, (/3, 4, 5/), 6]+[1, 2, 3, 4, 5]

   arr = [1, (/3, 4, 5/), &
          6]+[1, 2, 3, 4, 5]

   arr = [1, (/3, 4, 5/), &
          6]+ &
         [1, 2, 3, 4, 5]

   arr = [1, (/3, 4, &
               5/), &
          6]+ &
         [1, 2, 3, 4, 5]

! example 3.3
   l = test_function(1, 2, &
                     3, 4)

   l = test_function( &
       1, 2, 3, 4)

   arr = [1, 2, &
          3, 4, 5]
   arr = [ &
         1, 2, 3, 4, 5]

! 4) More complex formatting and tricky test cases !
!**************************************************!

! example 4.1
   l = 0
   DO r = 1, 10
      SELECT CASE (r)
      CASE (1)
         DO i = 1, 100; IF (i <= 2) THEN ! comment
               DO j = 1, 5
                  DO k = 1, 3
                     l = l+1
! unindented comment
                     ! indented comment
                  END DO; ENDDO
            ELSEIF (.NOT. j == 4) THEN
               my_integer = 4
            ELSE
               WRITE (*, *) "hello"
            ENDIF
         ENDDO
      CASE (2)
         l = i+j+k
      END SELECT
   ENDDO

! example 4.2
   IF ( &
      l == &
      111) &
      THEN
      DO k = 1, 2
         IF (k == 1) &
            l = test_function(1, &
                              test_function(r=4, i=5, &
                                            j=6, k=test_function(1, 2*(3*(1+1)), str_function(")a!(b['(;=dfe"), &
                                                                 9)+ &
                                            test_function(1, 2, 3, 4)), 9, 10) &
                ! test_function(1,2,3,4)),9,10) &
                ! +13*str_function('') + str_function('"')
                +13*str_function('')+str_function('"')
      END & ! comment
         ! comment
         DO
   ENDIF

! example 4.3
   arr = [1, (/3, 4, &
               5/), &
          6]+ &
         [1, 2, 3, 4, 5]; arr = [1, 2, &
 3, 4, 5]

! example 4.4
   ENDIF = 3
   IF (ENDIF == 2) THEN
      ENDIF = 5
   ELSE IF (ENDIF == 3) THEN
      WRITE (*, *) ENDIF
   ENDIF

! example 4.5
   DO i = 1, 2; IF (.TRUE.) THEN
         WRITE (*, *) "hello"
      ENDIF; ENDDO

END PROGRAM
"""

# EOF
