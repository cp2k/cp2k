!-------------------------------------------------------------------------------
! some utils to convert and do simple stuff
!-------------------------------------------------------------------------------
MODULE utils
  IMPLICIT NONE
CONTAINS
!-------------------------------------------------------------------------------
! upcases a string
!-------------------------------------------------------------------------------
   FUNCTION toupper(label) RESULT(RES)
     CHARACTER(LEN=*) :: label
     CHARACTER(LEN=LEN(label)) ::  RES
     INTEGER I,J
     DO I=1,LEN(label)
       J=INDEX("the quick brown fox jumps over the lazy dog",label(I:I))
       SELECT CASE(J)
       CASE (0)
         RES(I:I)=label(I:I)
       CASE DEFAULT
         RES(I:I)="THE QUICK BROWN FOX JUMPS OVER THE LAZY DOG"(J:J)
       END SELECT
     ENDDO
   END FUNCTION toupper
!-------------------------------------------------------------------------------
! folds the coordates into a box with edges and center (not origin) 
!-------------------------------------------------------------------------------
   FUNCTION fold(coords,center,edges)
     REAL, DIMENSION(3) :: fold,coords,center,edges
     fold=coords-edges*ANINT((1.0/edges)*(coords-center))
   END FUNCTION
END MODULE

