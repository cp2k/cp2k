MODULE hip_profiling

 INTERFACE 
     SUBROUTINE roctxMarkA(message) BIND(c, name="roctxMarkA")
       USE ISO_C_BINDING,   ONLY: C_CHAR
       IMPLICIT NONE
       CHARACTER(C_CHAR) :: message(*)
     END SUBROUTINE roctxMarkA

     FUNCTION roctxRangePushA(message) BIND(c, name="roctxRangePushA")
       USE ISO_C_BINDING,   ONLY: C_INT,&
                                  C_CHAR
       IMPLICIT NONE
       INTEGER(C_INT) :: roctxRangePushA
       CHARACTER(C_CHAR) :: message(*)
     END FUNCTION roctxRangePushA

     SUBROUTINE roctxRangePop() BIND(c, name="roctxRangePop")
       IMPLICIT NONE
     END SUBROUTINE roctxRangePop

 END INTERFACE

END MODULE hip_profiling

