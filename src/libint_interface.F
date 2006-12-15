!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/libint_interface [1.0] *
!!
!!   NAME
!!     libint_interface
!!
!!   FUNCTION
!!     Interface to the Libint-Library 
!!
!!   NOTES
!!     IMPORTANT NOTE : this file currently is for a libint configured for
!!                      LIBINT_MAX_AM 5
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!****************************************************************************
MODULE libint_interface

#if defined(__LIBINT)
  USE, INTRINSIC :: ISO_C_BINDING
#endif
  USE gamma,                           ONLY: fgamma
  USE mathconstants
  USE orbital_pointers
  USE kinds,                           ONLY: dp

#include "cp_common_uses.h"      
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC initialize_libint, evaluate_eri, lib_int, terminate_libint

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_interface'
 
#if defined(__LIBINT)
  INTEGER, PARAMETER :: libint_max_am = 5
  INTEGER, PARAMETER :: build_eri_size = libint_max_am-1
  INTEGER, PARAMETER :: prim_data_f_size = 17
  INTEGER, PARAMETER :: lib_int_vrr_classes_size = 9

  TYPE, BIND(C) :: prim_data
    REAL(C_DOUBLE) :: F(prim_data_f_size)
    REAL(C_DOUBLE) :: U(3,6)
    REAL(C_DOUBLE) :: twozeta_a
    REAL(C_DOUBLE) :: twozeta_b
    REAL(C_DOUBLE) :: twozeta_c
    REAL(C_DOUBLE) :: twozeta_d
    REAL(C_DOUBLE) :: oo2z
    REAL(C_DOUBLE) :: oo2n
    REAL(C_DOUBLE) :: oo2zn
    REAL(C_DOUBLE) :: poz
    REAL(C_DOUBLE) :: pon
    REAL(C_DOUBLE) :: oo2p
    REAL(C_DOUBLE) :: ss_r12_ss
  END TYPE prim_data
        
  TYPE, BIND(C) :: lib_int
    TYPE(C_PTR)     :: int_stack
    TYPE(C_PTR)     :: PrimQuartet
    REAL(C_DOUBLE)  :: AB(3)
    REAL(C_DOUBLE)  :: CD(3)
    TYPE(C_PTR)     :: vrr_classes(lib_int_vrr_classes_size,lib_int_vrr_classes_size)
    TYPE(C_PTR)     :: vrr_stack
  END TYPE lib_int

  !global vars
  TYPE(C_FUNPTR), DIMENSION(0:build_eri_size,0:build_eri_size,0:build_eri_size,0:build_eri_size), BIND(C) :: build_eri

  INTEGER                            :: global_m = 0
  REAL(KIND=C_DOUBLE), DIMENSION(17) :: global_Fm = 0.0_dp
  REAL(KIND=C_DOUBLE)                :: global_T  = 0.0_dp

  INTERFACE
    FUNCTION build(lib, np) BIND(C)
      IMPORT
      TYPE(C_PTR)                    :: build
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION build
  
    FUNCTION init_lib(lib, max_am, np) BIND(C, name="init_libint")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_lib
      TYPE(lib_int)                  :: lib
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT),VALUE      :: np
    END FUNCTION INIT_LIB
    
    SUBROUTINE init_base() BIND(C, name="init_libint_base")
    END SUBROUTINE init_base

    SUBROUTINE free_library(lib) BIND(C, name="free_libint")
      IMPORT
      TYPE(lib_int)                  :: lib
    END SUBROUTINE free_library
  END INTERFACE

#else

  ! a bogus lib_int type to be able to compile without libint being present

  TYPE lib_int
     INTEGER :: do_not_use_me
  END TYPE lib_int

#endif

  CONTAINS
  
  SUBROUTINE initialize_libint(lib,max_am,error)
    TYPE(lib_int)                            :: lib
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libint', &
      routineP = moduleN//':'//routineN

    LOGICAL                          :: failure


#if defined(__LIBINT)
    INTEGER(KIND=C_INT)              :: max_am_local
    INTEGER(KIND=C_INT)              :: lib_storage 
    INTEGER(KIND=C_INT)              :: max_prim

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    CALL init_base()
    lib_storage = init_lib(lib, max_am_local, max_prim)

    IF (lib_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libint ", &
                        error,failure)
    ENDIF
#else 
    failure = .FALSE.
    CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                           " libint required, but CP2K not linked against libint" ,error,failure)

#endif

  END SUBROUTINE initialize_libint

  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                    :: lib
#if defined(__LIBINT)
     CALL free_library(lib)
#endif
  END SUBROUTINE terminate_libint

  SUBROUTINE evaluate_eri(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                          n_a,n_b,n_c,n_d,work_prim,&
                          work_prim_perm,omega)
    
    TYPE(lib_int)                           :: lib
    REAL(dp), INTENT(IN),TARGET             :: A(3), B(3), C(3), D(3)
    REAL(dp), INTENT(IN),TARGET             :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                     :: n_a, n_b, n_c, n_d
    REAL(dp), DIMENSION(:), POINTER         :: work_prim, work_prim_perm
    REAL(dp), DIMENSION(:), POINTER         :: p_work
    REAL(dp), INTENT(IN), OPTIONAL          :: omega

#if defined(__LIBINT)
   INTEGER                                  :: i,j,k,l,p_i,p_j,p_k,p_l,index,p_index,temp
   TYPE(prim_data),TARGET                   :: prim
   INTEGER(KIND=C_INT)                      :: Config(4)
   INTEGER(KIND=C_INT)                      :: max_am, m_max
   INTEGER(KIND=C_INT)                      :: mysize
   INTEGER(KIND=C_INT)                      :: n1, n2, n3, n4, n
   
   REAL(dp),DIMENSION(:),POINTER            :: A1,B1,C1,D1,X
   REAL(dp), POINTER                        :: Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, Y
   
   PROCEDURE(build), POINTER                :: pbuild
   TYPE(C_PTR)                              :: pc_result 
   




   !permutation of configuration

   
   
   n1=n_a; n2=n_b; n3=n_c; n4=n_d
   A1 => A; B1 => B; C1 => C; D1 => D;
   Zeta_A1 => Zeta_A; Zeta_B1 => Zeta_B; Zeta_C1 => Zeta_C; Zeta_D1 => Zeta_D
    
   If(n_a<n_b) THEN
     n1 = n_b
     n2 = n_a
     A1 => B
     B1 => A
     Zeta_A1 => Zeta_B
     Zeta_B1 => Zeta_A
   END IF

   IF(n_c<n_d) THEN
     n3 = n_d
     n4 = n_c
     C1 => D
     D1 => C
     Zeta_C1 => Zeta_D
     Zeta_D1 => Zeta_C
   END IF
 
   IF( n_a+n_b > n_c+n_d) THEN
     n   = n1
     n1  = n3
     n3  = n
     n   = n2
     n2  = n4
     n4  = n
     
     X => A1
     A1 => C1
     C1 => X
     X => B1
     B1 => D1
     D1 => X
     
     Y => Zeta_A1
     Zeta_A1 => Zeta_C1
     Zeta_C1 => Y
     y => Zeta_B1
     Zeta_B1 => Zeta_D1
     Zeta_D1 => Y  
   END IF
   
   Config=(/n1,n2,n3,n4/)
   max_am=MAX(Config(1), Config(2), Config(3), Config(4))
   m_max=Config(1) + Config(2) + Config(3) + Config(4)

   IF(PRESENT(omega)) THEN
     prim=build_quartet_data(A1,B1,C1,D1,Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, m_max, omega=omega)
   ELSE
     prim=build_quartet_data(A1,B1,C1,D1,Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, m_max)
   END IF

   lib%PrimQuartet=C_LOC(prim%F(1))

   lib%AB=A1-B1
   lib%CD=C1-D1
  
   mysize = nco(n1)*nco(n2)*nco(n3)*nco(n4)
     
   IF(n_a+ n_b + n_c + n_d  /=0) THEN
     CALL C_F_PROCPOINTER(build_eri(config(4),config(3),config(2),config(1)),pbuild)
     pc_result = pbuild(lib,1)
     CALL C_F_POINTER(pc_result, p_work, (/mysize/))
     work_prim(1:mysize) = p_work(1:mysize)
   ELSE
     work_prim(1) = prim%F(1)
   END IF
   DO l = 1,nco(n_d)
     DO k = 1,nco(n_c)     
       DO j = 1,nco(n_b)
         DO i = 1,nco(n_a)
           p_l = l
           p_i = i
           p_j = j
           p_k = k 

           IF(n_a<n_b) THEN
             temp = p_i
             p_i = p_j
             p_j = temp
           END IF
 
           IF(n_c<n_d) THEN
             temp = p_k
             p_k = p_l
             p_l = temp 
           END IF
           IF( n_a+n_b > n_c+n_d) THEN
             temp = p_i
             p_i = p_k
             p_k = temp
             temp = p_j
             p_j = p_l
             p_l = temp
           END IF
   
           index = (((i-1)*nco(n_b)+(j-1))*nco(n_c)+(k-1))*nco(n_d)+(l-1) + 1
           p_index = (((p_i-1)*nco(n2)+(p_j-1))*nco(n3)+(p_k-1))*nco(n4)+(p_l-1) + 1
           work_prim_perm(index) = work_prim(p_index)
         END DO
       END DO
     END DO
   END DO
#endif
  END SUBROUTINE evaluate_eri

#if defined(__LIBINT)
  FUNCTION build_quartet_data(A,B,C,D,Zeta_A, Zeta_B, Zeta_C, Zeta_D, m_max, omega)
     
    TYPE(prim_data)                    :: build_quartet_data
    REAL(KIND=C_DOUBLE), INTENT(IN)    :: A(3), B(3), C(3), D(3)
    REAL(KIND=C_DOUBLE), INTENT(IN)    :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER(KIND=C_INT), INTENT(IN)    :: m_max
    REAL(KIND=C_DOUBLE), INTENT(IN),&
                         OPTIONAL      :: omega
     
    TYPE(prim_data)                    :: prim
    REAL(KIND=C_DOUBLE)                :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3)
    REAL(KIND=C_DOUBLE)                :: Zeta, Eta, Rho
    REAL(KIND=C_DOUBLE)                :: AB2, CD2, PQ2,factor
    REAL(KIND=C_DOUBLE)                :: S12, S34, T, omega_corr
    REAL(KIND=C_DOUBLE), DIMENSION(17) :: Fm
    
    Zeta = Zeta_A + Zeta_B
    Eta  = Zeta_C + Zeta_D
    Rho  = Zeta*Eta/(Zeta+Eta)
    
    P = (Zeta_A*A + Zeta_B*B)/Zeta
    Q = (Zeta_C*C + Zeta_D*D)/Eta
    W = (Zeta*P + Eta*Q)/(Zeta+Eta)
    AB = A-B
    CD = C-D
    PQ = P-Q
    AB2 = DOT_PRODUCT(AB,AB)
    CD2 = DOT_PRODUCT(CD,CD)
    PQ2 = DOT_PRODUCT(PQ,PQ)

    S12  = EXP(-Zeta_A*Zeta_B/Zeta*AB2)
    S34  = EXP(-Zeta_C*Zeta_D/Eta*CD2)
    T = Rho*PQ2

    IF( T == global_T .AND. m_max == global_m) THEN
      prim%F = global_Fm
    ELSE 
      CALL fgamma(m_max,T,prim%F)
      global_Fm = prim%F
      global_T = T
      global_m = m_max
    END IF
    IF(PRESENT(omega)) THEN
      omega_corr = omega**2.0_dp/(omega**2.0_dp+Rho)
      T = T*omega_corr
      CALL fgamma(m_max,T,Fm)
      prim%F = prim%F - 2.0_dp*Pi/Rho*SQRT(omega_corr)*Fm
    END IF
    factor = 2.0_dp*Pi/Rho 
    factor = factor*S12*S34*SQRT((Pi/(Eta+Zeta))**(3.0_dp))

    CALL dscal(17,factor,prim%F,1)
    prim%U(:,1) = P-A
    prim%U(:,3) = Q-C
    prim%U(:,5) = W-P
    prim%U(:,6) = W-Q
    prim%twozeta_a = 0.0_dp                      !Not used in libint
    prim%twozeta_b = 0.0_dp                      !Not used in libint
    prim%twozeta_c = 0.0_dp                      !Not used in libint 
    prim%twozeta_d = 0.0_dp                      !Not used in libint
    prim%oo2z      = 1.0_dp/(2.0_dp*Zeta)
    prim%oo2n      = 1.0_dp/(2.0_dp*Eta)
    prim%oo2zn     = 1.0_dp/(2.0_dp*(Zeta+Eta))
    prim%poz       = Rho/Zeta
    prim%pon       = Rho/Eta
    prim%oo2p      = 1.0_dp/(2.0_dp*Rho)
    prim%ss_r12_ss = 0.0_dp                      !Not used in libint

    build_quartet_data = prim
  END FUNCTION build_quartet_data
#endif 

END MODULE libint_interface
