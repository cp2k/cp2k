!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
!!                      LIBINT_MAX_AM 5 and LIBINT_MAX_AM1 4
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
  PUBLIC initialize_libint, evaluate_eri, lib_int, terminate_libint,&
         initialize_libderiv, lib_deriv, terminate_libderiv, evaluate_deriv_eri

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'libint_interface'
 
#if defined(__LIBINT)
  INTEGER, PARAMETER :: libint_max_am = 5
  INTEGER, PARAMETER :: build_eri_size = libint_max_am-1
  INTEGER, PARAMETER :: prim_data_f_size = 17
  INTEGER, PARAMETER :: lib_int_vrr_classes_size = 9
  INTEGER, PARAMETER :: libderiv_max_am1 = 4
  INTEGER, PARAMETER :: build_deriv1_eri_size = libderiv_max_am1 - 1 

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

  TYPE, BIND(C) :: lib_deriv
    TYPE(C_PTR)     :: int_stack
    TYPE(C_PTR)     :: PrimQuartet
    TYPE(C_PTR)     :: zero_stack
    TYPE(C_PTR)     :: ABCD(156)
    REAL(C_DOUBLE)  :: AB(3)
    REAL(C_DOUBLE)  :: CD(3)
    TYPE(C_PTR)     :: deriv_classes(12,7,7)
    TYPE(C_PTR)     :: deriv2_classes(144,7,7)
    TYPE(C_PTR)     :: dvrr_classes(7,7)
    TYPE(C_PTR)     :: dvtt_stack
  END TYPE lib_deriv

  !global vars
  TYPE(C_FUNPTR), DIMENSION(0:build_eri_size,0:build_eri_size,0:build_eri_size,0:build_eri_size), BIND(C) :: build_eri
  TYPE(C_FUNPTR), DIMENSION(0:build_deriv1_eri_size,0:build_deriv1_eri_size,&
                            0:build_deriv1_eri_size,0:build_deriv1_eri_size), BIND(C) :: build_deriv1_eri

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
      INTEGER(KIND=C_INT), VALUE     :: np
    END FUNCTION init_lib
    
    SUBROUTINE init_base() BIND(C, name="init_libint_base")
    END SUBROUTINE init_base

    SUBROUTINE free_lib_int(lib) BIND(C, name="free_libint")
      IMPORT
      TYPE(lib_int)                  :: lib
    END SUBROUTINE free_lib_int

    SUBROUTINE init_deriv_base() BIND(C, name="init_libderiv_base")
    END SUBROUTINE init_deriv_base

    FUNCTION init_deriv1(deriv, max_am, np, ccs) BIND(C, name="init_libderiv1")
      IMPORT
      INTEGER(KIND=C_INT)            :: init_deriv1
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT), VALUE     :: max_am
      INTEGER(KIND=C_INT), VALUE     :: np
      INTEGER(KIND=C_INT), VALUE     :: ccs
    END FUNCTION init_deriv1

    SUBROUTINE build_deriv1(deriv, np) BIND(C)
      IMPORT
      TYPE(lib_deriv)                :: deriv
      INTEGER(KIND=C_INT),VALUE      :: np
    END SUBROUTINE build_deriv1

    SUBROUTINE free_lib_deriv(deriv) BIND(C, name="free_libderiv")
      IMPORT
      TYPE(lib_deriv)                  :: deriv
    END SUBROUTINE free_lib_deriv

  END INTERFACE

#else

  ! a bogus lib_int type to be able to compile without libint being present

  TYPE lib_int
     INTEGER :: do_not_use_me
  END TYPE lib_int
  
  TYPE lib_deriv
     INTEGER :: do_not_use_me
  END TYPE lib_deriv

#endif

!***
!****************************************************************************

  CONTAINS

!!***** libint_interface/inititalize_libint [1.0] *
!!
!!   NAME
!!     initialize_libint
!!
!!   FUNCTION
!!     Allocates memory for using libint
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!!****************************************************************************


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
!***************************************************************************


!!***** libint_interface/inititalize_libderiv [1.0] *
!!
!!   NAME
!!     initialize_libderiv
!!
!!   FUNCTION
!!     Allocates memory for using libderiv
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE initialize_libderiv(deriv,max_am,error)
    TYPE(lib_deriv)                          :: deriv
    INTEGER                                  :: max_am
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_libderiv', &
      routineP = moduleN//':'//routineN

    LOGICAL                          :: failure


#if defined(__LIBINT)
    INTEGER(KIND=C_INT)              :: max_am_local
    INTEGER(KIND=C_INT)              :: lib_deriv_storage
    INTEGER(KIND=C_INT)              :: max_prim
    INTEGER(KIND=C_INT)              :: max_classes

    failure = .FALSE.
    max_am_local= max_am
    max_prim = 1
    max_classes = 10000!max_am**4

    CALL init_deriv_base()
    lib_deriv_storage = init_deriv1(deriv, max_am_local, max_prim, max_classes)

    IF (lib_deriv_storage<0) THEN
        CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                        " the angular momentum needed exceeds the value assumed when configuring libderiv ", &
                        error,failure)
    ENDIF
#else
    failure = .FALSE.
    CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                           " libint required, but CP2K not linked against libint" ,error,failure)

#endif

  END SUBROUTINE initialize_libderiv
!***************************************************************************

!!***** libint_interface/terminate_libint [1.0] *
!!
!!   NAME
!!     terminate_libint
!!
!!   FUNCTION
!!     Frees allocated memory in libint
!!
!!   NOTES
!!     Does not work at the moment
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE terminate_libint(lib)
    TYPE(lib_int)                    :: lib
#if defined(__LIBINT)
     CALL free_lib_int(lib)
#endif
  END SUBROUTINE terminate_libint
!***************************************************************************

!!***** libint_interface/terminate_libderiv [1.0] *
!!
!!   NAME
!!     terminate_libderiv
!!
!!   FUNCTION
!!     Frees allocated memory in libderiv
!!
!!   NOTES
!!     Does not work at the moment
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE terminate_libderiv(deriv)
    TYPE(lib_deriv)                    :: deriv
#if defined(__LIBINT)
     CALL free_lib_deriv(deriv)
#endif
  END SUBROUTINE terminate_libderiv
!***************************************************************************

!!***** libint_interface/evaluate_eri [1.0] *
!!
!!   NAME
!!     evaluate_eri
!!
!!   FUNCTION
!!     Evaluate electron repulsion integrals for a primitive quartet
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!!****************************************************************************


  SUBROUTINE evaluate_eri(lib,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                          n_a,n_b,n_c,n_d,work_prim,&
                          work_prim_perm,omega,deriv)
    
    TYPE(lib_int)                           :: lib
    REAL(dp), INTENT(IN),TARGET             :: A(3), B(3), C(3), D(3)
    REAL(dp), INTENT(IN),TARGET             :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                     :: n_a, n_b, n_c, n_d
    REAL(dp), DIMENSION(:), POINTER         :: work_prim, work_prim_perm
    REAL(dp), DIMENSION(:), POINTER         :: p_work
    REAL(dp), INTENT(IN), OPTIONAL          :: omega
    TYPE(lib_deriv), OPTIONAL               :: deriv

#if defined(__LIBINT)
    INTEGER                                  :: i,j,k,l,p_i,p_j,p_k,p_l,index,p_index,temp
    TYPE(prim_data),TARGET                   :: prim
    INTEGER(KIND=C_INT)                      :: Config(4)
    INTEGER(KIND=C_INT)                      :: m_max
    INTEGER(KIND=C_INT)                      :: mysize
    INTEGER(KIND=C_INT)                      :: n1, n2, n3, n4, n
   
    REAL(dp),DIMENSION(:),POINTER            :: A1,B1,C1,D1,X
    REAL(dp), POINTER                        :: Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, Y
   
    PROCEDURE(build), POINTER                :: pbuild
    TYPE(C_PTR)                              :: pc_result 
   
    !permutation of configuration

    n1 = n_a; n2 = n_b; n3 = n_c; n4 = n_d
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
   
    Config = (/n1,n2,n3,n4/)
    m_max = Config(1) + Config(2) + Config(3) + Config(4)


    IF(PRESENT(omega)) THEN
      prim = build_quartet_data(A1,B1,C1,D1,Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, m_max, omega=omega)
    ELSE
      prim = build_quartet_data(A1,B1,C1,D1,Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, m_max)
    END IF

    lib%PrimQuartet = C_LOC(prim)

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
!***************************************************************************

!!***** libint_interface/build_quartet_data [1.0] *
!!
!!   NAME
!!     build_quartet_data
!!
!!   FUNCTION
!!     Fill data structure used in libint and libderiv
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

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
    prim%U(:,2) = P-B                 !Not used in libint
    prim%U(:,3) = Q-C
    prim%U(:,4) = Q-D                 !Not used in libint
    prim%U(:,5) = W-P
    prim%U(:,6) = W-Q
    prim%twozeta_a = 2.0_dp*Zeta_A     !0.0_dp                      !Not used in libint
    prim%twozeta_b = 2.0_dp*Zeta_B     !0.0_dp                      !Not used in libint
    prim%twozeta_c = 2.0_dp*Zeta_C     !0.0_dp                      !Not used in libint 
    prim%twozeta_d = 2.0_dp*Zeta_D     !                      !Not used in libint
    prim%oo2z      = 1.0_dp/(2.0_dp*Zeta)
    prim%oo2n      = 1.0_dp/(2.0_dp*Eta)
    prim%oo2zn     = 1.0_dp/(2.0_dp*(Zeta+Eta))
    prim%poz       = Rho/Zeta
    prim%pon       = Rho/Eta
    prim%oo2p      = 1.0_dp/(2.0_dp*Rho)
    prim%ss_r12_ss = 0.0_dp                      !Not used in libint,libderiv

    build_quartet_data = prim
  END FUNCTION build_quartet_data
#endif 
!***************************************************************************

!!***** libint_interface/evaluate_deriv_eri [1.0] *
!!
!!   NAME
!!     evaluate_deriv_eri
!!
!!   FUNCTION
!!     Evaluates derivatives of electron repulsion integrals for a primitive quartet
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************


  SUBROUTINE evaluate_deriv_eri(deriv,A,B,C,D,Zeta_A,Zeta_B,Zeta_C,Zeta_D,&
                                n_a,n_b,n_c,n_d,work_forces,&
                                work_forces_perm,omega)
    
    TYPE(lib_deriv)                         :: deriv
    REAL(dp), INTENT(IN),TARGET             :: A(3), B(3), C(3), D(3)
    REAL(dp), INTENT(IN),TARGET             :: Zeta_A, Zeta_B, Zeta_C, Zeta_D
    INTEGER, INTENT(IN)                     :: n_a, n_b, n_c, n_d
    REAL(dp), DIMENSION(:,:), POINTER       :: work_forces, work_forces_perm
    REAL(dp), INTENT(IN), OPTIONAL          :: omega

#if defined(__LIBINT)
    INTEGER                                  :: i,j,k,l,p_i,p_j,p_k,p_l,index,p_index,temp
    TYPE(prim_data),TARGET                   :: prim
    INTEGER(KIND=C_INT)                      :: Config(4)
    INTEGER(KIND=C_INT)                      :: m_max
    INTEGER(KIND=C_INT)                      :: mysize
    INTEGER(KIND=C_INT)                      :: n1, n2, n3, n4, n
   
    REAL(dp),DIMENSION(:),POINTER            :: A1,B1,C1,D1,X
    REAL(dp), POINTER                        :: Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, Y
    
    TYPE(C_PTR)                              :: pc_result 
   
    PROCEDURE(build_deriv1), POINTER         :: pbuild_deriv1
    REAL(C_DOUBLE), DIMENSION(:), POINTER    :: tmp_data
   
    INTEGER                                  :: perm(4)
    !permutation of configuration

    perm=(/1,2,3,4/)
    n1 = n_a; n2 = n_b; n3 = n_c; n4 = n_d
    A1 => A; B1 => B; C1 => C; D1 => D;
    Zeta_A1 => Zeta_A; Zeta_B1 => Zeta_B; Zeta_C1 => Zeta_C; Zeta_D1 => Zeta_D
    
    If(n_a<n_b) THEN
      n1 = n_b
      n2 = n_a
      A1 => B
      B1 => A
      Zeta_A1 => Zeta_B
      Zeta_B1 => Zeta_A
      perm(1) = 2
      perm(2) = 1
    END IF

    IF(n_c<n_d) THEN
      n3 = n_d
      n4 = n_c
      C1 => D
      D1 => C
      Zeta_C1 => Zeta_D
      Zeta_D1 => Zeta_C
      perm(3)=4
      perm(4)=3
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
      n = perm(1)
      perm(1) = perm(3)
      perm(3) = n

      n = perm(2)
      perm(2) = perm(4)
      perm(4) = n

    END IF

   
    Config = (/n1,n2,n3,n4/)
    m_max = Config(1) + Config(2) + Config(3) + Config(4)

    m_max = m_max + 1
    IF(PRESENT(omega)) THEN
      prim = build_quartet_data(A1,B1,C1,D1,Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, m_max, omega=omega)
    ELSE
      prim = build_quartet_data(A1,B1,C1,D1,Zeta_A1, Zeta_B1, Zeta_C1, Zeta_D1, m_max)
    END IF

    mysize = nco(n1)*nco(n2)*nco(n3)*nco(n4)
  

    deriv%PrimQuartet = C_LOC(prim)
    deriv%AB = A1-B1
    deriv%CD = C1-D1
    CALL C_F_PROCPOINTER(build_deriv1_eri(config(4),config(3),config(2),config(1)),pbuild_deriv1)
    CALL pbuild_deriv1(deriv,1)

    DO k=1,12
      pc_result = deriv%ABCD(k)
      CALL C_F_POINTER(pc_result, tmp_data ,(/mysize/))
      IF (TRANSFER(pc_result,0_C_INTPTR_T).NE.0) THEN
        work_forces(1:mysize,k) = tmp_data
      ENDIF
    END DO

    DO k=4,6
     work_forces(1:mysize,k) = - 1.0_dp* (work_forces(1:mysize,k-3) + &
                                          work_forces(1:mysize,k+3) + &
                                          work_forces(1:mysize,k+6) )
    END DO
 
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
            DO n=1,4
              work_forces_perm(index,(perm(n)-1)*3+1:(perm(n)-1)*3+3) = work_forces(p_index,(n-1)*3+1:(n-1)*3+3)
            END DO
          END DO
        END DO
      END DO
    END DO

#endif
  END SUBROUTINE evaluate_deriv_eri
!***************************************************************************

END MODULE libint_interface
