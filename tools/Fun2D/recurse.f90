! * Copyright (c) 2008, Joost VandeVondele and Manuel Guidon
! * All rights reserved.
! *
! * Redistribution and use in source and binary forms, with or without
! * modification, are permitted provided that the following conditions are met:
! *     * Redistributions of source code must retain the above copyright
! *       notice, this list of conditions and the following disclaimer.
! *     * Redistributions in binary form must reproduce the above copyright
! *       notice, this list of conditions and the following disclaimer in the
! *       documentation and/or other materials provided with the distribution.
! *
! * THIS SOFTWARE IS PROVIDED BY Joost VandeVondele and Manuel Guidon ''AS IS'' AND ANY
! * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! * DISCLAIMED. IN NO EVENT SHALL Joost VandeVondele or Manuel Guidon BE LIABLE FOR ANY
! * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!   Joost VandeVondele and Manuel Guidon, Nov 2008.
!

!   recursively bisect the rectangular domain of a bivariate function
!   till the lagrange interpolation of specified degree through
!   the Padua points for it converges to a preset threshold
!   generate a Fortran module for the efficient evaluation
!  
!   This version computes the basic integrals for the
!   truncated coulomb operator
!
!   G_0(R,T)= ((2*erf(sqrt(t))+erf(R-sqrt(t))-erf(R+sqrt(t)))/sqrt(t))
!
!   and up to 21 derivatives with respect to T
!
!   (-1)**n d^n/dT^n G_0(R,T)
!
!
!   The function is only computed for values of R,T which fulfil
!
!   R**2 - 11.0_dp*R + 0.0_dp < T < R**2 + 11.0_dp*R + 50.0_dp
!   where R>=0 T>=0
! 
!   for T larger than the upper bound, 0 is returned
!   (which is accurate at least up to 1.0E-16)
!   while for T smaller than the lower bound, the caller is instructed 
!   to use the gamma function instead
!   This program requires the auxiliary program drvrec_mpfr.x
!
    MODULE RECURSE      
      INTEGER, SAVE :: find_count=0
      INTEGER, PARAMETER :: dp=KIND(0.0D0)
      REAL(KIND=dp) ERR
      integer :: functionid,digits
    CONTAINS

!====================================================================================================================
!
!     a function specific wrapper to the bisection routines
!
!     this function must open a unit 13 (the coefficients of the interpolation)
!     this function must open a unit 17 (some info on the bisected domains )
!
!====================================================================================================================
      SUBROUTINE GO()
         IMPLICIT NONE
         REAL(KIND=dp) A1,B1,A2,B2
         INTEGER :: DEG,LEVEL
         INTEGER :: nderiv


         ! these are the selected settings for CP2K
         ! some initial guess for the needed number of digits
         digits=128
         ! number of simultaneous function values to be evaluated for a given 2D point 
         nderiv=21
         ! degree of the interpolating polynomial
         DEG=13
         ! target error
         ERR=1.0E-9

         OPEN(UNIT=13,FILE="T_C_G.dat")
         LEVEL=0

         CALL write_copy_right()

         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   This module computes the basic integrals for the"
         WRITE(6,'(A)') "!   truncated coulomb operator"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   res(1) =G_0(R,T)= ((2*erf(sqrt(t))+erf(R-sqrt(t))-erf(R+sqrt(t)))/sqrt(t))"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   and up to 21 derivatives with respect to T"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   res(n+1)=(-1)**n d^n/dT^n G_0(R,T)"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   The function is only computed for values of R,T which fulfil"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   R**2 - 11.0_dp*R + 0.0_dp < T < R**2 + 11.0_dp*R + 50.0_dp"
         WRITE(6,'(A)') "!   where R>=0 T>=0"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   for T larger than the upper bound, 0 is returned"
         WRITE(6,'(A)') "!   (which is accurate at least up to 1.0E-16)"
         WRITE(6,'(A)') "!   while for T smaller than the lower bound, the caller is instructed"
         WRITE(6,'(A)') "!   to use the conventional gamma function instead"
         WRITE(6,'(A)') "!   (i.e. the limit of above expression for R to Infinity)"
         WRITE(6,'(A)') "! instead of a module one could use: INTEGER, PARAMETER ::  dp=KIND(0.0D0)"

         CALL write_module_head("T_C_G0",DEG,nderiv,err)

         WRITE(6,'(A,I0,A)') "SUBROUTINE T_C_G0_n(RES,use_gamma,R,T,NDERIV)"
         WRITE(6,'(A)')      " IMPLICIT NONE"
         WRITE(6,'(A)')      " REAL(KIND=dp), INTENT(OUT) :: RES(*)"
         WRITE(6,'(A)')      " LOGICAL, INTENT(OUT) :: use_gamma"
         WRITE(6,'(A)')      " REAL(KIND=dp),INTENT(IN) :: R,T"
         WRITE(6,'(A)')      " INTEGER, INTENT(IN)      :: NDERIV"
         WRITE(6,'(A)')      " REAL(KIND=dp)            :: upper,lower,X1,X2,TG1,TG2"

         WRITE(6,'(A)')      "  use_gamma=.FALSE."
         WRITE(6,'(A)')      "  upper=R**2 + 11.0_dp*R + 50.0_dp"
         WRITE(6,'(A)')      "  lower=R**2 - 11.0_dp*R + 0.0_dp"
         WRITE(6,'(A)')      "   IF (T>upper) THEN"
         WRITE(6,'(A)')      "      RES(1:NDERIV+1)=0.0_dp"
         WRITE(6,'(A)')      "      RETURN"
         WRITE(6,'(A)')      "   ENDIF"

         WRITE(6,'(A)')      " IF (R<=11.0_dp) THEN"
         WRITE(6,'(A)')      "  X2=R/11.0_dp"
         WRITE(6,'(A)')      "  upper=R**2 + 11.0_dp*R + 50.0_dp"
         WRITE(6,'(A)')      "  lower=0.0_dp"
         WRITE(6,'(A)')      "  X1=(T-lower)/(upper-lower)"

             functionid=2
             open(UNIT=17,FILE="nearfield_domains.dat")
             A1=0.0D0
             B1=1.0D0
             A2=0.0D0
             B2=1.0D0
             CALL FIND_DOMAINS(A1,B1,A2,B2,DEG,NDERIV,ERR,LEVEL)
             close(17)

         WRITE(6,'(A)') " ELSE"

         WRITE(6,'(A)') "   IF (T<lower) THEN"
         WRITE(6,'(A)') "      use_gamma=.TRUE."
         WRITE(6,'(A)') "      RETURN"
         WRITE(6,'(A)') "   ENDIF"
         WRITE(6,'(A)') "  X2=11.0_dp/R"
         WRITE(6,'(A)') "  X1=(T-lower)/(upper-lower)"

             functionid=1
             open(UNIT=17,FILE="farfield_domains.dat")
             A1=0.0D0
             B1=1.0D0
             A2=0.0D0
             B2=1.0D0
             CALL FIND_DOMAINS(A1,B1,A2,B2,DEG,NDERIV,ERR,LEVEL)
             close(17)

         WRITE(6,'(A)') " ENDIF"

         WRITE(6,'(A)') "END SUBROUTINE T_C_G0_n"
           
         CLOSE(UNIT=13)

         CALL WRITE_INIT()

         CALL WRITE_PD2VAL(DEG)

         WRITE(6,'(A)') "END MODULE T_C_G0"

      END SUBROUTINE GO
!====================================================================================================================
!
!     write the top of the module 
!
!====================================================================================================================
      SUBROUTINE WRITE_MODULE_HEAD(mod_name,DEG,nderiv,err)
         CHARACTER(LEN=*) :: mod_name
         INTEGER          :: nderiv,deg
         REAL(KIND=dp)    :: err

         WRITE(6,'(A,A)')    "MODULE ",mod_name
         WRITE(6,'(A)')      " USE kinds, ONLY: dp"
         WRITE(6,'(A)')      " IMPLICIT NONE"
         WRITE(6,'(A)')      " REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE,  SAVE :: C0"
         WRITE(6,'(A,I0)')   " INTEGER, PARAMETER :: degree=",DEG
         WRITE(6,'(A,E18.6)')" REAL(KIND=dp), PARAMETER :: target_error=",err
         WRITE(6,'(A,I0)')   " INTEGER, PARAMETER :: nderiv_max=",nderiv
         WRITE(6,'(A)')      " INTEGER            :: nderiv_init=-1"
         WRITE(6,'(A)')      " INTEGER, SAVE      :: patches=-1"
         WRITE(6,'(A)')      "CONTAINS"

      END SUBROUTINE WRITE_MODULE_HEAD
!====================================================================================================================
!
!     writes a BSD style copy right header
!
!====================================================================================================================
      SUBROUTINE WRITE_COPY_RIGHT()

         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "! * Copyright (c) 2008, Joost VandeVondele and Manuel Guidon"
         WRITE(6,'(A)') "! * All rights reserved."
         WRITE(6,'(A)') "! *"
         WRITE(6,'(A)') "! * Redistribution and use in source and binary forms, with or without"
         WRITE(6,'(A)') "! * modification, are permitted provided that the following conditions are met:"
         WRITE(6,'(A)') "! *     * Redistributions of source code must retain the above copyright"
         WRITE(6,'(A)') "! *       notice, this list of conditions and the following disclaimer."
         WRITE(6,'(A)') "! *     * Redistributions in binary form must reproduce the above copyright"
         WRITE(6,'(A)') "! *       notice, this list of conditions and the following disclaimer in the"
         WRITE(6,'(A)') "! *       documentation and/or other materials provided with the distribution."
         WRITE(6,'(A)') "! *"
         WRITE(6,'(A)') "! * THIS SOFTWARE IS PROVIDED BY Joost VandeVondele and Manuel Guidon AS IS AND ANY"
         WRITE(6,'(A)') "! * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED"
         WRITE(6,'(A)') "! * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE"
         WRITE(6,'(A)') "! * DISCLAIMED. IN NO EVENT SHALL Joost VandeVondele or Manuel Guidon BE LIABLE FOR ANY"
         WRITE(6,'(A)') "! * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES"
         WRITE(6,'(A)') "! * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;"
         WRITE(6,'(A)') "! * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND"
         WRITE(6,'(A)') "! * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT"
         WRITE(6,'(A)') "! * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS"
         WRITE(6,'(A)') "! * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!"
         WRITE(6,'(A)') "!   Joost VandeVondele and Manuel Guidon, Nov 2008."

      END SUBROUTINE WRITE_COPY_RIGHT
!====================================================================================================================
!
!     write the table initialization routines
!
!====================================================================================================================
      SUBROUTINE WRITE_INIT()

         WRITE(6,'(A)') "! iunit contains the data file to initialize the table"
         WRITE(6,'(A)') "! Nder is the number of derivatives that will actually be used"
         WRITE(6,'(A)') "SUBROUTINE INIT(Nder,iunit)"
         WRITE(6,'(A)') "  IMPLICIT NONE"
         WRITE(6,'(A)') "  INTEGER, INTENT(IN) :: Nder,iunit"
         WRITE(6,'(A)') "  INTEGER :: I"
         WRITE(6,'(A)') "  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: chunk"
         WRITE(6,'(A,I0)') "  patches=",find_count
         WRITE(6,'(A)') '  IF (Nder>nderiv_max) STOP "Reading data for initialization of C0 failed"'
         WRITE(6,'(A)') "  nderiv_init=Nder"
         WRITE(6,'(A)') "  IF(ALLOCATED(C0)) DEALLOCATE(C0)"
         WRITE(6,'(A)') "  ! round up to a multiple of 32 to give some generous alignment for each C0"
         WRITE(6,'(A)') "  ALLOCATE(C0(32*((31+(Nder+1)*(degree+1)*(degree+2)/2)/32),patches))"
         WRITE(6,'(A)') "  ALLOCATE(chunk((nderiv_max+1)*(degree+1)*(degree+2)/2))"
         WRITE(6,'(A)') "  DO I=1,patches"
         WRITE(6,'(A)') "     READ(iunit,*) chunk"
         WRITE(6,'(A)') "     C0(1:(Nder+1)*(degree+1)*(degree+2)/2,I)=chunk(1:(Nder+1)*(degree+1)*(degree+2)/2)"
         WRITE(6,'(A)') "  ENDDO"
         WRITE(6,'(A)') "  DEALLOCATE(chunk)"
         WRITE(6,'(A)') "END SUBROUTINE INIT"
         WRITE(6,'(A)') "SUBROUTINE FREE()"
         WRITE(6,'(A)') "  IF(ALLOCATED(C0)) DEALLOCATE(C0)"
         WRITE(6,'(A)') "  nderiv_init=-1"
         WRITE(6,'(A)') "END SUBROUTINE FREE"
      END SUBROUTINE WRITE_INIT
!====================================================================================================================
!
!     write the code to evaluate the interpolation
!
!====================================================================================================================
      SUBROUTINE WRITE_PD2VAL(DEG)

         INTEGER :: DEG

         INTEGER :: I,J

         WRITE(6,'(A)') "SUBROUTINE PD2VAL(RES,NDERIV,TG1,TG2,C0) "
         WRITE(6,'(A)') " IMPLICIT NONE"
         WRITE(6,'(A)') " REAL(KIND=dp), INTENT(OUT) :: res(*)"
         WRITE(6,'(A)') " INTEGER, INTENT(IN) :: NDERIV"
         WRITE(6,'(A)') " REAL(KIND=dp),INTENT(IN) :: TG1,TG2"
         WRITE(6,'(A,I0,A)') " REAL(KIND=dp) :: T1(0:",DEG,")"
         WRITE(6,'(A,I0,A)') " REAL(KIND=dp) :: T2(0:",DEG,")"
         WRITE(6,'(A,I0,A)') " REAL(KIND=dp), INTENT(IN) :: C0(",(DEG+1)*(DEG+2)/2,",*)"
         WRITE(6,'(A)') " INTEGER :: I,J,K"
         WRITE(6,'(A)') " REAL(KIND=dp), PARAMETER :: SQRT2=1.4142135623730950488016887242096980785696718753_dp"
         WRITE(6,'(A)') " T1(0)=1.0_dp"
         WRITE(6,'(A)') " T2(0)=1.0_dp"
         WRITE(6,'(A)') " T1(1)=SQRT2 * TG1"
         WRITE(6,'(A)') " T2(1)=SQRT2 * TG2"
         WRITE(6,'(A)') " T1(2)= 2 * TG1 * T1(1) - SQRT2 "
         WRITE(6,'(A)') " T2(2)= 2 * TG2 * T2(1) - SQRT2 "
         DO I=3,DEG
         WRITE(6,'(A,I0,A,I0,A,I0,A)') " T1(",I,") = 2 * TG1 * T1(",I-1,") - T1(",I-2,")" 
         WRITE(6,'(A,I0,A,I0,A,I0,A)') " T2(",I,") = 2 * TG2 * T2(",I-1,") - T2(",I-2,")" 
         ENDDO

         WRITE(6,'(A)')      " DO K=1,NDERIV+1"
         WRITE(6,'(A)')      "    RES(K) = 0.0_dp"
           J=1
           DO I=DEG,0,-1
              WRITE(6,'(A,I0,A,I0,A,I0,A,I0,A)') "  RES(K)=RES(K)+DOT_PRODUCT(T1(0:",I,"),C0(",&
                                                           J,":",J+I,",K))*T2(",DEG-I,")"
              J=J+I+1
           ENDDO
         WRITE(6,'(A)')      " ENDDO"

         WRITE(6,'(A)') "END SUBROUTINE PD2VAL"

      END SUBROUTINE WRITE_PD2VAL
!====================================================================================================================
!
!     The basic piece of the code, does the bisection till convergence
!
!====================================================================================================================
      RECURSIVE SUBROUTINE FIND_DOMAINS(A1,B1,A2,B2,DEG,NDERIV,ERR,LEVEL)

         IMPLICIT NONE
         REAL(KIND=dp) A1,B1,A2,B2,ERR,R1,R2
         INTEGER :: DEG,LEVEL,I,DIRECTION,NDERIV,K,I1,I2

         REAL(KIND=dp) ESTIMATED_ERR,MIDDLE,ESTIMATED_ERR1,ESTIMATED_ERR2,ESTIMATED_ERR3,ESTIMATED_ERR4
         REAL(KIND=dp) :: C0(0:DEG,0:DEG,0:NDERIV)

         IF (DEG<3) STOP "Not supported"
         IF (A1==B1 .AND. A2==B2) STOP "BUG 1"
         IF (LEVEL>20) STOP "Too many bisections: is the function sufficiently differentiable in the full domain?"

         ! compute the interpolating polynomial and its estimated error

         CALL drvrec_wrap(A1,B1,A2,B2,DEG,NDERIV,C0,ESTIMATED_ERR)

         IF (ESTIMATED_ERR<ERR) THEN

            ! we have converged ... turn this into a call to pd2val, and append to the C0 table
            find_count=find_count+1
            WRITE(6,'(A,E25.18,A,E25.18,A)') REPEAT(" ",LEVEL+3)//"TG1= ( 2 * X1 - ",(B1+A1),"_dp)*",1.0D0/(B1-A1),"_dp"
            WRITE(6,'(A,E25.18,A,E25.18,A)') REPEAT(" ",LEVEL+3)//"TG2= ( 2 * X2 - ",(B2+A2),"_dp)*",1.0D0/(B2-A2),"_dp"
          
            WRITE(6,'(A,I0,A)') REPEAT(" ",LEVEL+3)//"CALL PD2VAL(RES,NDERIV,TG1,TG2,C0(1,", find_count, "))"
            WRITE(17,'(2F20.16,I10,E16.8)') A1,A2,find_count,ESTIMATED_ERR
            WRITE(17,'(2F20.16,I10,E16.8)') A1,B2 ,find_count,ESTIMATED_ERR
            WRITE(17,*)
            WRITE(17,'(2F20.16,I10,E16.8)') A1,A2 ,find_count,ESTIMATED_ERR
            WRITE(17,'(2F20.16,I10,E16.8)') B1,A2 ,find_count,ESTIMATED_ERR
            WRITE(17,*)
            WRITE(17,'(2F20.16,I10,E16.8)') B1,B2,find_count,ESTIMATED_ERR
            WRITE(17,'(2F20.16,I10,E16.8)') A1,B2,find_count,ESTIMATED_ERR
            WRITE(17,*)
            WRITE(17,'(2F20.16,I10,E16.8)') B1,B2,find_count,ESTIMATED_ERR
            WRITE(17,'(2F20.16,I10,E16.8)') B1,A2,find_count,ESTIMATED_ERR
            WRITE(17,*)
            DO K=0,NDERIV
            DO I=DEG,0,-1
               write(13,*) C0(0:I,DEG-I,K) 
            ENDDO
            ENDDO

         ELSE

            ! decide in which dimension to bisection. 
            ! the easiest and fastest way
            ! DIRECTION=MOD(LEVEL,2)+1
            ! The way leading to more efficient code
            ! using a section that minimizes the maximum error is slightly more
            ! efficient than alternating bisection.
            MIDDLE=(A1+B1)/2
            CALL drvrec_wrap(A1,MIDDLE,A2,B2,DEG,NDERIV,C0,ESTIMATED_ERR1)
            CALL drvrec_wrap(MIDDLE,B1,A2,B2,DEG,NDERIV,C0,ESTIMATED_ERR2)
            MIDDLE=(A2+B2)/2
            CALL drvrec_wrap(A1,B1,A2,MIDDLE,DEG,NDERIV,C0,ESTIMATED_ERR3)
            CALL drvrec_wrap(A1,B1,MIDDLE,B2,DEG,NDERIV,C0,ESTIMATED_ERR4)
            IF (MAX(ESTIMATED_ERR1,ESTIMATED_ERR2)<MAX(ESTIMATED_ERR3,ESTIMATED_ERR4)) THEN
               DIRECTION=1
            ELSE
               DIRECTION=2
            ENDIF

            IF (DIRECTION==1) THEN
                MIDDLE=(A1+B1)/2
                WRITE(6,'(A,E25.18,A)') REPEAT(" ",LEVEL)//"  IF (X1<=",MIDDLE,"_dp) THEN "
                   CALL FIND_DOMAINS(A1,MIDDLE,A2,B2,DEG,NDERIV,ERR,LEVEL+1)
                WRITE(6,'(A)') REPEAT(" ",LEVEL)//"  ELSE"
                   CALL FIND_DOMAINS(MIDDLE,B1,A2,B2,DEG,NDERIV,ERR,LEVEL+1)
                WRITE(6,'(A)') REPEAT(" ",LEVEL)//"  ENDIF"
            ELSE
                MIDDLE=(A2+B2)/2
                WRITE(6,'(A,E25.18,A)') REPEAT(" ",LEVEL)//"  IF (X2<=",MIDDLE,"_dp) THEN "
                   CALL FIND_DOMAINS(A1,B1,A2,MIDDLE,DEG,NDERIV,ERR,LEVEL+1)
                WRITE(6,'(A)') REPEAT(" ",LEVEL)//"  ELSE"
                   CALL FIND_DOMAINS(A1,B1,MIDDLE,B2,DEG,NDERIV,ERR,LEVEL+1)
                WRITE(6,'(A)') REPEAT(" ",LEVEL)//"  ENDIF"
            ENDIF
         ENDIF

      END SUBROUTINE FIND_DOMAINS
!====================================================================================================================
!
!     evaluates C0 coefs increasing the number of digits, till the resulting C0 is converged.
!     This is a simple wrapper to the program drvrec_mpfr.x
!
!====================================================================================================================
      SUBROUTINE drvrec_wrap(A1,B1,A2,B2,DEG,NDERIV,C0,ESTIMATED_ERR)

         REAL(KIND=dp) :: a1,b1,a2,b2,estimated_err,estimated_err2,estimated_err3,esterr
         INTEGER :: deg,nderiv
         REAL(KIND=dp) :: C0(0:DEG,0:DEG,0:NDERIV), C02(0:DEG,0:DEG,0:NDERIV)
         REAL(KIND=dp) :: R,T,vals(0:nderiv)
         integer :: digits_local

         digits_local=digits
         C02 = 0.0_dp

         DO
            OPEN(UNIT=51,FILE="drvrec_mpfr.in")
            WRITE(51,*) digits_local
            WRITE(51,*) deg
            WRITE(51,*) nderiv
            WRITE(51,*) A1,B1
            WRITE(51,*) A2,B2
            WRITE(51,*) functionid
            CLOSE(51)
            CALL SYSTEM("./drvrec_mpfr.x")
            OPEN(UNIT=51,FILE="drvrec_mpfr.out")
            READ(51,*) ESTIMATED_ERR
            READ(51,*) C0
            CLOSE(51)

            !  do not allow for differences larger than ERR times a given factor
            IF(MAXVAL(ABS(C0-C02))>ERR*1.0E-5) THEN
               C02=C0
               digits_local=digits_local*2
            ELSE
               OPEN(UNIT=97,FILE="function_vals.dat")
               DO
                 READ(97,*,END=999) R,T,vals
                 WRITE(98,*) R,T,vals
               ENDDO
999            CONTINUE
               CLOSE(97)
               EXIT
            ENDIF
         ENDDO

      END SUBROUTINE
END MODULE Recurse

!====================================================================================================================
!
!     The main program ... small and easy
!
!     output files....
!     std out          :  a program that can evaluate the function in a given domain
!     history.dat      :  a list of function values, evaluated during the construction of the approximation,
!                         that can be used as a reference data.
!
!====================================================================================================================
PROGRAM Fun2D

 USE Recurse

 OPEN(UNIT=98,FILE="history.dat")
 CALL GO()
 CLOSE(98)

END PROGRAM Fun2D
