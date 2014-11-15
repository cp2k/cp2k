! *****************************************************************************
MODULE fparser
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! Fortran 90 function parser v1.1
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  !
  ! This public domain function parser module is intended for applications
  ! where a set of mathematical expressions is specified at runtime and is
  ! then evaluated for a large number of variable values. This is done by
  ! compiling the set of function strings into byte code, which is interpreted
  ! very efficiently for the various variable values.
  !
  ! The source code is available from:
  ! http://www.its.uni-karlsruhe.de/~schmehl/opensource/fparser-v1.1.tar.gz
  !
  ! Please send comments, corrections or questions to the author:
  ! Roland Schmehl <Roland.Schmehl@mach.uni-karlsruhe.de>
  !
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             rn => dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fparser'

  !------- -------- --------- --------- --------- --------- --------- --------- -------
  PUBLIC :: initf,&  ! Initialize function parser for n functions
            parsef,&  ! Parse single function string
            evalf,&  ! Evaluate single function
            finalizef,&  ! Finalize the function parser
            evalfd
  INTEGER, PUBLIC            :: EvalErrType ! =0: no error occured, >0: evaluation error
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  PRIVATE
  SAVE
  INTEGER, PARAMETER, PRIVATE :: is = SELECTED_INT_KIND(1) !Data type of bytecode
  INTEGER(is),                              PARAMETER :: cImmed   = 1,          &
                                                         cNeg     = 2,          &
                                                         cAdd     = 3,          &
                                                         cSub     = 4,          &
                                                         cMul     = 5,          &
                                                         cDiv     = 6,          &
                                                         cPow     = 7,          &
                                                         cAbs     = 8,          &
                                                         cExp     = 9,          &
                                                         cLog10   = 10,         &
                                                         cLog     = 11,         &
                                                         cSqrt    = 12,         &
                                                         cSinh    = 13,         &
                                                         cCosh    = 14,         &
                                                         cTanh    = 15,         &
                                                         cSin     = 16,         &
                                                         cCos     = 17,         &
                                                         cTan     = 18,         &
                                                         cAsin    = 19,         &
                                                         cAcos    = 20,         &
                                                         cAtan    = 21,         &
                                                         VarBegin = 22
  CHARACTER (LEN=1), DIMENSION(cAdd:cPow),  PARAMETER :: Ops      = (/ '+',     &
                                                                       '-',     &
                                                                       '*',     &
                                                                       '/',     &
                                                                       '^' /)
  CHARACTER (LEN=5), DIMENSION(cAbs:cAtan), PARAMETER :: Funcs    = (/ 'abs  ', &
                                                                       'exp  ', &
                                                                       'log10', &
                                                                       'log  ', &
                                                                       'sqrt ', &
                                                                       'sinh ', &
                                                                       'cosh ', &
                                                                       'tanh ', &
                                                                       'sin  ', &
                                                                       'cos  ', &
                                                                       'tan  ', &
                                                                       'asin ', &
                                                                       'acos ', &
                                                                       'atan ' /)
! *****************************************************************************
  TYPE tComp
     INTEGER(is), DIMENSION(:), POINTER :: ByteCode
     INTEGER                            :: ByteCodeSize
     REAL(rn),    DIMENSION(:), POINTER :: Immed
     INTEGER                            :: ImmedSize
     REAL(rn),    DIMENSION(:), POINTER :: Stack
     INTEGER                            :: StackSize, &
                                           StackPtr
  END TYPE tComp
  TYPE(tComp),   DIMENSION(:),  POINTER :: Comp  ! Bytecode
  INTEGER,   DIMENSION(:),  ALLOCATABLE :: ipos  ! Associates function strings
  !
CONTAINS
  !
! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE finalizef()
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Finalize function parser
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(len=*), PARAMETER :: routineN = 'finalizef', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat

!----- -------- --------- --------- --------- --------- --------- --------- -------

    DO i = 1, SIZE(Comp)
       IF (ASSOCIATED(Comp(i)%ByteCode)) THEN
          DEALLOCATE ( Comp(i)%ByteCode, stat = istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"Comp(i)%ByteCode")
       END IF
       IF (ASSOCIATED(Comp(i)%Immed)) THEN
          DEALLOCATE ( Comp(i)%Immed, stat=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"Comp(i)%Immed")
       END IF
       IF (ASSOCIATED(Comp(i)%Stack)) THEN
          DEALLOCATE ( Comp(i)%Stack, stat=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"Comp(i)%Stack")
       END IF
    END DO

    DEALLOCATE ( Comp, stat=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"Comp")

  END SUBROUTINE finalizef
  !
! *****************************************************************************
!> \brief ...
!> \param n ...
! *****************************************************************************
  SUBROUTINE initf (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Initialize function parser for n functions
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: n

    INTEGER                                  :: i

! Number of functions
!----- -------- --------- --------- --------- --------- --------- --------- -------

    ALLOCATE (Comp(n))
    DO i=1,n
       NULLIFY (Comp(i)%ByteCode,Comp(i)%Immed,Comp(i)%Stack)
    END DO
  END SUBROUTINE initf
  !
! *****************************************************************************
!> \brief ...
!> \param i ...
!> \param FuncStr ...
!> \param Var ...
! *****************************************************************************
  SUBROUTINE parsef (i, FuncStr, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Parse ith function string FuncStr and compile it into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: i
    CHARACTER(LEN=*), INTENT(in)             :: FuncStr
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: Var

    CHARACTER(len=*), PARAMETER :: routineN = 'parsef', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=LEN(FuncStr))              :: Func

! Function identifier
! Function string
! Array with variable names
! Function string, local use
!----- -------- --------- --------- --------- --------- --------- --------- -------

    IF (i < 1 .OR. i > SIZE(Comp)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "Function number is out of range")
    END IF
    ALLOCATE (ipos(LEN_TRIM(FuncStr)))                       ! Char. positions in orig. string
    Func = FuncStr                                           ! Local copy of function string
    CALL Replace ('**','^ ',Func)                            ! Exponent into 1-Char. format
    CALL RemoveSpaces (Func)                                 ! Condense function string
    CALL CheckSyntax (Func,FuncStr,Var)
    DEALLOCATE (ipos)
    CALL Compile (i,Func,Var)                                ! Compile into bytecode

  END SUBROUTINE parsef
  !
! *****************************************************************************
!> \brief ...
!> \param i ...
!> \param Val ...
!> \retval res ...
! *****************************************************************************
  FUNCTION evalf (i, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function for the values passed in array Val(:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: i
    REAL(rn), DIMENSION(:), INTENT(in)       :: Val
    REAL(rn)                                 :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'evalf', &
      routineP = moduleN//':'//routineN
    REAL(rn), PARAMETER                      :: zero = 0._rn

    INTEGER                                  :: DP, IP, ipow, SP

! Function identifier
! Variable values
! Result
! Instruction pointer
! Data pointer
! Stack pointer
!----- -------- --------- --------- --------- --------- --------- --------- -------

    DP = 1
    SP = 0
    DO IP=1,Comp(i)%ByteCodeSize
       SELECT CASE (Comp(i)%ByteCode(IP))

       CASE (cImmed); SP=SP+1; Comp(i)%Stack(SP)=Comp(i)%Immed(DP); DP=DP+1
       CASE   (cNeg); Comp(i)%Stack(SP)=-Comp(i)%Stack(SP)
       CASE   (cAdd); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)+Comp(i)%Stack(SP); SP=SP-1
       CASE   (cSub); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)-Comp(i)%Stack(SP); SP=SP-1
       CASE   (cMul); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)*Comp(i)%Stack(SP); SP=SP-1
       CASE   (cDiv); IF (Comp(i)%Stack(SP)==0._rn) THEN; EvalErrType=1; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)/Comp(i)%Stack(SP); SP=SP-1
       CASE   (cPow)
          ! Fixing for possible Negative floating-point value raised to a real power
          IF (Comp(i)%Stack(SP-1)<0.0_rn) THEN
             ipow = FLOOR(Comp(i)%Stack(SP))
             IF (MOD(Comp(i)%Stack(SP),REAL(ipow,KIND=rn))==0.0_rn) THEN
                Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)**ipow
             ELSE
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "Negative floating-point value raised to a real power!")
             END IF
          ELSE
             Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)**Comp(i)%Stack(SP)
          END IF
          SP=SP-1
       CASE   (cAbs); Comp(i)%Stack(SP)=ABS(Comp(i)%Stack(SP))
       CASE   (cExp); Comp(i)%Stack(SP)=EXP(Comp(i)%Stack(SP))
       CASE (cLog10); IF (Comp(i)%Stack(SP)<=0._rn) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=LOG10(Comp(i)%Stack(SP))
       CASE   (cLog); IF (Comp(i)%Stack(SP)<=0._rn) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=LOG(Comp(i)%Stack(SP))
       CASE  (cSqrt); IF (Comp(i)%Stack(SP)<0._rn) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=SQRT(Comp(i)%Stack(SP))
       CASE  (cSinh); Comp(i)%Stack(SP)=SINH(Comp(i)%Stack(SP))
       CASE  (cCosh); Comp(i)%Stack(SP)=COSH(Comp(i)%Stack(SP))
       CASE  (cTanh); Comp(i)%Stack(SP)=TANH(Comp(i)%Stack(SP))
       CASE   (cSin); Comp(i)%Stack(SP)=SIN(Comp(i)%Stack(SP))
       CASE   (cCos); Comp(i)%Stack(SP)=COS(Comp(i)%Stack(SP))
       CASE   (cTan); Comp(i)%Stack(SP)=TAN(Comp(i)%Stack(SP))
       CASE  (cAsin); IF ((Comp(i)%Stack(SP)<-1._rn).OR.(Comp(i)%Stack(SP)>1._rn)) THEN
                      EvalErrType=4; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=ASIN(Comp(i)%Stack(SP))
       CASE  (cAcos); IF ((Comp(i)%Stack(SP)<-1._rn).OR.(Comp(i)%Stack(SP)>1._rn)) THEN
                      EvalErrType=4; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=ACOS(Comp(i)%Stack(SP))
       CASE  (cAtan); Comp(i)%Stack(SP)=ATAN(Comp(i)%Stack(SP))
       CASE  DEFAULT; SP=SP+1; Comp(i)%Stack(SP)=Val(Comp(i)%ByteCode(IP)-VarBegin+1)
       END SELECT
    END DO
    EvalErrType = 0
    res = Comp(i)%Stack(1)
  END FUNCTION evalf
  !
! *****************************************************************************
!> \brief ...
!> \param Func ...
!> \param FuncStr ...
!> \param Var ...
! *****************************************************************************
  SUBROUTINE CheckSyntax (Func,FuncStr,Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check syntax of function string,  returns 0 if syntax is ok
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(in)             :: Func, FuncStr
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: Var

    INTEGER                                  :: ib, in, j, lFunc, ParCnt
    CHARACTER(LEN=1)                         :: c
    INTEGER(is)                              :: n
    LOGICAL                                  :: err
    REAL(rn)                                 :: r

! Function string without spaces
! Original function string
! Array with variable names
! Parenthesis counter
!----- -------- --------- --------- --------- --------- --------- --------- -------

    j = 1
    ParCnt = 0
    lFunc = LEN_TRIM(Func)
    step: DO
       IF (j > lFunc) CALL ParseErrMsg (j, FuncStr)
       c = Func(j:j)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Check for valid operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (c == '-' .OR. c == '+') THEN                      ! Check for leading - or +
          j = j+1
          IF (j > lFunc) CALL ParseErrMsg (j, FuncStr, 'Missing operand')
          c = Func(j:j)
          IF (ANY(c == Ops)) CALL ParseErrMsg (j, FuncStr, 'Multiple operators')
       END IF
       n = MathFunctionIndex (Func(j:))
       IF (n > 0) THEN                                       ! Check for math function
          j = j+LEN_TRIM(Funcs(n))
          IF (j > lFunc) CALL ParseErrMsg (j, FuncStr, 'Missing function argument')
          c = Func(j:j)
          IF (c /= '(') CALL ParseErrMsg (j, FuncStr, 'Missing opening parenthesis')
       END IF
       IF (c == '(') THEN                                    ! Check for opening parenthesis
          ParCnt = ParCnt+1
          j = j+1
          CYCLE step
       END IF
       IF (SCAN(c,'0123456789.') > 0) THEN                   ! Check for number
          r = RealNum (Func(j:),ib,in,err)
          IF (err) CALL ParseErrMsg (j, FuncStr, 'Invalid number format:  '//Func(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       ELSE                                                  ! Check for variable
          n = VariableIndex (Func(j:),Var,ib,in)
          IF (n == 0) CALL ParseErrMsg (j, FuncStr, 'Invalid element: '//Func(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       END IF
       DO WHILE (c == ')')                                   ! Check for closing parenthesis
          ParCnt = ParCnt-1
          IF (ParCnt < 0) CALL ParseErrMsg (j, FuncStr, 'Mismatched parenthesis')
          IF (Func(j-1:j-1) == '(') CALL ParseErrMsg (j-1, FuncStr, 'Empty parentheses')
          j = j+1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       END DO
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have a legal operand: A legal operator or end of string must follow
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (j > lFunc) EXIT
       IF (ANY(c == Ops)) THEN                               ! Check for multiple operators
          IF (j+1 > lFunc) CALL ParseErrMsg (j, FuncStr)
          IF (ANY(Func(j+1:j+1) == Ops)) CALL ParseErrMsg (j+1, FuncStr, 'Multiple operators')
       ELSE                                                  ! Check for next operand
          CALL ParseErrMsg (j, FuncStr, 'Missing operator')
       END IF
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have an operand and an operator: the next loop will check for another
       ! operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       j = j+1
    END DO step
    IF (ParCnt > 0) CALL ParseErrMsg (j, FuncStr, 'Missing )')
  END SUBROUTINE CheckSyntax
  !
! *****************************************************************************
!> \brief ...
!> \retval msg ...
! *****************************************************************************
  FUNCTION EvalErrMsg () RESULT (msg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return error message
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), DIMENSION(4), PARAMETER :: m = (/ &
      'Division by zero                ', 'Argument of SQRT negative       ', &
      'Argument of LOG negative        ', 'Argument of ASIN or ACOS illegal' &
      /)
    CHARACTER(LEN=LEN(m))                    :: msg

!----- -------- --------- --------- --------- --------- --------- --------- -------

    IF (EvalErrType < 1 .OR. EvalErrType > SIZE(m)) THEN
       msg = ''
    ELSE
       msg = m(EvalErrType)
    ENDIF
  END FUNCTION EvalErrMsg
  !
! *****************************************************************************
!> \brief ...
!> \param j ...
!> \param FuncStr ...
!> \param Msg ...
! *****************************************************************************
  SUBROUTINE ParseErrMsg (j, FuncStr, Msg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Print error message and terminate program
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: j
    CHARACTER(LEN=*), INTENT(in)             :: FuncStr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL   :: Msg

    CHARACTER(len=*), PARAMETER :: routineN = 'ParseErrMsg', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: message

! Original function string
!----- -------- --------- --------- --------- --------- --------- --------- -------

    IF (PRESENT(Msg)) THEN
       WRITE (UNIT=message,FMT="(A)") "Syntax error in function string: "//Msg
    ELSE
       WRITE (UNIT=message,FMT="(A)") "Syntax error in function string"
    END IF
    WRITE (*,'(/,T2,A)') TRIM(FuncStr)
    IF ((j > LBOUND(ipos,DIM=1)).AND.(j <= UBOUND(ipos,DIM=1))) THEN
       WRITE (*,'(A)') REPEAT(" ",ipos(j))//"?"
    ELSE
       WRITE (*,'(A)') REPEAT(" ",SIZE(ipos)+1)//"?"
    END IF
    CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))

  END SUBROUTINE ParseErrMsg
  !
! *****************************************************************************
!> \brief ...
!> \param c ...
!> \retval n ...
! *****************************************************************************
  FUNCTION OperatorIndex (c) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return operator index
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=1), INTENT(in)             :: c
    INTEGER(is)                              :: n

    INTEGER(is)                              :: j

!----- -------- --------- --------- --------- --------- --------- --------- -------

    n = 0
    DO j=cAdd,cPow
       IF (c == Ops(j)) THEN
          n = j
          EXIT
       END IF
    END DO
  END FUNCTION OperatorIndex
  !
! *****************************************************************************
!> \brief ...
!> \param str ...
!> \retval n ...
! *****************************************************************************
  FUNCTION MathFunctionIndex (str) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return index of math function beginnig at 1st position of string str
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(in)             :: str
    INTEGER(is)                              :: n

    CHARACTER(LEN=LEN(Funcs))                :: fun
    INTEGER                                  :: k
    INTEGER(is)                              :: j

!----- -------- --------- --------- --------- --------- --------- --------- -------

    n = 0
    DO j=cAbs,cAtan                                          ! Check all math functions
       k = MIN(LEN_TRIM(Funcs(j)), LEN(str))
       CALL LowCase (str(1:k), fun)
       IF (fun == Funcs(j)) THEN                             ! Compare lower case letters
          n = j                                              ! Found a matching function
          EXIT
       END IF
    END DO
  END FUNCTION MathFunctionIndex
  !
! *****************************************************************************
!> \brief ...
!> \param str ...
!> \param Var ...
!> \param ibegin ...
!> \param inext ...
!> \retval n ...
! *****************************************************************************
  FUNCTION VariableIndex (str, Var, ibegin, inext) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return index of variable at begin of string str (returns 0 if no variable found)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(in)             :: str
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: Var
    INTEGER, INTENT(out), OPTIONAL           :: ibegin, inext
    INTEGER(is)                              :: n

    INTEGER                                  :: ib, in, j, lstr

! String
! Array with variable names
! Index of variable
! Start position of variable name
! Position of character after name
!----- -------- --------- --------- --------- --------- --------- --------- -------

    n = 0
    lstr = LEN_TRIM(str)
    IF (lstr > 0) THEN
       DO ib=1,lstr                                          ! Search for first character in str
          IF (str(ib:ib) /= ' ') EXIT                        ! When lstr>0 at least 1 char in str
       END DO
       DO in=ib,lstr                                         ! Search for name terminators
          IF (SCAN(str(in:in),'+-*/^) ') > 0) EXIT
       END DO
       DO j=1,SIZE(Var)
          IF (str(ib:in-1) == Var(j)) THEN
             n = j                                           ! Variable name found
             EXIT
          END IF
       END DO
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
  END FUNCTION VariableIndex
  !
! *****************************************************************************
!> \brief ...
!> \param str ...
! *****************************************************************************
  SUBROUTINE RemoveSpaces (str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Remove Spaces from string, remember positions of characters in old string
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(inout)          :: str

    INTEGER                                  :: k, lstr

!----- -------- --------- --------- --------- --------- --------- --------- -------

    lstr = LEN_TRIM(str)
    ipos(:) = (/ (k,k=1,lstr) /)
    k = 1
    DO WHILE (str(k:lstr) /= ' ')
       IF (str(k:k) == ' ') THEN
          str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
          ipos(k:lstr) = (/ ipos(k+1:lstr), 0 /)             ! Move 1 element to left
          k = k-1
       END IF
       k = k+1
    END DO
  END SUBROUTINE RemoveSpaces
  !
! *****************************************************************************
!> \brief ...
!> \param ca ...
!> \param cb ...
!> \param str ...
! *****************************************************************************
  SUBROUTINE Replace (ca,cb,str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Replace ALL appearances of character set ca in string str by character set cb
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(in)             :: ca
    CHARACTER(LEN=LEN(ca)), INTENT(in)       :: cb
    CHARACTER(LEN=*), INTENT(inout)          :: str

    INTEGER                                  :: j, lca

! LEN(ca) must be LEN(cb)
!----- -------- --------- --------- --------- --------- --------- --------- -------

    lca = LEN(ca)
    DO j=1,LEN_TRIM(str)-lca+1
       IF (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    END DO
  END SUBROUTINE Replace
  !
! *****************************************************************************
!> \brief ...
!> \param i ...
!> \param F ...
!> \param Var ...
! *****************************************************************************
  SUBROUTINE Compile (i, F, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Compile i-th function string F into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: i
    CHARACTER(LEN=*), INTENT(in)             :: F
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: Var

    CHARACTER(len=*), PARAMETER :: routineN = 'Compile', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

! Function identifier
! Function string
! Array with variable names
!----- -------- --------- --------- --------- --------- --------- --------- -------

    IF (ASSOCIATED(Comp(i)%ByteCode)) DEALLOCATE ( Comp(i)%ByteCode, &
                                                   Comp(i)%Immed,    &
                                                   Comp(i)%Stack     )
    Comp(i)%ByteCodeSize = 0
    Comp(i)%ImmedSize    = 0
    Comp(i)%StackSize    = 0
    Comp(i)%StackPtr     = 0
    CALL CompileSubstr (i,F,1,LEN_TRIM(F),Var)               ! Compile string to determine size
    ALLOCATE ( Comp(i)%ByteCode(Comp(i)%ByteCodeSize), &
               Comp(i)%Immed(Comp(i)%ImmedSize),       &
               Comp(i)%Stack(Comp(i)%StackSize),       &
               STAT = istat                            )
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"Byte code")
    ELSE
       Comp(i)%ByteCodeSize = 0
       Comp(i)%ImmedSize    = 0
       Comp(i)%StackSize    = 0
       Comp(i)%StackPtr     = 0
       CALL CompileSubstr (i,F,1,LEN_TRIM(F),Var)            ! Compile string into bytecode
    END IF
    !
  END SUBROUTINE Compile
  !
! *****************************************************************************
!> \brief ...
!> \param i ...
!> \param b ...
! *****************************************************************************
  SUBROUTINE AddCompiledByte (i, b)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Add compiled byte to bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: i
    INTEGER(is), INTENT(in)                  :: b

! Function identifier
! Value of byte to be added
!----- -------- --------- --------- --------- --------- --------- --------- -------

    Comp(i)%ByteCodeSize = Comp(i)%ByteCodeSize + 1
    IF (ASSOCIATED(Comp(i)%ByteCode)) Comp(i)%ByteCode(Comp(i)%ByteCodeSize) = b
  END SUBROUTINE AddCompiledByte
  !
! *****************************************************************************
!> \brief ...
!> \param i ...
!> \param F ...
!> \param Var ...
!> \retval n ...
! *****************************************************************************
  FUNCTION MathItemIndex (i, F, Var) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return math item index, if item is real number, enter it into Comp-structure
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: i
    CHARACTER(LEN=*), INTENT(in)             :: F
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: Var
    INTEGER(is)                              :: n

! Function identifier
! Function substring
! Array with variable names
! Byte value of math item
!----- -------- --------- --------- --------- --------- --------- --------- -------

    n = 0
    IF (SCAN(F(1:1),'0123456789.') > 0) THEN                 ! Check for begin of a number
       Comp(i)%ImmedSize = Comp(i)%ImmedSize + 1
       IF (ASSOCIATED(Comp(i)%Immed)) Comp(i)%Immed(Comp(i)%ImmedSize) = RealNum (F)
       n = cImmed
    ELSE                                                     ! Check for a variable
       n = VariableIndex (F, Var)
       IF (n > 0) n = VarBegin+n-1_is
    END IF
  END FUNCTION MathItemIndex
  !
! *****************************************************************************
!> \brief ...
!> \param F ...
!> \param b ...
!> \param e ...
!> \retval res ...
! *****************************************************************************
  FUNCTION CompletelyEnclosed (F, b, e) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(in)             :: F
    INTEGER, INTENT(in)                      :: b, e
    LOGICAL                                  :: res

    INTEGER                                  :: j, k

! Function substring
! First and last pos. of substring
!----- -------- --------- --------- --------- --------- --------- --------- -------

    res=.FALSE.
    IF (F(b:b) == '(' .AND. F(e:e) == ')') THEN
       k = 0
       DO j=b+1,e-1
          IF     (F(j:j) == '(') THEN
             k = k+1
          ELSEIF (F(j:j) == ')') THEN
             k = k-1
          END IF
          IF (k < 0) EXIT
       END DO
       IF (k == 0) res=.TRUE.                                ! All opened parenthesis closed
    END IF
  END FUNCTION CompletelyEnclosed
  !
! *****************************************************************************
!> \brief ...
!> \param i ...
!> \param F ...
!> \param b ...
!> \param e ...
!> \param Var ...
! *****************************************************************************
  RECURSIVE SUBROUTINE CompileSubstr (i, F, b, e, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Compile i-th function string F into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: i
    CHARACTER(LEN=*), INTENT(in)             :: F
    INTEGER, INTENT(in)                      :: b, e
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in)                             :: Var

    CHARACTER(LEN=*), PARAMETER :: &
      calpha = 'abcdefghijklmnopqrstuvwxyz'// 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    INTEGER                                  :: b2, io, j, k
    INTEGER(is)                              :: n

! Function identifier
! Function substring
! Begin and end position substring
! Array with variable names
!----- -------- --------- --------- --------- --------- --------- --------- -------
! Check for special cases of substring
!----- -------- --------- --------- --------- --------- --------- --------- -------

    IF     (F(b:b) == '+') THEN                              ! Case 1: F(b:e) = '+...'
!      WRITE(*,*)'1. F(b:e) = "+..."'
       CALL CompileSubstr (i, F, b+1, e, Var)
       RETURN
    ELSEIF (CompletelyEnclosed (F, b, e)) THEN               ! Case 2: F(b:e) = '(...)'
!      WRITE(*,*)'2. F(b:e) = "(...)"'
       CALL CompileSubstr (i, F, b+1, e-1, Var)
       RETURN
    ELSEIF (SCAN(F(b:b),calpha) > 0) THEN
       n = MathFunctionIndex (F(b:e))
       IF (n > 0) THEN
          b2 = b+INDEX(F(b:e),'(')-1
          IF (CompletelyEnclosed(F, b2, e)) THEN             ! Case 3: F(b:e) = 'fcn(...)'
!            WRITE(*,*)'3. F(b:e) = "fcn(...)"'
             CALL CompileSubstr(i, F, b2+1, e-1, Var)
             CALL AddCompiledByte (i, n)
             RETURN
          END IF
       END IF
    ELSEIF (F(b:b) == '-') THEN
       IF (CompletelyEnclosed (F, b+1, e)) THEN              ! Case 4: F(b:e) = '-(...)'
!         WRITE(*,*)'4. F(b:e) = "-(...)"'
          CALL CompileSubstr (i, F, b+2, e-1, Var)
          CALL AddCompiledByte (i, cNeg)
          RETURN
       ELSEIF (SCAN(F(b+1:b+1),calpha) > 0) THEN
          n = MathFunctionIndex (F(b+1:e))
          IF (n > 0) THEN
             b2 = b+INDEX(F(b+1:e),'(')
             IF (CompletelyEnclosed(F, b2, e)) THEN          ! Case 5: F(b:e) = '-fcn(...)'
!               WRITE(*,*)'5. F(b:e) = "-fcn(...)"'
                CALL CompileSubstr(i, F, b2+1, e-1, Var)
                CALL AddCompiledByte (i, n)
                CALL AddCompiledByte (i, cNeg)
                RETURN
             END IF
          END IF
       ENDIF
    END IF
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    DO io=cAdd,cPow                                          ! Increasing priority +-*/^
       k = 0
       DO j=e,b,-1
          IF     (F(j:j) == ')') THEN
             k = k+1
          ELSEIF (F(j:j) == '(') THEN
             k = k-1
          END IF
          IF (k == 0 .AND. F(j:j) == Ops(io) .AND. IsBinaryOp (j, F)) THEN
             IF (ANY(F(j:j) == Ops(cMul:cPow)) .AND. F(b:b) == '-') THEN ! Case 6: F(b:e) = '-...Op...' with Op > -
!               WRITE(*,*)'6. F(b:e) = "-...Op..." with Op > -'
                CALL CompileSubstr (i, F, b+1, e, Var)
                CALL AddCompiledByte (i, cNeg)
                RETURN
             ELSE                                                        ! Case 7: F(b:e) = '...BinOp...'
!               WRITE(*,*)'7. Binary operator',F(j:j)
                CALL CompileSubstr (i, F, b, j-1, Var)
                CALL CompileSubstr (i, F, j+1, e, Var)
                CALL AddCompiledByte (i, OperatorIndex(Ops(io)))
                Comp(i)%StackPtr = Comp(i)%StackPtr - 1
                RETURN
             END IF
          END IF
       END DO
    END DO
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for remaining items, i.e. variables or explicit numbers
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    b2 = b
    IF (F(b:b) == '-') b2 = b2+1
    n = MathItemIndex(i, F(b2:e), Var)
!   WRITE(*,*)'8. AddCompiledByte ',n
    CALL AddCompiledByte (i, n)
    Comp(i)%StackPtr = Comp(i)%StackPtr + 1
    IF (Comp(i)%StackPtr > Comp(i)%StackSize) Comp(i)%StackSize = Comp(i)%StackSize + 1
    IF (b2 > b) CALL AddCompiledByte (i, cNeg)
  END SUBROUTINE CompileSubstr
  !
! *****************************************************************************
!> \brief ...
!> \param j ...
!> \param F ...
!> \retval res ...
! *****************************************************************************
  FUNCTION IsBinaryOp (j, F) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check if operator F(j:j) in string F is binary operator
    ! Special cases already covered elsewhere:              (that is corrected in v1.1)
    ! - operator character F(j:j) is first character of string (j=1)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    INTEGER, INTENT(in)                      :: j
    CHARACTER(LEN=*), INTENT(in)             :: F
    LOGICAL                                  :: res

    INTEGER                                  :: k
    LOGICAL                                  :: Dflag, Pflag

! Position of Operator
! String
! Result
!----- -------- --------- --------- --------- --------- --------- --------- -------

    res=.TRUE.
    IF (F(j:j) == '+' .OR. F(j:j) == '-') THEN               ! Plus or minus sign:
       IF (j == 1) THEN                                      ! - leading unary operator ?
          res = .FALSE.
       ELSEIF (SCAN(F(j-1:j-1),'+-*/^(') > 0) THEN           ! - other unary operator ?
          res = .FALSE.
       ELSEIF (SCAN(F(j+1:j+1),'0123456789') > 0 .AND. &     ! - in exponent of real number ?
               SCAN(F(j-1:j-1),'eEdD')       > 0) THEN
          Dflag=.FALSE.; Pflag=.FALSE.
          k = j-1
          DO WHILE (k > 1)                                   !   step to the left in mantissa
             k = k-1
             IF     (SCAN(F(k:k),'0123456789') > 0) THEN
                Dflag=.TRUE.
             ELSEIF (F(k:k) == '.') THEN
                IF (Pflag) THEN
                   EXIT                                      !   * EXIT: 2nd appearance of '.'
                ELSE
                   Pflag=.TRUE.                              !   * mark 1st appearance of '.'
                ENDIF
             ELSE
                EXIT                                         !   * all other characters
             END IF
          END DO
          IF (Dflag .AND. (k == 1 .OR. SCAN(F(k:k),'+-*/^(') > 0)) res = .FALSE.
       END IF
    END IF
  END FUNCTION IsBinaryOp
  !
! *****************************************************************************
!> \brief ...
!> \param str ...
!> \param ibegin ...
!> \param inext ...
!> \param error ...
!> \retval res ...
! *****************************************************************************
  FUNCTION RealNum (str, ibegin, inext, error) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(in)             :: str
    INTEGER, INTENT(out), OPTIONAL           :: ibegin, inext
    LOGICAL, INTENT(out), OPTIONAL           :: error
    REAL(rn)                                 :: res

    INTEGER                                  :: ib, in, istat
    LOGICAL                                  :: Bflag, DInExp, DInMan, Eflag, &
                                                err, InExp, InMan, Pflag

! String
! Real number
! Start position of real number
! 1st character after real number
! Error flag
! .T. at begin of number in str
! .T. in mantissa of number
! .T. after 1st '.' encountered
! .T. at exponent identifier 'eEdD'
! .T. in exponent of number
! .T. if at least 1 digit in mant.
! .T. if at least 1 digit in exp.
! Local error flag
!----- -------- --------- --------- --------- --------- --------- --------- -------

    Bflag=.TRUE.; InMan=.FALSE.; Pflag=.FALSE.; Eflag=.FALSE.; InExp=.FALSE.
    DInMan=.FALSE.; DInExp=.FALSE.
    ib   = 1
    in   = 1
    DO WHILE (in <= LEN_TRIM(str))
       SELECT CASE (str(in:in))
       CASE (' ')                                            ! Only leading blanks permitted
          ib = ib+1
          IF (InMan .OR. Eflag .OR. InExp) EXIT
       CASE ('+','-')                                        ! Permitted only
          IF     (Bflag) THEN
             InMan=.TRUE.; Bflag=.FALSE.                     ! - at beginning of mantissa
          ELSEIF (Eflag) THEN
             InExp=.TRUE.; Eflag=.FALSE.                     ! - at beginning of exponent
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE ('0':'9')                                        ! Mark
          IF     (Bflag) THEN
             InMan=.TRUE.; Bflag=.FALSE.                     ! - beginning of mantissa
          ELSEIF (Eflag) THEN
             InExp=.TRUE.; Eflag=.FALSE.                     ! - beginning of exponent
          ENDIF
          IF (InMan) DInMan=.TRUE.                           ! Mantissa contains digit
          IF (InExp) DInExp=.TRUE.                           ! Exponent contains digit
       CASE ('.')
          IF     (Bflag) THEN
             Pflag=.TRUE.                                    ! - mark 1st appearance of '.'
             InMan=.TRUE.; Bflag=.FALSE.                     !   mark beginning of mantissa
          ELSEIF (InMan .AND..NOT.Pflag) THEN
             Pflag=.TRUE.                                    ! - mark 1st appearance of '.'
          ELSE
             EXIT                                            ! - otherwise STOP
          END IF
       CASE ('e','E','d','D')                                ! Permitted only
          IF (InMan) THEN
             Eflag=.TRUE.; InMan=.FALSE.                     ! - following mantissa
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE DEFAULT
          EXIT                                               ! STOP at all other characters
       END SELECT
       in = in+1
    END DO
    err = (ib > in-1) .OR. (.NOT.DInMan) .OR. ((Eflag.OR.InExp).AND..NOT.DInExp)
    IF (err) THEN
       res = 0.0_rn
    ELSE
       READ(str(ib:in-1),*,IOSTAT=istat) res
       err = istat /= 0
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
    IF (PRESENT(error))  error  = err
  END FUNCTION RealNum
  !
! *****************************************************************************
!> \brief ...
!> \param str1 ...
!> \param str2 ...
! *****************************************************************************
  SUBROUTINE LowCase (str1, str2)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Transform upper case letters in str1 into lower case letters, result is str2
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    CHARACTER(LEN=*), INTENT(in)             :: str1
    CHARACTER(LEN=*), INTENT(out)            :: str2

    CHARACTER(LEN=*), PARAMETER :: lc = 'abcdefghijklmnopqrstuvwxyz', &
      uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    INTEGER                                  :: j, k

!----- -------- --------- --------- --------- --------- --------- --------- -------

    str2 = str1
    DO j=1,LEN_TRIM(str1)
       k = INDEX(uc,str1(j:j))
       IF (k > 0) str2(j:j) = lc(k:k)
    END DO
  END SUBROUTINE LowCase

! *****************************************************************************
!> \brief Evaluates derivatives
!> \param id_fun ...
!> \param ipar ...
!> \param vals ...
!> \param h ...
!> \param err ...
!> \retval derivative ...
!> \author Main algorithm from Numerical Recipes
!>      Ridders, C.J.F. 1982 - Advances in Engineering Software, Vol.4, no. 2, pp. 75-76
! *****************************************************************************
  FUNCTION evalfd(id_fun,ipar,vals,h,err) RESULT(derivative)
    INTEGER, INTENT(IN)                      :: id_fun, ipar
    REAL(KIND=rn), DIMENSION(:), &
      INTENT(INOUT)                          :: vals
    REAL(KIND=rn), INTENT(IN)                :: h
    REAL(KIND=rn), INTENT(OUT)               :: err
    REAL(KIND=rn)                            :: derivative

    CHARACTER(len=*), PARAMETER :: routineN = 'evalfd', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: ntab = 10
    REAL(KIND=rn), PARAMETER                 :: big_error = 1.0E30_rn, &
                                                con = 1.4_rn, con2 = con*con, &
                                                safe = 2.0_rn

    INTEGER                                  :: i, j
    REAL(KIND=rn)                            :: a(ntab,ntab), errt, fac, &
                                                funcm, funcp, hh, xval

    derivative = HUGE(0.0_rn)
    IF(h/=0._rn) THEN
       xval = vals(ipar)
       hh=h
       vals(ipar) = xval + hh
       funcp = evalf(id_fun, vals)
       vals(ipar) = xval - hh
       funcm = evalf(id_fun, vals)
       a(1,1)=(funcp-funcm)/(2.0_rn*hh)
       err=big_error
       DO i=2,ntab
          hh=hh/con
          vals(ipar) = xval + hh
          funcp = evalf(id_fun, vals)
          vals(ipar) = xval - hh
          funcm = evalf(id_fun, vals)
          a(1,i)=(funcp-funcm)/(2.0_rn*hh)
          fac=con2
          DO j=2,i
             a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.0_rn)
             fac=con2*fac
             errt=MAX(ABS(a(j,i)-a(j-1,i)),ABS(a(j,i)-a(j-1,i-1)))
             IF (errt.LE.err) THEN
                err=errt
                derivative=a(j,i)
             ENDIF
          END DO
          IF(ABS(a(i,i)-a(i-1,i-1)).GE.safe*err)RETURN
       END DO
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "DX provided equals zero!")
    END IF
    vals(ipar)=xval
  END FUNCTION evalfd

END MODULE fparser
