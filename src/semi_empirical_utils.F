!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 CP2K developers group                                  !
!-----------------------------------------------------------------------------!
!!****** cp2k/semi_empirical_utils [1.0] *
!!
!!   NAME
!!     semi_empirical_utils
!!
!!   FUNCTION
!!     Working with the  semi empirical parameter types.
!!
!!   AUTHOR
!!     JGH (14.08.2004)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE semi_empirical_utils

  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE basis_set_types,                 ONLY: gto_basis_set_type,&
                                             sto_basis_set_type,&
                                             allocate_sto_basis_set,&
                                             deallocate_sto_basis_set,&
                                             set_sto_basis_set,&
                                             create_gto_from_sto_basis
  USE global_types,                    ONLY: global_environment_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: evolt
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_parameters,       ONLY: eht_default_parameter
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "semi_empirical_utils"

  PUBLIC :: allocate_se_param,&
            deallocate_se_param,&
            get_se_param,&
            init_se_param,&
            read_se_param,&
            set_se_param,&
            write_se_param,&
            se_param_set_default

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_se_param(se_parameter)

    TYPE(semi_empirical_type), POINTER        :: se_parameter

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_se_param"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) CALL deallocate_se_param(se_parameter)

    ALLOCATE (se_parameter,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "se_parameter",0)

    NULLIFY ( se_parameter % sto_exponents )
    NULLIFY ( se_parameter % basis )

    se_parameter % defined = .FALSE.
    se_parameter % name = ""
    se_parameter % typ = "NONE"
    se_parameter % z = -1
    se_parameter % zeff = -1
    se_parameter % natorb = 0
    se_parameter % ngauss = 0

    ALLOCATE (se_parameter%beta(0:3),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "se_parameter%beta",0)
    se_parameter % beta(:) = 0._dp

    ALLOCATE (se_parameter%sto_exponents(0:3),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "se_parameter%states",0)
    se_parameter % sto_exponents(:) = -1

  END SUBROUTINE allocate_se_param

! *****************************************************************************

  SUBROUTINE deallocate_se_param(se_parameter)

    TYPE(semi_empirical_type), POINTER        :: se_parameter

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_se_param"
    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) THEN
      IF (ASSOCIATED(se_parameter%basis)) THEN
        CALL deallocate_sto_basis_set(se_parameter%basis)
      END IF
      DEALLOCATE (se_parameter%beta,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "se_parameter%beta")
      DEALLOCATE (se_parameter%sto_exponents,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "se_parameter%sto_exponents")
      DEALLOCATE (se_parameter,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "se_parameter")
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer se_parameter is not associated")
    END IF

  END SUBROUTINE deallocate_se_param

! *****************************************************************************

  SUBROUTINE get_se_param(se_parameter,name,typ,defined,z,zeff,natorb,&
                          beta,sto_exponents,ngauss)

    TYPE(semi_empirical_type), POINTER     :: se_parameter
    CHARACTER(LEN=default_string_length), &
                   INTENT(OUT), OPTIONAL   :: name
    CHARACTER(LEN=default_string_length), &
                   INTENT(OUT), OPTIONAL   :: typ
    LOGICAL, INTENT(OUT), OPTIONAL         :: defined
    INTEGER, INTENT(OUT), OPTIONAL         :: z,zeff
    INTEGER, INTENT(OUT), OPTIONAL         :: natorb
    REAL(KIND = dp), DIMENSION(:), POINTER, &
                   OPTIONAL                :: beta
    REAL(KIND = dp), DIMENSION(:), POINTER, &
                   OPTIONAL                :: sto_exponents
    INTEGER, INTENT(OUT), OPTIONAL         :: ngauss

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_se_param"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) THEN

      IF (PRESENT(name)) name = se_parameter%name
      IF (PRESENT(typ)) typ = se_parameter%typ
      IF (PRESENT(defined)) defined = se_parameter%defined
      IF (PRESENT(z)) z = se_parameter%z
      IF (PRESENT(zeff)) zeff = se_parameter%zeff
      IF (PRESENT(natorb)) natorb = se_parameter%natorb
      IF (PRESENT(beta)) beta => se_parameter%beta
      IF (PRESENT(sto_exponents)) sto_exponents => se_parameter%sto_exponents
      IF (PRESENT(ngauss)) ngauss = se_parameter%ngauss

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE get_se_param

! *****************************************************************************

  SUBROUTINE init_se_param(se_parameter,orb_basis_set)

    TYPE(semi_empirical_type), POINTER     :: se_parameter
    TYPE(gto_basis_set_type), POINTER      :: orb_basis_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_se_param"

    INTEGER                                :: nshell, l, istat, nl, dl
    CHARACTER(LEN=6), DIMENSION(:),POINTER :: symbol
    INTEGER, DIMENSION(:), POINTER         :: nq, lq
    REAL(KIND=dp), DIMENSION(:), POINTER   :: zet

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) THEN

      CALL allocate_sto_basis_set (se_parameter%basis)
      nshell = 0
      IF (se_parameter%natorb ==  1) nshell=1
      IF (se_parameter%natorb ==  4) nshell=2
      IF (se_parameter%natorb ==  9) nshell=3
      IF (se_parameter%natorb == 16) nshell=4
      ALLOCATE (nq(0:3),lq(0:3),zet(0:3),STAT=istat)
      ALLOCATE (symbol(0:3),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "nq, lq, zet, symbol",0)
      nq=0
      lq=0
      zet=0._dp
      DO l=0,nshell-1
        dl=2*l+1
        nl=ptable(se_parameter%z)%e_conv(l)
        IF ( nl == 0 ) THEN
          nq(l) = l+1
        ELSE
          nq(l)=(nl-1+2*dl)/(2*dl) + l
        END IF
        lq(l)=l
        zet(l)=se_parameter%sto_exponents(l)
        IF(l==0) WRITE(symbol(0),'(I1,A1)') nq(l),"S" 
        IF(l==1) WRITE(symbol(1),'(I1,A1)') nq(l),"P" 
        IF(l==2) WRITE(symbol(2),'(I1,A1)') nq(l),"D" 
        IF(l==3) WRITE(symbol(3),'(I1,A1)') nq(l),"F" 
      END DO

      IF (nshell > 0) THEN
        CALL set_sto_basis_set(se_parameter%basis,name=se_parameter%name,&
                               nshell=nshell,symbol=symbol,&
                               nq=nq,lq=lq,zet=zet)
        CALL create_gto_from_sto_basis(se_parameter%basis,orb_basis_set,&
                               se_parameter%ngauss)
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE init_se_param

! *****************************************************************************

  SUBROUTINE read_se_param()

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_se_param"

!   ---------------------------------------------------------------------------

  END SUBROUTINE read_se_param

! *****************************************************************************

  SUBROUTINE set_se_param(se_parameter,name,typ,defined,z,zeff,natorb,&
                          beta,sto_exponents,ngauss)

    TYPE(semi_empirical_type), POINTER     :: se_parameter
    CHARACTER(LEN=default_string_length), &
                   INTENT(IN), OPTIONAL   :: name
    CHARACTER(LEN=default_string_length), &
                   INTENT(IN), OPTIONAL   :: typ
    LOGICAL, INTENT(IN), OPTIONAL         :: defined
    INTEGER, INTENT(IN), OPTIONAL         :: z,zeff
    INTEGER, INTENT(IN), OPTIONAL         :: natorb
    INTEGER, DIMENSION(:), OPTIONAL       :: beta
    REAL(KIND = dp), DIMENSION(:), &
                   OPTIONAL               :: sto_exponents
    INTEGER, INTENT(IN), OPTIONAL         :: ngauss

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_se_param"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) THEN

      IF (PRESENT(name)) se_parameter%name = name
      IF (PRESENT(typ)) se_parameter%typ = typ
      IF (PRESENT(defined)) se_parameter%defined = defined
      IF (PRESENT(z)) se_parameter%z = z
      IF (PRESENT(zeff)) se_parameter%zeff = zeff
      IF (PRESENT(natorb)) se_parameter%natorb = natorb
      IF (PRESENT(beta)) se_parameter%beta = beta
      IF (PRESENT(sto_exponents)) se_parameter%sto_exponents = sto_exponents
      IF (PRESENT(ngauss)) se_parameter%ngauss = ngauss

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE set_se_param

! *****************************************************************************

  SUBROUTINE write_se_param(se_parameter,globenv)

    TYPE(semi_empirical_type), POINTER     :: se_parameter
    TYPE(global_environment_type), &
      INTENT(IN)                           :: globenv

    CHARACTER(LEN=default_string_length)   :: name
    CHARACTER(LEN=default_string_length)   :: typ
    LOGICAL                                :: defined
    INTEGER                                :: z,zeff
    INTEGER                                :: natorb
    REAL(KIND = dp), DIMENSION(:), POINTER :: beta
    REAL(KIND = dp), DIMENSION(:), POINTER :: sexp
    INTEGER                                :: ngauss
    INTEGER                                :: output_unit, l

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_se_param"

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(se_parameter).AND.globenv%print%basis_set) THEN
      CALL get_se_param(se_parameter,name=name,typ=typ,defined=defined,&
                        z=z,zeff=zeff,natorb=natorb,beta=beta,&
                        sto_exponents=sexp,ngauss=ngauss)
      output_unit = globenv%scr

      SELECT CASE (typ)
      CASE DEFAULT
        CALL stop_program(routine_name,module_name,__LINE__,&
                        "Typ unknown: "//typ)
      CASE ('EHT')
        WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
          " Semi empirical parameters: ","Extended Huckel ",&
          ADJUSTR(TRIM(name))
        IF ( defined ) THEN
          WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
            "Effective core charge:",zeff, &
            "Number of orbitals:",natorb, &
            "Basis set expansion (STO-NG)",ngauss
          DO l=0,3
            IF (ABS(beta(l)) > 0._dp) THEN
              WRITE (UNIT=output_unit,FMT="(T16,A,I2,T48,A,F17.4,/,T47,A,F17.4)") &
               "Shell: ",l,"Slater exponent:",sexp(l),&
               "H parameter [eV]:",beta(l)*evolt
            END IF
          END DO
        ELSE
          WRITE (UNIT=output_unit,FMT="(T55,A)")&
             "Parameters are not defined"
        END IF
      END SELECT

    END IF

  END SUBROUTINE write_se_param

! *****************************************************************************

  SUBROUTINE se_param_set_default(se_parameter,z,method)

    TYPE(semi_empirical_type), POINTER     :: se_parameter
    INTEGER, INTENT(IN)                    :: z
    CHARACTER(LEN=*)                       :: method

    CHARACTER(LEN=*), PARAMETER :: routine_name = "se_param_set_default"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) THEN
      IF ( z < 0 ) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                        "Atomic number < 0")
      END IF
      SELECT CASE (method)
        CASE ("EHT")
          CALL eht_default_parameter(se_parameter,z)
        CASE ("AM1")
        CASE ("PM3")
        CASE ("MNDO")
        CASE DEFAULT
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "Invalid method:"//method)
      END SELECT
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE se_param_set_default

! *****************************************************************************

END MODULE semi_empirical_utils
