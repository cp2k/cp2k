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

  USE basis_set_types,                 ONLY: allocate_sto_basis_set,&
                                             create_gto_from_sto_basis,&
                                             deallocate_sto_basis_set,&
                                             gto_basis_set_type,&
                                             set_sto_basis_set
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: angstrom,&
                                             evolt,&
                                             kcalmol
  USE semi_empirical_parameters,       ONLY: am1_default_parameter,&
                                             eht_default_parameter,&
                                             mndo_default_parameter,&
                                             pm3_default_parameter
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "semi_empirical_utils"

  PUBLIC :: allocate_se_param,&
            deallocate_se_param,&
            get_se_param,&
            get_se_type,&
            init_se_param,&
            read_se_param,&
            set_se_param,&
            write_se_param,&
            se_param_set_default,&
            expand,&
            contract

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_se_param(se_parameter)

    TYPE(semi_empirical_type), POINTER       :: se_parameter

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

    se_parameter % uss = 0._dp
    se_parameter % upp = 0._dp
    se_parameter % udd = 0._dp
    se_parameter % uff = 0._dp
    se_parameter % alp = 0._dp
    se_parameter % eisol = 0._dp

    se_parameter % gss = 0._dp
    se_parameter % gsp = 0._dp
    se_parameter % gpp = 0._dp
    se_parameter % gp2 = 0._dp
    se_parameter % hsp = 0._dp

    se_parameter % gsd = 0._dp
    se_parameter % gpd = 0._dp
    se_parameter % gdd = 0._dp

  END SUBROUTINE allocate_se_param

! *****************************************************************************

  SUBROUTINE deallocate_se_param(se_parameter)

    TYPE(semi_empirical_type), POINTER       :: se_parameter

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

  SUBROUTINE get_se_param(se_parameter,name,typ,defined,z,zeff,natorb,eheat,&
                          beta,sto_exponents,uss,upp,udd,uff,alp,eisol,&
                          gss,gsp,gpp,gp2,hsp,gsd,gpd,gdd,ngauss)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name, typ
    LOGICAL, INTENT(OUT), OPTIONAL           :: defined
    INTEGER, INTENT(OUT), OPTIONAL           :: z, zeff, natorb
    REAL(KIND=dp), OPTIONAL                  :: eheat
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: beta, sto_exponents
    REAL(KIND=dp), OPTIONAL                  :: uss, upp, udd, uff, alp, &
                                                eisol, gss, gsp, gpp, gp2, &
                                                hsp, gsd, gpd, gdd
    INTEGER, INTENT(OUT), OPTIONAL           :: ngauss

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "get_se_param"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) THEN

      IF (PRESENT(name)) name = se_parameter%name
      IF (PRESENT(typ)) typ = se_parameter%typ
      IF (PRESENT(defined)) defined = se_parameter%defined
      IF (PRESENT(z)) z = se_parameter%z
      IF (PRESENT(zeff)) zeff = se_parameter%zeff
      IF (PRESENT(natorb)) natorb = se_parameter%natorb
      IF (PRESENT(eheat)) eheat = se_parameter%eheat
      IF (PRESENT(beta)) beta => se_parameter%beta
      IF (PRESENT(sto_exponents)) sto_exponents => se_parameter%sto_exponents
      IF (PRESENT(ngauss)) ngauss = se_parameter%ngauss
      IF (PRESENT(uss)) uss = se_parameter%uss
      IF (PRESENT(upp)) upp = se_parameter%upp
      IF (PRESENT(udd)) udd = se_parameter%udd
      IF (PRESENT(uff)) uff = se_parameter%uff
      IF (PRESENT(alp)) alp = se_parameter%alp
      IF (PRESENT(eisol)) eisol = se_parameter%eisol
      IF (PRESENT(gss)) gss = se_parameter%gss
      IF (PRESENT(gsp)) gsp = se_parameter%gsp
      IF (PRESENT(gpp)) gpp = se_parameter%gpp
      IF (PRESENT(gp2)) gp2 = se_parameter%gp2
      IF (PRESENT(hsp)) hsp = se_parameter%hsp
      IF (PRESENT(gsd)) gsd = se_parameter%gsd
      IF (PRESENT(gpd)) gpd = se_parameter%gpd
      IF (PRESENT(gdd)) gdd = se_parameter%gdd

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE get_se_param

! *****************************************************************************

  SUBROUTINE init_se_param(se_parameter,orb_basis_set)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "init_se_param"

    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: symbol
    INTEGER                                  :: dl, istat, l, nl, nshell
    INTEGER, DIMENSION(:), POINTER           :: lq, nq
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet

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

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "read_se_param"

!   ---------------------------------------------------------------------------

  END SUBROUTINE read_se_param

! *****************************************************************************

  SUBROUTINE set_se_param(se_parameter,name,typ,defined,z,zeff,natorb,eheat,&
                          beta,sto_exponents,uss,upp,udd,uff,alp,eisol,&
                          gss,gsp,gpp,gp2,hsp,gsd,gpd,gdd,ngauss)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name, typ
    LOGICAL, INTENT(IN), OPTIONAL            :: defined
    INTEGER, INTENT(IN), OPTIONAL            :: z, zeff, natorb
    REAL(KIND=dp), OPTIONAL                  :: eheat
    REAL(dp), DIMENSION(0:), OPTIONAL        :: beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: sto_exponents
    REAL(KIND=dp), OPTIONAL                  :: uss, upp, udd, uff, alp, &
                                                eisol, gss, gsp, gpp, gp2, &
                                                hsp, gsd, gpd, gdd
    INTEGER, INTENT(IN), OPTIONAL            :: ngauss

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "set_se_param"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(se_parameter)) THEN

      IF (PRESENT(name)) se_parameter%name = name
      IF (PRESENT(typ)) se_parameter%typ = typ
      IF (PRESENT(defined)) se_parameter%defined = defined
      IF (PRESENT(z)) se_parameter%z = z
      IF (PRESENT(zeff)) se_parameter%zeff = zeff
      IF (PRESENT(natorb)) se_parameter%natorb = natorb
      IF (PRESENT(eheat)) se_parameter%eheat = eheat
      IF (PRESENT(beta)) se_parameter%beta = beta
      IF (PRESENT(sto_exponents)) se_parameter%sto_exponents = sto_exponents
      IF (PRESENT(ngauss)) se_parameter%ngauss = ngauss
      IF (PRESENT(uss)) se_parameter%uss = uss
      IF (PRESENT(upp)) se_parameter%upp = upp
      IF (PRESENT(udd)) se_parameter%udd = udd
      IF (PRESENT(uff)) se_parameter%uff = uff
      IF (PRESENT(alp)) se_parameter%alp = alp
      IF (PRESENT(eisol)) se_parameter%eisol = eisol
      IF (PRESENT(gss)) se_parameter%gss = gss
      IF (PRESENT(gsp)) se_parameter%gsp = gsp
      IF (PRESENT(gpp)) se_parameter%gpp = gpp
      IF (PRESENT(gp2)) se_parameter%gp2 = gp2
      IF (PRESENT(hsp)) se_parameter%hsp = hsp
      IF (PRESENT(gsd)) se_parameter%gsd = gsd
      IF (PRESENT(gpd)) se_parameter%gpd = gpd
      IF (PRESENT(gdd)) se_parameter%gdd = gdd

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE set_se_param

! *****************************************************************************

  SUBROUTINE write_se_param(se_parameter,globenv)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_se_param"

    CHARACTER(LEN=default_string_length)     :: name, typ
    INTEGER                                  :: l, lm, natorb, ngauss, &
                                                output_unit, z, zeff
    LOGICAL                                  :: defined
    REAL(KIND=dp)                            :: alp, eheat, eisol, gp2, gpp, &
                                                gsp, gss, hsp, udd, uff, upp, &
                                                uss
    REAL(KIND=dp), DIMENSION(0:3)            :: u
    REAL(KIND=dp), DIMENSION(:), POINTER     :: beta, sexp

!   ---------------------------------------------------------------------------
!   REAL(KIND = dp)                        :: dd,qq,am,ad,aq
!   REAL(KIND = dp),DIMENSION(4)           :: fn1,fn2,fn3

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(se_parameter).AND.globenv%print%basis_set) THEN
      CALL get_se_param(se_parameter,name=name,typ=typ,defined=defined,&
                        z=z,zeff=zeff,natorb=natorb,eheat=eheat,beta=beta,&
                        sto_exponents=sexp,uss=uss,upp=upp,udd=udd,uff=uff,&
                        alp=alp,eisol=eisol,&
                        gss=gss,gsp=gsp,gpp=gpp,gp2=gp2,hsp=hsp,&
                        ngauss=ngauss)
      output_unit = globenv%scr

      u(0) = uss
      u(1) = upp
      u(2) = udd
      u(3) = uff

      lm=-1

      SELECT CASE (typ)
        CASE DEFAULT
          CALL stop_program(routine_name,module_name,__LINE__,&
                          "Typ unknown: "//typ)
        CASE ('EHT')
          WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
            " Semi empirical parameters: ","Extended Huckel ",TRIM(name)
          IF ( defined ) THEN
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
              "Effective core charge:",zeff, &
              "Number of orbitals:",natorb, &
              "Basis set expansion (STO-NG)",ngauss
            WRITE (UNIT=output_unit,FMT="(T16,A,T66,F15.5)")&
              "Atomic heat of formation [kcal/mol]:",eheat*kcalmol
            DO l=0,3
              IF (ABS(beta(l)) > 0._dp) THEN
                WRITE (UNIT=output_unit,&
                 FMT="(T16,A,I2,T48,A,F17.4,/,T47,A,F17.4)") &
                 "Shell: ",l,"Slater exponent:",sexp(l),&
                 "H parameter [eV]:",beta(l)*evolt
              END IF
            END DO
          ELSE
            WRITE (UNIT=output_unit,FMT="(T55,A)")&
               "Parameters are not defined"
          END IF
        CASE ('AM1')
          WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
            " Semi empirical parameters: ","Austin Model 1 (AM1)",TRIM(name)
          IF ( defined ) THEN
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
              "Effective core charge:",zeff, &
              "Number of orbitals:",natorb, &
              "Basis set expansion (STO-NG)",ngauss
            WRITE (UNIT=output_unit,FMT="(T16,A,T66,F15.5)")&
              "Atomic heat of formation [kcal/mol]:",eheat*kcalmol
            DO l=0,3
              IF (ABS(beta(l)) > 0._dp) THEN
                WRITE (UNIT=output_unit,&
                 FMT="(T16,A,I2,T48,A,F17.4,/,T22,A,F12.4,T52,A,F12.4)") &
                 "Shell: ",l,"Slater exponent:",sexp(l),&
                 "U [eV]    :",u(l)*evolt,"Beta [eV]       :",beta(l)*evolt
                lm=l
              END IF
            END DO
            WRITE (UNIT=output_unit,FMT="(T16,A,F12.4,T52,A,F12.4)") &
               "alp [A^-1]      :",alp/angstrom,&
               "eisol [eV]      :",eisol*evolt
            WRITE (UNIT=output_unit,FMT="(T16,A,/,T16,5(A,F8.3))") &
               "One-centre-two-electron integrals [eV]",&
               " gss=",gss*evolt," gsp=",gsp*evolt,&
               " gpp=",gpp*evolt," gp2=",gp2*evolt," hsp=",hsp*evolt
          ELSE
            WRITE (UNIT=output_unit,FMT="(T55,A)")&
               "Parameters are not defined"
          END IF
        CASE ('PM3')
          WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
            " Semi empirical parameters: ","Parametric Method 3 (PM3) ",&
            TRIM(name)
          IF ( defined ) THEN
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
              "Effective core charge:",zeff, &
              "Number of orbitals:",natorb, &
              "Basis set expansion (STO-NG)",ngauss
            WRITE (UNIT=output_unit,FMT="(T16,A,T66,F15.5)")&
              "Atomic heat of formation [kcal/mol]:",eheat*kcalmol
            DO l=0,3
              IF (ABS(beta(l)) > 0._dp) THEN
                WRITE (UNIT=output_unit,&
                 FMT="(T16,A,I2,T48,A,F17.4,/,T22,A,F12.4,T52,A,F12.4)") &
                 "Shell: ",l,"Slater exponent:",sexp(l),&
                 "U [eV]    :",u(l)*evolt,"Beta [eV]       :",beta(l)*evolt
                lm=l
              END IF
            END DO
            WRITE (UNIT=output_unit,FMT="(T16,A,F12.4,T52,A,F12.4)") &
               "alp [A^-1]      :",alp/angstrom,&
               "eisol [eV]      :",eisol*evolt
            WRITE (UNIT=output_unit,FMT="(T16,A,/,T16,5(A,F8.3))") &
               "One-centre-two-electron integrals [eV]",&
               " gss=",gss*evolt," gsp=",gsp*evolt,&
               " gpp=",gpp*evolt," gp2=",gp2*evolt," hsp=",hsp*evolt
          ELSE
            WRITE (UNIT=output_unit,FMT="(T55,A)")&
               "Parameters are not defined"
          END IF
        CASE ('MNDO')
          WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
            " Semi empirical parameters: ","MNDO ",TRIM(name)
          IF ( defined ) THEN
            WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
              "Effective core charge:",zeff, &
              "Number of orbitals:",natorb, &
              "Basis set expansion (STO-NG)",ngauss
            WRITE (UNIT=output_unit,FMT="(T16,A,T66,F15.5)")&
              "Atomic heat of formation [kcal/mol]:",eheat*kcalmol
            DO l=0,3
              IF (ABS(beta(l)) > 0._dp) THEN
                WRITE (UNIT=output_unit,&
                 FMT="(T16,A,I2,T48,A,F17.4,/,T22,A,F12.4,T52,A,F12.4)") &
                 "Shell: ",l,"Slater exponent:",sexp(l),&
                 "U [eV]    :",u(l)*evolt,"Beta [eV]       :",beta(l)*evolt
                lm=l
              END IF
            END DO
            WRITE (UNIT=output_unit,FMT="(T16,A,F12.4,T52,A,F12.4)") &
               "alp [A^-1]      :",alp/angstrom,&
               "eisol [eV]      :",eisol*evolt
            WRITE (UNIT=output_unit,FMT="(T16,A,/,T16,5(A,F8.3))") &
               "One-centre-two-electron integrals [eV]",&
               " gss=",gss*evolt," gsp=",gsp*evolt,&
               " gpp=",gpp*evolt," gp2=",gp2*evolt," hsp=",hsp*evolt
          ELSE
            WRITE (UNIT=output_unit,FMT="(T55,A)")&
               "Parameters are not defined"
          END IF
      END SELECT

    END IF

  END SUBROUTINE write_se_param

! *****************************************************************************

  SUBROUTINE se_param_set_default(se_parameter,z,method)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    INTEGER, INTENT(IN)                      :: z
    CHARACTER(LEN=*)                         :: method

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
          CALL am1_default_parameter(se_parameter,z)
        CASE ("PM3")
          CALL pm3_default_parameter(se_parameter,z)
        CASE ("MNDO")
          CALL mndo_default_parameter(se_parameter,z)
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

  FUNCTION get_se_type(se_method) RESULT(se_type)

    CHARACTER(LEN=*)                         :: se_method
    INTEGER                                  :: se_type

    SELECT CASE (se_method)
       CASE DEFAULT
         se_type=0
       CASE ("EHT")
         se_type=0
       CASE ("AM1")
         se_type=2
       CASE ("PM3")
         se_type=3
       CASE ("MNDO")
         se_type=1
    END SELECT 

  END FUNCTION get_se_type

! *****************************************************************************

  SUBROUTINE contract ( amat, avec )
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: amat
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: avec

    avec( 1) = amat(1,1)
    avec( 2) = amat(4,1) + amat(1,4)
    avec( 3) = amat(4,4)
    avec( 4) = amat(2,1) + amat(1,2)
    avec( 5) = amat(2,4) + amat(4,2)
    avec( 6) = amat(2,2)
    avec( 7) = amat(3,1) + amat(1,3)
    avec( 8) = amat(3,4) + amat(4,3)
    avec( 9) = amat(3,2) + amat(2,3)
    avec(10) = amat(3,3)

  END SUBROUTINE contract

  SUBROUTINE expand ( amat, avec )
    REAL(dp), DIMENSION(:, :), INTENT(INOUT)   :: amat
    REAL(dp), DIMENSION(:), INTENT(IN)       :: avec

    amat(1,1) = amat(1,1) + avec( 1)
    amat(4,1) = amat(4,1) + avec( 2)
    amat(1,4) = amat(1,4) + avec( 2)
    amat(4,4) = amat(4,4) + avec( 3)
    amat(2,1) = amat(2,1) + avec( 4)
    amat(1,2) = amat(1,2) + avec( 4)
    amat(2,4) = amat(2,4) + avec( 5)
    amat(4,2) = amat(4,2) + avec( 5)
    amat(2,2) = amat(2,2) + avec( 6)
    amat(3,1) = amat(3,1) + avec( 7)
    amat(1,3) = amat(1,3) + avec( 7)
    amat(3,4) = amat(3,4) + avec( 8)
    amat(4,3) = amat(4,3) + avec( 8)
    amat(3,2) = amat(3,2) + avec( 9)
    amat(2,3) = amat(2,3) + avec( 9)
    amat(3,3) = amat(3,3) + avec(10)

  END SUBROUTINE expand

! *****************************************************************************

END MODULE semi_empirical_utils

! *****************************************************************************
