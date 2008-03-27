!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of the semi empirical parameter types.
!> \author JGH (14.08.2004)
! *****************************************************************************
MODULE semi_empirical_types

  USE basis_set_types,                 ONLY: deallocate_sto_basis_set,&
                                             sto_basis_set_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndo,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_rm1
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE physcon,                         ONLY: angstrom,&
                                             evolt,&
                                             kcalmol
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_types'

! *****************************************************************************
  TYPE semi_empirical_type
    INTEGER                                :: typ
    CHARACTER(LEN=default_string_length)   :: name
    LOGICAL                                :: defined
    INTEGER                                :: z
    REAL(KIND=dp)                          :: zeff
    INTEGER                                :: natorb
    REAL(KIND = dp), DIMENSION(:), POINTER :: beta
    REAL(KIND = dp), DIMENSION(:), POINTER :: sto_exponents
    TYPE(sto_basis_set_type), POINTER      :: basis
    INTEGER                                :: ngauss
    REAL(KIND = dp)                        :: eheat
    REAL(KIND = dp)                        :: uss,upp,udd,uff
    REAL(KIND = dp)                        :: alp
    REAL(KIND = dp)                        :: eisol
    REAL(KIND = dp)                        :: gss,gsp,gpp,gp2
    REAL(KIND = dp)                        :: gsd,gpd,gdd
    REAL(KIND = dp)                        :: hsp
    REAL(KIND = dp)                        :: dd,qq,am,ad,aq
    REAL(KIND = dp),DIMENSION(2)           :: pre, d
    REAL(KIND = dp),DIMENSION(4)           :: fn1,   fn2,  fn3
    REAL(KIND = dp),DIMENSION(4,4)         :: bfn1, bfn2, bfn3
    ! Specific PM6
    REAL(KIND = dp), DIMENSION(0:115)      :: xab
    REAL(KIND = dp), DIMENSION(0:115)      :: aab
    REAL(KIND = dp)                        :: a, b, c, rho
    REAL(KIND = dp)                        :: f0dd, f2dd, f4dd, f0sd, f0pd, f2pd,&
                                              g1pd, g2sd, g3pd
    REAL(KIND = dp), DIMENSION(:), POINTER :: zn
    REAL(KIND = dp), DIMENSION(9)          :: ko
    REAL(KIND = dp), DIMENSION(6)          :: cs
  END TYPE semi_empirical_type

  PUBLIC :: semi_empirical_type,&
            allocate_se_param,&
            deallocate_se_param,&
            get_se_param,&
            set_se_param,&
            write_se_param

CONTAINS

! *****************************************************************************
  SUBROUTINE allocate_se_param(se_parameter, error)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_se_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL deallocate_se_param(se_parameter, error)

    ALLOCATE (se_parameter,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY ( se_parameter % sto_exponents )
    NULLIFY ( se_parameter % zn  )
    NULLIFY ( se_parameter % basis )

    se_parameter % defined = .FALSE.
    se_parameter % name = ""
    se_parameter % typ = HUGE(0)
    se_parameter % z = -1
    se_parameter % zeff = -1.0_dp
    se_parameter % natorb = 0
    se_parameter % ngauss = 0

    ALLOCATE (se_parameter%beta(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    se_parameter % beta(:) = 0._dp

    ALLOCATE (se_parameter%sto_exponents(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    se_parameter % sto_exponents(:) = -1

    ALLOCATE (se_parameter%zn(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    se_parameter % zn(:) = 0.0_dp

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

    se_parameter % pre = 0._dp
    se_parameter % d = 0._dp

    se_parameter % xab = 0._dp
    se_parameter % aab = 0._dp
    se_parameter % a = 0._dp
    se_parameter % b = 0._dp
    se_parameter % c = 0._dp

  END SUBROUTINE allocate_se_param

! *****************************************************************************
  SUBROUTINE deallocate_se_param(se_parameter, error)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_se_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(se_parameter)) THEN
       CALL deallocate_sto_basis_set(se_parameter%basis, error)
       DEALLOCATE (se_parameter%beta,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (se_parameter%sto_exponents,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (se_parameter%zn,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (se_parameter,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE deallocate_se_param

! *****************************************************************************
  SUBROUTINE get_se_param(se_parameter,name,typ,defined,z,zeff,natorb,eheat,&
       beta,sto_exponents,uss,upp,udd,uff,alp,eisol,&
       gss,gsp,gpp,gp2,hsp,gsd,gpd,gdd,ppddg,dpddg,ngauss)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    INTEGER, INTENT(OUT), OPTIONAL           :: typ
    LOGICAL, INTENT(OUT), OPTIONAL           :: defined
    INTEGER, INTENT(OUT), OPTIONAL           :: z
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: zeff
    INTEGER, INTENT(OUT), OPTIONAL           :: natorb
    REAL(KIND=dp), OPTIONAL                  :: eheat
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: beta, sto_exponents
    REAL(KIND=dp), OPTIONAL                  :: uss, upp, udd, uff, alp, &
                                                eisol, gss, gsp, gpp, gp2, &
                                                hsp, gsd, gpd, gdd
    REAL(KIND=dp), DIMENSION(2), OPTIONAL    :: ppddg, dpddg
    INTEGER, INTENT(OUT), OPTIONAL           :: ngauss

    CHARACTER(len=*), PARAMETER :: routineN = 'get_se_param', &
      routineP = moduleN//':'//routineN

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
       IF (PRESENT(ppddg)) ppddg = se_parameter%pre
       IF (PRESENT(dpddg)) dpddg = se_parameter%d

    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE get_se_param

! *****************************************************************************
  SUBROUTINE set_se_param(se_parameter,name,typ,defined,z,zeff,natorb,eheat,&
       beta,sto_exponents,uss,upp,udd,uff,alp,eisol,&
       gss,gsp,gpp,gp2,hsp,gsd,gpd,gdd,ppddg,dpddg,ngauss)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    INTEGER, INTENT(IN), OPTIONAL            :: typ
    LOGICAL, INTENT(IN), OPTIONAL            :: defined
    INTEGER, INTENT(IN), OPTIONAL            :: z
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: zeff
    INTEGER, INTENT(IN), OPTIONAL            :: natorb
    REAL(KIND=dp), OPTIONAL                  :: eheat
    REAL(dp), DIMENSION(0:), OPTIONAL        :: beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: sto_exponents
    REAL(KIND=dp), OPTIONAL                  :: uss, upp, udd, uff, alp, &
                                                eisol, gss, gsp, gpp, gp2, &
                                                hsp, gsd, gpd, gdd
    REAL(dp), DIMENSION(2), OPTIONAL         :: ppddg, dpddg
    INTEGER, INTENT(IN), OPTIONAL            :: ngauss

    CHARACTER(len=*), PARAMETER :: routineN = 'set_se_param', &
      routineP = moduleN//':'//routineN

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
       IF (PRESENT(ppddg)) se_parameter%pre = ppddg
       IF (PRESENT(dpddg)) se_parameter%d = dpddg

    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer se_parameter is not associated")

    END IF

  END SUBROUTINE set_se_param

! *****************************************************************************
  SUBROUTINE write_se_param(se_parameter,subsys_section,error)

    TYPE(semi_empirical_type), POINTER       :: se_parameter
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_se_param', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: l, lm, natorb, ngauss, &
                                                output_unit, typ, z
    LOGICAL                                  :: defined
    REAL(KIND=dp)                            :: alp, eheat, eisol, gp2, gpp, &
                                                gsp, gss, hsp, udd, uff, upp, &
                                                uss, zeff
    REAL(KIND=dp), DIMENSION(0:3)            :: u
    REAL(KIND=dp), DIMENSION(2)              :: dpddg, ppddg
    REAL(KIND=dp), DIMENSION(:), POINTER     :: beta, sexp
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (ASSOCIATED(se_parameter).AND.BTEST(cp_print_key_should_output(logger%iter_info,subsys_section,&
         "PRINT%KINDS/BASIS_SET",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,subsys_section,"PRINT%KINDS",&
            extension=".Log",error=error)

       IF (output_unit >0) THEN
          CALL get_se_param(se_parameter,name=name,typ=typ,defined=defined,&
               z=z,zeff=zeff,natorb=natorb,eheat=eheat,beta=beta,&
               sto_exponents=sexp,uss=uss,upp=upp,udd=udd,uff=uff,&
               alp=alp,eisol=eisol,&
               gss=gss,gsp=gsp,gpp=gpp,gp2=gp2,hsp=hsp,&
               ppddg=ppddg,dpddg=dpddg,&
               ngauss=ngauss)

          u(0) = uss
          u(1) = upp
          u(2) = udd
          u(3) = uff

          lm=-1

          SELECT CASE (typ)
          CASE DEFAULT
             CALL stop_program(routineN,moduleN,__LINE__,"Semiempirical method unknown")
          CASE (do_method_am1,do_method_rm1)
             SELECT CASE (typ)
             CASE(do_method_am1)
                WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                     " Semi empirical parameters: ","Austin Model 1 (AM1)",TRIM(name)
             CASE(do_method_rm1)
                WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                     " Semi empirical parameters: ","Recife Model 1 (RM1)",TRIM(name)
             END SELECT
             IF ( defined ) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.2)")&
                     "Effective core charge:",zeff
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
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
          CASE (do_method_pm3)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Parametric Method 3 (PM3) ",&
                  TRIM(name)
             IF ( defined ) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.2)")&
                     "Effective core charge:",zeff
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
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
          CASE (do_method_pm6)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Parametric Method 6 (PM6) ",&
                  TRIM(name)
             IF ( defined ) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.2)")&
                     "Effective core charge:",zeff
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
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
          CASE (do_method_pdg)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","PDDG/PM3 ",&
                  TRIM(name)
             IF ( defined ) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.2)")&
                     "Effective core charge:",zeff
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
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
                WRITE (UNIT=output_unit,FMT="(T16,A,T52,F14.10,T67,F14.10)") &
                     "d_PDDG [A^-1]      :",dpddg/angstrom,&
                     "P_PDDG [eV]      :",ppddg*evolt
             ELSE
                WRITE (UNIT=output_unit,FMT="(T55,A)")&
                     "Parameters are not defined"
             END IF
          CASE (do_method_mndo)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","MNDO ",TRIM(name)
             IF ( defined ) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.2)")&
                     "Effective core charge:",zeff
                WRITE (UNIT=output_unit,FMT="(T16,A,T71,I10)")&
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
       CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
            "PRINT%KINDS",error=error)
    END IF

  END SUBROUTINE write_se_param

END MODULE semi_empirical_types
