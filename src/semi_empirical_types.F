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
                                             do_method_mndod,&
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
!> \brief Semi-empirical type
! *****************************************************************************
  TYPE semi_empirical_type
    INTEGER                                :: typ
    CHARACTER(LEN=default_string_length)   :: name
    LOGICAL                                :: defined, dorb
    INTEGER                                :: z
    REAL(KIND=dp)                          :: zeff
    INTEGER                                :: natorb
    REAL(KIND = dp), DIMENSION(:), POINTER :: beta
    REAL(KIND = dp), DIMENSION(:), POINTER :: sto_exponents
    REAL(KIND = dp), DIMENSION(:), POINTER :: zn
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
    REAL(KIND = dp)                        :: f0dd, f2dd, f4dd, f0sd, f0pd, f2pd,&
                                              g1pd, g2sd, g3pd
    REAL(KIND = dp), DIMENSION(9)          :: ko
    REAL(KIND = dp), DIMENSION(6)          :: cs
    REAL(KIND = dp), DIMENSION(52)         :: onec2el
    ! Specific for PM6
    REAL(KIND = dp), DIMENSION(0:115)      :: xab
    REAL(KIND = dp), DIMENSION(0:115)      :: aab
    REAL(KIND = dp)                        :: a, b, c, rho
  END TYPE semi_empirical_type

  PUBLIC :: semi_empirical_type,&
            allocate_se_param,&
            deallocate_se_param,&
            get_se_param,&
            set_se_param,&
            write_se_param

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical type
! *****************************************************************************
  SUBROUTINE allocate_se_param(sep, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_se_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL deallocate_se_param(sep, error)

    ALLOCATE (sep,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sep%beta(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sep%sto_exponents(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sep%zn(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(sep%basis)
    CALL zero_se_param(sep,error)

  END SUBROUTINE allocate_se_param

! *****************************************************************************
!> \brief Zero the whole semi-empirical type
! ***************************************************************************** 
  SUBROUTINE zero_se_param(sep, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'zero_se_param', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(sep),cp_failure_level,routineP,error,failure)
    sep%defined = .FALSE.
    sep%dorb    = .FALSE.
    sep%name    = ""
    sep%typ     = HUGE(0)
    sep%z       = HUGE(0)
    sep%zeff    = HUGE(0.0_dp)
    sep%natorb  = 0
    sep%ngauss  = 0
    sep%eheat   = HUGE(0.0_dp)

    sep%zn            = 0.0_dp
    sep%sto_exponents = 0.0_dp
    sep%beta          = 0.0_dp

    sep%uss    = 0.0_dp !eV
    sep%upp    = 0.0_dp !eV
    sep%udd    = 0.0_dp !eV
    sep%uff    = 0.0_dp
    sep%alp    = 0.0_dp
    sep%eisol  = 0.0_dp
    sep%gss    = 0.0_dp
    sep%gsp    = 0.0_dp
    sep%gpp    = 0.0_dp
    sep%gp2    = 0.0_dp
    sep%gsd    = 0.0_dp
    sep%gpd    = 0.0_dp
    sep%gdd    = 0.0_dp
    sep%hsp    = 0.0_dp
    sep%dd     = 0.0_dp
    sep%qq     = 0.0_dp
    sep%am     = 0.0_dp
    sep%ad     = 0.0_dp
    sep%aq     = 0.0_dp

    sep%fn1    = 0.0_dp
    sep%fn2    = 0.0_dp
    sep%fn3    = 0.0_dp
    sep%bfn1   = 0.0_dp
    sep%bfn2   = 0.0_dp
    sep%bfn3   = 0.0_dp

    sep%pre    = 0.0_dp
    sep%d      = 0.0_dp

    sep%xab    = 0.0_dp
    sep%aab    = 0.0_dp
    sep%a      = 0.0_dp
    sep%b      = 0.0_dp
    sep%c      = 0.0_dp
    sep%rho    = 0.0_dp

    sep%f0dd   = 0.0_dp
    sep%f2dd   = 0.0_dp
    sep%f4dd   = 0.0_dp
    sep%f0sd   = 0.0_dp
    sep%f0pd   = 0.0_dp
    sep%f2pd   = 0.0_dp
    sep%g1pd   = 0.0_dp
    sep%g2sd   = 0.0_dp
    sep%g3pd   = 0.0_dp
    sep%ko     = 0.0_dp
    sep%cs     = 0.0_dp
    sep%onec2el= 0.0_dp

  END SUBROUTINE zero_se_param

! *****************************************************************************
!> \brief Deallocate the semi-empirical type
! ***************************************************************************** 
  SUBROUTINE deallocate_se_param(sep, error)

    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_se_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(sep)) THEN
       CALL deallocate_sto_basis_set(sep%basis, error)
       DEALLOCATE (sep%beta,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (sep%sto_exponents,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (sep%zn,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (sep,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE deallocate_se_param

! *****************************************************************************
!> \brief Get info from the semi-empirical type
! ***************************************************************************** 
  SUBROUTINE get_se_param(sep,name,typ,defined,z,zeff,natorb,eheat,&
       beta,sto_exponents,uss,upp,udd,uff,alp,eisol,gss,gsp,gpp,gp2,&
       hsp,gsd,gpd,gdd,ppddg,dpddg,ngauss)

    TYPE(semi_empirical_type), POINTER       :: sep
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

    IF (ASSOCIATED(sep)) THEN
       IF (PRESENT(name)) name = sep%name
       IF (PRESENT(typ)) typ = sep%typ
       IF (PRESENT(defined)) defined = sep%defined
       IF (PRESENT(z)) z = sep%z
       IF (PRESENT(zeff)) zeff = sep%zeff
       IF (PRESENT(natorb)) natorb = sep%natorb
       IF (PRESENT(eheat)) eheat = sep%eheat
       IF (PRESENT(beta)) beta => sep%beta
       IF (PRESENT(sto_exponents)) sto_exponents => sep%sto_exponents
       IF (PRESENT(ngauss)) ngauss = sep%ngauss
       IF (PRESENT(uss)) uss = sep%uss
       IF (PRESENT(upp)) upp = sep%upp
       IF (PRESENT(udd)) udd = sep%udd
       IF (PRESENT(uff)) uff = sep%uff
       IF (PRESENT(alp)) alp = sep%alp
       IF (PRESENT(eisol)) eisol = sep%eisol
       IF (PRESENT(gss)) gss = sep%gss
       IF (PRESENT(gsp)) gsp = sep%gsp
       IF (PRESENT(gpp)) gpp = sep%gpp
       IF (PRESENT(gp2)) gp2 = sep%gp2
       IF (PRESENT(hsp)) hsp = sep%hsp
       IF (PRESENT(gsd)) gsd = sep%gsd
       IF (PRESENT(gpd)) gpd = sep%gpd
       IF (PRESENT(gdd)) gdd = sep%gdd
       IF (PRESENT(ppddg)) ppddg = sep%pre
       IF (PRESENT(dpddg)) dpddg = sep%d
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sep is not associated")
    END IF

  END SUBROUTINE get_se_param

! *****************************************************************************
!> \brief Set info from the semi-empirical type
! ***************************************************************************** 
  SUBROUTINE set_se_param(sep,name,typ,defined,z,zeff,natorb,eheat,&
       beta,sto_exponents,uss,upp,udd,uff,alp,eisol,gss,gsp,gpp,gp2,&
       hsp,gsd,gpd,gdd,ppddg,dpddg,ngauss)

    TYPE(semi_empirical_type), POINTER       :: sep
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

    IF (ASSOCIATED(sep)) THEN
       IF (PRESENT(name)) sep%name = name
       IF (PRESENT(typ)) sep%typ = typ
       IF (PRESENT(defined)) sep%defined = defined
       IF (PRESENT(z)) sep%z = z
       IF (PRESENT(zeff)) sep%zeff = zeff
       IF (PRESENT(natorb)) sep%natorb = natorb
       IF (PRESENT(eheat)) sep%eheat = eheat
       IF (PRESENT(beta)) sep%beta = beta
       IF (PRESENT(sto_exponents)) sep%sto_exponents = sto_exponents
       IF (PRESENT(ngauss)) sep%ngauss = ngauss
       IF (PRESENT(uss)) sep%uss = uss
       IF (PRESENT(upp)) sep%upp = upp
       IF (PRESENT(udd)) sep%udd = udd
       IF (PRESENT(uff)) sep%uff = uff
       IF (PRESENT(alp)) sep%alp = alp
       IF (PRESENT(eisol)) sep%eisol = eisol
       IF (PRESENT(gss)) sep%gss = gss
       IF (PRESENT(gsp)) sep%gsp = gsp
       IF (PRESENT(gpp)) sep%gpp = gpp
       IF (PRESENT(gp2)) sep%gp2 = gp2
       IF (PRESENT(hsp)) sep%hsp = hsp
       IF (PRESENT(gsd)) sep%gsd = gsd
       IF (PRESENT(gpd)) sep%gpd = gpd
       IF (PRESENT(gdd)) sep%gdd = gdd
       IF (PRESENT(ppddg)) sep%pre = ppddg
       IF (PRESENT(dpddg)) sep%d = dpddg
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sep is not associated")
    END IF

  END SUBROUTINE set_se_param

! *****************************************************************************
!> \brief Writes the semi-empirical type
! ***************************************************************************** 
  SUBROUTINE write_se_param(sep, subsys_section, error)

    TYPE(semi_empirical_type), POINTER       :: sep
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
    IF (ASSOCIATED(sep).AND.BTEST(cp_print_key_should_output(logger%iter_info,subsys_section,&
         "PRINT%KINDS/BASIS_SET",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,subsys_section,"PRINT%KINDS",&
            extension=".Log",error=error)

       IF (output_unit >0) THEN
          CALL get_se_param(sep,name=name,typ=typ,defined=defined,&
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
          CASE (do_method_mndo, do_method_mndod)
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
