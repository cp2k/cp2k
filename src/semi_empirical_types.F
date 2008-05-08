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
    LOGICAL                                :: defined, dorb, extended_basis_set
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
    ! One center - two electron integrals
    REAL(KIND = dp), DIMENSION(:,:),&
         POINTER                           :: w
  END TYPE semi_empirical_type

! *****************************************************************************
!> \brief Taper type
! *****************************************************************************
  TYPE taper_type
     REAL(KIND=dp)                         :: r0, al
  END TYPE taper_type

! *****************************************************************************
!> \brief Rotation Matrix Type
! *****************************************************************************
  TYPE rotmat_type
     ! Value of Rotation Matrices
     REAL(KIND=dp), DIMENSION (    3, 3)      :: sp
     REAL(KIND=dp), DIMENSION (    5, 5)      :: sd
     REAL(KIND=dp), DIMENSION ( 6, 3, 3)      :: pp
     REAL(KIND=dp), DIMENSION (15, 5, 3)      :: pd
     REAL(KIND=dp), DIMENSION (15, 5, 5)      :: dd
     ! Derivatives of Rotation Matrices
     REAL(KIND=dp), DIMENSION (3,     3, 3)   :: sp_d
     REAL(KIND=dp), DIMENSION (3,     5, 5)   :: sd_d
     REAL(KIND=dp), DIMENSION (3,  6, 3, 3)   :: pp_d
     REAL(KIND=dp), DIMENSION (3, 15, 5, 3)   :: pd_d
     REAL(KIND=dp), DIMENSION (3, 15, 5, 5)   :: dd_d
  END TYPE rotmat_type

  PUBLIC :: semi_empirical_type,&
            semi_empirical_create,&
            semi_empirical_release,&
            taper_type,&
            taper_create,&
            taper_release,&
            rotmat_type,&
            rotmat_create,&
            rotmat_release,&
            get_se_param,&
            set_se_param,&
            write_se_param

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical type
! *****************************************************************************
  SUBROUTINE semi_empirical_create(sep, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(sep),cp_failure_level,routineP,error,failure)
    ALLOCATE (sep,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sep%beta(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sep%sto_exponents(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sep%zn(0:3),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(sep%basis)
    NULLIFY(sep%w)
    CALL zero_se_param(sep,error)
    
  END SUBROUTINE semi_empirical_create

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
    sep%defined            = .FALSE.
    sep%dorb               = .FALSE.
    sep%extended_basis_set = .FALSE.
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
  SUBROUTINE semi_empirical_release(sep, error)

    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_release', &
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
       DEALLOCATE (sep%w,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (sep,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE semi_empirical_release

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
!> \brief Creates taper type
! ***************************************************************************** 
  SUBROUTINE taper_create(taper, rc, range, error)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rc, range
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(taper),cp_failure_level,routineP,error,failure)
    ALLOCATE (taper,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    taper%r0 = rc - 0.5_dp * range
    taper%al = 24._dp/range
    
  END SUBROUTINE taper_create

! *****************************************************************************
!> \brief Releases taper type
! ***************************************************************************** 
  SUBROUTINE taper_release(taper, error)
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(taper)) THEN
       DEALLOCATE (taper,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
   
  END SUBROUTINE taper_release

! *****************************************************************************
!> \brief Creates rotmat type
! ***************************************************************************** 
  SUBROUTINE rotmat_create(rotmat, error)
    TYPE(rotmat_type), POINTER               :: rotmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotmat_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(rotmat),cp_failure_level,routineP,error,failure)
    ALLOCATE (rotmat,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    rotmat%sp   = 0.0_dp
    rotmat%pp   = 0.0_dp
    rotmat%sd   = 0.0_dp
    rotmat%pd   = 0.0_dp
    rotmat%dd   = 0.0_dp
    rotmat%sp_d = 0.0_dp
    rotmat%pp_d = 0.0_dp
    rotmat%sd_d = 0.0_dp
    rotmat%pd_d = 0.0_dp
    rotmat%dd_d = 0.0_dp

  END SUBROUTINE rotmat_create

! *****************************************************************************
!> \brief Releases rotmat type
! ***************************************************************************** 
  SUBROUTINE rotmat_release(rotmat, error)
    TYPE(rotmat_type), POINTER               :: rotmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotmat_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(rotmat)) THEN
       DEALLOCATE (rotmat,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
   
  END SUBROUTINE rotmat_release

! *****************************************************************************
!> \brief Writes the semi-empirical type
!> \par History
!>        04.2008 Teodoro Laino [tlaino] - University of Zurich: rewriting with
!>                support for the whole set of parameters
! ***************************************************************************** 
  SUBROUTINE write_se_param(sep, subsys_section, error)

    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_se_param', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=1), DIMENSION(0:3), &
      PARAMETER                              :: orb_lab = (/"S","P","D","F"/)
    CHARACTER(LEN=2), DIMENSION(0:3), PARAMETER :: &
      z_lab = (/"ZS","ZP","ZD","ZF"/)
    CHARACTER(LEN=3), DIMENSION(0:3), PARAMETER :: &
      zeta_lab = (/"ZSN","ZPN","ZDN","ZFN"/)
    CHARACTER(LEN=5), DIMENSION(0:3), PARAMETER :: &
      beta_lab = (/"BETAS","BETAP","BETAD","BETAF"/)
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER                                  :: l, natorb, ngauss, &
                                                output_unit, typ, z
    LOGICAL                                  :: defined
    REAL(KIND=dp)                            :: alp, eheat, eisol, gp2, gpp, &
                                                gsp, gss, hsp, udd, uff, upp, &
                                                uss, zeff
    CHARACTER(LEN=3), DIMENSION(0:3), PARAMETER :: &
      u_lab = (/"USS","UPP","UDD","UFF"/)

    REAL(KIND=dp), DIMENSION(0:3)            :: u
    REAL(KIND=dp), DIMENSION(2)              :: dpddg, ppddg
    REAL(KIND=dp), DIMENSION(:), POINTER     :: beta, sexp
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (ASSOCIATED(sep).AND.BTEST(cp_print_key_should_output(logger%iter_info,subsys_section,&
         "PRINT%KINDS/SE_PARAMETERS",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,subsys_section,"PRINT%KINDS/SE_PARAMETERS",&
            extension=".Log",error=error)

       IF (output_unit >0) THEN
          CALL get_se_param(sep,name=name,typ=typ,defined=defined,&
               z=z,zeff=zeff,natorb=natorb,eheat=eheat,beta=beta,&
               sto_exponents=sexp,uss=uss,upp=upp,udd=udd,uff=uff,&
               alp=alp,eisol=eisol,gss=gss,gsp=gsp,gpp=gpp,gp2=gp2,&
               hsp=hsp,ppddg=ppddg,dpddg=dpddg,ngauss=ngauss)

          u(0) = uss
          u(1) = upp
          u(2) = udd
          u(3) = uff

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
          CASE (do_method_pm3)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Parametric Method 3 (PM3) ",&
                  TRIM(name)
          CASE (do_method_pm6)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Parametric Method 6 (PM6) ",&
                  TRIM(name)
          CASE (do_method_pdg)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","PDDG/PM3 ",&
                  TRIM(name)
          CASE (do_method_mndo)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","MNDO ",TRIM(name)
          CASE (do_method_mndod)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","MNDOD",TRIM(name)
          END SELECT

          ! If defined print all its semi-empirical parameters
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
                   WRITE (UNIT=output_unit,FMT="(T16,A,I2)") "Parameters for Shell: ",l
                   WRITE (UNIT=output_unit,FMT="(T22,A5,T30,A,T64,F17.4)")&
                         ADJUSTR(z_lab(l)),"- "//"Slater  Exponent for "//orb_lab(l)//"     [A]: ",sexp(l)
                   WRITE (UNIT=output_unit,FMT="(T22,A5,T30,A,T64,F17.4)")&
                         ADJUSTR(u_lab(l)),"- "//"One Center Energy for "//orb_lab(l)//"   [eV]: ",u(l)*evolt
                   WRITE (UNIT=output_unit,FMT="(T22,A5,T30,A,T64,F17.4)")&
                         ADJUSTR(beta_lab(l)),"- "//"Beta Parameter for "//orb_lab(l)//"      [eV]: ",beta(l)*evolt
                   WRITE (UNIT=output_unit,FMT="(T22,A5,T30,A,T64,F17.4)")&
                         ADJUSTR(zeta_lab(l)),"- "//"Internal Exponent for "//orb_lab(l)//" [a.u.]: ",sep%zn(l)
                END IF
             END DO
             WRITE (UNIT=output_unit,FMT="(/,T16,A)")"Additional Parameters (Derived or Fitted):"
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("ALP"),"- "//"Alpha Parameter for Core    [A^-1]: ",alp/angstrom
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("EISOL"),"- "//"Atomic Energy (Calculated)    [eV]: ",eisol*evolt
             ! One center Two electron Integrals
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("GSS"),"- "//"One Center Integral (SS ,SS ) [eV]: ",gss*evolt
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("GSP"),"- "//"One Center Integral (SS ,PP ) [eV]: ",gsp*evolt
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("GPP"),"- "//"One Center Integral (PP ,PP ) [eV]: ",gpp*evolt
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("GP2"),"- "//"One Center Integral (PP*,PP*) [eV]: ",gp2*evolt
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("HSP"),"- "//"One Center Integral (SP ,SP ) [eV]: ",hsp*evolt
             ! Slater Condon Parameters
             IF (sep%dorb) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("F0DD"),"- "//"Slater Condon Parameter F0DD  [eV]: ",sep%f0dd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("F2DD"),"- "//"Slater Condon Parameter F2DD  [eV]: ",sep%f2dd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("F4DD"),"- "//"Slater Condon Parameter F4DD  [eV]: ",sep%f4dd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("FOSD"),"- "//"Slater Condon Parameter FOSD  [eV]: ",sep%f0sd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("G2SD"),"- "//"Slater Condon Parameter G2SD  [eV]: ",sep%g2sd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("F0PD"),"- "//"Slater Condon Parameter F0PD  [eV]: ",sep%f0pd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("F2PD"),"- "//"Slater Condon Parameter F2PD  [eV]: ",sep%f2pd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("G1PD"),"- "//"Slater Condon Parameter G1PD  [eV]: ",sep%g1pd
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("G3PD"),"- "//"Slater Condon Parameter G3PD  [eV]: ",sep%g3pd
             END IF
             ! Charge Separation
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("DD2"),"- "//"Charge Separation  SP, L=1  [bohr]: ",sep%cs(2)
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("DD3"),"- "//"Charge Separation  PP, L=2  [bohr]: ",sep%cs(3)
             IF (sep%dorb) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("DD4"),"- "//"Charge Separation  SD, L=2  [bohr]: ",sep%cs(4)
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("DD5"),"- "//"Charge Separation  PD, L=1  [bohr]: ",sep%cs(5)
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("DD6"),"- "//"Charge Separation  DD, L=2  [bohr]: ",sep%cs(6)
             END IF
             ! Klopman-Ohno Terms
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("PO1"),"- "//"Klopman-Ohno term, SS, L=0  [bohr]: ",sep%ko(1)
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("PO2"),"- "//"Klopman-Ohno term, SP, L=1  [bohr]: ",sep%ko(2)
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("PO3"),"- "//"Klopman-Ohno term, PP, L=2  [bohr]: ",sep%ko(3)
             IF (sep%dorb) THEN
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("PO4"),"- "//"Klopman-Ohno term, SD, L=2  [bohr]: ",sep%ko(4)
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("PO5"),"- "//"Klopman-Ohno term, PD, L=1  [bohr]: ",sep%ko(5)
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("PO6"),"- "//"Klopman-Ohno term, DD, L=2  [bohr]: ",sep%ko(6)
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("PO7"),"- "//"Klopman-Ohno term, PP, L=0  [bohr]: ",sep%ko(7)
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                      ADJUSTR("PO8"),"- "//"Klopman-Ohno term, DD, L=0  [bohr]: ",sep%ko(8)
             END IF
             WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                   ADJUSTR("PO9"),"- "//"Klopman-Ohno term, CORE     [bohr]: ",sep%ko(9)
          ELSE
             WRITE (UNIT=output_unit,FMT="(T55,A)") "Parameters are not defined"
          END IF

          ! Additional Parameters not common to all semi-empirical methods
          SELECT CASE (typ)
          CASE (do_method_pdg)
                WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T52,F14.10,T67,F14.10)") &
                     ADJUSTR("d_PDDG"),"- "//"Exponent [A^-1]:",dpddg/angstrom,&
                     ADJUSTR("P_PDDG"),"- "//"Parameter  [eV]:",ppddg*evolt
          END SELECT
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
            "PRINT%KINDS/SE_PARAMETERS",error=error)
    END IF
  END SUBROUTINE write_se_param

END MODULE semi_empirical_types
