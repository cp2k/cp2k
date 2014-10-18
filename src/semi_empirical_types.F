!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE dg_types,                        ONLY: dg_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_mndo, do_method_mndod, do_method_pdg, &
       do_method_pm3, do_method_pm6, do_method_pnnl, do_method_rm1, &
       do_se_IS_kdso_d, do_se_IS_slater
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE multipole_types,                 ONLY: do_multipole_charge,&
                                             do_multipole_dipole,&
                                             do_multipole_none,&
                                             do_multipole_quadrupole
  USE physcon,                         ONLY: angstrom,&
                                             evolt,&
                                             kcalmol
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE semi_empirical_expns3_types,     ONLY: semi_empirical_expns3_p_type,&
                                             semi_empirical_expns3_release
  USE semi_empirical_mpole_types,      ONLY: semi_empirical_mpole_p_release,&
                                             semi_empirical_mpole_p_type
  USE taper_types,                     ONLY: taper_create,&
                                             taper_release,&
                                             taper_type
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_types'

! *****************************************************************************
!> \brief Semi-empirical type
! *****************************************************************************
  TYPE semi_empirical_type
    INTEGER                                :: typ
    INTEGER                                :: nr
    INTEGER                                :: core_size, atm_int_size
    CHARACTER(LEN=default_string_length)   :: name
    LOGICAL                                :: defined, dorb, extended_basis_set
    LOGICAL                                :: p_orbitals_on_h
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
    REAL(KIND = dp)                        :: ass,asp,app,de,acoul
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
    TYPE(semi_empirical_mpole_p_type),&
         POINTER, DIMENSION(:)             :: w_mpole
    ! 1/R^3 residual integral part
    TYPE(semi_empirical_expns3_p_type),&
         POINTER, DIMENSION(:)             :: expns3_int
  END TYPE semi_empirical_type

  TYPE semi_empirical_p_type
     TYPE(semi_empirical_type), POINTER    :: se_param
  END TYPE semi_empirical_p_type

! *****************************************************************************
!> \brief  Rotation Matrix Type
!> \author 05.2008 Teodoro Laino [tlaino] - University of Zurich
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

! *****************************************************************************
!> \brief  Ewald control type (for periodic SE)
!> \author Teodoro Laino [tlaino]  - 12.2008
! *****************************************************************************
  TYPE ewald_gks_type
     REAL(KIND=dp)                            :: alpha
     TYPE(dg_type), POINTER                   :: dg
     TYPE(pw_pool_type), POINTER              :: pw_pool
  END TYPE ewald_gks_type

  TYPE se_int_control_type
     LOGICAL                                  :: shortrange
     LOGICAL                                  :: do_ewald_r3
     LOGICAL                                  :: do_ewald_gks
     LOGICAL                                  :: pc_coulomb_int
     INTEGER                                  :: integral_screening
     INTEGER                                  :: max_multipole
     TYPE(ewald_gks_type)                     :: ewald_gks
  END TYPE se_int_control_type

! *****************************************************************************
!> \brief Store the value of the tapering function and possibly its derivative
!>        for screened integrals
! *****************************************************************************
  TYPE se_int_screen_type
     REAL(KIND=dp)                            :: ft, dft
  END TYPE se_int_screen_type

! *****************************************************************************
!> \brief Taper type use in semi-empirical calculations
! *****************************************************************************
  TYPE se_taper_type
     TYPE(taper_type), POINTER             :: taper
     TYPE(taper_type), POINTER             :: taper_cou
     TYPE(taper_type), POINTER             :: taper_exc
     TYPE(taper_type), POINTER             :: taper_lrc
     ! This taper is for KDSO-D integrals
     TYPE(taper_type), POINTER             :: taper_add
  END TYPE se_taper_type

  PUBLIC :: semi_empirical_type,&
            semi_empirical_p_type,&
            semi_empirical_create,&
            semi_empirical_release,&
            rotmat_type,&
            rotmat_create,&
            rotmat_release,&
            get_se_param,&
            write_se_param,&
            se_int_control_type,&
            setup_se_int_control_type,&
            se_int_screen_type,&
            se_taper_type,&
            se_taper_release,&
            se_taper_create

CONTAINS

! *****************************************************************************
!> \brief Allocate semi-empirical type
!> \param sep ...
!> \param error ...
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
    NULLIFY(sep%w_mpole)
    NULLIFY(sep%expns3_int)
    CALL zero_se_param(sep,error)

  END SUBROUTINE semi_empirical_create

! *****************************************************************************
!> \brief Deallocate the semi-empirical type
!> \param sep ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE semi_empirical_release(sep, error)

    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(sep)) THEN
       CALL deallocate_sto_basis_set(sep%basis, error)
       CALL semi_empirical_mpole_p_release(sep%w_mpole, error)
       IF (ASSOCIATED(sep%beta)) THEN
          DEALLOCATE (sep%beta,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(sep%sto_exponents)) THEN
          DEALLOCATE (sep%sto_exponents,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(sep%zn)) THEN
          DEALLOCATE (sep%zn,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(sep%w)) THEN
          DEALLOCATE (sep%w,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(sep%expns3_int)) THEN
          DO i = 1, SIZE(sep%expns3_int)
             CALL semi_empirical_expns3_release(sep%expns3_int(i)%expns3,error)
          END DO
          DEALLOCATE (sep%expns3_int,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE (sep,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE semi_empirical_release

! *****************************************************************************
!> \brief Zero the whole semi-empirical type
!> \param sep ...
!> \param error ...
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
    sep%p_orbitals_on_h    = .FALSE.
    sep%name          = ""
    sep%typ           = HUGE(0)
    sep%core_size     = HUGE(0)
    sep%atm_int_size  = HUGE(0)
    sep%z             = HUGE(0)
    sep%zeff          = HUGE(0.0_dp)
    sep%natorb        = 0
    sep%ngauss        = 0
    sep%eheat         = HUGE(0.0_dp)

    sep%zn            = 0.0_dp
    sep%sto_exponents = 0.0_dp
    sep%beta          = 0.0_dp

    sep%uss    = 0.0_dp !eV
    sep%upp    = 0.0_dp !eV
    sep%udd    = 0.0_dp !eV
    sep%uff    = 0.0_dp
    sep%alp    = 0.0_dp
    sep%eisol  = 0.0_dp
    sep%nr    = 1
    sep%acoul = 0.0_dp
    sep%de    = 0.0_dp
    sep%ass    = 0.0_dp
    sep%asp    = 0.0_dp
    sep%app    = 0.0_dp
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
!> \brief Get info from the semi-empirical type
!> \param sep ...
!> \param name ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param eheat ...
!> \param beta ...
!> \param sto_exponents ...
!> \param uss ...
!> \param upp ...
!> \param udd ...
!> \param uff ...
!> \param alp ...
!> \param eisol ...
!> \param gss ...
!> \param gsp ...
!> \param gpp ...
!> \param gp2 ...
!> \param acoul ...
!> \param nr ...
!> \param de ...
!> \param ass ...
!> \param asp ...
!> \param app ...
!> \param hsp ...
!> \param gsd ...
!> \param gpd ...
!> \param gdd ...
!> \param ppddg ...
!> \param dpddg ...
!> \param ngauss ...
! *****************************************************************************
  SUBROUTINE get_se_param(sep,name,typ,defined,z,zeff,natorb,eheat,&
       beta,sto_exponents,uss,upp,udd,uff,alp,eisol,gss,gsp,gpp,gp2,&
       acoul,nr,de,ass,asp,app,hsp,gsd,gpd,gdd,ppddg,dpddg,ngauss)

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
                                                acoul
    INTEGER, INTENT(OUT), OPTIONAL           :: nr
    REAL(KIND=dp), OPTIONAL                  :: de, ass, asp, app, hsp, gsd, &
                                                gpd, gdd
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
       IF (PRESENT(nr)) nr = sep%nr
       IF (PRESENT(acoul)) acoul = sep%acoul
       IF (PRESENT(de)) de = sep%de
       IF (PRESENT(ass)) ass = sep%ass
       IF (PRESENT(asp)) asp = sep%asp
       IF (PRESENT(app)) app = sep%app
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
!> \param sep ...
!> \param name ...
!> \param typ ...
!> \param defined ...
!> \param z ...
!> \param zeff ...
!> \param natorb ...
!> \param eheat ...
!> \param beta ...
!> \param sto_exponents ...
!> \param uss ...
!> \param upp ...
!> \param udd ...
!> \param uff ...
!> \param alp ...
!> \param eisol ...
!> \param gss ...
!> \param gsp ...
!> \param gpp ...
!> \param gp2 ...
!> \param acoul ...
!> \param nr ...
!> \param de ...
!> \param ass ...
!> \param asp ...
!> \param app ...
!> \param hsp ...
!> \param gsd ...
!> \param gpd ...
!> \param gdd ...
!> \param ppddg ...
!> \param dpddg ...
!> \param ngauss ...
! *****************************************************************************
  SUBROUTINE set_se_param(sep,name,typ,defined,z,zeff,natorb,eheat,&
       beta,sto_exponents,uss,upp,udd,uff,alp,eisol,gss,gsp,gpp,gp2,&
       acoul,nr,de,ass,asp,app,hsp,gsd,gpd,gdd,ppddg,dpddg,ngauss)

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
                                                acoul
    INTEGER, INTENT(IN), OPTIONAL            :: nr
    REAL(KIND=dp), OPTIONAL                  :: de, ass, asp, app, hsp, gsd, &
                                                gpd, gdd
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
       IF (PRESENT(acoul)) sep%acoul = acoul
       IF (PRESENT(nr)) sep%nr = nr
       IF (PRESENT(de)) sep%de = de
       IF (PRESENT(ass)) sep%ass = ass
       IF (PRESENT(asp)) sep%asp = asp
       IF (PRESENT(app)) sep%app = app
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
!> \brief Creates rotmat type
!> \param rotmat ...
!> \param error ...
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

  END SUBROUTINE rotmat_create

! *****************************************************************************
!> \brief Releases rotmat type
!> \param rotmat ...
!> \param error ...
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
!> \brief Setup the Semiempirical integral control type
!> \param se_int_control ...
!> \param shortrange ...
!> \param do_ewald_r3 ...
!> \param do_ewald_gks ...
!> \param integral_screening ...
!> \param max_multipole ...
!> \param pc_coulomb_int ...
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE setup_se_int_control_type(se_int_control, shortrange, do_ewald_r3,&
       do_ewald_gks, integral_screening, max_multipole, pc_coulomb_int)
    TYPE(se_int_control_type)                :: se_int_control
    LOGICAL, INTENT(IN)                      :: shortrange, do_ewald_r3, &
                                                do_ewald_gks
    INTEGER, INTENT(IN)                      :: integral_screening, &
                                                max_multipole
    LOGICAL, INTENT(IN)                      :: pc_coulomb_int

    se_int_control%shortrange         = shortrange
    se_int_control%do_ewald_r3        = do_ewald_r3
    se_int_control%integral_screening = integral_screening
    ! This makes the assignment independent of the value of the different constants
    SELECT CASE(max_multipole)
    CASE (do_multipole_none)
       se_int_control%max_multipole   = -1
    CASE (do_multipole_charge)
       se_int_control%max_multipole   =  0
    CASE (do_multipole_dipole)
       se_int_control%max_multipole   =  1
    CASE (do_multipole_quadrupole)
       se_int_control%max_multipole   =  2
    END SELECT

    se_int_control%do_ewald_gks       = do_ewald_gks
    se_int_control%pc_coulomb_int     = pc_coulomb_int
    NULLIFY(se_int_control%ewald_gks%dg,se_int_control%ewald_gks%pw_pool)

  END SUBROUTINE setup_se_int_control_type

! *****************************************************************************
!> \brief Creates the taper type used in SE calculations
!> \param se_taper ...
!> \param integral_screening ...
!> \param do_ewald ...
!> \param taper_cou ...
!> \param range_cou ...
!> \param taper_exc ...
!> \param range_exc ...
!> \param taper_scr ...
!> \param range_scr ...
!> \param taper_lrc ...
!> \param range_lrc ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 03.2009
! *****************************************************************************
  SUBROUTINE se_taper_create(se_taper, integral_screening, do_ewald, &
       taper_cou, range_cou, taper_exc, range_exc, taper_scr, range_scr,&
       taper_lrc, range_lrc, error)
    TYPE(se_taper_type), POINTER             :: se_taper
    INTEGER, INTENT(IN)                      :: integral_screening
    LOGICAL, INTENT(IN)                      :: do_ewald
    REAL(KIND=dp), INTENT(IN)                :: taper_cou, range_cou, &
                                                taper_exc, range_exc, &
                                                taper_scr, range_scr, &
                                                taper_lrc, range_lrc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_taper_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(se_taper),cp_failure_level,routineP,error,failure)
    ALLOCATE(se_taper, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(se_taper%taper)
    NULLIFY(se_taper%taper_cou)
    NULLIFY(se_taper%taper_exc)
    NULLIFY(se_taper%taper_lrc)
    NULLIFY(se_taper%taper_add)
    ! Create the sub-typo taper
    CALL taper_create(se_taper%taper_cou, taper_cou, range_cou, error)
    CALL taper_create(se_taper%taper_exc, taper_exc, range_exc, error)
    IF (integral_screening==do_se_IS_kdso_d) THEN
       CALL taper_create(se_taper%taper_add, taper_scr, range_scr, error)
    END IF
    IF ((integral_screening/=do_se_IS_slater).AND.do_ewald) THEN
       CALL taper_create(se_taper%taper_lrc, taper_lrc, range_lrc, error)
    END IF
  END SUBROUTINE se_taper_create

! *****************************************************************************
!> \brief Releases the taper type used in SE calculations
!> \param se_taper ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 03.2009
! *****************************************************************************
  SUBROUTINE se_taper_release(se_taper, error)
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_taper_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(se_taper)) THEN
       CALL taper_release(se_taper%taper_cou, error=error)
       CALL taper_release(se_taper%taper_exc, error=error)
       CALL taper_release(se_taper%taper_lrc, error=error)
       CALL taper_release(se_taper%taper_add, error=error)

       DEALLOCATE(se_taper, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE se_taper_release

! *****************************************************************************
!> \brief Writes the semi-empirical type
!> \param sep ...
!> \param subsys_section ...
!> \param error ...
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
    CHARACTER(LEN=default_string_length)     :: i_string, name
    INTEGER                                  :: i, l, natorb, ngauss, nr, &
                                                output_unit, typ, z
    LOGICAL                                  :: defined
    REAL(KIND=dp)                            :: acoul, alp, app, asp, ass, &
                                                de, eheat, eisol, gp2, gpp, &
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
               de=de,ass=ass,asp=asp,app=app,hsp=hsp,ppddg=ppddg, &
               acoul=acoul,nr=nr,dpddg=dpddg,ngauss=ngauss)

          u(0) = uss
          u(1) = upp
          u(2) = udd
          u(3) = uff

          SELECT CASE (typ)
          CASE DEFAULT
             CALL stop_program(routineN,moduleN,__LINE__,"Semiempirical method unknown")
          CASE(do_method_am1)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Austin Model 1 (AM1)",TRIM(name)
          CASE(do_method_rm1)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Recife Model 1 (RM1)",TRIM(name)
          CASE (do_method_pm3)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Parametric Method 3 (PM3) ",TRIM(name)
          CASE (do_method_pnnl)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","PNNL method ",TRIM(name)
          CASE (do_method_pm6)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","Parametric Method 6 (PM6) ",TRIM(name)
          CASE (do_method_pdg)
             WRITE (UNIT=output_unit,FMT="(/,A,T35,A,T67,A14)")&
                  " Semi empirical parameters: ","PDDG/PM3 ",TRIM(name)
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
             SELECT CASE (typ)
             CASE (do_method_am1, do_method_rm1, do_method_pm3, do_method_pdg, do_method_pnnl)
                IF(typ == do_method_pnnl) THEN
                   WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                         ADJUSTR("ASS"),"- "//" SS polarization [au]: ",sep%ass
                   WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                         ADJUSTR("ASP"),"- "//" SP polarization [au]: ",sep%asp
                   WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                         ADJUSTR("APP"),"- "//" PP polarization[au]: ",sep%app
                   WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                         ADJUSTR("DE"),"- "//" Dispersion Parameter [eV]: ",sep%de*evolt
                   WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                         ADJUSTR("ACOUL"),"- "//" Slater parameter: ",sep%acoul
                   WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,I12)")&
                         ADJUSTR("NR"),"- "//" Slater parameter: ",sep%nr
                ELSEIF((typ == do_method_am1 .OR. typ == do_method_rm1) .AND. sep%z == 5) THEN
                   ! Standard case
                   DO i = 1, SIZE(sep%bfn1,1)
                      i_string = cp_to_string(i)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN1"//TRIM(ADJUSTL(i_string))//"_ALL"),&
                           "- "//"Core-Core VDW, Multiplier   [a.u.]: ",sep%bfn1(i,1)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN2"//TRIM(ADJUSTL(i_string))//"_ALL"),&
                           "- "//"Core-Core VDW, Exponent     [a.u.]: ",sep%bfn2(i,1)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN3"//TRIM(ADJUSTL(i_string))//"_ALL"),&
                           "- "//"Core-Core VDW, Position     [a.u.]: ",sep%bfn3(i,1)
                   END DO
                   ! Special Case : Hydrogen
                   DO i = 1, SIZE(sep%bfn1,1)
                      i_string = cp_to_string(i)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN1"//TRIM(ADJUSTL(i_string))//"_H"),&
                           "- "//"Core-Core VDW, Multiplier   [a.u.]: ",sep%bfn1(i,2)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN2"//TRIM(ADJUSTL(i_string))//"_H"),&
                           "- "//"Core-Core VDW, Exponent     [a.u.]: ",sep%bfn2(i,2)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN3"//TRIM(ADJUSTL(i_string))//"_H"),&
                           "- "//"Core-Core VDW, Position     [a.u.]: ",sep%bfn3(i,2)
                   END DO
                   ! Special Case : Carbon
                   DO i = 1, SIZE(sep%bfn1,1)
                      i_string = cp_to_string(i)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN1"//TRIM(ADJUSTL(i_string))//"_C"),&
                           "- "//"Core-Core VDW, Multiplier   [a.u.]: ",sep%bfn1(i,3)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN2"//TRIM(ADJUSTL(i_string))//"_C"),&
                           "- "//"Core-Core VDW, Exponent     [a.u.]: ",sep%bfn2(i,3)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN3"//TRIM(ADJUSTL(i_string))//"_C"),&
                           "- "//"Core-Core VDW, Position     [a.u.]: ",sep%bfn3(i,3)
                   END DO
                   ! Special Case : Halogens
                   DO i = 1, SIZE(sep%bfn1,1)
                      i_string = cp_to_string(i)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN1"//TRIM(ADJUSTL(i_string))//"_HALO"),&
                           "- "//"Core-Core VDW, Multiplier   [a.u.]: ",sep%bfn1(i,4)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN2"//TRIM(ADJUSTL(i_string))//"_HALO"),&
                           "- "//"Core-Core VDW, Exponent     [a.u.]: ",sep%bfn2(i,4)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN3"//TRIM(ADJUSTL(i_string))//"_HALO"),&
                           "- "//"Core-Core VDW, Position     [a.u.]: ",sep%bfn3(i,4)
                   END DO
                ELSE
                   DO i = 1, SIZE(sep%fn1,1)
                      i_string = cp_to_string(i)
                      ! Skip the printing of params that are zero..
                      IF (sep%fn1(i)==0.0_dp.AND.sep%fn2(i)==0.0_dp.AND.sep%fn3(i)==0.0_dp) CYCLE
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN1"//TRIM(ADJUSTL(i_string))),&
                           "- "//"Core-Core VDW, Multiplier   [a.u.]: ",sep%fn1(i)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN2"//TRIM(ADJUSTL(i_string))),&
                           "- "//"Core-Core VDW, Exponent     [a.u.]: ",sep%fn2(i)
                      WRITE (UNIT=output_unit,FMT="(T16,A11,T30,A,T69,F12.4)")&
                           ADJUSTR("FN3"//TRIM(ADJUSTL(i_string))),&
                           "- "//"Core-Core VDW, Position     [a.u.]: ",sep%fn3(i)
                   END DO
                END IF
             END SELECT
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
