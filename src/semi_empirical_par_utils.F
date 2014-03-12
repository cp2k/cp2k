!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities to post-process semi-empirical parameters
!> \par History
!>         [tlaino] 03.2008 - Splitting from semi_empirical_parameters and
!>                            keeping there only the setting of the SE params
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
MODULE semi_empirical_par_utils
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fac
  USE mathlib,                         ONLY: binomial
  USE physcon,                         ONLY: bohr,&
                                             evolt
  USE semi_empirical_int_arrays,       ONLY: int_ij,&
                                             int_kl,&
                                             int_onec2el
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PRIVATE :: nelem = 106

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_par_utils'

!                STANDARD MOPAC PARAMETERS USED FOR AM1, RM1, MNDO, PM3, PM6
!
!      H                                                                      He
!      Li Be                                                 B  C  N  O  F    Ne
!      Na Mg                                                 Al Si P  S  Cl   Ar
!      K  Ca Sc                Ti V  Cr Mn Fe Co Ni Cu Zn    Ga Ge As Se Br   Kr
!      Rb Sr Y                 Zr Nb Mo Tc Ru Rh Pd Ag Cd    In Sn Sb Te I    Xe
!      Cs Ba La Ce-Lu          Hf Ta W  Re Os Ir Pt Au Hg    Tl Pb Bi Po At   Rn
!      Fr Ra Ac Th Pa U        Np Pu Am Cm Bk Cf Es Fm Md    No Lr Rf Ha 106

!                                      "s" shell
  INTEGER, DIMENSION(0:nelem), PRIVATE :: Nos =(/                           -1,&!    0
        1,                                                                   2,&!    2
        1, 2,                                                 2, 2, 2, 2, 2, 0,&!   10
        1, 2,                                                 2, 2, 2, 2, 2, 0,&!   18
        1, 2, 2,                2, 2, 1, 2, 2, 2, 2, 1, 2,    2, 2, 2, 2, 2, 0,&!   36
        1, 2, 2,                2, 1, 1, 2, 1, 1, 0, 1, 2,    2, 2, 2, 2, 2, 0,&!   54
        1, 2, 2, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,                     &
                                2, 2, 1, 2, 2, 2, 1, 1, 2,    2, 2, 2, 2, 2, 0,&!   86
        1, 1, 2, 4, 2, 2,       2, 2, 2, 2, 2, 2, 0, 3,-3,    1, 2, 1,-2,-1 /)

!                                      "p" shell
  INTEGER, DIMENSION(0:nelem), PRIVATE :: Nop =(/                           -1,&!    0
        0,                                                                   0,&!    2
        0, 0,                                                 1, 2, 3, 4, 5, 6,&!   10
        0, 0,                                                 1, 2, 3, 4, 5, 6,&!   18
        0, 0, 0,                0, 0, 0, 0, 0, 0, 0, 0, 0,    1, 2, 3, 4, 5, 6,&!   36
        0, 0, 0,                0, 0, 0, 0, 0, 0, 0, 0, 0,    1, 2, 3, 4, 5, 6,&!   54
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                     &
                                0, 0, 0, 0, 0, 0, 0, 0, 0,    1, 2, 3, 4, 5, 6,&!   86
        0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0 /)

!                                      "d" shell
  INTEGER, DIMENSION(0:nelem), PRIVATE :: Nod =(/                           -1,&!    0
        0,                                                                   0,&!    2
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   10
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   18
        0, 0, 1,                2, 3, 5, 5, 6, 7, 8,10, 0,    0, 0, 0, 0, 0, 0,&!   36
        0, 0, 1,                2, 4, 5, 5, 7, 8,10,10, 0,    0, 0, 0, 0, 0, 0,&!   54
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,                     &
                                2, 3, 5, 5, 6, 7, 9,10, 0,    0, 0, 0, 0, 0, 0,&!   86
        0, 0, 1, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0 /)

!      H          <Quantum Numbers for s, p, d and f orbitals>                He
!      Li Be                                                 B  C  N  O  F    Ne
!      Na Mg                                                 Al Si P  S  Cl   Ar
!      K  Ca Sc                Ti V  Cr Mn Fe Co Ni Cu Zn    Ga Ge As Se Br   Kr
!      Rb Sr Y                 Zr Nb Mo Tc Ru Rh Pd Ag Cd    In Sn Sb Te I    Xe
!      Cs Ba La Ce-Lu          Hf Ta W  Re Os Ir Pt Au Hg    Tl Pb Bi Po At   Rn
!      Fr Ra Ac Th Pa U        Np Pu Am Cm Bk Cf Es Fm Md    No Lr Rf Ha 106

  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqs = (/             -1,&!    0
        1,                                                                   1,&!    2
        2, 2,                                                 2, 2, 2, 2, 2, 2,&!   10
        3, 3,                                                 3, 3, 3, 3, 3, 3,&!   18
        4, 4, 4,                4, 4, 4, 4, 4, 4, 4, 4, 4,    4, 4, 4, 4, 4, 4,&!   36
        5, 5, 5,                5, 5, 5, 5, 5, 5, 5, 5, 5,    5, 5, 5, 5, 5, 5,&!   54
        6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 6,                     &
                                6, 6, 6, 6, 6, 6, 6, 6, 6,    6, 6, 6,-1,-1,-1,&!   86
       -1,-1,-1,-1,-1,-1,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqp = (/             -1,&!    0
       -1,                                                                  -1,&!    2
        2, 2,                                                 2, 2, 2, 2, 2, 2,&!   10
        3, 3,                                                 3, 3, 3, 3, 3, 3,&!   18
        4, 4, 4,                4, 4, 4, 4, 4, 4, 4, 4, 4,    4, 4, 4, 4, 4, 4,&!   36
        5, 5, 5,                5, 5, 5, 5, 5, 5, 5, 5, 5,    5, 5, 5, 5, 5, 5,&!   54
        6, 6, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 6,                     &
                                6, 6, 6, 6, 6, 6, 6, 6, 6,    6, 6, 6,-1,-1,-1,&!   86
        -1,-1,-1,-1,-1,-1,     -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqd = (/             -1,&!    0
       -1,                                                                  -1,&!    2
       -1,-1,                                                -1,-1,-1,-1,-1,-1,&!   10
       -1,-1,                                                 3, 3, 3, 3, 3,-1,&!   18
       -1,-1, 3,                3, 3, 3, 3, 3, 3, 3, 3, 3,    4, 4, 4, 4, 4,-1,&!   36
       -1,-1, 4,                4, 4, 4, 4, 4, 4, 4, 4, 4,    5, 5, 5, 5, 5,-1,&!   54
       -1,-1, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,                     &
                                5, 5, 5, 5, 5, 5, 5, 5, 5,   -1,-1,-1,-1,-1,-1,&!   86
       -1,-1,-1,-1,-1,-1,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE ::   nqf = (/             -1,&!    0
       -1,                                                                  -1,&!    2
       -1,-1,                                                -1,-1,-1,-1,-1,-1,&!   10
       -1,-1,                                                -1,-1,-1,-1,-1,-1,&!   18
       -1,-1,-1,               -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1,-1,&!   36
       -1,-1,-1,               -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1,-1,&!   54
       -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,                     &
                               -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1,-1,&!   86
       -1,-1,-1,-1,-1,-1,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  ! Element Valence
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: zval = (/              -1,&!    0
        1,                                                                   2,&!    2
        1, 2,                                                 3, 4, 5, 6, 7, 8,&!   10
        1, 2,                                                 3, 4, 5, 6, 7, 8,&!   18
        1, 2, 3,                4, 5, 6, 7, 8, 9,10,11, 2,    3, 4, 5, 6, 7, 8,&!   36
        1, 2, 3,                4, 5, 6, 7, 8, 9,10,11, 2,    3, 4, 5, 6, 7, 8,&!   54
        1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16, 3,                     &
                                4, 5, 6, 7, 8, 9,10,11, 2,    3, 4, 5, 6, 7,-1,&!   86
       -1,-1,-1, 4,-1, 6,      -1,-1,-1,-1,-1,-1,-1,-1,-1,   -1,-1,-1,-1,-1 /)

  ! Number of 1 center 2 electron integrals involving partially filled d shells
  ! r016:  <SS|DD>
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: ir016 =  (/             0,&!    0
        0,                                                                   0,&!    2
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   10
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   18
        0, 0, 2,                4, 6, 5,10,12,14,16,10, 0,    0, 0, 0, 0, 0, 0,&!   36
        0, 0, 4,                4, 4, 5,10, 7, 8, 0,10, 0,    0, 0, 0, 0, 0, 0,&!   54
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                     &
                                4, 6, 8,10,12,14, 9,10, 0,    0, 0, 0, 0, 0, 0,&!   86
        0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0 /)

  ! r066:  <DD|DD> "0" term
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: ir066 =  (/             0,&!    0
        0,                                                                   0,&!    2
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   10
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   18
        0, 0, 0,                1, 3,10,10,15,21,28,45, 0,    0, 0, 0, 0, 0, 0,&!   36
        0, 0, 0,                1, 6,10,10,21,28,45,45, 0,    0, 0, 0, 0, 0, 0,&!   54
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                     &
                                1, 3, 6,10,15,21,36,45, 0,    0, 0, 0, 0, 0, 0,&!   86
        0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0 /)

  ! r244:  <SD|SD>
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: ir244 =  (/             0,&!    0
        0,                                                                   0,&!    2
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   10
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   18
        0, 0, 1,                2, 3, 5, 5, 6, 7, 8, 5, 0,    0, 0, 0, 0, 0, 0,&!   36
        0, 0, 1,                2, 4, 5, 5, 5, 5, 0, 5, 0,    0, 0, 0, 0, 0, 0,&!   54
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                     &
                                2, 3, 4, 5, 6, 7, 5, 5, 0,    0, 0, 0, 0, 0, 0,&!   86
        0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0 /)

  ! r266:  <DD|DD> "2" term
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: ir266 =  (/             0,&!    0
        0,                                                                   0,&!    2
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   10
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   18
        0, 0, 0,                8,15,35,35,35,43,50,70, 0,    0, 0, 0, 0, 0, 0,&!   36
        0, 0, 0,                8,21,35,35,43,50,70,70, 0,    0, 0, 0, 0, 0, 0,&!   54
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                     &
                                8,15,21,35,35,43,56,70, 0,    0, 0, 0, 0, 0, 0,&!   86
        0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0 /)

  ! r466:  <DD|DD> "4" term
  INTEGER, DIMENSION(0:nelem), PARAMETER, PRIVATE :: ir466 =  (/             0,&!    0
        0,                                                                   0,&!    2
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   10
        0, 0,                                                 0, 0, 0, 0, 0, 0,&!   18
        0, 0, 0,                1, 8,35,35,35,36,43,70, 0,    0, 0, 0, 0, 0, 0,&!   36
        0, 0, 0,                1,21,35,35,36,43,70,70, 0,    0, 0, 0, 0, 0, 0,&!   54
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                     &
                                1, 8,21,35,35,36,56,70, 0,    0, 0, 0, 0, 0, 0,&!   86
        0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0 /)

  INTERFACE amn_l
     MODULE PROCEDURE amn_l1, amn_l2
  END INTERFACE

  PUBLIC :: convert_param_to_cp2k, calpar, valence_electrons, get_se_basis,&
            setup_1c_2el_int, amn_l

CONTAINS

! *****************************************************************************
!> \brief  Gives back the number of valence electrons for element z and also the
!>         number of atomic orbitals for that specific element
! *****************************************************************************
  SUBROUTINE valence_electrons(sep, extended_basis_set, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    LOGICAL, INTENT(IN)                      :: extended_basis_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'valence_electrons', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: natorb, z
    LOGICAL                                  :: check, failure, use_p_orbitals
    REAL(KIND=dp)                            :: zeff

    failure        = .FALSE.
    use_p_orbitals = .TRUE.
    z              = sep%z
    CPPostcondition(z>=0, cp_failure_level, routineP, error, failure)
    ! Special case for Hydrogen.. If requested allow p-orbitals on it..
    SELECT CASE (z)
    CASE (0,2)
       use_p_orbitals = .FALSE.
    CASE (1)
       use_p_orbitals = sep%p_orbitals_on_h
    CASE DEFAULT
       ! Nothing to do..
    END SELECT
    ! Determine the number of atomic orbitals
    natorb = 0
    IF (nqs(z)>0)                                        natorb = natorb+1
    IF((nqp(z)>0).OR.use_p_orbitals)                     natorb = natorb+3
    IF (extended_basis_set.AND.element_has_d(sep,error)) natorb = natorb+5
    IF (extended_basis_set.AND.element_has_f(sep,error)) natorb = natorb+7
    ! Check and assignemnt
    check = (natorb<=4).OR.(extended_basis_set)
    CPPostcondition(check, cp_failure_level, routineP, error, failure)
    sep%natorb             = natorb
    sep%extended_basis_set = extended_basis_set
    ! Determine the Z eff
    zeff = REAL(zval(z),KIND=dp)
    sep%zeff               = zeff
  END SUBROUTINE valence_electrons

! *****************************************************************************
!> \brief  Gives back the number of basis function for each l
! *****************************************************************************
  FUNCTION get_se_basis(sep,l,error) RESULT(n)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER, INTENT(IN)                      :: l
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'get_se_basis', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF ( sep%z < 0 .OR. sep%z > nelem ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Invalid atomic number !")
    ELSE
       IF     ( l==0 ) THEN
          n = nqs(sep%z)
       ELSEIF ( l==1 ) THEN
          ! Special case for Hydrogen.. If requested allow p-orbitals on it..
          IF ((sep%z==1).AND.sep%p_orbitals_on_h) THEN
             n = 1
          ELSE
             n = nqp(sep%z)
          END IF
       ELSEIF ( l==2 ) THEN
          n = nqd(sep%z)
       ELSEIF ( l==3 ) THEN
          n = nqf(sep%z)
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,"Invalid l quantum number !")
       END IF
       IF ( n < 0 ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,"Invalid n quantum number !")
       END IF
    END IF
  END FUNCTION get_se_basis

! *****************************************************************************
!> \brief  Converts parameter units to internal
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE convert_param_to_cp2k(sep,error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'convert_param_to_cp2k', &
      routineP = moduleN//':'//routineN

    sep%beta  = sep%beta/evolt
    sep%uss   = sep%uss/evolt
    sep%upp   = sep%upp/evolt
    sep%udd   = sep%udd/evolt
    sep%alp   = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss   = sep%gss/evolt
    sep%gsp   = sep%gsp/evolt
    sep%gpp   = sep%gpp/evolt
    sep%gp2   = sep%gp2/evolt
    sep%gsd   = sep%gsd/evolt
    sep%gpd   = sep%gpd/evolt
    sep%gdd   = sep%gdd/evolt
    sep%hsp   = sep%hsp/evolt
    sep%fn1   = sep%fn1*bohr/evolt
    sep%fn2   = sep%fn2/bohr/bohr
    sep%fn3   = sep%fn3*bohr
    sep%bfn1  = sep%bfn1*bohr/evolt
    sep%bfn2  = sep%bfn2/bohr/bohr
    sep%bfn3  = sep%bfn3*bohr
    sep%f0sd  = sep%f0sd
    sep%g2sd  = sep%g2sd
    sep%a     = sep%a*bohr/evolt
    sep%b     = sep%b/bohr/bohr
    sep%c     = sep%c*bohr
    sep%pre   = sep%pre/evolt
    sep%d     = sep%d/bohr

  END SUBROUTINE convert_param_to_cp2k

! *****************************************************************************
!> \brief  Calculates missing parameters
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE calpar(z,sep,error)
    INTEGER                                  :: z
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calpar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iod, iop, ios, j, jmax, k, l
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: ad, am, aq, d1, d2, d3, dd, df, eisol, gdd1, gp2, gp2c, &
      gpp, gppc, gqq, gsp, gspc, gss, gssc, hpp, hpp1, hpp2, hsp, hsp1, hsp2, &
      hspc, p, p4, q1, q2, q3, qf, qn, qq, udd, upp, uss, zp, zs

    failure = .FALSE.
    IF (.NOT.sep%defined) RETURN
    uss = sep%uss
    upp = sep%upp
    udd = sep%udd
    gss = sep%gss
    gpp = sep%gpp
    gsp = sep%gsp
    gp2 = sep%gp2
    hsp = sep%hsp
    zs  = sep%sto_exponents(0)
    zp  = sep%sto_exponents(1)
    ios = Nos(z)
    iop = Nop(z)
    iod = Nod(z)

    p  = 2.0_dp
    p4 = p**4
    !  GSSC is the number of two-electron terms of type <SS|SS>
    gssc  = REAL(MAX(ios-1,0),KIND=dp)
    k = iop
    !  GSPC is the number of two-electron terms of type <SS|PP>
    gspc = REAL(ios*k,KIND=dp)
    l = MIN(k,6-k)
    !  GP2C is the number of two-electron terms of type <PP|PP>
    !       plus 0.5 of the number of HPP integrals.
    !  (HPP is not used; instead it is replaced by 0.5(GPP-GP2))
    gp2c = REAL((k*(k - 1))/2,KIND=dp) + 0.5_dp*REAL((l*(l - 1))/2,KIND=dp)
    !  GPPC is minus 0.5 times the number of HPP integrals.
    gppc = -0.5_dp*REAL((l*(l - 1))/2,KIND=dp)
    !  HSPC is the number of two-electron terms of type <SP|SP>.
    !       (S and P must have the same spin.  In all cases, if
    !  P is non-zero, there are two S electrons)
    hspc = REAL(-k ,KIND=dp)
    !  Constraint the value of the STO exponent
    zp = MAX(0.3_dp,zp)
    !  Take into account constraints on the values of the integrals
    hpp = 0.5_dp*(gpp-gp2)
    hpp = MAX(0.1_dp,hpp)
    hsp = MAX(1.E-7_dp,hsp)

    ! Evaluation of EISOL
    eisol = uss*ios + upp*iop + udd*iod + gss*gssc + gpp*gppc + gsp*gspc + gp2*gp2c + hsp*hspc

    ! Principal quantum number
    qn = REAL(nqs(z),KIND=dp)
    CPPostcondition(qn>0, cp_failure_level, routineP, error, failure)

    ! Charge separation evaluation
    dd = (2.0_dp*qn + 1)*(4.0_dp*zs*zp)**(qn + 0.5_dp)/(zs+zp)**(2.0_dp*qn + 2)/SQRT(3.0_dp)
    qq = SQRT((4.0_dp*qn*qn + 6.0_dp*qn + 2.0_dp)/20.0_dp)/zp

    ! Calculation of the additive terms in atomic units
    jmax = 5
    gdd1 = (hsp/(evolt*dd**2))**(1.0_dp/3.0_dp)
    d1 = gdd1
    d2 = gdd1 + 0.04_dp
    DO j = 1, jmax
       df = d2 - d1
       hsp1 = 0.5_dp*d1 - 0.5_dp/SQRT(4.0_dp*dd**2+1.0_dp/d1**2)
       hsp2 = 0.5_dp*d2 - 0.5_dp/SQRT(4.0_dp*dd**2+1.0_dp/d2**2)
       IF (ABS(hsp2 - hsp1) < EPSILON(0.0_dp)) EXIT
       d3 = d1 + df*(hsp/evolt-hsp1)/(hsp2 - hsp1)
       d1 = d2
       d2 = d3
    END DO
    gqq = (p4*hpp/(evolt*48.0_dp*qq**4))**0.2_dp
    q1 = gqq
    q2 = gqq + 0.04_dp
    DO j = 1, jmax
       qf = q2 - q1
       hpp1 = 0.25_dp*q1 - 0.5_dp/SQRT(4.0_dp*qq**2+1.0_dp/q1**2) + 0.25_dp/SQRT(8.0_dp*qq**2+1.0_dp/q1**2)
       hpp2 = 0.25_dp*q2 - 0.5_dp/SQRT(4.0_dp*qq**2+1.0_dp/q2**2) + 0.25_dp/SQRT(8.0_dp*qq**2+1.0_dp/q2**2)
       IF (ABS(hpp2 - hpp1) < EPSILON(0.0_dp)) EXIT
       q3 = q1 + qf*(hpp/evolt - hpp1)/(hpp2 - hpp1)
       q1 = q2
       q2 = q3
    END DO
    am = gss/evolt
    ad = d2
    aq = q2
    IF (z==1) THEN
       ad = am
       aq = am
       dd = 0.0_dp
       qq = 0.0_dp
    END IF
    ! Overwrite these parameters if they were undefined.. otherwise keep the defined
    ! value
    IF (ABS(sep%eisol) < EPSILON(0.0_dp)) sep%eisol = eisol
    IF (ABS(sep%dd)    < EPSILON(0.0_dp)) sep%dd = dd
    IF (ABS(sep%qq)    < EPSILON(0.0_dp)) sep%qq = qq
    IF (ABS(sep%am)    < EPSILON(0.0_dp)) sep%am = am
    IF (ABS(sep%ad)    < EPSILON(0.0_dp)) sep%ad = ad
    IF (ABS(sep%aq)    < EPSILON(0.0_dp)) sep%aq = aq
    ! Proceed with d-orbitals and fill the Kolpman-Ohno and Charge Separation
    ! arrays
    CALL calpar_d(sep,error)
  END SUBROUTINE calpar

! *****************************************************************************
!> \brief  Finalize the initialization of parameters, defining additional
!>         parameters for d-orbitals
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE calpar_d(sep,error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calpar_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(6)              :: amn

    failure = .FALSE.
    ! Determine if this element owns d-orbitals (only if the parametrization
    ! supports the d-orbitals)
    IF (sep%extended_basis_set) sep%dorb = element_has_d(sep,error)
    IF (sep%dorb) THEN
       CALL amn_l(sep, amn, error)
       CALL eval_1c_2el_spd(sep, error)
       CALL eval_cs_ko(sep, amn, error)
    END IF
    IF (.NOT.sep%dorb) THEN
       ! Use the old integral module
       IF (ABS(sep%am) > EPSILON(0.0_dp)) THEN
          sep%ko(1) = 0.5_dp/sep%am
       END IF
       IF (ABS(sep%ad) > EPSILON(0.0_dp).AND.(sep%z/=1)) THEN
          sep%ko(2) = 0.5_dp/sep%ad
       END IF
       IF (ABS(sep%aq) > EPSILON(0.0_dp).AND.(sep%z/=1)) THEN
          sep%ko(3) = 0.5_dp/sep%aq
       END IF
       sep%ko(7) = sep%ko(1)
       sep%ko(9) = sep%ko(1)
       sep%cs(2) = sep%dd
       sep%cs(3) = sep%qq * SQRT(2.0_dp)
    ELSE
       ! Use the new integral module
       sep%ko(9) = sep%ko(1)
       sep%aq    = 0.5_dp/sep%ko(3)
    END IF
    ! In case the Klopman-Ohno CORE therm is provided let's overwrite the
    ! computed one
    IF ( ABS(sep%rho) > EPSILON(0.0_dp) ) THEN
       sep%ko(9) = sep%rho
    END IF
  END SUBROUTINE calpar_d

! *****************************************************************************
!> \brief  Determines if the elements has d-orbitals
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   05.2008 [tlaino]
! *****************************************************************************
  FUNCTION element_has_d(sep, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'element_has_d', &
      routineP = moduleN//':'//routineN

    res = (nqd(sep%z)>0.AND.sep%sto_exponents(2)>EPSILON(0.0_dp))
  END FUNCTION element_has_d

! *****************************************************************************
!> \brief  Determines if the elements has f-orbitals
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   05.2008 [tlaino]
! *****************************************************************************
  FUNCTION element_has_f(sep, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'element_has_f', &
      routineP = moduleN//':'//routineN

    res = (nqf(sep%z)>0.AND.sep%sto_exponents(3)>EPSILON(0.0_dp))
  END FUNCTION element_has_f

! *****************************************************************************
!> \brief  Computes the A^{\mu \nu}_l values for the evaluation of the two-center
!>          two-electron integrals. The term is the one reported in Eq.(7) of TCA
!>
!> \par    Notation Index: 1 (SS), 2 (SP), 3 (SD), 4 (PP), 5 (PD), 6 (DD)
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE amn_l1 (sep, amn, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(6), INTENT(OUT) :: amn
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'amn_l1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nd, nsp
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: z1, z2, z3

    failure = .FALSE.
    z1 = sep%sto_exponents(0)
    z2 = sep%sto_exponents(1)
    z3 = sep%sto_exponents(2)
    CALL cp_assert(z1>0.0_dp,cp_fatal_level,cp_assertion_failed,routineP,&
         "Trying to use s-orbitals, but the STO exponents is set to 0. "//&
         "Please check if your parameterization supports the usage of s orbitals! "//&
CPSourceFileRef)
    amn = 0.0_dp
    nsp = nqs(sep%z)
    IF (sep%natorb>=4) THEN
       CALL cp_assert(z2>0.0_dp,cp_fatal_level,cp_assertion_failed,routineP,&
            "Trying to use p-orbitals, but the STO exponents is set to 0. "//&
            "Please check if your parameterization supports the usage of p orbitals! "//&
CPSourceFileRef)
       amn(2) = amn_l_low(z1, z2, nsp, nsp, 1)
       amn(3) = amn_l_low(z2, z2, nsp, nsp, 2)
       IF (sep%dorb) THEN
          CALL cp_assert(z3>0.0_dp,cp_fatal_level,cp_assertion_failed,routineP,&
               "Trying to use d-orbitals, but the STO exponents is set to 0. "//&
               "Please check if your parameterization supports the usage of d orbitals! "//&
CPSourceFileRef)
          nd     = nqd(sep%z)
          amn(4) = amn_l_low(z1, z3, nsp, nd, 2)
          amn(5) = amn_l_low(z2, z3, nsp, nd, 1)
          amn(6) = amn_l_low(z3, z3,  nd, nd, 2)
       END IF
    END IF
  END SUBROUTINE amn_l1

! *****************************************************************************
!> \brief  Computes the A^{\mu \nu}_l values for the evaluation of the two-center
!>          two-electron integrals. The term is the one reported in Eq.(7) of TCA
!>
!> \par
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   09.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE amn_l2 (sep, amn, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(6, 0:2), &
      INTENT(OUT)                            :: amn
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'amn_l2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nd, nsp
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: z1, z2, z3

    failure = .FALSE.
    z1 = sep%sto_exponents(0)
    z2 = sep%sto_exponents(1)
    z3 = sep%sto_exponents(2)
    CPPostcondition(z1>0.0_dp, cp_failure_level, routineP, error, failure)
    amn = 0.0_dp
    nsp = nqs(sep%z)
    amn(1,0) = amn_l_low(z1, z1, nsp, nsp, 0)
    IF (sep%natorb>=4) THEN
       CPPostcondition(z2>0.0_dp, cp_failure_level, routineP, error, failure)
       amn(2,1) = amn_l_low(z1, z2, nsp, nsp, 1)
       amn(3,0) = amn_l_low(z2, z2, nsp, nsp, 0)
       amn(3,2) = amn_l_low(z2, z2, nsp, nsp, 2)
       IF (sep%dorb) THEN
          CPPostcondition(z3>0.0_dp, cp_failure_level, routineP, error, failure)
          nd     = nqd(sep%z)
          amn(4,2) = amn_l_low(z1, z3, nsp, nd, 2)
          amn(5,1) = amn_l_low(z2, z3, nsp, nd, 1)
          amn(6,0) = amn_l_low(z3, z3,  nd, nd, 0)
          amn(6,2) = amn_l_low(z3, z3,  nd, nd, 2)
       END IF
    END IF
  END SUBROUTINE amn_l2

! *****************************************************************************
!> \brief  Low level for computing Eq.(7) of TCA
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  FUNCTION amn_l_low(z1, z2, n1, n2, l) RESULT(amnl)
    REAL(KIND=dp), INTENT(IN)                :: z1, z2
    INTEGER, INTENT(IN)                      :: n1, n2, l
    REAL(KIND=dp)                            :: amnl

    amnl = fac(n1+n2+l)/SQRT(fac(2*n1)*fac(2*n2)) * (2.0_dp*z1/(z1+z2))**n1 * &
         (2.0_dp*z2/(z1+z2))**n2 * 2.0_dp * SQRT(z1*z2) / (z1+z2)**(l+1)

  END FUNCTION amn_l_low


! *****************************************************************************
!> \brief  Calculation of chare separations and additive terms used for computing
!>         the two-center two-electron integrals with d-orbitals
!> \par    Notation
!>         -) Charge separations [sep%cs(1:6)]  [see equations (12)-(16) of TCA]
!>         -) Additive terms of Klopman-Ohno terms [sep%ko(1:9)] [see equations
!>            (19)-(26) of TCA]
!>         -) Atomic core additive term stored in sep%ko(9): used in the calculation
!>            of the core-electron attractions and core-core repulsions
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE eval_cs_ko (sep,amn,error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(6), INTENT(IN)  :: amn
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_cs_ko', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: d, fg

    failure = .FALSE.
    ! SS term
    fg = sep%gss
    sep%ko(1) = ko_ij (0, 1.0_dp, fg, error)
    IF (sep%natorb>= 4) THEN
       ! Other terms for SP basis
       ! SP
       d  = amn(2)/SQRT(3.0_dp)
       fg = sep%hsp
       sep%cs(2) = d
       sep%ko(2) = ko_ij(1, d, fg, error)
       ! PP
       sep%ko(7) = sep%ko(1)
       d  = SQRT(amn(3)*2.0_dp/5.0_dp)
       fg = 0.5_dp*(sep%gpp-sep%gp2)
       sep%cs(3) = d
       sep%ko(3) = ko_ij(2, d, fg, error)
       ! Terms involving d-orbitals
       IF (sep%dorb) THEN
          ! SD
          d  = SQRT(amn(4)*2.0_dp/SQRT(15.0_dp))
          fg = sep%onec2el(19)
          sep%cs(4) = d
          sep%ko(4) = ko_ij(2, d, fg, error)
          ! PD
          d  = amn(5)/SQRT(5.0_dp)
          fg = sep%onec2el(23) - 1.8_dp * sep%onec2el(35)
          sep%cs(5) = d
          sep%ko(5) = ko_ij(1, d, fg, error)
          ! DD
          fg = 0.2_dp * (sep%onec2el(29)+2.0_dp*sep%onec2el(30)+2.0_dp*sep%onec2el(31))
          sep%ko(8) = ko_ij (0, 1.0_dp, fg, error)
          d  = SQRT(amn(6)*2.0_dp/7.0_dp)
          fg = sep%onec2el(44) - (20.0_dp/35.0_dp) * sep%onec2el(52)
          sep%cs(6) = d
          sep%ko(6) = ko_ij(2, d, fg, error)
       END IF
    END IF
  END SUBROUTINE eval_cs_ko

! *****************************************************************************
!> \brief  Computes the 1 center two-electrons integrals for a SPD basis
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE eval_1c_2el_spd(sep,error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_1c_2el_spd', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: r016, r036, r066, r125, r155, &
                                                r234, r236, r244, r246, r266, &
                                                r355, r466, s15, s3, s5

    failure = .FALSE.
    IF (sep%dorb) THEN
       s3  = SQRT(3.0_dp)
       s5  = SQRT(5.0_dp)
       s15 = SQRT(15.0_dp)

       ! We evaluate now the Slater-Condon parameters (Rlij)
       CALL sc_param(sep, r066, r266, r466, r016, r244, r036, r236, r155, r355, r125,&
            r234, r246, error)

       IF (ABS(sep%f0sd) > EPSILON(0.0_dp)) THEN
          r016 = sep%f0sd
       END IF
       IF (ABS(sep%g2sd) > EPSILON(0.0_dp)) THEN
          r244 = sep%g2sd
       END IF
       CALL eisol_corr(sep, r016, r066, r244, r266, r466, error)
       sep%onec2el(1)  =  r016
       sep%onec2el(2)  =  2.0_dp/(3.0_dp*s5)*r125
       sep%onec2el(3)  =  1.0_dp/s15*r125
       sep%onec2el(4)  =  2.0_dp/(5.0_dp*s5)*r234
       sep%onec2el(5)  =  r036 + 4.0_dp/35.0_dp*r236
       sep%onec2el(6)  =  r036 + 2.0_dp/35.0_dp*r236
       sep%onec2el(7)  =  r036 - 4.0_dp/35.0_dp*r236
       sep%onec2el(8)  = -1.0_dp/(3.0_dp*s5)*r125
       sep%onec2el(9)  =  SQRT(3.0_dp/125.0_dp)*r234
       sep%onec2el(10) =  s3/35.0_dp*r236
       sep%onec2el(11) =  3.0_dp/35.0_dp*r236
       sep%onec2el(12) = -1.0_dp/(5.0_dp*s5)*r234
       sep%onec2el(13) =  r036 - 2.0_dp/35.0_dp*r236
       sep%onec2el(14) = -2.0_dp*s3/35.0_dp*r236
       sep%onec2el(15) = -sep%onec2el(3)
       sep%onec2el(16) = -sep%onec2el(11)
       sep%onec2el(17) = -sep%onec2el(9)
       sep%onec2el(18) = -sep%onec2el(14)
       sep%onec2el(19) =  1.0_dp/5.0_dp*r244
       sep%onec2el(20) =  2.0_dp/(7.0_dp*s5)*r246
       sep%onec2el(21) =  sep%onec2el(20)/2.0_dp
       sep%onec2el(22) = -sep%onec2el(20)
       sep%onec2el(23) =  4.0_dp/15.0_dp*r155 + 27.0_dp/245.0_dp*r355
       sep%onec2el(24) =  2.0_dp*s3/15.0_dp*r155 - 9.0_dp*s3/245.0_dp*r355
       sep%onec2el(25) =  1.0_dp/15.0_dp*r155 + 18.0_dp/245.0_dp*r355
       sep%onec2el(26) = -s3/15.0_dp*r155 + 12.0_dp*s3/245.0_dp*r355
       sep%onec2el(27) = -s3/15.0_dp*r155 -  3.0_dp*s3/245.0_dp*r355
       sep%onec2el(28) = -sep%onec2el(27)
       sep%onec2el(29) =  r066 + 4.0_dp/49.0_dp* r266 +  4.0_dp/ 49.0_dp*r466
       sep%onec2el(30) =  r066 + 2.0_dp/49.0_dp* r266 - 24.0_dp/441.0_dp*r466
       sep%onec2el(31) =  r066 - 4.0_dp/49.0_dp* r266 +  6.0_dp/441.0_dp*r466
       sep%onec2el(32) =  SQRT(3.0_dp/245.0_dp)*r246
       sep%onec2el(33) =  1.0_dp/ 5.0_dp*r155 + 24.0_dp/245.0_dp*r355
       sep%onec2el(34) =  1.0_dp/ 5.0_dp*r155 -  6.0_dp/245.0_dp*r355
       sep%onec2el(35) =  3.0_dp/49.0_dp*r355
       sep%onec2el(36) =  1.0_dp/49.0_dp*r266 + 30.0_dp/441.0_dp*r466
       sep%onec2el(37) =  s3/49.0_dp*r266   - 5.0_dp*s3/441.0_dp*r466
       sep%onec2el(38) =  r066 - 2.0_dp/49.0_dp*r266 - 4.0_dp/441.0_dp*r466
       sep%onec2el(39) = -2.0_dp*s3/49.0_dp*r266 + 10.0_dp*s3/441.0_dp*r466
       sep%onec2el(40) = -sep%onec2el(32)
       sep%onec2el(41) = -sep%onec2el(34)
       sep%onec2el(42) = -sep%onec2el(35)
       sep%onec2el(43) = -sep%onec2el(37)
       sep%onec2el(44) =  3.0_dp/49.0_dp*r266 + 20.0_dp/441.0_dp*r466
       sep%onec2el(45) = -sep%onec2el(39)
       sep%onec2el(46) =  1.0_dp/5.0_dp*r155 - 3.0_dp/35.0_dp*r355
       sep%onec2el(47) = -sep%onec2el(46)
       sep%onec2el(48) =  4.0_dp/49.0_dp*r266 + 15.0_dp/441.0_dp*r466
       sep%onec2el(49) =  3.0_dp/49.0_dp*r266 -  5.0_dp/147.0_dp*r466
       sep%onec2el(50) = -sep%onec2el(49)
       sep%onec2el(51) =  r066 + 4.0_dp/49.0_dp*r266 - 34.0_dp/441.0_dp*r466
       sep%onec2el(52) =  35.0_dp/441.0_dp*r466
       sep%f0dd     =  r066
       sep%f2dd     =  r266
       sep%f4dd     =  r466
       sep%f0sd     =  r016
       sep%g2sd     =  r244
       sep%f0pd     =  r036
       sep%f2pd     =  r236
       sep%g1pd     =  r155
       sep%g3pd     =  r355
    END IF
  END SUBROUTINE eval_1c_2el_spd

! *****************************************************************************
!> \brief  Slater-Condon parameters for 1 center 2 electrons integrals
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE sc_param (sep, r066, r266, r466, r016, r244, r036, r236, r155, r355, &
       r125, r234, r246, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), INTENT(out)               :: r066, r266, r466, r016, r244, &
                                                r036, r236, r155, r355, r125, &
                                                r234, r246
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nd, ns
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ed, ep, es

    failure = .FALSE.
    ns = nqs(sep%z)
    nd = nqd(sep%z)
    es = sep%zn(0)
    ep = sep%zn(1)
    ed = sep%zn(2)
    r016 = sc_param_low (0, ns, es, ns, es, nd, ed, nd, ed, error)
    r036 = sc_param_low (0, ns, ep, ns, ep, nd, ed, nd, ed, error)
    r066 = sc_param_low (0, nd, ed, nd, ed, nd, ed, nd, ed, error)
    r155 = sc_param_low (1, ns, ep, nd, ed, ns, ep, nd, ed, error)
    r125 = sc_param_low (1, ns, es, ns, ep, ns, ep, nd, ed, error)
    r244 = sc_param_low (2, ns, es, nd, ed, ns, es, nd, ed, error)
    r236 = sc_param_low (2, ns, ep, ns, ep, nd, ed, nd, ed, error)
    r266 = sc_param_low (2, nd, ed, nd, ed, nd, ed, nd, ed, error)
    r234 = sc_param_low (2, ns, ep, ns, ep, ns, es, nd, ed, error)
    r246 = sc_param_low (2, ns, es, nd, ed, nd, ed, nd, ed, error)
    r355 = sc_param_low (3, ns, ep, nd, ed, ns, ep, nd, ed, error)
    r466 = sc_param_low (4, nd, ed, nd, ed, nd, ed, nd, ed, error)
  END SUBROUTINE sc_param

! *****************************************************************************
!> \brief  Slater-Condon parameters for 1 center 2 electrons integrals - Low level
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \par    Notation
!>         -) k:      Type of integral
!>         -) na,na:  Principle Quantum Number of AO,corresponding to electron 1
!>         -) ea,eb:  Exponents of AO,corresponding to electron 1
!>         -) nb,nc:  Principle Quantum Number of AO,corresponding to electron 2
!>         -) ec,ed:  Exponents of AO,corresponding to electron 2
!> \date   03.2008 [tlaino]
! *****************************************************************************
  FUNCTION sc_param_low (k, na, ea, nb, eb, nc, ec, nd, ed, error) RESULT(res)
    INTEGER, INTENT(in)                      :: k, na
    REAL(KIND=dp), INTENT(in)                :: ea
    INTEGER, INTENT(in)                      :: nb
    REAL(KIND=dp), INTENT(in)                :: eb
    INTEGER, INTENT(in)                      :: nc
    REAL(KIND=dp), INTENT(in)                :: ec
    INTEGER, INTENT(in)                      :: nd
    REAL(KIND=dp), INTENT(in)                :: ed
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'sc_param_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, m, m1, m2, n, nab, ncd
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a2, aab, acd, ae, aea, aeb, &
                                                aec, aed, c, e, eab, ecd, ff, &
                                                s0, s1, s2, s3, tmp

    failure = .FALSE.
    CPPostcondition(ea>0.0_dp, cp_failure_level, routineP, error, failure)
    CPPostcondition(eb>0.0_dp, cp_failure_level, routineP, error, failure)
    CPPostcondition(ec>0.0_dp, cp_failure_level, routineP, error, failure)
    CPPostcondition(ed>0.0_dp, cp_failure_level, routineP, error, failure)
    aea = LOG(ea)
    aeb = LOG(eb)
    aec = LOG(ec)
    aed = LOG(ed)
    nab = na + nb
    ncd = nc + nd
    ecd = ec + ed
    eab = ea + eb
    e   = ecd + eab
    n   = nab + ncd
    ae  = LOG(e)
    a2  = LOG(2.0_dp)
    acd = LOG(ecd)
    aab = LOG(eab)
    ff  = fac(n-1)/SQRT(fac(2*na)*fac(2*nb)*fac(2*nc)*fac(2*nd))
    tmp = na*aea + nb*aeb + nc*aec + nd*aed + 0.5_dp*(aea+aeb+aec+aed) + a2*(n+2) - ae*n
    c   = evolt * ff * EXP(tmp)
    s0  = 1.0_dp/e
    s1  = 0.0_dp
    s2  = 0.0_dp
    m = ncd - k
    DO i = 1, m
       s0 = s0*e/ecd
       s1 = s1 + s0 * (binomial(ncd-k-1, i-1)-binomial(ncd+k, i-1))/binomial(n-1, i-1)
    END DO
    m1 = m
    m2 = ncd + k
    DO i = m1, m2
       s0 = s0 * e / ecd
       s2 = s2 + s0 * binomial(m2, i)/binomial(n-1, i)
    END DO
    s3 = EXP(ae*n-acd*(m2+1)-aab*(nab-k))/binomial(n-1, m2)
    res = c * (s1-s2+s3)
  END FUNCTION sc_param_low

! *****************************************************************************
!> \brief  Corrects the EISOL fo the one-center terms coming from those atoms
!>         that have partially filled "d" shells
!> \par    Notation
!>         r016:  <SS|DD>
!>         r066:  <DD|DD> "0" term
!>         r244:  <SD|SD>
!>         r266:  <DD|DD> "2" term
!>         r466:  <DD|DD> "4" term
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE eisol_corr (sep, r016, r066, r244, r266, r466, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), INTENT(in)                :: r016, r066, r244, r266, r466
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eisol_corr', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    sep%eisol = sep%eisol + ir016(sep%z) * r016          +&
                            ir066(sep%z) * r066          -&
                            ir244(sep%z) * r244/5.0_dp   -&
                            ir266(sep%z) * r266/49.0_dp  -&
                            ir466(sep%z) * r466/49.0_dp
  END SUBROUTINE eisol_corr

! *****************************************************************************
!> \brief  Computes the Klopman-Ohno additive terms for 2-center 2-electron
!>         integrals requiring that the corresponding 1-center 2-electron integral
!>         is reproduced from the 2-center one for r->0
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   03.2008 [tlaino]
! *****************************************************************************
  FUNCTION ko_ij(l, d, fg, error) RESULT(res)
    INTEGER, INTENT(in)                      :: l
    REAL(KIND=dp), INTENT(in)                :: d, fg
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ko_ij', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: niter = 100
    REAL(KIND=dp), PARAMETER                 :: epsil = 1.0E-08_dp, &
                                                g1 = 0.382_dp, g2 = 0.618_dp

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a1, a2, delta, dsq, ev4, ev8, &
                                                f1, f2, y1, y2

    failure = .FALSE.
    CPPostcondition(fg/=0.0_dp, cp_failure_level, routineP, error, failure)
    ! Term for SS
    IF (l == 0) THEN
       res = 0.5_dp * evolt/fg
       RETURN
    END IF
    ! Term for Higher angular momentum
    dsq = d * d
    ev4 = evolt * 0.25_dp
    ev8 = evolt / 8.0_dp
    a1 = 0.1_dp
    a2 = 5.0_dp
    DO i = 1, niter
       delta = a2 - a1
       IF (delta < epsil) EXIT
       y1 = a1 + delta * g1
       y2 = a1 + delta * g2
       IF      (l == 1) THEN
          f1 = (ev4*(1/y1 - 1/SQRT(y1**2 + dsq)) - fg) ** 2
          f2 = (ev4*(1/y2 - 1/SQRT(y2**2 + dsq)) - fg) ** 2
       ELSE IF (l == 2) THEN
          f1 = (ev8*(1.0_dp/y1 - 2.0_dp/SQRT(y1**2 + dsq*0.5_dp) + 1.0_dp/SQRT(y1**2 + dsq)) - fg)** 2
          f2 = (ev8*(1/y2 - 2.0_dp/SQRT(y2**2 + dsq*0.5_dp) + 1.0_dp/SQRT(y2**2 + dsq)) - fg)** 2
       END IF
       IF (f1 < f2) THEN
          a2 = y2
       ELSE
          a1 = y1
       END IF
    END DO
    ! Convergence reached.. define additive terms
    IF (f1 >= f2) THEN
       res = a2
    ELSE
       res = a1
    END IF
  END FUNCTION ko_ij

! *****************************************************************************
!> \brief  Fills the 1 center 2 electron integrals for the construction of the
!>         one-electron fock matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE setup_1c_2el_int(sep, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_1c_2el_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, ij0, ind, ip, ipx, &
                                                ipy, ipz, isize, kl, natorb, &
                                                stat
    LOGICAL                                  :: defined, failure
    REAL(KIND=dp)                            :: gp2, gpp, gsp, gss, hsp

    failure = .FALSE.
    CALL get_se_param(sep, defined=defined, natorb=natorb,&
                      gss=gss, gsp=gsp, gpp=gpp, gp2=gp2, hsp=hsp)
    CPPostcondition(defined, cp_failure_level, routineP, error, failure)

    isize     = natorb*(natorb+1)/2
    ALLOCATE(sep%w(isize,isize),stat=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! Initialize array
    sep%w = 0.0_dp
    ! Fill the array
    IF (natorb>0) THEN
       ip            = 1
       sep%w(ip, ip) = gss
       IF (natorb > 2) THEN
          ipx = ip + 2
          ipy = ip + 5
          ipz = ip + 9
          sep%w(ipx,   ip ) = gsp
          sep%w(ipy,   ip ) = gsp
          sep%w(ipz,   ip ) = gsp
          sep%w(ip,   ipx ) = gsp
          sep%w(ip,   ipy ) = gsp
          sep%w(ip,   ipz ) = gsp
          sep%w(ipx,  ipx ) = gpp
          sep%w(ipy,  ipy ) = gpp
          sep%w(ipz,  ipz ) = gpp
          sep%w(ipy,  ipx ) = gp2
          sep%w(ipz,  ipx ) = gp2
          sep%w(ipz,  ipy ) = gp2
          sep%w(ipx,  ipy ) = gp2
          sep%w(ipx,  ipz ) = gp2
          sep%w(ipy,  ipz ) = gp2
          sep%w(ip+1, ip+1) = hsp
          sep%w(ip+3, ip+3) = hsp
          sep%w(ip+6, ip+6) = hsp
          sep%w(ip+4, ip+4) = 0.5_dp*(gpp-gp2)
          sep%w(ip+7, ip+7) = 0.5_dp*(gpp-gp2)
          sep%w(ip+8, ip+8) = 0.5_dp*(gpp-gp2)
          IF (sep%dorb) THEN
             ij0 = ip - 1
             DO i = 1, 243
                ij  = int_ij(i)
                kl  = int_kl(i)
                ind = int_onec2el(i)
                sep%w(ij+ij0, kl+ij0) = sep%onec2el(ind)/evolt
             END DO
          END IF
       END IF
    END IF
  END SUBROUTINE setup_1c_2el_int

END MODULE semi_empirical_par_utils

