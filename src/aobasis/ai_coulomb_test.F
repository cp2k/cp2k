!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Test of Electron Repulsion Routines (ERI)
!> \par History
!>      none
!> \author JGH (01.07.2009)
! *****************************************************************************
MODULE ai_coulomb_test

  USE ai_coulomb,                      ONLY: coulomb2
  USE ai_eri_debug,                    ONLY: init_os,&
                                             os
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE orbital_pointers,                ONLY: coset,&
                                             deallocate_orbital_pointers,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE physcon,                         ONLY: bohr
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_coulomb_test'

  REAL(KIND=dp), PARAMETER :: threshold=1.0E-6_dp

  PRIVATE

  PUBLIC :: eri_test
! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE eri_test (iw,error)

    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eri_test', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: lmax = 6

    CHARACTER(LEN=11), DIMENSION(0:lmax)     :: i2g
    CHARACTER(LEN=5), DIMENSION(0:lmax)      :: i2c
    CHARACTER(LEN=7), DIMENSION(0:lmax)      :: i2e
    CHARACTER(LEN=9), DIMENSION(0:lmax)      :: i2f
    INTEGER                                  :: i, ii, istat, l, la_max, &
                                                la_min, lc_max, lc_min, ll, &
                                                n, na, nb, npgfa, npgfb, &
                                                npgfc, npgfd
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: perf, rac2, t, tend, tstart
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: f
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: vac
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: v, za, zb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: iab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :, :)            :: iabc
    REAL(KIND=dp), DIMENSION(3)              :: ra, rb, rc, rd
    REAL(KIND=dp), DIMENSION(6)              :: a, b
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rpgf, zeta, zetb, zetc, zetd

    failure = .FALSE.
    IF ( iw>0 ) WRITE(iw,'(/,A)') " Test of Electron Repulsion Integrals (ERI) "

    CALL init_orbital_pointers(lmax)

    i2c(0) = "(s|s)"
    i2c(1) = "(p|p)"
    i2c(2) = "(d|d)"
    i2c(3) = "(f|f)"
    i2c(4) = "(g|g)"
    i2c(5) = "(h|h)"
    i2c(6) = "(i|i)"

    i2g(0) = "[(ss)|(ss)]"
    i2g(1) = "[(pp)|(pp)]"
    i2g(2) = "[(dd)|(dd)]"
    i2g(3) = "[(ff)|(ff)]"
    i2g(4) = "[(gg)|(gg)]"
    i2g(5) = "[(hh)|(hh)]"
    i2g(6) = "[(ii)|(ii)]"

    i2f(0) = "[ss|(ss)]"
    i2f(1) = "[pp|(pp)]"
    i2f(2) = "[dd|(dd)]"
    i2f(3) = "[ff|(ff)]"
    i2f(4) = "[gg|(gg)]"
    i2f(5) = "[hh|(hh)]"
    i2f(6) = "[ii|(ii)]"

    i2e(0) = "(ss|ss)"
    i2e(1) = "(pp|pp)"
    i2e(2) = "(dd|dd)"
    i2e(3) = "(ff|ff)"
    i2e(4) = "(gg|gg)"
    i2e(5) = "(hh|hh)"
    i2e(6) = "(ii|ii)"

    npgfa = 4
    npgfb = 2
    npgfc = 4
    npgfd = 1
    n = MAX(npgfa,npgfb,npgfc,npgfd)

    ALLOCATE(zeta(npgfa),zetb(npgfb),zetc(npgfc),zetd(npgfd),rpgf(n),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    zeta(1:npgfa) = 0.5_dp
    zetb(1:npgfb) = 0.4_dp
    zetc(1:npgfc) = 0.3_dp
    zetd(1:npgfd) = 0.2_dp

    ra = (/0.0_dp,0.0_dp,0.0_dp/)
    rb = (/1.0_dp,0.0_dp,0.0_dp/)
    rc = (/0.0_dp,0.3_dp,0.3_dp/)
    rd = (/0.7_dp,0.2_dp,0.1_dp/)

    rac2 = SUM((ra-rc)**2)
    rpgf = 1.e10_dp

    ! Performance test of coulomb2 routine
    IF ( iw>0 ) THEN

       WRITE(iw,'(//,A,/)') " Test of 2-Electron-2-Center Integrals (coulomb2) "
       DO l=0,lmax
         la_max=l
         la_min=l
         lc_max=l
         lc_min=l
         ll = ncoset(l)
         ALLOCATE(f(0:2*l+2),v(npgfa*ll,npgfc*ll,2*l+1),vac(npgfa*ll,npgfc*ll),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         vac = 0._dp
         ii = MAX(100/(l+1)**2,1)
         tstart = m_walltime ( )
         DO i=1,ii
           CALL coulomb2(la_max,npgfa,zeta,rpgf,la_min,lc_max,npgfc,zetc,rpgf,lc_min,rc,rac2,vac,v,f)
         END DO
         tend = m_walltime ( )
         t = tend - tstart + threshold
         perf = REAL(ii*nco(l)**2)*1.e-6_dp * REAL(npgfa*npgfc,dp)/t
         WRITE(iw,'(A,T40,A,T66,F15.3)') " Performance [Mintegrals/s] ",i2c(l),perf
         DEALLOCATE(f,v,vac,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO

    END IF

    DEALLOCATE(zeta,zetb,zetc,zetd,rpgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL deallocate_orbital_pointers()

  END SUBROUTINE eri_test

! *****************************************************************************

END MODULE ai_coulomb_test

