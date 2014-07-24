!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of procedures to modify the density matrix, in order to collocate
!>      on the real space grid specific functions of the primitives
!> \author MI 06.2006
! *****************************************************************************
MODULE  qs_modify_pab_block

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_modify_pab_block'

! *** Public subroutines ***

  PUBLIC :: prepare_dadb,prepare_adb_m_dab, prepare_dab_p_adb, prepare_ardb_m_darb,prepare_arb
  PUBLIC :: prepare_diadib, prepare_dijadijb,prepare_diiadiib
  PUBLIC :: FUNC_AB, FUNC_ARB, FUNC_DADB, FUNC_ADBmDAB, FUNC_ARDBmDARB, FUNC_DABpADB
  PUBLIC :: FUNC_DER,FUNC_DXDY,FUNC_DYDZ,FUNC_DZDX,FUNC_DX,FUNC_DY,FUNC_DZ
  PUBLIC :: FUNC_DXDX,FUNC_DYDY,FUNC_DZDZ

  ! Flags for the type of function to be collocated
  INTEGER, PARAMETER :: FUNC_AB = 401, FUNC_DADB = 402
  INTEGER, PARAMETER :: FUNC_ADBmDAB = 403, FUNC_ARDBmDARB = 404
  INTEGER, PARAMETER :: FUNC_DABpADB = 405, FUNC_ARB = 406
  INTEGER, DIMENSION(0:9), PARAMETER :: FUNC_DER = (/401,501,502,503,601,602,603,604,605,606/)
  INTEGER, PARAMETER :: FUNC_DXDY=601,FUNC_DYDZ= 602, FUNC_DZDX = 603
  INTEGER, PARAMETER :: FUNC_DXDX=604,FUNC_DYDY= 605, FUNC_DZDZ = 606
  INTEGER, PARAMETER :: FUNC_DX=501,FUNC_DY=502, FUNC_DZ = 503

CONTAINS

!

! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>      is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
!>      (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})
!> \param pab_local ...
!> \param pab ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
!> \param zeta ...
!> \param zetb ...
! *****************************************************************************
  SUBROUTINE prepare_dadb(pab_local,pab,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    INTEGER, INTENT(IN)                      :: lxa, lya, lza, lxb, lyb, lzb, &
                                                o1, o2
    REAL(dp), INTENT(IN)                     :: zeta, zetb

    INTEGER                                  :: ico, ico_l, jco, jco_l

! this element of pab results in 12 elements of pab_local

     ico=coset(lxa,lya,lza)
     jco=coset(lxb,lyb,lzb)
     ! x  (all safe if lxa = 0, as the spurious added terms have zero prefactor)

     ico_l=coset(MAX(lxa-1,0),lya,lza)
     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lxa* lxb*pab(o1+ico,o2+jco)
     ico_l=coset(MAX(lxa-1,0),lya,lza)
     jco_l=coset(   (lxb+1  ),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lxa*zetb*pab(o1+ico,o2+jco)
     ico_l=coset(   (lxa+1  ),lya,lza)
     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lxb*pab(o1+ico,o2+jco)
     ico_l=coset(   (lxa+1  ),lya,lza)
     jco_l=coset(   (lxb+1  ),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

     ! y

     ico_l=coset(lxa,MAX(lya-1,0),lza)
     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lya* lyb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,MAX(lya-1,0),lza)
     jco_l=coset(lxb,   (lyb+1  ),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lya*zetb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,   (lya+1  ),lza)
     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lyb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,   (lya+1  ),lza)
     jco_l=coset(lxb,   (lyb+1  ),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

     ! z

     ico_l=coset(lxa,lya,MAX(lza-1,0))
     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lza* lzb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,lya,MAX(lza-1,0))
     jco_l=coset(lxb,lyb,   (lzb+1  ))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lza*zetb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,lya,   (lza+1  ))
     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lzb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,lya,   (lza+1  ))
     jco_l=coset(lxb,lyb,   (lzb+1  ))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

  END SUBROUTINE prepare_dadb

!
! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>      is equivalent to mapping pab with  (ddi pgf_a) . (ddi pgf_b)
!>      (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})
!> \param pab_local ...
!> \param pab ...
!> \param ider ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
!> \param zeta ...
!> \param zetb ...
! *****************************************************************************
  SUBROUTINE prepare_diadib(pab_local,pab,ider,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

    REAL(dp), DIMENSION(:, :), INTENT(inout) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(in)    :: pab
    INTEGER, INTENT(in)                      :: ider, lxa, lya, lza, lxb, &
                                                lyb, lzb, o1, o2
    REAL(dp), INTENT(in)                     :: zeta, zetb

    INTEGER                                  :: ico, ico_l, jco, jco_l

! this element of pab results in 12 elements of pab_local

     ico=coset(lxa,lya,lza)
     jco=coset(lxb,lyb,lzb)
     IF(ider==1) THEN
     ! x  (all safe if lxa = 0, as the spurious added terms have zero prefactor)
     ico_l=coset(MAX(lxa-1,0),lya,lza)
     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lxa* lxb*pab(o1+ico,o2+jco)
     ico_l=coset(MAX(lxa-1,0),lya,lza)
     jco_l=coset(   (lxb+1  ),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lxa*zetb*pab(o1+ico,o2+jco)
     ico_l=coset(   (lxa+1  ),lya,lza)
     jco_l=coset(MAX(lxb-1,0),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lxb*pab(o1+ico,o2+jco)
     ico_l=coset(   (lxa+1  ),lya,lza)
     jco_l=coset(   (lxb+1  ),lyb,lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)

     ELSEIF(ider==2) THEN
     ! y

     ico_l=coset(lxa,MAX(lya-1,0),lza)
     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lya* lyb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,MAX(lya-1,0),lza)
     jco_l=coset(lxb,   (lyb+1  ),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lya*zetb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,   (lya+1  ),lza)
     jco_l=coset(lxb,MAX(lyb-1,0),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lyb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,   (lya+1  ),lza)
     jco_l=coset(lxb,   (lyb+1  ),lzb)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)
     ELSEIF(ider==3) THEN
     ! z
     ico_l=coset(lxa,lya,MAX(lza-1,0))
     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lza* lzb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,lya,MAX(lza-1,0))
     jco_l=coset(lxb,lyb,   (lzb+1  ))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp* lza*zetb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,lya,   (lza+1  ))
     jco_l=coset(lxb,lyb,MAX(lzb-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta* lzb*pab(o1+ico,o2+jco)
     ico_l=coset(lxa,lya,   (lza+1  ))
     jco_l=coset(lxb,lyb,   (lzb+1  ))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zeta*zetb*pab(o1+ico,o2+jco)
     END IF

  END SUBROUTINE prepare_diadib
!
! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>      is equivalent to mapping pab with  (ddidj pgf_a) . (ddidj pgf_b)
!>      (ddxdy pgf_a ) (ddxdy pgf_b) =
!> \param pab_local ...
!> \param pab ...
!> \param ider1 ...
!> \param ider2 ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
!> \param zeta ...
!> \param zetb ...
! *****************************************************************************
  SUBROUTINE prepare_dijadijb(pab_local,pab,ider1,ider2,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

    REAL(dp), DIMENSION(:, :), INTENT(inout) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(in)    :: pab
    INTEGER, INTENT(in)                      :: ider1, ider2, lxa, lya, lza, &
                                                lxb, lyb, lzb, o1, o2
    REAL(dp), INTENT(in)                     :: zeta, zetb

    INTEGER                                  :: ico, ico_l, jco
    REAL(dp)                                 :: func_a

! this element of pab results in 12 elements of pab_local

     ico=coset(lxa,lya,lza)
     jco=coset(lxb,lyb,lzb)

     IF((ider1==1 .AND. ider2==2) .OR. (ider1==2 .AND. ider2==1)) THEN  ! xy
       ico_l=coset(MAX(lxa-1,0),MAX(lya-1,0),lza)
       func_a = lxa*lya*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=1)
       ico_l=coset(lxa+1,MAX(lya-1,0),lza)
       func_a = -2.0_dp*zeta*lya*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=1)
       ico_l=coset(MAX(lxa-1,0),lya+1,lza)
       func_a = -2.0_dp*zeta*lxa*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=1)
       ico_l=coset(lxa+1,lya+1,lza)
       func_a = 4.0_dp*zeta*zeta*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=1)
     ELSEIF((ider1==2 .AND. ider2==3) .OR. (ider1==3 .AND. ider2==2)) THEN  ! yz
       ico_l=coset(lxa,MAX(lya-1,0),MAX(lza-1,0))
       func_a = lya*lza*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=2)
       ico_l=coset(lxa,lya+1,MAX(lza-1,0))
       func_a = -2.0_dp*zeta*lza*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=2)
       ico_l=coset(lxa,MAX(lya-1,0),lza+1)
       func_a = -2.0_dp*zeta*lya*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=2)
       ico_l=coset(lxa,lya+1,lza+1)
       func_a = 4.0_dp*zeta*zeta*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=2)
     ELSEIF((ider1==3 .AND. ider2==1) .OR. (ider1==1 .AND. ider2==3)) THEN  ! zx
       ico_l=coset(MAX(lxa-1,0),lya,MAX(lza-1,0))
       func_a = lza*lxa*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=3)
       ico_l=coset(MAX(lxa-1,0),lya,lza+1)
       func_a = -2.0_dp*zeta*lxa*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=3)
       ico_l=coset(lxa+1,lya,MAX(lza-1,0))
       func_a = -2.0_dp*zeta*lza*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=3)
       ico_l=coset(lxa+1,lya,lza+1)
       func_a = 4.0_dp*zeta*zeta*pab(o1+ico,o2+jco)
       CALL oneterm_dijdij(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=3)
     END IF

  END SUBROUTINE prepare_dijadijb

! *****************************************************************************
!> \brief ...
!> \param pab_local ...
!> \param func_a ...
!> \param ico_l ...
!> \param lx ...
!> \param ly ...
!> \param lz ...
!> \param zet ...
!> \param idir ...
! *****************************************************************************
  SUBROUTINE oneterm_dijdij(pab_local,func_a,ico_l,lx,ly,lz,zet,idir)

    REAL(dp), DIMENSION(:, :), INTENT(inout) :: pab_local
    REAL(dp), INTENT(in)                     :: func_a
    INTEGER, INTENT(in)                      :: ico_l, lx, ly, lz
    REAL(dp), INTENT(in)                     :: zet
    INTEGER, INTENT(in)                      :: idir

    INTEGER                                  :: jco_l, l1, l2

    IF(idir ==1) THEN
     l1 = lx
     l2 = ly
     jco_l=coset(MAX(lx-1,0),MAX(ly-1,0),lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        l1*  l2*func_a
     jco_l=coset(lx+1,MAX(ly-1,0),lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*  l2*func_a
     jco_l=coset(MAX(lx-1,0),ly+1,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*  l1*func_a
     jco_l=coset(lx+1,ly+1,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zet*zet*func_a
    ELSEIF(idir==2) THEN
     l1 = ly
     l2 = lz
     jco_l=coset(lx,MAX(ly-1,0),MAX(lz-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        l1*  l2*func_a
     jco_l=coset(lx,ly+1,MAX(lz-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*  l2*func_a
     jco_l=coset(lx,MAX(ly-1,0),lz+1)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*  l1*func_a
     jco_l=coset(lx,ly+1,lz+1)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zet*zet*func_a
    ELSEIF(idir==3) THEN
     l1 = lz
     l2 = lx
     jco_l=coset(MAX(lx-1,0),ly,MAX(lz-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        l1*  l2*func_a
     jco_l=coset(MAX(lx-1,0),ly,lz+1)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*  l2*func_a
     jco_l=coset(lx+1,ly,MAX(lz-1,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*  l1*func_a
     jco_l=coset(lx+1,ly,lz+1)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zet*zet*func_a

    END IF
  END SUBROUTINE oneterm_dijdij

!
! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>      is equivalent to mapping pab with  (ddidj pgf_a) . (ddidj pgf_b)
!>      (ddxdx pgf_a ) (ddxdx pgf_b) =
!> \param pab_local ...
!> \param pab ...
!> \param ider ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
!> \param zeta ...
!> \param zetb ...
! *****************************************************************************
  SUBROUTINE prepare_diiadiib(pab_local,pab,ider,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

    REAL(dp), DIMENSION(:, :), INTENT(inout) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(in)    :: pab
    INTEGER, INTENT(in)                      :: ider, lxa, lya, lza, lxb, &
                                                lyb, lzb, o1, o2
    REAL(dp), INTENT(in)                     :: zeta, zetb

    INTEGER                                  :: ico, ico_l, jco
    REAL(dp)                                 :: func_a

! this element of pab results in  9 elements of pab_local

     ico=coset(lxa,lya,lza)
     jco=coset(lxb,lyb,lzb)

     IF(ider==1 ) THEN  ! x
       ico_l=coset(MAX(lxa-2,0),lya,lza)
       func_a = lxa*(lxa-1)*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=1)
       ico_l=coset(lxa,lya,lza)
       func_a = -2.0_dp*zeta*(2*lxa+1)*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=1)
       ico_l=coset(lxa+2,lya,lza)
       func_a = 4.0_dp*zeta*zeta*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=1)
     ELSEIF(ider==2 ) THEN  ! y
       ico_l=coset(lxa,MAX(lya-2,0),lza)
       func_a = lya*(lya-1)*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=2)
       ico_l=coset(lxa,lya,lza)
       func_a = -2.0_dp*zeta*(2*lya+1)*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=2)
       ico_l=coset(lxa,lya+2,lza)
       func_a = 4.0_dp*zeta*zeta*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=2)
     ELSEIF(ider==3 ) THEN  ! z
       ico_l=coset(lxa,lya,MAX(lza-2,0))
       func_a = lza*(lza-1)*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=3)
       ico_l=coset(lxa,lya,lza)
       func_a = -2.0_dp*zeta*(2*lza+1)*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=3)
       ico_l=coset(lxa,lya,lza+2)
       func_a = 4.0_dp*zeta*zeta*pab(o1+ico,o2+jco)
       CALL oneterm_diidii(pab_local,func_a,ico_l,lxb,lyb,lzb,zetb,idir=3)
     END IF

  END SUBROUTINE prepare_diiadiib


! *****************************************************************************
!> \brief ...
!> \param pab_local ...
!> \param func_a ...
!> \param ico_l ...
!> \param lx ...
!> \param ly ...
!> \param lz ...
!> \param zet ...
!> \param idir ...
! *****************************************************************************
  SUBROUTINE oneterm_diidii(pab_local,func_a,ico_l,lx,ly,lz,zet,idir)
    REAL(dp), DIMENSION(:, :), INTENT(inout) :: pab_local
    REAL(dp), INTENT(in)                     :: func_a
    INTEGER, INTENT(in)                      :: ico_l, lx, ly, lz
    REAL(dp), INTENT(in)                     :: zet
    INTEGER, INTENT(in)                      :: idir

    INTEGER                                  :: jco_l, l1

    IF(idir ==1) THEN
     l1 = lx
     jco_l=coset(MAX(lx-2,0),ly,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        l1*(l1-1)*func_a
     jco_l=coset(lx,ly,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*(2*l1+1)*func_a
     jco_l=coset(lx+2,ly,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zet*zet*func_a
    ELSEIF(idir==2) THEN
     l1 = ly
     jco_l=coset(lx,MAX(ly-2,0),lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        l1*(l1-1)*func_a
     jco_l=coset(lx,ly,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*(2*l1+1)*func_a
     jco_l=coset(lx,ly+2,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zet*zet*func_a
    ELSEIF(idir==3) THEN
     l1 = lz
     jco_l=coset(lx,ly,MAX(lz-2,0))
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        l1*(l1-1)*func_a
     jco_l=coset(lx,ly,lz)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zet*(2*l1+1)*func_a
     jco_l=coset(lx,ly,lz+2)
     pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+4.0_dp*zet*zet*func_a
    END IF

  END SUBROUTINE oneterm_diidii

! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>      is equivalent to mapping pab with  pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
!>      ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) =
!>              pgf_a *(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x}) - (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b
!> \param pab_local ...
!> \param pab ...
!> \param idir ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
!> \param zeta ...
!> \param zetb ...
! *****************************************************************************
  SUBROUTINE prepare_adb_m_dab(pab_local,pab,idir,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    INTEGER, INTENT(IN)                      :: idir, lxa, lya, lza, lxb, &
                                                lyb, lzb, o1, o2
    REAL(dp), INTENT(IN)                     :: zeta, zetb

    INTEGER                                  :: ico, ico_l, jco, jco_l

! this element of pab results in 4 elements of pab_local

     ico=coset(lxa,lya,lza)
     jco=coset(lxb,lyb,lzb)

     IF(idir==1) THEN  ! x
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(MAX(lxb-1,0),lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lxb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,lza)
        jco_l=coset((lxb+1),lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zetb*pab(o1+ico,o2+jco)
        ico_l=coset(MAX(lxa-1,0),lya,lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-        lxa*pab(o1+ico,o2+jco)
        ico_l=coset((lxa+1),lya,lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+2.0_dp*zeta*pab(o1+ico,o2+jco)
     ELSEIF(idir==2) THEN ! y
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,MAX(lyb-1,0),lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lyb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,(lyb+1),lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zetb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,MAX(lya-1,0),lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-        lya*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,(lya+1),lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+2.0_dp*zeta*pab(o1+ico,o2+jco)
     ELSE  ! z
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,lyb,MAX(lzb-1,0))
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lzb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,lyb,   (lzb+1  ))
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zetb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,MAX(lza-1,0))
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-        lza*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,   (lza+1  ))
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+2.0_dp*zeta*pab(o1+ico,o2+jco)
     END IF

  END SUBROUTINE prepare_adb_m_dab

!
! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>      is equivalent to mapping pab with  pgf_a (nabla_{idir} pgf_b) + (nabla_{idir} pgf_a) pgf_b
!>      ( pgf_a ) (ddx pgf_b) + (ddx pgf_a)( pgf_b ) =
!>              pgf_a *(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x}) + (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b
!> \param pab_local ...
!> \param pab ...
!> \param idir ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
!> \param zeta ...
!> \param zetb ...
! *****************************************************************************
  SUBROUTINE prepare_dab_p_adb(pab_local,pab,idir,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    INTEGER, INTENT(IN)                      :: idir, lxa, lya, lza, lxb, &
                                                lyb, lzb, o1, o2
    REAL(dp), INTENT(IN)                     :: zeta, zetb

    INTEGER                                  :: ico, ico_l, jco, jco_l

! this element of pab results in 4 elements of pab_local

     ico=coset(lxa,lya,lza)
     jco=coset(lxb,lyb,lzb)

     IF(idir==1) THEN  ! x
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(MAX(lxb-1,0),lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lxb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,lza)
        jco_l=coset((lxb+1),lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zetb*pab(o1+ico,o2+jco)
        ico_l=coset(MAX(lxa-1,0),lya,lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lxa*pab(o1+ico,o2+jco)
        ico_l=coset((lxa+1),lya,lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta*pab(o1+ico,o2+jco)
     ELSEIF(idir==2) THEN ! y
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,MAX(lyb-1,0),lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lyb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,(lyb+1),lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zetb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,MAX(lya-1,0),lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lya*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,(lya+1),lza)
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta*pab(o1+ico,o2+jco)
     ELSE  ! z
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,lyb,MAX(lzb-1,0))
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lzb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,lza)
        jco_l=coset(lxb,lyb,   (lzb+1  ))
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zetb*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,MAX(lza-1,0))
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)+        lza*pab(o1+ico,o2+jco)
        ico_l=coset(lxa,lya,   (lza+1  ))
        jco_l=coset(lxb,lyb,lzb)
        pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l)-2.0_dp*zeta*pab(o1+ico,o2+jco)
     END IF

  END SUBROUTINE prepare_dab_p_adb

!

! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>      pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}  pgf_b
!>       ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b ) =
!>                     pgf_a *(lbx pgf_{b-1x+1ir} - 2*zetb*pgf_{b+1x+1ir}) -
!>                              (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}
!> \param pab_local ...
!> \param pab ...
!> \param idir ...
!> \param ir ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
!> \param zeta ...
!> \param zetb ...
! *****************************************************************************
  SUBROUTINE prepare_ardb_m_darb(pab_local,pab,idir,ir,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    INTEGER, INTENT(IN)                      :: idir, ir, lxa, lya, lza, lxb, &
                                                lyb, lzb, o1, o2
    REAL(dp), INTENT(IN)                     :: zeta, zetb

    INTEGER                                  :: ico, ico_l, jco, jco_l

! this element of pab results in 4 elements of pab_local

     ico=coset(lxa,lya,lza)
     jco=coset(lxb,lyb,lzb)

     IF(idir==1 .AND. ir==1) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lxb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb+2),lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(MAX(lxa-1,0),lya,lza)
       jco_l=coset((lxb+1),lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lxa*pab(o1+ico,o2+jco)

       ico_l=coset((lxa+1),lya,lza)
       jco_l=coset((lxb+1),lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==1 .AND. ir==2) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(MAX(lxb-1,0),(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lxb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb+1),(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(MAX(lxa-1,0),lya,lza)
       jco_l=coset(lxb,(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lxa*pab(o1+ico,o2+jco)

       ico_l=coset((lxa+1),lya,lza)
       jco_l=coset(lxb,(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==1 .AND. ir==3) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(MAX(lxb-1,0),lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lxb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb+1),lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(MAX(lxa-1,0),lya,lza)
       jco_l=coset(lxb,lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lxa*pab(o1+ico,o2+jco)

       ico_l=coset((lxa+1),lya,lza)
       jco_l=coset(lxb,lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==2 .AND. ir==1) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb+1),MAX(lyb-1,0),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lyb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb+1),(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,MAX(lya-1,0),lza)
       jco_l=coset((lxb+1),lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lya*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,(lya+1),lza)
       jco_l=coset((lxb+1),lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==2 .AND. ir==2) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lyb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,(lyb+2),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,MAX(lya-1,0),lza)
       jco_l=coset(lxb,(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lya*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,(lya+1),lza)
       jco_l=coset(lxb,(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==2 .AND. ir==3) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,MAX(lyb-1,0),(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lyb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,(lyb+1),(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,MAX(lya-1,0),lza)
       jco_l=coset(lxb,lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lya*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,(lya+1),lza)
       jco_l=coset(lxb,lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==3 .AND. ir==1) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb+1),lyb,MAX(lzb-1,0))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lzb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb+1),lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,MAX(lza-1,0))
       jco_l=coset((lxb+1),lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lza*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,(lza+1))
       jco_l=coset((lxb+1),lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==3 .AND. ir==2) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,(lyb+1),MAX(lzb-1,0))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lzb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,(lyb+1),(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,MAX(lza-1,0))
       jco_l=coset(lxb,(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lza*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,(lza+1))
       jco_l=coset(lxb,(lyb+1),lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     ELSEIF(idir==3 .AND. ir==3) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,lyb,lzb)
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) +         lzb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,lyb,(lzb+2))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) - 2.0_dp*zetb*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,MAX(lza-1,0))
       jco_l=coset(lxb,lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) -         lza*pab(o1+ico,o2+jco)

       ico_l=coset(lxa,lya,(lza+1))
       jco_l=coset(lxb,lyb,(lzb+1))
       pab_local(ico_l,jco_l)=pab_local(ico_l,jco_l) + 2.0_dp*zeta*pab(o1+ico,o2+jco)

     END IF

  END SUBROUTINE prepare_ardb_m_darb

! *****************************************************************************
!> \brief create a new pab_local so that mapping pab_local with pgf_a pgf_b
!>    pgf_a (r-Rb)_{ir} pgf_b = pgf_a * pgf_b{+1ir}
!> \param pab_local ...
!> \param pab ...
!> \param ir ...
!> \param lxa ...
!> \param lya ...
!> \param lza ...
!> \param lxb ...
!> \param lyb ...
!> \param lzb ...
!> \param o1 ...
!> \param o2 ...
! *****************************************************************************
  SUBROUTINE prepare_arb(pab_local,pab,ir,lxa,lya,lza,lxb,lyb,lzb,o1,o2)

    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: pab_local
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pab
    INTEGER, INTENT(IN)                      :: ir, lxa, lya, lza, lxb, lyb, &
                                                lzb, o1, o2

    INTEGER                                  :: ico, ico_l, jco, jco_l

    ico=coset(lxa,lya,lza)
    jco=coset(lxb,lyb,lzb)

    IF( ir == 1 ) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset((lxb + 1),lyb,lzb)
       pab_local(ico_l,jco_l) = pab_local(ico_l,jco_l) + pab(o1+ico,o2+jco)

    ELSEIF( ir == 2 ) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,(lyb + 1),lzb)
       pab_local(ico_l,jco_l) = pab_local(ico_l,jco_l) + pab(o1+ico,o2+jco)

    ELSEIF( ir == 3 ) THEN

       ico_l=coset(lxa,lya,lza)
       jco_l=coset(lxb,lyb,(lzb + 1))
       pab_local(ico_l,jco_l) = pab_local(ico_l,jco_l) + pab(o1+ico,o2+jco)

    END IF

  END SUBROUTINE prepare_arb



END MODULE qs_modify_pab_block

