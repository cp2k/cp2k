
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/qs_modify_pab_block [1.0] *
!!
!!   NAME
!!     qs_modify_pab_block
!!
!!   FUNCTION
!!     Set of procedures to modify the density matrix, in order to collocate
!!     on the real space grid specific functions of the primitives
!!
!!   AUTHOR
!!     MI 06.2006
!!
!
!!   SOURCE
!******************************************************************************

MODULE  qs_modify_pab_block

! *****************************************************************************

  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_modify_pab_block'

! *** Public subroutines ***

  PUBLIC :: prepare_dadb,prepare_adb_m_dab,prepare_ardb_m_darb
  PUBLIC :: FUNC_AB, FUNC_DADB, FUNC_ADBmDAB, FUNC_ARDBmDARB

  ! Flags for the type of function to be collocated
  INTEGER, PARAMETER :: FUNC_AB = 401, FUNC_DADB = 402
  INTEGER, PARAMETER :: FUNC_ADBmDAB = 403, FUNC_ARDBmDARB = 404

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
!
!!****f* qs_modify_pab_block/prepare_dadb *
!!
!!   NAME
!!     prepare_dadb
!!
!!   FUNCTION
!!     create a new pab_local so that mapping pab_local with pgf_a pgf_b
!!     is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
!!     (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})
!!
!!
!!*** **********************************************************************

  SUBROUTINE prepare_dadb(pab_local,pab,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

   REAL(dp),DIMENSION(:,:), INTENT(INOUT) :: pab_local
   REAL(dp),DIMENSION(:,:), INTENT(IN)    :: pab
   INTEGER, INTENT(IN)                    :: lxa, lya, lza
   INTEGER, INTENT(IN)                    :: lxb, lyb, lzb
   INTEGER, INTENT(IN)                    :: o1, o2
   REAL(dp), INTENT(IN)                   :: zeta, zetb

   INTEGER :: ico,ico_l,jco,jco_l

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

! *****************************************************************************
!
!!****f* qs_modify_pab_block/prepare_adb_m_dab *
!!
!!   NAME
!!     prepare_adb_m_dab
!!
!!   FUNCTION
!!     create a new pab_local so that mapping pab_local with pgf_a pgf_b
!!     is equivalent to mapping pab with  pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
!!     ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) = 
!!             pgf_a *(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x}) - (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b
!!
!!
!!*** **********************************************************************

  SUBROUTINE prepare_adb_m_dab(pab_local,pab,idir,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

   REAL(dp),DIMENSION(:,:), INTENT(INOUT) :: pab_local
   REAL(dp),DIMENSION(:,:), INTENT(IN)    :: pab
   INTEGER, INTENT(IN)                    :: idir
   INTEGER, INTENT(IN)                    :: lxa, lya, lza
   INTEGER, INTENT(IN)                    :: lxb, lyb, lzb
   INTEGER, INTENT(IN)                    :: o1, o2
   REAL(dp), INTENT(IN)                   :: zeta, zetb

   INTEGER :: ico,ico_l,jco,jco_l

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

! *****************************************************************************
!
!!****f* qs_modify_pab_block/prepare_ardb_m_darb
!!
!!   NAME
!!     prepare_ardb_m_darb
!!
!!   FUNCTION
!!     create a new pab_local so that mapping pab_local with pgf_a pgf_b
!!     pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}  pgf_b
!      ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b ) = 
!!                    pgf_a *(lbx pgf_{b-1x+1ir} - 2*zetb*pgf_{b+1x+1ir}) - 
!!                             (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}

!!
!!
!!*** **********************************************************************

  SUBROUTINE prepare_ardb_m_darb(pab_local,pab,idir,ir,lxa,lya,lza,lxb,lyb,lzb,o1,o2,zeta,zetb)

   REAL(dp),DIMENSION(:,:), INTENT(INOUT) :: pab_local
   REAL(dp),DIMENSION(:,:), INTENT(IN)    :: pab
   INTEGER, INTENT(IN)                    :: idir
   INTEGER, INTENT(IN)                    :: ir
   INTEGER, INTENT(IN)                    :: lxa, lya, lza
   INTEGER, INTENT(IN)                    :: lxb, lyb, lzb
   INTEGER, INTENT(IN)                    :: o1, o2
   REAL(dp), INTENT(IN)                   :: zeta, zetb

   INTEGER :: ico,ico_l,jco,jco_l

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

!  *****************************************************************************************
END MODULE qs_modify_pab_block

