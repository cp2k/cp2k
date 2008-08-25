!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the integrals of the Vxc potential calculated
!>      for the atomic code
! *****************************************************************************
MODULE atom_output

  USE atom_types,                      ONLY: atom_state,&
                                             atom_type
  USE atom_utils,                      ONLY: get_maxl_occ,&
                                             get_maxn_occ
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_output'

  PUBLIC :: atom_print_state, atom_print_energies, atom_print_iteration

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_print_state(state,iw,error)
    TYPE(atom_state)                         :: state
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_state', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=1), DIMENSION(0:3)         :: label = (/ "S","P","D","F" /)
    INTEGER                                  :: j, l, mc, mlc, mlo, mm(0:3), &
                                                mo
    LOGICAL                                  :: failure = .FALSE.

    WRITE(iw,'(//,A)') " Electronic structure "
    WRITE(iw,'(T5,A,T71,F10.2)') "Total number of core electrons",SUM(state%core)
    WRITE(iw,'(T5,A,T71,F10.2)') "Total number of valence electrons",SUM(state%occ)
    WRITE(iw,'(T5,A,T71,F10.2)') "Total number of electrons",SUM(state%occ+state%core)

    mlo = get_maxl_occ(state%occ)
    mlc = get_maxl_occ(state%core)
    mm  = get_maxn_occ(state%core)

    DO l=0,MAX(mlo,mlc)
      mo = state%maxn_occ(l)
      IF ( SUM(state%core(l,:)) == 0 ) THEN
         WRITE(iw,'(A5,T10,10F6.2)') label(l),(state%occ(l,j),j=1,mo)
      ELSE
         mc = mm(l)
         CPPrecondition(SUM(state%occ(l,1:mc))==0,cp_failure_level,routineP,error,failure)
         WRITE(iw,ADVANCE="no",FMT='(A5,T9,A1,10F6.2)') label(l),"[",(state%core(l,j),j=1,mc)
         WRITE(iw,FMT='(A1,F5.2,10F6.2)') "]",(state%occ(l,j),j=mc+1,mo)
      END IF
    END DO
    WRITE(iw,*)

  END SUBROUTINE atom_print_state
! *****************************************************************************
  SUBROUTINE atom_print_energies(atom,iw,error)
    TYPE(atom_type)                          :: atom
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_energies', &
      routineP = moduleN//':'//routineN

    WRITE(iw,'(/,A)') " Energy components [Hartree]"
    WRITE(iw,'(T10,A,T61,F20.12)') "    Total Energy ::",atom%energy%etot
    WRITE(iw,'(T10,A,T61,F20.12)') "     Band Energy ::",atom%energy%eband
    WRITE(iw,'(T10,A,T61,F20.12)') "  Kinetic Energy ::",atom%energy%ekin
    WRITE(iw,'(T10,A,T61,F20.12)') "     Core Energy ::",atom%energy%ecore
    WRITE(iw,'(T10,A,T61,F20.12)') "       XC Energy ::",atom%energy%exc
    WRITE(iw,'(T10,A,T61,F20.12)') "  Coulomb Energy ::",atom%energy%ecoulomb
    WRITE(iw,'(T10,A,T61,F20.12)') " Exchange Energy ::",atom%energy%eexchange
    IF ( atom%ppot_type /= 0 ) THEN
       WRITE(iw,'(T10,A,T61,F20.12)') "    Total Pseudopotential Energy ::",atom%energy%epseudo
       WRITE(iw,'(T10,A,T61,F20.12)') "    Local Pseudopotential Energy ::",atom%energy%eploc
       WRITE(iw,'(T10,A,T61,F20.12)') " Nonlocal Pseudopotential Energy ::",atom%energy%epnl
    END IF

  END SUBROUTINE atom_print_energies
! *****************************************************************************
  SUBROUTINE atom_print_iteration(iter,deps,etot,iw,error)
    INTEGER, INTENT(IN)                      :: iter
    REAL(dp), INTENT(IN)                     :: deps, etot
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_print_iteration', &
      routineP = moduleN//':'//routineN

    IF (iter==1) THEN
      WRITE(iw,'(/," ",79("*"),/,T19,A,T38,A,T70,A,/," ",79("*"))') &
           "Iteration","Convergence","Energy [au]"
    END IF
    WRITE(iw,'(T20,i8,T34,G14.6,T61,F20.12)') iter,deps,etot

  END SUBROUTINE atom_print_iteration
! *****************************************************************************

END MODULE atom_output

! *****************************************************************************
