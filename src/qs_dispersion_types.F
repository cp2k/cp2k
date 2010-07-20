!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition of disperson types for DFT calculations
!> \author JGH (20.10.2008)
! *****************************************************************************
MODULE qs_dispersion_types

  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_types'

! *****************************************************************************
  INTEGER, PARAMETER                       :: dftd2_pp = 1
  INTEGER, PARAMETER                       :: dftd3_pp = 2

  TYPE qs_dispersion_type
    INTEGER                                :: TYPE
    CHARACTER (LEN=default_string_length)  :: pp_type
    CHARACTER (LEN=default_string_length)  :: ref_functional
    REAL(KIND=dp)                          :: scaling
    REAL(KIND=dp)                          :: rc_disp
    REAL(KIND=dp)                          :: exp_pre
    LOGICAL                                :: verbose  !extended output
    !DFT-D3 global parameters
    INTEGER                                :: max_elem !elements parametrized
    INTEGER                                :: maxc     !max coordination number references per element
    REAL(KIND=dp)                          :: k1,k2,k3 !ad hoc parameters
    REAL(KIND=dp)                          :: alp      !ad hoc parameters
    REAL(KIND=dp)                          :: s6,s8,sr6!scaling parameters
    REAL(KIND=dp)                          :: eps_cn
    LOGICAL                                :: doabc    !neglect C9 terms
    LOGICAL                                :: c9cnst   !use constant c9 terms
    LOGICAL                                :: lrc      !calculate a long range correction
    REAL(KIND=dp),DIMENSION(:,:,:,:,:), POINTER &
                                           :: c6ab
    INTEGER,DIMENSION(:),POINTER           :: maxci
    REAL(KIND=dp),DIMENSION(:,:), POINTER  :: r0ab
    REAL(KIND=dp),DIMENSION(:), POINTER    :: rcov     !covalent radii
    REAL(KIND=dp),DIMENSION(:), POINTER    :: r2r4     !atomic <r^2>/<r^4> values
    REAL(KIND=dp),DIMENSION(:), POINTER    :: cn       !coordination numbers (defaults)
    TYPE(cn_kind_list),DIMENSION(:), POINTER &
                                           :: cnkind
    TYPE(cn_atom_list),DIMENSION(:), POINTER &
                                           :: cnlist
  END TYPE qs_dispersion_type

  TYPE qs_atom_dispersion_type
    INTEGER                                :: TYPE
    LOGICAL                                :: defined
    REAL(KIND=dp)                          :: vdw_radii  !van der Waals radii
    REAL(KIND=dp)                          :: c6         !c6 coefficients
  END TYPE qs_atom_dispersion_type

  TYPE cn_kind_list
    REAL(KIND=dp)                          :: cnum
    INTEGER                                :: kind
  END TYPE cn_kind_list
  TYPE cn_atom_list
    REAL(KIND=dp)                          :: cnum
    INTEGER                                :: natom
    INTEGER,DIMENSION(:), POINTER          :: atom
  END TYPE cn_atom_list

! *****************************************************************************

  PUBLIC :: qs_atom_dispersion_type, qs_dispersion_type, cn_atom_list, cn_kind_list
  PUBLIC :: dftd2_pp, dftd3_pp

! *****************************************************************************

END MODULE qs_dispersion_types

