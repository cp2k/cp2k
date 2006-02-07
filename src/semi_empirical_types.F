!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 CP2K developers group                                  !
!-----------------------------------------------------------------------------!
!!****** cp2k/semi_empirical_types [1.0] *
!!
!!   NAME
!!     semi_empirical_types
!!
!!   FUNCTION
!!     Definition of the semi empirical parameter types.
!!
!!   AUTHOR
!!     JGH (14.08.2004)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE semi_empirical_types

  USE basis_set_types,                 ONLY: sto_basis_set_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_types'

  TYPE semi_empirical_type
    ! PRIVATE
    CHARACTER(LEN=default_string_length)   :: typ
    CHARACTER(LEN=default_string_length)   :: name
    LOGICAL                                :: defined
    INTEGER                                :: z
    REAL(KIND=dp)                          :: zeff
    INTEGER                                :: natorb
    REAL(KIND = dp), DIMENSION(:), POINTER :: beta
    REAL(KIND = dp), DIMENSION(:), &
                                   POINTER :: sto_exponents
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
    REAL(KIND = dp),DIMENSION(4)           :: fn1,fn2,fn3
    REAL(KIND = dp),DIMENSION(4,4)         :: bfn1,bfn2,bfn3
  END TYPE semi_empirical_type

! *** Public data types ***

  PUBLIC :: semi_empirical_type

!!***
! *****************************************************************************

END MODULE semi_empirical_types
