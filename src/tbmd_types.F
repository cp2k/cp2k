!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****s* cp2k/tbmd_types [1.0] *
!!
!!   NAME
!!     tbmd_types
!!
!!   FUNCTION
!!     Define basic types for Tight-Binding Calculations
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

    MODULE tbmd_types

  USE kinds,                           ONLY: dp
!


      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: tbatom_param_type, gsp_param_type, tb_hopping_type
!
!!    Maximum number of orbitals of each l-quantum number
      INTEGER, PARAMETER :: mx_orb = 3
!
      TYPE tbatom_param_type
        CHARACTER (len=10) :: atom_name
        INTEGER :: atomic_number
        REAL(KIND=dp) :: valence_charge
        INTEGER :: natoms
        INTEGER :: nbsf
        INTEGER, DIMENSION (0:3) :: bsf
        REAL(KIND=dp), DIMENSION (0:3,1:mx_orb) :: orbital_energy
        REAL(KIND=dp), DIMENSION (:), POINTER :: h_diag
        INTEGER, DIMENSION (:,:), POINTER :: spoint
      END TYPE tbatom_param_type
!
      TYPE gsp_param_type !Goodwin parameter set
        REAL(KIND=dp) :: er0
        REAL(KIND=dp) :: r0, rt
        REAL(KIND=dp) :: na, nb
      END TYPE gsp_param_type
!
      TYPE tb_hopping_type
        CHARACTER (len=10) :: htype
        REAL(KIND=dp) :: rcutsq
        INTEGER :: nint
        TYPE (gsp_param_type), POINTER :: gsp(:)
        INTEGER, POINTER :: lls(:,:)
        INTEGER, POINTER :: llpoint(:,:)
      END TYPE tb_hopping_type
!
!------------------------------------------------------------------------------!
    END MODULE tbmd_types
!------------------------------------------------------------------------------!
!!*****

