!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/external_potential_types [1.0] *
!!
!!   NAME
!!     external_potential_types
!!
!!   FUNCTION
!!     Definition of the atomic potential types.
!!
!!   AUTHOR
!!     Matthias Krack (04.07.2000)
!!
!!   MODIFICATION HISTORY
!!     GT, 22.09.2002: added elp_potential_types
!!
!!   SOURCE
!******************************************************************************

MODULE external_potential_types

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: dfac,&
                                             fac,&
                                             pi
  USE mathlib,                         ONLY: symmetrize_matrix
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: co,&
                                             coset,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset,&
                                             nso
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE periodic_table,                  ONLY: ptable
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             test_object
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "external_potential_types"

! *** Define the all-electron potential type ***

! *** Literature: M. Krack and M. Parrinello,
!                 Phys. Chem. Chem. Phys. 2, 2105 (2000)

  TYPE all_potential_type
    PRIVATE
    CHARACTER(LEN=default_string_length)   :: name
    CHARACTER(LEN=default_string_length),&
      DIMENSION(2)                         :: description
    REAL(KIND = dp)                        :: alpha_core_charge,&
                                              ccore_charge,&
                                              core_charge_radius
    INTEGER                                :: z,zeff
    INTEGER, DIMENSION(:), POINTER         :: elec_conf
  END TYPE all_potential_type

! *** Define the GTH potential type ***

! *** Literature: - S. Goedecker, M. Teter and J. Hutter,
!                   Phys. Rev. B 54, 1703 (1996)
!                 - C. Hartwigsen, S. Goedecker and J. Hutter,
!                   Phys. Rev. B 58, 3641 (1998)

  TYPE gth_potential_type
    PRIVATE
    CHARACTER(LEN=default_string_length)       :: name
    CHARACTER(LEN=default_string_length),&
      DIMENSION(3)                             :: description
    REAL(KIND = dp)                            :: alpha_core_charge,&
                                                  alpha_ppl,ccore_charge,&
                                                  cerf_ppl,&
                                                  core_charge_radius,&
                                                  ppl_radius,ppnl_radius
    INTEGER                                    :: lppnl,lprj_ppnl_max,&
                                                  nexp_ppl,nppnl,&
                                                  nprj_ppnl_max,z,zeff
    REAL(KIND = dp), DIMENSION(:), POINTER     :: alpha_ppnl,cexp_ppl
    INTEGER, DIMENSION(:), POINTER             :: elec_conf,nprj_ppnl
    REAL(KIND = dp), DIMENSION(:,:), POINTER   :: cprj,cprj_ppnl,vprj_ppnl
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER :: hprj_ppnl
  END TYPE gth_potential_type

! *** Define the effective point charge potential type ***

  TYPE elp_potential_type
    PRIVATE
    CHARACTER(LEN=default_string_length)   :: name
    CHARACTER(LEN=default_string_length),&
      DIMENSION(2) :: description
    REAL(KIND = dp)                        :: alpha_core_charge,&
                                              ccore_charge,&
                                              cerf_ppl,&
                                              core_charge_radius,&
                                              ppl_radius,qeff 
    INTEGER                                :: nalpha_ppl
    INTEGER,DIMENSION (:), POINTER         :: nexp_ppl               
    REAL(KIND = dp),DIMENSION(:), POINTER  :: alpha_ppl
    REAL(KIND = dp),DIMENSION(:,:),POINTER :: cexp_ppl   
  END TYPE elp_potential_type

! *** Public subroutines ***

  PUBLIC :: allocate_potential,&
            deallocate_potential,&
            get_potential,&
            init_potential,&
            read_potential,&
            set_potential,&
            set_default_all_potential,&
            write_potential

! *** Public data types ***

  PUBLIC :: all_potential_type,&
            elp_potential_type,&
            gth_potential_type

!!***
! *****************************************************************************

  INTERFACE allocate_potential
    MODULE PROCEDURE allocate_all_potential,&
                     allocate_elp_potential,&
                     allocate_gth_potential
  END INTERFACE

  INTERFACE deallocate_potential
    MODULE PROCEDURE deallocate_all_potential,&
                     deallocate_elp_potential,&
                     deallocate_gth_potential
  END INTERFACE

  INTERFACE get_potential
    MODULE PROCEDURE get_all_potential,&
                     get_elp_potential,&
                     get_gth_potential
  END INTERFACE

  INTERFACE init_potential
    MODULE PROCEDURE init_all_potential,&
                     init_gth_potential
  END INTERFACE

  INTERFACE read_potential
    MODULE PROCEDURE read_all_potential,&
                     read_elp_potential,&
                     read_gth_potential
  END INTERFACE

  INTERFACE set_potential
    MODULE PROCEDURE set_all_potential,&
                     set_elp_potential,&
                     set_gth_potential
  END INTERFACE

  INTERFACE write_potential
    MODULE PROCEDURE write_all_potential,&
                     write_elp_potential,&
                     write_gth_potential
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_all_potential(potential)

!   Purpose: Allocate an atomic all-electron potential data set.

!   History: Creation (25.07.2000,MK)

!   ***************************************************************************

    TYPE(all_potential_type), POINTER        :: potential

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_all_potential"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

    ALLOCATE (potential,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "potential",0)

    NULLIFY (potential%elec_conf)

    potential%description(1) = "All-Electron potential"
    potential%description(2) =&
      "- Krack, Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)"

  END SUBROUTINE allocate_all_potential

! *****************************************************************************

  SUBROUTINE allocate_elp_potential(potential)

!   Purpose: Allocate an effective charge potential data set.
!   History: Creation (22.09.2002, gt)

!   ***************************************************************************

    TYPE(elp_potential_type), POINTER        :: potential

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_elp_potential"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

    ALLOCATE (potential,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "potential",0)
    NULLIFY (potential%alpha_ppl)
    NULLIFY (potential%cexp_ppl)
    NULLIFY (potential%nexp_ppl)  

    potential%description(1) = "Local pseudopotential for KG-MD"

  END SUBROUTINE allocate_elp_potential

! *****************************************************************************

  SUBROUTINE allocate_gth_potential(potential)

!   Purpose: Allocate an atomic GTH potential data set.

!   History: Creation (25.07.2000,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_gth_potential"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) CALL deallocate_potential(potential)

    ALLOCATE (potential,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "potential",0)

    NULLIFY (potential%alpha_ppnl)
    NULLIFY (potential%cexp_ppl)
    NULLIFY (potential%elec_conf)
    NULLIFY (potential%nprj_ppnl)
    NULLIFY (potential%cprj)
    NULLIFY (potential%cprj_ppnl)
    NULLIFY (potential%vprj_ppnl)
    NULLIFY (potential%hprj_ppnl)

    potential%description(1) = "Goedecker-Teter-Hutter pseudopotential"
    potential%description(2) =&
      "- Goedecker, Teter, Hutter, PRB 54, 1703 (1996)"
    potential%description(3) =&
      "- Hartwigsen, Goedecker, Hutter, PRB 58, 3641 (1998)"

  END SUBROUTINE allocate_gth_potential

! *****************************************************************************

  SUBROUTINE deallocate_all_potential(potential)

!   Purpose: Deallocate an atomic all-electron potential data set.

!   History: Creation (03.11.2000,MK)

!   ***************************************************************************

    TYPE(all_potential_type), POINTER        :: potential

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_all_potential"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      DEALLOCATE (potential%elec_conf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "potential%elec_conf")

      DEALLOCATE (potential,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "potential")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE deallocate_all_potential

! *****************************************************************************

  SUBROUTINE deallocate_elp_potential(potential)

!   Purpose: deallocate an effective charge potential data set.
!   History: Creation (22.09.2002, gt)

!   ***************************************************************************

    TYPE(elp_potential_type), POINTER        :: potential

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_elp_potential"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      IF (ASSOCIATED(potential%nexp_ppl)) THEN
        DEALLOCATE (potential%nexp_ppl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%nexp_ppl")
      END IF

      IF (ASSOCIATED(potential%cexp_ppl)) THEN
        DEALLOCATE (potential%cexp_ppl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%cexp_ppl")
      END IF

      IF (ASSOCIATED(potential%alpha_ppl)) THEN
        DEALLOCATE (potential%alpha_ppl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%alpha_ppl")
      END IF

      DEALLOCATE (potential,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "potential")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE deallocate_elp_potential

! *****************************************************************************

  SUBROUTINE deallocate_gth_potential(potential)

!   Purpose: Deallocate an atomic GTH potential data set.

!   History: Creation (03.11.2000,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_gth_potential"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      DEALLOCATE (potential%elec_conf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "potential%elec_conf")

!     *** Deallocate the parameters of the local part ***

      IF (potential%nexp_ppl > 0) THEN
        DEALLOCATE (potential%cexp_ppl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%cexp_ppl")
      END IF

!     *** Deallocate the parameters of the non-local part ***

      IF (potential%nppnl > 0) THEN

        DEALLOCATE (potential%alpha_ppnl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%alpha_ppnl")

        DEALLOCATE (potential%cprj,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%cprj")

        DEALLOCATE (potential%cprj_ppnl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%cprj_ppnl")

        DEALLOCATE (potential%hprj_ppnl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%hprj_ppnl")

        DEALLOCATE (potential%nprj_ppnl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%nprj_ppnl")

        DEALLOCATE (potential%vprj_ppnl,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "potential%vprj_ppnl")

      END IF

      DEALLOCATE (potential,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "potential")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE deallocate_gth_potential

! *****************************************************************************

  SUBROUTINE get_all_potential(potential,name,alpha_core_charge,&
                               ccore_charge,core_charge_radius,z,zeff,&
                               elec_conf)

!   Purpose: Get informations about an all-electron potential data set.

!   History: Creation (11.01.2002,MK)

!   ***************************************************************************

    TYPE(all_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius
    INTEGER, INTENT(OUT), OPTIONAL           :: z, zeff
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_all_potential"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      IF (PRESENT(name)) name = potential%name
      IF (PRESENT(alpha_core_charge))&
        alpha_core_charge = potential%alpha_core_charge
      IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
      IF (PRESENT(core_charge_radius))&
        core_charge_radius = potential%core_charge_radius
      IF (PRESENT(z)) z = potential%z
      IF (PRESENT(zeff)) zeff = potential%zeff
      IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE get_all_potential

! *****************************************************************************

  SUBROUTINE get_elp_potential(potential,name,alpha_core_charge,&
                               alpha_ppl,cerf_ppl,cexp_ppl,nalpha_ppl,&
                               ccore_charge,core_charge_radius,&
                               nexp_ppl,ppl_radius,ppl_present,qeff)

!   Purpose: Get informations about effective pint charge potential.
!   History: 22-09-2002 GT 
!            28-11-2003 gt: supports a local gth-like analytical form

!   ***************************************************************************

    TYPE(elp_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha_core_charge
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppl
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: cerf_ppl
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cexp_ppl
    INTEGER, OPTIONAL                        :: nalpha_ppl
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: ccore_charge, &
                                                core_charge_radius
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nexp_ppl
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: ppl_radius
    LOGICAL, INTENT(OUT), OPTIONAL           :: ppl_present
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: qeff

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_elp_potential"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      IF (PRESENT(name)) name = potential%name
      IF (PRESENT(qeff)) qeff = potential%qeff
      IF (PRESENT(nalpha_ppl)) nalpha_ppl=potential%nalpha_ppl
      IF (PRESENT(alpha_ppl)) alpha_ppl => potential%alpha_ppl
      IF (PRESENT(alpha_core_charge))&
        alpha_core_charge = potential%alpha_core_charge
      IF (PRESENT(cerf_ppl)) cerf_ppl=potential%cerf_ppl
      IF (PRESENT(cexp_ppl)) cexp_ppl => potential%cexp_ppl
      IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
      IF (PRESENT(core_charge_radius))&
        core_charge_radius = potential%core_charge_radius
      IF (PRESENT(ppl_radius)) ppl_radius=potential%ppl_radius
      IF (PRESENT(nexp_ppl)) nexp_ppl => potential%nexp_ppl
      IF (PRESENT(ppl_present)) ppl_present = (ASSOCIATED(potential%nexp_ppl) )

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE get_elp_potential

! *****************************************************************************

  SUBROUTINE get_gth_potential(potential,name,alpha_core_charge,&
                               alpha_ppl,ccore_charge,cerf_ppl,&
                               core_charge_radius,ppl_radius,ppnl_radius,&
                               lppnl,lprj_ppnl_max,nexp_ppl,nppnl,&
                               nprj_ppnl_max,z,zeff,ppl_present,ppnl_present,&
                               alpha_ppnl,cexp_ppl,elec_conf,nprj_ppnl,cprj,&
                               cprj_ppnl,vprj_ppnl,hprj_ppnl)

!   Purpose: Get informations about a GTH potential data set.

!   History: Creation (11.01.2002,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: alpha_core_charge, alpha_ppl, &
      ccore_charge, cerf_ppl, core_charge_radius, ppl_radius, ppnl_radius
    INTEGER, INTENT(OUT), OPTIONAL           :: lppnl, lprj_ppnl_max, &
                                                nexp_ppl, nppnl, &
                                                nprj_ppnl_max, z, zeff
    LOGICAL, INTENT(OUT), OPTIONAL           :: ppl_present, ppnl_present
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppnl, cexp_ppl
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf, nprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cprj, cprj_ppnl, vprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: hprj_ppnl

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_gth_potential"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      IF (PRESENT(name)) name = potential%name
      IF (PRESENT(alpha_core_charge))&
        alpha_core_charge = potential%alpha_core_charge
      IF (PRESENT(alpha_ppl)) alpha_ppl = potential%alpha_ppl
      IF (PRESENT(ccore_charge)) ccore_charge = potential%ccore_charge
      IF (PRESENT(cerf_ppl)) cerf_ppl = potential%cerf_ppl
      IF (PRESENT(core_charge_radius))&
        core_charge_radius = potential%core_charge_radius
      IF (PRESENT(ppl_radius)) ppl_radius = potential%ppl_radius
      IF (PRESENT(ppnl_radius)) ppnl_radius = potential%ppnl_radius
      IF (PRESENT(lppnl)) lppnl = potential%lppnl
      IF (PRESENT(lprj_ppnl_max)) lprj_ppnl_max = potential%lprj_ppnl_max
      IF (PRESENT(nexp_ppl)) nexp_ppl = potential%nexp_ppl
      IF (PRESENT(nppnl)) nppnl = potential%nppnl
      IF (PRESENT(nprj_ppnl_max)) nprj_ppnl_max = potential%nprj_ppnl_max
      IF (PRESENT(z)) z = potential%z
      IF (PRESENT(zeff)) zeff = potential%zeff
      IF (PRESENT(ppl_present)) ppl_present = (potential%nexp_ppl > 0)
      IF (PRESENT(ppnl_present)) ppnl_present = (potential%nppnl > 0)
      IF (PRESENT(alpha_ppnl)) alpha_ppnl => potential%alpha_ppnl
      IF (PRESENT(cexp_ppl)) cexp_ppl => potential%cexp_ppl
      IF (PRESENT(elec_conf)) elec_conf => potential%elec_conf
      IF (PRESENT(nprj_ppnl)) nprj_ppnl => potential%nprj_ppnl
      IF (PRESENT(cprj)) cprj => potential%cprj
      IF (PRESENT(cprj_ppnl)) cprj_ppnl => potential%cprj_ppnl
      IF (PRESENT(vprj_ppnl)) vprj_ppnl => potential%vprj_ppnl
      IF (PRESENT(hprj_ppnl)) hprj_ppnl => potential%hprj_ppnl

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE get_gth_potential

! *****************************************************************************

  SUBROUTINE init_cprj_ppnl(potential)

!   Purpose: Initialise the coefficients of the projectors of the non-local
!            part of the GTH pseudopotential and the transformation matrices
!            for Cartesian overlap integrals between the orbital basis
!            functions and the projector functions.

!   History: - Creation (16.10.2000,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential

    INTEGER                                  :: cpx, cpy, cpz, cx, cy, cz, &
                                                ico, iprj, iprj_ppnl, l, lp, &
                                                lprj_ppnl, nprj, px, py, pz
    REAL(KIND=dp)                            :: alpha_ppnl, cp

!   ---------------------------------------------------------------------------

    nprj = 0

    DO l=0,potential%lppnl
      alpha_ppnl = potential%alpha_ppnl(l)
      DO iprj_ppnl=1,potential%nprj_ppnl(l)
        lp = iprj_ppnl - 1
        lprj_ppnl = l + 2*lp
        cp = SQRT(2.0_dp**(2.0_dp*REAL(lprj_ppnl,dp) + 3.5_dp)*&
                  alpha_ppnl**(REAL(lprj_ppnl,dp) + 1.5_dp)/&
                  (SQRT(pi)*dfac(2*lprj_ppnl + 1)))
        potential%cprj_ppnl(iprj_ppnl,l) = cp
        DO cx=0,l
          DO cy=0,l-cx
            cz = l - cx - cy
            iprj = nprj + co(cx,cy,cz)
            DO px=0,lp
              DO py=0,lp-px
                pz = lp - px - py
                cpx = cx + 2*px
                cpy = cy + 2*py
                cpz = cz + 2*pz
                ico = coset(cpx,cpy,cpz)
                potential%cprj(ico,iprj) = cp*fac(lp)/(fac(px)*fac(py)*fac(pz))
              END DO
            END DO
          END DO
        END DO
        nprj = nprj + nco(l)
      END DO
    END DO

  END SUBROUTINE init_cprj_ppnl

! *****************************************************************************

  SUBROUTINE init_gth_potential(potential)

!   Purpose: Initialise a GTH potential data set structure.

!   History: Creation (27.10.2000,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(potential)) RETURN

    IF (potential%nppnl > 0) THEN

!     *** Initialise the projector coefficients of the    ***
!     *** non-local part of the GTH pseudopotential and   ***
!     *** the transformation matrices "pgf" -> "prj_ppnl" ***

      CALL init_cprj_ppnl(potential)

!     *** Initialise the h(i,j) projector coefficients of ***
!     *** the non-local part of the GTH pseudopotential   ***

      CALL init_vprj_ppnl(potential)

    END IF

  END SUBROUTINE init_gth_potential

! *****************************************************************************

  SUBROUTINE init_vprj_ppnl(potential)

!   Purpose: Initialise the h(i,j) projector coefficients of the non-local part
!            of the GTH pseudopotential.

!   History: - Creation (24.10.2000,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential

    INTEGER                                  :: i, ico, iprj, iprj_ppnl, iso, &
                                                j, jco, jprj, jprj_ppnl, l, &
                                                nprj

!   ---------------------------------------------------------------------------

    nprj = 0

    DO l=0,potential%lppnl
      DO iprj_ppnl=1,potential%nprj_ppnl(l)
        iprj = nprj + (iprj_ppnl - 1)*nco(l)
        DO jprj_ppnl=1,potential%nprj_ppnl(l)
          jprj = nprj + (jprj_ppnl - 1)*nco(l)
          DO ico=1,nco(l)
            i = iprj + ico
            DO jco=1,nco(l)
              j = jprj + jco
              DO iso=1,nso(l)
                potential%vprj_ppnl(i,j) = potential%vprj_ppnl(i,j) +&
                                           orbtramat(l)%slm(iso,ico)*&
                                           potential%hprj_ppnl(iprj_ppnl,&
                                                               jprj_ppnl,l)*&
                                           orbtramat(l)%slm(iso,jco)
              END DO
            END DO
          END DO
        END DO
      END DO
      nprj = nprj + potential%nprj_ppnl(l)*nco(l)
    END DO

  END SUBROUTINE init_vprj_ppnl

! *****************************************************************************

  SUBROUTINE init_all_potential(potential,itype,zeff)

    TYPE(all_potential_type), POINTER        :: potential
    CHARACTER (LEN=*), OPTIONAL              :: itype
    INTEGER, OPTIONAL                        :: zeff

    INTEGER   :: dz

!   ---------------------------------------------------------------------------

    IF (.NOT.ASSOCIATED(potential)) RETURN

    IF ( PRESENT (zeff) ) THEN
      potential%zeff = zeff
    END IF
    dz = potential%z - potential%zeff
    SELECT CASE (dz)
      CASE DEFAULT
      CASE (2)
        potential%elec_conf(0) = potential%elec_conf(0) - 2
      CASE (10)
        potential%elec_conf(0) = potential%elec_conf(0) - 4
        potential%elec_conf(1) = potential%elec_conf(1) - 6
      CASE (18)
        potential%elec_conf(0) = potential%elec_conf(0) - 6
        potential%elec_conf(1) = potential%elec_conf(1) - 12
      CASE (28)
        potential%elec_conf(0) = potential%elec_conf(0) - 6
        potential%elec_conf(1) = potential%elec_conf(1) - 12
        potential%elec_conf(2) = potential%elec_conf(2) - 10
      CASE (30)
        potential%elec_conf(0) = potential%elec_conf(0) - 8
        potential%elec_conf(1) = potential%elec_conf(1) - 12
        potential%elec_conf(2) = potential%elec_conf(2) - 10
      CASE (36)
        potential%elec_conf(0) = potential%elec_conf(0) - 8
        potential%elec_conf(1) = potential%elec_conf(1) - 18
        potential%elec_conf(2) = potential%elec_conf(2) - 10
      CASE (46)
        potential%elec_conf(0) = potential%elec_conf(0) - 8
        potential%elec_conf(1) = potential%elec_conf(1) - 18
        potential%elec_conf(2) = potential%elec_conf(2) - 20
      CASE (48)
        potential%elec_conf(0) = potential%elec_conf(0) - 10
        potential%elec_conf(1) = potential%elec_conf(1) - 18
        potential%elec_conf(2) = potential%elec_conf(2) - 20
      CASE (54)
        potential%elec_conf(0) = potential%elec_conf(0) - 10
        potential%elec_conf(1) = potential%elec_conf(1) - 24
        potential%elec_conf(2) = potential%elec_conf(2) - 20
      CASE (68)
        potential%elec_conf(0) = potential%elec_conf(0) - 10
        potential%elec_conf(1) = potential%elec_conf(1) - 24
        potential%elec_conf(2) = potential%elec_conf(2) - 20
        potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (78)
        potential%elec_conf(0) = potential%elec_conf(0) - 10
        potential%elec_conf(1) = potential%elec_conf(1) - 24
        potential%elec_conf(2) = potential%elec_conf(2) - 30
        potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (80)
        potential%elec_conf(0) = potential%elec_conf(0) - 12
        potential%elec_conf(1) = potential%elec_conf(1) - 24
        potential%elec_conf(2) = potential%elec_conf(2) - 30
        potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (86)
        potential%elec_conf(0) = potential%elec_conf(0) - 12
        potential%elec_conf(1) = potential%elec_conf(1) - 30
        potential%elec_conf(2) = potential%elec_conf(2) - 30
        potential%elec_conf(3) = potential%elec_conf(3) - 14
      CASE (100)
        potential%elec_conf(0) = potential%elec_conf(0) - 12
        potential%elec_conf(1) = potential%elec_conf(1) - 30
        potential%elec_conf(2) = potential%elec_conf(2) - 30
        potential%elec_conf(3) = potential%elec_conf(3) - 28
    END SELECT

    IF ( PRESENT (itype) ) THEN
      IF ( itype == "BARE" ) THEN
        potential%description(1) = "Bare Coulomb Potential"
        IF ( dz > 0 ) THEN
          potential%description(2) = "Valence charge only"
        ELSE
          potential%description(2) = "Full atomic charge"
        END IF
      ENDIF
    END IF

  END SUBROUTINE init_all_potential

! *****************************************************************************

  SUBROUTINE read_all_potential(element_symbol,potential_name,potential,&
                                globenv)

!   Purpose: Read an atomic all-electron potential data set.

!   History: - Creation (14.05.2000,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, potential_name
    TYPE(all_potential_type), POINTER        :: potential
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_all_potential"

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    CHARACTER(LEN=LEN(potential_name))       :: apname
    CHARACTER(LEN=LEN(potential_name)+2)     :: apname2
    INTEGER                                  :: istat, l, strlen1, strlen2
    INTEGER, DIMENSION(:), POINTER           :: elec_conf
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: alpha, r

!   ---------------------------------------------------------------------------

    line2 = ""
    symbol2 = ""
    apname2 = ""

    potential%name = potential_name

    CALL start_parser(globenv%potential_file_name,globenv)

!   *** Search for the requested potential in the potential file   ***
!   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol

    search_loop: DO

      CALL search_object(apname,.TRUE.,found,line)

      IF (found) THEN

        CALL uppercase(line)
        CALL uppercase(symbol)
        CALL uppercase(apname)

!       *** Check both the element symbol and the atomic potential name ***

        line2 = " "//line//" "
        symbol2 = " "//TRIM(symbol)//" "
        apname2 = " "//TRIM(apname)//" "
        strlen1 = LEN_TRIM(symbol2) + 1
        strlen2 = LEN_TRIM(apname2) + 1

        IF ((INDEX(line2,symbol2(:strlen1)) > 0).AND.&
            (INDEX(line2,apname2(:strlen2)) > 0)) THEN

!         *** Read the electronic configuration ***

          NULLIFY (elec_conf)

          l = 0
          CALL reallocate(elec_conf,0,l)
          CALL read_object(elec_conf(l),newline=.TRUE.)
          DO WHILE (test_object() == "INT")
            l = l + 1
            CALL reallocate(elec_conf,0,l)
            CALL read_object(elec_conf(l))
          END DO

          CALL reallocate(potential%elec_conf,0,l)
          potential%elec_conf(:) = elec_conf(:)

          potential%zeff = SUM(elec_conf)

          DEALLOCATE (elec_conf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "elec_conf")

!         *** Read r(loc) to define the exponent of the core charge    ***
!         *** distribution and calculate the corresponding coefficient ***

          CALL read_object(r,newline=.TRUE.)

          alpha = 1.0_dp/(2.0_dp*r**2)

          potential%alpha_core_charge = alpha
          potential%ccore_charge = REAL(potential%zeff,dp)*SQRT((alpha/pi)**3)

          EXIT search_loop

        END IF

      ELSE

!       *** Stop program, if the end of file is reached ***

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The requested atomic potential <"//&
                          TRIM(potential_name)//&
                          "> was not found in the potential file <"//&
                          TRIM(globenv%potential_file_name)//">")

      END IF

    END DO search_loop

    CALL finish_parser()

  END SUBROUTINE read_all_potential

! *****************************************************************************

  SUBROUTINE read_elp_potential(element_symbol,potential_name,potential,&
                                globenv)

!   Purpose: Read an effective charge potential data set.
!   History: creation 22-09-2002, GT

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, potential_name
    TYPE(elp_potential_type), POINTER        :: potential
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_elp_potential"

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    CHARACTER(LEN=LEN(potential_name))       :: apname
    CHARACTER(LEN=LEN(potential_name)+2)     :: apname2
    INTEGER                                  :: i, j, maxlppl, n, nalpha_ppl, &
                                                nmax, strlen1, strlen2
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: alpha, qeff

!   ---------------------------------------------------------------------------

    line2 = ""
    symbol2 = ""
    apname2 = ""
    nalpha_ppl = 0
    nmax = 0
    

    potential%name = potential_name

    CALL start_parser(globenv%potential_file_name,globenv)

!   *** Search for the requested potential in the potential file   ***
!   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol

    search_loop: DO

      CALL search_object(TRIM(apname),.TRUE.,found,line)

      IF (found) THEN

        CALL uppercase(line)
        CALL uppercase(symbol)
        CALL uppercase(apname)

!       *** Check both the element symbol and the atomic potential name ***

        line2 = " "//line//" "
        symbol2 = " "//TRIM(symbol)//" "
        apname2 = " "//TRIM(apname)//" "
        strlen1 = LEN_TRIM(symbol2) + 1
        strlen2 = LEN_TRIM(apname2) + 1

        IF ((INDEX(line2,symbol2(:strlen1)) > 0).AND.&
            (INDEX(line2,apname2(:strlen2)) > 0)) THEN

!         *** Read the effective charge qeff ***

          CALL read_object(qeff,newline=.TRUE.)

          potential%qeff = qeff

!         *** read alpha and calculate the corresponding coefficient ***

          CALL read_object(alpha,newline=.TRUE.)

          potential%alpha_core_charge = alpha
          potential%ccore_charge = potential%qeff*SQRT((alpha/pi)**3)
          potential%cerf_ppl = potential%qeff*SQRT((alpha/pi)**3)

          CALL read_object(nalpha_ppl,newline=.TRUE.)

          potential%nalpha_ppl = nalpha_ppl

          IF (nalpha_ppl > 0) THEN
            CALL reallocate(potential%alpha_ppl,1,nalpha_ppl)
            CALL reallocate(potential%nexp_ppl,1,nalpha_ppl)
            DO i=1,nalpha_ppl
              CALL read_object(potential%alpha_ppl(i),newline=.TRUE.)
              CALL read_object(n)
              potential%nexp_ppl(i)=n
              IF (n > nmax) THEN
                nmax = n
                CALL reallocate(potential%cexp_ppl,1,nalpha_ppl,1,nmax)
              END IF
              DO j=1,n
                 CALL read_object(potential%cexp_ppl(i,j))
                 potential%cexp_ppl(i,j) = (2.0_dp*potential%alpha_ppl(i))**(j - 1)*&
                                           potential%cexp_ppl(i,j)
              END DO
            END DO

!           *** Check if alpha_core_charge = alpha_ppl(1) ***
 
            IF (potential%alpha_core_charge /= potential%alpha_ppl(1)) THEN
              CALL stop_program(routine_name,module_name,__LINE__,&
                                "exponent of the core charge and "//&
                                "alpha_ppl(1) must be equal")
            END IF

            maxlppl = 2*(nmax - 1)

            IF (maxlppl > -1)  CALL init_orbital_pointers(maxlppl)

          END IF 

          EXIT search_loop

        END IF

      ELSE

!       *** Stop program, if the end of file is reached ***

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The requested atomic potential <"//&
                          TRIM(potential_name)//&
                          "> was not found in the potential file <"//&
                          TRIM(globenv%potential_file_name)//">",globenv)

      END IF

    END DO search_loop

    CALL finish_parser()

  END SUBROUTINE read_elp_potential

! *****************************************************************************

  SUBROUTINE read_gth_potential(element_symbol,potential_name,potential,&
                                globenv)

!   Purpose: Read an atomic GTH potential data set.

!   Literature: - S. Goedecker, M. Teter and J. Hutter,
!                 Phys. Rev. B 54, 1703 (1996)
!               - C. Hartwigsen, S. Goedecker and J. Hutter,
!                 Phys. Rev. B 58, 3641 (1998)

!   History: - Creation (14.05.2000,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: element_symbol, potential_name
    TYPE(gth_potential_type), POINTER        :: potential
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_gth_potential"

    CHARACTER(LEN=240)                       :: line
    CHARACTER(LEN=242)                       :: line2
    CHARACTER(LEN=LEN(element_symbol))       :: symbol
    CHARACTER(LEN=LEN(element_symbol)+2)     :: symbol2
    CHARACTER(LEN=LEN(potential_name))       :: apname
    CHARACTER(LEN=LEN(potential_name)+2)     :: apname2
    INTEGER                                  :: i, istat, j, l, lppnl, &
                                                lprj_ppnl_max, maxlppl, n, &
                                                nppnl, nprj_ppnl, &
                                                nprj_ppnl_max, strlen1, &
                                                strlen2
    INTEGER, DIMENSION(:), POINTER           :: elec_conf
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: alpha, r
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: hprj_ppnl

!   ---------------------------------------------------------------------------

    line2 = ""
    symbol2 = ""
    apname2 = ""

    potential%name = potential_name

    CALL start_parser(globenv%potential_file_name,globenv)

!   *** Search for the requested potential in the potential file   ***
!   *** until the potential is found or the end of file is reached ***

    apname = potential_name
    symbol = element_symbol

    search_loop: DO

      CALL search_object(TRIM(apname),.TRUE.,found,line)

      IF (found) THEN

        CALL uppercase(line)
        CALL uppercase(symbol)
        CALL uppercase(apname)

!       *** Check both the element symbol and the atomic potential name ***

        line2 = " "//line//" "
        symbol2 = " "//TRIM(symbol)//" "
        apname2 = " "//TRIM(apname)//" "
        strlen1 = LEN_TRIM(symbol2) + 1
        strlen2 = LEN_TRIM(apname2) + 1

        IF ((INDEX(line2,symbol2(:strlen1)) > 0).AND.&
            (INDEX(line2,apname2(:strlen2)) > 0)) THEN

!         *** Read the electronic configuration ***

          NULLIFY (elec_conf)
          l = 0
          CALL reallocate(elec_conf,0,l)
          CALL read_object(elec_conf(l),newline=.TRUE.)
          DO WHILE (test_object() == "INT")
            l = l + 1
            CALL reallocate(elec_conf,0,l)
            CALL read_object(elec_conf(l))
          END DO

          CALL reallocate(potential%elec_conf,0,l)
          potential%elec_conf(:) = elec_conf(:)

          potential%zeff = SUM(elec_conf)

          DEALLOCATE (elec_conf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "elec_conf")

!         *** Read r(loc) to define the exponent of the core charge    ***
!         *** distribution and calculate the corresponding coefficient ***

          CALL read_object(r,newline=.TRUE.)

          alpha = 1.0_dp/(2.0_dp*r**2)

          potential%alpha_core_charge = alpha
          potential%ccore_charge = REAL(potential%zeff,dp)*SQRT((alpha/pi)**3)

          potential%alpha_ppl = alpha
          potential%cerf_ppl = REAL(potential%zeff,dp)*SQRT((alpha/pi)**3)

!         *** Read the parameters for the local part ***
!         *** of the GTH pseudopotential (ppl)       ***

          CALL read_object(n)
          potential%nexp_ppl = n
          CALL reallocate(potential%cexp_ppl,1,n)

          DO i=1,n
            CALL read_object(potential%cexp_ppl(i))
            potential%cexp_ppl(i) = (2.0_dp*potential%alpha_ppl)**(i - 1)*&
                                    potential%cexp_ppl(i)
          END DO

          maxlppl = 2*(n - 1)

          IF (maxlppl > -1) CALL init_orbital_pointers(maxlppl)

!         *** Read the parameters for the non-local  ***
!         *** part of the GTH pseudopotential (ppnl) ***

          CALL read_object(n,newline=.TRUE.)

          potential%lppnl = n - 1
          potential%nppnl = 0

          potential%lprj_ppnl_max = n - 1
          potential%nprj_ppnl_max = 0

          IF (n > 0) THEN

            lppnl = potential%lppnl
            nppnl = potential%nppnl

            CALL init_orbital_pointers(lppnl)

            NULLIFY (hprj_ppnl)

!           *** Load the parameter for n non-local projectors ***

            CALL reallocate(potential%alpha_ppnl,0,lppnl)
            CALL reallocate(potential%nprj_ppnl,0,lppnl)

            lprj_ppnl_max = -1
            nprj_ppnl_max = 0

            DO l=0,lppnl
              CALL read_object(r,newline=.TRUE.)
              potential%alpha_ppnl(l) = 1.0_dp/(2.0_dp*r**2)
              CALL read_object(nprj_ppnl)
              potential%nprj_ppnl(l) = nprj_ppnl
              nppnl = nppnl + nprj_ppnl*nco(l)
              IF (nprj_ppnl > nprj_ppnl_max) THEN
                nprj_ppnl_max = nprj_ppnl
                CALL reallocate(hprj_ppnl,1,nprj_ppnl_max,&
                                          1,nprj_ppnl_max,&
                                          0,lppnl)
              END IF
              DO i=1,nprj_ppnl
                IF (i == 1) THEN
                  CALL read_object(hprj_ppnl(i,i,l))
                ELSE
                  CALL read_object(hprj_ppnl(i,i,l),newline=.TRUE.)
                END IF
                DO j=i+1,nprj_ppnl
                  CALL read_object(hprj_ppnl(i,j,l))
                END DO
              END DO
              IF (nprj_ppnl > 1) THEN
                CALL symmetrize_matrix(hprj_ppnl(:,:,l),"upper_to_lower")
              END IF
              lprj_ppnl_max = MAX(lprj_ppnl_max,l + 2*(nprj_ppnl - 1))
            END DO

            potential%nppnl = nppnl

            CALL init_orbital_pointers(lprj_ppnl_max)

            potential%lprj_ppnl_max = lprj_ppnl_max
            potential%nprj_ppnl_max = nprj_ppnl_max

            CALL reallocate(potential%hprj_ppnl,1,nprj_ppnl_max,&
                                                1,nprj_ppnl_max,&
                                                0,lppnl)
            potential%hprj_ppnl(:,:,:) = hprj_ppnl(:,:,:)

            CALL reallocate(potential%cprj,1,ncoset(lprj_ppnl_max),1,nppnl)
            CALL reallocate(potential%cprj_ppnl,1,nprj_ppnl_max,0,lppnl)
            CALL reallocate(potential%vprj_ppnl,1,nppnl,1,nppnl)

            DEALLOCATE (hprj_ppnl,STAT=istat)
            IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                             "hprj_ppnl")

          END IF

          EXIT search_loop

        END IF

      ELSE

!       *** Stop program, if the end of file is reached ***

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The requested atomic potential <"//&
                          TRIM(potential_name)//&
                          "> was not found in the potential file <"//&
                          TRIM(globenv%potential_file_name)//">",globenv)

      END IF

    END DO search_loop

    CALL finish_parser()

  END SUBROUTINE read_gth_potential

! *****************************************************************************

  SUBROUTINE set_default_all_potential(potential,z)

    TYPE(all_potential_type), POINTER        :: potential
    INTEGER, INTENT(IN)                      :: z

    CHARACTER(LEN=default_string_length)     :: name
    REAL(KIND=dp)                            :: r, alpha
    REAL(KIND=dp)                            :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius
    INTEGER                                  :: zeff
    INTEGER, DIMENSION(:),POINTER            :: elec_conf

    ALLOCATE ( elec_conf(0:3) )
    elec_conf(0:3) = ptable(z)%e_conv(0:3)
    zeff = SUM ( elec_conf )
    name = ptable(z)%name

    r = ptable(z)%covalent_radius * 0.5_dp
    r = MAX ( r, 0.2_dp )
    r = MIN ( r, 1.0_dp )
    alpha = 1.0_dp/(2.0_dp*r**2)

    core_charge_radius = r
    alpha_core_charge = alpha
    ccore_charge = REAL(zeff,dp)*SQRT((alpha/pi)**3)

    CALL set_all_potential(potential,&
                           name=name,&
                           alpha_core_charge=alpha_core_charge,&
                           ccore_charge=ccore_charge,&
                           core_charge_radius=core_charge_radius,&
                           z=z,&
                           zeff=zeff,&
                           elec_conf=elec_conf)
    
    DEALLOCATE ( elec_conf )

  END SUBROUTINE set_default_all_potential

! *****************************************************************************

  SUBROUTINE set_all_potential(potential,name,alpha_core_charge,&
                               ccore_charge,core_charge_radius,z,zeff,&
                               elec_conf)

!   Purpose: Set the components of an all-electron potential data set.

!   History: Creation (11.01.2002,MK)

!   ***************************************************************************

    TYPE(all_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius
    INTEGER, INTENT(IN), OPTIONAL            :: z, zeff
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_all_potential"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      IF (PRESENT(name)) potential%name = name
      IF (PRESENT(alpha_core_charge))&
        potential%alpha_core_charge = alpha_core_charge
      IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
      IF (PRESENT(core_charge_radius))&
        potential%core_charge_radius = core_charge_radius
      IF (PRESENT(z)) potential%z = z
      IF (PRESENT(zeff)) potential%zeff = zeff
      IF (PRESENT(elec_conf)) THEN
        IF ( .NOT. ASSOCIATED( potential%elec_conf ) ) THEN
          CALL reallocate(potential%elec_conf,0,SIZE(elec_conf)-1)
        END IF
        potential%elec_conf(:) = elec_conf(:)
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE set_all_potential

! *****************************************************************************

  SUBROUTINE set_elp_potential(potential,name,alpha_core_charge,ccore_charge,&
                               core_charge_radius,cexp_ppl,cerf_ppl,qeff,&
                               ppl_radius,alpha_ppl,nexp_ppl,nalpha_ppl)

!   Purpose: Set the effective charge in the elp potential data set.
!   History: creation (22-09-2002, GT)

!   ***************************************************************************

    TYPE(elp_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cexp_ppl
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: cerf_ppl, qeff, ppl_radius
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppl
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: nexp_ppl
    INTEGER, OPTIONAL                        :: nalpha_ppl

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_elp_potential"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      IF (PRESENT(name)) potential%name = name
      IF (PRESENT(nalpha_ppl)) potential%nalpha_ppl = nalpha_ppl
      IF (PRESENT(alpha_core_charge))&
        potential%alpha_core_charge = alpha_core_charge
      IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
      IF (PRESENT(core_charge_radius))&
        potential%core_charge_radius = core_charge_radius
      IF (PRESENT(qeff)) potential%qeff = qeff
      IF (PRESENT(alpha_ppl)) potential%alpha_ppl => alpha_ppl
      IF (PRESENT(cexp_ppl)) potential%cexp_ppl => cexp_ppl
      IF (PRESENT(nexp_ppl)) potential%nexp_ppl => nexp_ppl
      IF (PRESENT(cerf_ppl)) potential%cerf_ppl = cerf_ppl
      IF (PRESENT(ppl_radius)) potential%ppl_radius = ppl_radius

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE set_elp_potential

! *****************************************************************************

  SUBROUTINE set_gth_potential(potential,name,alpha_core_charge,alpha_ppl,&
                               ccore_charge,cerf_ppl,core_charge_radius,&
                               ppl_radius,ppnl_radius,lppnl,lprj_ppnl_max,&
                               nexp_ppl,nppnl,nprj_ppnl_max,z,zeff,alpha_ppnl,&
                               cexp_ppl,elec_conf,nprj_ppnl,cprj,cprj_ppnl,&
                               vprj_ppnl,hprj_ppnl)

!   Purpose: Get the components of a GTH potential data set.

!   History: Creation (11.01.2002,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential
    CHARACTER(LEN=default_string_length), &
      INTENT(IN), OPTIONAL                   :: name
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: alpha_core_charge, alpha_ppl, &
      ccore_charge, cerf_ppl, core_charge_radius, ppl_radius, ppnl_radius
    INTEGER, INTENT(IN), OPTIONAL            :: lppnl, lprj_ppnl_max, &
                                                nexp_ppl, nppnl, &
                                                nprj_ppnl_max, z, zeff
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: alpha_ppnl, cexp_ppl
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf, nprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: cprj, cprj_ppnl, vprj_ppnl
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: hprj_ppnl

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_gth_potential"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(potential)) THEN

      IF (PRESENT(name)) potential%name = name
      IF (PRESENT(alpha_core_charge))&
        potential%alpha_core_charge = alpha_core_charge
      IF (PRESENT(alpha_ppl)) potential%alpha_ppl = alpha_ppl
      IF (PRESENT(ccore_charge)) potential%ccore_charge = ccore_charge
      IF (PRESENT(cerf_ppl)) potential%cerf_ppl = cerf_ppl
      IF (PRESENT(core_charge_radius))&
        potential%core_charge_radius = core_charge_radius
      IF (PRESENT(ppl_radius)) potential%ppl_radius = ppl_radius
      IF (PRESENT(ppnl_radius)) potential%ppnl_radius = ppnl_radius
      IF (PRESENT(lppnl)) potential%lppnl = lppnl
      IF (PRESENT(lprj_ppnl_max)) potential%lprj_ppnl_max = lprj_ppnl_max
      IF (PRESENT(nexp_ppl)) potential%nexp_ppl = nexp_ppl
      IF (PRESENT(nppnl)) potential%nppnl = nppnl
      IF (PRESENT(nprj_ppnl_max)) potential%nprj_ppnl_max = nprj_ppnl_max
      IF (PRESENT(z)) potential%z = z
      IF (PRESENT(zeff)) potential%zeff = zeff
      IF (PRESENT(alpha_ppnl)) potential%alpha_ppnl => alpha_ppnl
      IF (PRESENT(cexp_ppl)) potential%cexp_ppl => cexp_ppl
      IF (PRESENT(elec_conf)) potential%elec_conf => elec_conf
      IF (PRESENT(nprj_ppnl)) potential%nprj_ppnl => nprj_ppnl
      IF (PRESENT(cprj)) potential%cprj => cprj
      IF (PRESENT(cprj_ppnl)) potential%cprj_ppnl => cprj_ppnl
      IF (PRESENT(vprj_ppnl)) potential%vprj_ppnl => vprj_ppnl
      IF (PRESENT(hprj_ppnl)) potential%hprj_ppnl => hprj_ppnl

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer potential is not associated")

    END IF

  END SUBROUTINE set_gth_potential

! *****************************************************************************

  SUBROUTINE write_all_potential(potential,globenv)

!   Purpose: Write an atomic all-electron potential data set to the output
!            unit.

!   History: - Creation (09.02.2002,MK)

!   ***************************************************************************

    TYPE(all_potential_type), POINTER        :: potential
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(potential).AND.globenv%print%potential) THEN

      output_unit = globenv%scr

      WRITE (UNIT=output_unit,FMT="(/,T4,A,T57,A24)")&
        "- Potential information",TRIM(potential%name)

      WRITE (UNIT=output_unit,FMT="(T12,A,T41,A40)")&
          "Description: ",TRIM(potential%description(1)),&
          "             ",TRIM(potential%description(2))
      WRITE (UNIT=output_unit,FMT="(T12,A,T69,F12.6)")&
          "Gaussian exponent of the core charge distribution: ",&
          potential%alpha_core_charge
      WRITE (UNIT=output_unit,FMT="(T12,A,T65,4I4)")&
        "Electronic configuration (s p d ...):",potential%elec_conf

    END IF

  END SUBROUTINE write_all_potential

! *****************************************************************************

  SUBROUTINE write_elp_potential(potential,globenv)

!   Purpose: Write an effective-point-charge potential data set to the output
!            unit.

!   History: - Creation (24.09.2002,GT)

!   ***************************************************************************

    TYPE(elp_potential_type), POINTER        :: potential
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, j, output_unit
    REAL(KIND=dp)                            :: r

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(potential).AND.globenv%print%potential) THEN

      output_unit = globenv%scr

      WRITE (UNIT=output_unit,FMT="(/,/,T4,A,/,/,T6,A)")&
        "- Potential information",&
        "Potential name: "//TRIM(potential%name)

      WRITE (UNIT=output_unit,FMT="(/,(T6,A))")&
        "Description: "//TRIM(potential%description(1))
      WRITE (UNIT=output_unit,FMT="(/,T6,A,F12.6)")&
        "Atomic Effective Charge:",potential%qeff
      WRITE (UNIT=output_unit,FMT="(/,T6,A,F12.6)")&
        "Gaussian exponent of the core charge distribution: ",&
        potential%alpha_core_charge

!      IF (potential%nalpha_ppl>0) THEN
!        DO I=1,potential%nalpha_ppl
!           r=sqrt(1.0_dp/(2.0_dp*potential%alpha_ppl(i)))
!           WRITE (UNIT=output_unit,FMT="(/,T12,A,/,/,T8,A,/,T8,6F12.6)")&
!           "Parameters of the local pseudopotential:",&
!           "alpha       rloc        C1          C2          C3          C4",&
!           potential%alpha_ppl(i), r, &
!           (potential%cexp_ppl(i,j)*r**(2*(j-1)),j=1,potential%nexp_ppl(i))
!        END DO
!      END IF

    END IF

  END SUBROUTINE write_elp_potential

! *****************************************************************************

  SUBROUTINE write_gth_potential(potential,globenv)

!   Purpose: Write an atomic GTH potential data set to the output unit.

!   History: - Creation (09.02.2002,MK)

!   ***************************************************************************

    TYPE(gth_potential_type), POINTER        :: potential
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: i, j, l, output_unit
    REAL(KIND=dp)                            :: r

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(potential).AND.globenv%print%potential) THEN

      output_unit = globenv%scr

      WRITE (UNIT=output_unit,FMT="(/,/,T4,A,/,/,T6,A)")&
        "- Potential information",&
        "Potential name: "//TRIM(potential%name)

      WRITE (UNIT=output_unit,FMT="(/,(T6,A))")&
        "Description: "//TRIM(potential%description(1)),&
        "             "//TRIM(potential%description(2)),&
        "             "//TRIM(potential%description(3))
      WRITE (UNIT=output_unit,FMT="(/,T6,A,F12.6)")&
        "Gaussian exponent of the core charge distribution: ",&
        potential%alpha_core_charge
      WRITE (UNIT=output_unit,FMT="(/,T6,A,10I3)")&
        "Electronic configuration (s p d ...):",potential%elec_conf
        r = 1.0_dp/SQRT(2.0_dp*potential%alpha_ppl)

      WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T8,A,/,T8,4F12.6)")&
        "Parameters of the local part of the GTH pseudopotential:",&
        "      rloc        C1          C2          C3          C4",&
        r,(potential%cexp_ppl(i)*r**(2*(i-1)),i=1,potential%nexp_ppl)

      IF (potential%lppnl > -1) THEN
        WRITE (UNIT=output_unit,FMT="(/,T6,A,/,/,T6,A,/)")&
          "Parameters of the non-local part of the GTH pseudopotential:",&
          " l      r(l)      h(i,j,l)"
        DO l=0,potential%lppnl
          r = SQRT(0.5_dp/potential%alpha_ppnl(l))
          WRITE (UNIT=output_unit,FMT="(T6,I2,F12.6,4F12.6)")&
            l,r,(potential%hprj_ppnl(1,j,l),j=1,potential%nprj_ppnl(l))
          DO i=2,potential%nprj_ppnl(l)
            WRITE (UNIT=output_unit,FMT="(T20,4F12.6)")&
              (potential%hprj_ppnl(i,j,l),j=1,potential%nprj_ppnl(l))
          END DO
        END DO
      END IF

    END IF

  END SUBROUTINE write_gth_potential

! *****************************************************************************

END MODULE external_potential_types
