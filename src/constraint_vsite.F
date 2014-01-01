!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to handle the virtual site constraint/restraint
!> \par History
!>      Teodoro Laino [tlaino] 12.2008 - Preparing for VIRTUAL SITE constraints
!>                                       (patch by Marcel Baer)
! *****************************************************************************
MODULE constraint_vsite
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             vsite_constraint_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: shake_vsite_int,&
            force_vsite_int,&
            force_vsite_ext,&
            shake_vsite_ext,&
            vsite_force_control

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_vsite'

CONTAINS

! *****************************************************************************
!> \brief control force distribution for virtual sites
!> \par History
!>      - none
!> \author Marcel Baer
!> \date 12.2008
! *****************************************************************************
  SUBROUTINE vsite_force_control( force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, ikind, imol, nconstraint, &
                                                nkind, nmol_per_kind, &
                                                nvsitecon
    LOGICAL                                  :: do_ext_constraint
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(gci, subsys, local_molecules, local_particles,&
         molecule_kinds)

    CALL force_env_get(force_env=force_env, subsys=subsys,error=error)

    CALL cp_subsys_get(subsys=subsys,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,&
         molecule_kinds_new=molecule_kinds,gci=gci,molecules_new=molecules,&
         error=error)

    molecule_kind_set => molecule_kinds%els
    molecule_set => molecules%els
    particle_set => particles%els
    nkind = SIZE ( molecule_kind_set )
    ! Intermolecular Virtual Site Constraints
    do_ext_constraint = .FALSE.
    IF (ASSOCIATED(gci)) THEN
       do_ext_constraint = (gci%ntot/=0)
    END IF
    ! Intramolecular Virtual Site Constraints
    MOL:  DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i )
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind,nconstraint=nconstraint,nvsite=nvsitecon)
          IF (nconstraint==0) CYCLE
          IF (nvsitecon/=0) &
               CALL force_vsite_int( molecule, particle_set)
       END DO
    END DO MOL
    ! Intermolecular Virtual Site Constraints
    IF (do_ext_constraint) THEN
       IF ( gci%nvsite /= 0 ) &
            CALL force_vsite_ext( gci, particle_set)
    END IF

  END SUBROUTINE vsite_force_control

! *****************************************************************************
!> \brief Intramolecular virtual site
!> \par History
!>      12.2008 Marcel Baer
! *****************************************************************************
  SUBROUTINE shake_vsite_int( molecule, pos, error)
    TYPE(molecule_type), POINTER             :: molecule
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: first_atom, nvsite
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(vsite_constraint_type), POINTER     :: vsite_list( : )

    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, nvsite = nvsite, vsite_list = vsite_list )
    CALL get_molecule ( molecule, first_atom = first_atom)
    ! Real Shake
    CALL shake_vsite_low( vsite_list, nvsite, first_atom, pos, error)

  END SUBROUTINE shake_vsite_int

! *****************************************************************************
!> \brief Intramolecular virtual site
!> \par History
!>      12.2008 Marcel Baer
! *****************************************************************************
  SUBROUTINE shake_vsite_ext( gci, pos, error)

    TYPE(global_constraint_type), POINTER    :: gci
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: first_atom, nvsite
    TYPE(vsite_constraint_type), POINTER     :: vsite_list( : )

    first_atom = 1
    nvsite = gci%nvsite
    vsite_list => gci%vsite_list
    ! Real Shake
    CALL shake_vsite_low( vsite_list, nvsite, first_atom, pos, error)

  END SUBROUTINE shake_vsite_ext

! *****************************************************************************
!> \par History
!>      12.2008 Marcel Bear
! *****************************************************************************
  SUBROUTINE shake_vsite_low(vsite_list, nvsite, first_atom, pos, error)
    TYPE(vsite_constraint_type)              :: vsite_list( : )
    INTEGER, INTENT(IN)                      :: nvsite, first_atom
    REAL(KIND=dp), INTENT(INOUT)             :: pos( :, : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c, index_d
    REAL(KIND=dp), DIMENSION(3)              :: r1, r2

    DO iconst = 1, nvsite
       IF (vsite_list(iconst)%restraint%active) CYCLE
       index_a = vsite_list(iconst)%a + first_atom -1
       index_b = vsite_list(iconst)%b + first_atom -1
       index_c = vsite_list(iconst)%c + first_atom -1
       index_d = vsite_list(iconst)%d + first_atom -1

       r1(:) = pos(:,index_b) - pos(:,index_c )
       r2(:) = pos(:,index_d) - pos(:,index_c )
       pos(:,index_a)=pos(:,index_c)+vsite_list(iconst)%wbc*r1(:)+&
                                     vsite_list(iconst)%wdc*r2(:)
    END DO
  END SUBROUTINE shake_vsite_low

! *****************************************************************************
!> \brief Intramolecular virtual site
!> \par History
!>      12.2008 Marcel Bear
! *****************************************************************************
  SUBROUTINE force_vsite_int( molecule, particle_set)
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, index_d, &
                                                nvsite
    REAL(KIND=dp)                            :: wb, wc, wd
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(vsite_constraint_type), POINTER     :: vsite_list( : )

    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, nvsite = nvsite, vsite_list = vsite_list )
    CALL get_molecule ( molecule, first_atom = first_atom)

    DO iconst = 1, nvsite
       IF (vsite_list(iconst)%restraint%active) CYCLE
       index_a = vsite_list(iconst)%a + first_atom -1
       index_b = vsite_list(iconst)%b + first_atom -1
       index_c = vsite_list(iconst)%c + first_atom -1
       index_d = vsite_list(iconst)%d + first_atom -1

       wb=vsite_list(iconst)%wbc
       wd=vsite_list(iconst)%wdc
       wc=1.0_dp-vsite_list(iconst)%wbc-vsite_list(iconst)%wdc

       particle_set(index_b)%f(:)=particle_set(index_b)%f(:) + wb*particle_set(index_a)%f(:)
       particle_set(index_c)%f(:)=particle_set(index_c)%f(:) + wc*particle_set(index_a)%f(:)
       particle_set(index_d)%f(:)=particle_set(index_d)%f(:) + wd*particle_set(index_a)%f(:)
       particle_set(index_a)%f(:)=0.0_dp
    END DO

  END SUBROUTINE force_vsite_int

! *****************************************************************************
!> \brief Intramolecular virtual site
!> \par History
!>      12.2008 Marcel Bear
! *****************************************************************************
  SUBROUTINE force_vsite_ext( gci, particle_set)
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, index_d, &
                                                nvsite
    REAL(KIND=dp)                            :: wb, wc, wd
    TYPE(vsite_constraint_type), POINTER     :: vsite_list( : )

    first_atom = 1
    nvsite = gci%nvsite
    vsite_list => gci%vsite_list
    ! Real Shake

    DO iconst = 1, nvsite
       IF (vsite_list(iconst)%restraint%active) CYCLE
       index_a = vsite_list(iconst)%a + first_atom -1
       index_b = vsite_list(iconst)%b + first_atom -1
       index_c = vsite_list(iconst)%c + first_atom -1
       index_d = vsite_list(iconst)%d + first_atom -1

       wb=vsite_list(iconst)%wbc
       wd=vsite_list(iconst)%wdc
       wc=1.0_dp-vsite_list(iconst)%wbc-vsite_list(iconst)%wdc

       particle_set(index_b)%f(:)=particle_set(index_b)%f(:) + wb*particle_set(index_a)%f(:)
       particle_set(index_c)%f(:)=particle_set(index_c)%f(:) + wc*particle_set(index_a)%f(:)
       particle_set(index_d)%f(:)=particle_set(index_d)%f(:) + wd*particle_set(index_a)%f(:)
       particle_set(index_a)%f(:)=0.0_dp
    END DO
  END SUBROUTINE force_vsite_ext

END MODULE constraint_vsite
