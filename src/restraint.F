!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Handles all possible kinds of restraints in CP2K
!> \author Teodoro Laino 08.2006 [tlaino] - University of Zurich
!> \par    History
!>         Teodoro Laino [tlaino] - 11.2008 : Improved the fixd_list restraints
! *****************************************************************************
MODULE restraint
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE colvar_types,                    ONLY: colvar_counters,&
                                             diff_colvar
  USE constraint_fxd,                  ONLY: create_local_fixd_list,&
                                             release_local_fixd_list
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set,&
                                             force_env_type
  USE input_constants,                 ONLY: use_perd_x,&
                                             use_perd_xy,&
                                             use_perd_xyz,&
                                             use_perd_xz,&
                                             use_perd_y,&
                                             use_perd_yz,&
                                             use_perd_z
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             local_fixd_constraint_type,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             local_colvar_constraint_type,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             update_particle_set
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: restraint_control
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'restraint'

CONTAINS

! *****************************************************************************
!> \brief Computes restraints
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_control( force_env, error )

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ifixd, ii, ikind, imol, iparticle, &
      n3x3con_restraint, n4x6con_restraint, n_restraint, nkind, &
      nmol_per_kind, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: energy_3x3, energy_4x6, energy_colv, energy_fixd, &
      extended_energies, k0, rab(3), rab2, targ(3)
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: force
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_counters)                    :: ncolv
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(local_fixd_constraint_type), &
      POINTER                                :: lfixd_list(:)
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind, &
                                                molecule_kind_set(:)
    TYPE(molecule_type), POINTER             :: molecule, molecule_set(:)
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), POINTER             :: particle_set(:)

    failure = .FALSE.
    NULLIFY(cell, subsys, local_molecules, local_particles, fixd_list, molecule_kinds,&
         molecules, molecule_kind, molecule_kind_set, molecule, molecule_set, particles,&
         particle_set, gci, lfixd_list)
    CALL timeset(routineN,handle)
    IF (.NOT. failure) THEN
       CALL force_env_get(force_env=force_env, subsys=subsys, cell=cell,error=error)
       energy_4x6  = 0.0_dp
       energy_3x3  = 0.0_dp
       energy_colv = 0.0_dp
       energy_fixd = 0.0_dp
       n_restraint=0
       CALL cp_subsys_get(subsys=subsys,particles=particles, molecules_new=molecules,&
            local_particles=local_particles, local_molecules_new=local_molecules,&
            gci=gci, molecule_kinds_new=molecule_kinds ,error=error)

       nkind             =  molecule_kinds%n_els
       particle_set      => particles%els
       molecule_set      => molecules%els
       molecule_kind_set => molecule_kinds%els

       ! Intramolecular Restraints
       ALLOCATE(force(3,SIZE(particle_set)),stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       force = 0.0_dp

       ! Create the list of locally fixed atoms
       CALL create_local_fixd_list(lfixd_list, nkind, molecule_kind_set, local_particles,&
            error)

       DO ifixd = 1, SIZE(lfixd_list)
          ikind = lfixd_list(ifixd)%ikind
          ii    = lfixd_list(ifixd)%ifixd_index
          molecule_kind => molecule_kind_set(ikind)
          CALL get_molecule_kind(molecule_kind, fixd_list=fixd_list)
          IF (fixd_list(ii)%restraint%active) THEN
             n_restraint = n_restraint+1
             iparticle   = fixd_list(ii)%fixd
             k0          = fixd_list(ii)%restraint%k0
             targ        = fixd_list(ii)%coord
             rab         = 0.0_dp
             SELECT CASE(fixd_list(ii)%itype)
             CASE(use_perd_x)
                rab(1)  = particle_set(iparticle)%r(1)-targ(1)
             CASE(use_perd_y)
                rab(2)  = particle_set(iparticle)%r(2)-targ(2)
             CASE(use_perd_z)
                rab(3)  = particle_set(iparticle)%r(3)-targ(3)
             CASE(use_perd_xy)
                rab(1)  = particle_set(iparticle)%r(1)-targ(1)
                rab(2)  = particle_set(iparticle)%r(2)-targ(2)
             CASE(use_perd_xz)
                rab(1)  = particle_set(iparticle)%r(1)-targ(1)
                rab(3)  = particle_set(iparticle)%r(3)-targ(3)
             CASE(use_perd_yz)
                rab(2)  = particle_set(iparticle)%r(2)-targ(2)
                rab(3)  = particle_set(iparticle)%r(3)-targ(3)
             CASE(use_perd_xyz)
                rab  = particle_set(iparticle)%r-targ
             END SELECT
             rab2 = DOT_PRODUCT(rab,rab)
             ! Energy
             energy_fixd = energy_fixd + k0 * rab2
             ! Forces
             force(:,iparticle) = force(:,iparticle) - 2.0_dp * k0 * rab
          END IF
       END DO
       CALL release_local_fixd_list(lfixd_list, error)

       ! Loop over other kind of Restraints
       MOL:  DO ikind = 1, nkind
          molecule_kind => molecule_kind_set(ikind)
          nmol_per_kind = local_molecules % n_el ( ikind )
          DO imol = 1, nmol_per_kind
             i = local_molecules % list ( ikind ) % array ( imol )
             molecule => molecule_set ( i )
             molecule_kind => molecule % molecule_kind

             CALL get_molecule_kind ( molecule_kind,&
                  ncolv = ncolv,  &
                  ng3x3_restraint = n3x3con_restraint,&
                  ng4x6_restraint = n4x6con_restraint )
             ! 3x3
             IF ( n3x3con_restraint /= 0 ) THEN
                n_restraint=n_restraint+n3x3con_restraint
                CALL restraint_3x3_int( molecule, particle_set, energy_3x3, force, error )
             ENDIF
             ! 4x6
             IF ( n4x6con_restraint /= 0 ) THEN
                n_restraint=n_restraint+n4x6con_restraint
                CALL restraint_4x6_int( molecule, particle_set, energy_4x6, force, error )
             ENDIF
             ! collective variables
             IF ( ncolv%nrestraint /= 0 ) THEN
                n_restraint=n_restraint+ncolv%nrestraint
                CALL restraint_colv_int( molecule, particle_set, cell, energy_colv, force,&
                     error )
             ENDIF
          END DO
       END DO MOL
       CALL mp_sum(n_restraint,force_env%para_env%group)
       IF (n_restraint>0) THEN
          CALL mp_sum(energy_fixd,force_env%para_env%group)
          CALL mp_sum(energy_3x3,force_env%para_env%group)
          CALL mp_sum(energy_4x6,force_env%para_env%group)
          CALL mp_sum(energy_colv,force_env%para_env%group)
          CALL update_particle_set(particle_set,force_env%para_env%group,for=force,add=.TRUE.,&
               error=error)
          force = 0.0_dp
          n_restraint = 0
       ENDIF
       ! Intermolecular Restraints
       IF (ASSOCIATED(gci)) THEN
          IF (gci%nrestraint>0) THEN
             ! 3x3
             IF ( gci%ng3x3_restraint /= 0 ) THEN
                n_restraint=n_restraint+gci%ng3x3_restraint
                CALL restraint_3x3_ext( gci, particle_set, energy_3x3, force, error )
             ENDIF
             ! 4x6
             IF ( gci%ng4x6_restraint /= 0 ) THEN
                n_restraint=n_restraint+gci%ng4x6_restraint
                CALL restraint_4x6_ext( gci, particle_set, energy_4x6, force, error )
             ENDIF
             ! collective variables
             IF ( gci%ncolv%nrestraint /= 0 ) THEN
                n_restraint=n_restraint+gci%ncolv%nrestraint
                CALL restraint_colv_ext( gci, particle_set, cell, energy_colv, force,&
                     error )
             ENDIF
             DO iparticle = 1, SIZE(particle_set)
                particle_set(iparticle)%f = particle_set(iparticle)%f + force(:,iparticle)
             END DO
          END IF
       END IF
       DEALLOCATE(force,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Store restraint energies
       CALL force_env_get(force_env=force_env, additional_potential=extended_energies, error=error)
       extended_energies=extended_energies+ energy_3x3  +&
                                            energy_fixd +&
                                            energy_4x6  +&
                                            energy_colv
       CALL force_env_set(force_env=force_env, additional_potential=extended_energies, error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE restraint_control

! *****************************************************************************
!> \brief Computes restraints 3x3 - Intramolecular
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_3x3_int( molecule, particle_set, energy, force, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_3x3_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, ng3x3
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list,&
         fixd_list = fixd_list)
    CALL get_molecule ( molecule, first_atom = first_atom )
    CALL restraint_3x3_low( ng3x3, g3x3_list, fixd_list, first_atom, particle_set,&
         energy, force, error )

  END SUBROUTINE restraint_3x3_int

! *****************************************************************************
!> \brief Computes restraints 4x6 - Intramolecular
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_4x6_int( molecule, particle_set, energy, force, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_4x6_int', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, ng4x6
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list,&
         fixd_list=fixd_list )
    CALL get_molecule ( molecule, first_atom = first_atom )
    CALL restraint_4x6_low( ng4x6, g4x6_list, fixd_list, first_atom, particle_set,&
         energy, force, error )

  END SUBROUTINE restraint_4x6_int

! *****************************************************************************
!> \brief Computes restraints colv - Intramolecular
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_colv_int( molecule, particle_set, cell, energy, force, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_colv_int', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)

    failure = .FALSE.
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, lcolv=lcolv )
    CALL restraint_colv_low( colv_list, fixd_list, lcolv, particle_set,&
         cell, energy, force, error )

  END SUBROUTINE restraint_colv_int

! *****************************************************************************
!> \brief Computes restraints 3x3 - Intermolecular
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_3x3_ext( gci, particle_set, energy, force, error )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_3x3_ext', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, ng3x3
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )

    first_atom = 1
    ng3x3     =  gci%ng3x3
    g3x3_list => gci%g3x3_list
    fixd_list => gci%fixd_list
    CALL restraint_3x3_low( ng3x3, g3x3_list, fixd_list, first_atom, particle_set,&
         energy, force, error )

  END SUBROUTINE restraint_3x3_ext

! *****************************************************************************
!> \brief Computes restraints 4x6 - Intermolecular
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_4x6_ext( gci, particle_set, energy, force, error )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_4x6_ext', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, ng4x6
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )

    first_atom = 1
    ng4x6     =  gci%ng4x6
    g4x6_list => gci%g4x6_list
    fixd_list => gci%fixd_list
    CALL restraint_4x6_low( ng4x6, g4x6_list, fixd_list, first_atom, particle_set,&
         energy, force, error )

  END SUBROUTINE restraint_4x6_ext

! *****************************************************************************
!> \brief Computes restraints colv - Intermolecular
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_colv_ext( gci, particle_set, cell, energy, force, error )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_colv_ext', &
      routineP = moduleN//':'//routineN

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )

    colv_list => gci%colv_list
    fixd_list => gci%fixd_list
    lcolv => gci%lcolv
    CALL restraint_colv_low( colv_list, fixd_list, lcolv, particle_set,&
         cell, energy, force, error )

  END SUBROUTINE restraint_colv_ext

! *****************************************************************************
!> \brief Computes restraints 3x3 - Real 3x3 restraints
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_3x3_low( ng3x3, g3x3_list, fixd_list, first_atom,&
       particle_set, energy, force, error )

    INTEGER                                  :: ng3x3
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    INTEGER, INTENT(IN)                      :: first_atom
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_3x3_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c
    REAL(KIND=dp)                            :: k, rab, rac, rbc, tab, tac, &
                                                tbc
    REAL(KIND=dp), DIMENSION(3)              :: r0_12, r0_13, r0_23

    DO iconst = 1, ng3x3
       IF (.NOT.g3x3_list(iconst)%restraint%active) CYCLE
       index_a = g3x3_list ( iconst ) % a + first_atom -1
       index_b = g3x3_list ( iconst ) % b + first_atom -1
       index_c = g3x3_list ( iconst ) % c + first_atom -1
       r0_12 ( : ) = particle_set ( index_a )%r - particle_set ( index_b )%r
       r0_13 ( : ) = particle_set ( index_a )%r - particle_set ( index_c )%r
       r0_23 ( : ) = particle_set ( index_b )%r - particle_set ( index_c )%r

       rab = SQRT(DOT_PRODUCT(r0_12,r0_12))
       rac = SQRT(DOT_PRODUCT(r0_13,r0_13))
       rbc = SQRT(DOT_PRODUCT(r0_23,r0_23))
       tab = rab-g3x3_list(ng3x3)%dab
       tac = rac-g3x3_list(ng3x3)%dac
       tbc = rbc-g3x3_list(ng3x3)%dbc
       k   = g3x3_list(iconst)%restraint%k0
       ! Update Energy
       energy = energy + k *(tab**2+tac**2+tbc**2)
       ! Update Forces
       force(:,index_a) = force(:,index_a) - 2.0_dp*k*( r0_12/rab*tab+r0_13/rac*tac)
       force(:,index_b) = force(:,index_b) - 2.0_dp*k*(-r0_12/rab*tab+r0_23/rbc*tbc)
       force(:,index_c) = force(:,index_c) - 2.0_dp*k*(-r0_13/rac*tac-r0_23/rbc*tbc)
       ! Fixed atoms
       IF (ASSOCIATED(fixd_list)) THEN
          IF (SIZE(fixd_list)>0) THEN
             IF (ANY(fixd_list(:)%fixd==index_a)) force(:,index_a) = 0.0_dp
             IF (ANY(fixd_list(:)%fixd==index_b)) force(:,index_b) = 0.0_dp
             IF (ANY(fixd_list(:)%fixd==index_c)) force(:,index_c) = 0.0_dp
          END IF
       END IF
    END DO

  END SUBROUTINE restraint_3x3_low

! *****************************************************************************
!> \brief Computes restraints 4x6 - Real 4x6 restraints
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_4x6_low( ng4x6, g4x6_list, fixd_list, first_atom,&
       particle_set, energy, force, error )

    INTEGER                                  :: ng4x6
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    INTEGER, INTENT(IN)                      :: first_atom
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_4x6_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c, index_d
    REAL(KIND=dp)                            :: k, rab, rac, rad, rbc, rbd, &
                                                rcd, tab, tac, tad, tbc, tbd, &
                                                tcd
    REAL(KIND=dp), DIMENSION(3)              :: r0_12, r0_13, r0_14, r0_23, &
                                                r0_24, r0_34

    DO iconst = 1, ng4x6
       IF (.NOT.g4x6_list(iconst)%restraint%active) CYCLE
       index_a = g4x6_list ( iconst ) % a + first_atom -1
       index_b = g4x6_list ( iconst ) % b + first_atom -1
       index_c = g4x6_list ( iconst ) % c + first_atom -1
       index_d = g4x6_list ( iconst ) % d + first_atom -1
       r0_12 ( : ) = particle_set( index_a )%r - particle_set( index_b )%r
       r0_13 ( : ) = particle_set( index_a )%r - particle_set( index_c )%r
       r0_14 ( : ) = particle_set( index_a )%r - particle_set( index_d )%r
       r0_23 ( : ) = particle_set( index_b )%r - particle_set( index_c )%r
       r0_24 ( : ) = particle_set( index_b )%r - particle_set( index_d )%r
       r0_34 ( : ) = particle_set( index_c )%r - particle_set( index_d )%r

       rab = SQRT(DOT_PRODUCT(r0_12,r0_12))
       rac = SQRT(DOT_PRODUCT(r0_13,r0_13))
       rad = SQRT(DOT_PRODUCT(r0_14,r0_14))
       rbc = SQRT(DOT_PRODUCT(r0_23,r0_23))
       rbd = SQRT(DOT_PRODUCT(r0_24,r0_24))
       rcd = SQRT(DOT_PRODUCT(r0_34,r0_34))

       tab = rab - g4x6_list(ng4x6)%dab
       tac = rac - g4x6_list(ng4x6)%dac
       tad = rad - g4x6_list(ng4x6)%dad
       tbc = rbc - g4x6_list(ng4x6)%dbc
       tbd = rbd - g4x6_list(ng4x6)%dbd
       tcd = rcd - g4x6_list(ng4x6)%dcd

       k   = g4x6_list(iconst)%restraint%k0
       ! Update Energy
       energy = energy + k *(tab**2+tac**2+tad**2+tbc**2+tbd**2+tcd**2)
       ! Update Forces
       force(:,index_a) = force(:,index_a)-2.0_dp*k*( r0_12/rab*tab+r0_13/rac*tac+r0_14/rad*tad)
       force(:,index_b) = force(:,index_b)-2.0_dp*k*(-r0_12/rab*tab+r0_23/rbc*tbc+r0_24/rbd*tbd)
       force(:,index_c) = force(:,index_c)-2.0_dp*k*(-r0_13/rac*tac-r0_23/rbc*tbc+r0_34/rcd*tcd)
       force(:,index_d) = force(:,index_d)-2.0_dp*k*(-r0_14/rad*tad-r0_24/rbd*tbd-r0_34/rcd*tcd)
       ! Fixed atoms
       IF (ASSOCIATED(fixd_list)) THEN
          IF (SIZE(fixd_list)>0) THEN
             IF (ANY(fixd_list(:)%fixd==index_a)) force(:,index_a) = 0.0_dp
             IF (ANY(fixd_list(:)%fixd==index_b)) force(:,index_b) = 0.0_dp
             IF (ANY(fixd_list(:)%fixd==index_c)) force(:,index_c) = 0.0_dp
             IF (ANY(fixd_list(:)%fixd==index_d)) force(:,index_d) = 0.0_dp
          END IF
       END IF
    END DO

  END SUBROUTINE restraint_4x6_low

! *****************************************************************************
!> \brief Computes restraints colv - Real COLVAR restraints
!> \author Teodoro Laino 08.2006 [tlaino]
! *****************************************************************************
  SUBROUTINE restraint_colv_low( colv_list, fixd_list, lcolv, &
       particle_set, cell, energy, force, error )

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_colv_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatm, iconst, ind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: k, tab, targ

    failure = .FALSE.
    DO iconst = 1, SIZE(colv_list)
       IF (.NOT.colv_list(iconst)%restraint%active) CYCLE
       ! Update colvar
       CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, &
            particles=particle_set, fixd_list=fixd_list, error=error)

       k    = colv_list ( iconst )%restraint%k0
       targ = colv_list ( iconst )%expected_value
       tab  = diff_colvar(lcolv(iconst)%colvar,targ)
       ! Update Energy
       energy = energy + k * tab**2
       ! Update Forces
       DO iatm = 1, SIZE( lcolv (iconst) % colvar % i_atom )
          ind = lcolv(iconst) % colvar % i_atom (iatm)
          force(:,ind) = force(:,ind)-2.0_dp*k*tab*lcolv(iconst) % colvar % dsdr(:,iatm)
       END DO
    END DO

  END SUBROUTINE restraint_colv_low

END MODULE restraint
