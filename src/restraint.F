!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****** cp2k/restraint [1.0] *
!!
!!   NAME
!!     restraint
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
MODULE restraint
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE colvar_types,                    ONLY: colvar_counters
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: extended_energy_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             local_colvar_constraint_type,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: restraint_control
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'restraint'

CONTAINS

!!****** restraint/restraint_control [1.0] *
!!
!!   NAME
!!     restraint_control
!!
!!   FUNCTION
!!     Computes restraints
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006 [tlaino]
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE restraint_control( force_env, error )

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ii, ikind, imol, iparticle, iparticle_local, &
      isubsys, k, n3x3con_restraint, n4x6con_restraint, nfixd_restraint, &
      nkind, nmol_per_kind, nparticle_local, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: energy_3x3, energy_4x6, &
                                                energy_colv, energy_fixd, &
                                                rab(3), rab2, targ(3)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_counters)                    :: ncolv
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(extended_energy_type), POINTER      :: extended_energies
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind, &
                                                molecule_kind_set(:)
    TYPE(molecule_type), POINTER             :: molecule, molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), POINTER             :: particle_set( : )

    failure = .FALSE.
    CALL timeset ( routineN, 'I', ' ', handle )
    IF (.NOT. failure) THEN
       CALL force_env_get(force_env=force_env, subsys=subsys, cell=cell)
       energy_4x6  = 0.0_dp
       energy_3x3  = 0.0_dp
       energy_colv = 0.0_dp
       energy_fixd = 0.0_dp
       DO isubsys = 1, SIZE(subsys)
          CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
               particles=particles, &
               molecules_new=molecules, &
               local_particles=local_particles,&
               local_molecules_new=local_molecules,&
               molecule_kinds_new=molecule_kinds )

          nkind             =  molecule_kinds % n_els
          particle_set      => particles%els
          molecule_set      => molecules % els
          molecule_kind_set => molecule_kinds % els
          MOL:  DO ikind = 1, nkind
             ! Fixed Atoms restraint
             
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind ( molecule_kind, nfixd=nfixd_restraint,&
                  fixd_list=fixd_list )

             IF (nfixd_restraint /=0) THEN
                DO i = 1, SIZE(local_particles%n_el)
                   nparticle_local   = local_particles%n_el(i)
                   DO iparticle_local=1,nparticle_local
                      iparticle = local_particles%list(i)%array(iparticle_local)
                      DO ii = 1, SIZE(fixd_list)
                         IF (fixd_list(ii)%fixd == iparticle) THEN
                            IF (fixd_list(ii)%restraint%active) THEN
                               k    = fixd_list(ii)%restraint%k0
                               targ = fixd_list(ii)%coord
                               rab  = particle_set(iparticle)%r-targ
                               rab2 = DOT_PRODUCT(rab,rab)
                               ! Energy
                               energy_fixd = energy_fixd + k * rab2
                               ! Forces
                               particle_set(iparticle)%f(:) =  -2.0_dp * k* rab
                            END IF
                            EXIT
                         END IF
                      END DO
                   END DO
                END DO
             END IF
             ! Other Restraints
             nmol_per_kind = local_molecules % n_el ( ikind )
             DO imol = 1, nmol_per_kind
                i = local_molecules % list ( ikind ) % array ( imol )
                molecule => molecule_set ( i ) 
                molecule_kind => molecule % molecule_kind

                CALL get_molecule_kind ( molecule_kind,&
                                         ncolv = ncolv,  &
                                         ng3x3_restraint = n3x3con_restraint,&
                                         ng4x6_restraint=n4x6con_restraint )
                
                ! 3x3
                IF ( n3x3con_restraint /= 0 ) CALL restraint_3x3( molecule, particle_set,&
                     energy_3x3, error )
                
                ! 4x6
                IF ( n4x6con_restraint /= 0 ) CALL restraint_4x6( molecule, particle_set,&
                     energy_4x6, error )
                
                ! collective variables
                IF ( ncolv%nrestraint /= 0 ) CALL restraint_colv( molecule, particle_set,&
                     cell, energy_colv, error )
             END DO
          END DO MOL
       END DO
       ! Store restraint energies
       NULLIFY(extended_energies)
       ALLOCATE(extended_energies,stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       extended_energies%energy_g3x3_restraint=energy_3x3
       extended_energies%energy_fixd_restraint=energy_fixd
       extended_energies%energy_g4x6_restraint=energy_4x6
       extended_energies%energy_colv_restraint=energy_colv
       CALL force_env_set(force_env=force_env, extended_energies=extended_energies, error=error)
    END IF
    CALL timestop ( 0.0_dp, handle )

  END SUBROUTINE restraint_control

!!****** restraint/restraint_3x3 [1.0] *
!!
!!   NAME
!!     restraint_3x3
!!
!!   FUNCTION
!!     Computes restraints 3x3
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006 [tlaino]
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE restraint_3x3( molecule, particle_set, energy, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_3x3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, ng3x3
    REAL(KIND=dp)                            :: k, rab, rac, rbc, tab, tac, &
                                                tbc
    REAL(KIND=dp), DIMENSION(3)              :: r0_12, r0_13, r0_23
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, g3x3_list = g3x3_list )
    CALL get_molecule ( molecule, first_atom = first_atom )
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
       particle_set ( index_a )%f = particle_set ( index_a )%f -&
            2.0_dp*k*( r0_12/rab*tab+r0_13/rac*tac)
       particle_set ( index_b )%f = particle_set ( index_b )%f -&
            2.0_dp*k*(-r0_12/rab*tab+r0_23/rbc*tbc)
       particle_set ( index_c )%f = particle_set ( index_c )%f -&
            2.0_dp*k*(-r0_13/rac*tac-r0_23/rbc*tbc)
    END DO

  END SUBROUTINE restraint_3x3

!!****** restraint/restraint_4x6 [1.0] *
!!
!!   NAME
!!     restraint_4x6
!!
!!   FUNCTION
!!     Computes restraints 4x6
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006 [tlaino]
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE restraint_4x6( molecule, particle_set, energy, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_4x6', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, iconst, index_a, &
                                                index_b, index_c, index_d, &
                                                ng4x6
    REAL(KIND=dp)                            :: k, rab, rac, rad, rbc, rbd, &
                                                rcd, tab, tac, tad, tbc, tbd, &
                                                tcd
    REAL(KIND=dp), DIMENSION(3)              :: r0_12, r0_13, r0_14, r0_23, &
                                                r0_24, r0_34
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, ng4x6 = ng4x6, g4x6_list = g4x6_list )
    CALL get_molecule ( molecule, first_atom = first_atom )
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
       particle_set ( index_a )%f = particle_set ( index_a )%f-&
            2.0_dp*k*( r0_12/rab*tab+r0_13/rac*tac+r0_14/rad*tad)
       particle_set ( index_b )%f = particle_set ( index_b )%f-&
            2.0_dp*k*(-r0_12/rab*tab+r0_23/rbc*tbc+r0_24/rbd*tbd)
       particle_set ( index_c )%f = particle_set ( index_c )%f-&
            2.0_dp*k*(-r0_13/rac*tac-r0_23/rbc*tbc+r0_34/rcd*tcd)
       particle_set ( index_d )%f = particle_set ( index_d )%f-&
            2.0_dp*k*(-r0_14/rad*tad-r0_24/rbd*tbd-r0_34/rcd*tcd)       
    END DO

  END SUBROUTINE restraint_4x6

!!****** restraint/restraint_colv [1.0] *
!!
!!   NAME
!!     restraint_colv
!!
!!   FUNCTION
!!     Computes restraints colv
!!
!!   AUTHOR
!!     Teodoro Laino 08.2006 [tlaino]
!!
!!   MODIFICATION HISTORY
!!     
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE restraint_colv( molecule, particle_set, cell, energy, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restraint_colv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, iatm, iconst, ind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: k, tab, targ
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    failure = .FALSE.
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list = colv_list )
    CALL get_molecule ( molecule, first_atom = first_atom, lcolv=lcolv )
    DO iconst = 1, SIZE(colv_list)
       IF (.NOT.colv_list(iconst)%restraint%active) CYCLE
       ! Update colvar
       CALL colvar_eval_mol_f( lcolv ( iconst ) % colvar, cell, &
            particles=particle_set, error=error)

       k    = colv_list ( iconst )%restraint%k0
       targ = colv_list ( iconst )%expected_value
       tab  = lcolv ( iconst ) % colvar % ss - targ
       ! Update Energy
       energy = energy + k * tab**2
       ! Update Forces
       DO iatm = 1, SIZE( lcolv (iconst) % colvar % i_atom )
          ind = lcolv(iconst) % colvar % i_atom (iatm)
          particle_set(ind)%f = particle_set(ind)%f &
               -2.0_dp*k*tab*lcolv(iconst) % colvar % dsdr(:,iatm)
       END DO
    END DO

  END SUBROUTINE restraint_colv

END MODULE restraint
