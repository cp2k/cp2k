!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of routines handling the localization for molecular properties
! *****************************************************************************
MODULE qs_loc_molecules
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: pbc
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_minloc,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qs_loc_types,                    ONLY: qs_loc_env_new_type
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

  IMPLICIT NONE

  ! *** Public ***
  PUBLIC :: wfc_to_molecule

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_molecules'

CONTAINS

! *****************************************************************************
!> \brief maps wfc's to molecules and also prints molecular dipoles
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE wfc_to_molecule ( qs_loc_env, loc_print_key, center, molecule_set, error )
    TYPE(qs_loc_env_new_type), INTENT(IN)    :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_print_key
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfc_to_molecule', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: mol_name
    INTEGER :: counter, first_atom, i, iatom, ikind, imol, imol_now, isos, &
      istate, local_location, natom, natom_loc, natom_max, nkind, nmol, &
      nstate, output_unit
    INTEGER, POINTER                         :: wfc_to_atom_map( : )
    REAL(KIND=dp)                            :: dipole(3), dr(3), mydist(2), &
                                                r_pbc( 3 ), zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: charge_set
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dipole_set
    REAL(KIND=dp), POINTER                   :: distance( : ), r( :, : )
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_type), POINTER             :: particle_set( : )

    logger => cp_error_get_logger(error)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    particle_set => qs_loc_env % particle_set
    para_env => qs_loc_env % para_env
    local_molecules => qs_loc_env % local_molecules
    nstate = SIZE ( center, 2 )
    ALLOCATE(wfc_to_atom_map(nstate),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nkind = SIZE (local_molecules%n_el)
    natom = 0
    natom_max = 0
    DO ikind=1,nkind
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol
          i = local_molecules%list(ikind)%array(imol)
          molecule_kind => molecule_set(i)%molecule_kind
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               natom=natom)
          natom_max = natom_max + natom
          molecule_set(i)%lmi%nstates=0
          IF(ASSOCIATED(molecule_set(i)%lmi%states)) THEN
             DEALLOCATE(molecule_set(i)%lmi%states,STAT=isos)
             IF(isos/=0) CALL stop_memory('wfc_to_molecule','states')
          END IF
       END DO
    END DO
    natom_loc = natom_max
    natom = natom_max

    CALL mp_max(natom_max,para_env%group)

    ALLOCATE(r(3,natom_max),STAT=isos)
    IF(isos/=0) CALL stop_memory('centers_and_spreads','r')

    ALLOCATE(distance(natom_max),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')

    !Zero all the stuff
    r(:,:) = 0.0_dp
    distance(:) = 1.E10_dp

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    counter = 0
    nkind = SIZE (local_molecules%n_el)
    DO ikind=1,nkind
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol
          i = local_molecules%list(ikind)%array(imol)
          molecule_kind => molecule_set(i)%molecule_kind
          first_atom = molecule_set(i)%first_atom
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               natom=natom,&
               name=mol_name,&
               atom_list=atom_list)
          DO iatom=1,natom
             counter = counter + 1
             r ( :, counter ) = particle_set(first_atom+iatom-1) % r ( : )
          END DO
       END DO
    END DO

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    DO istate=1,nstate
       distance(:) = 1.E10_dp
       DO iatom=1,natom_loc
          dr(1) = r(1,iatom) - center(1,istate)
          dr(2) = r(2,iatom) - center(2,istate)
          dr(3) = r(3,iatom) - center(3,istate)
          r_pbc = pbc ( dr, qs_loc_env%cell )
          distance(iatom) = SQRT ( DOT_PRODUCT ( r_pbc, r_pbc ) )
       END DO

       !combine distance() from all procs
       local_location = MINLOC(distance,DIM=1)

       mydist(1) = distance(local_location)
       mydist(2) = para_env%mepos

       CALL mp_minloc(mydist,para_env%group)

       IF(mydist(2)==para_env%mepos) THEN
          wfc_to_atom_map(istate)  =  local_location
       ELSE
          wfc_to_atom_map(istate)  =  0
       END IF
    END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    IF(natom_loc/=0) THEN
       DO istate = 1,nstate
          iatom = wfc_to_atom_map(istate)
          IF(iatom/=0) THEN
             counter = 0
             nkind = SIZE (local_molecules%n_el)
             DO ikind=1,nkind
                nmol = SIZE ( local_molecules%list(ikind)%array )
                DO imol=1,nmol
                   imol_now = local_molecules%list(ikind)%array(imol)
                   molecule_kind => molecule_set(imol_now)%molecule_kind
                   CALL get_molecule_kind(molecule_kind=molecule_kind,&
                        natom=natom)
                   counter = counter + natom
                   IF(counter>=iatom) EXIT
                END DO
                IF(counter>=iatom) EXIT
             END DO
             i = molecule_set(imol_now)%lmi%nstates
             i = i + 1
             molecule_set(imol_now)%lmi%nstates = i
             CALL reallocate(molecule_set(imol_now)%lmi%states,1,i)
             molecule_set(imol_now)%lmi%states(i) = istate
          END IF
       END DO
    END IF

    !---------------------------------------------------------------------------
    ! figure out dipole of the molecule.
    ! Assume it is neutral (so origin doesn't matter, we take the first atom)
    ! if the molecule is larger than half the box there might be a problem folding stuff
    ! FIXME this should be moved out of this place
    !---------------------------------------------------------------------------
    particle_set => qs_loc_env % particle_set
    para_env => qs_loc_env % para_env
    local_molecules => qs_loc_env % local_molecules
    nstate = SIZE ( center, 2 )
    ALLOCATE(dipole_set(3,SIZE(molecule_set)))
    ALLOCATE(charge_set(SIZE(molecule_set)))
    dipole_set=0.0_dp
    charge_set=0.0_dp
    nkind = SIZE (local_molecules%n_el)
    DO ikind=1,nkind ! loop over different molecules
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol ! all the molecules of the kind
          imol_now = local_molecules%list(ikind)%array(imol) ! index in the global array
          IF (.NOT.ASSOCIATED(molecule_set(imol_now)%lmi%states)) CYCLE
          molecule_kind => molecule_set(imol_now)%molecule_kind
          first_atom = molecule_set(imol_now)%first_atom
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               natom=natom,atom_list=atom_list)
          dipole=0.0_dp
          ! charges of the atoms involved
          DO iatom=1,natom
             i=first_atom+iatom-1
             atomic_kind=>particle_set(i)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=zeff)
             r_pbc = particle_set(i) % r - particle_set(first_atom) % r
             r_pbc = pbc ( r_pbc , qs_loc_env%cell )
             dipole=dipole + zeff*r_pbc
             charge_set(imol_now)=charge_set(imol_now)+zeff
          END DO
          ! charges of the wfc involved
          DO istate=1,SIZE(molecule_set(imol_now)%lmi%states)
             i = molecule_set(imol_now)%lmi%states(istate)
             r_pbc = center(1:3,i) - particle_set(first_atom) % r
             r_pbc = pbc ( r_pbc , qs_loc_env%cell )
             ! XXXXXXXXXXXXXXXXXXXX double occupation FIXME XXXXXXXXX (where are the spins?)
             dipole=dipole - 2.0_dp * r_pbc 
             charge_set(imol_now)=charge_set(imol_now)-2.0_dp
          ENDDO
          dipole_set(:,imol_now)=dipole ! a.u.
       ENDDO
    ENDDO
    CALL mp_sum(dipole_set,para_env%group)
    CALL mp_sum(charge_set,para_env%group)
    
    output_unit=cp_print_key_unit_nr(logger,loc_print_key,"/MOLECULAR_DIPOLES",&
         extension=".data",middle_name="MOLECULAR_DIPOLES",error=error)
    IF (output_unit>0) THEN
       WRITE(UNIT=output_unit,FMT='(A80)') "# molecule nr, charge, dipole vector, dipole (Debye) "
       dipole_set=dipole_set*debye ! Debye
       DO I=1,SIZE(dipole_set,2)
          WRITE(UNIT=output_unit,FMT='(I6,5F12.6)') I,charge_set(I),dipole_set(1:3,I), &
               SQRT(DOT_PRODUCT(dipole_set(1:3,I),dipole_set(1:3,I)))
       ENDDO
    ENDIF
    CALL cp_print_key_finished_output(output_unit,logger,loc_print_key,&
         "/MOLECULAR_DIPOLES",error=error)
    DEALLOCATE(dipole_set,charge_set)
    !---------------------------------------------------------------------------
    ! end of molecular dipole calculation
    !---------------------------------------------------------------------------

    DEALLOCATE(distance,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')
    DEALLOCATE(r,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','r')

    DEALLOCATE(wfc_to_atom_map,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')

  END SUBROUTINE wfc_to_molecule
  !------------------------------------------------------------------------------

END MODULE qs_loc_molecules

