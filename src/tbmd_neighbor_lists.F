!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_neighbor_lists [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_lists
!!
!!   FUNCTION
!!     Just another neighbor list, but this time for tight-binding and
!!     by making heavily use of the fist neighbor lists.
!!
!!   NOTES
!!     Linear scaling linked-list (cell-list) neighbor-list with verlet
!!     update.
!!
!!   REFERENCES
!!     - Hockney, R.W. & Eastwood, J.W.: Computer Simulation using Particles;
!!                                       McGraw-Hill, New York, 1981
!!     - Verlet, L.: Phys. Rev. 1967, 159, 88
!!     - Auerbach, D.J., et.al.: A Special Purpose Parallel Computer for 
!!                               Molecular Dynamics: Motivation, Design, 
!!                               Implementation, and Application;
!!                               J. Phys. Chem., 1987, 91, 4881-4890
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_neighbor_lists
  USE atomic_kind_types,                      ONLY: atomic_kind_type, &
                                                    get_atomic_kind_set, &
                                                    get_atomic_kind
  USE cell_types,                             ONLY: cell_type, &
                                                    real_to_scaled, &
                                                    scaled_to_real, &
                                                    pbc
  USE cp_output_handling,                     ONLY: cp_print_key_unit_nr, &
                                                    cp_print_key_should_output, &
                                                    cp_p_file, &
                                                    cp_print_key_finished_output
  USE distribution_1d_types,                  ONLY: distribution_1d_type
  USE fist_neighbor_list_types,               ONLY: neighbor_list_set_p_type
  USE fist_neighbor_lists,                    ONLY: build_neighbor_lists, &
                                                    write_neighbor_lists, &
                                                    write_distribution
  USE global_types,                           ONLY: global_environment_type
  USE input_section_types,                    ONLY: section_vals_type
  USE kinds,                                  ONLY: dp, &
                                                    dp_size, &
                                                    int_size
  USE particle_types,                         ONLY: particle_type
  USE qs_neighbor_lists,                      ONLY: local_atoms_type
  USE tbmd_neighbor_list_env_types,           ONLY: tbmd_neighbor_list_env_type, &
                                                    tbmd_neighbor_list_env_get, &
                                                    tbmd_neighbor_list_env_set
  USE termination,                            ONLY: stop_memory
  USE timings,                                ONLY: timeset, &
                                                    timestop

#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: moduleN = "tbmd_neighbor_lists"

! *** Public subroutines ***

  PUBLIC :: tbmd_neighbor_list_control, &
            build_tbmd_neighbor_lists

!!***
! *****************************************************************************

CONTAINS

!!****f* tbmd_neighbor_lists/tbmd_neighbor_list_control [1.0] *
!!
!!   NAME
!!     tbmd_neighbor_list_control
!!
!!   FUNCTION
!!     Verlet update logic, which decides if an update of the neighbor lists 
!!     has to be performed or not
!!
!!   NOTES
!!     An update of the neighbor list is performed if any particle has been
!!     moved by more than 0.5*verlet_skin from it's previous update position.
!!
!!   INPUTS
!!     - atomic_kind_set: The set of all involvolved atomic kinds
!!     - particle_set: The set of particles
!!     - local_particles: Particles which are local on this particluar node
!!     - cell: The simulation cell
!!     - tbmd_neighbor_list_env: The TBMD neighbor-list environment
!!     - force_env_section: Pointer to the force environment input section
!!     - globenv: The global environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_neighbor_list_control(atomic_kind_set, particle_set, &
                                        local_particles, cell, &
                                        tbmd_neighbor_list_env, globenv, &
                                        force_env_section, error)

    TYPE(atomic_kind_type), DIMENSION(:), POINTER     :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER        :: particle_set
    TYPE(distribution_1d_type), POINTER               :: local_particles
    TYPE(cell_type), POINTER                          :: cell
    TYPE(tbmd_neighbor_list_env_type), POINTER        :: tbmd_neighbor_list_env
    TYPE(global_environment_type), &
      pointer                                         :: globenv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                                     :: error
    TYPE(section_vals_type), POINTER :: force_env_section

    CHARACTER(LEN=*), PARAMETER :: routineN = "tbmd_neighbor_list_control"

    REAL(KIND=dp)                                     :: dr2, &
                                                         dr2_max, &
                                                         verlet_skin, &
                                                         aup, &
                                                         lup
    REAL(KIND=dp), DIMENSION(3)                       :: dr
    REAL(KIND=dp), DIMENSION(:, :), POINTER           :: r_last_update, &
                                                         r_cutoff
    INTEGER                                           :: istat, &
                                                         iparticle, &
                                                         iparticle_kind, &
                                                         iparticle_local, &
                                                         nparticle, &
                                                         nparticle_kind, &
                                                         nparticle_local, &
                                                         handle, &
                                                         counter, &
                                                         num_update, &
                                                         last_update
    LOGICAL                                           :: build_from_scratch, &
                                                         update_neighbor_lists
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                           :: tbmd_neighbor_list

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(r_last_update, tbmd_neighbor_list)

! *** Assigning local pointers ***
  CALL tbmd_neighbor_list_env_get(tbmd_neighbor_list_env, &
                                  tbmd_neighbor_list=tbmd_neighbor_list, &
                                  r_cutoff=r_cutoff, &
                                  aup=aup, lup=lup, &
                                  counter=counter, &
                                  r_last_update=r_last_update, &
                                  num_update=num_update, &
                                  verlet_skin=verlet_skin,&
                                  last_update=last_update)

    nparticle_kind = SIZE(atomic_kind_set)
    nparticle = SIZE(particle_set)

! *** Check, if the neighbor lists have to be built or updated ***
    update_neighbor_lists = .FALSE.
    build_from_scratch = .FALSE.

    IF (ASSOCIATED(r_last_update)) THEN
      dr2_max = 0.0_dp

      DO iparticle_kind = 1,nparticle_kind
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local = 1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          dr(:) = scaled_to_real(r_last_update(:,iparticle),cell) - &
                  particle_set(iparticle)%r(:)
          dr2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
          dr2_max = MAX(dr2_max,dr2)
        END DO
      END DO

#if defined(__parallel)
      CALL mp_max(dr2_max, globenv%group)
#endif

      IF (dr2_max > 0.25_dp*verlet_skin**2) THEN
        DO iparticle = 1,nparticle
          r_last_update(:,iparticle ) = real_to_scaled(particle_set(iparticle)%r(:), &
                                                       cell)
        END DO
        update_neighbor_lists = .TRUE.
      END IF
    ELSE
      ALLOCATE(r_last_update(3,nparticle), STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN, moduleN, __LINE__, "r_last_update",&
                         3*nparticle*dp_size)
      END IF

      DO iparticle = 1,nparticle
        r_last_update(:,iparticle ) = real_to_scaled(particle_set(iparticle)%r(:), &
                                                     cell)
      END DO

      update_neighbor_lists = .TRUE.
      build_from_scratch = .TRUE.
    END IF

    IF (update_neighbor_lists) THEN
      CALL build_tbmd_neighbor_lists(atomic_kind_set, particle_set, &
                                     local_particles, cell, r_max=r_cutoff, &
                                     tbmd_neighbor_list=tbmd_neighbor_list, &
                                     build_from_scratch=build_from_scratch, &
                                     force_env_section=force_env_section, &
                                     globenv=globenv, error=error)

      IF ( counter > 0 ) THEN
        num_update = num_update + 1
        lup = counter + 1 - last_update
        last_update = counter + 1
        aup = aup + (lup - aup)/REAL(num_update,KIND=dp)
      ELSE
        num_update = 0
        lup = 0
        last_update = 1
        aup = 0.0_dp
      END IF

      CALL tbmd_neighbor_list_env_set(tbmd_neighbor_list_env, lup=lup, aup=aup, &
                                      r_last_update=r_last_update, &
                                      tbmd_neighbor_list=tbmd_neighbor_list, &
                                      num_update=num_update, &
                                      last_update=last_update, error=error)

      IF ((globenv%print_level >= 2).AND.(globenv%ionode)) THEN
        WRITE (UNIT=globenv%scr, &
               FMT="(/,T2,A,/,T52,A,/,A,T31,A,T49,2(1X,F15.2),/,T2,A,/)")&
          REPEAT("*",79),"INSTANTANEOUS        AVERAGES",&
          " LIST UPDATES[steps]","= ", lup, aup, REPEAT("*",79)
      END IF
    END IF

    counter = counter + 1
    CALL tbmd_neighbor_list_env_set(tbmd_neighbor_list_env, counter=counter)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_neighbor_list_control
! *****************************************************************************

!!****f* tbmd_neighbor_lists/build_tbmd_neighbor_lists [1.0] *
!!
!!   NAME
!!     build_tbmd_neighbor_lists
!!
!!   FUNCTION
!!     Builds all the required neighbor lists for TBMD
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - atomic_kind_set: The set of all involvolved atomic kinds
!!     - particle_set: The set of particles
!!     - local_particles: Particles which are local on this particluar node
!!     - cell: The simulation cell
!!     - r_max: Maximum interaction radius
!!     - tbmd_neighbor_list: The TBMD neighbor list
!!     - build_from_scratch: Logical which indicates if the neighbor list can
!!                           be updated or has to be constructed from scratch
!!     - force_env_section: Pointer to the force environment input section
!!     - globenv: The global environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE build_tbmd_neighbor_lists(atomic_kind_set, particle_set, &
                                       local_particles, cell, r_max, &
                                       tbmd_neighbor_list, globenv, &
                                       build_from_scratch, force_env_section, &
                                       error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                   :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                   :: particle_set
    TYPE(distribution_1d_type), POINTER         :: local_particles
    TYPE(cell_type), POINTER                    :: cell
    REAL(dp), DIMENSION(:, :), INTENT(IN)       :: r_max
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                     :: tbmd_neighbor_list
    TYPE(global_environment_type), &
      POINTER                                   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                                  :: error
    TYPE(section_vals_type), POINTER            :: force_env_section

    LOGICAL, INTENT(IN)                         :: build_from_scratch

    CHARACTER(LEN=*), PARAMETER :: routineN = "build_tbmd_neighbor_lists"

    LOGICAL                                     :: print_subcell_grid
    INTEGER                                     :: handle, &
                                                   output_unit, &
                                                   maxatom, &
                                                   maxatom_local, &
                                                   nkind, &
                                                   istat, &
                                                   ikind, &
                                                   natom_a, &
                                                   natom_local_a, &
                                                   iatom, &
                                                   atom_a, &
                                                   iatom_local, &
                                                   iw
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                              :: atom
    TYPE(atomic_kind_type), POINTER             :: atomic_kind
    REAL(dp), DIMENSION(3)                      :: r_pbc

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(logger, atomic_kind)
    logger => cp_error_get_logger(error)

    print_subcell_grid = .FALSE.
    output_unit = cp_print_key_unit_nr(logger, force_env_section, &
                                       "SUBSYS%PRINT%SUBCELL", &
                                       extension=".Log", error=error)
    IF (output_unit>0)  print_subcell_grid = .TRUE.

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxatom=maxatom)

!   *** Allocate work storage ***
    nkind = SIZE(atomic_kind_set)

    ALLOCATE (atom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "atom", nkind*int_size)

!   *** Initialize the local data structures ***
    maxatom_local = 0

    DO ikind = 1,nkind
      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom(ikind)%list)
      NULLIFY (atom(ikind)%list_local_a_index)
      NULLIFY (atom(ikind)%r_pbc)
      NULLIFY (atom(ikind)%s_pbc)

      CALL get_atomic_kind(atomic_kind=atomic_kind, &
                           atom_list=atom(ikind)%list)

      natom_a = SIZE(atom(ikind)%list)
      natom_local_a = local_particles%n_el(ikind)

      maxatom_local = MAX(maxatom_local,natom_local_a)

      IF (natom_local_a > 0) THEN
        ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a), STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                         "atom(ikind)%list_local_a_index", &
                                         natom_local_a*int_size)

!       *** Build index vector for mapping ***
        DO iatom_local = 1,natom_local_a
          atom_a = local_particles%list(ikind)%array(iatom_local)
          DO iatom = 1,natom_a
            IF (atom_a == atom(ikind)%list(iatom)) THEN
              atom(ikind)%list_local_a_index(iatom_local) = iatom
              EXIT
            END IF
          END DO
        END DO
      END IF

      ALLOCATE(atom(ikind)%r_pbc(3,natom_a), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                       "atom(ikind)%r_pbc", 3*natom_a*dp_size)

      ALLOCATE(atom(ikind)%s_pbc(3,natom_a), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                       "atom(ikind)%s_pbc", 3*natom_a*dp_size)

!     *** Calculate PBC coordinates ***
      DO iatom = 1,natom_a
        atom_a = atom(ikind)%list(iatom)
        r_pbc(:) = pbc(particle_set(atom_a)%r(:), cell)
        atom(ikind)%r_pbc(:,iatom) = r_pbc(:)
        atom(ikind)%s_pbc(:,iatom) = real_to_scaled(r_pbc(:), cell)
      END DO

    END DO

!   *** Build the TBMD neighbor lists ***
    CALL build_neighbor_lists(tbmd_neighbor_list, particle_set, atom, cell, &
                              maxatom_local, print_subcell_grid, &
                              output_unit, r_max, build_from_scratch)

    IF (BTEST(cp_print_key_should_output(logger%iter_info, &
         force_env_section, "TBMD%PRINT%NEIGHBOR_LISTS", &
         error=error), cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,force_env_section, &
                                 "TBMD%PRINT%NEIGHBOR_LISTS", &
                                 extension=".mmLog", error=error)
       CALL write_neighbor_lists(tbmd_neighbor_list, particle_set, cell, &
                                 globenv, iw, "NONBONDED NEIGHBOR LISTS", &
                                 error=error)
       CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                         "TBMD%PRINT%NEIGHBOR_LISTS", &
                                         error=error)
    END IF

!   *** Print particle distribution ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info, &
         force_env_section, "PRINT%DISTRIBUTION", error=error), &
         cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger, force_env_section, &
                                 "PRINT%DISTRIBUTION", extension=".mmLog", &
                                 error=error)
       CALL write_distribution(tbmd_neighbor_list, atomic_kind_set, iw, &
                               globenv, error=error)
       CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                         "PRINT%DISTRIBUTION", error=error)
    END IF

!   *** Release work storage ***
    DO ikind = 1,nkind
      NULLIFY(atom(ikind)%list)
      IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
        DEALLOCATE(atom(ikind)%list_local_a_index, STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                         "atom(ikind)%list_local_a_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%r_pbc)) THEN
        DEALLOCATE(atom(ikind)%r_pbc, STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                         "atom(ikind)%r_pbc")
      END IF
      IF (ASSOCIATED(atom(ikind)%s_pbc)) THEN
        DEALLOCATE (atom(ikind)%s_pbc, STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                         "atom(ikind)%s_pbc")
      END IF
    END DO

    DEALLOCATE(atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, "atom")

    CALL cp_print_key_finished_output(output_unit, logger, &
                                      force_env_section, &
                                      "SUBSYS%PRINT%SUBCELL", error=error)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE build_tbmd_neighbor_lists

! *****************************************************************************

END MODULE tbmd_neighbor_lists
