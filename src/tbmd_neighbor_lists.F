!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             plane_distance,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE fist_neighbor_list_types,        ONLY: add_neighbor_list,&
                                             add_neighbor_node,&
                                             allocate_neighbor_list_set,&
                                             deallocate_neighbor_list_set,&
                                             init_neighbor_list_set,&
                                             neighbor_list_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type
  USE fist_neighbor_lists,             ONLY: write_distribution,&
                                             write_neighbor_lists
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_max
  USE particle_types,                  ONLY: particle_type
  USE qs_neighbor_lists,               ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             local_atoms_type,&
                                             subcell_type
  USE tbmd_neighbor_list_env_types,    ONLY: tbmd_neighbor_list_env_get,&
                                             tbmd_neighbor_list_env_set,&
                                             tbmd_neighbor_list_env_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_neighbor_lists'

  ! Public subroutines

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

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(tbmd_neighbor_list_env_type), &
      POINTER                                :: tbmd_neighbor_list_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tbmd_neighbor_list_control'

    INTEGER :: counter, handle, iparticle, iparticle_kind, iparticle_local, &
      istat, last_update, nparticle, nparticle_kind, nparticle_local, &
      num_update
    LOGICAL                                  :: build_from_scratch, failure, &
                                                update_neighbor_lists
    REAL(KIND=dp)                            :: aup, dr2, dr2_max, lup, &
                                                verlet_skin
    REAL(KIND=dp), DIMENSION(3)              :: dr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: r_cutoff, r_last_update, r_max
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: tbmd_neighbor_list

! -------------------------------------------------------------------------

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

    ALLOCATE(r_max(nparticle_kind,nparticle_kind), STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineN, error, failure)
    r_max(:,:) = r_cutoff(:,:) + verlet_skin

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

      CALL mp_max(dr2_max, globenv%group)

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
                                     local_particles, cell, r_max=r_max, &
                                     tbmd_neighbor_list=tbmd_neighbor_list, &
                                     build_from_scratch=build_from_scratch, &
                                     force_env_section=force_env_section, &
                                     globenv=globenv, error=error)

      IF ( counter > 0 ) THEN
        num_update = num_update + 1
        lup = counter + 1 - last_update
        last_update = counter + 1
        aup = aup + (lup - aup)/REAL(num_update, KIND=dp)
      ELSE
        num_update = 0
        lup = 0.0_dp
        last_update = 1
        aup = 0.0_dp
      END IF

      CALL tbmd_neighbor_list_env_set(tbmd_neighbor_list_env, lup=lup, aup=aup, &
                                      r_last_update=r_last_update, &
                                      tbmd_neighbor_list=tbmd_neighbor_list, &
                                      num_update=num_update, &
                                      last_update=last_update, error=error)

      IF ((globenv%print_level >= 2).AND.(globenv%ionode)) THEN
        WRITE(UNIT=globenv%scr, &
              FMT="(/,T2,A,/,T52,A,/,A,T31,A,T49,2(1X,F15.2),/,T2,A,/)")&
          REPEAT("*",79),"INSTANTANEOUS        AVERAGES",&
          " LIST UPDATES[steps]","= ", lup, aup, REPEAT("*",79)
      END IF
    END IF

    counter = counter + 1
    CALL tbmd_neighbor_list_env_set(tbmd_neighbor_list_env, counter=counter)

    DEALLOCATE(r_max, STAT=istat)
    CPPostcondition(istat==0, cp_failure_level, routineN, error, failure)

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
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: tbmd_neighbor_list
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN)                      :: build_from_scratch
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_tbmd_neighbor_lists'

    INTEGER :: atom_a, handle, iatom, iatom_local, ikind, istat, iw, maxatom, &
      maxatom_local, natom_a, natom_local_a, nkind, output_unit
    LOGICAL                                  :: print_subcell_grid
    REAL(dp), DIMENSION(3)                   :: r_pbc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom

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

    ALLOCATE (atom(nkind), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                     "atom", nkind*int_size)

!   *** Initialize the local data structures ***
    maxatom_local = 0

    DO ikind = 1,nkind
      atomic_kind => atomic_kind_set(ikind)

      NULLIFY(atom(ikind)%list)
      NULLIFY(atom(ikind)%list_local_a_index)
      NULLIFY(atom(ikind)%r_pbc)
      NULLIFY(atom(ikind)%s_pbc)

      CALL get_atomic_kind(atomic_kind=atomic_kind, &
                           atom_list=atom(ikind)%list)

      natom_a = SIZE(atom(ikind)%list)
      natom_local_a = local_particles%n_el(ikind)
  !    natom_local_a = distribution_2d%n_local_rows(ikind)
  !    natom_local_b = distribution_2d%n_local_cols(ikind)

     maxatom_local = MAX(maxatom_local, natom_local_a)
  !   maxatom_local = MAX(maxatom_local, natom_local_a, natom_local_b)

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


!     *** Block columns ***
   !   IF (natom_local_b > 0) THEN

   !     ALLOCATE(atom(ikind)%list_local_b_index(natom_local_b), STAT=istat)
   !     IF (istat /=0 ) CALL stop_memory(routineN, moduleN, __LINE__, &
   !                                      "atom(ikind)%list_local_b_index", &
   !                                      natom_local_b*int_size)

!       *** Build index vector for mapping ***
   !     DO iatom_local = 1,natom_local_b
   !       atom_b = local_particles%list(ikind)%array(iatom_local)
   !       DO iatom = 1,natom_a
   !         IF (atom_b == atom(ikind)%list(iatom)) THEN
   !           atom(ikind)%list_local_b_index(iatom_local) = iatom
   !           EXIT
   !         END IF
   !       END DO
   !     END DO

   !   END IF


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
       CALL write_neighbor_lists(tbmd_neighbor_list,particle_set,cell,&
                                 globenv,iw,"TIGHT-BINDING NEIGHBOR LISTS",&
                                 error)
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
      
   !   IF (ASSOCIATED(atom(ikind)%list_local_b_index)) THEN
   !     DEALLOCATE(atom(ikind)%list_local_b_index, STAT=istat)
   !     IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
   !                                      "atom(ikind)%list_local_b_index")
   !   END IF

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

!!****f* tbmd_neighbor_lists/build_neighbor_lists [1.0] *
!!
!!   NAME
!!     build_neighbor_lists
!!
!!   FUNCTION
!!     Builds all the required neighbor lists for TBMD
!!
!!   NOTES
!!     Pretty much Matthias' build_neighbor_lists@fist_neighbor_lists.F
!!
!!   INPUTS
!!     - ab: The pointer to the neighbor-list which should be built.
!!     - particle_set: The particle set of matter
!!     - atom: The local atom type
!!     - cell: The simulation cell
!!     - maxatom_local: The maximum of atoms local to each other
!!     - print_subcell_grid: Indicates if the subcell grid should be printed or not
!!     - output_unit: The output unit
!!     - r_max: The maximum interaction radius
!!     - build_from_scratch: Logical which indicates if an update can be performed,
!!                           or if the neighbor-list has to be built-up from scratch.
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     02.2006 initial create [tdk]
!!
!!*** **********************************************************************

  SUBROUTINE build_neighbor_lists(ab, particle_set, atom, cell, maxatom_local, &
                                  print_subcell_grid, output_unit, r_max,&
                                  build_from_scratch)

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: r_max
    LOGICAL, INTENT(IN)                      :: build_from_scratch

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_neighbor_lists'

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_b, handle, i, iab, iatom, iatom_local, &
      iatom_subcell, icell, ikind, istat, j, jatom_local, jcell, jkind, k, &
      kcell, natom_local_a, natom_local_b, nkind
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b, ncell, &
                                                nsubcell, periodic
    LOGICAL                                  :: cell000, include_ab
    REAL(dp)                                 :: rab2, rab2_max, rab_max, &
                                                subcells, unit_of_length
    REAL(dp), DIMENSION(3)                   :: abc, r, rab, rb, s, sab, &
                                                sab_max, sb, sb_max, sb_min, &
                                                sb_pbc
    TYPE(neighbor_list_p_type), &
      DIMENSION(maxatom_local)               :: kind_a
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    nkind = SIZE(atom)

    IF (build_from_scratch) THEN

!     *** Deallocate the old neighbor list structure ***
      IF (ASSOCIATED(ab)) THEN
        DO iab = 1,SIZE(ab)
          CALL deallocate_neighbor_list_set(ab(iab)%neighbor_list_set)
        END DO
        DEALLOCATE(ab, STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, "ab")
      END IF

!     *** Allocate and initialize the new neighbor list structure ***
      ALLOCATE(ab(nkind*nkind), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                       "ab", nkind*nkind*int_size)

      DO iab=1,SIZE(ab)
        NULLIFY (ab(iab)%neighbor_list_set)
      END DO

    ELSE

      IF (.NOT.ASSOCIATED(ab)) THEN
        CALL stop_program(routineN, moduleN, __LINE__, &
                          "ab is not associated and cannot be updated")
      END IF

    END IF

!   *** Load informations about the simulation cell ***
    CALL get_cell(cell=cell, abc=abc, periodic=periodic, &
                  subcells=subcells, &
                  unit_of_length=unit_of_length, &
                  unit_of_length_name=unit_of_length_name)

!   *** Print headline ***
    IF (print_subcell_grid) THEN
      WRITE(UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE TIGHT-BINDING NEIGHBOR LISTS", &
        "Atomic kind pair","Grid size", &
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE

      natom_local_a = SIZE(atom(ikind)%list_local_a_index)

      DO jkind = 1,nkind

        iab = ikind + nkind*(jkind - 1)

        natom_local_b = SIZE(atom(jkind)%list)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = r_max(ikind, jkind)
        rab2_max = rab_max*rab_max

        sab_max(1) = rab_max/plane_distance(1,0,0,cell)
        sab_max(2) = rab_max/plane_distance(0,1,0,cell)
        sab_max(3) = rab_max/plane_distance(0,0,1,cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        nsubcell(:) = MIN(MAX(1,NINT(0.5_dp*subcells/sab_max(:))),20)

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_subcell_grid) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),KIND=dp)
          r(:) = s(:)*abc(:)
          WRITE(UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind, jkind, nsubcell(1:3), r(1:3)/unit_of_length
        END IF

        IF (build_from_scratch) THEN
          CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          ab(iab)%neighbor_list_set, &
                                          r_max=rab_max)
          neighbor_list_set => ab(iab)%neighbor_list_set
        ELSE
          neighbor_list_set => ab(iab)%neighbor_list_set
          IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE
          CALL init_neighbor_list_set(neighbor_list_set=neighbor_list_set, &
                                      r_max=rab_max)
        END IF

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(nsubcell) == 3) THEN

!         *** Case 1: No subcell grid is needed ***

          DO iatom_local = 1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set, &
                                   atom=atom_a, cell=cell_a, &
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local = 1,natom_local_b

            atom_b = atom(jkind)%list(jatom_local)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom_local)

            loop1_kcell: DO kcell = -ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              sb_max(3) = sb(3) + sab_max(3)
              IF (periodic(3) /= 0) THEN
                IF (sb_min(3) >= 0.5_dp) EXIT loop1_kcell
                IF (sb_max(3) < -0.5_dp) CYCLE loop1_kcell
              END IF
              cell_b(3) = kcell

              loop1_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                sb_max(2) = sb(2) + sab_max(2)
                IF (periodic(2) /= 0) THEN
                  IF (sb_min(2) >= 0.5_dp) EXIT loop1_jcell
                  IF (sb_max(2) < -0.5_dp) CYCLE loop1_jcell
                END IF
                cell_b(2) = jcell

                loop1_icell: DO icell = -ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (periodic(1) /= 0) THEN
                    IF (sb_min(1) >= 0.5_dp) EXIT loop1_icell
                    IF (sb_max(1) < -0.5_dp) CYCLE loop1_icell
                  END IF
                  cell_b(1) = icell

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  rb(:) = scaled_to_real(sb(:),cell)

                  DO iatom_local = 1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)
                    IF (cell000) THEN
                      IF (atom_a > atom_b) THEN
                        include_ab = (MODULO(atom_a + atom_b,2) == 0)
                      ELSE
                        include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                      END IF
                    ELSE
                      IF (atom_a == atom_b) THEN
                        include_ab = .TRUE.
                      ELSE IF (atom_a > atom_b) THEN
                        include_ab = (MODULO(atom_a + atom_b,2) == 0)
                      ELSE
                        include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                      END IF
                    END IF
                    IF (include_ab) THEN
                      rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                      IF (rab2 < rab2_max) THEN
                        sab(:) = real_to_scaled(rab(:),cell)
                        CALL add_neighbor_node(&
                          neighbor_list=kind_a(iatom_local)%neighbor_list, &
                          neighbor=atom_b, r=sab, &   ! r=rab
                          exclusion_list = particle_set(atom_a)%lpi%&
                                           exclusion(jkind)%list_vdw)
                      END IF
                    END IF
                  END DO

                END DO loop1_icell
              END DO loop1_jcell
            END DO loop1_kcell

          END DO

        ELSE

!         *** Case 2: A subcell grid is needed ***

          CALL allocate_subcell(subcell, nsubcell, maxatom_local)

          DO iatom_local = 1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = INT(s(1)*REAL(nsubcell(1),KIND=dp)) + 1
            j = INT(s(2)*REAL(nsubcell(2),KIND=dp)) + 1
            k = INT(s(3)*REAL(nsubcell(3),KIND=dp)) + 1
            i = MIN(MAX(i,1),nsubcell(1))
            j = MIN(MAX(j,1),nsubcell(2))
            k = MIN(MAX(k,1),nsubcell(3))
            subcell(i,j,k)%natom  = subcell(i,j,k)%natom + 1
            subcell(i,j,k)%atom_list(subcell(i,j,k)%natom) = iatom_local
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set, &
                                   atom=atom_a, cell=cell_a, &
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local = 1,natom_local_b

            atom_b = atom(jkind)%list(jatom_local)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom_local)

            loop2_kcell: DO kcell = -ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              sb_max(3) = sb(3) + sab_max(3)
              IF (periodic(3) /= 0) THEN
                IF (sb_min(3) >= 0.5_dp) EXIT loop2_kcell
                IF (sb_max(3) < -0.5_dp) CYCLE loop2_kcell
              END IF
              cell_b(3) = kcell

              loop2_jcell: DO jcell = -ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                sb_max(2) = sb(2) + sab_max(2)
                IF (periodic(2) /= 0) THEN
                  IF (sb_min(2) >= 0.5_dp) EXIT loop2_jcell
                  IF (sb_max(2) < -0.5_dp) CYCLE loop2_jcell
                END IF
                cell_b(2) = jcell

                loop2_icell: DO icell = -ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (periodic(1) /= 0) THEN
                    IF (sb_min(1) >= 0.5_dp) EXIT loop2_icell
                    IF (sb_max(1) < -0.5_dp) CYCLE loop2_icell
                  END IF
                  cell_b(1) = icell

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  rb(:) = scaled_to_real(sb(:),cell)

                  loop_k: DO k = 1,nsubcell(3)
                    loop_j: DO j = 1,nsubcell(2)
                      loop_i: DO i = 1,nsubcell(1)

                        IF (periodic(3) /= 0) THEN
                          IF (sb_max(3) < subcell(i,j,k)%s_min(3)) EXIT loop_k
                          IF (sb_min(3) >= subcell(i,j,k)%s_max(3)) CYCLE loop_k
                        END IF

                        IF (periodic(2) /= 0) THEN
                          IF (sb_max(2) < subcell(i,j,k)%s_min(2)) EXIT loop_j
                          IF (sb_min(2) >= subcell(i,j,k)%s_max(2)) CYCLE loop_j
                        END IF

                        IF (periodic(1) /= 0) THEN
                          IF (sb_max(1) < subcell(i,j,k)%s_min(1)) EXIT loop_i
                          IF (sb_min(1) >= subcell(i,j,k)%s_max(1)) CYCLE loop_i
                        END IF

                        IF (subcell(i,j,k)%natom == 0) CYCLE
 
                        DO iatom_subcell = 1,subcell(i,j,k)%natom
                          iatom_local = subcell(i,j,k)%atom_list(iatom_subcell)
                          iatom = atom(ikind)%list_local_a_index(iatom_local)
                          atom_a = atom(ikind)%list(iatom)
                          IF (cell000) THEN
                            IF (atom_a > atom_b) THEN
                              include_ab = (MODULO(atom_a + atom_b,2) == 0)
                            ELSE
                              include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                            END IF
                          ELSE
                            IF (atom_a == atom_b) THEN
                              include_ab = .TRUE.
                            ELSE IF (atom_a > atom_b) THEN
                              include_ab = (MODULO(atom_a + atom_b,2) == 0)
                            ELSE
                              include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                            END IF
                          END IF
                          IF (include_ab) THEN
                            rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                            IF (rab2 < rab2_max) THEN
                              sab(:) = real_to_scaled(rab(:),cell)
                              CALL add_neighbor_node(&
                                neighbor_list=kind_a(iatom_local)%neighbor_list, &
                                neighbor=atom_b, r=sab, &   ! r=rab
                                exclusion_list=particle_set(atom_a)%lpi%&
                                               exclusion(jkind)%list_vdw)
                            END IF
                          END IF
                        END DO

                      END DO loop_i
                    END DO loop_j
                  END DO loop_k

                END DO loop2_icell
              END DO loop2_jcell
            END DO loop2_kcell

          END DO

          CALL deallocate_subcell(subcell)

        END IF

      END DO
    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_neighbor_lists

  ! ***************************************************************************

END MODULE tbmd_neighbor_lists
