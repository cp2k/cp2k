!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE qs_fb_distribution_methods

  USE cell_types,                      ONLY: cell_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_get_info, cp_dbcsr_nblkcols_total, cp_dbcsr_p_type, &
       cp_dbcsr_type, dbcsr_distribution_mp, dbcsr_distribution_obj, &
       dbcsr_mp_get_process, dbcsr_mp_npcols, dbcsr_mp_nprows, dbcsr_mp_obj
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_fb_atomic_halo_types,         ONLY: &
       fb_atomic_halo_build_halo_atoms, fb_atomic_halo_cost, &
       fb_atomic_halo_create, fb_atomic_halo_init, fb_atomic_halo_nullify, &
       fb_atomic_halo_obj, fb_atomic_halo_release, fb_atomic_halo_set, &
       fb_build_pair_radii
  USE qs_fb_env_types,                 ONLY: fb_env_get,&
                                             fb_env_obj,&
                                             fb_env_set
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fb_distribution_build

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_distribution_methods'

! *****************************************************************************
!> \brief derived type containing cost data used for process distribution
!> \param id               : global atomic index
!> \param cost             : computational cost for the atomic matrix associated
!>                           to this atom
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_distribution_element
     INTEGER :: id
     REAL(KIND=dp) :: cost
  END TYPE fb_distribution_element

! *****************************************************************************
!> \brief derived type containing the list of atoms currently allocated to a
!>        processor
!> \param list             : list of atoms and their associated costs
!> \param cost             : total cost of the list
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_distribution_list
     TYPE(fb_distribution_element), DIMENSION(:), POINTER :: list => NULL()
     INTEGER :: nelements
     REAL(KIND=dp) :: cost
  END TYPE fb_distribution_list

! *****************************************************************************
!> \brief In filter matrix algorithm, each atomic matrix contributes to a
!>        column in the filter matrix, which is stored in DBCSR format.
!>        When distributing the atoms (and hence the atomic matrics) to the
!>        processors, we want the processors to have atoms that would be
!>        correspond to the block columns in the DBCSR format local to them.
!>        This derived type stores this information. For each atom, it
!>        corresponds to a DBCSR block column, and the list of processors
!>        in the 2D processor grid responsible for this column will be the
!>        preferred processors for this atom.
!> \param list             : list of preferred processors for an atom
!>                           note that here the processors are indexed from
!>                           1, i.e. = MPI_RANK+1
!> \param nprocs           : number of processors in the list
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  TYPE fb_preferred_procs_list
     INTEGER, DIMENSION(:), POINTER :: list => NULL()
     INTEGER :: nprocs
  END TYPE fb_preferred_procs_list

! Parameters related to automatic resizing of the hash_table:
! Resize by EXPAND_FACTOR if total no. slots / no. of filled slots < ENLARGE_RATIO
  INTEGER, PARAMETER, PRIVATE :: ENLARGE_RATIO = 1
  INTEGER, PARAMETER, PRIVATE :: REDUCE_RATIO = 3
  INTEGER, PARAMETER, PRIVATE :: EXPAND_FACTOR = 2
  INTEGER, PARAMETER, PRIVATE :: SHRINK_FACTOR = 2

  INTERFACE fb_distribution_remove
     MODULE PROCEDURE fb_distribution_remove_ind, &
                      fb_distribution_remove_el
  END INTERFACE fb_distribution_remove

  INTERFACE fb_distribution_move
     MODULE PROCEDURE fb_distribution_move_ind, &
                      fb_distribution_move_el
  END INTERFACE fb_distribution_move

CONTAINS

! *****************************************************************************
!> \brief Build local atoms associated to filter matrix algorithm for each
!>        MPI process, trying to balance the load for calculating the
!>        filter matrix
!> \param fb_env : the filter matrix environment
!> \param qs_env : quickstep environment
!> \param scf_section : SCF input section
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_build(fb_env, qs_env, scf_section)
    TYPE(fb_env_obj), INTENT(INOUT)          :: fb_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_build', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i_common_set, iatom, ii, ipe, lb, lowest_cost_ind, &
      my_pe, n_common_sets, natoms, nhalo_atoms, nkinds, nprocs, &
      owner_id_in_halo, pref_pe, stat, ub
    INTEGER, ALLOCATABLE, DIMENSION(:) :: common_set_ids, local_atoms_all, &
      local_atoms_sizes, local_atoms_starts, pe, pos_in_preferred_list
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms, local_atoms
    LOGICAL                                  :: acceptable_move, failure, &
                                                move_happened
    REAL(KIND=dp)                            :: average_cost
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cost_per_atom, cost_per_proc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pair_radii
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rcut
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: mat_ks
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fb_atomic_halo_obj)                 :: atomic_halo
    TYPE(fb_distribution_element)            :: element
    TYPE(fb_distribution_list), &
      ALLOCATABLE, DIMENSION(:)              :: dist
    TYPE(fb_preferred_procs_list), &
      ALLOCATABLE, DIMENSION(:)              :: preferred_procs_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(mat_ks, rcut, cell, para_env, particle_set, qs_kind_set, &
            halo_atoms, local_atoms)
    CALL fb_atomic_halo_nullify(atomic_halo)

    ! obtain relevant data from fb_env, qs_env
    CALL fb_env_get(fb_env=fb_env, &
                    rcut=rcut)
    CALL get_qs_env(qs_env=qs_env, &
                    natom=natoms, &
                    particle_set=particle_set, &
                    qs_kind_set=qs_kind_set, &
                    nkind=nkinds, &
                    cell=cell, &
                    para_env=para_env, &
                    matrix_ks=mat_ks)
    nprocs = para_env%num_pe
    my_pe = para_env%mepos + 1 ! counting from 1

    ! for each global atom, build atomic halo and get the associated cost
    ALLOCATE(pair_radii(nkinds,nkinds), STAT=stat)
    CPASSERT(stat==0)
    CALL fb_build_pair_radii(rcut, nkinds, pair_radii)
    CALL fb_atomic_halo_create(atomic_halo)
    ALLOCATE(cost_per_atom(natoms), STAT=stat)
    CPASSERT(stat==0)
    DO iatom = 1, natoms
       CALL fb_atomic_halo_init(atomic_halo)
       CALL fb_atomic_halo_build_halo_atoms(iatom, &
                                            particle_set, &
                                            cell, &
                                            pair_radii, &
                                            halo_atoms, &
                                            nhalo_atoms, &
                                            owner_id_in_halo)
       CALL fb_atomic_halo_set(atomic_halo=atomic_halo, &
                               owner_atom=iatom, &
                               natoms=nhalo_atoms, &
                               halo_atoms=halo_atoms)
       NULLIFY(halo_atoms)
       cost_per_atom(iatom) = fb_atomic_halo_cost(atomic_halo, particle_set, qs_kind_set)
    END DO
    DEALLOCATE(pair_radii, STAT=stat)
    CPASSERT(stat==0)
    CALL fb_atomic_halo_release(atomic_halo)

    ! build the preferred_procs_set according to DBCSR mat H
    ALLOCATE(preferred_procs_set(natoms), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(common_set_ids(natoms), STAT=stat)
    CPASSERT(stat==0)
    CALL fb_build_preferred_procs(mat_ks(1)%matrix, &
                                  natoms, &
                                  preferred_procs_set, &
                                  common_set_ids, &
                                  n_common_sets)

    ! for each atomic halo, construct distribution_element, and assign
    ! the element to a processors using preferred_procs_set in a
    ! round-robin manner
    ALLOCATE(dist(nprocs), STAT=stat)
    CPASSERT(stat==0)
    DO ipe = 1, nprocs
       CALL fb_distribution_init(dist=dist(ipe))
    END DO
    ALLOCATE(pos_in_preferred_list(n_common_sets), STAT=stat)
    CPASSERT(stat==0)
    pos_in_preferred_list(:) = 0
    DO iatom = 1, natoms
       element%id = iatom
       element%cost = cost_per_atom(iatom)
       i_common_set = common_set_ids(iatom)
       pos_in_preferred_list(i_common_set) = &
            MOD(pos_in_preferred_list(i_common_set), &
                preferred_procs_set(iatom)%nprocs) + 1
       ipe = preferred_procs_set(iatom)%list(pos_in_preferred_list(i_common_set))
       CALL fb_distribution_add(dist(ipe), element)
    END DO

    DEALLOCATE(pos_in_preferred_list, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(common_set_ids, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(cost_per_atom, STAT=stat)
    CPASSERT(stat==0)

    ! sort processors according to the overall cost of their assigned
    ! corresponding distribution
    ALLOCATE(cost_per_proc(nprocs), STAT=stat)
    CPASSERT(stat==0)
    DO ipe = 1, nprocs
       cost_per_proc(ipe) = dist(ipe)%cost
    END DO
    ALLOCATE(pe(nprocs), STAT=stat)
    CPASSERT(stat==0)
    CALL sort(cost_per_proc, nprocs, pe)
    ! now that cost_per_proc is sorted, ipe's no longer give mpi
    ! ranks, the correct one to use should be pe(ipe)

    ! work out the ideal average cost per proc if work load is evenly
    ! distributed
    average_cost = SUM(cost_per_proc) / REAL(nprocs, dp)

    DEALLOCATE(cost_per_proc, STAT=stat)
    CPASSERT(stat==0)

    ! loop over the processors, starting with the highest cost, move
    ! atoms one by one:
    !   1. FIRST to the next processor in the preferred list that has
    !      cost below average. IF no such proc is found, THEN
    !   2. to the next procesor in the overall list that has cost
    !      below average.
    ! repeat until the cost on this processor is less than or equal
    ! to the average cost
    lowest_cost_ind = 1
    DO ipe = nprocs, 1, -1
       redistribute: DO WHILE (dist(pe(ipe))%cost .GT. average_cost)
          iatom = dist(pe(ipe))%list(lowest_cost_ind)%id
          move_happened = .FALSE.
          ! first try to move to a preferred process
          preferred: DO ii = 1, preferred_procs_set(iatom)%nprocs
             pref_pe = preferred_procs_set(iatom)%list(ii)
             acceptable_move = &
                  fb_distribution_acceptable_move(dist(pe(ipe)), &
                                                  dist(pe(ipe))%list(lowest_cost_ind), &
                                                  dist(pref_pe), &
                                                  average_cost)
             IF ((pref_pe .NE. pe(ipe)) .AND. acceptable_move) THEN
                CALL fb_distribution_move(dist(pe(ipe)), &
                                          lowest_cost_ind, &
                                          dist(pref_pe))
                move_happened = .TRUE.
                EXIT preferred
             END IF
          END DO preferred
          ! if no preferred process is avaliable, move to a proc in
          ! the sorted list that has cost less than average.  remember
          ! that some of the proc may have already taken redistributed
          ! atoms, and thus may beome unavaliable (full)
          IF (.NOT. move_happened) THEN
             ! searching from the proc with the least initial cost
             next_in_line: DO ii = 1, nprocs
                acceptable_move = &
                  fb_distribution_acceptable_move(dist(pe(ipe)), &
                                                  dist(pe(ipe))%list(lowest_cost_ind), &
                                                  dist(pe(ii)), &
                                                  average_cost)
                IF ((pe(ii) .NE. pe(ipe)) .AND. acceptable_move) THEN
                   CALL fb_distribution_move(dist(pe(ipe)), &
                                             lowest_cost_ind, &
                                             dist(pe(ii)))
                   move_happened = .TRUE.
                   EXIT next_in_line
                END IF
             END DO next_in_line
          END IF
          ! if the atom cannot be moved, then this means it is too
          ! costly for all other processes to accept. When this
          ! happens we must stop the redistribution process for this
          ! processor---as all other of its atoms will be even more
          ! costly
          IF (.NOT. move_happened) THEN
             EXIT redistribute
          END IF
       END DO redistribute ! while
    END DO ! ipe

    DEALLOCATE(pe, STAT=stat)
    CPASSERT(stat==0)
    DO ii = 1, SIZE(preferred_procs_set)
       CALL fb_preferred_procs_list_release(preferred_procs_set(ii))
    END DO
    DEALLOCATE(preferred_procs_set, STAT=stat)
    CPASSERT(stat==0)

    ! generate local atoms from dist
    ALLOCATE(local_atoms_all(natoms), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(local_atoms_starts(nprocs), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(local_atoms_sizes(nprocs), STAT=stat)
    CPASSERT(stat==0)
    CALL fb_distribution_to_local_atoms(dist, &
                                        local_atoms_all, &
                                        local_atoms_starts, &
                                        local_atoms_sizes)
    ALLOCATE(local_atoms(local_atoms_sizes(my_pe)), STAT=stat)
    CPASSERT(stat==0)
    lb = local_atoms_starts(my_pe)
    ub = local_atoms_starts(my_pe) + local_atoms_sizes(my_pe) - 1
    local_atoms(1:local_atoms_sizes(my_pe)) = local_atoms_all(lb:ub)
    CALL fb_env_set(fb_env=fb_env, &
                    local_atoms=local_atoms, &
                    nlocal_atoms=local_atoms_sizes(my_pe))

    ! write out info
    CALL fb_distribution_write_info(dist, scf_section)

    DEALLOCATE(local_atoms_all, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(local_atoms_starts, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(local_atoms_sizes, STAT=stat)
    CPASSERT(stat==0)
    DO ipe = 1, SIZE(dist)
       CALL fb_distribution_release(dist(ipe))
    END DO
    DEALLOCATE(dist, STAT=stat)
    CPASSERT(stat==0)

    CALL timestop(handle)

  END SUBROUTINE fb_distribution_build

! *****************************************************************************
!> \brief Checks if moving an element from one distribution to another is
!>        allowed in mind of load balancing.
!> \param dist_from : the source distribution
!> \param element   : the element in source distribution considered for the
!>                    move
!> \param dist_to   : the destination distribution
!> \param threshold ...
!> \retval res : TRUE or FALSE
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  PURE FUNCTION fb_distribution_acceptable_move(dist_from, &
                                                element, &
                                                dist_to, &
                                                threshold) &
       RESULT(acceptable)
    TYPE(fb_distribution_list), INTENT(IN)   :: dist_from
    TYPE(fb_distribution_element), &
      INTENT(IN)                             :: element
    TYPE(fb_distribution_list), INTENT(IN)   :: dist_to
    REAL(KIND=dp), INTENT(IN)                :: threshold
    LOGICAL                                  :: acceptable

    acceptable = (dist_to%cost + element%cost .LT. dist_from%cost) .AND. &
                 (dist_to%cost .LT. threshold)
  END FUNCTION fb_distribution_acceptable_move

! *****************************************************************************
!> \brief Write out information on the load distribution on processors
!> \param dist_set    : set of distributions for the processors
!> \param scf_section : SCF input section
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_write_info(dist_set, scf_section)
    TYPE(fb_distribution_list), &
      DIMENSION(:), INTENT(IN)               :: dist_set
    TYPE(section_vals_type), POINTER         :: scf_section

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_write_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, max_natoms, min_natoms, &
                                                natoms, nprocs, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ave_cost, ave_natoms, &
                                                max_cost, min_cost, total_cost
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.

    nprocs = SIZE(dist_set)
    natoms = 0
    total_cost = 0.0_dp
    DO ii = 1, nprocs
       natoms = natoms + dist_set(ii)%nelements
       total_cost = total_cost + dist_set(ii)%cost
    END DO
    ave_natoms = REAL(natoms, dp) / REAL(nprocs, dp)
    ave_cost = total_cost / REAL(nprocs, dp)
    max_natoms = 0
    max_cost = 0._dp
    DO ii = 1, nprocs
       max_natoms = MAX(max_natoms, dist_set(ii)%nelements)
       max_cost = MAX(max_cost, dist_set(ii)%cost)
    END DO
    min_natoms = natoms
    min_cost = total_cost
    DO ii = 1, nprocs
       min_natoms = MIN(min_natoms, dist_set(ii)%nelements)
       min_cost = MIN(min_cost, dist_set(ii)%cost)
    END DO

    logger => cp_get_default_logger()
    unit_nr = cp_print_key_unit_nr(logger, scf_section, &
                                   "PRINT%FILTER_MATRIX", &
                                   extension="")

    IF (unit_nr > 0) THEN
       WRITE (UNIT=unit_nr, FMT="(/,A,I6,A)") &
            " FILTER_MAT_DIAG| Load distribution across ", nprocs, " processors:"
       WRITE (UNIT=unit_nr, &
              FMT="(A,T40,A,T55,A,T70,A,T85,A)") &
            " FILTER_MAT_DIAG| ", "Total", "Average", "Max", "Min"
       WRITE (UNIT=unit_nr, &
              FMT="(A,T40,I12,T55,F12.1,T70,I12,T85,I10)") &
            " FILTER_MAT_DIAG|   Atomic Matrices", &
            natoms, ave_natoms, max_natoms, min_natoms
       WRITE (UNIT=unit_nr, &
              FMT="(A,T40,D12.7,T55,D12.7,T70,D12.7,T85,D12.7)") &
            " FILTER_MAT_DIAG|   Cost*", &
            total_cost, ave_cost, max_cost, min_cost
       WRITE (UNIT=unit_nr, FMT="(A)") &
            " FILTER_MAT_DIAG| (* cost is calculated as sum of cube of atomic matrix sizes)"
    END IF
    CALL cp_print_key_finished_output(unit_nr, logger, scf_section, &
                                      "PRINT%FILTER_MATRIX")
  END SUBROUTINE fb_distribution_write_info

! *****************************************************************************
!> \brief Build the preferred list of processors for atoms
!> \param dbcsr_mat   : the reference DBCSR matrix, from which the local block
!>                      cols and the processor maps are obtained
!> \param natoms      : total number of atoms globally
!> \param preferred_procs_set : set of preferred procs list for each atom
!> \param common_set_ids : atoms (block cols) local to the same processor grid
!>                         col will have the same preferred list. This list
!>                         maps each atom to their corresponding group
!> \param n_common_sets  : number of unique preferred lists (groups)
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_build_preferred_procs(dbcsr_mat, &
                                      natoms, &
                                      preferred_procs_set, &
                                      common_set_ids, &
                                      n_common_sets)
    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_mat
    INTEGER, INTENT(IN)                      :: natoms
    TYPE(fb_preferred_procs_list), &
      DIMENSION(:), INTENT(INOUT)            :: preferred_procs_set
    INTEGER, DIMENSION(:), INTENT(OUT)       :: common_set_ids
    INTEGER, INTENT(OUT)                     :: n_common_sets

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_build_preferred_procs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icol, nprows, pcol, prow, stat
    INTEGER, DIMENSION(:), POINTER           :: col_dist
    LOGICAL                                  :: check_ok, failure
    TYPE(dbcsr_distribution_obj)             :: dbcsr_dist
    TYPE(dbcsr_mp_obj)                       :: dbcsr_mp

    failure = .FALSE.
    check_ok = natoms .LE. cp_dbcsr_nblkcols_total(dbcsr_mat)
    CPASSERT(check_ok)
    check_ok = SIZE(preferred_procs_set) .GE. natoms
    CPASSERT(check_ok)
    check_ok = SIZE(common_set_ids) .GE. natoms
    CPASSERT(check_ok)

    CALL cp_dbcsr_get_info(dbcsr_mat, distribution=dbcsr_dist, proc_col_dist=col_dist)
    dbcsr_mp = dbcsr_distribution_mp(dbcsr_dist)
    nprows = dbcsr_mp_nprows(dbcsr_mp)

    DO icol = 1, natoms
       IF (ASSOCIATED(preferred_procs_set(icol)%list)) THEN
          DEALLOCATE(preferred_procs_set(icol)%list, STAT=stat)
          CPASSERT(stat==0)
       END IF
       ALLOCATE(preferred_procs_set(icol)%list(nprows), STAT=stat)
       CPASSERT(stat==0)
       pcol = col_dist(icol)
       ! dbcsr prow and pcol counts from 0
       DO prow = 0, nprows-1
          ! here, we count processes from 1, so +1 from mpirank
          preferred_procs_set(icol)%list(prow+1) = &
               dbcsr_mp_get_process(dbcsr_mp, prow, pcol) + 1
       END DO
       preferred_procs_set(icol)%nprocs = nprows
    END DO

    common_set_ids(:) = 0
    common_set_ids(1:natoms) = col_dist(1:natoms) + 1
    n_common_sets = dbcsr_mp_npcols(dbcsr_mp)

  END SUBROUTINE fb_build_preferred_procs

! *****************************************************************************
!> \brief Release a preferred_procs_list
!> \param preferred_procs_list  : the preferred procs list in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_preferred_procs_list_release(preferred_procs_list)
    TYPE(fb_preferred_procs_list), &
      INTENT(INOUT)                          :: preferred_procs_list

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'fb_preferred_procs_list_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(preferred_procs_list%list)) THEN
       DEALLOCATE(preferred_procs_list%list, STAT=stat)
       CPASSERT(stat==0)
    END IF
  END SUBROUTINE fb_preferred_procs_list_release

! *****************************************************************************
!> \brief Convert distribution data to 1D array containing information of
!>        which atoms are distributed to which processor
!> \param dist_set    : set of distributions for the processors
!> \param local_atoms : continuous array of atoms arranged in order
!>                      corresponding their allocated processors
!> \param local_atoms_starts : starting position in local_atoms array for
!>                             each processor
!> \param local_atoms_sizes  : number of atoms local to each processor
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_to_local_atoms(dist_set, &
                                            local_atoms, &
                                            local_atoms_starts, &
                                            local_atoms_sizes)
    TYPE(fb_distribution_list), &
      DIMENSION(:), INTENT(IN)               :: dist_set
    INTEGER, DIMENSION(:), INTENT(OUT)       :: local_atoms, &
                                                local_atoms_starts, &
                                                local_atoms_sizes

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'fb_distribution_to_local_atoms', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ipe, n_procs, pos
    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.

    n_procs = SIZE(dist_set)

    check_ok = SIZE(local_atoms_starts) .GE. n_procs
    CPASSERT(check_ok)
    check_ok = SIZE(local_atoms_sizes) .GE. n_procs
    CPASSERT(check_ok)

    local_atoms(:) = 0
    local_atoms_starts(:) = 0
    local_atoms_sizes(:) = 0

    pos = 1
    DO ipe = 1, n_procs
       local_atoms_starts(ipe) = pos
       DO iatom = 1, dist_set(ipe)%nelements
          local_atoms(pos) = dist_set(ipe)%list(iatom)%id
          pos = pos + 1
          local_atoms_sizes(ipe) = local_atoms_sizes(ipe) + 1
       END DO
    END DO
  END SUBROUTINE fb_distribution_to_local_atoms

! *****************************************************************************
!> \brief Initialise a distribution
!> \param dist        : the distribution in question
!> \param nmax        : [OPTIONAL] size of the list array to be allocated
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_init(dist, nmax)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist
    INTEGER, INTENT(IN), OPTIONAL            :: nmax

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, my_nmax, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    my_nmax = 0
    IF (PRESENT(nmax)) my_nmax = nmax
    IF (ASSOCIATED(dist%list)) THEN
       DEALLOCATE(dist%list, STAT=stat)
       CPASSERT(stat==0)
    END IF
    NULLIFY(dist%list)
    IF (my_nmax .GT. 0) THEN
       ALLOCATE(dist%list(my_nmax), STAT=stat)
       CPASSERT(stat==0)
       DO ii = 1, SIZE(dist%list)
          dist%list(ii)%id = 0
          dist%list(ii)%cost = 0.0_dp
       END DO
    END IF
    dist%nelements = 0
    dist%cost = 0.0_dp
  END SUBROUTINE fb_distribution_init

! *****************************************************************************
!> \brief Resize the list array in a distribution
!> \param dist        : The distribution in question
!> \param nmax        : new size of the list array
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_resize(dist, nmax)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist
    INTEGER, INTENT(IN)                      :: nmax

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_resize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, my_nmax, stat
    LOGICAL                                  :: failure
    TYPE(fb_distribution_element), &
      DIMENSION(:), POINTER                  :: new_list

    failure = .FALSE.
    IF (.NOT. ASSOCIATED(dist%list)) THEN
       my_nmax = MAX(nmax, 1)
       ALLOCATE(dist%list(my_nmax), STAT=stat)
       CPASSERT(stat==0)
    ELSE
       my_nmax = MAX(nmax, dist%nelements)
       ALLOCATE(new_list(my_nmax), STAT=stat)
       CPASSERT(stat==0)
       DO ii = 1, SIZE(new_list)
          new_list(ii)%id = 0
          new_list(ii)%cost = 0.0_dp
       END DO
       DO ii = 1, dist%nelements
          new_list(ii) = dist%list(ii)
       END DO
       DEALLOCATE(dist%list, STAT=stat)
       CPASSERT(stat==0)
       dist%list => new_list
    END IF
  END SUBROUTINE fb_distribution_resize

! *****************************************************************************
!> \brief Add an atom (element) to a distribution
!> \param dist        : the distribution in question
!> \param element     : the element to be added
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_add(dist, element)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist
    TYPE(fb_distribution_element), &
      INTENT(IN)                             :: element

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, new_nelements, pos
    LOGICAL                                  :: failure

    failure = .FALSE.
    new_nelements = dist%nelements + 1

    ! resize list if necessary
    IF (.NOT. ASSOCIATED(dist%list)) THEN
       CALL fb_distribution_resize(dist, new_nelements)
    ELSE IF (new_nelements*ENLARGE_RATIO .GT. SIZE(dist%list)) THEN
       CALL fb_distribution_resize(dist, SIZE(dist%list)*EXPAND_FACTOR)
    END IF
    ! assuming the list of elements is always sorted with respect to cost
    ! slot the new element into the appropriate spot
    IF (new_nelements == 1) THEN
       dist%list(1) = element
    ELSE
       pos = fb_distribution_find_slot(dist, element)
       DO ii = dist%nelements, pos, -1
          dist%list(ii + 1) = dist%list(ii)
       END DO
       dist%list(pos) = element
    END IF
    dist%nelements = new_nelements
    dist%cost = dist%cost + element%cost
  END SUBROUTINE fb_distribution_add

! *****************************************************************************
!> \brief Find the correct slot in the list array to add a new element, so that
!>        the list will always be ordered with respect to cost
!> \param dist        : the distribution in question
!> \param element     : element to be added
!> \retval pos        : the correct position to add the new element
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  PURE FUNCTION fb_distribution_find_slot(dist, element) RESULT(pos)
    TYPE(fb_distribution_list), INTENT(IN)   :: dist
    TYPE(fb_distribution_element), &
      INTENT(IN)                             :: element
    INTEGER                                  :: pos

    INTEGER                                  :: lower, middle, N, upper

    N = dist%nelements
    IF (element%cost .LT. dist%list(1)%cost) THEN
       pos = 1
       RETURN
    END IF
    IF (element%cost .GE. dist%list(N)%cost) THEN
       pos = N + 1
       RETURN
    END IF
    lower = 1
    upper = N
    DO WHILE ((upper - lower) .GT. 1)
       middle = (lower + upper) / 2
       IF (element%cost .LT. dist%list(middle)%cost) THEN
          upper = middle
       ELSE
          lower = middle
       END IF
    END DO
    pos = upper
  END FUNCTION fb_distribution_find_slot

! *****************************************************************************
!> \brief Remove the pos-th element from a distribution
!> \param dist        : the distribution in question
!> \param pos         : index of the element in the list array
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_remove_ind(dist, pos)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist
    INTEGER, INTENT(IN)                      :: pos

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_remove_ind', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii
    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = pos .GT. 0
    CPASSERT(check_ok)
    IF (pos .LE. dist%nelements) THEN
       dist%cost = dist%cost - dist%list(pos)%cost
       DO ii = pos, dist%nelements-1
          dist%list(ii) = dist%list(ii+1)
       END DO
       dist%list(dist%nelements)%id = 0
       dist%list(dist%nelements)%cost = 0.0_dp
       dist%nelements = dist%nelements - 1
       ! auto resize if required
       IF (dist%nelements*REDUCE_RATIO .LT. SIZE(dist%list)) THEN
          CALL fb_distribution_resize(dist, dist%nelements/SHRINK_FACTOR)
       END IF
    END IF
  END SUBROUTINE fb_distribution_remove_ind

! *****************************************************************************
!> \brief Remove a given element from a distribution
!> \param dist        : the distribution in question
!> \param element     : the element in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_remove_el(dist, element)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist
    TYPE(fb_distribution_element), &
      INTENT(IN)                             :: element

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_remove_el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, pos

    pos = dist%nelements + 1
    DO ii = 1, dist%nelements
       IF (element%id == dist%list(ii)%id) THEN
          pos = ii
          EXIT
       END IF
    END DO
    CALL fb_distribution_remove_ind(dist, pos)
  END SUBROUTINE fb_distribution_remove_el

! *****************************************************************************
!> \brief Move the pos-th element from a distribution to another
!> \param dist_from   : the source distribution
!> \param pos         : index of the element in the source distribution
!> \param dist_to     : the destination distribution
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_move_ind(dist_from, pos, dist_to)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist_from
    INTEGER, INTENT(IN)                      :: pos
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist_to

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_move_ind', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check_ok, failure
    TYPE(fb_distribution_element)            :: element

    failure = .FALSE.
    check_ok = ASSOCIATED(dist_from%list)
    CPASSERT(check_ok)
    check_ok = pos .LE. dist_from%nelements
    CPASSERT(check_ok)
    element = dist_from%list(pos)
    CALL fb_distribution_add(dist_to, element)
    CALL fb_distribution_remove(dist_from, pos)
  END SUBROUTINE fb_distribution_move_ind

! *****************************************************************************
!> \brief Move a given element from a distribution to another
!> \param dist_from   : the source distribution
!> \param element     : the element in question
!> \param dist_to     : the destination distribution
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_move_el(dist_from, element, dist_to)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist_from
    TYPE(fb_distribution_element), &
      INTENT(IN)                             :: element
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist_to

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_move_el', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.
    check_ok = ASSOCIATED(dist_from%list)
    CPASSERT(check_ok)
    CALL fb_distribution_add(dist_to, element)
    CALL fb_distribution_remove(dist_from, element)
  END SUBROUTINE fb_distribution_move_el

! *****************************************************************************
!> \brief Release a distribution
!> \param dist  : the distribution in question
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_distribution_release(dist)
    TYPE(fb_distribution_list), &
      INTENT(INOUT)                          :: dist

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_distribution_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(dist%list)) THEN
       DEALLOCATE(dist%list, STAT=stat)
       CPASSERT(stat==0)
    END IF
  END SUBROUTINE fb_distribution_release

END MODULE qs_fb_distribution_methods
