! *****************************************************************************
!> \brief Routines used for force-mixing QM/MM calculations
!> \par History
!>      2.2012 created [noam]
!> \author Noam Bernstein
! *****************************************************************************
MODULE qmmm_force_mixing
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_deallocate,&
                                             fist_neighbor_type
  USE fist_neighbor_lists,             ONLY: build_fist_neighbor_lists
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_section_types,             ONLY: &
       section_vals_add_values, section_vals_duplicate, section_vals_get, &
       section_vals_get_subs_vals, section_vals_get_subs_vals3, &
       section_vals_release, section_vals_remove_values, &
       section_vals_set_subs_vals, section_vals_type, section_vals_val_get, &
       section_vals_val_set, section_vals_write
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE qmmm_types,                      ONLY: force_mixing_label_QM_core,&
                                             force_mixing_label_QM_core_list,&
                                             force_mixing_label_QM_dynamics,&
                                             force_mixing_label_QM_dynamics_list,&
                                             force_mixing_label_buffer,&
                                             force_mixing_label_buffer_list,&
                                             force_mixing_label_none,&
                                             force_mixing_label_termination
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_force_mixing'

  PUBLIC :: setup_force_mixing_qmmm_sections, update_force_mixing_labels


!***
CONTAINS

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param subsys ...
!> \param qmmm_section ...
!> \param labels_changed ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2012 created [noam]
!> \author Noam Bernstein
! *****************************************************************************
   SUBROUTINE update_force_mixing_labels(force_env, subsys, qmmm_section, labels_changed, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: qmmm_section
    LOGICAL, OPTIONAL                        :: labels_changed
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_force_mixing_labels', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      POINTER                                :: adaptive_exclude_molecules(:)
    INTEGER :: i_rep_section, i_rep_val, im, ip, max_n_qm, n_new, &
      n_rep_exclude, n_rep_section, n_rep_val, natoms, &
      QM_extended_seed_min_label_val
    INTEGER, ALLOCATABLE                     :: new_full_labels(:), &
                                                orig_full_labels(:)
    INTEGER, POINTER :: breakable_bonds(:,:), broken_bonds(:), &
      cur_indices(:), cur_labels(:), mm_index_entry(:), molecule_of_atom(:), &
      new_indices(:), new_labels(:)
    LOGICAL                                  :: explicit, &
                                                QM_extended_seed_is_core_list
    REAL(dp), ALLOCATABLE                    :: nearest_dist(:)
    REAL(dp), POINTER                        :: r_buf(:), r_core(:), r_qm(:)
    TYPE(cell_type), POINTER                 :: cell
    TYPE(fist_neighbor_type), POINTER        :: nlist
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: force_mixing_section, &
                                                non_adaptive_section, &
                                                qm_kind_section, &
                                                restart_section

IF (debug_this_module) PRINT *, "BOB starting update_force_mixing_labels"
    ! get cur indices, labels
    force_mixing_section => section_vals_get_subs_vals3(qmmm_section,"FORCE_MIXING",error=error)
    CALL get_force_mixing_indices(force_mixing_section, cur_indices, cur_labels, error=error)
IF (debug_this_module) PRINT *, "BOB got cur_indices ",SIZE(cur_indices)
IF (debug_this_module) PRINT *, "BOB got cur_labels ",SIZE(cur_labels)

    ! read from input 
    ![NB] breakable bonds will come from here, too
    NULLIFY(r_core, r_qm, r_buf, adaptive_exclude_molecules,broken_bonds)
    CALL section_vals_val_get(force_mixing_section,"R_CORE",r_vals=r_core,error=error)
    CALL section_vals_val_get(force_mixing_section,"R_QM",r_vals=r_qm,error=error)
    CALL section_vals_val_get(force_mixing_section,"QM_EXTENDED_SEED_IS_ONLY_CORE_LIST",&
      l_val=QM_extended_seed_is_core_list,error=error)
    CALL section_vals_val_get(force_mixing_section,"R_BUF",r_vals=r_buf,error=error)
    CALL section_vals_val_get(force_mixing_section,"MAX_N_QM",i_val=max_n_qm,error=error)

    CALL section_vals_val_get(force_mixing_section,"ADAPTIVE_EXCLUDE_MOLECULES",n_rep_val=n_rep_exclude,error=error)
    IF (n_rep_exclude > 0) THEN
       CALL section_vals_val_get(force_mixing_section,"ADAPTIVE_EXCLUDE_MOLECULES",c_vals=adaptive_exclude_molecules,error=error)
    ENDIF
    ![NB] need to read real list from input
    ! should be 2xN_bb integer arrays, with (1,:) indices of inside atoms, and (2,:) indices of outside atoms
    ! maybe also breakable_bond_types, with _atomic numbers_ of inside/outside atoms?
    ! separate lists for core/buffer?
    ALLOCATE(breakable_bonds(2,1))
    breakable_bonds(1:2,1) = 0

    ! get particles, molecules
    NULLIFY(particles,molecules)
    CALL cp_subsys_get(subsys=subsys, particles=particles,molecules_new=molecules, error=error)
    particle_set => particles%els
    molecule_set => molecules%els

    natoms = SIZE(particle_set)

    ! initialize new indices, labels, and new_full_labels
    NULLIFY(new_indices, new_labels)
    CALL reallocate(new_indices,1,SIZE(cur_indices))
    CALL reallocate(new_labels,1,SIZE(cur_labels))
    new_indices = 0
    new_labels = force_mixing_label_none
    ALLOCATE(new_full_labels(natoms))
    new_full_labels = force_mixing_label_none

    ! neighbor list for various hysteretic distance calls
    NULLIFY(cell)
    CALL force_env_get(force_env, cell=cell, error=error)
    NULLIFY(nlist)
    CALL make_neighbor_list(force_mixing_section, subsys, cell, MAX(r_core(2), r_qm(2), r_buf(2)), nlist, error)

    ! create labels for core_list from QM_KIND
    NULLIFY(mm_index_entry)
    qm_kind_section => section_vals_get_subs_vals3(qmmm_section,"QM_KIND",error=error)
    CALL section_vals_get(qm_kind_section,n_repetition=n_rep_section,error=error)
    n_new = 0
    DO i_rep_section=1,n_rep_section
      CALL section_vals_val_get(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section,n_rep_val=n_rep_val,error=error)
      DO i_rep_val=1,n_rep_val
        CALL section_vals_val_get(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section,i_rep_val=i_rep_val, &
          i_vals=mm_index_entry, error=error)
        DO ip=1, SIZE(mm_index_entry)
          CALL add_new_label(mm_index_entry(ip), force_mixing_label_QM_core_list, n_new, new_indices, new_labels, &
                             new_full_labels, max_n_qm)
        END DO ! ip
      END DO ! i_rep_val
    END DO ! i_rep_section

IF (debug_this_module) PRINT *, "BOB core_list new_indices ",new_indices(1:n_new)
IF (debug_this_module) PRINT *, "BOB core_list new_labels ",new_labels(1:n_new)

    ! create labels for non adaptive QM and buffer regions from *_NON_ADAPTIVE&QM_KIND sections
    non_adaptive_section => section_vals_get_subs_vals(qmmm_section,"FORCE_MIXING%QM_NON_ADAPTIVE",&
      can_return_null=.TRUE.,error=error)
    IF (ASSOCIATED(non_adaptive_section)) THEN
       qm_kind_section => section_vals_get_subs_vals3(non_adaptive_section,"QM_KIND",error=error)
       CALL section_vals_get(qm_kind_section,n_repetition=n_rep_section,error=error)
       DO i_rep_section=1,n_rep_section
         CALL section_vals_val_get(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section,n_rep_val=n_rep_val,error=error)
         DO i_rep_val=1,n_rep_val
           CALL section_vals_val_get(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section,i_rep_val=i_rep_val, &
             i_vals=mm_index_entry, error=error)
           DO ip=1, SIZE(mm_index_entry)
             CALL add_new_label(mm_index_entry(ip), force_mixing_label_QM_dynamics_list, n_new, new_indices, new_labels, &
                                new_full_labels, max_n_qm)
           END DO ! ip
         END DO ! i_rep_val
       END DO ! i_rep_section
    ENDIF
IF (debug_this_module) PRINT *, "BOB core_list + non adaptive QM new_indices ",new_indices(1:n_new)
IF (debug_this_module) PRINT *, "BOB core_list + non adaptive QM new_labels ",new_labels(1:n_new)
    non_adaptive_section => section_vals_get_subs_vals(qmmm_section,"FORCE_MIXING%BUFFER_NON_ADAPTIVE",&
      can_return_null=.TRUE.,error=error)
    IF (ASSOCIATED(non_adaptive_section)) THEN
       qm_kind_section => section_vals_get_subs_vals3(non_adaptive_section,"QM_KIND",error=error)
       CALL section_vals_get(qm_kind_section,n_repetition=n_rep_section,error=error)
       DO i_rep_section=1,n_rep_section
         CALL section_vals_val_get(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section,n_rep_val=n_rep_val,error=error)
         DO i_rep_val=1,n_rep_val
           CALL section_vals_val_get(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section,i_rep_val=i_rep_val, &
             i_vals=mm_index_entry, error=error)
           DO ip=1, SIZE(mm_index_entry)
             CALL add_new_label(mm_index_entry(ip), force_mixing_label_buffer_list, n_new, new_indices, new_labels, &
                                new_full_labels, max_n_qm)
           END DO ! ip
         END DO ! i_rep_val
       END DO ! i_rep_section
    ENDIF

IF (debug_this_module) PRINT *, "BOB core_list + non adaptive QM+buffer new_indices ",new_indices(1:n_new)
IF (debug_this_module) PRINT *, "BOB core_list + non adaptive QM+buffer new_labels ",new_labels(1:n_new)

    ! allocate and initialize full atom set labels for hysteretic loops
    ALLOCATE(nearest_dist(natoms))

    ! orig_full_labels is full array (natoms) with orig labels
    ALLOCATE(orig_full_labels(natoms))
    orig_full_labels = force_mixing_label_none
    orig_full_labels(cur_indices(:)) = cur_labels(:)

    ! hysteretically set QM core from QM_core_list and radii, whole molecule
    ![NB] need to replace all the whole molecule stuff with pad to breakable bonds. not quite done
    ! (need intra molecule bond info, which isn't available for QM molecules yet)

    ! find molecule of atom
    ALLOCATE(molecule_of_atom(natoms))
    DO im=1, SIZE(molecule_set)
      molecule_of_atom(molecule_set(im)%first_atom:molecule_set(im)%last_atom) = im
    END DO

    ! add core using hysteretic selection(core_list, r_core) + unbreakable bonds
    CALL add_layer_hysteretically(nlist, particle_set, cell, nearest_dist, &
        orig_full_labels, new_full_labels, n_new, new_indices, new_labels, &
        force_mixing_label_QM_core_list, force_mixing_label_QM_core_list, force_mixing_label_QM_core, r_core, &
        max_n_qm, adaptive_exclude_molecules, .TRUE., breakable_bonds, molecule_of_atom, molecule_set, broken_bonds)
    ![NB] should actually pass this back for making link sections?
    DEALLOCATE(broken_bonds)

IF (debug_this_module) PRINT *, "BOB core new_indices ",new_indices(1:n_new)
IF (debug_this_module) PRINT *, "BOB core new_labels ",new_labels(1:n_new)

    ![NB] need more sophisticated QM extended, buffer rules

    ! add QM using hysteretic selection (core_list, r_qm) + unbreakable bonds
IF (debug_this_module) PRINT *, "BOB QM_extended_seed_is_core_list ",QM_extended_seed_is_core_list
    IF (QM_extended_seed_is_core_list) THEN
      QM_extended_seed_min_label_val = force_mixing_label_QM_core_list
    ELSE ! QM region seed is all of core, not just core list + unbreakable bonds
      QM_extended_seed_min_label_val = force_mixing_label_QM_core
    ENDIF
    CALL add_layer_hysteretically(nlist, particle_set, cell, nearest_dist, &
        orig_full_labels, new_full_labels, n_new, new_indices, new_labels, &
        QM_extended_seed_min_label_val, force_mixing_label_QM_core_list, &
        force_mixing_label_QM_dynamics, r_qm, &
        max_n_qm, adaptive_exclude_molecules, .TRUE., breakable_bonds, molecule_of_atom, molecule_set)

IF (debug_this_module) PRINT *, "BOB extended new_indices ",new_indices(1:n_new)
IF (debug_this_module) PRINT *, "BOB extended new_labels ",new_labels(1:n_new)

    ! add buffer using hysteretic selection (>= QM extended, r_buf) + unbreakable bonds
    CALL add_layer_hysteretically(nlist, particle_set, cell, nearest_dist, &
        orig_full_labels, new_full_labels, n_new, new_indices, new_labels, &
        force_mixing_label_QM_dynamics, force_mixing_label_QM_core_list, force_mixing_label_buffer, r_buf, &
        max_n_qm, adaptive_exclude_molecules, .TRUE., breakable_bonds, molecule_of_atom, molecule_set, broken_bonds)
    ![NB] should actually pass this back for making link sections?
    DEALLOCATE(broken_bonds)

IF (debug_this_module) PRINT *, "BOB buffer new_indices ",new_indices(1:n_new)
IF (debug_this_module) PRINT *, "BOB buffer new_labels ",new_labels(1:n_new)

     DEALLOCATE(molecule_of_atom, breakable_bonds)
     DEALLOCATE(nearest_dist)

     IF (PRESENT(labels_changed)) labels_changed = ANY (new_full_labels /= orig_full_labels)

     DEALLOCATE(orig_full_labels)
     DEALLOCATE(new_full_labels)

     ! reduce new indices, labels to actually used size
     CALL reallocate(new_indices,1,n_new)
     CALL reallocate(new_labels,1,n_new)

     ! save info in input structure
     restart_section => section_vals_get_subs_vals(qmmm_section, "FORCE_MIXING%RESTART_INFO", error=error)
     CALL section_vals_get(restart_section,explicit=explicit,error=error)
     IF (explicit) CALL section_vals_remove_values(restart_section,error)
     CALL section_vals_val_set(restart_section,"INDICES",i_vals_ptr=new_indices,error=error)
     CALL section_vals_val_set(restart_section,"LABELS",i_vals_ptr=new_labels,error=error)

     DEALLOCATE(cur_indices, cur_labels)
     CALL fist_neighbor_deallocate(nlist,error)

     ![NB] perhap be controlled by some &PRINT section?
     IF (force_env%para_env%ionode) THEN
       WRITE(unit=*, fmt='(A,A,I6,A,I5,A,I5,A,I5)') &
          "QMMM FORCE MIXING final count (not including links): ",&
          " N_QM core_list ", COUNT(new_labels == force_mixing_label_QM_core_list), &
          " N_QM core ",      COUNT(new_labels == force_mixing_label_QM_core), &
          " N_QM extended ",  COUNT(new_labels == force_mixing_label_QM_dynamics .OR. &
                                    new_labels == force_mixing_label_QM_dynamics_list), &
          " N_QM buffered ",  COUNT(new_labels == force_mixing_label_buffer .OR. &
                                    new_labels == force_mixing_label_buffer_list)
     ENDIF

   END SUBROUTINE update_force_mixing_labels

! *****************************************************************************
!> \brief ...
!> \param ip ...
!> \param label ...
!> \param n_new ...
!> \param new_indices ...
!> \param new_labels ...
!> \param new_full_labels ...
!> \param max_n_qm ...
! *****************************************************************************
  SUBROUTINE add_new_label(ip, label, n_new, new_indices, new_labels, new_full_labels, max_n_qm)
    INTEGER                                  :: ip, label, n_new
    INTEGER, POINTER                         :: new_indices(:), new_labels(:)
    INTEGER                                  :: new_full_labels(:), max_n_qm

    CHARACTER(len=*), PARAMETER :: routineN = 'add_new_label', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, old_index

    IF (new_full_labels(ip) > force_mixing_label_none) THEN ! already marked, just change mark
         old_index = -1
         DO i=1, n_new
            IF (new_indices(i) == ip) THEN
               old_index = i
               EXIT
            ENDIF
         END DO
         CALL cp_assert(old_index > 0, cp_failure_level,cp_assertion_failed,routineP,&
              "add_new_label found atom with a label already set, but not in new_indices array"//&
CPSourceFileRef)
         new_labels(old_index) = label
    ELSE
         n_new = n_new + 1
         CALL cp_assert(n_new <= max_n_qm, cp_failure_level,cp_assertion_failed,routineP,&
              "add_new_label tried to add more atoms than allowed by &FORCE_MIXING&MAX_N_QM!"//&
CPSourceFileRef)
         IF (n_new > SIZE(new_indices)) CALL reallocate(new_indices,1,n_new+9)
         IF (n_new > SIZE(new_labels))  CALL reallocate(new_labels,1,n_new+9)
         new_indices(n_new) = ip
         new_labels(n_new) = label
    ENDIF
    new_full_labels(ip) = label
  END SUBROUTINE add_new_label


! *****************************************************************************
!> \brief ...
!> \param nlist ...
!> \param particle_set ...
!> \param cell ...
!> \param nearest_dist ...
!> \param orig_full_labels ...
!> \param new_full_labels ...
!> \param n_new ...
!> \param new_indices ...
!> \param new_labels ...
!> \param seed_min_label_val ...
!> \param seed_max_label_val ...
!> \param set_label_val ...
!> \param r_inout ...
!> \param max_n_qm ...
!> \param adaptive_exclude_molecules ...
!> \param pad_to_breakable_bonds ...
!> \param breakable_bonds ...
!> \param molecule_of_atom ...
!> \param molecule_set ...
!> \param broken_bonds ...
! *****************************************************************************
  SUBROUTINE add_layer_hysteretically(nlist, particle_set, cell, nearest_dist, &
    orig_full_labels, new_full_labels, n_new, new_indices, new_labels, &
    seed_min_label_val, seed_max_label_val, set_label_val, r_inout, max_n_qm, &
    adaptive_exclude_molecules, &
    pad_to_breakable_bonds, breakable_bonds, molecule_of_atom, molecule_set, &
    broken_bonds)
    TYPE(fist_neighbor_type), POINTER        :: nlist
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp)                                 :: nearest_dist(:)
    INTEGER                                  :: orig_full_labels(:), &
                                                new_full_labels(:), n_new
    INTEGER, POINTER                         :: new_indices(:), new_labels(:)
    INTEGER                                  :: seed_min_label_val, &
                                                seed_max_label_val, &
                                                set_label_val
    REAL(dp)                                 :: r_inout(2)
    INTEGER                                  :: max_n_qm
    CHARACTER(len=*), POINTER                :: adaptive_exclude_molecules(:)
    LOGICAL                                  :: pad_to_breakable_bonds
    INTEGER, OPTIONAL, POINTER               :: breakable_bonds(:,:), &
                                                molecule_of_atom(:)
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    INTEGER, OPTIONAL, POINTER               :: broken_bonds(:)

    CHARACTER(len=*), PARAMETER :: routineN = 'add_layer_hysteretically', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_ind, im, im_exclude, ip, &
                                                ipair, ipairkind, j_ind
    LOGICAL :: adaptive_exclude, i_in_new_seed, i_outside_new_seed, &
      j_in_new_seed, j_outside_new_seed, molec_in_inner, molec_in_outer
    REAL(dp)                                 :: r_ij(3), r_ij_mag

IF (debug_this_module) PRINT *,"BOB adding hysteretically seed ",seed_min_label_val,seed_max_label_val,&
" set ",set_label_val," r ",r_inout
    ! calculate nearest dist from each atom outside of new seed to nearest atom inside of new seed
    nearest_dist=HUGE(1.0_dp)
    ! loop over pairs of all kinds in random order
    DO ipairkind=1,SIZE(nlist%neighbor_kind_pairs)
      DO ipair=1,nlist%neighbor_kind_pairs(ipairkind)%npairs

        i_ind = nlist%neighbor_kind_pairs(ipairkind)%list(1,ipair)
        j_ind = nlist%neighbor_kind_pairs(ipairkind)%list(2,ipair)

        i_in_new_seed = (new_full_labels(i_ind) >= seed_min_label_val .AND. new_full_labels(i_ind) <= seed_max_label_val)
        i_outside_new_seed = (new_full_labels(i_ind) < seed_min_label_val)
        j_in_new_seed = (new_full_labels(j_ind) >= seed_min_label_val .AND. new_full_labels(j_ind) <= seed_max_label_val)
        j_outside_new_seed = (new_full_labels(j_ind) < seed_min_label_val)

        IF ((i_in_new_seed .AND. j_outside_new_seed) .OR. (j_in_new_seed .AND. i_outside_new_seed)) THEN
          r_ij = pbc(particle_set(i_ind)%r - particle_set(j_ind)%r, cell)
          r_ij_mag = SQRT(SUM(r_ij**2))
          IF (i_in_new_seed .AND. j_outside_new_seed .AND. (r_ij_mag < nearest_dist(j_ind))) THEN
            nearest_dist(j_ind) = r_ij_mag
          ENDIF
          IF (j_in_new_seed .AND. i_outside_new_seed .AND. (r_ij_mag < nearest_dist(i_ind))) THEN
            nearest_dist(i_ind) = r_ij_mag
          ENDIF
        ENDIF

      END DO
    END DO
! IF (debug_this_module) PRINT *,"BOB nearest_dist ", nearest_dist

![NB] this is whole molecule.  Should be replaced with labeling of individual atoms + 
! pad_to_breakable_bonds (below), but QM molecule bond information isn't available yet
DO im=1, SIZE(molecule_set)
! IF (debug_this_module) PRINT *, "BOB check molecule ",im,trim(molecule_set(im)%molecule_kind%name)," atoms ",&
! molecule_set(im)%first_atom,molecule_set(im)%last_atom
  IF (ASSOCIATED(adaptive_exclude_molecules)) THEN
     adaptive_exclude = .FALSE.
     DO im_exclude=1, SIZE(adaptive_exclude_molecules)
        IF ( TRIM(molecule_set(im)%molecule_kind%name) == TRIM(adaptive_exclude_molecules(im_exclude)) .OR. &
             TRIM(molecule_set(im)%molecule_kind%name) == '_QM_'//TRIM(adaptive_exclude_molecules(im_exclude)) ) &
          adaptive_exclude = .TRUE.
     END DO
     IF (adaptive_exclude) CYCLE
  ENDIF
  molec_in_inner = ANY(nearest_dist(molecule_set(im)%first_atom:molecule_set(im)%last_atom) <= r_inout(1))
  molec_in_outer = ANY(nearest_dist(molecule_set(im)%first_atom:molecule_set(im)%last_atom) <= r_inout(2))
! IF (debug_this_module) PRINT *, "BOB got molec_in_inner ", molec_in_inner, " outer ", molec_in_outer
  IF (molec_in_inner) THEN
    DO ip=molecule_set(im)%first_atom, molecule_set(im)%last_atom
      ! labels are being rebuild from scratch, so never overwrite new label that's higher level
      IF (new_full_labels(ip) < set_label_val) &
        CALL add_new_label(ip, set_label_val, n_new, new_indices, new_labels, new_full_labels, max_n_qm)
    END DO
  ELSE IF (molec_in_outer) THEN
    IF (ANY(orig_full_labels(molecule_set(im)%first_atom:molecule_set(im)%last_atom) >= set_label_val)) THEN
      DO ip=molecule_set(im)%first_atom, molecule_set(im)%last_atom
        ! labels are being rebuild from scratch, so never overwrite new label that's higher level
        IF (new_full_labels(ip) < set_label_val) &
          CALL add_new_label(ip, set_label_val, n_new, new_indices, new_labels, new_full_labels, max_n_qm)
      END DO
    ENDIF
  ENDIF
END DO
IF (PRESENT(broken_bonds)) CALL reallocate(broken_bonds,1,0)

!!![NB] enable this (and remove previous bit) once intra-QM-molecule bond information
!!! is available
!!!    ! label new layer
!!!    ! loop over every atom
!!!    DO ip=1, SIZE(new_full_labels)
!!!      ! labels are being rebuild from scratch, so never overwrite new label that's higher level
!!!      if (new_full_labels(ip) > set_label_val) cycle
!!!
!!!if (debug_this_module) print *,"BOB atom ",ip," new_full_label ",new_full_labels(ip), &
!!!" is not in seed, check for for adding"
!!!      ![NB] is this logic right?  maybe it is
!!!      IF ( nearest_dist(ip) <= r_inout(1) .OR. &
!!!          (nearest_dist(ip) <= r_inout(2) .AND. orig_full_labels(ip) >= set_label_val) ) THEN 
!!!          ! atom in inner radius
!!!          ! or
!!!          ! atom in outer radius and was previously set to set_label_val or higher 
!!!if (debug_this_module) print *,"BOB add label for atom ", ip
!!!          call add_new_label(ip, set_label_val, n_new, new_indices, new_labels, new_full_labels, max_n_qm)
!!!      ENDIF
!!!    END DO ! ip
!!!
!!!    if (present(broken_bonds)) n_broken_bonds = 0
!!!    if (pad_to_breakable_bonds) then
!!!if (debug_this_module) print *,"BOB pad to breakable bonds"
!!!      CALL cp_assert(PRESENT(molecule_of_atom),cp_failure_level,cp_assertion_failed,routineP,&
!!!           "pad_to_breakable_bonds is true, but molecule_of_atom not present!"//&
!!!CPSourceFileRef)
!!!      CALL cp_assert(PRESENT(molecule_set),cp_failure_level,cp_assertion_failed,routineP,&
!!!           "pad_to_breakable_bonds is true, but molecule_set not present!"//&
!!!CPSourceFileRef)
!!!      CALL cp_assert(PRESENT(breakable_bonds),cp_failure_level,cp_assertion_failed,routineP,&
!!!           "whole_molecule is true, but breakable_bonds not present!"//&
!!!CPSourceFileRef)
!!!
!!!      orig_n_new = n_new
!!!      ! loop over all recently marked atoms
!!!      do ip=1, orig_n_new
!!!        if (new_labels(ip) /= set_label_val) cycle
!!!
!!!if (debug_this_module) then
!!!  print *, "BOB molecule set"
!!!  do im=1, size(molecule_set)
!!!    print *, "BOB molecule im ",im, " #bonds ",size(molecule_set(im)%molecule_kind%bond_list), " name ", &
!!!      trim(molecule_set(im)%molecule_kind%name), molecule_set(im)%molecule_kind%kind_number
!!!    do ib=1, size(molecule_set(im)%molecule_kind%bond_list)
!!!      print *, "BOB    ",molecule_set(im)%molecule_kind%bond_list(ib)%a + molecule_set(im)%first_atom - 1, &
!!!                     molecule_set(im)%molecule_kind%bond_list(ib)%b + molecule_set(im)%first_atom - 1
!!!    end do
!!!  end do
!!!endif
!!!
!!!if (debug_this_module) print *,"BOB check atom local index ",ip," full index ", new_indices(ip)
!!!        ! atom has just been set to set_label_val, check the rest of its molecule
!!!        im = molecule_of_atom(new_indices(ip))
!!!
!!!if (debug_this_module) print *,"in molecule ",im
!!!        any_additions = .true.
!!!        do while (any_additions)
!!!if (debug_this_module) print *,"start loop checking for additions"
!!!          any_additions = .false.
!!!          ! loop over bonds in this molecule in random order
!!!          do ib=1, size(molecule_set(im)%molecule_kind%bond_list)
!!!            bond_atom_a = molecule_set(im)%molecule_kind%bond_list(ib)%a + molecule_set(im)%first_atom - 1
!!!            bond_atom_b = molecule_set(im)%molecule_kind%bond_list(ib)%b + molecule_set(im)%first_atom - 1
!!!if (debug_this_module) print *,"bond ", ib, " between a,b ", bond_atom_a,bond_atom_b
!!!
!!!            if (new_full_labels(bond_atom_a) == set_label_val .and. new_full_labels(bond_atom_b) < set_label_val) then
!!!              ! a is inside, b is outside
!!!              bond_atom_inside = bond_atom_a
!!!              bond_atom_outside = bond_atom_b
!!!if (debug_this_module) print *,"  a is inside, b is outside"
!!!            else if (new_full_labels(bond_atom_b) == set_label_val .and. new_full_labels(bond_atom_a) < set_label_val) then
!!!              ! b is inside, a is outside
!!!              bond_atom_inside = bond_atom_b
!!!              bond_atom_outside = bond_atom_a
!!!if (debug_this_module) print *,"  b is inside, a is outside"
!!!            else
!!!              ! inside-inside or outside-outside
!!!              cycle
!!!            endif
!!!
!!!            ! check for breakable bonds
!!!            breakable = .false.
!!!            do ibb=1,size(breakable_bonds,2)
!!!              if ((breakable_bonds(1,ibb) == bond_atom_inside) .and. breakable_bonds(2,ibb) == bond_atom_outside) then
!!!                ! this bond is OK to break.
!!!                breakable = .true.
!!!                if (set_label_val == force_mixing_label_buffer) then ! set termination mark on outside atom b
!!!                  call add_new_label(bond_atom_outside, force_mixing_label_termination, n_new, new_indices, new_labels, &
!!!                    new_full_labels, max_n_qm)
!!!                endif
!!!                if (present(broken_bonds)) then ! add a new one
!!!                  n_broken_bonds = n_broken_bonds + 1
!!!                  CALL reallocate(broken_bonds,1,n_broken_bonds)
!!!                  broken_bonds(n_broken_bonds) = ibb
!!!                endif
!!!              endif
!!!            end do ! ibb
!!!if (debug_this_module) print *,"  breakable ", breakable
!!!
!!!            if (.not. breakable) then ! add correct atom
!!!              call add_new_label(bond_atom_outside, set_label_val, n_new, new_indices, new_labels, new_full_labels, max_n_qm)
!!!              any_additions = .true.
!!!if (debug_this_module) print *, " not breakable, adding outside atom ", bond_atom_outside
!!!            endif
!!!
!!!          end do ! ib
!!!        end do ! any_additions
!!!
!!!      end do ! ip
!!!
!!!      if (present(broken_bonds)) CALL reallocate(broken_bonds,1,n_broken_bonds)
!!!    endif ! pad_to_breakable_bonds

  END SUBROUTINE add_layer_hysteretically

! *****************************************************************************
!> \brief ...
!> \param force_mixing_section ...
!> \param subsys ...
!> \param cell ...
!> \param r_max ...
!> \param nlist ...
!> \param error ...
! *****************************************************************************
   SUBROUTINE make_neighbor_list(force_mixing_section, subsys, cell, r_max, nlist, error)
    TYPE(section_vals_type), POINTER         :: force_mixing_section
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp)                                 :: r_max
    TYPE(fist_neighbor_type), POINTER        :: nlist
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_neighbor_list', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kind_name
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: kind_name_a(:)
    INTEGER                                  :: ik
    LOGICAL                                  :: skip_kind
    REAL(dp), ALLOCATABLE                    :: r_max_a(:,:), r_minsq_a(:,:)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

     ALLOCATE(r_max_a(SIZE(subsys%atomic_kinds%els),SIZE(subsys%atomic_kinds%els)))
     ALLOCATE(r_minsq_a(SIZE(subsys%atomic_kinds%els),SIZE(subsys%atomic_kinds%els)))
     r_max_a = r_max
     r_minsq_a = EPSILON(1.0_dp)

      ! save kind names
      ALLOCATE(kind_name_a(SIZE(subsys%atomic_kinds%els)))
      DO ik=1, SIZE(subsys%atomic_kinds%els)
         atomic_kind => subsys%atomic_kinds%els(ik)
         CALL get_atomic_kind(atomic_kind=atomic_kind,name=kind_name)
         kind_name_a(ik) = kind_name
      END DO

     ! overwrite kind names so that none are QM, and so excluding QM-QM interactions
     ! (which is not what we want) will not happen
     DO ik=1, SIZE(subsys%atomic_kinds%els)
       atomic_kind => subsys%atomic_kinds%els(ik)
       CALL get_atomic_kind(atomic_kind=atomic_kind,name=kind_name)
       ! when atom is QM atom, kind_name is replaced with original
       ! mm kind name, and return status is logical .TRUE.
       skip_kind=qmmm_ff_precond_only_qm(kind_name)
       CALL set_atomic_kind(atomic_kind=atomic_kind,name=kind_name)
     END DO

     NULLIFY(nlist)
     CALL build_fist_neighbor_lists(subsys%atomic_kinds%els, subsys%particles%els, &
          cell=cell, r_max=r_max_a, r_minsq=r_minsq_a, &
          ei_scale14=1.0_dp, vdw_scale14=1.0_dp, nonbonded=nlist, &
          para_env=subsys%para_env, build_from_scratch=.TRUE., geo_check=.FALSE., &
          mm_section=force_mixing_section, error=error)

      DEALLOCATE(r_max_a, r_minsq_a)

      ! restore kind names
      DO ik=1, SIZE(subsys%atomic_kinds%els)
         CALL set_atomic_kind(atomic_kind=atomic_kind,name=kind_name_a(ik))
      END DO
     DEALLOCATE(kind_name_a)

   END SUBROUTINE make_neighbor_list

! *****************************************************************************
!> \brief ...
!> \param subsys ...
!> \param qmmm_section ...
!> \param qmmm_core_section ...
!> \param qmmm_extended_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2012 created [noam]
!> \author Noam Bernstein
! *****************************************************************************
   SUBROUTINE setup_force_mixing_qmmm_sections(subsys, qmmm_section, qmmm_core_section, qmmm_extended_section, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: qmmm_section, &
                                                qmmm_core_section, &
                                                qmmm_extended_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'setup_force_mixing_qmmm_sections', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length), &
      POINTER                                :: elem_mapping(:,:), &
                                                elem_mapping_entry(:)
    INTEGER :: delta_charge, i_rep_section_core, i_rep_section_extended, &
      i_rep_val_core, i_rep_val_extended, ielem, ip, n_elements
    INTEGER, POINTER                         :: cur_indices(:), cur_labels(:)
    LOGICAL                                  :: mapped, new_element_core, &
                                                new_element_extended
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    TYPE(section_vals_type), POINTER :: buffer_non_adaptive_section, &
      dup_link_section, force_mixing_section, link_section, qm_kind_section

!
! create new qmmm sections for core and extended

    NULLIFY(qmmm_core_section, qmmm_extended_section)
    CALL section_vals_duplicate(qmmm_section, qmmm_core_section, error=error)
    CALL section_vals_duplicate(qmmm_section, qmmm_extended_section, error=error)

    ! remove LINKs (specified by user for core) from extended
    link_section => section_vals_get_subs_vals(qmmm_extended_section, "LINK", can_return_null=.TRUE., error=error)
    IF (ASSOCIATED(link_section)) THEN
      CALL section_vals_remove_values(link_section,error)
    END IF
    ! for LINKs to be added to extended
    buffer_non_adaptive_section => section_vals_get_subs_vals(qmmm_extended_section, "FORCE_MIXING%BUFFER_NON_ADAPTIVE", &
      can_return_null=.TRUE., error=error)
    link_section => section_vals_get_subs_vals(buffer_non_adaptive_section, "LINK", can_return_null=.TRUE., error=error)
    IF (ASSOCIATED(link_section)) THEN
       NULLIFY(dup_link_section)
       CALL section_vals_duplicate(link_section, dup_link_section, error=error)
       CALL section_vals_set_subs_vals(qmmm_extended_section, "LINK", dup_link_section, error=error)
       CALL section_vals_release(dup_link_section, error=error)
    ENDIF

IF (debug_this_module) THEN
link_section => section_vals_get_subs_vals(qmmm_core_section, "LINK", can_return_null=.TRUE., error=error)
PRINT *, "core section has LINKs ",ASSOCIATED(link_section)
CALL section_vals_write(link_section,unit_nr=6,error=error)
link_section => section_vals_get_subs_vals(qmmm_extended_section, "LINK", can_return_null=.TRUE., error=error)
PRINT *, "extended section has LINKs ",ASSOCIATED(link_section)
CALL section_vals_write(link_section,unit_nr=6,error=error)
ENDIF

    force_mixing_section => section_vals_get_subs_vals(qmmm_section,"FORCE_MIXING",error=error)

    ! get QM_KIND_ELEMENT_MAPPING
    CALL section_vals_val_get(force_mixing_section,"QM_KIND_ELEMENT_MAPPING",n_rep_val=n_elements,error=error)
    ALLOCATE(elem_mapping(2,n_elements))
    DO ielem=1,n_elements
      CALL section_vals_val_get(force_mixing_section,"QM_KIND_ELEMENT_MAPPING",i_rep_val=ielem,c_vals=elem_mapping_entry,&
        error=error)
      elem_mapping(1:2,ielem) = elem_mapping_entry(1:2)
    END DO

    ! get CUR_INDICES, CUR_LABELS
    CALL get_force_mixing_indices(force_mixing_section, cur_indices, cur_labels, error=error)
    CALL cp_assert(SIZE(cur_indices) > 0,cp_failure_level,cp_assertion_failed,routineP,&
         "cur_indices is empty, found no QM atoms"//&
CPSourceFileRef)

IF (debug_this_module) THEN
PRINT *, "cur_indices ",cur_indices
PRINT *, "cur_labels ",cur_labels
ENDIF

    ! loop through elements and atoms, and set up new QM_KIND sections
    particles => subsys%particles%els

    DO ip=1, SIZE(cur_indices)
      IF (cur_labels(ip) > force_mixing_label_none .AND. cur_labels(ip) < force_mixing_label_QM_core_list .AND. &
          cur_labels(ip) /= force_mixing_label_termination) THEN
          mapped = .FALSE.
          DO ielem=1,n_elements
              IF (TRIM(particles(cur_indices(ip))%atomic_kind%element_symbol) == TRIM(elem_mapping(1,ielem))) THEN
                 mapped = .TRUE.
                 EXIT
              ENDIF
          END DO
          CALL cp_assert(mapped,cp_failure_level,cp_assertion_failed,routineP,&
               "Force-mixing failed to find QM_KIND mapping for atom of type "//&
               TRIM(particles(cur_indices(ip))%atomic_kind%element_symbol)//"! "//&
CPSourceFileRef)
      ENDIF
    END DO

    ! pre-existing QM_KIND section specifies list of core atom
    qm_kind_section => section_vals_get_subs_vals3(qmmm_section,"QM_KIND",error=error)
    CALL section_vals_get(qm_kind_section,n_repetition=i_rep_section_core,error=error)
    CALL cp_assert(i_rep_section_core > 0,cp_failure_level,cp_assertion_failed,routineP,&
         "Force-mixing QM didn't find any QM_KIND sections, so no core specified!"//&
CPSourceFileRef)
    i_rep_section_extended = i_rep_section_core
    DO ielem=1,n_elements
      new_element_core = .TRUE.
      new_element_extended = .TRUE.
      DO ip=1, SIZE(cur_indices) ! particles with label
        IF (TRIM(particles(cur_indices(ip))%atomic_kind%element_symbol) /= TRIM(elem_mapping(1,ielem))) CYCLE
        ! extended
        ! if current particle is some sort of QM atom, and not in core list
        ! (those the user gave explicit QM_KIND sections for), and not a
        ! termination atom, need to make a QM_KIND section for it
        IF (cur_labels(ip) > force_mixing_label_none .AND. &
            cur_labels(ip) /= force_mixing_label_QM_core_list .AND. &
            cur_labels(ip) /= force_mixing_label_termination) THEN
            qm_kind_section => section_vals_get_subs_vals3(qmmm_extended_section,"QM_KIND",error=error)
            IF (new_element_extended) THEN ! add new QM_KIND section for this element
              i_rep_section_extended = i_rep_section_extended + 1
              CALL section_vals_add_values(qm_kind_section,error=error)
              CALL section_vals_val_set(qm_kind_section,"_SECTION_PARAMETERS_",i_rep_section=i_rep_section_extended, &
                c_val=elem_mapping(2,ielem),error=error)
              i_rep_val_extended = 0
              new_element_extended = .FALSE.
            ENDIF
            i_rep_val_extended = i_rep_val_extended + 1
            CALL section_vals_val_set(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section_extended, &
              i_rep_val=i_rep_val_extended, i_val=cur_indices(ip),error=error)
        ENDIF ! is a non-termination QM atom

        ! core
        ! if current particle is a core QM atom, and not in core list (those the user
        ! gave explicit QM_KIND sections for, need to make a QM_KIND section for it
        IF (cur_labels(ip) == force_mixing_label_QM_core) THEN
            qm_kind_section => section_vals_get_subs_vals3(qmmm_core_section,"QM_KIND",error=error)
            IF (new_element_core) THEN ! add new QM_KIND section for this element
              i_rep_section_core = i_rep_section_core + 1
              CALL section_vals_add_values(qm_kind_section,error=error)
              CALL section_vals_val_set(qm_kind_section,"_SECTION_PARAMETERS_",i_rep_section=i_rep_section_core, &
                c_val=elem_mapping(2,ielem),error=error)
              i_rep_val_core = 0
              new_element_core = .FALSE.
            ENDIF
            i_rep_val_core = i_rep_val_core + 1
            CALL section_vals_val_set(qm_kind_section,"MM_INDEX",i_rep_section=i_rep_section_core, &
              i_rep_val=i_rep_val_core, i_val=cur_indices(ip),error=error)
        ENDIF ! is a non-termination QM atom

      END DO ! atom index ip
    END DO ! element index ielem

    CALL section_vals_val_get(force_mixing_section,"EXTENDED_DELTA_CHARGE",i_val=delta_charge,error=error)
    CALL section_vals_val_set(qmmm_extended_section,"DELTA_CHARGE",i_val=delta_charge,error=error)

    ![NB] check
    DEALLOCATE(elem_mapping, cur_indices, cur_labels)

IF (debug_this_module) THEN
IF (subsys%para_env%ionode) THEN
PRINT *, "qmmm_core_section"
CALL section_vals_write(qmmm_core_section,unit_nr=6,error=error)
PRINT *, "qmmm_extended_section"
CALL section_vals_write(qmmm_extended_section,unit_nr=6,error=error)
ENDIF
ENDIF

   END SUBROUTINE setup_force_mixing_qmmm_sections

! *****************************************************************************
!> \brief ...
!> \param force_mixing_section ...
!> \param indices ...
!> \param labels ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_force_mixing_indices(force_mixing_section, indices, labels, error)
    TYPE(section_vals_type), POINTER         :: force_mixing_section
    INTEGER, POINTER                         :: indices(:), labels(:)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_force_mixing_indices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_rep_val, n_indices, &
                                                n_labels, n_reps
    INTEGER, POINTER                         :: indices_entry(:), &
                                                labels_entry(:)
    LOGICAL                                  :: explicit
    TYPE(section_vals_type), POINTER         :: restart_section

    NULLIFY(indices, labels)
    restart_section => section_vals_get_subs_vals(force_mixing_section, "RESTART_INFO", error=error)
    CALL section_vals_get(restart_section,explicit=explicit,error=error)
    IF (.NOT. explicit) THEN ! no old indices, labels, return empty arrays
      ALLOCATE(indices(0))
      ALLOCATE(labels(0))
      RETURN
    ENDIF

    ![NB] maybe switch to reallocatable array
    CALL section_vals_val_get(restart_section,"INDICES",n_rep_val=n_reps,error=error)
    n_indices=0
    DO i_rep_val = 1,n_reps
      CALL section_vals_val_get(restart_section,"INDICES",&
        i_rep_val=i_rep_val,i_vals=indices_entry,error=error)
      n_indices = n_indices + SIZE(indices_entry)
    END DO
    ALLOCATE(indices(n_indices))
    n_indices=0
    DO i_rep_val = 1,n_reps
      CALL section_vals_val_get(restart_section,"INDICES",&
        i_rep_val=i_rep_val,i_vals=indices_entry,error=error)
      indices(n_indices+1:n_indices+SIZE(indices_entry)) = indices_entry
      n_indices = n_indices + SIZE(indices_entry)
    END DO

    CALL section_vals_val_get(restart_section,"LABELS",n_rep_val=n_reps,error=error)
    n_labels=0
    DO i_rep_val = 1,n_reps
      CALL section_vals_val_get(restart_section,"LABELS",&
        i_rep_val=i_rep_val,i_vals=labels_entry,error=error)
      n_labels = n_labels + SIZE(labels_entry)
    END DO
    ALLOCATE(labels(n_labels))
    n_labels=0
    DO i_rep_val = 1,n_reps
      CALL section_vals_val_get(restart_section,"LABELS",&
        i_rep_val=i_rep_val,i_vals=labels_entry,error=error)
      labels(n_labels+1:n_labels+SIZE(labels_entry)) = labels_entry
      n_labels = n_labels + SIZE(labels_entry)
    END DO

    CALL cp_assert(n_indices == n_labels, cp_failure_level,cp_assertion_failed,routineP,&
         "got unequal numbers of force_mixing indices and labels!"//&
CPSourceFileRef)
  END SUBROUTINE get_force_mixing_indices

END MODULE qmmm_force_mixing
