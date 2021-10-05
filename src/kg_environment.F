!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for a Kim-Gordon-like partitioning into molecular subunits
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! **************************************************************************************************
MODULE kg_environment
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Andermatt2016,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE external_potential_types,        ONLY: get_potential,&
                                              local_potential_type
   USE input_constants,                 ONLY: kg_tnadd_atomic,&
                                              kg_tnadd_embed,&
                                              kg_tnadd_embed_ri,&
                                              kg_tnadd_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE integration_grid_types,          ONLY: allocate_intgrid,&
                                              integration_grid_type
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kg_vertex_coloring_methods,      ONLY: kg_vertex_coloring
   USE kinds,                           ONLY: dp,&
                                              int_4,&
                                              int_4_size,&
                                              int_8
   USE lri_environment_init,            ONLY: lri_env_basis,&
                                              lri_env_init
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_gather,&
                                              mp_max
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: initialize_atomic_grid
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup,&
                                              write_neighbor_lists
   USE string_utilities,                ONLY: uppercase
   USE task_list_types,                 ONLY: deallocate_task_list
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_environment'

   PUBLIC :: kg_env_create, kg_build_neighborlist, kg_build_subsets

CONTAINS

! **************************************************************************************************
!> \brief Allocates and intitializes kg_env
!> \param qs_env ...
!> \param kg_env the object to create
!> \param qs_kind_set ...
!> \param input ...
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! **************************************************************************************************
   SUBROUTINE kg_env_create(qs_env, kg_env, qs_kind_set, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input

      ALLOCATE (kg_env)
      CALL init_kg_env(qs_env, kg_env, qs_kind_set, input)
   END SUBROUTINE kg_env_create

! **************************************************************************************************
!> \brief Initializes kg_env
!> \param qs_env ...
!> \param kg_env ...
!> \param qs_kind_set ...
!> \param input ...
!> \par History
!>       2012.07 created [Martin Haeufel]
!>       2018.01 TNADD correction {JGH]
!> \author Martin Haeufel
! **************************************************************************************************
   SUBROUTINE init_kg_env(qs_env, kg_env, qs_kind_set, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_kg_env'

      CHARACTER(LEN=10)                                  :: intgrid
      INTEGER                                            :: handle, i, iatom, ib, ikind, iunit, n, &
                                                            na, natom, nbatch, nkind, np, nr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bid
      REAL(KIND=dp)                                      :: load, radb, rmax
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(gto_basis_set_type), POINTER                  :: lri_aux_basis
      TYPE(integration_grid_type), POINTER               :: ig_full, ig_mol
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(section_vals_type), POINTER                   :: lri_section

      CALL timeset(routineN, handle)

      CALL cite_reference(Andermatt2016)

      NULLIFY (para_env)
      NULLIFY (kg_env%sab_orb_full)
      NULLIFY (kg_env%sac_kin)
      NULLIFY (kg_env%subset_of_mol)
      NULLIFY (kg_env%subset)
      NULLIFY (kg_env%tnadd_mat)
      NULLIFY (kg_env%lri_env)
      NULLIFY (kg_env%lri_env1)
      NULLIFY (kg_env%int_grid_atom)
      NULLIFY (kg_env%int_grid_molecules)
      NULLIFY (kg_env%int_grid_full)
      NULLIFY (kg_env%lri_density)
      NULLIFY (kg_env%lri_rho1)

      kg_env%nsubsets = 0

      ! get coloring method settings
      CALL section_vals_val_get(input, "DFT%KG_METHOD%COLORING_METHOD", i_val=kg_env%coloring_method)
      ! get method for nonadditive kinetic energy embedding potential
      CALL section_vals_val_get(input, "DFT%KG_METHOD%TNADD_METHOD", i_val=kg_env%tnadd_method)
      !
      SELECT CASE (kg_env%tnadd_method)
      CASE (kg_tnadd_embed, kg_tnadd_embed_ri)
         ! kinetic energy functional
         kg_env%xc_section_kg => section_vals_get_subs_vals(input, "DFT%KG_METHOD%XC")
         IF (.NOT. ASSOCIATED(kg_env%xc_section_kg)) THEN
            CALL cp_abort(__LOCATION__, &
                          "KG runs require a kinetic energy functional set in &KG_METHOD")
         END IF
      CASE (kg_tnadd_atomic, kg_tnadd_none)
         NULLIFY (kg_env%xc_section_kg)
      CASE DEFAULT
         CPABORT("KG:TNADD METHOD")
      END SELECT

      IF (kg_env%tnadd_method == kg_tnadd_embed_ri) THEN
         ! initialize the LRI environment
         ! Check if LRI_AUX basis is available
         rmax = 0.0_dp
         nkind = SIZE(qs_kind_set)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            NULLIFY (lri_aux_basis)
            CALL get_qs_kind(qs_kind, basis_set=lri_aux_basis, basis_type="LRI_AUX")
            CPASSERT(ASSOCIATED(lri_aux_basis))
            CALL get_gto_basis_set(gto_basis_set=lri_aux_basis, kind_radius=radb)
            rmax = MAX(rmax, radb)
         END DO
         rmax = 1.25_dp*rmax
         lri_section => section_vals_get_subs_vals(input, "DFT%KG_METHOD%LRIGPW")
         CALL lri_env_init(kg_env%lri_env, lri_section)
         CALL lri_env_basis("LRI", qs_env, kg_env%lri_env, qs_kind_set)
         !
         ! if energy correction is performed with force calculation,
         ! then response calculation requires
         ! perturbation density for kernel calculation
         CALL lri_env_init(kg_env%lri_env1, lri_section)
         CALL lri_env_basis("LRI", qs_env, kg_env%lri_env1, qs_kind_set)
         !
         ! integration grid
         !
         CALL section_vals_val_get(input, "DFT%KG_METHOD%INTEGRATION_GRID", c_val=intgrid)
         CALL uppercase(intgrid)
         SELECT CASE (intgrid)
         CASE ("SMALL")
            nr = 50
            na = 38
         CASE ("MEDIUM")
            nr = 100
            na = 110
         CASE ("LARGE")
            nr = 200
            na = 302
         CASE ("HUGE")
            nr = 400
            na = 590
         CASE DEFAULT
            CPABORT("KG:INTEGRATION_GRID")
         END SELECT
         NULLIFY (logger)
         logger => cp_get_default_logger()
         iunit = cp_logger_get_default_io_unit(logger)
         CALL initialize_atomic_grid(kg_env%int_grid_atom, nr, na, rmax, iunit=iunit)
         ! load balancing
         CALL get_qs_env(qs_env=qs_env, natom=natom, para_env=para_env)
         np = para_env%num_pe
         load = REAL(natom, KIND=dp)*kg_env%int_grid_atom%ntot/REAL(np, KIND=dp)
         !
         CALL allocate_intgrid(kg_env%int_grid_full)
         ig_full => kg_env%int_grid_full
         CALL allocate_intgrid(kg_env%int_grid_molecules)
         ig_mol => kg_env%int_grid_molecules
         nbatch = (natom*kg_env%int_grid_atom%nbatch)/np
         nbatch = NINT((nbatch + 1)*1.2_dp)
         ALLOCATE (bid(2, nbatch))
         nbatch = 0
         DO iatom = 1, natom
            DO ib = 1, kg_env%int_grid_atom%nbatch
               IF (para_env%mepos == MOD(iatom + ib, np)) THEN
                  nbatch = nbatch + 1
                  CPASSERT(nbatch <= SIZE(bid, 2))
                  bid(1, nbatch) = iatom
                  bid(2, nbatch) = ib
               END IF
            END DO
         END DO
         !
         ig_full%nbatch = nbatch
         ALLOCATE (ig_full%grid_batch(nbatch))
         !
         ig_mol%nbatch = nbatch
         ALLOCATE (ig_mol%grid_batch(nbatch))
         !
         DO i = 1, nbatch
            iatom = bid(1, i)
            ib = bid(2, i)
            !
            ig_full%grid_batch(i)%ref_atom = iatom
            ig_full%grid_batch(i)%ibatch = ib
            ig_full%grid_batch(i)%np = kg_env%int_grid_atom%batch(ib)%np
            ig_full%grid_batch(i)%radius = kg_env%int_grid_atom%batch(ib)%rad
            ig_full%grid_batch(i)%rcenter(1:3) = kg_env%int_grid_atom%batch(ib)%rcenter(1:3)
            n = ig_full%grid_batch(i)%np
            ALLOCATE (ig_full%grid_batch(i)%rco(3, n))
            ALLOCATE (ig_full%grid_batch(i)%weight(n))
            ALLOCATE (ig_full%grid_batch(i)%wref(n))
            ALLOCATE (ig_full%grid_batch(i)%wsum(n))
            ig_full%grid_batch(i)%weight(:) = kg_env%int_grid_atom%batch(ib)%weight(:)
            !
            ig_mol%grid_batch(i)%ref_atom = iatom
            ig_mol%grid_batch(i)%ibatch = ib
            ig_mol%grid_batch(i)%np = kg_env%int_grid_atom%batch(ib)%np
            ig_mol%grid_batch(i)%radius = kg_env%int_grid_atom%batch(ib)%rad
            ig_mol%grid_batch(i)%rcenter(1:3) = kg_env%int_grid_atom%batch(ib)%rcenter(1:3)
            n = ig_mol%grid_batch(i)%np
            ALLOCATE (ig_mol%grid_batch(i)%rco(3, n))
            ALLOCATE (ig_mol%grid_batch(i)%weight(n))
            ALLOCATE (ig_mol%grid_batch(i)%wref(n))
            ALLOCATE (ig_mol%grid_batch(i)%wsum(n))
            ig_mol%grid_batch(i)%weight(:) = kg_env%int_grid_atom%batch(ib)%weight(:)
         END DO
         !
         DEALLOCATE (bid)
      END IF

      CALL timestop(handle)

   END SUBROUTINE init_kg_env

! **************************************************************************************************
!> \brief builds either the full neighborlist or neighborlists of molecular
!> \brief subsets, depending on parameter values
!> \param qs_env ...
!> \param sab_orb the return type, a neighborlist
!> \param sac_kin ...
!> \param molecular if false, the full neighborlist is build
!> \param subset_of_mol the molecular subsets
!> \param current_subset the subset of which the neighborlist is to be build
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! **************************************************************************************************
   SUBROUTINE kg_build_neighborlist(qs_env, sab_orb, sac_kin, &
                                    molecular, subset_of_mol, current_subset)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb, sac_kin
      LOGICAL, OPTIONAL                                  :: molecular
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: subset_of_mol
      INTEGER, OPTIONAL                                  :: current_subset

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_build_neighborlist'

      INTEGER                                            :: handle, ikind, nkind
      LOGICAL                                            :: mic, molecule_only
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present, tpot_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius, tpot_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(local_potential_type), POINTER                :: tnadd_potential
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: neighbor_list_section

      CALL timeset(routineN, handle)
      NULLIFY (para_env)

      ! restrict lists to molecular subgroups
      molecule_only = .FALSE.
      IF (PRESENT(molecular)) molecule_only = molecular
      ! enforce minimum image convention if we use molecules
      mic = molecule_only

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      distribution_2d=distribution_2d, &
                      molecule_set=molecule_set, &
                      local_particles=distribution_1d, &
                      particle_set=particle_set, &
                      para_env=para_env)

      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ! Allocate work storage
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (orb_radius(nkind), tpot_radius(nkind))
      orb_radius(:) = 0.0_dp
      tpot_radius(:) = 0.0_dp
      ALLOCATE (orb_present(nkind), tpot_present(nkind))
      ALLOCATE (pair_radius(nkind, nkind))
      ALLOCATE (atom2d(nkind))

      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, molecule_only, particle_set=particle_set)

      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom2d(ikind)%list)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            orb_present(ikind) = .TRUE.
            IF (PRESENT(subset_of_mol)) THEN
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, kind_radius=orb_radius(ikind))
            ELSE
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, short_kind_radius=orb_radius(ikind))
            END IF
         ELSE
            orb_present(ikind) = .FALSE.
            orb_radius(ikind) = 0.0_dp
         END IF
      END DO

      IF (PRESENT(sab_orb)) THEN
         ! Build the orbital-orbital overlap neighbor list
         CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)
         IF (PRESENT(subset_of_mol)) THEN
            CALL build_neighbor_lists(sab_orb, particle_set, atom2d, cell, pair_radius, &
                                      mic=mic, subcells=subcells, molecular=molecule_only, subset_of_mol=subset_of_mol, &
                                      current_subset=current_subset, nlname="sab_orb")
         ELSE
            CALL build_neighbor_lists(sab_orb, particle_set, atom2d, cell, pair_radius, &
                                      mic=mic, subcells=subcells, molecular=molecule_only, nlname="sab_orb")
         END IF
         ! Print out the neighborlist
         neighbor_list_section => section_vals_get_subs_vals(qs_env%input, "DFT%KG_METHOD%PRINT%NEIGHBOR_LISTS")
         IF (molecule_only) THEN
            CALL write_neighbor_lists(sab_orb, particle_set, cell, para_env, neighbor_list_section, &
                                      "/SAB_ORB_MOLECULAR", "sab_orb", "MOLECULAR SUBSET NEIGHBORLIST")
         ELSE
            CALL write_neighbor_lists(sab_orb, particle_set, cell, para_env, neighbor_list_section, &
                                      "/SAB_ORB_FULL", "sab_orb", "FULL NEIGHBORLIST")
         END IF
      END IF

      IF (PRESENT(sac_kin)) THEN
         DO ikind = 1, nkind
            tpot_present(ikind) = .FALSE.
            CALL get_qs_kind(qs_kind_set(ikind), tnadd_potential=tnadd_potential)
            IF (ASSOCIATED(tnadd_potential)) THEN
               CALL get_potential(potential=tnadd_potential, radius=tpot_radius(ikind))
               tpot_present(ikind) = .TRUE.
            END IF
         END DO
         CALL pair_radius_setup(orb_present, tpot_present, orb_radius, tpot_radius, pair_radius)
         CALL build_neighbor_lists(sac_kin, particle_set, atom2d, cell, pair_radius, &
                                   subcells=subcells, operator_type="ABC", nlname="sac_kin")
         neighbor_list_section => section_vals_get_subs_vals(qs_env%input, &
                                                             "DFT%KG_METHOD%PRINT%NEIGHBOR_LISTS")
         CALL write_neighbor_lists(sac_kin, particle_set, cell, para_env, neighbor_list_section, &
                                   "/SAC_KIN", "sac_kin", "ORBITAL kin energy potential")
      END IF

      ! Release work storage
      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d)
      DEALLOCATE (orb_present, tpot_present)
      DEALLOCATE (orb_radius, tpot_radius)
      DEALLOCATE (pair_radius)

      CALL timestop(handle)

   END SUBROUTINE kg_build_neighborlist

! **************************************************************************************************
!> \brief Removes all replicated pairs from a 2d integer buffer array
!> \param pairs_buffer the array, assumed to have the shape (2,:)
!> \param n number of pairs (in), number of disjunct pairs (out)
!> \par History
!>       2012.07 created [Martin Haeufel]
!>       2014.11 simplified [Ole Schuett]
!> \author Martin Haeufel
! **************************************************************************************************
   SUBROUTINE kg_remove_duplicates(pairs_buffer, n)
      INTEGER(KIND=int_4), DIMENSION(:, :), &
         INTENT(INOUT)                                   :: pairs_buffer
      INTEGER, INTENT(INOUT)                             :: n

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_remove_duplicates'

      INTEGER                                            :: handle, i, npairs
      INTEGER, DIMENSION(n)                              :: ind
      INTEGER(KIND=int_8), DIMENSION(n)                  :: sort_keys
      INTEGER(KIND=int_4), DIMENSION(2, n)               :: work

      CALL timeset(routineN, handle)

      IF (n > 0) THEN
         ! represent a pair of int_4 as a single int_8 number, simplifies sorting.
         sort_keys(1:n) = ISHFT(INT(pairs_buffer(1, 1:n), KIND=int_8), 8*int_4_size)
         sort_keys(1:n) = sort_keys(1:n) + pairs_buffer(2, 1:n) !upper + lower bytes
         CALL sort(sort_keys, n, ind)

         ! add first pair, the case npairs==0 was excluded above
         npairs = 1
         work(:, 1) = pairs_buffer(:, ind(1))

         ! remove duplicates from the sorted list
         DO i = 2, n
            IF (sort_keys(i) /= sort_keys(i - 1)) THEN
               npairs = npairs + 1
               work(:, npairs) = pairs_buffer(:, ind(i))
            END IF
         END DO

         n = npairs
         pairs_buffer(:, :n) = work(:, :n)
      END IF

      CALL timestop(handle)

   END SUBROUTINE kg_remove_duplicates

! **************************************************************************************************
!> \brief writes the graph to file using the DIMACS standard format
!>        for a definition of the file format see
!>        mat.gsia.cmu.edu?COLOR/general/ccformat.ps
!>        c comment line
!>        p edge NODES EDGES
!>        with NODES - number of nodes
!>        EDGES - numer of edges
!>        e W V
!>        ...
!>        there is one edge descriptor line for each edge in the graph
!>        for an edge (w,v) the fields W and V specify its endpoints
!> \param pairs ...
!> \param nnodes the total number of nodes
! **************************************************************************************************
   SUBROUTINE write_to_file(pairs, nnodes)
      INTEGER(KIND=int_4), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(IN)                     :: pairs
      INTEGER, INTENT(IN)                                :: nnodes

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_to_file'

      INTEGER                                            :: handle, i, imol, jmol, npairs, unit_nr
      INTEGER(KIND=int_4), ALLOCATABLE, DIMENSION(:, :)  :: sorted_pairs

      CALL timeset(routineN, handle)

      ! get the number of disjunct pairs
      npairs = SIZE(pairs, 2)

      ALLOCATE (sorted_pairs(2, npairs))

      ! reorder pairs such that pairs(1,*) < pairs(2,*)
      DO i = 1, npairs
         ! get molecular ids
         imol = pairs(1, i)
         jmol = pairs(2, i)
         IF (imol > jmol) THEN
            ! switch pair and store
            sorted_pairs(1, i) = jmol
            sorted_pairs(2, i) = imol
         ELSE
            ! keep ordering just copy
            sorted_pairs(1, i) = imol
            sorted_pairs(2, i) = jmol
         END IF
      END DO

      ! remove duplicates and get the number of disjunct pairs (number of edges)
      CALL kg_remove_duplicates(sorted_pairs, npairs)

      ! should now be half as much pairs as before
      CPASSERT(npairs == SIZE(pairs, 2)/2)

      CALL open_file(unit_number=unit_nr, file_name="graph.col")

      WRITE (unit_nr, '(A6,1X,I8,1X,I8)') "p edge", nnodes, npairs

      ! only write out the first npairs entries
      DO i = 1, npairs
         WRITE (unit_nr, '(A1,1X,I8,1X,I8)') "e", sorted_pairs(1, i), sorted_pairs(2, i)
      END DO

      CALL close_file(unit_nr)

      DEALLOCATE (sorted_pairs)

      CALL timestop(handle)

   END SUBROUTINE write_to_file

! **************************************************************************************************
!> \brief ...
!> \param kg_env ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE kg_build_subsets(kg_env, para_env)
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kg_build_subsets'

      INTEGER                                            :: color, handle, i, iatom, imol, isub, &
                                                            jatom, jmol, nmol, npairs, npairs_local
      INTEGER(KIND=int_4)                                :: ncolors
      INTEGER(KIND=int_4), ALLOCATABLE, DIMENSION(:)     :: color_of_node
      INTEGER(KIND=int_4), ALLOCATABLE, DIMENSION(:, :)  :: msg_gather, pairs, pairs_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nnodes_of_color
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)

      ! first: get a (local) list of pairs from the (local) neighbor list data
      nmol = SIZE(kg_env%molecule_set)

      npairs = 0
      CALL neighbor_list_iterator_create(nl_iterator, kg_env%sab_orb_full)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom)

         imol = kg_env%atom_to_molecule(iatom)
         jmol = kg_env%atom_to_molecule(jatom)

         !IF (imol<jmol) THEN
         IF (imol .NE. jmol) THEN

            npairs = npairs + 2

         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ALLOCATE (pairs_buffer(2, npairs))

      npairs = 0
      CALL neighbor_list_iterator_create(nl_iterator, kg_env%sab_orb_full)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom)

         imol = kg_env%atom_to_molecule(iatom)
         jmol = kg_env%atom_to_molecule(jatom)

         IF (imol .NE. jmol) THEN

            ! add pair to the local list

            ! add both orderings - makes it easier to build the neighborlist
            npairs = npairs + 1

            pairs_buffer(1, npairs) = imol
            pairs_buffer(2, npairs) = jmol

            npairs = npairs + 1

            pairs_buffer(2, npairs) = imol
            pairs_buffer(1, npairs) = jmol

         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! remove duplicates
      CALL kg_remove_duplicates(pairs_buffer, npairs)

      ! get the maximum number of local pairs on all nodes (size of the mssg)
      ! remember how many pairs we have local
      npairs_local = npairs
      CALL mp_max(npairs, para_env%group)

      ! allocate message
      ALLOCATE (pairs(2, npairs))

      pairs(:, 1:npairs_local) = pairs_buffer(:, 1:npairs_local)
      pairs(:, npairs_local + 1:) = 0

      DEALLOCATE (pairs_buffer)

      ! second: gather all data on the master node
      ! better would be a scheme where duplicates are removed in a tree-like reduction scheme.
      ! this step can be needlessly memory intensive in the current implementation.

      IF (para_env%source .EQ. para_env%mepos) THEN
         ALLOCATE (msg_gather(2, npairs*para_env%num_pe))
      ELSE
         ALLOCATE (msg_gather(2, 1))
      ENDIF

      msg_gather = 0

      CALL mp_gather(pairs, msg_gather, para_env%source, para_env%group)

      DEALLOCATE (pairs)

      IF (para_env%source .EQ. para_env%mepos) THEN

         ! shift all non-zero entries to the beginning of the array and count the number of actual pairs
         npairs = 0

         DO i = 1, SIZE(msg_gather, 2)
            IF (msg_gather(1, i) .NE. 0) THEN
               npairs = npairs + 1
               msg_gather(:, npairs) = msg_gather(:, i)
            END IF
         END DO

         ! remove duplicates
         CALL kg_remove_duplicates(msg_gather, npairs)

         ALLOCATE (pairs(2, npairs))

         pairs(:, 1:npairs) = msg_gather(:, 1:npairs)

         DEALLOCATE (msg_gather)

         !WRITE(*,'(A48,5X,I10,4X,A2,1X,I10)') " KG| Total number of overlapping molecular pairs",npairs/2,"of",nmol*(nmol-1)/2

         ! write to file, nnodes = number of molecules
         IF (.FALSE.) THEN
            CALL write_to_file(pairs, SIZE(kg_env%molecule_set))
         ENDIF

         ! vertex coloring algorithm
         CALL kg_vertex_coloring(kg_env, pairs, ncolors, color_of_node)

         DEALLOCATE (pairs)

      ELSE

         DEALLOCATE (msg_gather)

      END IF

      !WRITE(*,'(A27,40X,I6,1X,A6)') ' KG| Vertex coloring result', ncolors, 'colors'

      ! broadcast the number of colors to all nodes
      CALL mp_bcast(ncolors, para_env%source, para_env%group)

      IF (.NOT. ALLOCATED(color_of_node)) ALLOCATE (color_of_node(nmol))

      ! broadcast the resulting coloring to all nodes.....
      CALL mp_bcast(color_of_node, para_env%source, para_env%group)

      IF ((kg_env%nsubsets .NE. 0) .AND. (ncolors .NE. kg_env%nsubsets)) THEN
         ! number of subsets has changed

         ! deallocate stuff if necessary
         IF (ASSOCIATED(kg_env%subset)) THEN
            DO isub = 1, kg_env%nsubsets
               CALL release_neighbor_list_sets(kg_env%subset(isub)%sab_orb)
               CALL deallocate_task_list(kg_env%subset(isub)%task_list)
            END DO
            DEALLOCATE (kg_env%subset)
            NULLIFY (kg_env%subset)
         END IF

      END IF

      ! allocate and nullify some stuff
      IF (.NOT. ASSOCIATED(kg_env%subset)) THEN

         ALLOCATE (kg_env%subset(ncolors))

         DO i = 1, ncolors
            NULLIFY (kg_env%subset(i)%sab_orb)
            NULLIFY (kg_env%subset(i)%task_list)
         END DO
      END IF

      ! set the number of subsets
      kg_env%nsubsets = ncolors

      ! counting loop
      ALLOCATE (nnodes_of_color(ncolors))
      nnodes_of_color = 0
      DO i = 1, nmol ! nmol=nnodes
         color = color_of_node(i)
         kg_env%subset_of_mol(i) = color
         nnodes_of_color(color) = nnodes_of_color(color) + 1
      END DO

      DEALLOCATE (nnodes_of_color)
      DEALLOCATE (color_of_node)

      CALL timestop(handle)

   END SUBROUTINE kg_build_subsets

END MODULE kg_environment
