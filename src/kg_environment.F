!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a Kim-Gordon-like partitioning into molecular subunits 
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************
MODULE kg_environment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE external_potential_types,        ONLY: get_potential,&
                                             local_potential_type
  USE f77_blas
  USE input_constants,                 ONLY: kg_tnadd_atomic,&
                                             kg_tnadd_embed,&
                                             xc_funct_no_shortcut
  USE input_section_types,             ONLY: &
       section_vals_create, section_vals_duplicate, section_vals_get, &
       section_vals_get_subs_vals, section_vals_release, &
       section_vals_set_subs_vals, section_vals_type, section_vals_val_get, &
       section_vals_val_set
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kg_vertex_coloring_methods,      ONLY: kg_vertex_coloring
  USE kinds,                           ONLY: dp,&
                                             int_4,&
                                             int_8
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_gather,&
                                             mp_max
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                             atom2d_cleanup,&
                                             build_neighbor_lists,&
                                             local_atoms_type,&
                                             pair_radius_setup,&
                                             write_neighbor_lists
  USE task_list_types,                 ONLY: deallocate_task_list
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_environment'

  PUBLIC :: kg_env_create, kg_build_neighborlist, kg_build_subsets

CONTAINS

! *****************************************************************************
!> \brief Allocates and intitializes kg_env
!> \param kg_env the object to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************
  SUBROUTINE kg_env_create(kg_env, input, error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(kg_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL init_kg_env(kg_env, input, error=error)
    END IF
  END SUBROUTINE kg_env_create

! *****************************************************************************
!> \brief Initializes kg_env
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************
  SUBROUTINE init_kg_env(kg_env, input, error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_kg_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ke_sections
    LOGICAL                                  :: failure, is_set
    TYPE(section_vals_type), POINTER         :: kef_section, xc_fun_section, &
                                                xc_fun_section_kg, &
                                                xc_section, xc_section_kg

    failure=.FALSE.

    CALL timeset(routineN, handle)

    NULLIFY(kg_env%sab_orb_full)
    NULLIFY(kg_env%sac_kin)
    NULLIFY(kg_env%subset_of_mol)
    NULLIFY(kg_env%subset)
    NULLIFY(kg_env%tnadd_mat)

    kg_env%nsubsets=0

    ! get coloring method settings
    CALL section_vals_val_get(input,"DFT%KG_METHOD%COLORING_METHOD",i_val=kg_env%coloring_method,error=error)
    ! get method for nonadditive kinetic energy embedding potential
    CALL section_vals_val_get(input,"DFT%KG_METHOD%TNADD_METHOD",i_val=kg_env%tnadd_method,error=error)

    IF (kg_env%tnadd_method == kg_tnadd_embed) THEN
       ! generate a new XC section with only the KE functional
       NULLIFY(xc_fun_section_kg, xc_section_kg)

       xc_section => section_vals_get_subs_vals(input, "DFT%XC", error=error)
       xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL",error=error)

       CALL section_vals_duplicate(xc_section, xc_section_kg, error=error)

       ! test for TF, TFW, KE_GGA, and KE_LIBXC sections in turn
       ke_sections = 0
       NULLIFY(kef_section)
       kef_section => section_vals_get_subs_vals(xc_fun_section,"TF",error=error)
       CALL section_vals_get(kef_section, explicit=is_set, error=error)
       IF(is_set) THEN
         NULLIFY(xc_fun_section_kg)
         CALL section_vals_create(xc_fun_section_kg,xc_fun_section%section,&
              error=error)
         CALL section_vals_val_set(xc_fun_section_kg,"_SECTION_PARAMETERS_",&
              i_val=xc_funct_no_shortcut,error=error)
         CALL section_vals_set_subs_vals(xc_fun_section_kg,"TF",&
              kef_section,error=error)
         CALL section_vals_set_subs_vals(xc_section_kg,"XC_FUNCTIONAL",&
              xc_fun_section_kg,error=error)
         ke_sections = ke_sections + 1
       END IF

       NULLIFY(kef_section)
       kef_section => section_vals_get_subs_vals(xc_fun_section,"TFW",error=error)
       CALL section_vals_get(kef_section, explicit=is_set, error=error)
       IF(is_set) THEN
         NULLIFY(xc_fun_section_kg)
         CALL section_vals_create(xc_fun_section_kg,xc_fun_section%section,&
              error=error)
         CALL section_vals_val_set(xc_fun_section_kg,"_SECTION_PARAMETERS_",&
              i_val=xc_funct_no_shortcut,error=error)
         CALL section_vals_set_subs_vals(xc_fun_section_kg,"TFW",&
              kef_section,error=error)
         CALL section_vals_set_subs_vals(xc_section_kg,"XC_FUNCTIONAL",&
              xc_fun_section_kg,error=error)
         ke_sections = ke_sections + 1
       END IF

       NULLIFY(kef_section)
       kef_section => section_vals_get_subs_vals(xc_fun_section,"KE_GGA",error=error)
       CALL section_vals_get(kef_section, explicit=is_set, error=error)
       IF(is_set) THEN
         NULLIFY(xc_fun_section_kg)
         CALL section_vals_create(xc_fun_section_kg,xc_fun_section%section,&
              error=error)
         CALL section_vals_val_set(xc_fun_section_kg,"_SECTION_PARAMETERS_",&
              i_val=xc_funct_no_shortcut,error=error)
         CALL section_vals_set_subs_vals(xc_fun_section_kg,"KE_GGA",&
              kef_section,error=error)
         CALL section_vals_set_subs_vals(xc_section_kg,"XC_FUNCTIONAL",&
              xc_fun_section_kg,error=error)
         ke_sections = ke_sections + 1
       END IF

       NULLIFY(kef_section)
       kef_section => section_vals_get_subs_vals(xc_fun_section,"KE_LIBXC",error=error)
       CALL section_vals_get(kef_section, explicit=is_set, error=error)
       IF(is_set) THEN
         NULLIFY(xc_fun_section_kg)
         CALL section_vals_create(xc_fun_section_kg,xc_fun_section%section,&
              error=error)
         CALL section_vals_val_set(xc_fun_section_kg,"_SECTION_PARAMETERS_",&
              i_val=xc_funct_no_shortcut,error=error)
         CALL section_vals_set_subs_vals(xc_fun_section_kg,"KE_LIBXC",&
              kef_section,error=error)
         CALL section_vals_set_subs_vals(xc_section_kg,"XC_FUNCTIONAL",&
              xc_fun_section_kg,error=error)
         ke_sections = ke_sections + 1
       END IF

       ! stop if there is more than one kinetic energy functional section present
       !CPPostcondition(ke_sections==1,cp_failure_level,routineP,error,failure)
       CALL cp_assert(ke_sections==1, cp_fatal_level, cp_assertion_failed, routineP, &
          "KG runs require a kinetic energy functional. Exactly one of the following "// &
          "subsections of XC_FUNCTIONAL needs to be present: KE_GGA, TF, TFW, or KE_LIBXC.") 

       kg_env%xc_section_kg => xc_section_kg

       NULLIFY(kef_section, xc_fun_section)

       CALL section_vals_release(xc_fun_section_kg, error=error)
    
    ELSEIF(kg_env%tnadd_method == kg_tnadd_atomic) THEN
       NULLIFY(kg_env%xc_section_kg)
    ELSE
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE init_kg_env

! *****************************************************************************
!> \brief builds either the full neighborlist or neighborlists of molecular
!> \brief subsets, depending on parameter values
!> \param sab_orb: the return type, a neighborlist
!> \param molecular: if false, the full neighborlist is build
!> \param subset: the molecular subsets
!> \param current_subset: the subset of which the neighborlist is to be build
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************
  SUBROUTINE kg_build_neighborlist(kg_env, qs_env, sab_orb, sac_kin,&
             molecular, subset_of_mol, current_subset, error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_orb, sac_kin
    LOGICAL, OPTIONAL                        :: molecular
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: subset_of_mol
    INTEGER, OPTIONAL                        :: current_subset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_build_neighborlist', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, nkind, stat
    LOGICAL                                  :: failure, mic, molecule_only
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: orb_present, tpot_present
    REAL(dp)                                 :: subcells
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: orb_radius, tpot_radius
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: pair_radius
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom2d
    TYPE(local_potential_type), POINTER      :: tnadd_potential
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: neighbor_list_section

    CALL timeset(routineN, handle)

    failure = .FALSE.

    molecule_only = .FALSE.
    IF (PRESENT(molecular)) molecule_only = molecular

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    distribution_2d=distribution_2d,&
                    molecule_set=molecule_set,&
                    local_particles=distribution_1d,&
                    particle_set=particle_set,&
                    error=error)

    ! Allocate work storage
    nkind = SIZE(atomic_kind_set)
    ALLOCATE (orb_radius(nkind),tpot_radius(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    orb_radius(:) = 0.0_dp
    tpot_radius(:) = 0.0_dp
    ALLOCATE (orb_present(nkind),tpot_present(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (pair_radius(nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (atom2d(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL section_vals_val_get(qs_env%input,"DFT%SUBCELLS",r_val=subcells,error=error)

    IF (PRESENT(subset_of_mol)) THEN
      CALL atom2d_build(atom2d,orb_radius,orb_present,distribution_1d,distribution_2d,&
                     atomic_kind_set,molecule_set,molecule_only,kg=.FALSE.,&
                     dftb=.FALSE.,particle_set=particle_set,error=error)
    ELSE
      CALL atom2d_build(atom2d,orb_radius,orb_present,distribution_1d,distribution_2d,&
                     atomic_kind_set,molecule_set,molecule_only,kg=.TRUE.,&
                     dftb=.FALSE.,particle_set=particle_set,error=error)
    END IF

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind, atom_list=atom2d(ikind)%list)
    END DO

    IF (PRESENT(sab_orb)) THEN

       ! Build the orbital-orbital overlap neighbor list
       CALL pair_radius_setup(orb_present,orb_present,orb_radius,orb_radius,pair_radius,error)

       IF (PRESENT(subset_of_mol)) THEN
         CALL build_neighbor_lists(sab_orb,particle_set,atom2d,cell,pair_radius,&
                                mic=.FALSE.,subcells=subcells,molecular=molecule_only,subset_of_mol=subset_of_mol,&
                                current_subset=current_subset,name="sab_orb",error=error)
       ELSE
         CALL build_neighbor_lists(sab_orb,particle_set,atom2d,cell,pair_radius,&
                                mic=.FALSE.,subcells=subcells,molecular=molecule_only,name="sab_orb",error=error)
       END IF

       ! Print out the neighborlist
       neighbor_list_section => section_vals_get_subs_vals(qs_env%input,"DFT%KG_METHOD%PRINT%NEIGHBOR_LISTS",error=error) 

       IF (molecule_only) THEN
         CALL write_neighbor_lists(sab_orb,particle_set,cell,qs_env%para_env,neighbor_list_section,&
              "/SAB_ORB_MOLECULAR","sab_orb","MOLECULAR SUBSET NEIGHBORLIST",error)
       ELSE
         CALL write_neighbor_lists(sab_orb,particle_set,cell,qs_env%para_env,neighbor_list_section,&
              "/SAB_ORB_FULL","sab_orb","FULL NEIGHBORLIST",error)
       END IF

    END IF

    IF (PRESENT(sac_kin)) THEN
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          tpot_present(ikind) = .FALSE.
          CALL get_atomic_kind(atomic_kind=atomic_kind,tnadd_potential=tnadd_potential)
          IF (ASSOCIATED(tnadd_potential)) THEN
             CALL get_potential(potential=tnadd_potential,radius=tpot_radius(ikind))
             tpot_present(ikind) = .TRUE.
          END IF
       END DO
       CALL pair_radius_setup(orb_present,tpot_present,orb_radius,tpot_radius,pair_radius,error)
       CALL build_neighbor_lists(sac_kin,particle_set,atom2d,cell,pair_radius,&
            subcells=subcells,operator_type="ABC",name="sac_kin",error=error)
       neighbor_list_section => section_vals_get_subs_vals(qs_env%input,&
                                         "DFT%KG_METHOD%PRINT%NEIGHBOR_LISTS",error=error) 
       CALL write_neighbor_lists(sac_kin,particle_set,cell,qs_env%para_env,neighbor_list_section,&
              "/SAC_KIN","sac_kin","ORBITAL kin energy potential",error)
    END IF

    ! Release work storage
    CALL atom2d_cleanup(atom2d,error)
    DEALLOCATE (atom2d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (orb_present,tpot_present,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (orb_radius,tpot_radius,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (pair_radius,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE kg_build_neighborlist

! *****************************************************************************
!> \brief Removes all replicated pairs from a 2d integer buffer array
!> \param pairs_buffer: the array, assumed to have the shape (2,:)
!> \param n: number of pairs (in), number of disjunct pairs (out)
!> \par History
!>       2012.07 created [Martin Haeufel]
!> \author Martin Haeufel
! *****************************************************************************
  SUBROUTINE kg_remove_duplicates(pairs_buffer, n)
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:, :), INTENT(INOUT)         :: pairs_buffer
    INTEGER, INTENT(INOUT)                   :: n

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_remove_duplicates', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, npairs
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: work
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ind

    CALL timeset(routineN, handle)

    IF (n.GT.0) THEN

      ! allocate work storage
      ALLOCATE(work(n))
  
      ! copy the 2d int_4 array into a 1d int_8 array
      work(1:n) = INT(pairs_buffer(1,1:n), KIND=int_8)*HUGE(pairs_buffer(1,1)) &
                                       + pairs_buffer(2,1:n)
      ! sort
      ALLOCATE(ind(n))
      CALL sort(work, n, ind)
      DEALLOCATE(ind)

      npairs=1 ! n is always greater than 0
  
      ! remove duplicates from the sorted list
      DO i=2, n
        IF (work(i).NE.work(i-1)) THEN
          npairs=npairs+1
          work(npairs)=work(i)
        END IF
      END DO
  
      ! copy back to a 2d int_4 array
      pairs_buffer(1,1:npairs) = work(1:npairs)/HUGE(pairs_buffer(1,1))
      pairs_buffer(2,1:npairs) = work(1:npairs) &
                                         - INT(pairs_buffer(1,1:npairs), KIND=int_8)*HUGE(pairs_buffer(1,1)) 
  
      ! deallocate work memory
      DEALLOCATE(work)
  
      n = npairs

    END IF

    CALL timestop(handle)

  END SUBROUTINE

  ! writes the graph to file using the DIMACS standard format
  ! for a definition of the file format see 
  ! mat.gsia.cmu.edu?COLOR/general/ccformat.ps
  ! in brief:
  ! c comment line
  ! p edge NODES EDGES
  ! with NODES - number of nodes
  !      EDGES - numer of edges
  ! e W V
  ! ...
  ! there is one edge descriptor line for each edge in the graph
  ! for an edge (w,v) the fields W and V specify its endpoints
  SUBROUTINE write_to_file(pairs, nnodes, error)
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:, :), INTENT(IN)            :: pairs
    INTEGER, INTENT(IN)                      :: nnodes
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_to_file', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, imol, jmol, npairs
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sorted_pairs
    LOGICAL                                  :: failure

    CALL timeset(routineN, handle)

    failure=.FALSE.

    ! get the number of disjunct pairs
    npairs=SIZE(pairs,2)

    ALLOCATE(sorted_pairs(2,npairs))

    ! reorder pairs such that pairs(1,*) < pairs(2,*) 
    DO i=1,npairs
      ! get molecular ids
      imol = pairs(1,i)
      jmol = pairs(2,i)
      IF (imol>jmol) THEN
        ! switch pair and store
        sorted_pairs(1,i)=jmol
        sorted_pairs(2,i)=imol
      ELSE
        ! keep ordering just copy
        sorted_pairs(1,i)=imol
        sorted_pairs(2,i)=jmol
      END IF
    END DO

    ! remove duplicates and get the number of disjunct pairs (number of edges)
    CALL kg_remove_duplicates(sorted_pairs, npairs)

    ! should now be half as much pairs as before
    CPPostcondition(npairs==SIZE(pairs,2)/2,cp_failure_level,routineP,error,failure)

    OPEN(33, FILE="graph.col")

    WRITE(33, '(A6,1X,I8,1X,I8)') "p edge", nnodes, npairs

    ! only write out the first npairs entries
    DO i=1,npairs
      WRITE(33, '(A1,1X,I8,1X,I8)') "e", sorted_pairs(1,i), sorted_pairs(2,i)
    END DO

    CLOSE(33)

    DEALLOCATE(sorted_pairs)

    CALL timestop(handle)

  END SUBROUTINE

  SUBROUTINE kg_build_subsets(kg_env, para_env, error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_build_subsets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: color, handle, i, iab, iatom, &
                                                imol, isub, jatom, jmol, &
                                                nmol, npairs, npairs_local
    INTEGER(KIND=int_4)                      :: ncolors
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:)                           :: color_of_node
    INTEGER(KIND=int_4), ALLOCATABLE, &
      DIMENSION(:, :)                        :: msg_gather, pairs, &
                                                pairs_buffer
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index_of_color, &
                                                nnodes_of_color
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

! first: get a (local) list of pairs from the (local) neighbor list data

    CALL timeset(routineN, handle)

    nmol = SIZE(kg_env%molecule_set)


    npairs=0
    CALL neighbor_list_iterator_create(nl_iterator, kg_env%sab_orb_full) 
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
      CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom)
                 
      imol = kg_env%atom_to_molecule(iatom)
      jmol = kg_env%atom_to_molecule(jatom)

      !IF (imol<jmol) THEN
      IF (imol.NE.jmol) THEN

        npairs=npairs+2

      END IF

    END DO
    CALL neighbor_list_iterator_release(nl_iterator) 

    ALLOCATE(pairs_buffer(2,npairs))

    npairs=0
    CALL neighbor_list_iterator_create(nl_iterator, kg_env%sab_orb_full) 
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
      CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom)

      imol = kg_env%atom_to_molecule(iatom)
      jmol = kg_env%atom_to_molecule(jatom)

      IF (imol.NE.jmol) THEN

        ! add pair to the local list

        ! add both orderings - makes it easier to build the neighborlist
        npairs=npairs+1

        pairs_buffer(1,npairs) = imol
        pairs_buffer(2,npairs) = jmol

        npairs=npairs+1

        pairs_buffer(2,npairs) = imol
        pairs_buffer(1,npairs) = jmol

      END IF

    END DO
    CALL neighbor_list_iterator_release(nl_iterator) 

    ! remove duplicates
    CALL kg_remove_duplicates(pairs_buffer, npairs)

    ! get the maximum number of local pairs on all nodes (size of the mssg)
    ! remember how many pairs we have local
    npairs_local=npairs
    CALL mp_max(npairs, para_env%group)

    ! allocate message
    ALLOCATE(pairs(2,npairs))

    pairs(:,1:npairs_local)  = pairs_buffer(:,1:npairs_local)
    pairs(:,npairs_local+1:) = 0

    DEALLOCATE(pairs_buffer)

    ! second: gather all data on the master node

    IF (para_env%source.EQ.para_env%mepos) THEN
       ALLOCATE(msg_gather(2,npairs*para_env%num_pe))
    ELSE
       ALLOCATE(msg_gather(2,1))
    ENDIF

    msg_gather=0

    CALL mp_gather(pairs, msg_gather, para_env%source, para_env%group) 

    DEALLOCATE(pairs)

    IF (para_env%source.EQ.para_env%mepos) THEN

      ! shift all non-zero entries to the beginning of the array and count the number of actual pairs
      npairs=0
  
      DO i=1, SIZE(msg_gather,2)
        IF(msg_gather(1,i).NE.0) THEN
          npairs = npairs+1
          msg_gather(:,npairs) = msg_gather(:,i) 
        END IF
      END DO
  
      ! remove duplicates
      CALL kg_remove_duplicates(msg_gather, npairs)    
 
      ALLOCATE(pairs(2,npairs))

      pairs(:,1:npairs) = msg_gather(:,1:npairs)

      DEALLOCATE(msg_gather)

      !WRITE(*,'(A48,5X,I10,4X,A2,1X,I10)') " KG| Total number of overlapping molecular pairs",npairs/2,"of",nmol*(nmol-1)/2

      ! write to file, nnodes = number of molecules
      !CALL write_to_file(pairs, SIZE(kg_env%molecule_set), error)

      ! vertex coloring algorithm
      CALL kg_vertex_coloring(kg_env, pairs, ncolors, color_of_node, error)

      DEALLOCATE(pairs)

    ELSE

      DEALLOCATE(msg_gather)

    END IF

    !WRITE(*,'(A27,40X,I6,1X,A6)') ' KG| Vertex coloring result', ncolors, 'colors'


    ! broadcast the number of colors to all nodes
    CALL mp_bcast(ncolors, para_env%source, para_env%group)

    IF (.NOT.ALLOCATED(color_of_node)) ALLOCATE(color_of_node(nmol))

    ! broadcast the resulting coloring to all nodes.....
    CALL mp_bcast(color_of_node, para_env%source, para_env%group)

    IF ((kg_env%nsubsets.NE.0).AND.(ncolors.NE.kg_env%nsubsets)) THEN
      ! number of subsets has changed

      ! deallocate stuff if necessary
      IF (ASSOCIATED(kg_env%subset)) THEN

        DO isub=1,kg_env%nsubsets

          DO iab=1,SIZE(kg_env%subset(isub)%sab_orb)

            CALL deallocate_neighbor_list_set(kg_env%subset(isub)%sab_orb(iab)%neighbor_list_set)

          END DO

          DEALLOCATE(kg_env%subset(isub)%sab_orb)
     
          CALL deallocate_task_list(kg_env%subset(isub)%task_list, error)

        END DO

        DEALLOCATE(kg_env%subset)

        NULLIFY(kg_env%subset)
     
      END IF

    END IF

    ! allocate and nullify some stuff
    IF (.NOT.ASSOCIATED(kg_env%subset)) THEN

      ALLOCATE(kg_env%subset(ncolors))

        DO i=1,ncolors
          NULLIFY(kg_env%subset(i)%sab_orb)
          NULLIFY(kg_env%subset(i)%task_list)
        END DO
    END IF

    ! set the number of subsets
    kg_env%nsubsets = ncolors

    ! counting loop
    ALLOCATE(nnodes_of_color(ncolors))
    nnodes_of_color = 0
    DO i=1, nmol ! nmol=nnodes
      color=color_of_node(i)
      kg_env%subset_of_mol(i)=color
      nnodes_of_color(color)=nnodes_of_color(color)+1
    END DO

    DEALLOCATE(nnodes_of_color) 

    ! store the subset information
    ALLOCATE(index_of_color(ncolors))

    index_of_color=0
    DO i=1,nmol
      color=color_of_node(i)
      index_of_color(color)=index_of_color(color)+1
    END DO 
    DEALLOCATE(index_of_color)

    DEALLOCATE(color_of_node)

    CALL timestop(handle)

  END SUBROUTINE

END MODULE kg_environment
