!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_neighbor_lists [1.0] *
!!
!!   NAME
!!     qs_neighbor_lists
!!
!!   FUNCTION
!!     Generate the atomic neighbor lists.
!!
!!   AUTHOR
!!     Matthias Krack (08.10.1999,26.03.2002,16.06.2003)
!!
!!   MODIFICATION HISTORY
!!     - List rebuild for sab_orb neighbor list (10.09.2002,MK)
!!     - List rebuild for all lists (25.09.2002,MK)
!!     - Row-wise parallelized version (16.06.2003,MK)
!!     - Row- and column-wise parallelized version (19.07.2003,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_neighbor_lists

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: LOW,&
                                             MEDIUM,&
                                             global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: exclusion_list_type,&
                                             particle_type
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_local_rho_types,              ONLY: rhoz_type
  USE qs_neighbor_list_types,          ONLY: add_neighbor_list,&
                                             add_neighbor_node,&
                                             allocate_neighbor_list_set,&
                                             deallocate_neighbor_list_set,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Type definitions ***

  TYPE local_atoms_type
    INTEGER, DIMENSION(:), POINTER                   :: list,&
                                                        list_local_a_index,&
                                                        list_local_b_index,&
                                                        list_a_mol,&
                                                        list_b_mol
    REAL(dp), DIMENSION(:,:), POINTER                :: r_pbc,&
                                                        s_pbc
  END TYPE local_atoms_type

  TYPE subcell_type
    INTEGER                        :: natom
    REAL(dp), DIMENSION(3)         :: s_max,s_min
    INTEGER, DIMENSION(:), POINTER :: atom_list
  END TYPE subcell_type

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_neighbor_lists"

! *** Public data types ***

  PUBLIC :: local_atoms_type,&
            subcell_type

! *** Public subroutines ***

  PUBLIC :: allocate_subcell,&
            build_gth_ppl_neighbor_lists,&
            build_qs_neighbor_lists,&
            deallocate_subcell,&
            write_distribution,&
            write_neighbor_lists

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_qs_neighbor_lists(qs_env,globenv,molecular)

!   Purpose: Build all the required neighbor lists for Quickstep.

!   History: - Creation (28.08.2000,MK)
!            - Parallized version (16.06.2003,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, OPTIONAL                        :: molecular

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_qs_neighbor_lists"

    INTEGER :: atom_a, handle, iatom, iatom_local, ikind, istat, maxatom, &
      maxatom_local, natom_a, natom_local_a, natom_local_b, nkind, &
      output_unit, imol
    LOGICAL                                  :: all_potential_present,&
                                                gth_potential_present, &
                                                print_subcell_grid,&
                                                paw_atom_present,paw_atom,&
                                                rho0_present
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: all_present, orb_present,&
                                                oce_present, ppl_present, &
                                                ppnl_present, Qlm_present
    REAL(dp)                                 :: r_skin,rcprj
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: orb_rad_3c, orb_radius, &
                                                oce_radius, ppl_radius, &
                                                ppnl_radius, V0hmv0s_radius
    REAL(dp), DIMENSION(3)                   :: r_pbc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(paw_proj_set_type), POINTER         :: paw_proj
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_ppl, sac_ppnl, &
                                                sbc_ppnl,sac_oce,sbc_oce,&
                                                sab_2c,sac_3c,sbc_3c
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    INTEGER, DIMENSION(:), POINTER           :: local_rows_array, &
                                                local_cols_array
    LOGICAL                                  :: molecule_only
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(rhoz_type), DIMENSION(:), POINTER   :: rhoz_set
    CHARACTER(20)                            :: name

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I"," ",handle)

    NULLIFY (atomic_kind)
    NULLIFY (atomic_kind_set)
    NULLIFY (cell)
    NULLIFY (distribution_2d)
    NULLIFY (distribution_1d)
    NULLIFY (gth_potential)
    NULLIFY (orb_basis_set)
    NULLIFY (particle_set)
    NULLIFY (molecule_set)
    NULLIFY (rho0_mpole)
    NULLIFY (rhoz_set)
    NULLIFY (sab_orb)
    NULLIFY (sac_ppl)
    NULLIFY (sac_ppnl)
    NULLIFY (sbc_ppnl)
    NULLIFY (sac_oce)
    NULLIFY (sbc_oce)
    NULLIFY (sab_2c)
    NULLIFY (sac_3c)
    NULLIFY (sbc_3c)
!    NULLIFY (sac_all)
!    NULLIFY (sbc_all)

    output_unit = globenv%scr

    IF ((globenv%ionode.AND.globenv%print%cell_parameters).AND.&
        (globenv%print%level > LOW)) THEN
      print_subcell_grid = .TRUE.
    ELSE
      print_subcell_grid = .FALSE.
    END IF

    molecule_only = .FALSE.
    IF ( PRESENT ( molecular ) ) molecule_only = molecular

!   *** Set Verlet skin ***

    r_skin = 0.0_dp

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    distribution_2d=distribution_2d,&
                    local_particles=distribution_1d,&
                    particle_set=particle_set,&
                    molecule_set=molecule_set,&
                    rho0_mpole=rho0_mpole,&
                    rhoz_set=rhoz_set,&
                    sab_orb=sab_orb,&
                    sac_ppl=sac_ppl,&
                    sac_ppnl=sac_ppnl,&
                    sbc_ppnl=sbc_ppnl,&
                    sac_oce=sac_oce,&
                    sbc_oce=sbc_oce,&
                    sab_2c=sab_2c,&
                    sac_3c=sac_3c,&
                    sbc_3c=sbc_3c)
!                    sac_all=sac_all,&
!                    sbc_all=sbc_all)

    rho0_present=ASSOCIATED(rho0_mpole)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             paw_atom_present=paw_atom_present,&
                             all_potential_present=all_potential_present,&
                             maxatom=maxatom)

!   *** Allocate work storage ***

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (atom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "atom",nkind*int_size)

    ALLOCATE (orb_present(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "orb_present",nkind*int_size)

    ALLOCATE (orb_radius(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "orb_radius",nkind*dp_size)
    orb_radius(:) = 0.0_dp

    IF (gth_potential_present) THEN

      ALLOCATE (ppl_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppl_present",nkind*int_size)

      ALLOCATE (ppl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppl_radius",nkind*dp_size)
      ppl_radius = 0.0_dp

      ALLOCATE (ppnl_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppnl_present",nkind*int_size)

      ALLOCATE (ppnl_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppnl_radius",nkind*dp_size)
      ppnl_radius = 0.0_dp

    END IF
 
    IF (paw_atom_present) THEN
      ALLOCATE (oce_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "oce_present",nkind*int_size)

      ALLOCATE (oce_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "oce_radius",nkind*dp_size)
      oce_radius = 0.0_dp
    END IF

    IF (rho0_present) THEN

      IF(.NOT. ASSOCIATED(rhoz_set)) &
            CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer rhoz is not associated") 

      ALLOCATE (Qlm_present(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "Qlm_present",nkind*int_size)

      ALLOCATE (V0hmV0s_radius(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "V0hmV0s_radius",nkind*dp_size)
      V0hmV0s_radius = 0.0_dp

      ALLOCATE (orb_rad_3c(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "orb_rad_3c",nkind*dp_size)
      orb_rad_3c = 0.0_dp

    END IF

!   *** Initialize the local data structures ***

    maxatom_local = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom(ikind)%list)
      NULLIFY (atom(ikind)%list_local_a_index)
      NULLIFY (atom(ikind)%list_local_b_index)
      NULLIFY (atom(ikind)%list_a_mol)
      NULLIFY (atom(ikind)%list_b_mol)
      NULLIFY (atom(ikind)%r_pbc)
      NULLIFY (atom(ikind)%s_pbc)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom(ikind)%list,&
!                           all_potential=all_potential,&
                           gth_potential=gth_potential,&
                           orb_basis_set=orb_basis_set,&
                           paw_proj_set=paw_proj,&
                           paw_atom=paw_atom)

      natom_a = SIZE(atom(ikind)%list)
      IF ( molecule_only ) THEN
        natom_local_a = distribution_1d%n_el(ikind)
        natom_local_b = distribution_1d%n_el(ikind)
        local_rows_array => distribution_1d%list(ikind)%array
        local_cols_array => distribution_1d%list(ikind)%array
      ELSE
        natom_local_a = distribution_2d%n_local_rows(ikind)
        natom_local_b = distribution_2d%n_local_cols(ikind)
        local_rows_array => distribution_2d%local_rows(ikind)%array
        local_cols_array => distribution_2d%local_cols(ikind)%array
      END IF

      maxatom_local = MAX(maxatom_local,natom_local_a,natom_local_b)

      IF (ASSOCIATED(orb_basis_set)) THEN

        orb_present(ikind) = .TRUE.

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               kind_radius=orb_radius(ikind))

!       *** Block rows ***

        IF (natom_local_a > 0) THEN

          ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "atom(ikind)%list_local_a_index",&
                                           natom_local_a*int_size)
          ALLOCATE (atom(ikind)%list_a_mol(natom_local_a),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "atom(ikind)%list_a_mol",&
                                           natom_local_a*int_size)
          atom(ikind)%list_a_mol(:) = 0

!         *** Build index vector for mapping ***

          DO iatom_local=1,natom_local_a
            atom_a = local_rows_array(iatom_local)
            DO iatom=1,natom_a
              IF (atom_a == atom(ikind)%list(iatom)) THEN
                atom(ikind)%list_local_a_index(iatom_local) = iatom
                EXIT
              END IF
            END DO
            IF ( molecule_only ) THEN
              DO imol=1,SIZE(molecule_set)
                IF ( atom_a >= molecule_set(imol)%first_atom .AND. &
                     atom_a <= molecule_set(imol)%last_atom ) THEN
                  atom(ikind)%list_a_mol(iatom_local) = imol
                  EXIT
                END IF
              END DO
            END IF
          END DO

        END IF

!       *** Block columns ***

        IF (natom_local_b > 0) THEN

          ALLOCATE (atom(ikind)%list_local_b_index(natom_local_b),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "atom(ikind)%list_local_b_index",&
                                           natom_local_b*int_size)
          ALLOCATE (atom(ikind)%list_b_mol(natom_local_b),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "atom(ikind)%list_b_mol",&
                                           natom_local_b*int_size)
          atom(ikind)%list_b_mol(:) = 0

!         *** Build index vector for mapping ***

          DO iatom_local=1,natom_local_b
            atom_a = local_cols_array(iatom_local)
            DO iatom=1,natom_a
              IF (atom_a == atom(ikind)%list(iatom)) THEN
                atom(ikind)%list_local_b_index(iatom_local) = iatom
                EXIT
              END IF
            END DO
            IF ( molecule_only ) THEN
              DO imol=1,SIZE(molecule_set)
                IF ( atom_a >= molecule_set(imol)%first_atom .AND. &
                     atom_a <= molecule_set(imol)%last_atom ) THEN
                  atom(ikind)%list_b_mol(iatom_local) = imol
                  EXIT
                END IF
              END DO
            END IF
          END DO

        END IF

        ALLOCATE (atom(ikind)%r_pbc(3,natom_a),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%r_pbc",3*natom_a*dp_size)

        ALLOCATE (atom(ikind)%s_pbc(3,natom_a),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%s_pbc",3*natom_a*dp_size)

!       *** Calculate PBC coordinates ***

        DO iatom=1,natom_a
          atom_a = atom(ikind)%list(iatom)
          r_pbc(:) = pbc(particle_set(atom_a)%r(:),cell)
          atom(ikind)%r_pbc(:,iatom) = r_pbc(:)
          atom(ikind)%s_pbc(:,iatom) = real_to_scaled(r_pbc(:),cell)
        END DO

      ELSE

        orb_present(ikind) = .FALSE.

      END IF

      IF (gth_potential_present) THEN
        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             ppl_present=ppl_present(ikind),&
                             ppl_radius=ppl_radius(ikind),&
                             ppnl_present=ppnl_present(ikind),&
                             ppnl_radius=ppnl_radius(ikind))
        ELSE
          ppl_present(ikind) = .FALSE.
          ppnl_present(ikind) = .FALSE.
        END IF
      END IF

!     *** Check the presence of rho1_hard and rho1_soft ***
      IF (paw_atom_present) THEN
        IF (paw_atom) THEN
          oce_present(ikind) = .TRUE.
          CALL get_paw_proj_set(paw_proj_set=paw_proj,&
                                rcprj=oce_radius(ikind))
        ELSE
          oce_present(ikind) = .FALSE.
        END IF
      END IF

!     *** Check the presence of rho0_hard and rho0_soft ***
      IF(rho0_present) THEN
        IF(ABS(rhoz_set(ikind)%one_atom) > 1.E-8_dp) THEN
          Qlm_present(ikind) = .TRUE.
          CALL get_rho0_mpole(rho0_mpole=rho0_mpole, ikind=ikind,&
                              rpgf0_s=V0hmV0s_radius(ikind))
!                              r_V0hmV0s=V0hmV0s_radius(ikind))
        END IF
        orb_rad_3c(ikind) = rho0_mpole%orb_radius_3c(ikind)
      END IF

    END DO

!   *** Build the orbital-orbital overlap neighbor lists ***

    IF ( molecule_only ) THEN
      CALL build_restricted_neighbor_lists(sab_orb,atom,cell,r_skin,&
                 maxatom_local,print_subcell_grid,output_unit,&
                 orb_present,orb_present,orb_radius,orb_radius)
    ELSE
      CALL build_neighbor_lists(sab_orb,atom,cell,r_skin,maxatom_local,&
                 print_subcell_grid,output_unit,&
                 orb_present,orb_present,orb_radius,orb_radius)
    END IF

!   *** Put the generated neighbor lists into the QS environment ***

    CALL set_qs_env(qs_env=qs_env,sab_orb=sab_orb)

    IF (globenv%print%sab_orb_neighbor_lists) THEN
      CALL write_neighbor_lists(sab_orb,particle_set,cell,globenv,&
                                "ORBITAL ORBITAL NEIGHBOR LISTS")
    END IF

!   *** Build orbital GTH-PPL operator overlap list ***

    IF (gth_potential_present) THEN

      IF (ANY(ppl_present)) THEN

        CALL build_gth_ppl_neighbor_lists(sac_ppl,atom,cell,r_skin,&
                                  maxatom_local,print_subcell_grid,&
                                  output_unit,orb_present,ppl_present,&
                                  orb_radius,ppl_radius)

!       *** Put the generated neighbor lists into the QS environment ***

        CALL set_qs_env(qs_env=qs_env,sac_ppl=sac_ppl)

        IF (globenv%print%sac_ppl_neighbor_lists) THEN
          CALL write_neighbor_lists(sac_ppl,particle_set,cell,globenv,&
                                    "ORBITAL GTH-PPL NEIGHBOR LISTS")
        END IF

      END IF

      IF (ANY(ppnl_present)) THEN
        name = " ORBITAL GTH-PPNL"
        CALL build_gth_ppnl_neighbor_lists(sac_ppnl,sbc_ppnl,atom,cell,r_skin,&
                                      maxatom,maxatom_local,print_subcell_grid,&
                                      output_unit,orb_present,ppnl_present,&
                                      orb_radius,ppnl_radius,name)

!       *** Put the generated neighbor lists into the QS environment ***

        CALL set_qs_env(qs_env=qs_env,sac_ppnl=sac_ppnl,sbc_ppnl=sbc_ppnl)

        IF (globenv%print%sac_ppnl_neighbor_lists) THEN
          CALL write_neighbor_lists(sac_ppnl,particle_set,cell,globenv,&
                                    "ORBITAL(A) GTH-PPNL NEIGHBOR LISTS")
        END IF

        IF (globenv%print%sbc_ppnl_neighbor_lists) THEN
          CALL write_neighbor_lists(sbc_ppnl,particle_set,cell,globenv,&
                                    "ORBITAL(B) GTH-PPNL NEIGHBOR LISTS")
        END IF

      END IF

    END IF

!   *** Build orbital orbital-gapw projector  overlap list ***

    IF(paw_atom_present) THEN
      IF(ANY(oce_present)) THEN
        name = " ORBITAL GTH-PAW PRJ"
        CALL build_gth_ppnl_neighbor_lists(sac_oce,sbc_oce,atom,cell,r_skin,&
                                 maxatom,maxatom_local,print_subcell_grid,&
                                 output_unit,orb_present,oce_present,&
                                 orb_radius,oce_radius,name)

!       *** Put the generated neighbor lists into the QS environment ***

        CALL set_qs_env(qs_env=qs_env,sac_oce=sac_oce,sbc_oce=sbc_oce)

        IF (globenv%print%sac_oce_neighbor_lists) THEN
          CALL write_neighbor_lists(sac_oce,particle_set,cell,globenv,&
                                   "ORBITAL(A) PAW-PRJ NEIGHBOR LISTS")
        END IF
        IF (globenv%print%sbc_oce_neighbor_lists) THEN
          CALL write_neighbor_lists(sbc_oce,particle_set,cell,globenv,&
                                   "ORBITAL(B) PAW-PRJ NEIGHBOR LISTS")
        END IF
      ENDIF
    END IF
    IF (rho0_present) THEN

!     *** Build the neighbor lists for the 2-centers Hartree_local terms ***
!     ***              Vh(rho0_h)rho0_h - Vh(rho0_s)rho0_s               ***

      CALL build_neighbor_lists(sab_2c,atom,cell,r_skin,maxatom_local,&
                          .FALSE.,output_unit,Qlm_present,Qlm_present,&
                          V0hmV0s_radius,V0hmV0s_radius)

!     *** Put the generated neighbor lists into the QS environment ***

      CALL set_qs_env(qs_env=qs_env,sab_2c=sab_2c)

      IF (globenv%print%sac_oce_neighbor_lists) THEN
         CALL write_neighbor_lists(sab_2c,particle_set,cell,globenv,&
                                 "2C VH(RHO0) NEIGHBOR LISTS")
      END IF

      name = " ORBITAL GTH-3C VH" 

      CALL build_gth_ppnl_neighbor_lists(sac_3c,sbc_3c,atom,cell,r_skin,&
                               maxatom,maxatom_local,print_subcell_grid,&
                               output_unit,orb_present,Qlm_present,&
                               orb_rad_3c,V0hmV0s_radius,name)

!     *** Put the generated neighbor lists into the QS environment ***

      CALL set_qs_env(qs_env=qs_env,sac_3c=sac_3c,sbc_3c=sbc_3c)

      IF (globenv%print%sac_oce_neighbor_lists) THEN
         CALL write_neighbor_lists(sac_3c,particle_set,cell,globenv,&
                             "ORBITAL(A) 3-CENTERS-VH NEIGHBOR LISTS")
      END IF
      IF (globenv%print%sbc_oce_neighbor_lists) THEN
        CALL write_neighbor_lists(sbc_3c,particle_set,cell,globenv,&
                             "ORBITAL(B) 3-CENTERS-VH NEIGHBOR LISTS")
      END IF

    END IF

!   *** Print particle distribution ***

    CALL write_distribution(sab_orb,atomic_kind_set,globenv)

!   *** Release work storage ***

    DO ikind=1,nkind
      NULLIFY (atom(ikind)%list)
      IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_a_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%list_local_a_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%list_local_b_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_b_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%list_local_b_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%list_a_mol)) THEN
        DEALLOCATE (atom(ikind)%list_a_mol,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%list_a_mol")
      END IF
      IF (ASSOCIATED(atom(ikind)%list_b_mol)) THEN
        DEALLOCATE (atom(ikind)%list_b_mol,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%list_b_mol")
      END IF
      IF (ASSOCIATED(atom(ikind)%r_pbc)) THEN
        DEALLOCATE (atom(ikind)%r_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%r_pbc")
      END IF
      IF (ASSOCIATED(atom(ikind)%s_pbc)) THEN
        DEALLOCATE (atom(ikind)%s_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atom(ikind)%s_pbc")
      END IF
    END DO

    DEALLOCATE (atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "atom")

    DEALLOCATE (orb_present,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "orb_present")

    DEALLOCATE (orb_radius,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "orb_radius")

    IF (gth_potential_present) THEN

      DEALLOCATE (ppl_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppl_present")

      DEALLOCATE (ppl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppl_radius")

      DEALLOCATE (ppnl_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppnl_present")

      DEALLOCATE (ppnl_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "ppnl_radius")

    END IF

    IF (paw_atom_present) THEN

      DEALLOCATE (oce_present,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "oce_present")

      DEALLOCATE (oce_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "oce_radius")

      DEALLOCATE (V0hmV0s_radius,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "V0hmV0s_radius")

      DEALLOCATE(orb_rad_3c,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "hard_radius")
    ENDIF

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_qs_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_neighbor_lists(ab,atom,cell,r_skin,maxatom_local,&
                                  print_subcell_grid,output_unit,&
                                  present_a,present_b,radius_a,radius_b)

!   Purpose: Build simple pair neighbor lists.

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(IN)                     :: r_skin
    INTEGER, INTENT(IN)                      :: maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, DIMENSION(:), INTENT(IN)        :: present_a, present_b
    REAL(dp), DIMENSION(:), INTENT(IN)       :: radius_a, radius_b

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_neighbor_lists"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_b, handle, i, iab, iatom, iatom_local, &
      iatom_subcell, icell, ikind, istat, j, jatom, jatom_local, jcell, &
      jkind, k, kcell, natom_local_a, natom_local_b, nkind
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b, ncell, &
                                                nsubcell, periodic
    LOGICAL                                  :: include_ab
    REAL(dp)                                 :: rab2, rab2_max, rab_max, &
                                                subcells, unit_of_length
    REAL(dp), DIMENSION(3)                   :: r, rab, rb, s, sab_max, sb, &
                                                sb_max, sb_min, sb_pbc
    TYPE(neighbor_list_p_type), &
      DIMENSION(maxatom_local)               :: kind_a
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I"," ",handle)

!   *** Deallocate the old neighbor list structure ***

    IF (ASSOCIATED(ab)) THEN
      DO iab=1,SIZE(ab)
        CALL deallocate_neighbor_list_set(ab(iab)%neighbor_list_set)
      END DO
      DEALLOCATE (ab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"ab")
    END IF

!   *** Allocate and initialize the new neighbor list structure ***

    nkind = SIZE(atom)

    ALLOCATE (ab(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ab",nkind*nkind*int_size)

    DO iab=1,SIZE(ab)
      NULLIFY (ab(iab)%neighbor_list_set)
    END DO

!   *** Load informations about the simulation cell ***

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Print headline ***

    IF (print_subcell_grid) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE ORBITAL ORBITAL NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.present_a(ikind)) CYCLE

      IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE

      natom_local_a = SIZE(atom(ikind)%list_local_a_index)

      DO jkind=1,nkind

        IF (.NOT.present_b(jkind)) CYCLE

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(atom(jkind)%list_local_b_index)) CYCLE

        natom_local_b = SIZE(atom(jkind)%list_local_b_index)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius_a(ikind) + radius_b(jkind) + r_skin
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
        nsubcell(:) = MAX(1,NINT(0.5_dp*subcells/sab_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_subcell_grid) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),dp)
          r(:) = scaled_to_real(s(:),cell)
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,nsubcell(1:3),r(1:3)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                          ab(iab)%neighbor_list_set,&
                                        r_max=rab_max)
        neighbor_list_set => ab(iab)%neighbor_list_set

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(nsubcell) == 3) THEN

!         *** Case 1: No subcell grid is needed ***

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local=1,natom_local_b

            jatom = atom(jkind)%list_local_b_index(jatom_local)
            atom_b = atom(jkind)%list(jatom)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom)

            loop1_kcell: DO kcell=-ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              IF (sb_min(3) >= 0.5_dp) EXIT loop1_kcell
              sb_max(3) = sb(3) + sab_max(3)
              IF (sb_max(3) < -0.5_dp) CYCLE loop1_kcell
              cell_b(3) = kcell

              loop1_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                IF (sb_min(2) >= 0.5_dp) EXIT loop1_jcell
                sb_max(2) = sb(2) + sab_max(2)
                IF (sb_max(2) < -0.5_dp) CYCLE loop1_jcell
                cell_b(2) = jcell

                loop1_icell: DO icell=-ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  IF (sb_min(1) >= 0.5_dp) EXIT loop1_icell
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (sb_max(1) < -0.5_dp) CYCLE loop1_icell
                  cell_b(1) = icell

                  rb(:) = scaled_to_real(sb(:),cell)

                  DO iatom_local=1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)
                    IF (atom_a > atom_b) THEN
                      include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                    ELSE
                      include_ab = (MODULO(atom_a + atom_b,2) == 0)
                    END IF
                    IF (include_ab) THEN
                      rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                      IF (rab2 < rab2_max) THEN
                        CALL add_neighbor_node(&
                          neighbor_list=kind_a(iatom_local)%neighbor_list,&
                          neighbor=atom_b,&
                          cell=cell_b,&
                          r=rab)
                      END IF
                    END IF
                  END DO

                END DO loop1_icell
              END DO loop1_jcell
            END DO loop1_kcell

          END DO

        ELSE

!         *** Case 2: A subcell grid is needed ***

          CALL allocate_subcell(subcell,nsubcell,maxatom_local)

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = INT(s(1)*REAL(nsubcell(1),dp) + 1.0_dp)
            j = INT(s(2)*REAL(nsubcell(2),dp) + 1.0_dp)
            k = INT(s(3)*REAL(nsubcell(3),dp) + 1.0_dp)
            i = MIN(MAX(i,1),nsubcell(1))
            j = MIN(MAX(j,1),nsubcell(2))
            k = MIN(MAX(k,1),nsubcell(3))
            subcell(i,j,k)%natom =  subcell(i,j,k)%natom + 1
            subcell(i,j,k)%atom_list(subcell(i,j,k)%natom) = iatom_local
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local=1,natom_local_b

            jatom = atom(jkind)%list_local_b_index(jatom_local)
            atom_b = atom(jkind)%list(jatom)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom)

            loop2_kcell: DO kcell=-ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              IF (sb_min(3) >= 0.5_dp) EXIT loop2_kcell
              sb_max(3) = sb(3) + sab_max(3)
              IF (sb_max(3) < -0.5_dp) CYCLE loop2_kcell
              cell_b(3) = kcell

              loop2_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                IF (sb_min(2) >= 0.5_dp) EXIT loop2_jcell
                sb_max(2) = sb(2) + sab_max(2)
                IF (sb_max(2) < -0.5_dp) CYCLE loop2_jcell
                cell_b(2) = jcell

                loop2_icell: DO icell=-ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  IF (sb_min(1) >= 0.5_dp) EXIT loop2_icell
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (sb_max(1) < -0.5_dp) CYCLE loop2_icell
                  cell_b(1) = icell

                  rb(:) = scaled_to_real(sb(:),cell)

                  loop_k: DO k=1,nsubcell(3)
                    loop_j: DO j=1,nsubcell(2)
                      loop_i: DO i=1,nsubcell(1)

                        IF (subcell(i,j,k)%natom == 0) CYCLE

                        IF (sb_max(3) < subcell(i,j,k)%s_min(3)) EXIT loop_k
                        IF (sb_min(3) >= subcell(i,j,k)%s_max(3)) CYCLE loop_k

                        IF (sb_max(2) < subcell(i,j,k)%s_min(2)) EXIT loop_j
                        IF (sb_min(2) >= subcell(i,j,k)%s_max(2)) CYCLE loop_j

                        IF (sb_max(1) < subcell(i,j,k)%s_min(1)) EXIT loop_i
                        IF (sb_min(1) >= subcell(i,j,k)%s_max(1)) CYCLE loop_i

                        DO iatom_subcell=1,subcell(i,j,k)%natom
                          iatom_local = subcell(i,j,k)%atom_list(iatom_subcell)
                          iatom = atom(ikind)%list_local_a_index(iatom_local)
                          atom_a = atom(ikind)%list(iatom)
                          IF (atom_a > atom_b) THEN
                            include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                          ELSE
                            include_ab = (MODULO(atom_a + atom_b,2) == 0)
                          END IF
                          IF (include_ab) THEN
                            rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                            IF (rab2 < rab2_max) THEN
                              CALL add_neighbor_node(&
                                neighbor_list=kind_a(iatom_local)%neighbor_list,&
                                neighbor=atom_b,&
                                cell=cell_b,&
                                r=rab)
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

! *****************************************************************************

  SUBROUTINE build_restricted_neighbor_lists(ab,atom,cell,r_skin,maxatom_local,&
                                  print_subcell_grid,output_unit,&
                                  present_a,present_b,radius_a,radius_b)

!   Purpose: Build simple pair neighbor lists restricted to molecules.

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(IN)                     :: r_skin
    INTEGER, INTENT(IN)                      :: maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, DIMENSION(:), INTENT(IN)        :: present_a, present_b
    REAL(dp), DIMENSION(:), INTENT(IN)       :: radius_a, radius_b

    CHARACTER(LEN=*), PARAMETER :: routine_name = &
                                   "build_restricted_neighbor_lists"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_b, handle, i, iab, iatom, iatom_local, &
      iatom_subcell, icell, ikind, istat, j, jatom, jatom_local, jcell, &
      jkind, k, kcell, natom_local_a, natom_local_b, nkind, mol_a, mol_b
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b, ncell, &
                                                nsubcell, periodic
    LOGICAL                                  :: include_ab
    REAL(dp)                                 :: rab2, rab2_max, rab_max, &
                                                subcells, unit_of_length
    REAL(dp), DIMENSION(3)                   :: r, rab, rb, s, sab_max, sb, &
                                                sb_max, sb_min, sb_pbc
    TYPE(neighbor_list_p_type), &
      DIMENSION(maxatom_local)               :: kind_a
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I"," ",handle)

!   *** Deallocate the old neighbor list structure ***

    IF (ASSOCIATED(ab)) THEN
      DO iab=1,SIZE(ab)
        CALL deallocate_neighbor_list_set(ab(iab)%neighbor_list_set)
      END DO
      DEALLOCATE (ab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"ab")
    END IF

!   *** Allocate and initialize the new neighbor list structure ***

    nkind = SIZE(atom)

    ALLOCATE (ab(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ab",nkind*nkind*int_size)

    DO iab=1,SIZE(ab)
      NULLIFY (ab(iab)%neighbor_list_set)
    END DO

!   *** Load informations about the simulation cell ***

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Print headline ***

    IF (print_subcell_grid) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE ORBITAL ORBITAL NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.present_a(ikind)) CYCLE

      IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE

      natom_local_a = SIZE(atom(ikind)%list_local_a_index)

      DO jkind=1,nkind

        IF (.NOT.present_b(jkind)) CYCLE

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(atom(jkind)%list_local_b_index)) CYCLE

        natom_local_b = SIZE(atom(jkind)%list_local_b_index)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius_a(ikind) + radius_b(jkind) + r_skin
        rab2_max = rab_max*rab_max

        r(:) = rab_max
        sab_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = 0
        nsubcell(:) = MAX(1,NINT(0.5_dp*subcells/sab_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_subcell_grid) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),dp)
          r(:) = scaled_to_real(s(:),cell)
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,jkind,nsubcell(1:3),r(1:3)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                        ab(iab)%neighbor_list_set,r_max=rab_max)
        neighbor_list_set => ab(iab)%neighbor_list_set

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(nsubcell) == 3) THEN

!         *** Case 1: No subcell grid is needed ***

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local=1,natom_local_b

            jatom = atom(jkind)%list_local_b_index(jatom_local)
            atom_b = atom(jkind)%list(jatom)
            mol_b = atom(jkind)%list_b_mol(jatom_local)
            sb(:) = atom(jkind)%s_pbc(:,jatom)
            rb(:) = scaled_to_real(sb(:),cell)
            cell_b = 0

            DO iatom_local=1,natom_local_a
              iatom = atom(ikind)%list_local_a_index(iatom_local)
              atom_a = atom(ikind)%list(iatom)
              mol_a = atom(ikind)%list_a_mol(iatom_local)
              IF ( mol_a /= mol_b ) CYCLE
              IF (atom_a > atom_b) THEN
                include_ab = (MODULO(atom_a + atom_b,2) /= 0)
              ELSE
                include_ab = (MODULO(atom_a + atom_b,2) == 0)
              END IF
              IF (include_ab) THEN
                rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                r(:) = pbc(rab(:),cell)
                rab2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
                IF (rab2 < rab2_max) THEN
                  CALL add_neighbor_node(&
                    neighbor_list=kind_a(iatom_local)%neighbor_list,&
                    neighbor=atom_b,cell=cell_b,r=rab)
                END IF
              END IF
            END DO

          END DO

        ELSE

!         *** Case 2: A subcell grid is needed ***

          CALL allocate_subcell(subcell,nsubcell,maxatom_local)

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = INT(s(1)*REAL(nsubcell(1),dp) + 1.0_dp)
            j = INT(s(2)*REAL(nsubcell(2),dp) + 1.0_dp)
            k = INT(s(3)*REAL(nsubcell(3),dp) + 1.0_dp)
            i = MIN(MAX(i,1),nsubcell(1))
            j = MIN(MAX(j,1),nsubcell(2))
            k = MIN(MAX(k,1),nsubcell(3))
            subcell(i,j,k)%natom =  subcell(i,j,k)%natom + 1
            subcell(i,j,k)%atom_list(subcell(i,j,k)%natom) = iatom_local
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO jatom_local=1,natom_local_b

            jatom = atom(jkind)%list_local_b_index(jatom_local)
            atom_b = atom(jkind)%list(jatom)
            mol_b = atom(jkind)%list_b_mol(jatom_local)
            sb(:) = atom(jkind)%s_pbc(:,jatom)
            rb(:) = scaled_to_real(sb(:),cell)
            cell_b(:) = 0
            sb_min(:) = sb(:) - sab_max(:)
            sb_max(:) = sb(:) + sab_max(:)

            loop_k: DO k=1,nsubcell(3)
              loop_j: DO j=1,nsubcell(2)
                loop_i: DO i=1,nsubcell(1)

                  IF (subcell(i,j,k)%natom == 0) CYCLE

                  IF (sb_max(3) < subcell(i,j,k)%s_min(3)) EXIT loop_k
                  IF (sb_min(3) >= subcell(i,j,k)%s_max(3)) CYCLE loop_k

                  IF (sb_max(2) < subcell(i,j,k)%s_min(2)) EXIT loop_j
                  IF (sb_min(2) >= subcell(i,j,k)%s_max(2)) CYCLE loop_j

                  IF (sb_max(1) < subcell(i,j,k)%s_min(1)) EXIT loop_i
                  IF (sb_min(1) >= subcell(i,j,k)%s_max(1)) CYCLE loop_i

                  DO iatom_subcell=1,subcell(i,j,k)%natom
                    iatom_local = subcell(i,j,k)%atom_list(iatom_subcell)
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)
                    mol_a = atom(ikind)%list_a_mol(iatom_local)
                    IF ( mol_a /= mol_b ) CYCLE
                    IF (atom_a > atom_b) THEN
                      include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                    ELSE
                      include_ab = (MODULO(atom_a + atom_b,2) == 0)
                    END IF
                    IF (include_ab) THEN
                      rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                      r(:) = pbc(rab(:),cell)
                      rab2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
                      IF (rab2 < rab2_max) THEN
                        CALL add_neighbor_node(&
                          neighbor_list=kind_a(iatom_local)%neighbor_list,&
                          neighbor=atom_b,&
                          cell=cell_b,&
                          r=rab)
                      END IF
                    END IF
                  END DO

                END DO loop_i
              END DO loop_j
            END DO loop_k

          END DO

          CALL deallocate_subcell(subcell)

        END IF

      END DO
    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_restricted_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_gth_ppl_neighbor_lists(ac,atom,cell,r_skin,maxatom_local,&
                                          print_subcell_grid,output_unit,&
                                          present_a,present_c,&
                                          radius_a,radius_c, particle_set)

!   Purpose: Build the neighbor lists for the calculation of the local part of
!            the GTH pseudo potential (PPL). Atom c is a PPL operator atom,
!            atom a is an atom in the cell(0,0,0).

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ac
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(IN)                     :: r_skin
    INTEGER, INTENT(IN)                      :: maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, DIMENSION(:), INTENT(IN)        :: present_a, present_c
    REAL(dp), DIMENSION(:), INTENT(IN)       :: radius_a, radius_c
    TYPE(particle_type), DIMENSION(:), &
      POINTER, OPTIONAL                      :: particle_set


    CHARACTER(LEN=*), PARAMETER :: &
      routine_name = "build_gth_ppl_neighbor_lists"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_c, handle, i, iac, iatom, iatom_local, &
      iatom_subcell, icell, iexl, ikind, istat, j, jcell, jkind, k,&
      katom, kcell, kkind, nexl, &
      natom_c, natom_local_a, nkind
    INTEGER, DIMENSION(3)                    :: cell_a, cell_c, ncell, &
                                                nsubcell, periodic
    INTEGER, DIMENSION ( : ), POINTER        :: list
    LOGICAL                                  :: cell000, exclude
    REAL(dp)                                 :: rac2, rac2_max, rac_max, &
                                                subcells, unit_of_length
    REAL(dp), DIMENSION(3)                   :: r, rac, rc, s, sac_max, sc, &
                                                sc_max, sc_min, sc_pbc
    TYPE(neighbor_list_p_type), &
      DIMENSION(maxatom_local)               :: kind_a
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I"," ",handle)

!   *** Deallocate the old neighbor list structure ***

    IF (ASSOCIATED(ac)) THEN
      DO iac=1,SIZE(ac)
        CALL deallocate_neighbor_list_set(ac(iac)%neighbor_list_set)
      END DO
      DEALLOCATE (ac,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"ac")
    END IF

!   *** Allocate and initialize the new neighbor list structure ***

    nkind = SIZE(atom)

    ALLOCATE (ac(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ac",nkind*nkind*int_size)

    DO iac=1,SIZE(ac)
      NULLIFY (ac(iac)%neighbor_list_set)
    END DO

!   *** Load informations about the simulation cell ***

    exclude=.FALSE.
    NULLIFY (list)                 
   
!   *** Load informations about the simulation cell ***

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Print headline ***

    IF (print_subcell_grid) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE ORBITAL GTH-PPL NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.present_a(ikind)) CYCLE

      IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE

      natom_local_a = SIZE(atom(ikind)%list_local_a_index)

      DO kkind=1,nkind

        IF (.NOT.present_c(kkind)) CYCLE

        iac = ikind + nkind*(kkind - 1)

        natom_c = SIZE(atom(kkind)%list)

!       *** Calculate the square of the maximum interaction distance ***

        rac_max = radius_a(ikind) + radius_c(kkind) + r_skin
        rac2_max = rac_max*rac_max

        r(:) = rac_max
        sac_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sac_max(:)) + 1)*periodic(:)
        nsubcell(:) = MAX(1,NINT(0.5_dp*subcells/sac_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_subcell_grid) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),dp)
          r(:) = scaled_to_real(s(:),cell)
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,kkind,nsubcell(1:3),r(1:3)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=&
                                        ac(iac)%neighbor_list_set,r_max=rac_max)
        neighbor_list_set => ac(iac)%neighbor_list_set

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(nsubcell) == 3) THEN

!         *** Case 1: No subcell grid is needed ***

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO katom=1,natom_c

            atom_c = atom(kkind)%list(katom)
            sc_pbc(:) = atom(kkind)%s_pbc(:,katom)

            loop1_kcell: DO kcell=-ncell(3),ncell(3)

              sc(3) = sc_pbc(3) + REAL(kcell,dp)
              sc_min(3) = sc(3) - sac_max(3)
              IF (sc_min(3) >= 0.5_dp) EXIT loop1_kcell
              sc_max(3) = sc(3) + sac_max(3)
              IF (sc_max(3) < -0.5_dp) CYCLE loop1_kcell
              cell_c(3) = kcell

              loop1_jcell: DO jcell=-ncell(2),ncell(2)

                sc(2) = sc_pbc(2) + REAL(jcell,dp)
                sc_min(2) = sc(2) - sac_max(2)
                IF (sc_min(2) >= 0.5_dp) EXIT loop1_jcell
                sc_max(2) = sc(2) + sac_max(2)
                IF (sc_max(2) < -0.5_dp) CYCLE loop1_jcell
                cell_c(2) = jcell

                loop1_icell: DO icell=-ncell(1),ncell(1)

                  sc(1) = sc_pbc(1) + REAL(icell,dp)
                  sc_min(1) = sc(1) - sac_max(1)
                  IF (sc_min(1) >= 0.5_dp) EXIT loop1_icell
                  sc_max(1) = sc(1) + sac_max(1)
                  IF (sc_max(1) < -0.5_dp) CYCLE loop1_icell
                  cell_c(1) = icell

                  rc(:) = scaled_to_real(sc(:),cell)

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  DO iatom_local=1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)

! *** check for exclusion ***                          

                    IF (PRESENT(particle_set)) THEN
                      loop_out1: DO jkind = 1, nkind
                        list => particle_set(atom_a)%lpi%exclusion(jkind) % list
                        IF (.NOT.ASSOCIATED(list)) CYCLE
                        nexl = SIZE ( list )
                        loop_in1: DO iexl = 1, nexl
                          IF ((list(iexl)==atom_c).AND.cell000) THEN
                            exclude = .TRUE.
                            EXIT loop_out1
                          ELSE
                            exclude = .FALSE.
                          END IF
                        END DO loop_in1
                      END DO loop_out1
                      IF(exclude) CYCLE
                    END IF
!
                    rac(:) = rc(:) - atom(ikind)%r_pbc(:,iatom)
                    rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                    IF (rac2 < rac2_max) THEN
                      CALL add_neighbor_node(&
                        neighbor_list=kind_a(iatom_local)%neighbor_list,&
                        neighbor=atom_c,&
                        cell=cell_c,&
                        r=rac)
                    END IF
                  END DO

                END DO loop1_icell
              END DO loop1_jcell
            END DO loop1_kcell

          END DO

        ELSE

!         *** Case 2: A subcell grid is needed ***

          CALL allocate_subcell(subcell,nsubcell,maxatom_local)

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = INT(s(1)*REAL(nsubcell(1),dp) + 1.0_dp)
            j = INT(s(2)*REAL(nsubcell(2),dp) + 1.0_dp)
            k = INT(s(3)*REAL(nsubcell(3),dp) + 1.0_dp)
            i = MIN(MAX(i,1),nsubcell(1))
            j = MIN(MAX(j,1),nsubcell(2))
            k = MIN(MAX(k,1),nsubcell(3))
            subcell(i,j,k)%natom =  subcell(i,j,k)%natom + 1
            subcell(i,j,k)%atom_list(subcell(i,j,k)%natom) = iatom_local
            CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO katom=1,natom_c

            atom_c = atom(kkind)%list(katom)
            sc_pbc(:) = atom(kkind)%s_pbc(:,katom)

            loop2_kcell: DO kcell=-ncell(3),ncell(3)

              sc(3) = sc_pbc(3) + REAL(kcell,dp)
              sc_min(3) = sc(3) - sac_max(3)
              IF (sc_min(3) >= 0.5_dp) EXIT loop2_kcell
              sc_max(3) = sc(3) + sac_max(3)
              IF (sc_max(3) < -0.5_dp) CYCLE loop2_kcell
              cell_c(3) = kcell

              loop2_jcell: DO jcell=-ncell(2),ncell(2)

                sc(2) = sc_pbc(2) + REAL(jcell,dp)
                sc_min(2) = sc(2) - sac_max(2)
                IF (sc_min(2) >= 0.5_dp) EXIT loop2_jcell
                sc_max(2) = sc(2) + sac_max(2)
                IF (sc_max(2) < -0.5_dp) CYCLE loop2_jcell
                cell_c(2) = jcell

                loop2_icell: DO icell=-ncell(1),ncell(1)

                  sc(1) = sc_pbc(1) + REAL(icell,dp)
                  sc_min(1) = sc(1) - sac_max(1)
                  IF (sc_min(1) >= 0.5_dp) EXIT loop2_icell
                  sc_max(1) = sc(1) + sac_max(1)
                  IF (sc_max(1) < -0.5_dp) CYCLE loop2_icell
                  cell_c(1) = icell

                  rc(:) = scaled_to_real(sc(:),cell)

                  cell000 = ((icell == 0).AND.&
                             (jcell == 0).AND.&
                             (kcell == 0))

                  loop_k: DO k=1,nsubcell(3)
                    loop_j: DO j=1,nsubcell(2)
                      loop_i: DO i=1,nsubcell(1)

                        IF (subcell(i,j,k)%natom == 0) CYCLE

                        IF (sc_max(3) < subcell(i,j,k)%s_min(3)) EXIT loop_k
                        IF (sc_min(3) >= subcell(i,j,k)%s_max(3)) CYCLE loop_k

                        IF (sc_max(2) < subcell(i,j,k)%s_min(2)) EXIT loop_j
                        IF (sc_min(2) >= subcell(i,j,k)%s_max(2)) CYCLE loop_j

                        IF (sc_max(1) < subcell(i,j,k)%s_min(1)) EXIT loop_i
                        IF (sc_min(1) >= subcell(i,j,k)%s_max(1)) CYCLE loop_i

                        DO iatom_subcell=1,subcell(i,j,k)%natom
                          iatom_local = subcell(i,j,k)%atom_list(iatom_subcell)
                          iatom = atom(ikind)%list_local_a_index(iatom_local)
                          atom_a = atom(ikind)%list(iatom)

! *** check for exclusion ***

                          IF (PRESENT(particle_set)) THEN
                            loop_out2: DO jkind = 1, nkind
                              list => particle_set(atom_a)%lpi%exclusion(jkind)%list
                              IF (.NOT.ASSOCIATED(list)) CYCLE
                              nexl = SIZE (list)
                              loop_in2: DO iexl = 1, nexl
                                IF ((list(iexl)==atom_c).AND.cell000) THEN
                                  exclude = .TRUE.
                                  EXIT loop_out2
                                ELSE
                                  exclude = .FALSE.
                                END IF
                              END DO loop_in2
                            END DO loop_out2
                            IF(exclude) CYCLE
                          END IF

                          rac(:) = rc(:) - atom(ikind)%r_pbc(:,iatom)
                          rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                          IF (rac2 < rac2_max) THEN
                            CALL add_neighbor_node(&
                              neighbor_list=kind_a(iatom_local)%neighbor_list,&
                              neighbor=atom_c,&
                              cell=cell_c,&
                              r=rac)
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

  END SUBROUTINE build_gth_ppl_neighbor_lists

! *****************************************************************************

  SUBROUTINE build_gth_ppnl_neighbor_lists(ac,bc,atom,cell,r_skin,&
                                           maxatom,maxatom_local,&
                                           print_subcell_grid,output_unit,&
                                           present_a,present_c,&
                                           radius_a,radius_c,name)

!   Purpose: Build the neighbor lists for the calculation of the non-local part
!            of the GTH pseudo potential (PPNL). Atom c is a PPNL operator
!            atom, atom a is an atom in the cell(0,0,0), and atom b is an atom
!            in any cell(i,j,k).

!   History: - Creation (30.05.2003,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ac, bc
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), INTENT(IN)                     :: r_skin
    INTEGER, INTENT(IN)                      :: maxatom, maxatom_local
    LOGICAL, INTENT(IN)                      :: print_subcell_grid
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, DIMENSION(:), INTENT(IN)        :: present_a, present_c
    REAL(dp), DIMENSION(:), INTENT(IN)       :: radius_a, radius_c
    CHARACTER(20)                            :: name

    CHARACTER(LEN=*), PARAMETER :: &
      routine_name = "build_gth_ppnl_neighbor_lists"

    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER :: atom_a, atom_c, handle, i, iac, iatom, iatom_local, &
      iatom_subcell, icell, icell_b, ikind, istat, j, jcell, jcell_b, k, &
      katom, kcell, kcell_b, kkind, natom_c, natom_local_a, natom_local_b, &
      nkind
    INTEGER, DIMENSION(3)                    :: cell_a, cell_c, ncell, &
                                                nsubcell, periodic
    REAL(dp)                                 :: rac2, rac2_max, rac_max, &
                                                subcells, unit_of_length
    REAL(dp), DIMENSION(3)                   :: r, rac, rc, s, s_max, s_min, &
                                                sac_max, sc, sc_max, sc_min, &
                                                sc_pbc
    TYPE(neighbor_list_p_type), &
      DIMENSION(maxatom_local)               :: kind_a, kind_b
    TYPE(neighbor_list_set_type), POINTER    :: ac_neighbor_list_set, &
                                                bc_neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell_a, subcell_b

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I"," ",handle)

!   *** Deallocate the old neighbor list structures ***

    IF (ASSOCIATED(ac)) THEN
      DO iac=1,SIZE(ac)
        CALL deallocate_neighbor_list_set(ac(iac)%neighbor_list_set)
      END DO
      DEALLOCATE (ac,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"ac")
    END IF

    IF (ASSOCIATED(bc)) THEN
      DO iac=1,SIZE(bc)
        CALL deallocate_neighbor_list_set(bc(iac)%neighbor_list_set)
      END DO
      DEALLOCATE (bc,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"bc")
    END IF

!   *** Allocate and initialize the new neighbor list structures ***

    nkind = SIZE(atom)

    ALLOCATE (ac(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ac",nkind*nkind*int_size)

    ALLOCATE (bc(nkind*nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "bc",nkind*nkind*int_size)

    DO iac=1,SIZE(ac)
      NULLIFY (ac(iac)%neighbor_list_set)
      NULLIFY (bc(iac)%neighbor_list_set)
    END DO

!   *** Load informations about the simulation cell ***

    CALL get_cell(cell=cell,&
                  periodic=periodic,&
                  subcells=subcells,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

!   *** Print headline ***

    IF (print_subcell_grid) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/,/,T3,A,T29,A,T54,A)")&
        "SUBCELL GRID FOR THE"//TRIM(name)//" ORBITAL GTH-PPNL NEIGHBOR LISTS",&
        "Atomic kind pair","Grid size",&
        "Subcell size in "//unit_of_length_name
    END IF

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.present_a(ikind)) CYCLE

      IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
        natom_local_a = SIZE(atom(ikind)%list_local_a_index)
      ELSE
        natom_local_a = 0
      END IF

      IF (ASSOCIATED(atom(ikind)%list_local_b_index)) THEN
        natom_local_b = SIZE(atom(ikind)%list_local_b_index)
      ELSE
        natom_local_b = 0
      END IF

      DO kkind=1,nkind

        IF (.NOT.present_c(kkind)) CYCLE

        iac = ikind + nkind*(kkind - 1)

        natom_c = SIZE(atom(kkind)%list)

!       *** Calculate the square of the maximum interaction distance ***

        rac_max = radius_a(ikind) + radius_c(kkind) + r_skin
        rac2_max = rac_max*rac_max

        r(:) = rac_max
        sac_max(:) = real_to_scaled(r(:),cell)

        ncell(:) = (INT(sac_max(:)) + 1)*periodic(:)
        nsubcell(:) = MAX(1,NINT(0.5_dp*subcells/sac_max(:)))

!       *** Print subcell information for the current atomic kind pair ***

        IF (print_subcell_grid) THEN
          s(:) = 1.0_dp/REAL(nsubcell(:),dp)
          r(:) = scaled_to_real(s(:),cell)
          WRITE (UNIT=output_unit,FMT="(T3,2I8,4X,3I5,6X,3F12.6)")&
            ikind,kkind,nsubcell(1:3),r(1:3)/unit_of_length
        END IF

        CALL allocate_neighbor_list_set(neighbor_list_set=ac(iac)%neighbor_list_set,&
                                        r_max=rac_max)
        ac_neighbor_list_set => ac(iac)%neighbor_list_set

        CALL allocate_neighbor_list_set(neighbor_list_set=bc(iac)%neighbor_list_set,&
                                        r_max=rac_max)
        bc_neighbor_list_set => bc(iac)%neighbor_list_set

        cell_a = (/0,0,0/)

!       *** Check, if we have to consider a subcell grid ***

        IF (SUM(nsubcell) == 3) THEN

!         *** Case 1: No subcell grid is needed ***

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=ac_neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          DO iatom_local=1,natom_local_b
            iatom = atom(ikind)%list_local_b_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            CALL add_neighbor_list(neighbor_list_set=bc_neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_b(iatom_local)%neighbor_list)
          END DO

          DO katom=1,natom_c

            atom_c = atom(kkind)%list(katom)
            sc_pbc(:) = atom(kkind)%s_pbc(:,katom)

            loop1_kcell: DO kcell=-ncell(3),ncell(3)

              sc(3) = sc_pbc(3) + REAL(kcell,dp)
              sc_min(3) = sc(3) - sac_max(3)
              IF (sc_min(3) >= 0.5_dp) EXIT loop1_kcell
              sc_max(3) = sc(3) + sac_max(3)
              IF (sc_max(3) < -0.5_dp) CYCLE loop1_kcell
              cell_c(3) = kcell

              loop1_jcell: DO jcell=-ncell(2),ncell(2)

                sc(2) = sc_pbc(2) + REAL(jcell,dp)
                sc_min(2) = sc(2) - sac_max(2)
                IF (sc_min(2) >= 0.5_dp) EXIT loop1_jcell
                sc_max(2) = sc(2) + sac_max(2)
                IF (sc_max(2) < -0.5_dp) CYCLE loop1_jcell
                cell_c(2) = jcell

                loop1_icell: DO icell=-ncell(1),ncell(1)

                  sc(1) = sc_pbc(1) + REAL(icell,dp)
                  sc_min(1) = sc(1) - sac_max(1)
                  IF (sc_min(1) >= 0.5_dp) EXIT loop1_icell
                  sc_max(1) = sc(1) + sac_max(1)
                  IF (sc_max(1) < -0.5_dp) CYCLE loop1_icell
                  cell_c(1) = icell

                  rc(:) = scaled_to_real(sc(:),cell)

                  DO iatom_local=1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    rac(:) = rc(:) - atom(ikind)%r_pbc(:,iatom)
                    rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                    IF (rac2 < rac2_max) THEN
                      CALL add_neighbor_node(&
                        neighbor_list=kind_a(iatom_local)%neighbor_list,&
                        neighbor=atom_c,&
                        cell=cell_c,&
                        r=rac)
                    END IF
                  END DO

                  DO icell_b=cell_c(1)-ncell(1),cell_c(1)+ncell(1)

                    s(1) = REAL(icell_b,dp)

                    DO jcell_b=cell_c(2)-ncell(2),cell_c(2)+ncell(2)

                      s(2) = REAL(jcell_b,dp)

                      DO kcell_b=cell_c(3)-ncell(3),cell_c(3)+ncell(3)

                        s(3) = REAL(kcell_b,dp)

                        r(:) = rc(:) - scaled_to_real(s(:),cell)

                        DO iatom_local=1,natom_local_b
                          iatom = atom(ikind)%list_local_b_index(iatom_local)
                          rac(:) = r(:) - atom(ikind)%r_pbc(:,iatom)
                          rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                          IF (rac2 < rac2_max) THEN
                            CALL add_neighbor_node(&
                              neighbor_list=kind_b(iatom_local)%neighbor_list,&
                              neighbor=atom_c,&
                              cell=cell_c,&
                              r=rac)
                          END IF
                        END DO

                      END DO
                    END DO
                  END DO

                END DO loop1_icell
              END DO loop1_jcell
            END DO loop1_kcell

          END DO

        ELSE

          CALL allocate_subcell(subcell_a,nsubcell,maxatom_local)

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = INT(s(1)*REAL(nsubcell(1),dp) + 1.0_dp)
            j = INT(s(2)*REAL(nsubcell(2),dp) + 1.0_dp)
            k = INT(s(3)*REAL(nsubcell(3),dp) + 1.0_dp)
            i = MIN(MAX(i,1),nsubcell(1))
            j = MIN(MAX(j,1),nsubcell(2))
            k = MIN(MAX(k,1),nsubcell(3))
            subcell_a(i,j,k)%natom =  subcell_a(i,j,k)%natom + 1
            subcell_a(i,j,k)%atom_list(subcell_a(i,j,k)%natom) = iatom_local
            CALL add_neighbor_list(neighbor_list_set=ac_neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_a(iatom_local)%neighbor_list)
          END DO

          CALL allocate_subcell(subcell_b,nsubcell,maxatom_local)

          DO iatom_local=1,natom_local_b
            iatom = atom(ikind)%list_local_b_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
            s(:) = atom(ikind)%s_pbc(:,iatom) + 0.5_dp
            i = INT(s(1)*REAL(nsubcell(1),dp) + 1.0_dp)
            j = INT(s(2)*REAL(nsubcell(2),dp) + 1.0_dp)
            k = INT(s(3)*REAL(nsubcell(3),dp) + 1.0_dp)
            i = MIN(MAX(i,1),nsubcell(1))
            j = MIN(MAX(j,1),nsubcell(2))
            k = MIN(MAX(k,1),nsubcell(3))
            subcell_b(i,j,k)%natom =  subcell_b(i,j,k)%natom + 1
            subcell_b(i,j,k)%atom_list(subcell_b(i,j,k)%natom) = iatom_local
            CALL add_neighbor_list(neighbor_list_set=bc_neighbor_list_set,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   neighbor_list=kind_b(iatom_local)%neighbor_list)
          END DO

          DO katom=1,natom_c

            atom_c = atom(kkind)%list(katom)
            sc_pbc(:) = atom(kkind)%s_pbc(:,katom)

            loop2_kcell: DO kcell=-ncell(3),ncell(3)

              sc(3) = sc_pbc(3) + REAL(kcell,dp)
              sc_min(3) = sc(3) - sac_max(3)
              IF (sc_min(3) >= 0.5_dp) EXIT loop2_kcell
              sc_max(3) = sc(3) + sac_max(3)
              IF (sc_max(3) < -0.5_dp) CYCLE loop2_kcell
              cell_c(3) = kcell

              loop2_jcell: DO jcell=-ncell(2),ncell(2)

                sc(2) = sc_pbc(2) + REAL(jcell,dp)
                sc_min(2) = sc(2) - sac_max(2)
                IF (sc_min(2) >= 0.5_dp) EXIT loop2_jcell
                sc_max(2) = sc(2) + sac_max(2)
                IF (sc_max(2) < -0.5_dp) CYCLE loop2_jcell
                cell_c(2) = jcell

                loop2_icell: DO icell=-ncell(1),ncell(1)

                  sc(1) = sc_pbc(1) + REAL(icell,dp)
                  sc_min(1) = sc(1) - sac_max(1)
                  IF (sc_min(1) >= 0.5_dp) EXIT loop2_icell
                  sc_max(1) = sc(1) + sac_max(1)
                  IF (sc_max(1) < -0.5_dp) CYCLE loop2_icell
                  cell_c(1) = icell

                  rc(:) = scaled_to_real(sc(:),cell)

                  loop1_k: DO k=1,nsubcell(3)
                    loop1_j: DO j=1,nsubcell(2)
                      loop1_i: DO i=1,nsubcell(1)

                        IF (subcell_a(i,j,k)%natom == 0) CYCLE

                        IF (sc_max(3) < subcell_a(i,j,k)%s_min(3)) EXIT loop1_k
                        IF (sc_min(3) >= subcell_a(i,j,k)%s_max(3)) CYCLE loop1_k

                        IF (sc_max(2) < subcell_a(i,j,k)%s_min(2)) EXIT loop1_j
                        IF (sc_min(2) >= subcell_a(i,j,k)%s_max(2)) CYCLE loop1_j

                        IF (sc_max(1) < subcell_a(i,j,k)%s_min(1)) EXIT loop1_i
                        IF (sc_min(1) >= subcell_a(i,j,k)%s_max(1)) CYCLE loop1_i

                        DO iatom_subcell=1,subcell_a(i,j,k)%natom
                          iatom_local = subcell_a(i,j,k)%atom_list(iatom_subcell)
                          iatom = atom(ikind)%list_local_a_index(iatom_local)
                          rac(:) = rc(:) - atom(ikind)%r_pbc(:,iatom)
                          rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                          IF (rac2 < rac2_max) THEN
                            CALL add_neighbor_node(&
                              neighbor_list=kind_a(iatom_local)%neighbor_list,&
                              neighbor=atom_c,&
                              cell=cell_c,&
                              r=rac)
                          END IF
                        END DO

                      END DO loop1_i
                    END DO loop1_j
                  END DO loop1_k

                  loop_kcell_b: DO kcell_b=cell_c(3)-ncell(3),cell_c(3)+ncell(3)

                    s(3) = REAL(kcell_b,dp)
                    s_min(3) = sc_min(3) - s(3)
                    IF (s_min(3) >= 0.5_dp) CYCLE loop_kcell_b
                    s_max(3) = sc_max(3) - s(3)
                    IF (s_max(3) < -0.5_dp) EXIT loop_kcell_b

                    loop_jcell_b: DO jcell_b=cell_c(2)-ncell(2),cell_c(2)+ncell(2)

                      s(2) = REAL(jcell_b,dp)
                      s_min(2) = sc_min(2) - s(2)
                      IF (s_min(2) >= 0.5_dp) CYCLE loop_jcell_b
                      s_max(2) = sc_max(2) - s(2)
                      IF (s_max(2) < -0.5_dp) EXIT loop_jcell_b

                      loop_icell_b: DO icell_b=cell_c(1)-ncell(1),cell_c(1)+ncell(1)

                        s(1) = REAL(icell_b,dp)
                        s_min(1) = sc_min(1) - s(1)
                        IF (s_min(1) >= 0.5_dp) CYCLE loop_icell_b
                        s_max(1) = sc_max(1) - s(1)
                        IF (s_max(1) < -0.5_dp) EXIT loop_icell_b

                        r(:) = rc(:) - scaled_to_real(s(:),cell)

                        loop2_k: DO k=1,nsubcell(3)
                          loop2_j: DO j=1,nsubcell(2)
                            loop2_i: DO i=1,nsubcell(1)

                              IF (subcell_b(i,j,k)%natom == 0) CYCLE

                              IF (s_max(3) < subcell_b(i,j,k)%s_min(3)) EXIT loop2_k
                              IF (s_min(3) >= subcell_b(i,j,k)%s_max(3)) CYCLE loop2_k

                              IF (s_max(2) < subcell_b(i,j,k)%s_min(2)) EXIT loop2_j
                              IF (s_min(2) >= subcell_b(i,j,k)%s_max(2)) CYCLE loop2_j

                              IF (s_max(1) < subcell_b(i,j,k)%s_min(1)) EXIT loop2_i
                              IF (s_min(1) >= subcell_b(i,j,k)%s_max(1)) CYCLE loop2_i

                              DO iatom_subcell=1,subcell_b(i,j,k)%natom
                                iatom_local = subcell_b(i,j,k)%atom_list(iatom_subcell)
                                iatom = atom(ikind)%list_local_b_index(iatom_local)
                                rac(:) = r(:) - atom(ikind)%r_pbc(:,iatom)
                                rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                                IF (rac2 < rac2_max) THEN
                                  CALL add_neighbor_node(&
                                    neighbor_list=kind_b(iatom_local)%neighbor_list,&
                                    neighbor=atom_c,&
                                    cell=cell_c,&
                                    r=rac)
                                END IF
                              END DO

                            END DO loop2_i
                          END DO loop2_j
                        END DO loop2_k

                      END DO loop_icell_b
                    END DO loop_jcell_b
                  END DO loop_kcell_b

                END DO loop2_icell
              END DO loop2_jcell
            END DO loop2_kcell

          END DO

          CALL deallocate_subcell(subcell_a)
          CALL deallocate_subcell(subcell_b)

        END IF

      END DO

    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE build_gth_ppnl_neighbor_lists

! *****************************************************************************

  SUBROUTINE allocate_subcell(subcell,nsubcell,maxatom)

!   Purpose: Allocate and initialize a subcell grid structure for the atomic
!            neighbor search.

!   History: - Creation (12.06.2003,MK)

!   ***************************************************************************

    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell
    INTEGER, DIMENSION(3), INTENT(IN)        :: nsubcell
    INTEGER, INTENT(IN)                      :: maxatom

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_subcell"

    INTEGER                                  :: i, istat, j, k, na, nb, nc
    REAL(dp)                                 :: a_max, a_min, b_max, b_min, &
                                                c_max, c_min, delta_a, &
                                                delta_b, delta_c

!   ---------------------------------------------------------------------------

    na = nsubcell(1)
    nb = nsubcell(2)
    nc = nsubcell(3)

    ALLOCATE (subcell(na,nb,nc),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "subcell",na*nb*nc*int_size)

    delta_a = 1.0_dp/REAL(na,dp)
    delta_b = 1.0_dp/REAL(nb,dp)
    delta_c = 1.0_dp/REAL(nc,dp)

    c_min = -0.5_dp

    DO k=1,nc
      c_max = c_min + delta_c
      b_min = -0.5_dp
      DO j=1,nb
        b_max = b_min + delta_b
        a_min = -0.5_dp
        DO i=1,na
          a_max = a_min + delta_a
          subcell(i,j,k)%s_min(1) = a_min
          subcell(i,j,k)%s_min(2) = b_min
          subcell(i,j,k)%s_min(3) = c_min
          subcell(i,j,k)%s_max(1) = a_max
          subcell(i,j,k)%s_max(2) = b_max
          subcell(i,j,k)%s_max(3) = c_max
          subcell(i,j,k)%natom = 0
          ALLOCATE (subcell(i,j,k)%atom_list(maxatom),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "subcell(i,j,k)%atom_list",&
                                           maxatom*int_size)
          subcell(i,j,k)%atom_list(:) = 0
          a_min = a_max
        END DO
        b_min = b_max
      END DO
      c_min = c_max
    END DO

  END SUBROUTINE allocate_subcell

! *****************************************************************************

  SUBROUTINE deallocate_subcell(subcell)

!   Purpose: Deallocate a subcell grid structure.

!   History: - Creation (16.06.2003,MK)

!   ***************************************************************************

    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_subcell"

    INTEGER                                  :: i, istat, j, k

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(subcell)) THEN

      DO k=1,SIZE(subcell,3)
        DO j=1,SIZE(subcell,2)
          DO i=1,SIZE(subcell,1)
            DEALLOCATE (subcell(i,j,k)%atom_list,STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine_name,module_name,__LINE__,&
                               "subcell(i,j,k)%atom_list")
            END IF
          END DO
        END DO
      END DO

      DEALLOCATE (subcell,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "subcell")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer subcell is not associated")

    END IF

  END SUBROUTINE deallocate_subcell

! *****************************************************************************

  SUBROUTINE write_distribution(ab,atomic_kind_set,globenv)

!   Purpose: Print the distribution of the simple pair neighbor list.

!   History: - Creation (19.06.2003,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_distribution"

    INTEGER :: group, iab, ikind, ilist, ipe, istat, jkind, mype, n, nkind, &
      nlist, nnode, npe, nsgf_a, nsgf_b, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nblock, nelement
    LOGICAL                                  :: ionode
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

!   ---------------------------------------------------------------------------
!   *** Quick return, if no printing is requested ***

    IF (.NOT.globenv%print%distribution) RETURN

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos + 1
    npe = globenv%num_pe
    output_unit = globenv%scr

!   *** Allocate work storage ***

    ALLOCATE (nblock(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nblock",npe*int_size)
    nblock(:) = 0

    ALLOCATE (nelement(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nelement",npe*int_size)
    nelement(:) = 0

    nkind = SIZE(atomic_kind_set)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (ASSOCIATED(orb_basis_set)) THEN
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nsgf_a)
      ELSE
        nsgf_a = 1
      END IF

      DO jkind=1,nkind

        iab = ikind + nkind*(jkind - 1)

        neighbor_list_set => ab(iab)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (ASSOCIATED(orb_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nsgf_b)
        ELSE
          nsgf_b = 1
        END IF

        n = nsgf_a*nsgf_b

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 nnode=nnode)

          nblock(mype) = nblock(mype) + nnode
          nelement(mype) = nelement(mype) + n*nnode

          neighbor_list => next(neighbor_list)

        END DO

      END DO

    END DO

    CALL mp_sum(nblock,group)
    CALL mp_sum(nelement,group)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I6,T27,I10,T55,I10))")&
        "DISTRIBUTION OF THE NEIGHBOR LISTS",&
        "Process   Number of particle pairs   Number of matrix elements",&
        (ipe-1,nblock(ipe),nelement(ipe),ipe=1,npe)
      WRITE (UNIT=output_unit,FMT="(/,T7,A3,T27,I10,T55,I10)")&
        "Sum",SUM(nblock),SUM(nelement)
    END IF

!   *** Release work storage ***

    DEALLOCATE (nblock,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nblock")

    DEALLOCATE (nelement,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nelement")

  END SUBROUTINE write_distribution

! *****************************************************************************

  SUBROUTINE write_neighbor_lists(ab,particle_set,cell,globenv,name)

!   Purpose: Write a set of neighbor lists to the output unit.

!   History: - Creation (04.03.2002,MK)
!            - Adapted to the new parallelized neighbor list version
!              (26.06.2003,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), INTENT(IN)               :: ab
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    CHARACTER(LEN=*), INTENT(IN)             :: name

    CHARACTER(LEN=8)                         :: unit_of_length_name
    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: atom_a, atom_b, group, iab, &
                                                ilist, inode, ipe, mype, &
                                                nlist, nneighbor, nnode, npe, &
                                                output_unit
    INTEGER, DIMENSION(3)                    :: cell_a, cell_b
    LOGICAL                                  :: ionode, print_headline
    REAL(dp)                                 :: dab, unit_of_length
    REAL(dp), DIMENSION(3)                   :: ra, rab, rb
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos + 1
    npe = globenv%num_pe
    output_unit = globenv%scr

    CALL get_cell(cell=cell,&
                  unit_of_length=unit_of_length,&
                  unit_of_length_name=unit_of_length_name)

    CALL mp_sync(group)
    CALL m_flush(output_unit)
    CALL mp_sync(group)

!   *** Create the headline ***

    string = ""

    IF (npe == 1) THEN
      string = TRIM(name)//" IN "//TRIM(unit_of_length_name)
    ELSE
      WRITE (UNIT=string,FMT="(A,I5,A)")&
        TRIM(name)//" IN "//TRIM(unit_of_length_name)//" (PROCESS",mype-1,")"
        CALL compress(string)
    END IF

    print_headline = .TRUE.

!   *** All processes print now their neighbor lists portion consecutively ***

    nneighbor = 0

    DO ipe=1,npe

      IF (ipe == mype) THEN

!       *** Print first part of headline ***

        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(string)

        DO iab=1,SIZE(ab)

          neighbor_list_set => ab(iab)%neighbor_list_set

          IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

!         *** Loop over all atoms and their corresponding neighbor lists ***

          CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                     nlist=nlist)

          neighbor_list => first_list(neighbor_list_set)

          DO ilist=1,nlist

            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a,&
                                   cell=cell_a,&
                                   nnode=nnode)
            nneighbor = nneighbor + nnode

            IF (globenv%print%level > LOW) THEN

!             *** Print second part of headline ***

              IF (print_headline) THEN
                WRITE (UNIT=output_unit,FMT="(/,T3,A,7X,A,2(11X,A),10X,A)")&
                  "Atom  Neighbors  Cell(i,j,k)","X","Y","Z","Distance"
                print_headline = .FALSE.
              END IF

              ra(:) = pbc(particle_set(atom_a)%r,cell,cell_a)

              WRITE (UNIT=output_unit,FMT="(/,T2,I5,3X,I6,3X,3I4,3F12.6)")&
                atom_a,nnode,cell_a(1:3),ra(1:3)/unit_of_length

!             *** Direct the work pointer to the  ***
!             *** start point of the current list ***

              neighbor_node => first_node(neighbor_list)

!             *** Traverse the neighbor list of the current ***
!             *** atom and print the stored information     ***

              IF (globenv%print%level > MEDIUM) THEN

                DO inode=1,nnode

                  CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                         neighbor=atom_b,&
                                         cell=cell_b,&
                                         r=rab)

                  rb(:) = ra(:) + rab(:)
                  dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

                  WRITE (UNIT=output_unit,&
                         FMT="(T10,I6,3X,3I4,3F12.6,2X,F12.6)")&
                    atom_b,cell_b(1:3),rb(1:3)/unit_of_length,dab/unit_of_length

                  neighbor_node => next(neighbor_node)

                END DO

              END IF

            END IF

            neighbor_list => next(neighbor_list)

          END DO

        END DO

        string = ""

        IF (npe == 1) THEN
          WRITE (UNIT=string,FMT="(A,I12)")&
            "Total number of neighbor interactions:",nneighbor
        ELSE
          WRITE (UNIT=string,FMT="(A,I12,A,I12)")&
            "Total number of neighbor interactions for process",mype-1,":",&
            nneighbor
        END IF
        CALL compress(string)
        WRITE (UNIT=output_unit,FMT="(/,T2,A)") TRIM(string)

      END IF

      CALL mp_sync(group)
      CALL m_flush(output_unit)
      CALL mp_sync(group)

    END DO

  END SUBROUTINE write_neighbor_lists

! *****************************************************************************

END MODULE qs_neighbor_lists
