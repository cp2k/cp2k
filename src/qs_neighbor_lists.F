!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Generate the atomic neighbor lists.
!> \par History
!>      - List rebuild for sab_orb neighbor list (10.09.2002,MK)
!>      - List rebuild for all lists (25.09.2002,MK)
!>      - Row-wise parallelized version (16.06.2003,MK)
!>      - Row- and column-wise parallelized version (19.07.2003,MK)
!>      - bug fix for non-periodic case (23.02.06,MK)
!>      - major refactoring (25.07.10,jhu)
!> \author Matthias Krack (08.10.1999,26.03.2002,16.06.2003)
! *****************************************************************************
MODULE qs_neighbor_lists
  USE almo_scf_types,                  ONLY: almo_max_cutoff_multiplier
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             plane_distance,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE input_constants,                 ONLY: dispersion_uff,&
                                             do_method_lrigpw,&
                                             do_se_IS_slater,&
                                             vdw_pairpot_dftd3,&
                                             vdw_pairpot_dftd3bj,&
                                             xc_vdw_fun_pairpot
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE kpoint_types,                    ONLY: kpoint_type
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: bohr
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_dispersion_types,             ONLY: qs_dispersion_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_env_type,&
                                             set_ks_env
  USE qs_neighbor_list_types,          ONLY: &
       add_neighbor_list, add_neighbor_node, allocate_neighbor_list_set, &
       deallocate_neighbor_list_set, get_iterator_info, &
       neighbor_list_iterate, neighbor_list_iterator_create, &
       neighbor_list_iterator_p_type, neighbor_list_iterator_release, &
       neighbor_list_p_type, neighbor_list_set_p_type, neighbor_list_set_type
  USE scptb_types,                     ONLY: get_scptb_parameter,&
                                             scptb_parameter_type
  USE string_utilities,                ONLY: compress
  USE subcell_types,                   ONLY: allocate_subcell,&
                                             deallocate_subcell,&
                                             give_ijk_subcell,&
                                             subcell_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate,&
                                             sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE local_atoms_type
    INTEGER, DIMENSION(:), POINTER                   :: list,&
                                                        list_local_a_index,&
                                                        list_local_b_index,&
                                                        list_1d,&
                                                        list_a_mol,&
                                                        list_b_mol
  END TYPE local_atoms_type
! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_neighbor_lists'

  ! private counter, used to version qs neighbor lists
  INTEGER, SAVE, PRIVATE :: last_qs_neighbor_list_id_nr=0

  ! Public subroutines
  PUBLIC :: build_qs_neighbor_lists, local_atoms_type, atom2d_cleanup, &
            atom2d_build, build_neighbor_lists, pair_radius_setup, &
            write_neighbor_lists
CONTAINS

! *****************************************************************************
!> \brief   free the internals of atom2d
!> \param atom2d ...
!> \param error ...
!> \param
! *****************************************************************************
  SUBROUTINE atom2d_cleanup(atom2d,error)
    TYPE(local_atoms_type), DIMENSION(:)     :: atom2d
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom2d_cleanup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
    DO ikind=1,SIZE(atom2d)
       NULLIFY (atom2d(ikind)%list)
       IF (ASSOCIATED(atom2d(ikind)%list_local_a_index)) THEN
          DEALLOCATE (atom2d(ikind)%list_local_a_index,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(atom2d(ikind)%list_local_b_index)) THEN
          DEALLOCATE (atom2d(ikind)%list_local_b_index,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(atom2d(ikind)%list_a_mol)) THEN
          DEALLOCATE (atom2d(ikind)%list_a_mol,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(atom2d(ikind)%list_b_mol)) THEN
          DEALLOCATE (atom2d(ikind)%list_b_mol,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(atom2d(ikind)%list_1d)) THEN
          DEALLOCATE (atom2d(ikind)%list_1d,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO
    CALL timestop(handle)

  END SUBROUTINE atom2d_cleanup

! *****************************************************************************
!> \brief   Build some distribution structure of atoms, refactored from build_qs_neighbor_lists
!> \param atom2d output
!> \param orb_radius output
!> \param orb_present output
!> \param distribution_1d ...
!> \param distribution_2d ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param molecule_set ...
!> \param molecule_only ...
!> \param kg ...
!> \param dftb ...
!> \param particle_set ...
!> \param error ...
!> \author  JH
! *****************************************************************************
  SUBROUTINE atom2d_build(atom2d,orb_radius,orb_present,distribution_1d,distribution_2d,&
                          atomic_kind_set,qs_kind_set,molecule_set,molecule_only,kg,dftb,particle_set,error)
    TYPE(local_atoms_type), DIMENSION(:)     :: atom2d
    REAL(dp), DIMENSION(:)                   :: orb_radius
    LOGICAL, DIMENSION(:)                    :: orb_present
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    LOGICAL                                  :: molecule_only
    LOGICAL, OPTIONAL                        :: kg
    LOGICAL                                  :: dftb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom2d_build', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, handle, ia, iat, iatom, iatom_local, ikind, imol, &
      natom, natom_a, natom_local_a, natom_local_b, nel, nkind, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom2mol, listindex, listsort
    INTEGER, DIMENSION(:), POINTER           :: atom_of_kind, &
                                                local_cols_array, &
                                                local_rows_array
    LOGICAL                                  :: do_kg, failure
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_atom

    CALL timeset(routineN,handle)
    failure=.FALSE.

    do_kg=.FALSE.
    IF(PRESENT(kg)) do_kg=kg

    nkind=SIZE(atomic_kind_set)
    natom=SIZE(particle_set)
    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind)

    IF ( molecule_only ) THEN
       ALLOCATE(atom2mol(natom))            
       DO imol=1,SIZE(molecule_set)
          DO iat=molecule_set(imol)%first_atom,molecule_set(imol)%last_atom
             atom2mol(iat)=imol
          ENDDO
       END DO
    ENDIF

    DO ikind=1,nkind
       NULLIFY (atom2d(ikind)%list)
       NULLIFY (atom2d(ikind)%list_local_a_index)
       NULLIFY (atom2d(ikind)%list_local_b_index)
       NULLIFY (atom2d(ikind)%list_1d)
       NULLIFY (atom2d(ikind)%list_a_mol)
       NULLIFY (atom2d(ikind)%list_b_mol)

       CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom2d(ikind)%list)
       CALL get_qs_kind(qs_kind_set(ikind), orb_basis_set=orb_basis_set)

       natom_a = SIZE(atom2d(ikind)%list)


       natom_local_a = distribution_2d%n_local_rows(ikind)
       natom_local_b = distribution_2d%n_local_cols(ikind)
       local_rows_array => distribution_2d%local_rows(ikind)%array
       local_cols_array => distribution_2d%local_cols(ikind)%array

       nel = distribution_1d%n_el(ikind)
       ALLOCATE (atom2d(ikind)%list_1d(nel),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iat=1,nel
         ia = distribution_1d%list(ikind)%array(iat)
         atom2d(ikind)%list_1d(iat) = atom_of_kind(ia)
       END DO

       IF (ASSOCIATED(orb_basis_set)) THEN
          orb_present(ikind) = .TRUE.
          IF (do_kg) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,short_kind_radius=orb_radius(ikind))
          ELSE
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,kind_radius=orb_radius(ikind))
          END IF
       ELSE IF (dftb) THEN
          ! Set the interaction radius for the neighbor lists (DFTB case)
          ! This includes all interactions (orbitals and short range pair potential) except vdW
          CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom)
          CALL get_dftb_atom_param(dftb_parameter=dftb_atom,&
                                   cutoff=orb_radius(ikind),&
                                   defined=orb_present(ikind),error=error)
       ELSE
          orb_present(ikind) = .FALSE.
       END IF

       IF ( orb_present(ikind) ) THEN
          ALLOCATE (listsort(natom_a),listindex(natom_a),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          listsort(1:natom_a) = atom2d(ikind)%list(1:natom_a)
          CALL sort(listsort,natom_a,listindex)
          ! Block rows
          IF (natom_local_a > 0) THEN
             ALLOCATE (atom2d(ikind)%list_local_a_index(natom_local_a),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE (atom2d(ikind)%list_a_mol(natom_local_a),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             atom2d(ikind)%list_a_mol(:) = 0


             ! Build index vector for mapping
             DO iatom_local=1,natom_local_a
                atom_a = local_rows_array(iatom_local)
                iatom = locate(listsort,atom_a)
                atom2d(ikind)%list_local_a_index(iatom_local) = listindex(iatom)
                IF ( molecule_only ) atom2d(ikind)%list_a_mol(iatom_local) = atom2mol(atom_a)
             END DO

          END IF

          ! Block columns
          IF (natom_local_b > 0) THEN

             ALLOCATE (atom2d(ikind)%list_local_b_index(natom_local_b),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE (atom2d(ikind)%list_b_mol(natom_local_b),STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             atom2d(ikind)%list_b_mol(:) = 0

             ! Build index vector for mapping
             DO iatom_local=1,natom_local_b
                atom_a = local_cols_array(iatom_local)
                iatom = locate(listsort,atom_a)
                atom2d(ikind)%list_local_b_index(iatom_local) = listindex(iatom)
                IF (molecule_only) atom2d(ikind)%list_b_mol(iatom_local) = atom2mol(atom_a)
             END DO

          END IF

          DEALLOCATE (listsort,listindex,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       END IF

    ENDDO

    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE atom2d_build

! *****************************************************************************
!> \brief   Build all the required neighbor lists for Quickstep.
!> \param qs_env ...
!> \param para_env ...
!> \param molecular ...
!> \param force_env_section ...
!> \param error ...
!> \date    28.08.2000
!> \par History
!>          - Major refactoring (25.07.2010,jhu)
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_qs_neighbor_lists(qs_env,para_env,molecular,force_env_section,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, OPTIONAL                        :: molecular
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_qs_neighbor_lists', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: print_key_path
    INTEGER                                  :: handle, ikind, iw, maxatom, &
                                                nkind, stat, zat
    LOGICAL :: all_potential_present, almo, dftb, dokp, failure, &
      gth_potential_present, lri_optbas, lrigpw, mic, molecule_only, nddo, &
      paw_atom, paw_atom_present, scptb
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: all_present, aux_fit_present, &
      aux_present, core_present, default_present, oce_present, orb_present, &
      ppl_present, ppnl_present
    REAL(dp)                                 :: almo_rcov, almo_rvdw, alpha, &
                                                subcells
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: all_pot_rad, aux_fit_radius, &
      c_radius, calpha, core_radius, oce_radius, orb_radius, ppl_radius, &
      ppnl_radius, zeff
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: pair_radius
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set, &
                                                aux_fit_basis_set, &
                                                orb_basis_set
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom2d
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_all, &
      sab_almo, sab_aux_fit, sab_aux_fit_asymm, sab_aux_fit_vs_orb, sab_cn, &
      sab_core, sab_lrc, sab_orb, sab_scp, sab_se, sab_tbe, sab_vdw, sac_ae, &
      sac_ppl, sap_oce, sap_ppnl, soo_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(paw_proj_set_type), POINTER         :: paw_proj
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_atom
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(scptb_parameter_type), POINTER      :: scptb_kind
    TYPE(section_vals_type), POINTER         :: neighbor_list_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY (atomic_kind_set, qs_kind_set, cell, neighbor_list_section,&
             distribution_1d, distribution_2d, gth_potential, orb_basis_set,&
             particle_set, molecule_set, dft_control,ks_env)

    NULLIFY (sab_orb)
    NULLIFY (sac_ae)
    NULLIFY (sac_ppl)
    NULLIFY (sap_ppnl)
    NULLIFY (sap_oce)
    NULLIFY (sab_se)
    NULLIFY (sab_lrc)
    NULLIFY (sab_tbe)
    NULLIFY (sab_core)
    NULLIFY (sab_all)
    NULLIFY (sab_vdw)
    NULLIFY (sab_cn)
    NULLIFY (sab_aux_fit)
    NULLIFY (sab_aux_fit_vs_orb)
    NULLIFY (sab_scp)
    NULLIFY (soo_list)
    NULLIFY (sab_almo)

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    cell=cell,&
                    kpoints=kpoints,&
                    distribution_2d=distribution_2d,&
                    local_particles=distribution_1d,&
                    particle_set=particle_set,&
                    molecule_set=molecule_set,&
                    dft_control=dft_control,&
                    error=error)

    neighbor_list_section => section_vals_get_subs_vals(force_env_section,"DFT%PRINT%NEIGHBOR_LISTS",error=error)

    molecule_only = .FALSE.
    IF (PRESENT(molecular)) molecule_only = molecular

    ! This sets the id number of the qs neighbor lists, new lists, means new version
    ! new version implies new sparsity of the matrices

    last_qs_neighbor_list_id_nr=last_qs_neighbor_list_id_nr+1
    CALL set_ks_env(ks_env=ks_env,neighbor_list_id=last_qs_neighbor_list_id_nr,error=error)

    CALL get_ks_env(ks_env=ks_env,&
                    sab_orb=sab_orb,&
                    sab_aux_fit=sab_aux_fit,&
                    sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                    sab_aux_fit_asymm=sab_aux_fit_asymm,&
                    sac_ae=sac_ae,&
                    sac_ppl=sac_ppl,&
                    sab_vdw=sab_vdw,&
                    sap_ppnl=sap_ppnl,&
                    sap_oce=sap_oce,&
                    sab_se=sab_se,&
                    sab_lrc=sab_lrc,&
                    sab_tbe=sab_tbe,&
                    sab_core=sab_core,&
                    sab_all=sab_all, &
                    sab_scp=sab_scp, &
                    sab_almo=sab_almo, &
                    error=error)

    dokp = (kpoints%nkp > 0)
    nddo = dft_control%qs_control%semi_empirical
    dftb = dft_control%qs_control%dftb
    scptb = dft_control%qs_control%scptb
    almo = dft_control%qs_control%do_almo_scf
    lrigpw = (dft_control%qs_control%method_id == do_method_lrigpw)
    lri_optbas = dft_control%qs_control%lri_optbas

    CALL get_atomic_kind_set(atomic_kind_set, maxatom=maxatom)
    CALL get_qs_kind_set(qs_kind_set, paw_atom_present=paw_atom_present,&
                         gth_potential_present=gth_potential_present,&
                         all_potential_present=all_potential_present)

    CALL section_vals_val_get(qs_env%input,"DFT%SUBCELLS",r_val=subcells,error=error)

    ! Allocate work storage

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (orb_present(nkind),aux_fit_present(nkind),aux_present(nkind),&
              default_present(nkind),core_present(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (orb_radius(nkind),aux_fit_radius(nkind),c_radius(nkind),&
              core_radius(nkind),calpha(nkind),zeff(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    orb_radius(:) = 0.0_dp
    aux_fit_radius(:) = 0.0_dp
    c_radius(:) = 0.0_dp
    core_radius(:) = 0.0_dp
    calpha(:) = 0.0_dp
    zeff(:) = 0.0_dp

    ALLOCATE (pair_radius(nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (gth_potential_present) THEN
       ALLOCATE (ppl_present(nkind),ppl_radius(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ppl_radius = 0.0_dp
       ALLOCATE (ppnl_present(nkind),ppnl_radius(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ppnl_radius = 0.0_dp
    END IF

    IF (paw_atom_present) THEN
       ALLOCATE (oce_present(nkind),oce_radius(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       oce_radius = 0.0_dp
    END IF

    IF (all_potential_present) THEN
       ALLOCATE (all_present(nkind),all_pot_rad(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       all_pot_rad = 0.0_dp
    END IF

    ! Initialize the local data structures
    ALLOCATE (atom2d(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL atom2d_build(atom2d,orb_radius,orb_present,distribution_1d,distribution_2d,&
                      atomic_kind_set,qs_kind_set,molecule_set,molecule_only,&
                      dftb=dftb,particle_set=particle_set,error=error)

    DO ikind=1,nkind

       CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom2d(ikind)%list)

       CALL get_qs_kind(qs_kind_set(ikind),&
                        orb_basis_set=orb_basis_set,&
                        aux_basis_set=aux_basis_set,&
                        aux_fit_basis_set=aux_fit_basis_set,&
                        paw_proj_set=paw_proj,&
                        paw_atom=paw_atom,&
                        all_potential=all_potential,&
                        gth_potential=gth_potential)

       ! SCP
       IF (ASSOCIATED(aux_basis_set)) THEN
          aux_present(ikind) = .TRUE.
       ELSE
          aux_present(ikind) = .FALSE.
       ENDIF

       IF (ASSOCIATED(aux_fit_basis_set)) THEN
          aux_fit_present(ikind) = .TRUE.
          CALL get_gto_basis_set(gto_basis_set=aux_fit_basis_set,kind_radius=aux_fit_radius(ikind))
        ELSE
          aux_fit_present(ikind) = .FALSE.
       END IF

       ! core overlap
       CALL get_qs_kind(qs_kind_set(ikind),&
                        alpha_core_charge=calpha(ikind),&
                        core_charge_radius=core_radius(ikind),&
                        zeff=zeff(ikind))
       IF(zeff(ikind) /= 0._dp .AND. calpha(ikind) /= 0._dp) THEN
         core_present(ikind) = .TRUE.
       ELSE
         core_present(ikind) = .FALSE.
       END IF

       ! Pseudopotentials
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

       ! GAPW
       IF (paw_atom_present) THEN
          IF (paw_atom) THEN
             oce_present(ikind) = .TRUE.
             CALL get_paw_proj_set(paw_proj_set=paw_proj,rcprj=oce_radius(ikind))
          ELSE
             oce_present(ikind) = .FALSE.
          END IF
       END IF

       ! Check the presence of an all electron potential
       IF (all_potential_present) THEN
          IF (ASSOCIATED(all_potential)) THEN
             all_present(ikind) = .TRUE.
             CALL get_potential(potential=all_potential,core_charge_radius=all_pot_rad(ikind))
          ELSE
             all_present(ikind) = .FALSE.
             all_pot_rad(ikind) = 0.0_dp
          END IF
       END IF

    END DO

    IF (dokp) THEN
       ! no MIC for kpoints
       mic = .FALSE.
    ELSE
       ! enforce MIC for interaction lists in SE and DFTB
       IF (nddo .OR. dftb .OR. scptb) THEN
          mic = .TRUE.
       ELSE
          mic = .FALSE.
       END IF
    END IF

    ! Build the orbital-orbital overlap neighbor lists
    CALL pair_radius_setup(orb_present,orb_present,orb_radius,orb_radius,pair_radius,error)
    CALL build_neighbor_lists(sab_orb,particle_set,atom2d,cell,pair_radius,&
                              mic=mic,subcells=subcells,molecular=molecule_only,name="sab_orb",error=error)
    CALL set_ks_env(ks_env=ks_env,sab_orb=sab_orb,error=error)
    CALL write_neighbor_lists(sab_orb,particle_set,cell,para_env,neighbor_list_section,&
         "/SAB_ORB","sab_orb","ORBITAL ORBITAL",error)

    ! Build orbital-orbital list containing all the pairs, to be used with
    ! non-symmetric operators. Beware: the cutoff of the orbital-orbital overlap
    ! might be not optimal. It should be verified for each operator.
    ! right now only used for linear response (e.g. nmr). Should be computed more sparingly
    IF (.NOT. (nddo .OR. dftb .OR. scptb)) THEN
       CALL build_neighbor_lists(sab_all,particle_set,atom2d,cell,pair_radius,&
            mic=mic,symmetric=.FALSE.,subcells=subcells,molecular=molecule_only,name="sab_all",error=error)
       CALL set_ks_env(ks_env=ks_env,sab_all=sab_all,error=error)
    ENDIF

    ! Build the core-core overlap neighbor lists
    IF (.NOT. (nddo .OR. dftb .OR. scptb)) THEN
       CALL pair_radius_setup(core_present,core_present,core_radius,core_radius,pair_radius,error)
       CALL build_neighbor_lists(sab_core,particle_set,atom2d,cell,pair_radius,mic=mic,subcells=subcells,&
                                 operator_type="PP",name="sab_core",error=error)
       CALL set_ks_env(ks_env=ks_env,sab_core=sab_core,error=error)
       CALL write_neighbor_lists(sab_core,particle_set,cell,para_env,neighbor_list_section,&
            "/SAB_CORE","sab_core","CORE CORE",error)
    ENDIF


    IF(dft_control%do_admm) THEN
       CALL pair_radius_setup(aux_fit_present,aux_fit_present,aux_fit_radius,aux_fit_radius,pair_radius,error)
       CALL build_neighbor_lists(sab_aux_fit,particle_set,atom2d,cell,pair_radius,&
                                 mic=mic,molecular=molecule_only,subcells=subcells,name="sab_aux_fit",error=error)
       CALL build_neighbor_lists(sab_aux_fit_asymm,particle_set,atom2d,cell,pair_radius,&
                                 mic=mic,symmetric=.FALSE.,molecular=molecule_only,subcells=subcells,&
                                 name="sab_aux_fit_asymm",error=error)
       CALL pair_radius_setup(aux_fit_present,orb_present,aux_fit_radius,orb_radius,pair_radius,error)
       CALL build_neighbor_lists(sab_aux_fit_vs_orb,particle_set,atom2d,cell,pair_radius,&
                                 mic=mic,symmetric=.FALSE.,molecular=molecule_only,subcells=subcells,&
                                 name="sab_aux_fit_vs_orb",error=error)

       CALL set_ks_env(ks_env=ks_env,sab_aux_fit=sab_aux_fit,error=error)
       CALL set_ks_env(ks_env=ks_env,sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,error=error)
       CALL set_ks_env(ks_env=ks_env,sab_aux_fit_asymm=sab_aux_fit_asymm,error=error)

       CALL write_neighbor_lists(sab_aux_fit,particle_set,cell,para_env,neighbor_list_section,&
            "/SAB_AUX_FIT","sab_aux_fit","AUX_FIT_ORBITAL AUX_FIT_ORBITAL",error)
       CALL write_neighbor_lists(sab_aux_fit_vs_orb,particle_set,cell,para_env,neighbor_list_section,&
            "/SAB_AUX_FIT_VS_ORB","sab_aux_fit_vs_orb","ORBITAL AUX_FIT_ORBITAL",error)
    END IF

    ! Build orbital GTH-PPL operator overlap list
    IF (gth_potential_present) THEN
       IF (ANY(ppl_present)) THEN
          CALL pair_radius_setup(orb_present,ppl_present,orb_radius,ppl_radius,pair_radius,error)
          CALL build_neighbor_lists(sac_ppl,particle_set,atom2d,cell,pair_radius,&
               subcells=subcells,operator_type="ABC",name="sac_ppl",error=error)
          CALL set_ks_env(ks_env=ks_env,sac_ppl=sac_ppl,error=error)
          CALL write_neighbor_lists(sac_ppl,particle_set,cell,para_env,neighbor_list_section,&
               "/SAC_PPL","sac_ppl","ORBITAL GTH-PPL",error)
       END IF

       IF (ANY(ppnl_present)) THEN
          CALL pair_radius_setup(orb_present,ppnl_present,orb_radius,ppnl_radius,pair_radius,error)
          CALL build_neighbor_lists(sap_ppnl,particle_set,atom2d,cell,pair_radius,&
               subcells=subcells,operator_type="ABBA",name="sap_ppnl",error=error)
          CALL set_ks_env(ks_env=ks_env,sap_ppnl=sap_ppnl,error=error)
          CALL write_neighbor_lists(sap_ppnl,particle_set,cell,para_env,neighbor_list_section,&
               "/SAP_PPNL","sap_ppnl","ORBITAL GTH-PPNL",error)
       END IF
    END IF


    IF (paw_atom_present) THEN
       ! Build orbital-GAPW projector overlap list
       IF (ANY(oce_present)) THEN
          CALL pair_radius_setup(orb_present,oce_present,orb_radius,oce_radius,pair_radius,error)
          CALL build_neighbor_lists(sap_oce,particle_set,atom2d,cell,pair_radius,&
               subcells=subcells,operator_type="ABBA",name="sap_oce",error=error)
          CALL set_ks_env(ks_env=ks_env,sap_oce=sap_oce,error=error)
          CALL write_neighbor_lists(sap_oce,particle_set,cell,para_env,neighbor_list_section,&
               "/SAP_OCE","sap_oce","ORBITAL(A) PAW-PRJ",error)
       END IF
    END IF

    ! Build orbital-ERFC potential list
    IF (all_potential_present .AND. .NOT. (nddo .OR. dftb .OR. scptb)) THEN
       CALL pair_radius_setup(orb_present,all_present,orb_radius,all_pot_rad,pair_radius,error)
       CALL build_neighbor_lists(sac_ae,particle_set,atom2d,cell,pair_radius,&
            subcells=subcells,operator_type="ABC",name="sac_ae",error=error)
       CALL set_ks_env(ks_env=ks_env,sac_ae=sac_ae,error=error)
       CALL write_neighbor_lists(sac_ae,particle_set,cell,para_env,neighbor_list_section,&
               "/SAC_AE","sac_ae","ORBITAL ERFC POTENTIAL",error)
    END IF

    IF (nddo) THEN
       ! Semi-empirical neighbor lists
       default_present = .TRUE.
       c_radius = dft_control%qs_control%se_control%cutoff_cou
       ! Build the neighbor lists for the Hartree terms
       CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
       IF (dft_control%qs_control%se_control%do_ewald_gks)  THEN
         ! Use MIC for the periodic code of GKS
         CALL build_neighbor_lists(sab_se,particle_set,atom2d,cell,pair_radius,mic=mic,&
              subcells=subcells,name="sab_se",error=error)
       ELSE
         CALL build_neighbor_lists(sab_se,particle_set,atom2d,cell,pair_radius,&
              subcells=subcells,name="sab_se",error=error)
       END IF
       CALL set_ks_env(ks_env=ks_env,sab_se=sab_se,error=error)
       CALL write_neighbor_lists(sab_se,particle_set,cell,para_env,neighbor_list_section,&
            "/SAB_SE","sab_se","HARTREE INTERACTIONS",error)

       ! If requested build the SE long-range correction neighbor list
       IF ((dft_control%qs_control%se_control%do_ewald).AND.&
           (dft_control%qs_control%se_control%integral_screening/=do_se_IS_slater)) THEN
          c_radius = dft_control%qs_control%se_control%cutoff_lrc
          CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
          CALL build_neighbor_lists(sab_lrc,particle_set,atom2d,cell,pair_radius,&
               subcells=subcells,name="sab_lrc",error=error)
          CALL set_ks_env(ks_env=ks_env,sab_lrc=sab_lrc,error=error)
          CALL write_neighbor_lists(sab_lrc,particle_set,cell,para_env,neighbor_list_section,&
               "/SAB_LRC","sab_lrc","SE LONG-RANGE CORRECTION",error)
       END IF
    END IF

    IF (dftb) THEN
       ! Build the neighbor lists for the DFTB Ewald methods
       IF (dft_control%qs_control%dftb_control%do_ewald) THEN
          CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,error=error)
          CALL ewald_env_get ( ewald_env, alpha=alpha )
          c_radius = 0.5_dp*SQRT(-LOG(3.5_dp*alpha**3*1.e-12_dp))/alpha
          CALL pair_radius_setup(orb_present,orb_present,c_radius,c_radius,pair_radius,error)
          CALL build_neighbor_lists(sab_tbe,particle_set,atom2d,cell,pair_radius,mic=mic,&
               subcells=subcells,name="sab_tbe",error=error)
          CALL set_ks_env(ks_env=ks_env,sab_tbe=sab_tbe,error=error)
       END IF

       ! Build the neighbor lists for the DFTB vdW pair potential
       IF (dft_control%qs_control%dftb_control%dispersion) THEN
          IF (dft_control%qs_control%dftb_control%dispersion_type == dispersion_uff) THEN
             DO ikind = 1, nkind
                CALL get_qs_kind(qs_kind_set(ikind),dftb_parameter=dftb_atom)
                CALL get_dftb_atom_param(dftb_parameter=dftb_atom,&
                                         rcdisp=c_radius(ikind),error=error)
             END DO
             default_present=.TRUE.
             CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
             CALL build_neighbor_lists(sab_vdw,particle_set,atom2d,cell,pair_radius,&
                  subcells=subcells,name="sab_vdw",error=error)
             CALL set_ks_env(ks_env=ks_env,sab_vdw=sab_vdw,error=error)
          END IF
       END IF
    END IF

    ! Build the neighbor lists for the vdW pair potential
    CALL get_qs_env(qs_env=qs_env,dispersion_env=dispersion_env,error=error)
    sab_vdw => dispersion_env%sab_vdw
    sab_cn => dispersion_env%sab_cn
    IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN
       c_radius(:) = dispersion_env%rc_disp
       default_present=.TRUE. !include all atoms in vdW (even without basis)
       CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
       CALL build_neighbor_lists(sab_vdw,particle_set,atom2d,cell,pair_radius,&
            subcells=subcells,operator_type="PP",name="sab_vdw",error=error)
       dispersion_env%sab_vdw => sab_vdw

       IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 .OR. &
            dispersion_env%pp_type == vdw_pairpot_dftd3bj ) THEN
          ! Build the neighbor lists for coordination numbers as needed by the DFT-D3 method
          DO ikind = 1, nkind
             CALL get_atomic_kind(atomic_kind_set(ikind), z=zat)
             c_radius(ikind) = 4._dp*ptable(zat)%covalent_radius*bohr
          END DO
          CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
          CALL build_neighbor_lists(sab_cn,particle_set,atom2d,cell,pair_radius,&
               subcells=subcells,operator_type="PP",name="sab_cn",error=error)
          dispersion_env%sab_cn => sab_cn
       END IF
    END IF

    IF (scptb) THEN
       ! Build the neighbor lists for the Core term
       default_present = .FALSE.
       DO ikind = 1, nkind
          CALL get_qs_kind(qs_kind_set(ikind),scptb_parameter=scptb_kind)
          CALL get_scptb_parameter(scptb_kind,rcpair=c_radius(ikind),defined=default_present(ikind))
       END DO
       CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
       CALL build_neighbor_lists(sab_core,particle_set,atom2d,cell,pair_radius,mic=mic,&
            subcells=subcells,name="sab_core",error=error)
       CALL set_ks_env(ks_env=ks_env,sab_core=sab_core,error=error)
       CALL write_neighbor_lists(sab_core,particle_set,cell,para_env,neighbor_list_section,&
            "/SAB_CORE","sab_core","CORE INTERACTIONS",error)
       ! Build the neighbor lists for the SCPTB Ewald methods
       IF (dft_control%qs_control%scptb_control%do_ewald) THEN
          ! Ewald with alpha
          CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,error=error)
          CALL ewald_env_get ( ewald_env, alpha=alpha )
          c_radius = 0.5_dp*SQRT(-LOG(3.5_dp*alpha**3*1.e-12_dp))/alpha
          CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
          CALL build_neighbor_lists(sab_tbe,particle_set,atom2d,cell,pair_radius,mic=mic,&
               subcells=subcells,name="sab_tbe",error=error)
          CALL set_ks_env(ks_env=ks_env,sab_tbe=sab_tbe,error=error)
          ! Correction term for Gaussians
          DO ikind = 1, nkind
             CALL get_qs_kind(qs_kind_set(ikind),scptb_parameter=scptb_kind)
             CALL get_scptb_parameter(scptb_kind,ag=alpha)
             alpha = SQRT(alpha)
             c_radius(ikind) = 0.5_dp*SQRT(-LOG(3.5_dp*alpha**3*1.e-12_dp))/alpha
          END DO
          CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
          CALL build_neighbor_lists(sab_scp,particle_set,atom2d,cell,pair_radius,mic=mic,&
               subcells=subcells,name="sab_scp",error=error)
          CALL set_ks_env(ks_env=ks_env,sab_scp=sab_scp,error=error)
       END IF
    END IF
    IF(lrigpw .OR. lri_optbas) THEN
       ! set neighborlists in lri_env environment
       CALL pair_radius_setup(orb_present,orb_present,orb_radius,orb_radius,pair_radius,error)
       soo_list => qs_env%lri_env%soo_list
       CALL build_neighbor_lists(soo_list,particle_set,atom2d,cell,pair_radius,&
            subcells=subcells,name="soo_list",error=error)
       qs_env%lri_env%soo_list => soo_list
       CALL write_neighbor_lists(soo_list,particle_set,cell,para_env,neighbor_list_section,&
            "/SOO_LIST","soo_list","ORBITAL ORBITAL (RI)",error)
    END IF

    ! Build the neighbor lists for the ALMO delocalization
    IF ( almo ) THEN
       DO ikind = 1, nkind
          CALL get_atomic_kind(atomic_kind_set(ikind),rcov=almo_rcov,rvdw=almo_rvdw)
          ! multiply the radius by some hard-coded number
          c_radius(ikind) = MAX(almo_rcov,almo_rvdw)*bohr*&
                            almo_max_cutoff_multiplier
       END DO
       default_present=.TRUE. !include all atoms (even without basis)
       CALL pair_radius_setup(default_present,default_present,c_radius,c_radius,pair_radius,error)
       CALL build_neighbor_lists(sab_almo,particle_set,atom2d,cell,pair_radius,&
            subcells=subcells,operator_type="PP",name="sab_almo",error=error)
       CALL set_ks_env(ks_env=ks_env,sab_almo=sab_almo,error=error)
    ENDIF

    ! Print particle distribution
    print_key_path = "PRINT%DISTRIBUTION"
    IF (BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
         print_key_path,error=error),&
         cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger=logger,&
                                 basis_section=force_env_section,&
                                 print_key_path=print_key_path,&
                                 extension=".out",&
                                 error=error)
       CALL write_neighbor_distribution(sab_orb,qs_kind_set,iw,para_env,error)
       CALL cp_print_key_finished_output(unit_nr=iw,&
                                         logger=logger,&
                                         basis_section=force_env_section,&
                                         print_key_path=print_key_path,&
                                         error=error)
    END IF

    ! Release work storage
    CALL atom2d_cleanup(atom2d,error)

    DEALLOCATE (atom2d,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (orb_present,default_present,core_present,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (orb_radius,aux_fit_radius,c_radius,core_radius,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (calpha,zeff,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (pair_radius,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (gth_potential_present) THEN
       DEALLOCATE (ppl_present,ppl_radius,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (ppnl_present,ppnl_radius,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (paw_atom_present) THEN
       DEALLOCATE (oce_present,oce_radius,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (all_potential_present) THEN
       DEALLOCATE (all_present,all_pot_rad,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_qs_neighbor_lists

! *****************************************************************************
!> \brief   Build simple pair neighbor lists.
!> \param ab_list ...
!> \param particle_set ...
!> \param atom ...
!> \param cell ...
!> \param pair_radius ...
!> \param subcells ...
!> \param mic ...
!> \param symmetric ...
!> \param molecular ...
!> \param subset_of_mol ...
!> \param current_subset ...
!> \param operator_type ...
!> \param name ...
!> \param error ...
!> \date    20.03.2002
!> \par History
!>          - Major refactoring (25.07.2010,jhu)
!> \author  MK
!> \version 2.0
! *****************************************************************************
  SUBROUTINE build_neighbor_lists(ab_list,particle_set,atom,cell,pair_radius,subcells,&
                                  mic,symmetric,molecular,subset_of_mol,current_subset,&
                                  operator_type,name,error)

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab_list
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: pair_radius
    REAL(dp), INTENT(IN)                     :: subcells
    LOGICAL, INTENT(IN), OPTIONAL            :: mic, symmetric, molecular
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: subset_of_mol
    INTEGER, OPTIONAL                        :: current_subset
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: operator_type
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_neighbor_lists', &
      routineP = moduleN//':'//routineN

    TYPE local_lists
      INTEGER,DIMENSION(:),POINTER           :: list
    END TYPE local_lists
    INTEGER :: atom_a, atom_b, handle, i, iab, iatom, iatom_local, &
      iatom_subcell, icell, ikind, j, jatom, jatom_local, jcell, jkind, k, &
      kcell, maxat, mol_a, mol_b, nkind, otype, stat, natom
    INTEGER, DIMENSION(3)                    :: cell_b, ncell, &
                                                nsubcell, periodic
    INTEGER, DIMENSION(:), POINTER           :: index_list
    LOGICAL                                  :: failure, include_ab, my_mic, &
                                                my_molecular, my_symmetric
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: pres_a, pres_b
    REAL(dp)                                 :: rab2, rab2_max, rab_max, rabm, deth, subcell_scale
    REAL(dp), DIMENSION(3)                   :: r, rab, ra, rb,sab_max, sb, &
                                                sb_pbc, sb_min, sb_max, rab_pbc, pd, sab_max_guard
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nlista, nlistb
    TYPE(local_lists), DIMENSION(:), POINTER :: lista, listb
    TYPE(neighbor_list_p_type), &
      ALLOCATABLE, DIMENSION(:)              :: kind_a
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(subcell_type), DIMENSION(:, :, :), &
      POINTER                                :: subcell
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: r_pbc

    CALL timeset(routineN//"_"//TRIM(name),handle)
    failure = .FALSE.

    ! input options
    my_mic = .FALSE.
    IF (PRESENT(mic)) my_mic = mic
    my_symmetric = .TRUE.
    IF (PRESENT(symmetric)) my_symmetric = symmetric
    my_molecular = .FALSE.
    ! if we have a molecular NL, MIC has to be used
    IF (PRESENT(molecular)) my_molecular = molecular
    ! check for operator types
    IF (PRESENT(operator_type)) THEN
      SELECT CASE (operator_type)
        CASE ("AB")
           otype = 1 ! simple overlap
        CASE ("ABC")
           otype = 2 ! for three center operators
           CPAssert(.NOT.my_molecular,cp_failure_level,routineP,error,failure)
           my_symmetric = .FALSE.
        CASE ("ABBA")
           otype = 3 ! for separable nonlocal operators
           my_symmetric = .FALSE.
        CASE ("PP")
           otype = 4 ! simple atomic pair potential list
        CASE default
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT
    ELSE
      ! default is a simple AB neighbor list
      otype = 1
    END IF

    ! Deallocate the old neighbor list structure
    IF (ASSOCIATED(ab_list)) THEN
       DO iab=1,SIZE(ab_list)
          CALL deallocate_neighbor_list_set(ab_list(iab)%neighbor_list_set)
       END DO
       DEALLOCATE (ab_list,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    nkind = SIZE(atom)
    ! Allocate and initialize the new neighbor list structure
    ALLOCATE (ab_list(nkind*nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO iab=1,SIZE(ab_list)
      NULLIFY (ab_list(iab)%neighbor_list_set)
    END DO

    ! Allocate and initialize the kind availability
    ALLOCATE (pres_a(nkind),pres_b(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      pres_a(ikind) = ANY(pair_radius(ikind,:) > 0._dp)
      pres_b(ikind) = ANY(pair_radius(:,ikind) > 0._dp)
    END DO

    ! create a copy of the pbc'ed coordinates
    natom=SIZE(particle_set)
    ALLOCATE(r_pbc(3,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,natom
       r_pbc(1:3,i)=pbc(particle_set(i)%r(1:3),cell)
    ENDDO

    ! setup the local lists of atoms
    maxat = 0
    DO ikind=1,nkind
      maxat=MAX(maxat,SIZE(atom(ikind)%list))
    END DO
    ALLOCATE (index_list(maxat),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,maxat
      index_list(i) = i
    END DO
    ALLOCATE (lista(nkind),listb(nkind),nlista(nkind),nlistb(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    nlista = 0
    nlistb = 0
    DO ikind=1,nkind
      NULLIFY(lista(ikind)%list,listb(ikind)%list)
      SELECT CASE (otype)
        CASE (1)
          IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
             lista(ikind)%list => atom(ikind)%list_local_a_index
             nlista(ikind) = SIZE(lista(ikind)%list)
          END IF
          IF (ASSOCIATED(atom(ikind)%list_local_b_index)) THEN
             listb(ikind)%list => atom(ikind)%list_local_b_index
             nlistb(ikind) = SIZE(listb(ikind)%list)
          END IF
        CASE (2)
          IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
             lista(ikind)%list => atom(ikind)%list_local_a_index
             nlista(ikind) = SIZE(lista(ikind)%list)
          END IF
          nlistb(ikind) = SIZE(atom(ikind)%list)
          listb(ikind)%list => index_list
        CASE (3)
          CALL combine_lists(lista(ikind)%list,nlista(ikind),ikind,atom,error)
          nlistb(ikind) = SIZE(atom(ikind)%list)
          listb(ikind)%list => index_list
        CASE (4)
          nlista(ikind) = SIZE(atom(ikind)%list_1d)
          lista(ikind)%list => atom(ikind)%list_1d
          nlistb(ikind) = SIZE(atom(ikind)%list)
          listb(ikind)%list => index_list
        CASE default
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT
    END DO

    ! Determine max. number of local atoms
    maxat = 0
    DO ikind=1,nkind
      maxat=MAX(maxat,nlista(ikind),nlistb(ikind))
    END DO
    ALLOCATE (kind_a(2*maxat),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Load informations about the simulation cell
    CALL get_cell(cell=cell,periodic=periodic,deth=deth)

    ! Loop over all atomic kind pairs
    DO ikind=1,nkind
      IF (.NOT.pres_a(ikind)) CYCLE

      DO jkind=1,nkind
        IF (.NOT.pres_b(jkind)) CYCLE

        iab = ikind + nkind*(jkind - 1)

        ! Calculate the square of the maximum interaction distance
        IF ( pair_radius(ikind,jkind) <= 0._dp ) CYCLE
        rab_max = pair_radius(ikind,jkind)
        IF (otype==3) THEN
          ! Calculate the square of the maximum interaction distance
          ! for sac_max / ncell this must be the maximum over all kinds
          ! to be correct for three center terms involving different kinds
          rabm = MAXVAL(pair_radius(:,jkind))
        ELSE
          rabm = rab_max
        END IF
        rab2_max = rabm*rabm

        pd(1) = plane_distance(1,0,0,cell)
        pd(2) = plane_distance(0,1,0,cell)
        pd(3) = plane_distance(0,0,1,cell)

        sab_max = rabm / pd
        sab_max_guard = 15.0_dp / pd

        ! It makes sense to have fewer subcells for larger systems
        subcell_scale=((125.0_dp**3)/deth)**(1.0_dp/6.0_dp)

        ! guess the number of subcells for optimal performance, guard against crazy stuff triggered by very small rabm
        nsubcell(:) = INT(MAX(1.0_dp,MIN(0.5_dp*subcells*subcell_scale/sab_max(:),0.5_dp*subcells*subcell_scale/sab_max_guard(:))))

        ! number of image cells to be considered
        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)

        CALL allocate_neighbor_list_set(neighbor_list_set=ab_list(iab)%neighbor_list_set,&
             r_max=rabm,symmetric=my_symmetric,mic=my_mic,molecular=my_molecular)
        neighbor_list_set => ab_list(iab)%neighbor_list_set

        DO iatom_local=1,nlista(ikind)
          iatom = lista(ikind)%list(iatom_local)
          atom_a = atom(ikind)%list(iatom)
          CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                                 atom=atom_a,&
                                 neighbor_list=kind_a(iatom_local)%neighbor_list)
        END DO

        CALL allocate_subcell(subcell,nsubcell,error=error)
        DO iatom_local=1,nlista(ikind)
          iatom = lista(ikind)%list(iatom_local)
          atom_a = atom(ikind)%list(iatom)
          r = r_pbc(:,atom_a)
          CALL give_ijk_subcell(r,i,j,k,cell,nsubcell)
          subcell(i,j,k)%natom =  subcell(i,j,k)%natom + 1
        END DO
        DO k = 1,nsubcell(3)
           DO j = 1,nsubcell(2)
              DO i = 1,nsubcell(1)
                 maxat = subcell(i,j,k)%natom + subcell(i,j,k)%natom/10
                 ALLOCATE(subcell(i,j,k)%atom_list(maxat),stat=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                 subcell(i,j,k)%natom   = 0
              END DO
           END DO
        END DO
        DO iatom_local=1,nlista(ikind)
          iatom = lista(ikind)%list(iatom_local)
          atom_a = atom(ikind)%list(iatom)
          r = r_pbc(:,atom_a)
          CALL give_ijk_subcell(r,i,j,k,cell,nsubcell)
          subcell(i,j,k)%natom =  subcell(i,j,k)%natom + 1
          subcell(i,j,k)%atom_list(subcell(i,j,k)%natom) = iatom_local
        END DO

        DO jatom_local=1,nlistb(jkind)
          jatom = listb(jkind)%list(jatom_local)
          atom_b = atom(jkind)%list(jatom)
          IF (my_molecular) THEN
            mol_b = atom(jkind)%list_b_mol(jatom_local)
            IF (PRESENT(subset_of_mol)) THEN 
              IF (subset_of_mol(mol_b).NE.current_subset) CYCLE
            END IF
          END IF
          r = r_pbc(:,atom_b)
          CALL real_to_scaled(sb_pbc(:),r(:),cell)

          loop2_kcell: DO kcell=-ncell(3),ncell(3)
            sb(3) = sb_pbc(3) + REAL(kcell,dp)
            sb_min(3) = sb(3) - sab_max(3)
            sb_max(3) = sb(3) + sab_max(3)
            IF (periodic(3) /= 0) THEN
              IF (sb_min(3) >= 0.5_dp) EXIT loop2_kcell
              IF (sb_max(3) < -0.5_dp) CYCLE loop2_kcell
            END IF
            cell_b(3) = kcell

            loop2_jcell: DO jcell=-ncell(2),ncell(2)
              sb(2) = sb_pbc(2) + REAL(jcell,dp)
              sb_min(2) = sb(2) - sab_max(2)
              sb_max(2) = sb(2) + sab_max(2)
              IF (periodic(2) /= 0) THEN
                IF (sb_min(2) >= 0.5_dp) EXIT loop2_jcell
                IF (sb_max(2) < -0.5_dp) CYCLE loop2_jcell
              END IF
              cell_b(2) = jcell

              loop2_icell: DO icell=-ncell(1),ncell(1)
                sb(1) = sb_pbc(1) + REAL(icell,dp)
                sb_min(1) = sb(1) - sab_max(1)
                sb_max(1) = sb(1) + sab_max(1)
                IF (periodic(1) /= 0) THEN
                  IF (sb_min(1) >= 0.5_dp) EXIT loop2_icell
                  IF (sb_max(1) < -0.5_dp) CYCLE loop2_icell
                END IF
                cell_b(1) = icell

                CALL scaled_to_real(rb,sb,cell)

                loop_k: DO k=1,nsubcell(3)
                  loop_j: DO j=1,nsubcell(2)
                    loop_i: DO i=1,nsubcell(1)

                      ! FIXME for non-periodic systems, the whole subcell trick is skipped
                      ! yielding a Natom**2 pair list build.
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

                      DO iatom_subcell=1,subcell(i,j,k)%natom
                        iatom_local = subcell(i,j,k)%atom_list(iatom_subcell)
                        iatom = lista(ikind)%list(iatom_local)
                        atom_a = atom(ikind)%list(iatom)
                        IF (my_molecular) THEN
                          mol_a = atom(ikind)%list_a_mol(iatom_local)
                          IF (mol_a /= mol_b) CYCLE
                        END IF
                        IF (my_symmetric) THEN
                          IF (atom_a > atom_b) THEN
                            include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                          ELSE
                            include_ab = (MODULO(atom_a + atom_b,2) == 0)
                          END IF
                        ELSE
                          include_ab = .TRUE.
                        END IF
                        IF (include_ab) THEN
                          ra(:) = r_pbc(:,atom_a)
                          rab(:) = rb(:) - ra(:)
                          rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                          IF (rab2 < rab2_max) THEN
                            include_ab=.TRUE.
                            IF (my_mic) THEN
                               ! only if rab is minimum image the pair will be included
                               ! ideally the range of the pair list is < L/2 so that this never triggers
                               rab_pbc(:) = pbc(rab(:),cell)
                               IF (SUM((rab_pbc-rab)**2) > EPSILON(1.0_dp)) THEN
                                  include_ab = .FALSE.
                               ENDIF
                            ENDIF
                            IF (include_ab) THEN
                               CALL add_neighbor_node(&
                                    neighbor_list=kind_a(iatom_local)%neighbor_list,&
                                    neighbor=atom_b,&
                                    cell=cell_b,&
                                    r=rab,&
                                    nkind=nkind)
                            ENDIF
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

        CALL deallocate_subcell(subcell,error=error)

      END DO
    END DO

    SELECT CASE (otype)
      CASE (1:2,4)
      CASE (3)
         DO ikind=1,nkind
           DEALLOCATE(lista(ikind)%list,STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         END DO
      CASE default
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    DEALLOCATE (kind_a,pres_a,pres_b,lista,listb,nlista,nlistb,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (index_list,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(r_pbc,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE build_neighbor_lists

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param list ...
!> \param n ...
!> \param ikind ...
!> \param atom ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE combine_lists(list,n,ikind,atom,error)
    INTEGER, DIMENSION(:), POINTER           :: list
    INTEGER, INTENT(OUT)                     :: n
    INTEGER, INTENT(IN)                      :: ikind
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'combine_lists', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ib, ierr, na, nb
    INTEGER, DIMENSION(:), POINTER           :: lista, listb
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(.NOT.ASSOCIATED(list),cp_failure_level,routineP,error,failure)

    lista => atom(ikind)%list_local_a_index
    listb => atom(ikind)%list_local_b_index

    IF (ASSOCIATED(lista)) THEN
       na = SIZE(lista)
    ELSE
       na = 0
    END IF

    IF (ASSOCIATED(listb)) THEN
       nb = SIZE(listb)
    ELSE
       nb = 0
    END IF

    ALLOCATE(list(na+nb),STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    n = na
    IF(na.GT.0) list(1:na) = lista(1:na)
    IF(nb.GT.0) THEN
       loopb: DO ib = 1, nb
          DO i = 1, na
             IF(listb(ib) == list(i)) CYCLE loopb
          END DO
          n = n+1
          list(n) = listb(ib)
       END DO loopb
    ENDIF
  END SUBROUTINE combine_lists

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param present_a ...
!> \param present_b ...
!> \param radius_a ...
!> \param radius_b ...
!> \param pair_radius ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pair_radius_setup(present_a,present_b,radius_a,radius_b,pair_radius,error)
    LOGICAL, DIMENSION(:), INTENT(IN)        :: present_a, present_b
    REAL(dp), DIMENSION(:), INTENT(IN)       :: radius_a, radius_b
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: pair_radius
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pair_radius_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nkind

    nkind = SIZE(present_a)

    pair_radius = 0._dp

    DO i=1,nkind
      IF (.NOT. present_a(i)) CYCLE
      DO j=1,nkind
        IF (.NOT. present_b(j)) CYCLE
        pair_radius(i,j) = radius_a(i) + radius_b(j)
      END DO
    END DO

  END SUBROUTINE pair_radius_setup

! *****************************************************************************
!> \brief   Print the distribution of the simple pair neighbor list.
!> \param ab ...
!> \param qs_kind_set ...
!> \param output_unit ...
!> \param para_env ...
!> \param error ...
!> \date    19.06.2003
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_neighbor_distribution(ab,qs_kind_set,output_unit,para_env,error)
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER, INTENT(in)                      :: output_unit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_neighbor_distribution', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: full_output = .FALSE.

    INTEGER                                  :: group, handle, ikind, inode, &
                                                ipe, jkind, mype, n, nkind, &
                                                nnode, npe, stat
    INTEGER(int_8)                           :: nblock_max, nblock_sum, &
                                                nelement_max, nelement_sum, &
                                                tmp(2)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nblock, nelement, nnsgf
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    CALL timeset(routineN,handle)
    failure = .FALSE.
    group = para_env%group
    mype = para_env%mepos + 1
    npe = para_env%num_pe

    ! Allocate work storage
    ALLOCATE (nblock(npe),nelement(npe),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    nblock(:) = 0
    nelement(:) = 0
    nkind = SIZE(qs_kind_set)
    ALLOCATE (nnsgf(nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    nnsgf = 1
    DO ikind=1,nkind
      CALL get_qs_kind(qs_kind_set(ikind),orb_basis_set=orb_basis_set)
      IF (ASSOCIATED(orb_basis_set)) THEN
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nnsgf(ikind))
      END IF
    END DO


    CALL neighbor_list_iterator_create(nl_iterator,ab)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,nnode=nnode)
       IF (inode==1) THEN
          n = nnsgf(ikind)*nnsgf(jkind)
          nblock(mype) = nblock(mype) + nnode
          nelement(mype) = nelement(mype) + n*nnode
       END IF
    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    IF (full_output) THEN
       ! XXXXXXXX should gather/scatter this on ionode
       CALL mp_sum(nblock,group)
       CALL mp_sum(nelement,group)

       nblock_sum = SUM(INT(nblock,KIND=int_8))
       nelement_sum = SUM(INT(nelement,KIND=int_8))
    ELSE
       nblock_sum = nblock(mype)
       nblock_max = nblock(mype)
       nelement_sum = nelement(mype)
       nelement_max = nelement(mype)
       tmp=(/nblock_sum,nelement_sum/)
       CALL mp_sum(tmp,group)
       nblock_sum=tmp(1) ; nelement_sum=tmp(2)
       tmp=(/nblock_max,nelement_max/)
       CALL mp_max(tmp,group)
       nblock_max=tmp(1) ; nelement_max=tmp(2)
    ENDIF

    IF (output_unit > 0) THEN
      IF (full_output) THEN
         WRITE (UNIT=output_unit,&
                FMT="(/,/,T2,A,/,/,T3,A,/,/,(T4,I6,T27,I10,T55,I10))")&
           "DISTRIBUTION OF THE NEIGHBOR LISTS",&
           "Process   Number of particle pairs   Number of matrix elements",&
           (ipe-1,nblock(ipe),nelement(ipe),ipe=1,npe)
         WRITE (UNIT=output_unit,FMT="(/,T7,A3,T27,I10,T55,I10)")&
           "Sum",SUM(nblock),SUM(nelement)
      ELSE
         WRITE (UNIT=output_unit,FMT="(/,T2,A)") "DISTRIBUTION OF THE NEIGHBOR LISTS"
         WRITE (UNIT=output_unit,FMT="(T15,A,T68,I13)") "Total number of particle pairs:",nblock_sum
         WRITE (UNIT=output_unit,FMT="(T15,A,T68,I13)") "Total number of matrix elements:",nelement_sum
         WRITE (UNIT=output_unit,FMT="(T15,A,T68,I13)") "Average number of particle pairs:",(nblock_sum+npe-1)/npe
         WRITE (UNIT=output_unit,FMT="(T15,A,T68,I13)") "Maximum number of particle pairs:",nblock_max
         WRITE (UNIT=output_unit,FMT="(T15,A,T68,I13)") "Average number of matrix element:",(nelement_sum+npe-1)/npe
         WRITE (UNIT=output_unit,FMT="(T15,A,T68,I13)") "Maximum number of matrix elements:",nelement_max
      ENDIF
    END IF

    ! Release work storage

    DEALLOCATE (nblock,nelement,nnsgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE write_neighbor_distribution

! *****************************************************************************
!> \brief   Write a set of neighbor lists to the output unit.
!> \param ab ...
!> \param particle_set ...
!> \param cell ...
!> \param para_env ...
!> \param neighbor_list_section ...
!> \param nl_type ...
!> \param middle_name ...
!> \param name ...
!> \param error ...
!> \date    04.03.2002
!> \par History
!>       - Adapted to the new parallelized neighbor list version
!>         (26.06.2003,MK)
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_neighbor_lists(ab,particle_set,cell,para_env,neighbor_list_section,&
    nl_type,middle_name,name,error)

    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: ab
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: neighbor_list_section
    CHARACTER(LEN=*), INTENT(IN)             :: nl_type, middle_name, name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=default_string_length)     :: string, unit_str
    INTEGER                                  :: iatom, inode, iw, jatom, &
                                                mype, nneighbor, nnode
    INTEGER, DIMENSION(3)                    :: cell_b
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dab, unit_conv
    REAL(dp), DIMENSION(3)                   :: ra, rab, rb
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,neighbor_list_section,&
                                         TRIM(nl_type),error=error),&
              cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger=logger,&
                                 basis_section=neighbor_list_section,&
                                 print_key_path=TRIM(nl_type),&
                                 extension=".out",&
                                 middle_name=TRIM(middle_name),&
                                 local=.TRUE.,&
                                 log_filename=.FALSE.,&
                                 file_position="REWIND",&
                                 error=error)
       mype = para_env%mepos
       CALL section_vals_val_get(neighbor_list_section,"UNIT",c_val=unit_str,error=error)
       unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

       ! Print headline
       string = ""
       WRITE (UNIT=string,FMT="(A,I5,A)")&
         TRIM(name)//" IN "//TRIM(unit_str)//" (PROCESS",mype,")"
       CALL compress(string)
       IF (iw > 0) WRITE (UNIT=iw,FMT="(/,/,T2,A)") TRIM(string)

       nneighbor = 0

       CALL neighbor_list_iterator_create(nl_iterator,ab)
       DO WHILE (neighbor_list_iterate(nl_iterator)==0)
         CALL get_iterator_info(nl_iterator,inode=inode,nnode=nnode,&
                                iatom=iatom,jatom=jatom,cell=cell_b,r=rab)
         nneighbor = nneighbor + 1
         ra(:) = pbc(particle_set(iatom)%r,cell)
         rb(:) = ra(:) + rab(:)
         dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
         IF (iw > 0) THEN
           IF(inode==1) THEN
             WRITE (UNIT=iw,FMT="(/,T2,I5,3X,I6,3X,3F12.6)")&
               iatom,nnode,ra(1:3)*unit_conv
           END IF
           WRITE (UNIT=iw,FMT="(T10,I6,3X,3I4,3F12.6,2X,F12.6)")&
                 jatom,cell_b(1:3),rb(1:3)*unit_conv,dab*unit_conv
         END IF
       END DO
       CALL neighbor_list_iterator_release(nl_iterator)

       string = ""
       WRITE (UNIT=string,FMT="(A,I12,A,I12)")&
         "Total number of neighbor interactions for process",mype,":",&
         nneighbor
       CALL compress(string)
       IF (iw > 0) WRITE (UNIT=iw,FMT="(/,T2,A)") TRIM(string)
       CALL cp_print_key_finished_output(unit_nr=iw,&
                                         logger=logger,&
                                         basis_section=neighbor_list_section,&
                                         print_key_path=TRIM(nl_type),&
                                         local=.TRUE.,&
                                         error=error)
    END IF

  END SUBROUTINE write_neighbor_lists

END MODULE qs_neighbor_lists
