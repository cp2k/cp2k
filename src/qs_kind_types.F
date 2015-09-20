!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Define the quickstep kind type and their sub types
!> \author  Ole Schuett
!>
!> <b>Modification history:</b>
!> - 01.2002 creation [MK]
!> - 04.2002 added pao [fawzi]
!> - 09.2002 adapted for POL/KG use [GT]
!> - 02.2004 flexible normalization of basis sets [jgh]
!> - 03.2004 attach/detach routines [jgh]
!> - 10.2004 removed pao [fawzi]
!> - 08.2014 separated qs-related stuff from atomic_kind_types.F [Ole Schuett]
!> - 07.2015 new container for basis sets [jgh]
! *****************************************************************************
MODULE qs_kind_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_container_types,       ONLY: add_basis_set_to_container,&
                                             basis_set_container_type,&
                                             get_basis_from_container,&
                                             remove_basis_from_container,&
                                             remove_basis_set_container
  USE basis_set_types,                 ONLY: &
       allocate_gto_basis_set, combine_basis_sets, deallocate_gto_basis_set, &
       get_gto_basis_set, gto_basis_set_type, init_aux_basis_set, &
       init_orb_basis_set, read_gto_basis_set, write_gto_basis_set, &
       write_orb_basis_set
  USE cp_control_types,                ONLY: dft_control_type,&
                                             qs_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE external_potential_types,        ONLY: &
       all_potential_type, allocate_potential, deallocate_potential, &
       get_potential, gth_potential_type, init_potential, &
       local_potential_type, read_potential, set_default_all_potential, &
       set_potential, write_potential
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_dftb, do_method_mndo, do_method_mndod, &
       do_method_pdg, do_method_pm3, do_method_pm6, do_method_pnnl, &
       do_method_rm1, do_method_scptb, do_qs
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE paw_proj_set_types,              ONLY: allocate_paw_proj_set,&
                                             deallocate_paw_proj_set,&
                                             get_paw_proj_set,&
                                             paw_proj_set_type,&
                                             projectors
  USE periodic_table,                  ONLY: get_ptable_info,&
                                             ptable
  USE physcon,                         ONLY: bohr,&
                                             evolt
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: deallocate_dftb_atom_param,&
                                             get_dftb_atom_param,&
                                             write_dftb_atom_param
  USE qs_dispersion_types,             ONLY: qs_atom_dispersion_type
  USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                             deallocate_grid_atom,&
                                             grid_atom_type
  USE qs_harmonics_atom,               ONLY: allocate_harmonics_atom,&
                                             deallocate_harmonics_atom,&
                                             harmonics_atom_type
  USE scptb_types,                     ONLY: deallocate_scptb_parameter,&
                                             get_scptb_parameter,&
                                             scptb_parameter_type,&
                                             write_scptb_parameter
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_create,&
                                             semi_empirical_release,&
                                             semi_empirical_type,&
                                             write_se_param
  USE semi_empirical_utils,            ONLY: init_se_param,&
                                             se_param_set_default
  USE soft_basis_set,                  ONLY: create_soft_basis
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kind_types'

!> Input parameters for the DFT+U method
! *****************************************************************************
  TYPE dft_plus_u_type
     INTEGER                                :: l,max_scf
     REAL(KIND=dp)                          :: eps_u_ramping,&
                                               eps_scf,&
                                               u_minus_j_target,&
                                               u_minus_j,&
                                               u_ramping
     INTEGER, DIMENSION(:), POINTER         :: orbitals
     LOGICAL                                :: init_u_ramping_each_scf,&
                                               smear
  END TYPE dft_plus_u_type

!> Provides all information about a quickstep kind
! *****************************************************************************
  TYPE qs_kind_type
     CHARACTER(LEN=default_string_length)   :: name = ""
     CHARACTER(LEN=2)                       :: element_symbol = ""
     INTEGER                                :: natom = -1
     TYPE(all_potential_type), POINTER      :: all_potential     => Null()
     TYPE(local_potential_type), POINTER    :: tnadd_potential   => Null()
     TYPE(gth_potential_type), POINTER      :: gth_potential     => Null()
     TYPE(semi_empirical_type), POINTER     :: se_parameter      => Null()
     TYPE(qs_dftb_atom_type), POINTER       :: dftb_parameter    => Null()
     TYPE(scptb_parameter_type), POINTER    :: scptb_parameter   => Null()
     !
     TYPE(gto_basis_set_type), POINTER      :: soft_basis_set    => Null()
     TYPE(gto_basis_set_type), POINTER      :: hard_basis_set    => Null()
     !
     TYPE(basis_set_container_type), &
       DIMENSION(20)                        :: basis_sets
     !
     TYPE(paw_proj_set_type),  POINTER      :: paw_proj_set      => Null()
     !
     REAL(dp)                               :: hard_radius  = 0.8_dp*bohr    ! for hard and soft exp
     REAL(dp)                               :: hard0_radius = 0.8_dp*bohr    ! for hard exp of rho0
     REAL(dp)                               :: max_rad_local = 13.2_dp*bohr  ! max GTO radius used in GAPW
     LOGICAL                                :: paw_atom = .FALSE.        ! needs atomic rho1
     LOGICAL                                :: gpw_type_forced = .FALSE. ! gpw atom even if with hard exponents
     LOGICAL                                :: ghost = .FALSE.
     REAL(KIND = dp)                        :: dudq_dftb3 = 0.0_dp
     INTEGER, DIMENSION(:,:), POINTER       :: addel  => Null()
     INTEGER, DIMENSION(:,:), POINTER       :: laddel => Null()
     INTEGER, DIMENSION(:,:), POINTER       :: naddel => Null()
     TYPE(harmonics_atom_type), POINTER     :: harmonics => Null()
     TYPE(grid_atom_type), POINTER          :: grid_atom => Null()
     INTEGER                                :: ngrid_rad = 50
     INTEGER                                :: ngrid_ang = 50
     INTEGER                                :: lmax_rho0 = 0
     INTEGER                                :: mao = -1
     INTEGER, DIMENSION(:), POINTER         :: elec_conf => Null() ! used to set up the initial atomic guess
     LOGICAL                                :: bs_occupation = .FALSE.
     TYPE(dft_plus_u_type), POINTER         :: dft_plus_u => Null()
     LOGICAL                                :: no_optimize = .TRUE.
     !
     TYPE(qs_atom_dispersion_type), POINTER :: dispersion => Null()
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: reltmat => Null()
  END TYPE qs_kind_type

!> Provides a vector of pointers of type qs_kind_type
! *****************************************************************************
  TYPE qs_kind_p_type
     TYPE(qs_kind_type),DIMENSION(:),&
          POINTER                             :: qs_kind_set
  END TYPE qs_kind_p_type

  ! Public subroutines

  PUBLIC :: check_qs_kind_set,&
            deallocate_qs_kind_set,&
            get_qs_kind,&
            get_qs_kind_set,&
            init_qs_kind_set,&
            init_gapw_basis_set,&
            create_qs_kind_set,&
            set_qs_kind,&
            write_qs_kind_set,&
            write_gto_basis_sets,&
            init_atom_electronic_state

  ! Public data types
  PUBLIC :: qs_kind_type

CONTAINS

! *****************************************************************************
!> \brief   Destructor routine for a set of qs kinds
!> \param qs_kind_set ...
!> \date    02.01.2002
!> \author  Matthias Krack (MK)
!> \version 2.0
! *****************************************************************************
  SUBROUTINE deallocate_qs_kind_set(qs_kind_set)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_qs_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(qs_kind_set)) THEN

       nkind = SIZE(qs_kind_set)

       DO ikind=1,nkind
          IF (ASSOCIATED(qs_kind_set(ikind)%all_potential)) THEN
             CALL deallocate_potential(qs_kind_set(ikind)%all_potential)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%tnadd_potential)) THEN
             CALL deallocate_potential(qs_kind_set(ikind)%tnadd_potential)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%gth_potential)) THEN
             CALL deallocate_potential(qs_kind_set(ikind)%gth_potential)
          ENDIF
          IF (ASSOCIATED(qs_kind_set(ikind)%se_parameter)) THEN
             CALL semi_empirical_release(qs_kind_set(ikind)%se_parameter)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%dftb_parameter)) THEN
             CALL deallocate_dftb_atom_param(qs_kind_set(ikind)%dftb_parameter)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%scptb_parameter)) THEN
             CALL deallocate_scptb_parameter(qs_kind_set(ikind)%scptb_parameter)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%soft_basis_set).AND.&
               qs_kind_set(ikind)%paw_atom) THEN
             CALL deallocate_gto_basis_set(qs_kind_set(ikind)%soft_basis_set)
          ELSEIF(ASSOCIATED(qs_kind_set(ikind)%soft_basis_set).AND.&
               (.NOT.qs_kind_set(ikind)%paw_atom)) THEN
             NULLIFY (qs_kind_set(ikind)%soft_basis_set)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%paw_proj_set)) THEN
             CALL deallocate_paw_proj_set(qs_kind_set(ikind)%paw_proj_set)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%harmonics)) THEN
             CALL deallocate_harmonics_atom(qs_kind_set(ikind)%harmonics)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%grid_atom)) THEN
             CALL deallocate_grid_atom(qs_kind_set(ikind)%grid_atom)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%elec_conf)) THEN
             DEALLOCATE (qs_kind_set(ikind)%elec_conf,STAT=stat)
             CPASSERT(stat==0)
          END IF

          IF (ASSOCIATED(qs_kind_set(ikind)%dft_plus_u)) THEN
             IF (ASSOCIATED(qs_kind_set(ikind)%dft_plus_u%orbitals)) THEN
                DEALLOCATE (qs_kind_set(ikind)%dft_plus_u%orbitals,STAT=stat)
                CPASSERT(stat==0)
             END IF
             DEALLOCATE (qs_kind_set(ikind)%dft_plus_u,STAT=stat)
             CPASSERT(stat==0)
          END IF

          IF (ASSOCIATED(qs_kind_set(ikind)%dispersion)) THEN
             DEALLOCATE (qs_kind_set(ikind)%dispersion,STAT=stat)
             CPASSERT(stat==0)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%addel)) THEN
             DEALLOCATE (qs_kind_set(ikind)%addel,STAT=stat)
             CPASSERT(stat==0)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%naddel)) THEN
             DEALLOCATE (qs_kind_set(ikind)%naddel,STAT=stat)
             CPASSERT(stat==0)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%laddel)) THEN
             DEALLOCATE (qs_kind_set(ikind)%laddel,STAT=stat)
             CPASSERT(stat==0)
          END IF
          IF (ASSOCIATED(qs_kind_set(ikind)%reltmat)) THEN
             DEALLOCATE (qs_kind_set(ikind)%reltmat,STAT=stat)
             CPASSERT(stat==0)
          END IF

          CALL remove_basis_set_container(qs_kind_set(ikind)%basis_sets)

       END DO
       DEALLOCATE (qs_kind_set,STAT=stat)
       CPASSERT(stat==0)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer qs_kind_set is not associated and "//&
            "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_qs_kind_set

! *****************************************************************************
!> \brief Get attributes of an atomic kind.
!> \param qs_kind ...
!> \param basis_set ...
!> \param basis_type ...
!> \param ncgf ...
!> \param nsgf ...
!> \param all_potential ...
!> \param tnadd_potential ...
!> \param gth_potential ...
!> \param se_parameter ...
!> \param dftb_parameter ...
!> \param scptb_parameter ...
!> \param dftb3_param ...
!> \param zeff ...
!> \param elec_conf ...
!> \param mao ...
!> \param alpha_core_charge ...
!> \param ccore_charge ...
!> \param core_charge ...
!> \param core_charge_radius ...
!> \param soft_basis_set ...
!> \param hard_basis_set ...
!> \param paw_proj_set ...
!> \param softb ...
!> \param paw_atom ...
!> \param hard_radius ...
!> \param hard0_radius ...
!> \param max_rad_local ...
!> \param gpw_type_forced ...
!> \param harmonics ...
!> \param max_iso_not0 ...
!> \param max_s_harm ...
!> \param grid_atom ...
!> \param ngrid_ang ...
!> \param ngrid_rad ...
!> \param lmax_rho0 ...
!> \param dft_plus_u_atom ...
!> \param l_of_dft_plus_u ...
!> \param u_minus_j ...
!> \param dispersion ...
!> \param bs_occupation ...
!> \param no_optimize ...
!> \param addel ...
!> \param laddel ...
!> \param naddel ...
!> \param orbitals ...
!> \param max_scf ...
!> \param eps_scf ...
!> \param smear ...
!> \param u_ramping ...
!> \param u_minus_j_target ...
!> \param eps_u_ramping ...
!> \param init_u_ramping_each_scf ...
!> \param reltmat ...
!> \param ghost ...
!> \param name ...
!> \param element_symbol ...
! *****************************************************************************
  SUBROUTINE get_qs_kind(qs_kind, &
       basis_set, basis_type, ncgf, nsgf, &
       all_potential, tnadd_potential, gth_potential, &
       se_parameter, dftb_parameter, scptb_parameter, &
       dftb3_param, zeff, elec_conf, mao,&
       alpha_core_charge, ccore_charge, core_charge, core_charge_radius,&
       soft_basis_set, hard_basis_set, paw_proj_set, softb, &
       paw_atom, hard_radius, hard0_radius, max_rad_local, &
       gpw_type_forced, harmonics, max_iso_not0, max_s_harm, grid_atom, &
       ngrid_ang, ngrid_rad, lmax_rho0, &
       dft_plus_u_atom, l_of_dft_plus_u, u_minus_j, dispersion, &
       bs_occupation, no_optimize, addel, laddel, naddel, orbitals, &
       max_scf, eps_scf, smear, u_ramping, u_minus_j_target, eps_u_ramping, &
       init_u_ramping_each_scf, reltmat, ghost, name, element_symbol)

    TYPE(qs_kind_type)                       :: qs_kind
    TYPE(gto_basis_set_type), OPTIONAL, &
      POINTER                                :: basis_set
    CHARACTER(len=*), OPTIONAL               :: basis_type
    INTEGER, INTENT(OUT), OPTIONAL           :: ncgf, nsgf
    TYPE(all_potential_type), OPTIONAL, &
      POINTER                                :: all_potential
    TYPE(local_potential_type), OPTIONAL, &
      POINTER                                :: tnadd_potential
    TYPE(gth_potential_type), OPTIONAL, &
      POINTER                                :: gth_potential
    TYPE(semi_empirical_type), OPTIONAL, &
      POINTER                                :: se_parameter
    TYPE(qs_dftb_atom_type), OPTIONAL, &
      POINTER                                :: dftb_parameter
    TYPE(scptb_parameter_type), OPTIONAL, &
      POINTER                                :: scptb_parameter
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: dftb3_param, zeff
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: elec_conf
    INTEGER, INTENT(OUT), OPTIONAL           :: mao
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha_core_charge, &
                                                ccore_charge, core_charge, &
                                                core_charge_radius
    TYPE(gto_basis_set_type), OPTIONAL, &
      POINTER                                :: soft_basis_set, hard_basis_set
    TYPE(paw_proj_set_type), OPTIONAL, &
      POINTER                                :: paw_proj_set
    LOGICAL, INTENT(IN), OPTIONAL            :: softb
    LOGICAL, INTENT(OUT), OPTIONAL           :: paw_atom
    REAL(dp), INTENT(OUT), OPTIONAL          :: hard_radius, hard0_radius, &
                                                max_rad_local
    LOGICAL, INTENT(OUT), OPTIONAL           :: gpw_type_forced
    TYPE(harmonics_atom_type), OPTIONAL, &
      POINTER                                :: harmonics
    INTEGER, INTENT(OUT), OPTIONAL           :: max_iso_not0, max_s_harm
    TYPE(grid_atom_type), OPTIONAL, POINTER  :: grid_atom
    INTEGER, INTENT(OUT), OPTIONAL           :: ngrid_ang, ngrid_rad, &
                                                lmax_rho0
    LOGICAL, INTENT(OUT), OPTIONAL           :: dft_plus_u_atom
    INTEGER, INTENT(OUT), OPTIONAL           :: l_of_dft_plus_u
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: u_minus_j
    TYPE(qs_atom_dispersion_type), &
      OPTIONAL, POINTER                      :: dispersion
    LOGICAL, INTENT(OUT), OPTIONAL           :: bs_occupation, no_optimize
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: addel, laddel, naddel
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: orbitals
    INTEGER, OPTIONAL                        :: max_scf
    REAL(KIND=dp), OPTIONAL                  :: eps_scf
    LOGICAL, OPTIONAL                        :: smear
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: u_ramping, u_minus_j_target, &
                                                eps_u_ramping
    LOGICAL, OPTIONAL                        :: init_u_ramping_each_scf
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: reltmat
    LOGICAL, OPTIONAL                        :: ghost
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    CHARACTER(LEN=2), INTENT(OUT), OPTIONAL  :: element_symbol

    CHARACTER(len=*), PARAMETER :: routineN = 'get_qs_kind', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_basis_type
    INTEGER                                  :: l
    LOGICAL                                  :: softb_local
    TYPE(gto_basis_set_type), POINTER        :: tmp_basis_set

! Retrieve basis set from the kind container

    softb_local = .FALSE.
    IF (PRESENT(softb)) softb_local = softb

    IF ( PRESENT(basis_type) ) THEN
       my_basis_type = basis_type
    ELSE
       my_basis_type = "ORB"
    END IF

    IF ( PRESENT(basis_set) ) THEN
       IF(softb_local) THEN
          basis_set => qs_kind%soft_basis_set
       ELSE
          CALL get_basis_from_container(qs_kind%basis_sets,basis_set=basis_set,&
               basis_type=my_basis_type)
       END IF
    END IF

    IF (PRESENT(ncgf)) THEN
       CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis_set,&
            basis_type=my_basis_type)
       IF(ASSOCIATED(tmp_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,ncgf=ncgf)
       ELSE IF(ASSOCIATED(qs_kind%dftb_parameter)) THEN
          l = qs_kind%dftb_parameter%lmax
          ncgf = ((l+1)*(l+2)*(l+3))/6
       ELSE
          ncgf = 0
       END IF
    END IF

    IF (PRESENT(nsgf)) THEN
       CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis_set,&
            basis_type=my_basis_type)
       IF(ASSOCIATED(tmp_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,nsgf=nsgf)
       ELSE IF(ASSOCIATED(qs_kind%dftb_parameter)) THEN
          nsgf = qs_kind%dftb_parameter%natorb
       ELSE
          nsgf = 0
       END IF
    END IF

    IF (PRESENT(all_potential)) all_potential => qs_kind%all_potential
    IF (PRESENT(tnadd_potential)) tnadd_potential => qs_kind%tnadd_potential
    IF (PRESENT(gth_potential)) gth_potential => qs_kind%gth_potential
    IF (PRESENT(se_parameter))  se_parameter => qs_kind%se_parameter
    IF (PRESENT(dftb_parameter))  dftb_parameter => qs_kind%dftb_parameter
    IF (PRESENT(scptb_parameter))  scptb_parameter => qs_kind%scptb_parameter

    IF (PRESENT(element_symbol)) element_symbol = qs_kind%element_symbol
    IF (PRESENT(name)) name = qs_kind%name
    IF (PRESENT(dftb3_param)) dftb3_param = qs_kind%dudq_dftb3
    IF (PRESENT(elec_conf)) elec_conf => qs_kind%elec_conf
    IF (PRESENT(alpha_core_charge)) THEN
       IF (ASSOCIATED(qs_kind%all_potential)) THEN
          CALL get_potential(potential=qs_kind%all_potential,&
               alpha_core_charge=alpha_core_charge)
       ELSE IF (ASSOCIATED(qs_kind%gth_potential)) THEN
          CALL get_potential(potential=qs_kind%gth_potential,&
               alpha_core_charge=alpha_core_charge)
       ELSE
          alpha_core_charge = 1.0_dp
       END IF
    END IF
    IF (PRESENT(ccore_charge)) THEN
       IF (ASSOCIATED(qs_kind%all_potential)) THEN
          CALL get_potential(potential=qs_kind%all_potential,&
               ccore_charge=ccore_charge)
       ELSE IF (ASSOCIATED(qs_kind%gth_potential)) THEN
          CALL get_potential(potential=qs_kind%gth_potential,&
               ccore_charge=ccore_charge)
       ELSE
          ccore_charge = 0.0_dp
       END IF
    END IF
    IF (PRESENT(core_charge_radius)) THEN
       IF (ASSOCIATED(qs_kind%all_potential)) THEN
          CALL get_potential(potential=qs_kind%all_potential,&
               core_charge_radius=core_charge_radius)
       ELSE IF (ASSOCIATED(qs_kind%gth_potential)) THEN
          CALL get_potential(potential=qs_kind%gth_potential,&
               core_charge_radius=core_charge_radius)
       ELSE
          core_charge_radius = 0.0_dp
       END IF
    END IF
    IF (PRESENT(core_charge)) THEN
       IF (ASSOCIATED(qs_kind%all_potential)) THEN
          CALL get_potential(potential=qs_kind%all_potential,&
               zeff=core_charge)
       ELSE IF (ASSOCIATED(qs_kind%gth_potential)) THEN
          CALL get_potential(potential=qs_kind%gth_potential,&
               zeff=core_charge)
       ELSE
          core_charge = 0.0_dp
       END IF
    END IF

    IF (PRESENT(zeff)) THEN
       IF (ASSOCIATED(qs_kind%all_potential)) THEN
          CALL get_potential(potential=qs_kind%all_potential,zeff=zeff)
       ELSE IF (ASSOCIATED(qs_kind%gth_potential)) THEN
          CALL get_potential(potential=qs_kind%gth_potential,zeff=zeff)
       ELSE
          zeff = 0.0_dp
       END IF
    END IF

    IF (PRESENT(soft_basis_set)) soft_basis_set => qs_kind%soft_basis_set
    IF (PRESENT(hard_basis_set)) hard_basis_set => qs_kind%hard_basis_set

    IF (PRESENT(paw_proj_set))   paw_proj_set => qs_kind%paw_proj_set
    IF (PRESENT(paw_atom)) paw_atom = qs_kind%paw_atom
    IF (PRESENT(gpw_type_forced)) gpw_type_forced = qs_kind%gpw_type_forced
    IF (PRESENT(hard_radius)) hard_radius = qs_kind%hard_radius
    IF (PRESENT(hard0_radius)) hard0_radius = qs_kind%hard0_radius
    IF (PRESENT(max_rad_local)) max_rad_local = qs_kind%max_rad_local
    IF (PRESENT(harmonics))  harmonics => qs_kind%harmonics
    IF (PRESENT(max_s_harm)) THEN
       IF(ASSOCIATED(qs_kind%harmonics)) THEN
          max_s_harm = qs_kind%harmonics%max_s_harm
       ELSE
          max_s_harm = 0
       END IF
    END IF
    IF (PRESENT(max_iso_not0)) THEN
       IF(ASSOCIATED(qs_kind%harmonics)) THEN
          max_iso_not0 = qs_kind%harmonics%max_iso_not0
       ELSE
          max_iso_not0 = 0
       END IF
    END IF
    IF (PRESENT(grid_atom)) grid_atom => qs_kind%grid_atom
    IF (PRESENT(ngrid_ang)) ngrid_ang = qs_kind%ngrid_ang
    IF (PRESENT(ngrid_rad)) ngrid_rad = qs_kind%ngrid_rad
    IF (PRESENT(lmax_rho0)) lmax_rho0 = qs_kind%lmax_rho0
    IF (PRESENT(ghost)) ghost = qs_kind%ghost
    IF (PRESENT(dft_plus_u_atom)) dft_plus_u_atom = ASSOCIATED(qs_kind%dft_plus_u)
    IF (PRESENT(l_of_dft_plus_u)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          l_of_dft_plus_u = qs_kind%dft_plus_u%l
       ELSE
          l_of_dft_plus_u = -1
       END IF
    END IF
    IF (PRESENT(u_minus_j)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          u_minus_j = qs_kind%dft_plus_u%u_minus_j
       ELSE
          u_minus_j = 0.0_dp
       END IF
    END IF
    IF (PRESENT(u_minus_j_target)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          u_minus_j_target = qs_kind%dft_plus_u%u_minus_j_target
       ELSE
          u_minus_j_target = 0.0_dp
       END IF
    END IF
    IF (PRESENT(init_u_ramping_each_scf)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          init_u_ramping_each_scf = qs_kind%dft_plus_u%init_u_ramping_each_scf
       ELSE
          init_u_ramping_each_scf = .FALSE.
       END IF
    END IF
    IF (PRESENT(u_ramping)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          u_ramping = qs_kind%dft_plus_u%u_ramping
       ELSE
          u_ramping = 0.0_dp
       END IF
    END IF
    IF (PRESENT(eps_u_ramping)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          eps_u_ramping = qs_kind%dft_plus_u%eps_u_ramping
       ELSE
          eps_u_ramping = 1.0E-5_dp
       END IF
    END IF
    IF (PRESENT(orbitals)) THEN
       NULLIFY (orbitals)
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          IF (ASSOCIATED(qs_kind%dft_plus_u%orbitals)) THEN
             orbitals => qs_kind%dft_plus_u%orbitals
          END IF
       END IF
    END IF
    IF (PRESENT(eps_scf)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          eps_scf = qs_kind%dft_plus_u%eps_scf
       ELSE
          eps_scf = 1.0E30_dp
       END IF
    END IF
    IF (PRESENT(max_scf)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          max_scf = qs_kind%dft_plus_u%max_scf
       ELSE
          max_scf = -1
       END IF
    END IF
    IF (PRESENT(smear)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          smear = qs_kind%dft_plus_u%smear
       ELSE
          smear = .FALSE.
       END IF
    END IF
    IF (PRESENT(dispersion)) dispersion => qs_kind%dispersion
    IF (PRESENT(bs_occupation)) bs_occupation=qs_kind%bs_occupation
    IF (PRESENT(addel))  addel  => qs_kind%addel
    IF (PRESENT(laddel)) laddel => qs_kind%laddel
    IF (PRESENT(naddel)) naddel => qs_kind%naddel

    IF (PRESENT(no_optimize)) no_optimize=qs_kind%no_optimize

    IF (PRESENT(reltmat)) reltmat => qs_kind%reltmat

    IF (PRESENT(mao))  mao  = qs_kind%mao

  END SUBROUTINE get_qs_kind

! *****************************************************************************
!> \brief Get attributes of an atomic kind set.
!> \param qs_kind_set ...
!> \param all_potential_present ...
!> \param tnadd_potential_present ...
!> \param gth_potential_present ...
!> \param paw_atom_present ...
!> \param dft_plus_u_atom_present ...
!> \param maxcgf ...
!> \param maxsgf ...
!> \param maxco ...
!> \param maxco_proj ...
!> \param maxgtops ...
!> \param maxlgto ...
!> \param maxlprj ...
!> \param maxnset ...
!> \param maxsgf_set ...
!> \param ncgf ...
!> \param npgf ...
!> \param nset ...
!> \param nsgf ...
!> \param nshell ...
!> \param maxpol ...
!> \param maxlppl ...
!> \param maxlppnl ...
!> \param maxppnl ...
!> \param nelectron ...
!> \param maxder ...
!> \param max_ngrid_rad ...
!> \param max_sph_harm ...
!> \param maxg_iso_not0 ...
!> \param lmax_rho0 ...
!> \param basis_type ...
! *****************************************************************************
  SUBROUTINE get_qs_kind_set(qs_kind_set,&
             all_potential_present, tnadd_potential_present, gth_potential_present,&
             paw_atom_present, dft_plus_u_atom_present,&
             maxcgf, maxsgf, maxco, maxco_proj, maxgtops, maxlgto, maxlprj, maxnset, maxsgf_set,&
             ncgf, npgf, nset, nsgf, nshell, maxpol, maxlppl, maxlppnl, maxppnl,&
             nelectron, maxder, max_ngrid_rad, max_sph_harm, maxg_iso_not0, lmax_rho0,&
             basis_type)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    LOGICAL, INTENT(OUT), OPTIONAL :: all_potential_present, &
      tnadd_potential_present, gth_potential_present, paw_atom_present, &
      dft_plus_u_atom_present
    INTEGER, INTENT(OUT), OPTIONAL :: maxcgf, maxsgf, maxco, maxco_proj, &
      maxgtops, maxlgto, maxlprj, maxnset, maxsgf_set, ncgf, npgf, nset, &
      nsgf, nshell, maxpol, maxlppl, maxlppnl, maxppnl, nelectron
    INTEGER, INTENT(IN), OPTIONAL            :: maxder
    INTEGER, INTENT(OUT), OPTIONAL           :: max_ngrid_rad, max_sph_harm, &
                                                maxg_iso_not0, lmax_rho0
    CHARACTER(len=*), OPTIONAL               :: basis_type

    CHARACTER(len=*), PARAMETER :: routineN = 'get_qs_kind_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_basis_type
    INTEGER                                  :: ikind, imax, lmax_rho0_kind, &
                                                max_iso_not0, max_s_harm, n, &
                                                ngrid_rad, nkind
    LOGICAL                                  :: dft_plus_u_atom, paw_atom
    REAL(KIND=dp)                            :: zeff, zeff_correction
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: tmp_basis_set
    TYPE(local_potential_type), POINTER      :: tnadd_potential
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    TYPE(qs_kind_type), POINTER              :: qs_kind

    IF( PRESENT(basis_type) ) THEN
      my_basis_type = basis_type
    ELSE
      my_basis_type = "ORB"
    END IF

    IF (ASSOCIATED(qs_kind_set)) THEN

       IF (PRESENT(maxcgf)) maxcgf = 0
       IF (PRESENT(maxco)) maxco = 0
       IF (PRESENT(maxco_proj)) maxco_proj = 0
       IF (PRESENT(maxg_iso_not0)) maxg_iso_not0 = 0
       IF (PRESENT(maxgtops)) maxgtops = 0
       IF (PRESENT(maxlgto)) maxlgto = -1
       IF (PRESENT(maxlppl)) maxlppl = -1
       IF (PRESENT(maxlppnl)) maxlppnl = -1
       IF (PRESENT(maxpol)) maxpol = -1
       IF (PRESENT(maxlprj)) maxlprj = -1
       IF (PRESENT(maxnset)) maxnset = 0
       IF (PRESENT(maxppnl)) maxppnl = 0
       IF (PRESENT(maxsgf)) maxsgf = 0
       IF (PRESENT(maxsgf_set)) maxsgf_set = 0
       IF (PRESENT(ncgf)) ncgf = 0
       IF (PRESENT(nelectron)) nelectron = 0
       IF (PRESENT(npgf)) npgf = 0
       IF (PRESENT(nset)) nset = 0
       IF (PRESENT(nsgf)) nsgf = 0
       IF (PRESENT(nshell)) nshell = 0
       IF (PRESENT(all_potential_present)) all_potential_present = .FALSE.
       IF (PRESENT(tnadd_potential_present)) tnadd_potential_present = .FALSE.
       IF (PRESENT(gth_potential_present)) gth_potential_present = .FALSE.
       IF (PRESENT(paw_atom_present)) paw_atom_present = .FALSE.
       IF (PRESENT(max_ngrid_rad))  max_ngrid_rad = 0
       IF (PRESENT(max_sph_harm))   max_sph_harm = 0
       IF (PRESENT(lmax_rho0))      lmax_rho0 = 0

       nkind = SIZE(qs_kind_set)
       DO ikind=1,nkind
          qs_kind => qs_kind_set(ikind)
          CALL get_qs_kind(qs_kind=qs_kind,&
               all_potential=all_potential,&
               tnadd_potential=tnadd_potential,&
               gth_potential=gth_potential,&
               paw_proj_set=paw_proj_set,&
               dftb_parameter=dftb_parameter,&
               ngrid_rad=ngrid_rad,&
               max_s_harm=max_s_harm, &
               max_iso_not0=max_iso_not0,&
               paw_atom=paw_atom,&
               dft_plus_u_atom=dft_plus_u_atom,&
               lmax_rho0=lmax_rho0_kind)

          IF (PRESENT(maxlppl).AND.ASSOCIATED(gth_potential)) THEN
             CALL get_potential(potential=gth_potential,nexp_ppl=n)
             maxlppl = MAX(maxlppl,2*(n - 1))
          END IF

          IF (PRESENT(maxlppnl).AND.ASSOCIATED(gth_potential)) THEN
             CALL get_potential(potential=gth_potential,lprj_ppnl_max=imax)
             maxlppnl = MAX(maxlppnl,imax)
          END IF

          IF (PRESENT(maxpol).AND.ASSOCIATED(tnadd_potential)) THEN
             CALL get_potential(potential=tnadd_potential,npol=n)
             maxpol = MAX(maxpol,2*(n - 1))
          END IF

          IF (PRESENT(maxco_proj).AND.ASSOCIATED(paw_proj_set)) THEN
             CALL get_paw_proj_set(paw_proj_set=paw_proj_set,ncgauprj=imax)
             maxco_proj = MAX(maxco_proj,imax)
          ENDIF

          IF (PRESENT(maxlprj).AND.ASSOCIATED(paw_proj_set)) THEN
             CALL get_paw_proj_set(paw_proj_set=paw_proj_set,maxl=imax)
             maxlprj = MAX(maxlprj,imax)
          ENDIF

          IF (PRESENT(maxppnl).AND.ASSOCIATED(gth_potential)) THEN
             CALL get_potential(potential=gth_potential,nppnl=imax)
             maxppnl = MAX(maxppnl,imax)
          END IF

          CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis_set,&
               basis_type=my_basis_type)

          IF (PRESENT(maxcgf)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,ncgf=imax)
                maxcgf = MAX(maxcgf,imax)
             ELSE IF(ASSOCIATED(qs_kind%dftb_parameter)) THEN
                CALL get_dftb_atom_param(dftb_parameter=dftb_parameter,lmax=imax)
                imax=((imax+1)*(imax+2)*(imax+3))/6
                maxcgf = MAX(maxcgf,imax)
             END IF
          END IF

          IF (PRESENT(maxco)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                IF (PRESENT(maxder)) THEN
                   CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,&
                        maxco=imax,maxder=maxder)
                ELSE
                   CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,maxco=imax)
                END IF
                maxco = MAX(maxco,imax)
             END IF
             IF (ASSOCIATED(gth_potential)) THEN
                CALL get_potential(potential=gth_potential,lprj_ppnl_max=imax)
                maxco = MAX(maxco,ncoset(imax))
             END IF
          END IF

          IF (PRESENT(maxgtops)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,maxso=imax,nset=n)
                maxgtops = MAX(maxgtops, n*imax)
             END IF
          END IF

          IF (PRESENT(maxlgto)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,maxl=imax)
                maxlgto = MAX(maxlgto,imax)
             ELSE IF(ASSOCIATED(qs_kind%dftb_parameter)) THEN
                CALL get_dftb_atom_param(dftb_parameter=dftb_parameter,lmax=imax)
                maxlgto = MAX(maxlgto,imax)
             END IF
          END IF

          IF (PRESENT(maxnset)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,nset=n)
                maxnset = MAX(maxnset, n)
             END IF
          END IF

          IF (PRESENT(maxsgf)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,nsgf=imax)
                maxsgf = MAX(maxsgf,imax)
             END IF
          END IF

          IF (PRESENT(maxsgf_set)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,maxsgf_set=imax)
                maxsgf_set = MAX(maxsgf_set,imax)
             END IF
          END IF

          IF (PRESENT(ncgf)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,ncgf=n)
                ncgf = ncgf + n*qs_kind_set(ikind)%natom
             ELSE IF(ASSOCIATED(qs_kind%dftb_parameter)) THEN
                CALL get_dftb_atom_param(dftb_parameter=dftb_parameter,lmax=imax)
                n=((imax+1)*(imax+2)*(imax+3))/6
                ncgf = ncgf + n*qs_kind_set(ikind)%natom
             END IF
          END IF

          IF (PRESENT(npgf)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,npgf_sum=n)
                npgf = npgf + n*qs_kind_set(ikind)%natom
             END IF
          END IF

          IF (PRESENT(nset)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,nset=n)
                nset = nset + n*qs_kind_set(ikind)%natom
             END IF
          END IF

          IF (PRESENT(nsgf)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,nsgf=n)
                nsgf = nsgf + n*qs_kind_set(ikind)%natom
             ELSE IF(ASSOCIATED(qs_kind%dftb_parameter)) THEN
                CALL get_dftb_atom_param(dftb_parameter=dftb_parameter,natorb=n)
                nsgf = nsgf + n*qs_kind_set(ikind)%natom
             END IF
          END IF

          IF (PRESENT(nshell)) THEN
             IF(ASSOCIATED(tmp_basis_set)) THEN
                CALL get_gto_basis_set(gto_basis_set=tmp_basis_set,nshell_sum=n)
                nshell = nshell + n*qs_kind_set(ikind)%natom
             ELSE IF(ASSOCIATED(qs_kind%dftb_parameter)) THEN
                CALL get_dftb_atom_param(dftb_parameter=dftb_parameter,lmax=n)
                nshell = nshell + (n+1)*qs_kind_set(ikind)%natom
             END IF
          END IF

          IF (PRESENT(nelectron)) THEN
             IF (ASSOCIATED(qs_kind%all_potential)) THEN
                CALL get_potential(potential=qs_kind%all_potential,&
                     zeff=zeff,zeff_correction=zeff_correction)
             ELSE IF (ASSOCIATED(qs_kind%gth_potential)) THEN
                CALL get_potential(potential=qs_kind%gth_potential,&
                     zeff=zeff,zeff_correction=zeff_correction)
             ELSE
                zeff = 0.0_dp
                zeff_correction = 0.0_dp
             END IF
             nelectron = nelectron + qs_kind_set(ikind)%natom*NINT(zeff-zeff_correction)
          END IF

          IF (PRESENT(all_potential_present)) THEN
             IF (ASSOCIATED(all_potential)) THEN
                all_potential_present = .TRUE.
             END IF
          END IF

          IF (PRESENT(tnadd_potential_present)) THEN
             IF (ASSOCIATED(tnadd_potential)) THEN
                tnadd_potential_present = .TRUE.
             END IF
          END IF

          IF (PRESENT(gth_potential_present)) THEN
             IF (ASSOCIATED(gth_potential)) THEN
                gth_potential_present = .TRUE.
             END IF
          END IF

          IF (PRESENT(paw_atom_present)) THEN
             IF (paw_atom) THEN
                paw_atom_present = .TRUE.
             END IF
          END IF

          IF (PRESENT(dft_plus_u_atom_present)) THEN
             IF (dft_plus_u_atom) THEN
                dft_plus_u_atom_present = .TRUE.
             END IF
          END IF

          IF(PRESENT(max_ngrid_rad)) THEN
             max_ngrid_rad = MAX(max_ngrid_rad,ngrid_rad)
          ENDIF

          IF(PRESENT(max_sph_harm)) THEN
             max_sph_harm = MAX(max_sph_harm,max_s_harm)
          ENDIF

          IF(PRESENT(maxg_iso_not0)) THEN
             maxg_iso_not0 = MAX(maxg_iso_not0,max_iso_not0)
          END IF

          IF( PRESENT(lmax_rho0) ) THEN
             lmax_rho0 = MAX(lmax_rho0,lmax_rho0_kind)
          END IF

       END DO
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer qs_kind_set is not associated")
    END IF

  END SUBROUTINE get_qs_kind_set

! *****************************************************************************
!> \brief Initialise an atomic kind data set.
!> \param qs_kind ...
!> \author Creation (11.01.2002,MK)
!>                20.09.2002 adapted for pol/kg use, gtb
! *****************************************************************************
  SUBROUTINE init_qs_kind(qs_kind)
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CHARACTER(len=*), PARAMETER :: routineN = 'init_qs_kind', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: basis_type
    INTEGER                                  :: handle, i
    TYPE(gto_basis_set_type), POINTER        :: tmp_basis_set

    CALL timeset(routineN,handle)

    CPASSERT(ASSOCIATED(qs_kind))

    IF (ASSOCIATED(qs_kind%gth_potential)) THEN
       CALL init_potential(qs_kind%gth_potential)
    END IF

    DO i=1,SIZE(qs_kind%basis_sets,1)
       NULLIFY(tmp_basis_set)
       CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis_set,&
            inumbas=i,basis_type=basis_type)
       IF(basis_type=="") CYCLE
       IF(basis_type=="AUX") THEN
          IF(tmp_basis_set%norm_type < 0) tmp_basis_set%norm_type = 1
          CALL init_aux_basis_set(tmp_basis_set)
       ELSE
          IF(tmp_basis_set%norm_type < 0) tmp_basis_set%norm_type = 2
          CALL init_orb_basis_set(tmp_basis_set)
       END IF
    END DO

    CALL timestop(handle)

  END SUBROUTINE init_qs_kind

! *****************************************************************************
!> \brief Initialise an atomic kind set data set.
!> \param qs_kind_set ...
!> \author - Creation (17.01.2002,MK)
!>      - 20.09.2002 para_env passed (gt)
! *****************************************************************************
  SUBROUTINE init_qs_kind_set(qs_kind_set)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CHARACTER(len=*), PARAMETER :: routineN = 'init_qs_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    IF (.NOT. ASSOCIATED(qs_kind_set)) &
       STOP "init_qs_kind_set: The pointer qs_kind_set is not associated"

    DO ikind=1, SIZE(qs_kind_set)
       qs_kind => qs_kind_set(ikind)
       CALL init_qs_kind(qs_kind)
    END DO

    CALL timestop(handle)

  END SUBROUTINE init_qs_kind_set

! *****************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param qs_control ...
!> \param force_env_section ...
! *****************************************************************************
  SUBROUTINE init_gapw_basis_set(qs_kind_set,qs_control,force_env_section)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_control_type), POINTER           :: qs_control
    TYPE(section_vals_type), POINTER         :: force_env_section

    CHARACTER(len=*), PARAMETER :: routineN = 'init_gapw_basis_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, nkind
    LOGICAL                                  :: gpw, paw_atom
    REAL(dp)                                 :: max_rad_local_type, rc
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(paw_proj_set_type), POINTER         :: paw_proj
    TYPE(qs_kind_type), POINTER              :: qs_kind

    IF (ASSOCIATED(qs_kind_set)) THEN

       qs_control%gapw_control%non_paw_atoms = .FALSE.
       nkind = SIZE(qs_kind_set)

       DO ikind=1,nkind

          qs_kind => qs_kind_set(ikind)

          CALL allocate_gto_basis_set(qs_kind%soft_basis_set)

          ! The hard_basis points to the orb_basis
          CALL get_qs_kind(qs_kind=qs_kind,basis_set=orb_basis)
          qs_kind%hard_basis_set => orb_basis

          CALL get_qs_kind(qs_kind=qs_kind,hard_radius=rc,&
               max_rad_local=max_rad_local_type,gpw_type_forced=gpw)

          CALL create_soft_basis(orb_basis,qs_kind%soft_basis_set,&
               qs_control%gapw_control%eps_fit,rc,paw_atom,&
               qs_control%gapw_control%force_paw,gpw)

          CALL set_qs_kind(qs_kind=qs_kind,paw_atom=paw_atom)

          IF (paw_atom) THEN
             CALL allocate_paw_proj_set(qs_kind%paw_proj_set)
             CALL get_qs_kind(qs_kind=qs_kind,paw_proj_set=paw_proj)

             CALL projectors(paw_proj,orb_basis,rc,qs_control,max_rad_local_type,&
                  force_env_section)
          ELSE
             qs_control%gapw_control%non_paw_atoms = .TRUE.
          END IF

          ! grid_atom and harmonics are allocated even if NOT PAW_ATOM
          NULLIFY(qs_kind%grid_atom,qs_kind%harmonics)
          CALL allocate_grid_atom(qs_kind%grid_atom)
          CALL allocate_harmonics_atom(qs_kind%harmonics)
       END DO

       IF(qs_control%gapw_control%non_paw_atoms) THEN
          qs_control%gapw_control%nopaw_as_gpw = .TRUE.
       ELSE
          qs_control%gapw_control%nopaw_as_gpw = .FALSE.
       END IF
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer qs_kind_set is not associated")
    END IF

  END SUBROUTINE init_gapw_basis_set

! *****************************************************************************
!> \brief Read an atomic kind data set from the input file.
!> \param qs_kind ...
!> \param kind_section ...
!> \param para_env ...
!> \param force_env_section ...
!> \param no_fail ...
!> \param method_id ...
!> \par History
!>      - Creation (09.02.2002,MK)
!>      - 20.09.2002,gt: adapted for POL/KG use (elp_potential)
!>      - 05.03.2010: split elp_potential into fist_potential and kg_potential
! *****************************************************************************
  SUBROUTINE read_qs_kind(qs_kind,kind_section,para_env,force_env_section,&
       no_fail, method_id)

    TYPE(qs_kind_type), INTENT(INOUT)        :: qs_kind
    TYPE(section_vals_type), POINTER         :: kind_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    LOGICAL, INTENT(IN)                      :: no_fail
    INTEGER, INTENT(IN)                      :: method_id

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_qs_kind', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: maxbas = 20

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_string_length) :: akind_name, aux_basis_set_name, &
      aux_fit_basis_set_name, basis_type, keyword, kgpot_name, kgpot_type, &
      lri_basis_set_name, orb_basis_set_name, potential_name, potential_type, &
      ri_aux_basis_set_name, tmp
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    CHARACTER(LEN=default_string_length), &
      DIMENSION(maxbas)                      :: basis_set_name, basis_set_type
    INTEGER                                  :: handle, i, i_rep, ipos, j, &
                                                k_rep, l, m, n_rep, nb_rep, &
                                                nexp, ngauss, norbitals, &
                                                stat, z
    INTEGER, DIMENSION(:), POINTER           :: add_el, elec_conf, orbitals
    LOGICAL :: check, explicit, explicit_basis, explicit_kgpot, &
      explicit_potential, failure, nobasis, section_enabled, &
      subsection_enabled
    REAL(KIND=dp)                            :: zeff_correction
    TYPE(gto_basis_set_type), POINTER        :: sup_basis_set, tmp_basis_set
    TYPE(section_vals_type), POINTER :: basis_section, bs_section, &
      dft_plus_u_section, dft_section, enforce_occupation_section, &
      kgpot_section, potential_section, spin_section

    CALL timeset(routineN,handle)

    NULLIFY (elec_conf)

    failure = .FALSE.
    basis_set_name(:) = ""
    basis_set_type(:) = ""
    aux_basis_set_name = ""
    ri_aux_basis_set_name = ""
    orb_basis_set_name = ""
    potential_name = ""
    potential_type = ""
    kgpot_name = ""
    kgpot_type = ""
    z = -1
    zeff_correction = 0.0_dp
    dft_section => section_vals_get_subs_vals(force_env_section,"DFT")
    CALL section_vals_get(kind_section,n_repetition=n_rep)
    k_rep=-1
    akind_name=qs_kind%name
    CALL uppercase(akind_name)
    ! First we use the atom_name to find out the proper KIND section
    DO i_rep=1,n_rep
       CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
            c_val=keyword,i_rep_section=i_rep)
       CALL uppercase(keyword)
       IF (keyword==akind_name) THEN
          k_rep=i_rep
          EXIT
       END IF
    END DO
    ! The search for the KIND section failed.. check for a QM/MM link atom
    IF (k_rep<1) THEN
       ipos = INDEX(qs_kind%name,"_")
       IF (((ipos == 2).OR.(ipos == 3)).AND.(INDEX(qs_kind%name,"_ghost")==0)) THEN
          ! If the atm_name could not match any KIND section it maybe be a QM/MM link atom.
          ! ghost atoms will be treated differently.
          akind_name = qs_kind%name(1:ipos-1)
          CALL uppercase(akind_name)
          DO i_rep=1,n_rep
             CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
                  c_val=keyword,i_rep_section=i_rep)
             CALL uppercase(keyword)
             IF (keyword==akind_name) THEN
                k_rep=i_rep
                EXIT
             END IF
          END DO
       END IF
    END IF
    ! The search for the KIND section failed.. check element_symbol
    IF (k_rep<1) THEN
       ! If it's not a link atom let's check for the element and map
       ! the KIND section to the element.
       element_symbol = qs_kind%element_symbol(1:2)
       CALL uppercase(element_symbol)
       DO i_rep=1,n_rep
          CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
               c_val=keyword,i_rep_section=i_rep)
          CALL uppercase(keyword)
          IF (keyword==element_symbol) THEN
             k_rep=i_rep
             EXIT
          END IF
       END DO
    END IF
    ! In case it should not really match any possible KIND section
    ! let's look if a default one is defined..
    IF (k_rep<1) THEN
       DO i_rep=1,n_rep
          CALL section_vals_val_get(kind_section,"_SECTION_PARAMETERS_",&
               c_val=keyword,i_rep_section=i_rep)
          CALL uppercase(keyword)
          IF (keyword=="DEFAULT") THEN
             k_rep=i_rep
             EXIT
          END IF
       END DO
    END IF
    IF (k_rep<0.AND.(.NOT.no_fail)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "No &KIND section was possible to associate to the atomic kind <"//&
            TRIM(akind_name)//">. The KIND section were also scanned for the"//&
            " corresponding element <"//TRIM(qs_kind%element_symbol)//"> "//&
            " and for the DEFAULT section but no match was found. Check your input file!",para_env)
    END IF
    ! Retrieve information on element
    CALL get_ptable_info(qs_kind%element_symbol, ielement=z)

    ! Normal parsing of the KIND section
    IF (k_rep>0) THEN
       ! new style basis set input
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
                             keyword_name="BASIS_SET",&
                             explicit=explicit,&
                             n_rep_val=nb_rep)
       IF(.NOT.explicit) nb_rep = 0
       CPASSERT(nb_rep<=maxbas)
       DO i=1,nb_rep
          CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
               keyword_name="BASIS_SET",i_rep_val=i,c_vals=tmpstringlist)
          IF(SIZE(tmpstringlist)==1) THEN
             ! default is orbital type
             basis_set_type(i) = "ORB"
             basis_set_name(i) = tmpstringlist(1)
          ELSEIF(SIZE(tmpstringlist)==2) THEN
             basis_set_type(i) = tmpstringlist(1)
             basis_set_name(i) = tmpstringlist(2)
          ELSE
             CPASSERT(.FALSE.)
          END IF
       END DO

       ! Deprecated keywords
       DO i=1,nb_rep
          IF(basis_set_type(i)=="ORB") THEN
             orb_basis_set_name = basis_set_name(i)
             EXIT
          END IF
       END DO
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="AUX_FIT_BASIS_SET", c_val=aux_fit_basis_set_name)
       IF(aux_fit_basis_set_name /= "") THEN
          nb_rep = nb_rep + 1
          CPASSERT(nb_rep<=maxbas)
          basis_set_type(nb_rep) = "AUX_FIT"
          basis_set_name(nb_rep) = aux_fit_basis_set_name
       END IF
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="AUX_BASIS_SET", c_val=aux_basis_set_name)
       IF(aux_basis_set_name /= "") THEN
          nb_rep = nb_rep + 1
          CPASSERT(nb_rep<=maxbas)
          basis_set_type(nb_rep) = "AUX"
          basis_set_name(nb_rep) = aux_basis_set_name
       END IF
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="LRI_BASIS_SET", c_val=lri_basis_set_name)
       IF(lri_basis_set_name /= "") THEN
          nb_rep = nb_rep + 1
          CPASSERT(nb_rep<=maxbas)
          basis_set_type(nb_rep) = "LRI"
          basis_set_name(nb_rep) = lri_basis_set_name
       END IF
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="RI_AUX_BASIS_SET", c_val=ri_aux_basis_set_name)
       IF(ri_aux_basis_set_name /= "") THEN
          nb_rep = nb_rep + 1
          CPASSERT(nb_rep<=maxbas)
          basis_set_type(nb_rep) = "RI_AUX"
          basis_set_name(nb_rep) = ri_aux_basis_set_name
       END IF
       ! end of deprecated input

       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="ELEC_CONF",n_rep_val=i)
       IF (i>0) THEN
          CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
               keyword_name="ELEC_CONF",i_vals=elec_conf)
          CALL set_qs_kind(qs_kind,elec_conf=elec_conf)
       ENDIF
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="CORE_CORRECTION",r_val=zeff_correction)
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="POTENTIAL",c_val=potential_name)
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="KG_POTENTIAL",c_val=kgpot_name)
       ! assign atom dependent defaults, only H special case
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,n_rep_val=i,&
            keyword_name="HARD_EXP_RADIUS")
       IF (i==0) THEN
          IF (z==1) THEN
             qs_kind%hard_radius=1.2_dp
          ELSE
             qs_kind%hard_radius=0.8_dp*bohr
          ENDIF
       ELSE
          CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
               keyword_name="HARD_EXP_RADIUS",r_val=qs_kind%hard_radius)
       ENDIF

       ! assign atom dependent defaults, only H special case
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,n_rep_val=i,&
            keyword_name="RHO0_EXP_RADIUS")
       IF (i==0) THEN
          qs_kind%hard0_radius=qs_kind%hard_radius
       ELSE
          CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
               keyword_name="RHO0_EXP_RADIUS",r_val=qs_kind%hard0_radius)
       ENDIF
       CALL cp_assert(qs_kind%hard_radius>=qs_kind%hard0_radius,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "rc0 should be <= rc",failure)

       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="MAX_RAD_LOCAL",r_val=qs_kind%max_rad_local)
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="LEBEDEV_GRID",i_val=qs_kind%ngrid_ang)
       CALL cp_assert(qs_kind%ngrid_ang > 0,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "# point lebedev grid < 0",failure)
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="RADIAL_GRID",i_val=qs_kind%ngrid_rad)
       CALL cp_assert(qs_kind%ngrid_rad > 0,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "# point radial grid < 0",failure)
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="GPW_TYPE",l_val=qs_kind%gpw_type_forced)
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="GHOST",l_val=qs_kind%ghost)
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
            keyword_name="NO_OPTIMIZE",l_val=qs_kind%no_optimize)

       ! DFTB3 param
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
                keyword_name="DFTB3_PARAM",r_val=qs_kind%dudq_dftb3)

       ! MAOS
       CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
                keyword_name="MAO",i_val=qs_kind%mao)

       ! Read the BS subsection of the current atomic kind, if enabled
       NULLIFY(bs_section)
       bs_section => section_vals_get_subs_vals(kind_section,"BS",&
            i_rep_section=k_rep)
       section_enabled = .FALSE.
       CALL section_vals_val_get(bs_section,"_SECTION_PARAMETERS_",&
            l_val=section_enabled)
       IF (section_enabled) THEN
          qs_kind%bs_occupation = .TRUE.
          !Alpha spin
          NULLIFY(spin_section)
          spin_section => section_vals_get_subs_vals(bs_section,"ALPHA")
          CALL section_vals_get(spin_section, explicit=explicit)
          IF(explicit) THEN
              NULLIFY(add_el)
              CALL section_vals_val_get(spin_section,&
                   keyword_name="NEL",i_vals=add_el)
              CPASSERT(ASSOCIATED(add_el))
              ALLOCATE(qs_kind%addel(SIZE(add_el),2),STAT=stat)
              qs_kind%addel =0
              qs_kind%addel(1:SIZE(add_el),1) = add_el(1:SIZE(add_el))
              NULLIFY(add_el)
              CALL section_vals_val_get(spin_section,&
                   keyword_name="L",i_vals=add_el)
              CPASSERT(ASSOCIATED(add_el))
              CPASSERT(SIZE(add_el)==SIZE( qs_kind%addel,1))
              ALLOCATE(qs_kind%laddel(SIZE(add_el),2),STAT=stat)
              qs_kind%laddel =0
              qs_kind%laddel(1:SIZE(add_el),1) = add_el(1:SIZE(add_el))
              ALLOCATE(qs_kind%naddel(SIZE(add_el),2),STAT=stat)
              qs_kind%naddel =0
              NULLIFY(add_el)
              CALL section_vals_val_get(spin_section,&
                 keyword_name="N",n_rep_val=i)
              IF(i>0) THEN
                CALL section_vals_val_get(spin_section,&
                     keyword_name="N",i_vals=add_el)
                IF(SIZE(add_el)==SIZE( qs_kind%addel,1)) THEN
                  qs_kind%naddel(1:SIZE(add_el),1) = add_el(1:SIZE(add_el))
                END IF
              END IF
          END IF
          !Beta spin
          NULLIFY(spin_section)
          spin_section => section_vals_get_subs_vals(bs_section,"BETA")
          CALL section_vals_get(spin_section, explicit=explicit)
          IF(explicit) THEN

            NULLIFY(add_el)
            CALL section_vals_val_get(spin_section,&
                 keyword_name="NEL",i_vals=add_el)
            CPASSERT(SIZE(add_el)==SIZE( qs_kind%addel,1))
            qs_kind%addel(1:SIZE(add_el),2) = add_el(1:SIZE(add_el))
            qs_kind%addel(:,:) = qs_kind%addel(:,:)
            NULLIFY(add_el)
            CALL section_vals_val_get(spin_section,&
                 keyword_name="L",i_vals=add_el)
            CPASSERT(SIZE(add_el)==SIZE( qs_kind%addel,1))
            qs_kind%laddel(1:SIZE(add_el),2) = add_el(1:SIZE(add_el))

            CALL section_vals_val_get(spin_section,&
                keyword_name="N",n_rep_val=i)
            IF(i>0) THEN
              NULLIFY(add_el)
              CALL section_vals_val_get(spin_section,&
                    keyword_name="N",i_vals=add_el)
              IF(SIZE(add_el)==SIZE( qs_kind%addel,1)) THEN
                qs_kind%naddel(1:SIZE(add_el),2) = add_el(1:SIZE(add_el))
              END IF
            END IF
          END IF
       END IF

       ! Read the DFT+U subsection of the current atomic kind, if enabled

       NULLIFY (dft_plus_u_section)
       dft_plus_u_section => section_vals_get_subs_vals(kind_section,&
                                                        subsection_name="DFT_PLUS_U",&
                                                        i_rep_section=k_rep)
       section_enabled = .FALSE.
       CALL section_vals_val_get(dft_plus_u_section,&
                                 keyword_name="_SECTION_PARAMETERS_",&
                                 l_val=section_enabled)
       IF (section_enabled) THEN
          ALLOCATE (qs_kind%dft_plus_u,STAT=stat)
          CPASSERT(stat==0)
          NULLIFY (qs_kind%dft_plus_u%orbitals)
          CALL section_vals_val_get(dft_plus_u_section,&
                                    keyword_name="L",&
                                    i_val=l)
          qs_kind%dft_plus_u%l = l
          CALL cp_assert((l >= 0),cp_failure_level,&
                         cp_assertion_failed,routineP,&
                         "DFT+U| Invalid orbital angular momentum quantum number specified: l < 0",&
                         failure)
          CALL section_vals_val_get(dft_plus_u_section,&
                                    keyword_name="U_MINUS_J",&
                                    r_val=qs_kind%dft_plus_u%u_minus_j_target)
          CALL section_vals_val_get(dft_plus_u_section,&
                                    keyword_name="U_RAMPING",&
                                    r_val=qs_kind%dft_plus_u%u_ramping)
          CALL section_vals_val_get(dft_plus_u_section,&
                                    keyword_name="INIT_U_RAMPING_EACH_SCF",&
                                    l_val=qs_kind%dft_plus_u%init_u_ramping_each_scf)
          IF (qs_kind%dft_plus_u%u_ramping > 0.0_dp) THEN
             qs_kind%dft_plus_u%u_minus_j = 0.0_dp
          ELSE
             qs_kind%dft_plus_u%u_minus_j = qs_kind%dft_plus_u%u_minus_j_target
          END IF
          CALL section_vals_val_get(dft_plus_u_section,&
                                    keyword_name="EPS_U_RAMPING",&
                                    r_val=qs_kind%dft_plus_u%eps_u_ramping)

          NULLIFY (enforce_occupation_section)
          enforce_occupation_section => section_vals_get_subs_vals(dft_plus_u_section,&
                                                                   subsection_name="ENFORCE_OCCUPATION")
          subsection_enabled = .FALSE.
          CALL section_vals_val_get(enforce_occupation_section,&
                                    keyword_name="_SECTION_PARAMETERS_",&
                                    l_val=subsection_enabled)
          IF (subsection_enabled) THEN
             NULLIFY (orbitals)
             CALL section_vals_val_get(enforce_occupation_section,&
                                       keyword_name="ORBITALS",&
                                       i_vals=orbitals)
             norbitals = SIZE(orbitals)
             CALL cp_assert(((norbitals > 0).AND.(norbitals <= 2*l+1)),cp_failure_level,cp_assertion_failed,routineP,&
                            "DFT+U| Invalid number of ORBITALS specified: 1 to 2*L+1 integer numbers are expected",&
                            failure)
             ALLOCATE (qs_kind%dft_plus_u%orbitals(norbitals),STAT=stat)
             CPASSERT(stat==0)
             qs_kind%dft_plus_u%orbitals(:) = orbitals(:)
             NULLIFY (orbitals)
             DO m=1,norbitals
                CALL cp_assert((qs_kind%dft_plus_u%orbitals(m) <= l),cp_failure_level,&
                               cp_assertion_failed,routineP,&
                               "DFT+U| Invalid orbital magnetic quantum number specified: m > l",&
                               failure)
                CALL cp_assert((qs_kind%dft_plus_u%orbitals(m) >= -l),cp_failure_level,&
                               cp_assertion_failed,routineP,&
                               "DFT+U| Invalid orbital magnetic quantum number specified: m < -l",&
                               failure)
                DO j=1,norbitals
                   IF (j /= m) THEN
                      CALL cp_assert((qs_kind%dft_plus_u%orbitals(j) /= qs_kind%dft_plus_u%orbitals(m)),&
                                     cp_failure_level,cp_assertion_failed,routineP,&
                                     "DFT+U| An orbital magnetic quantum number was specified twice",&
                                     failure)
                   END IF
                END DO
             END DO
             CALL section_vals_val_get(enforce_occupation_section,&
                                       keyword_name="EPS_SCF",&
                                       r_val=qs_kind%dft_plus_u%eps_scf)
             CALL section_vals_val_get(enforce_occupation_section,&
                                       keyword_name="MAX_SCF",&
                                       i_val=i)
             qs_kind%dft_plus_u%max_scf = MAX(-1,i)
             CALL section_vals_val_get(enforce_occupation_section,&
                                       keyword_name="SMEAR",&
                                       l_val=qs_kind%dft_plus_u%smear)
          END IF ! subsection enabled
       END IF ! section enabled

    END IF

    ! Allocate and initialise the orbital basis set data set structure
    CALL init_orbital_pointers(5) ! debug the SUN optimizer

    ! BASIS  and POTENTIAL read only when strictly necessary otherwise, even if not used
    ! we just print misleading informations
    explicit_basis = .FALSE.
    IF (k_rep>0) THEN
       basis_section  => section_vals_get_subs_vals(kind_section,"BASIS",i_rep_section=k_rep,&
            can_return_null=.TRUE.)
       CALL section_vals_get(basis_section,explicit=explicit_basis)
    END IF

    explicit_potential = .FALSE.
    IF (k_rep>0) THEN
       potential_section  => section_vals_get_subs_vals(kind_section,"POTENTIAL",&
            i_rep_section=k_rep,can_return_null=.TRUE.)
       CALL section_vals_get(potential_section, explicit=explicit_potential)
    END IF

    explicit_kgpot = .FALSE.
    IF (k_rep>0) THEN
       kgpot_section  => section_vals_get_subs_vals(kind_section,"KG_POTENTIAL",&
            i_rep_section=k_rep,can_return_null=.TRUE.)
       CALL section_vals_get(kgpot_section, explicit=explicit_kgpot)
    END IF

    SELECT CASE(method_id)
    CASE(do_method_dftb, do_method_scptb, do_method_rm1, do_method_am1, do_method_mndo, &
         do_method_pdg, do_method_pm3, do_method_pm6, do_method_mndod,&
         do_method_pnnl )
       ! Allocate all_potential
       CALL allocate_potential(qs_kind%all_potential)
       CALL set_default_all_potential(qs_kind%all_potential,z,zeff_correction)
       CALL get_qs_kind(qs_kind,elec_conf=elec_conf)
       IF (.NOT.ASSOCIATED(elec_conf)) THEN
          CALL get_potential(potential=qs_kind%all_potential,elec_conf=elec_conf)
          CALL set_qs_kind(qs_kind,elec_conf=elec_conf)
       END IF
       IF (qs_kind%ghost) THEN
          CALL get_qs_kind(qs_kind=qs_kind,elec_conf=elec_conf)
          elec_conf(:) = 0
          CALL get_potential(potential=qs_kind%all_potential,&
               elec_conf=elec_conf)
          elec_conf(:) = 0
          CALL set_potential(potential=qs_kind%all_potential,&
               zeff=0.0_dp,&
               zeff_correction=0.0_dp)
       END IF

       ! Basis set (Parameters)
       SELECT CASE(method_id)
       CASE (do_method_rm1, do_method_am1, do_method_mndo, do_method_pdg,&
             do_method_pm3, do_method_pm6, do_method_mndod, do_method_pnnl )
          ! Setup proper semiempirical parameters
          check = .NOT.ASSOCIATED(qs_kind%se_parameter)
          CPASSERT(check)
          CALL semi_empirical_create(qs_kind%se_parameter)
          ! Check if we allow p-orbitals on H
          SELECT CASE(z)
          CASE (1)
             IF (k_rep>0) THEN
                CALL section_vals_val_get(kind_section,i_rep_section=k_rep,&
                     keyword_name="SE_P_ORBITALS_ON_H",l_val=qs_kind%se_parameter%p_orbitals_on_h)
             END IF
          CASE DEFAULT
             ! No special cases for other elements..
          END SELECT
          ! Set default parameters
          CALL section_vals_val_get(dft_section,"QS%SE%STO_NG",i_val=ngauss)
          CALL se_param_set_default(qs_kind%se_parameter,z,method_id)
          NULLIFY(tmp_basis_set)
          CALL init_se_param(qs_kind%se_parameter,tmp_basis_set,ngauss)
          CALL add_basis_set_to_container(qs_kind%basis_sets,tmp_basis_set,"ORB")
          CALL init_potential ( qs_kind%all_potential, itype="BARE", &
               zeff=qs_kind%se_parameter%zeff,zeff_correction=zeff_correction)
          qs_kind%se_parameter%zeff=qs_kind%se_parameter%zeff-zeff_correction
       CASE (do_method_dftb,do_method_scptb)
          ! Do nothing at this stage.. The initialization of the parameters will be done
          ! later
       CASE DEFAULT
          CPASSERT(.FALSE.)
       END SELECT

       check = ((potential_name /= '').OR.explicit_potential)
       CALL cp_assert(.NOT.check,cp_warning_level,cp_assertion_failed,routineP,&
            "Information provided in the input file regarding POTENTIAL for KIND <"//&
            TRIM(qs_kind%name)//"> will be ignored!"//&
            CPSourceFileRef,&
            only_ionode=.TRUE.)

       check = ((orb_basis_set_name /= '').OR.explicit_basis)
       CALL cp_assert(.NOT.check,cp_warning_level,cp_assertion_failed,routineP,&
            "Information provided in the input file regarding BASIS for KIND <"//&
            TRIM(qs_kind%name)//"> will be ignored!"//&
            CPSourceFileRef,&
            only_ionode=.TRUE.)

    CASE DEFAULT

       ! Allocate and initialise the basis set data set structure
       ! first external basis sets
       DO i=1,nb_rep
          NULLIFY(tmp_basis_set)
          CALL allocate_gto_basis_set(tmp_basis_set)
          CALL read_gto_basis_set(qs_kind%element_symbol,basis_set_name(i),&
               tmp_basis_set,para_env,dft_section)
          tmp=basis_set_type(i)
          CALL uppercase(tmp)
          CALL add_basis_set_to_container(qs_kind%basis_sets,tmp_basis_set,tmp)
       END DO
       ! now explicit basis sets
       IF(explicit_basis) THEN
          CALL section_vals_get(basis_section,n_repetition=nexp)
          DO i=1,nexp
             NULLIFY(tmp_basis_set)
             CALL allocate_gto_basis_set(tmp_basis_set)
             CALL read_gto_basis_set(qs_kind%element_symbol,basis_type,&
                  tmp_basis_set,basis_section,i)
             tmp=basis_type
             CALL uppercase(tmp)
             CALL add_basis_set_to_container(qs_kind%basis_sets,tmp_basis_set,tmp)
          END DO
       END IF
       ! combine multiple basis sets
       DO i=1,SIZE(qs_kind%basis_sets)
          NULLIFY(tmp_basis_set)
          CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis_set,&
               inumbas=i,basis_type=basis_type)
          IF(basis_type=="") CYCLE
          DO j=i+1,SIZE(qs_kind%basis_sets)
             NULLIFY(sup_basis_set)
             CALL get_basis_from_container(qs_kind%basis_sets,basis_set=sup_basis_set,&
                  inumbas=j,basis_type=tmp)
             IF(basis_type==tmp) THEN
                ! we found a match, combine the basis sets and delete the second
                CALL combine_basis_sets(tmp_basis_set,sup_basis_set)
                CALL remove_basis_from_container(qs_kind%basis_sets,j)
             END IF
          END DO
          NULLIFY(sup_basis_set)
       END DO
       NULLIFY(tmp_basis_set)

       ! check that we have an orbital basis set
       nobasis = .TRUE.
       DO i=1,SIZE(qs_kind%basis_sets)
          NULLIFY(tmp_basis_set)
          CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis_set,&
               inumbas=i,basis_type=basis_type)
          IF(basis_type == "ORB") nobasis = .FALSE.
       END DO
       IF (nobasis) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "No basis set type was defined for the "//&
               "atomic kind <"//TRIM(qs_kind%name)//">",para_env)
       END IF

       ! If Ghost atom we don't need to allocate/initialize anything connected to POTENTIAL
       IF (qs_kind%ghost) THEN
          IF (ASSOCIATED(qs_kind%elec_conf)) qs_kind%elec_conf = 0
       ELSE
          ! Allocate and initialise the potential data set structure
          IF ((potential_name /= '').OR.explicit_potential) THEN
             ipos = INDEX(potential_name,"-")
             IF (ipos > 1) THEN
                potential_type = potential_name(:ipos-1)
             ELSE
                potential_type = potential_name
             END IF
             CALL uppercase(potential_type)

             SELECT CASE (TRIM(potential_type))
             CASE ("ALL")
                CALL allocate_potential(qs_kind%all_potential)
                CALL read_potential(qs_kind%element_symbol,potential_name,&
                     qs_kind%all_potential,zeff_correction,para_env,&
                     dft_section,potential_section)
                CALL set_potential(qs_kind%all_potential,z=z)
                CALL get_qs_kind(qs_kind,elec_conf=elec_conf)
                IF (.NOT.ASSOCIATED(elec_conf)) THEN
                   CALL get_potential(potential=qs_kind%all_potential,elec_conf=elec_conf)
                   CALL set_qs_kind(qs_kind,elec_conf=elec_conf)
                END IF
             CASE ("GTH")
                CALL allocate_potential(qs_kind%gth_potential)
                CALL read_potential(qs_kind%element_symbol,potential_name,&
                     qs_kind%gth_potential,zeff_correction,para_env,&
                     dft_section,potential_section)
                CALL set_potential(qs_kind%gth_potential,z=z)
                CALL get_qs_kind(qs_kind,elec_conf=elec_conf)
                IF (.NOT.ASSOCIATED(elec_conf)) THEN
                   CALL get_potential(potential=qs_kind%gth_potential,elec_conf=elec_conf)
                   CALL set_qs_kind(qs_kind,elec_conf=elec_conf)
                END IF
             CASE DEFAULT
                CALL stop_program(routineN,moduleN,__LINE__,&
                     "An invalid potential type <"//&
                     TRIM(potential_name)//"> was specified "//&
                     "for the atomic kind <"//&
                     TRIM(qs_kind%name),para_env)
             END SELECT
          ELSE
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "No potential type was defined for the "//&
                  "atomic kind <"//TRIM(qs_kind%name)//">",para_env)
          END IF
          ! Allocate and initialise the potential data set structure
          IF ((kgpot_name /= '').OR.explicit_kgpot) THEN
             ipos = INDEX(kgpot_name,"-")
             IF (ipos > 1) THEN
                kgpot_type = kgpot_name(:ipos-1)
             ELSE
                kgpot_type = kgpot_name
             END IF
             CALL uppercase(kgpot_type)

             SELECT CASE (TRIM(kgpot_type))
             CASE ("TNADD")
                CALL allocate_potential(qs_kind%tnadd_potential)
                CALL read_potential(qs_kind%element_symbol,kgpot_name,&
                     qs_kind%tnadd_potential,para_env,&
                     dft_section,kgpot_section)
             CASE ("NONE")
                NULLIFY(qs_kind%tnadd_potential)
             CASE DEFAULT
                CALL stop_program(routineN,moduleN,__LINE__,&
                     "An invalid kg_potential type <"//&
                     TRIM(potential_name)//"> was specified "//&
                     "for the atomic kind <"//&
                     TRIM(qs_kind%name),para_env)
             END SELECT
          END IF
       END IF
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE read_qs_kind

! *****************************************************************************
!> \brief Read an atomic kind set data set from the input file.
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param kind_section ...
!> \param para_env ...
!> \param force_env_section ...
! *****************************************************************************
  SUBROUTINE create_qs_kind_set(qs_kind_set, atomic_kind_set,kind_section,para_env,force_env_section)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: kind_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: force_env_section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qs_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, method, nkind, &
                                                qs_method
    LOGICAL                                  :: no_fail

    CALL timeset(routineN,handle)

    IF(ASSOCIATED(qs_kind_set)) STOP "create_qs_kind_set: qs_kind_set already associated"
    IF(.NOT.ASSOCIATED(atomic_kind_set)) STOP "create_qs_kind_set: atomic_kind_set not associated"

    no_fail = .FALSE.

    ! Between all methods only SE and DFTB may not need a KIND section.
    CALL section_vals_val_get(force_env_section,"METHOD",i_val=method)
    IF (method==do_qs) THEN
       CALL section_vals_val_get(force_env_section,"DFT%QS%METHOD",i_val=qs_method)
       SELECT CASE (qs_method)
       CASE (do_method_mndo, do_method_am1, do_method_pm3, do_method_pm6,&
             do_method_pdg,  do_method_rm1, do_method_mndod, do_method_pnnl )
          no_fail = .TRUE.
       CASE (do_method_dftb, do_method_scptb)
          no_fail = .TRUE.
       END SELECT
    ELSE
       qs_method = method
    END IF

    nkind = SIZE(atomic_kind_set)
    ALLOCATE(qs_kind_set(nkind))

    DO ikind=1,nkind
       qs_kind_set(ikind)%name = atomic_kind_set(ikind)%name
       qs_kind_set(ikind)%element_symbol = atomic_kind_set(ikind)%element_symbol
       qs_kind_set(ikind)%natom = atomic_kind_set(ikind)%natom
       CALL read_qs_kind(qs_kind_set(ikind), kind_section, para_env, force_env_section, no_fail, qs_method)
    END DO

    CALL timestop(handle)

  END SUBROUTINE create_qs_kind_set

! *****************************************************************************
!> \brief This routines should perform only checks. no settings are allowed at
!>     this level anymore..
!> \param qs_kind ...
!> \param dft_control ...
!> \param subsys_section ...
! *****************************************************************************
  SUBROUTINE check_qs_kind(qs_kind,dft_control,subsys_section)

    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_qs_kind', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: defined
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter
    TYPE(scptb_parameter_type), POINTER      :: scptb_parameter
    TYPE(semi_empirical_type), POINTER       :: se_parameter

    IF (dft_control%qs_control%semi_empirical) THEN
       CALL get_qs_kind(qs_kind,se_parameter=se_parameter)
       CPASSERT(ASSOCIATED(se_parameter))
       CALL get_se_param(se_parameter,defined=defined)
       CPASSERT(defined)
       CALL write_se_param(se_parameter,subsys_section)
    ELSE IF (dft_control%qs_control%dftb) THEN
       CALL get_qs_kind(qs_kind,dftb_parameter=dftb_parameter)
       CPASSERT(ASSOCIATED(dftb_parameter))
       CALL get_dftb_atom_param(dftb_parameter,defined=defined)
       CPASSERT(defined)
       CALL write_dftb_atom_param(dftb_parameter,subsys_section)
    ELSE IF (dft_control%qs_control%scptb) THEN
       CALL get_qs_kind(qs_kind,scptb_parameter=scptb_parameter)
       CPASSERT(ASSOCIATED(scptb_parameter))
       CALL get_scptb_parameter(scptb_parameter,defined=defined)
       CPASSERT(defined)
       CALL write_scptb_parameter(scptb_parameter,subsys_section)
    END IF

  END SUBROUTINE check_qs_kind

! *****************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param dft_control ...
!> \param para_env ...
!> \param subsys_section ...
! *****************************************************************************
  SUBROUTINE check_qs_kind_set(qs_kind_set,dft_control,para_env,subsys_section)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(dft_control_type), INTENT(IN)       :: dft_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(len=*), PARAMETER :: routineN = 'check_qs_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, nkind
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)
    IF (ASSOCIATED(qs_kind_set)) THEN
       nkind = SIZE(qs_kind_set)
       DO ikind=1,nkind
          qs_kind => qs_kind_set(ikind)
          CALL check_qs_kind(qs_kind,dft_control,subsys_section)
       END DO
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer qs_kind_set is not associated",&
            para_env)
    END IF
    CALL timestop(handle)
  END SUBROUTINE check_qs_kind_set

! *****************************************************************************
!> \brief Set the components of an atomic kind data set.
!> \param qs_kind ...
!> \param paw_atom ...
!> \param ghost ...
!> \param hard_radius ...
!> \param hard0_radius ...
!> \param soft_basis_set ...
!> \param hard_basis_set ...
!> \param lmax_rho0 ...
!> \param zeff ...
!> \param no_optimize ...
!> \param dispersion ...
!> \param u_minus_j ...
!> \param reltmat ...
!> \param dftb_parameter ...
!> \param scptb_parameter ...
!> \param elec_conf ...
! *****************************************************************************
  SUBROUTINE set_qs_kind(qs_kind,paw_atom,ghost,hard_radius,hard0_radius,&
       soft_basis_set,hard_basis_set,lmax_rho0,zeff,&
       no_optimize,dispersion,u_minus_j,reltmat,&
       dftb_parameter, scptb_parameter,&
       elec_conf)

    TYPE(qs_kind_type), INTENT(INOUT)        :: qs_kind
    LOGICAL, INTENT(IN), OPTIONAL            :: paw_atom, ghost
    REAL(dp), INTENT(IN), OPTIONAL           :: hard_radius, hard0_radius
    TYPE(gto_basis_set_type), OPTIONAL, &
      POINTER                                :: soft_basis_set, hard_basis_set
    INTEGER, INTENT(IN), OPTIONAL            :: lmax_rho0
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: zeff
    LOGICAL, INTENT(IN), OPTIONAL            :: no_optimize
    TYPE(qs_atom_dispersion_type), &
      OPTIONAL, POINTER                      :: dispersion
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: u_minus_j
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: reltmat
    TYPE(qs_dftb_atom_type), OPTIONAL, &
      POINTER                                :: dftb_parameter
    TYPE(scptb_parameter_type), OPTIONAL, &
      POINTER                                :: scptb_parameter
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: elec_conf

    CHARACTER(len=*), PARAMETER :: routineN = 'set_qs_kind', &
      routineP = moduleN//':'//routineN

    IF (PRESENT(dftb_parameter))  qs_kind%dftb_parameter  => dftb_parameter
    IF (PRESENT(scptb_parameter)) qs_kind%scptb_parameter => scptb_parameter
    IF (PRESENT(elec_conf)) THEN
       IF (ASSOCIATED(qs_kind%elec_conf)) THEN
          DEALLOCATE(qs_kind%elec_conf)
       ENDIF
       ALLOCATE(qs_kind%elec_conf(0:SIZE(elec_conf)-1))
       qs_kind%elec_conf(:)=elec_conf(:)
    ENDIF
    IF (PRESENT(paw_atom)) qs_kind%paw_atom = paw_atom
    IF (PRESENT(hard_radius)) qs_kind%hard_radius = hard_radius
    IF (PRESENT(hard0_radius)) qs_kind%hard0_radius = hard0_radius
    IF (PRESENT(soft_basis_set)) qs_kind%soft_basis_set => soft_basis_set
    IF (PRESENT(hard_basis_set)) qs_kind%hard_basis_set => hard_basis_set
    IF (PRESENT(lmax_rho0)) qs_kind%lmax_rho0 = lmax_rho0
    IF (PRESENT(zeff)) THEN
       IF (ASSOCIATED(qs_kind%all_potential)) THEN
          CALL set_potential(potential=qs_kind%all_potential,zeff=zeff)
       ELSE IF (ASSOCIATED(qs_kind%gth_potential)) THEN
          CALL set_potential(potential=qs_kind%gth_potential,zeff=zeff)
       END IF
    END IF
    IF (PRESENT(ghost)) qs_kind%ghost = ghost

    IF (PRESENT(no_optimize)) qs_kind%no_optimize=no_optimize

    IF (PRESENT(dispersion)) qs_kind%dispersion => dispersion

    IF (PRESENT(u_minus_j)) THEN
       IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
          qs_kind%dft_plus_u%u_minus_j = u_minus_j
       END IF
    END IF

    IF (PRESENT(reltmat)) qs_kind%reltmat => reltmat

  END SUBROUTINE set_qs_kind

! *****************************************************************************
!> \brief Write an atomic kind data set to the output unit.
!> \param qs_kind ...
!> \param kind_number ...
!> \param output_unit ...
!> \par History
!>      Creation (09.02.2002,MK)
! *****************************************************************************
  SUBROUTINE write_qs_kind(qs_kind,kind_number,output_unit)

    TYPE(qs_kind_type), POINTER              :: qs_kind
    INTEGER, INTENT(in)                      :: kind_number, output_unit

    CHARACTER(len=*), PARAMETER :: routineN = 'write_qs_kind', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: yon
    CHARACTER(LEN=default_string_length)     :: basis_type, bstring
    INTEGER                                  :: ibas
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: tmp_basis

    failure = .FALSE.

    IF (output_unit > 0) THEN

       IF (ASSOCIATED(qs_kind)) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,I2,A,T57,A,T75,I6)")&
               kind_number,". Atomic kind: "//TRIM(qs_kind%name),&
               "Number of atoms: ",qs_kind%natom

          DO ibas=1,SIZE(qs_kind%basis_sets,1)
             NULLIFY(tmp_basis)
             CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis,&
                  inumbas=ibas,basis_type=basis_type)
             IF(basis_type=="") CYCLE
             SELECT CASE (basis_type)
                CASE DEFAULT
                   bstring = "Basis Set"
                CASE ("ORB")
                   bstring = "Orbital Basis Set"
                CASE ("AUX")
                   bstring = "Auxiliary Basis Set"
                CASE ("RI_AUX")
                   bstring = "RI Auxiliary Basis Set"
                CASE ("AUX_FIT")
                   bstring = "Auxiliary Fit Basis Set"
                CASE ("LRI")
                   bstring = "LRI Basis Set" 
             END SELECT

             CALL write_orb_basis_set(tmp_basis,output_unit,bstring)

          END DO

          IF (qs_kind%ghost) THEN
             WRITE (UNIT=output_unit,FMT="(/,T6,A)")&
                  "The atoms of this atomic kind are GHOST atoms!"
          END IF
          IF (qs_kind%paw_atom) THEN
             WRITE (UNIT=output_unit,FMT="(/,T6,A)")&
                  "The atoms of this atomic kind are PAW atoms (GAPW):"
             WRITE (UNIT=output_unit,FMT="(T8,A,T71,F10.3)")&
                  "Hard Gaussian function radius:",qs_kind%hard_radius,&
                  "Rho0 radius:",qs_kind%hard0_radius,&
                  "Maximum GTO radius used for PAW projector construction:",&
                  qs_kind%max_rad_local
             CALL write_orb_basis_set(qs_kind%soft_basis_set,output_unit,&
                  "GAPW Soft Basis Set")
          END IF
          ! Potentials
          CALL write_potential(qs_kind%all_potential,output_unit)
          CALL write_potential(qs_kind%gth_potential,output_unit)
          CALL write_potential(qs_kind%tnadd_potential,output_unit)
          IF (ASSOCIATED(qs_kind%dft_plus_u)) THEN
             WRITE (UNIT=output_unit,FMT="(/,T6,A,/,T8,A,T76,I5,/,T8,A,T73,F8.3)")&
               "A DFT+U correction is applied to atoms of this atomic kind:",&
               "Angular quantum momentum number L:",qs_kind%dft_plus_u%l,&
               "U(eff) = (U - J) value in [eV]:",qs_kind%dft_plus_u%u_minus_j_target*evolt
             IF (qs_kind%dft_plus_u%u_ramping > 0.0_dp) THEN
                IF (qs_kind%dft_plus_u%init_u_ramping_each_scf) THEN
                   yon = "YES"
                ELSE
                   yon = " NO"
                END IF
                WRITE (UNIT=output_unit,FMT="(T8,A,T73,F8.3,/,T8,A,T73,ES8.1,/,T8,A,T78,A3)")&
                  "Increment for U ramping in [eV]:",qs_kind%dft_plus_u%u_ramping*evolt,&
                  "SCF threshold value for U ramping:",qs_kind%dft_plus_u%eps_u_ramping,&
                  "Set U ramping value to zero before each wavefunction optimisation:",yon
             END IF
             IF (ASSOCIATED(qs_kind%dft_plus_u%orbitals)) THEN
                WRITE (UNIT=output_unit,FMT="(T8,A)")&
                  "An initial orbital occupation is requested:"
                WRITE (UNIT=output_unit,FMT="(T9,A,(T78,I3))")&
                  "Preferred (initial) orbital occupation order (orbital M values):",&
                  qs_kind%dft_plus_u%orbitals(:)
                WRITE (UNIT=output_unit,FMT="(T9,A,T71,ES10.3,/,T9,A,T76,I5)")&
                  "Threshold value for the SCF convergence criterion:",&
                  qs_kind%dft_plus_u%eps_scf,&
                  "Number of initial SCF iterations:",&
                  qs_kind%dft_plus_u%max_scf
                IF (qs_kind%dft_plus_u%smear) THEN
                   WRITE (UNIT=output_unit,FMT="(T9,A)")&
                     "A smearing of the orbital occupations will be performed"
                END IF
             END IF
          END IF
       ELSE
          CPASSERT(.FALSE.)
       END IF

    END IF

  END SUBROUTINE write_qs_kind

! *****************************************************************************
!> \brief Write an atomic kind set data set to the output unit.
!> \param qs_kind_set ...
!> \param subsys_section ...
!> \par History
!>      Creation (09.02.2002,MK)
! *****************************************************************************
  SUBROUTINE write_qs_kind_set(qs_kind_set,subsys_section)
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(len=*), PARAMETER :: routineN = 'write_qs_kind_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, nkind, &
                                                output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (logger)
    logger => cp_get_default_logger()
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%KINDS",extension=".Log")
    IF (output_unit > 0) THEN
       IF (ASSOCIATED(qs_kind_set)) THEN
          WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "ATOMIC KIND INFORMATION"
          nkind = SIZE(qs_kind_set)
          DO ikind=1,nkind
             qs_kind => qs_kind_set(ikind)
             CALL write_qs_kind(qs_kind,ikind,output_unit)
          END DO
       ELSE
          CPASSERT(.FALSE.)
       END IF
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%KINDS")

    CALL timestop(handle)

  END SUBROUTINE write_qs_kind_set

! *****************************************************************************
!> \brief Write all the GTO basis sets of an atomic kind set to the output
!>     unit (for the printing of the unnormalized basis sets as read from
!>           database).
!> \param qs_kind_set ...
!> \param subsys_section ...
!> \par History
!>      Creation (17.01.2002,MK)
! *****************************************************************************
  SUBROUTINE write_gto_basis_sets(qs_kind_set,subsys_section)

    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: subsys_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_gto_basis_sets', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: basis_type, bstring
    INTEGER                                  :: handle, ibas, ikind, nkind, &
                                                output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: tmp_basis
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (logger)
    logger => cp_get_default_logger()
    output_unit = cp_print_key_unit_nr(logger,subsys_section,&
         "PRINT%KINDS/BASIS_SET",&
         extension=".Log")
    IF (output_unit > 0) THEN
       IF (ASSOCIATED(qs_kind_set)) THEN
          WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
               "BASIS SET INFORMATION (Unnormalised Gaussian-type functions)"
          nkind = SIZE(qs_kind_set)
          DO ikind=1,nkind
             qs_kind => qs_kind_set(ikind)
             WRITE (UNIT=output_unit,FMT="(/,T2,I2,A)")&
                  ikind,". Atomic kind: "//TRIM(qs_kind%name)

             DO ibas=1,SIZE(qs_kind%basis_sets,1)
                NULLIFY(tmp_basis)
                CALL get_basis_from_container(qs_kind%basis_sets,basis_set=tmp_basis,&
                     inumbas=ibas,basis_type=basis_type)
                IF(basis_type=="") CYCLE
                SELECT CASE (basis_type)
                   CASE DEFAULT
                      bstring = "Basis Set"
                   CASE ("ORB")
                      bstring = "Orbital Basis Set"
                   CASE ("AUX")
                      bstring = "Auxiliary Basis Set"
                   CASE ("RI_AUX")
                      bstring = "RI Auxiliary Basis Set"
                   CASE ("AUX_FIT")
                      bstring = "Auxiliary Fit Basis Set"
                   CASE ("LRI")
                      bstring = "LRI Basis Set"
                END SELECT

                CALL write_gto_basis_set(tmp_basis,output_unit,bstring)
   
             END DO

             CALL write_gto_basis_set(qs_kind%soft_basis_set,&
                  output_unit,"GAPW Soft Basis Set")

          END DO
       ELSE
          CPASSERT(.FALSE.)
       END IF
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,subsys_section,&
         "PRINT%KINDS/BASIS_SET")

    CALL timestop(handle)

  END SUBROUTINE write_gto_basis_sets

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param atomic_kind ...
!> \param qs_kind ...
!> \param ispin ...
!> \param ncalc ...
!> \param ncore ...
!> \param nelem ...
! *****************************************************************************
  SUBROUTINE init_atom_electronic_state (atomic_kind,qs_kind,ispin,ncalc,ncore,nelem)
    
    TYPE(atomic_kind_type), INTENT(IN)       :: atomic_kind
    TYPE(qs_kind_type), INTENT(IN)           :: qs_kind
    INTEGER, INTENT(IN)                      :: ispin
    INTEGER, DIMENSION(0:3, 10), INTENT(OUT) :: ncalc, ncore, nelem

    CHARACTER(len=*), PARAMETER :: routineN = 'init_atom_electronic_state', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, l, ll, ne, nn, z
    INTEGER, DIMENSION(0:3)                  :: econfx
    INTEGER, DIMENSION(:), POINTER           :: econf
    INTEGER, DIMENSION(:, :), POINTER        :: addel, laddel, naddel
    LOGICAL                                  :: bs_occupation
    TYPE(gth_potential_type), POINTER        :: gth_potential

    CALL get_atomic_kind(atomic_kind,z=z)
    NULLIFY(gth_potential)
    CALL get_qs_kind(qs_kind,&
                     gth_potential=gth_potential,&
                     bs_occupation=bs_occupation,&
                     addel=addel,laddel=laddel,naddel=naddel)

    ! electronic state
    nelem = 0
    ncore = 0
    ncalc = 0
    NULLIFY(econf)
    IF ( ASSOCIATED(gth_potential) ) THEN
      CALL get_potential(gth_potential,elec_conf=econf)
      econfx = 0
      econfx(0:SIZE(econf)-1) = econf
      IF ( SUM(econf) >= 0 ) THEN
        DO l=0,3
          ll = 2*(2*l+1)
          nn = ptable(z)%e_conv(l)-econfx(l)
          ii = 0
          DO
            ii = ii + 1
            IF(nn <= ll) THEN
              ncore(l,ii) = nn
              EXIT
            ELSE
              ncore(l,ii) = ll
              nn = nn - ll
            END IF
          END DO
        END DO
        DO l=0,3
          ll = 2*(2*l+1)
          nn = ptable(z)%e_conv(l)
          ii = 0
          DO
            ii = ii + 1
            IF(nn <= ll) THEN
              nelem(l,ii) = nn
              EXIT
            ELSE
              nelem(l,ii) = ll
              nn = nn - ll
            END IF
          END DO
        END DO
        ncalc = nelem - ncore
      ELSE
        ncore = 0
        ncalc = 0
        DO l=0,3
          ll = 2*(2*l+1)
          nn = ABS(econfx(l))
          ii = 0
          DO
            ii = ii + 1
            IF(nn <= ll) THEN
              ncalc(l,ii) = -nn
              EXIT
            ELSE
              ncalc(l,ii) = -ll
              nn = nn - ll
            END IF
          END DO
        END DO
      END IF
    ELSE
      DO l=0,3
        ll = 2*(2*l+1)
        nn = ptable(z)%e_conv(l)
        ii = 0
        DO
          ii = ii + 1
          IF(nn <= ll) THEN
            nelem(l,ii) = nn
            EXIT
          ELSE
            nelem(l,ii) = ll
            nn = nn - ll
          END IF
        END DO
      END DO
      ncalc = nelem - ncore
    END IF

    ! readjust the occupation number of the orbitals as requested by user
    ! this is done to break symmetry (bs) and bias the intial guess 
    ! to the pre-defined multiplicity/charge state of the atom
    IF(bs_occupation) THEN
      DO i=1,SIZE(addel,1)
        ne=addel(i,ispin)
        l=laddel(i,ispin)
        nn=naddel(i,ispin)-l
        IF(ne/=0) THEN
          IF(nn==0) THEN
            DO ii=SIZE(nelem,2),1,-1
              IF(ncalc(l,ii)>0) THEN
                 IF((ncalc(l,ii)+ne) < 2*(2*l+1)+1) THEN
                   ncalc(l,ii) = ncalc(l,ii)+ne
                   nn = ii
                 ELSE
                   ncalc(l,ii+1) = ncalc(l,ii+1)+ne
                   nn = ii + 1
                 END IF
                 EXIT
              ELSE IF (ii==1) THEN
                 ncalc(l,ii)=ncalc(l,ii)+ne
                 nn = ii
              END IF
            END DO
          ELSE
            ncalc(l,nn) = ncalc(l,nn) + ne
          END IF
          IF( ncalc(l,nn)<0 ) THEN
             ncalc(l,nn) = 0
          END IF
        END IF
      END DO
    END IF

    IF ( qs_kind%ghost ) THEN
      nelem = 0
      ncore = 0
      ncalc = 0
    END IF

  END SUBROUTINE init_atom_electronic_state 

! *****************************************************************************

END MODULE qs_kind_types
