!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 -2004 CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/atomic_kind_types [1.0] *
!!
!!   NAME
!!     atomic_kind_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MK (02.01.2002)
!!
!!   MODIFICATION HISTORY
!!     - 01.2002 creation [MK]
!!     - 04.2002 added pao [fawzi]
!!     - 09.2002 adapted for POL/KG use [GT]
!!     - 02.2004 flexible normalization of basis sets [jgh]
!!     - 03.2004 attach/detach routines [jgh]
!!     - 10.2004 removed pao [fawzi]
!!   SOURCE
!******************************************************************************

MODULE atomic_kind_types

  USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                             deallocate_gto_basis_set,&
                                             get_gto_basis_set,&
                                             gto_basis_set_type,&
                                             init_aux_basis_set,&
                                             init_orb_basis_set,&
                                             read_gto_basis_set,&
                                             set_gto_basis_set,&
                                             write_aux_basis_set,&
                                             write_gto_basis_set,&
                                             write_orb_basis_set
  USE cp_control_types,                ONLY: dft_control_type, &
                                             gapw_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_linked_list_int,              ONLY: cp_sll_int_create,&
                                             cp_sll_int_dealloc,&
                                             cp_sll_int_get_empty,&
                                             cp_sll_int_insert_ordered,&
                                             cp_sll_int_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE external_potential_types,        ONLY: all_potential_type,&
                                             allocate_potential,&
                                             deallocate_potential,&
                                             elp_potential_type,&
                                             get_potential,&
                                             gth_potential_type,&
                                             init_potential,&
                                             read_potential,&
                                             set_default_all_potential,&
                                             set_potential,&
                                             write_potential
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE paw_proj_set_types,              ONLY: allocate_paw_proj_set,&
                                             deallocate_paw_proj_set,&
                                             get_paw_proj_set,&
                                             paw_proj_set_type,&
                                             projectors
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: bohr,massunit
  USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                             deallocate_grid_atom,&
                                             grid_atom_type
  USE qs_harmonics_atom,               ONLY: allocate_harmonics_atom,&
                                             deallocate_harmonics_atom,&
                                             harmonics_atom_type
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: allocate_se_param,&
                                             deallocate_se_param,&
                                             init_se_param,&
                                             write_se_param,&
                                             get_se_param,&
                                             se_param_set_default
  USE soft_basis_set,                  ONLY: create_soft_basis
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "atomic_kind_types"

! *** Define the atomic kind types ***

  TYPE atomic_kind_type
    PRIVATE
    TYPE(all_potential_type), POINTER      :: all_potential
    TYPE(elp_potential_type), POINTER      :: elp_potential
    TYPE(gth_potential_type), POINTER      :: gth_potential
    TYPE(semi_empirical_type), POINTER     :: se_parameter
    TYPE(gto_basis_set_type), POINTER      :: aux_basis_set
    TYPE(gto_basis_set_type), POINTER      :: orb_basis_set
    TYPE(gto_basis_set_type), POINTER      :: soft_basis_set
    TYPE(gto_basis_set_type), POINTER      :: hard_basis_set
    TYPE(paw_proj_set_type),  POINTER      :: paw_proj_set
    REAL(dp)                               :: hard_radius  ! for hard and soft exp
    REAL(dp)                               :: hard0_radius ! for hard exp of rho0
    LOGICAL                                :: paw_atom     ! needs atomic rho1
    LOGICAL                                :: qm_atom    ! needs atomic rho1
    CHARACTER(LEN=default_string_length)   :: name
    CHARACTER(LEN=2)                       :: element_symbol
    REAL(KIND = dp)                        :: mass
    INTEGER                                :: kind_number,&
                                              natom,&
                                              number_of_grid_points
    INTEGER, DIMENSION(:), POINTER         :: atom_list
    REAL(KIND=dp), DIMENSION(4)            :: ldaCoreRadius
    TYPE(harmonics_atom_type), POINTER     :: harmonics
    TYPE(grid_atom_type), POINTER          :: grid_atom
    INTEGER                                :: ngrid_rad,ngrid_ang 
    INTEGER                                :: lmax_rho0
  END TYPE atomic_kind_type

  TYPE atomic_kind_p_type
    TYPE(atomic_kind_type),DIMENSION(:),&
       POINTER                             :: atomic_kind_set
  END TYPE atomic_kind_p_type

! *** Public subroutines ***

  PUBLIC :: allocate_atomic_kind_set,&
            check_atomic_kind_set,&
            deallocate_atomic_kind_set,&
            get_atomic_kind,&
            get_atomic_kind_set,&
            init_atomic_kind_set,&
            init_gapw_basis_set,&
            read_atomic_kind_set,&
            set_atomic_kind,&
            write_atomic_kind,&
            write_atomic_kind_set,&
            write_gto_basis_sets

! *** Public data types ***

  PUBLIC :: atomic_kind_type,&
            atomic_kind_p_type

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_atomic_kind_set(atomic_kind_set,nkind)

!   Purpose: Allocate an atomic kind set.

!   History: - Creation (08.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_atomic_kind_set"

    INTEGER                                  :: ikind, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN
      CALL deallocate_atomic_kind_set(atomic_kind_set)
    END IF

    ALLOCATE (atomic_kind_set(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "atomic_kind_set",nkind*int_size)

    DO ikind=1,nkind
      NULLIFY (atomic_kind_set(ikind)%all_potential)
      NULLIFY (atomic_kind_set(ikind)%elp_potential)
      NULLIFY (atomic_kind_set(ikind)%gth_potential)
      NULLIFY (atomic_kind_set(ikind)%se_parameter)
      NULLIFY (atomic_kind_set(ikind)%aux_basis_set)
      NULLIFY (atomic_kind_set(ikind)%orb_basis_set)
      NULLIFY (atomic_kind_set(ikind)%atom_list)
      NULLIFY (atomic_kind_set(ikind)%hard_basis_set)
      NULLIFY (atomic_kind_set(ikind)%soft_basis_set)
      NULLIFY (atomic_kind_set(ikind)%paw_proj_set)
      NULLIFY (atomic_kind_set(ikind)%harmonics)
      NULLIFY (atomic_kind_set(ikind)%grid_atom)
      atomic_kind_set(ikind)%name = ""
      atomic_kind_set(ikind)%element_symbol = ""
      atomic_kind_set(ikind)%mass = 0.0_dp
      atomic_kind_set(ikind)%kind_number = 0
      atomic_kind_set(ikind)%natom = 0
      atomic_kind_set(ikind)%number_of_grid_points = 0
      atomic_kind_set(ikind)%ldaCoreRadius = 0.0_dp
      atomic_kind_set(ikind)%paw_atom = .FALSE.
      atomic_kind_set(ikind)%qm_atom = .FALSE.
      atomic_kind_set(ikind)%hard_radius = 0.8_dp*bohr
      atomic_kind_set(ikind)%hard0_radius = 0.8_dp*bohr
      atomic_kind_set(ikind)%ngrid_rad = 50
      atomic_kind_set(ikind)%ngrid_ang = 50 
      atomic_kind_set(ikind)%lmax_rho0 = 0 
    END DO

  END SUBROUTINE allocate_atomic_kind_set

! *****************************************************************************

  SUBROUTINE deallocate_atomic_kind_set(atomic_kind_set)

!   Purpose: Deallocate an atomic kind set.

!   History: - Creation (09.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_atomic_kind_set"

    INTEGER                                  :: ikind, istat, nkind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN

      nkind = SIZE(atomic_kind_set)

      DO ikind=1,nkind
        IF (ASSOCIATED(atomic_kind_set(ikind)%all_potential)) THEN
          CALL deallocate_potential(atomic_kind_set(ikind)%all_potential)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%elp_potential)) THEN
          CALL deallocate_potential(atomic_kind_set(ikind)%elp_potential)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%gth_potential)) THEN
          CALL deallocate_potential(atomic_kind_set(ikind)%gth_potential)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%se_parameter)) THEN
          CALL deallocate_se_param(atomic_kind_set(ikind)%se_parameter)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%aux_basis_set)) THEN
          CALL deallocate_gto_basis_set(atomic_kind_set(ikind)%aux_basis_set)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%soft_basis_set).AND.&
                                         atomic_kind_set(ikind)%paw_atom) THEN
          CALL deallocate_gto_basis_set(atomic_kind_set(ikind)%soft_basis_set)
        ELSEIF(ASSOCIATED(atomic_kind_set(ikind)%soft_basis_set).AND.&
                                    .NOT.  atomic_kind_set(ikind)%paw_atom) THEN
          NULLIFY (atomic_kind_set(ikind)%soft_basis_set)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%orb_basis_set)) THEN
          CALL deallocate_gto_basis_set(atomic_kind_set(ikind)%orb_basis_set)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%paw_proj_set)) THEN
          CALL deallocate_paw_proj_set(atomic_kind_set(ikind)%paw_proj_set)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%harmonics)) THEN
          CALL deallocate_harmonics_atom(atomic_kind_set(ikind)%harmonics)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%grid_atom)) THEN
          CALL deallocate_grid_atom(atomic_kind_set(ikind)%grid_atom)
        END IF

        DEALLOCATE (atomic_kind_set(ikind)%atom_list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "atomic_kind_set(ikind)%atom_list")
      END DO

      DEALLOCATE (atomic_kind_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atomic_kind_set")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind_set is not associated and "//&
                        "cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_atomic_kind_set

! *****************************************************************************

  SUBROUTINE get_atomic_kind(atomic_kind,all_potential,elp_potential,&
                             gth_potential,se_parameter,&
                             aux_basis_set,orb_basis_set,&
                             element_symbol,name,mass,kind_number,natom,&
                             atom_list,alpha_core_charge,ccore_charge,&
                             core_charge_radius,ncgf,nsgf,z,zeff,qeff,&
                             ldaCoreRadius,soft_basis_set,&
                             hard_basis_set,paw_proj_set,softb,&
                             paw_atom,qm_atom,hard_radius,hard0_radius,&
                             harmonics, max_iso_not0, max_s_harm,&
                             grid_atom, ngrid_ang, ngrid_rad, lmax_rho0)

!   Purpose: Get informations about an atomic kind.

!   History: - Creation (10.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(all_potential_type), OPTIONAL, &
      POINTER                                :: all_potential
    TYPE(elp_potential_type), OPTIONAL, &
      POINTER                                :: elp_potential
    TYPE(gth_potential_type), OPTIONAL, &
      POINTER                                :: gth_potential
    TYPE(semi_empirical_type), OPTIONAL, &
      POINTER                                :: se_parameter
    TYPE(gto_basis_set_type), OPTIONAL, &
      POINTER                                :: aux_basis_set,&
                                                orb_basis_set,&
                                                soft_basis_set,&
                                                hard_basis_set
    LOGICAL, INTENT(IN), OPTIONAL            :: softb
    TYPE(paw_proj_set_type),  OPTIONAL, &
      POINTER                                :: paw_proj_set
    CHARACTER(LEN=2), INTENT(OUT), OPTIONAL  :: element_symbol
    CHARACTER(LEN=default_string_length), &
      INTENT(OUT), OPTIONAL                  :: name
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: mass
    INTEGER, INTENT(OUT), OPTIONAL           :: kind_number, natom
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: atom_list
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: alpha_core_charge, &
                                                ccore_charge, &
                                                core_charge_radius
    INTEGER, INTENT(OUT), OPTIONAL           :: lmax_rho0, ncgf,&
                                                ngrid_ang,&
                                                ngrid_rad, nsgf,&
                                                z, zeff
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: qeff
    REAL(kind=dp), DIMENSION(4), &
      INTENT(out), OPTIONAL                  :: ldaCoreRadius
    LOGICAL, OPTIONAL,  INTENT(OUT)          :: paw_atom
    LOGICAL, OPTIONAL,  INTENT(OUT)          :: qm_atom
    REAL(dp), OPTIONAL, INTENT(OUT)          :: hard_radius
    REAL(dp), OPTIONAL, INTENT(OUT)          :: hard0_radius
    TYPE(harmonics_atom_type), OPTIONAL, &
                                 POINTER     :: harmonics
    INTEGER, INTENT(OUT), OPTIONAL           :: max_iso_not0, &
                                                max_s_harm
    TYPE(grid_atom_type), OPTIONAL, &
                            POINTER          :: grid_atom


    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_atomic_kind"

    LOGICAL                                  :: softb_local

!   ---------------------------------------------------------------------------

    softb_local = .FALSE.

    IF (ASSOCIATED(atomic_kind)) THEN

      IF (PRESENT(softb)) softb_local = softb
      IF (PRESENT(all_potential)) all_potential => atomic_kind%all_potential
      IF (PRESENT(elp_potential)) elp_potential => atomic_kind%elp_potential
      IF (PRESENT(gth_potential)) gth_potential => atomic_kind%gth_potential
      IF (PRESENT(se_parameter))  se_parameter  => atomic_kind%se_parameter
      IF (PRESENT(aux_basis_set)) aux_basis_set => atomic_kind%aux_basis_set
      IF (PRESENT(orb_basis_set)) THEN
        IF(softb_local) THEN
          orb_basis_set => atomic_kind%soft_basis_set
        ELSE
          orb_basis_set => atomic_kind%orb_basis_set
        END IF
      END IF
      IF (PRESENT(element_symbol)) element_symbol = atomic_kind%element_symbol
      IF (PRESENT(name)) name = atomic_kind%name
      IF (PRESENT(mass)) mass = atomic_kind%mass
      IF (PRESENT(kind_number)) kind_number = atomic_kind%kind_number
      IF (PRESENT(natom)) natom = atomic_kind%natom
      IF (PRESENT(atom_list)) atom_list => atomic_kind%atom_list
      IF (PRESENT(alpha_core_charge)) THEN
        IF (ASSOCIATED(atomic_kind%all_potential)) THEN
          CALL get_potential(potential=atomic_kind%all_potential,&
                             alpha_core_charge=alpha_core_charge)
        ELSE IF (ASSOCIATED(atomic_kind%elp_potential)) THEN
          CALL get_potential(potential=atomic_kind%elp_potential,&
                             alpha_core_charge=alpha_core_charge)
        ELSE IF (ASSOCIATED(atomic_kind%gth_potential)) THEN
          CALL get_potential(potential=atomic_kind%gth_potential,&
                             alpha_core_charge=alpha_core_charge)
        ELSE
          alpha_core_charge = 0.0_dp
        END IF
      END IF
      IF (PRESENT(ccore_charge)) THEN
        IF (ASSOCIATED(atomic_kind%all_potential)) THEN
          CALL get_potential(potential=atomic_kind%all_potential,&
                             ccore_charge=ccore_charge)
        ELSE IF (ASSOCIATED(atomic_kind%elp_potential)) THEN
          CALL get_potential(potential=atomic_kind%elp_potential,&
                             ccore_charge=ccore_charge)
        ELSE IF (ASSOCIATED(atomic_kind%gth_potential)) THEN
          CALL get_potential(potential=atomic_kind%gth_potential,&
                             ccore_charge=ccore_charge)
        ELSE
          ccore_charge = 0.0_dp
        END IF
      END IF
      IF (PRESENT(core_charge_radius)) THEN
        IF (ASSOCIATED(atomic_kind%all_potential)) THEN
          CALL get_potential(potential=atomic_kind%all_potential,&
                             core_charge_radius=core_charge_radius)
        ELSE IF (ASSOCIATED(atomic_kind%elp_potential)) THEN
          CALL get_potential(potential=atomic_kind%elp_potential,&
                             core_charge_radius=core_charge_radius)
        ELSE IF (ASSOCIATED(atomic_kind%gth_potential)) THEN
          CALL get_potential(potential=atomic_kind%gth_potential,&
                             core_charge_radius=core_charge_radius)
        ELSE
          core_charge_radius = 0.0_dp
        END IF
      END IF
      IF (PRESENT(ncgf)) THEN
        IF (ASSOCIATED(atomic_kind%orb_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=atomic_kind%orb_basis_set,&
                                 ncgf=ncgf)
        ELSE
          ncgf = 0
        END IF
      END IF
      IF (PRESENT(nsgf)) THEN
        IF (ASSOCIATED(atomic_kind%orb_basis_set)) THEN
          CALL get_gto_basis_set(gto_basis_set=atomic_kind%orb_basis_set,&
                                 nsgf=nsgf)
        ELSE
          nsgf = 0
        END IF
      END IF
      IF (PRESENT(z)) THEN
        IF (ASSOCIATED(atomic_kind%all_potential)) THEN
          CALL get_potential(potential=atomic_kind%all_potential,z=z)
        ELSE IF (ASSOCIATED(atomic_kind%gth_potential)) THEN
          CALL get_potential(potential=atomic_kind%gth_potential,z=z)
        ELSE
          z = 0
        END IF
      END IF
      IF (PRESENT(zeff)) THEN
        IF (ASSOCIATED(atomic_kind%all_potential)) THEN
          CALL get_potential(potential=atomic_kind%all_potential,zeff=zeff)
        ELSE IF (ASSOCIATED(atomic_kind%gth_potential)) THEN
          CALL get_potential(potential=atomic_kind%gth_potential,zeff=zeff)
        ELSE
          zeff = 0
        END IF
      END IF
      IF (PRESENT(qeff)) THEN
        IF (ASSOCIATED(atomic_kind%elp_potential)) THEN
          CALL get_potential(potential=atomic_kind%elp_potential,qeff=qeff)
        ELSE
          qeff = 0
        END IF
      END IF
      IF (PRESENT(ldaCoreRadius)) ldaCoreRadius = atomic_kind%ldaCoreRadius
      IF (PRESENT(soft_basis_set)) soft_basis_set => atomic_kind%soft_basis_set
      IF (PRESENT(hard_basis_set)) hard_basis_set => atomic_kind%hard_basis_set
      IF (PRESENT(paw_proj_set))   paw_proj_set => atomic_kind%paw_proj_set
      IF (PRESENT(paw_atom)) paw_atom = atomic_kind%paw_atom
      IF (PRESENT(qm_atom)) qm_atom = atomic_kind%qm_atom
      IF (PRESENT(hard_radius)) hard_radius = atomic_kind%hard_radius
      IF (PRESENT(hard0_radius)) hard0_radius = atomic_kind%hard0_radius
      IF (PRESENT(harmonics))  harmonics => atomic_kind%harmonics
      IF (PRESENT(max_s_harm)) THEN 
         IF(ASSOCIATED(atomic_kind%harmonics)) THEN 
            max_s_harm = atomic_kind%harmonics%max_s_harm
         ELSE
            max_s_harm = 0
         END IF
      END IF
      IF (PRESENT(max_iso_not0)) THEN 
         IF(ASSOCIATED(atomic_kind%harmonics)) THEN 
            max_iso_not0 = atomic_kind%harmonics%max_iso_not0
         ELSE
            max_iso_not0 = 0
         END IF
      END IF
      IF (PRESENT(grid_atom))  grid_atom => atomic_kind%grid_atom
      IF (PRESENT(ngrid_ang))  ngrid_ang = atomic_kind%ngrid_ang
      IF (PRESENT(ngrid_rad))  ngrid_rad = atomic_kind%ngrid_rad
      IF (PRESENT(lmax_rho0))  lmax_rho0 = atomic_kind%lmax_rho0
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind is not associated")

    END IF

  END SUBROUTINE get_atomic_kind

! *****************************************************************************

  SUBROUTINE get_atomic_kind_set(atomic_kind_set,maxatom,maxcgf,&
                                 maxco,maxco_proj,maxgtop,maxgtops,maxlgto,&
                                 maxlppl,maxlppnl,maxnset,maxpgf,&
                                 maxppnl,maxsgf,maxsgf_set,maxshell,&
                                 natom,nelectron,ncgf,ncgf_aux,&
                                 npgf,nset,nsgf,nshell,&
                                 all_potential_present,gth_potential_present,&
                                 elp_potential_present,paw_atom_present,&
                                 qm_atom_present,&
                                 atom_of_kind,kind_of,natom_of_kind,maxder,&
                                 max_ldaCoreRadius,max_ngrid_rad,&
                                 max_sph_harm,maxg_iso_not0,lmax_rho0,&
                                 zetsoft_max)

!   Purpose: Get informations about an atomic kind set.

!   History: - Creation (14.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(OUT), OPTIONAL :: maxatom, maxcgf, maxco, maxco_proj, &
      maxgtop,maxgtops, maxlgto, maxlppl, maxlppnl, maxnset,maxpgf,  &
      maxppnl, maxsgf, maxsgf_set, maxshell, natom, nelectron, &
      ncgf, ncgf_aux, npgf, nset, nsgf, nshell, &
      max_ngrid_rad, max_sph_harm, maxg_iso_not0
    LOGICAL, INTENT(OUT), OPTIONAL           :: all_potential_present, &
                                                gth_potential_present, &
                                                elp_potential_present, &
                                                paw_atom_present, &
                                                qm_atom_present
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: atom_of_kind, kind_of, &
                                                natom_of_kind
    INTEGER, INTENT(IN), OPTIONAL            :: maxder
    REAL(kind=dp), DIMENSION(4), &
      INTENT(out), OPTIONAL                  :: max_ldaCoreRadius
    REAL(kind=dp),  &
      INTENT(out), OPTIONAL                  :: zetsoft_max
    INTEGER, INTENT(out), OPTIONAL           :: lmax_rho0

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_atomic_kind_set"

    INTEGER :: atom_a, iatom, ikind, imax, ipgf, iset, &
      lmax_rho0_kind, max_iso_not0, &
      max_s_harm, n, ngrid_rad, nkind, p_excl_bas, &
      p_excl_bas_max, p_full_bas, p_full_bas_max, p_full_bas_tot, &
      p_l_angles_max, p_min_bas, p_min_bas_max, p_min_bas_tot, zeff
    INTEGER, DIMENSION(:), POINTER           :: n_ppl
    INTEGER, DIMENSION(:), POINTER           :: npgf_kind
    LOGICAL                                  :: paw_atom
    LOGICAL                                  :: qm_atom
    REAL(kind=dp), DIMENSION(4)              :: ldaCoreRadius
    REAL(kind=dp), DIMENSION(:,:), POINTER   :: zet_kind
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set, orb_basis_set
    TYPE(gto_basis_set_type), POINTER        :: soft_basis_set
    TYPE(paw_proj_set_type),  POINTER        :: paw_proj_set

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN

      IF (PRESENT(maxatom)) maxatom = 0
      IF (PRESENT(maxcgf)) maxcgf = 0
      IF (PRESENT(maxco)) maxco = 0
      IF (PRESENT(maxco_proj)) maxco_proj = 0
      IF (PRESENT(maxg_iso_not0)) maxg_iso_not0 = 0
      IF (PRESENT(maxgtop))  maxgtop = 0
      IF (PRESENT(maxgtops)) maxgtops = 0
      IF (PRESENT(maxlgto)) maxlgto = -1
      IF (PRESENT(maxlppl)) maxlppl = -1
      IF (PRESENT(maxlppnl)) maxlppnl = -1
      IF (PRESENT(maxnset)) maxnset = 0
      IF (PRESENT(maxpgf)) maxpgf = 0
      IF (PRESENT(maxppnl)) maxppnl = 0
      IF (PRESENT(maxsgf)) maxsgf = 0
      IF (PRESENT(maxsgf_set)) maxsgf_set = 0
      IF (PRESENT(maxshell)) maxshell = 0
      IF (PRESENT(natom)) natom = 0
      IF (PRESENT(ncgf)) ncgf = 0
      IF (PRESENT(ncgf_aux)) ncgf_aux = 0
      IF (PRESENT(nelectron)) nelectron = 0
      IF (PRESENT(npgf)) npgf = 0
      IF (PRESENT(nset)) nset = 0
      IF (PRESENT(nsgf)) nsgf = 0
      IF (PRESENT(nshell)) nshell = 0
      IF (PRESENT(all_potential_present)) all_potential_present = .FALSE.
      IF (PRESENT(elp_potential_present)) elp_potential_present = .FALSE.
      IF (PRESENT(gth_potential_present)) gth_potential_present = .FALSE.
      IF (PRESENT(paw_atom_present)) paw_atom_present = .FALSE.
      IF (PRESENT(qm_atom_present)) qm_atom_present = .FALSE.
      IF (PRESENT(atom_of_kind)) atom_of_kind(:) = 0
      IF (PRESENT(kind_of)) kind_of(:) = 0
      IF (PRESENT(natom_of_kind)) natom_of_kind(:) = 0
      IF (PRESENT(zetsoft_max))  zetsoft_max = 0.0_dp
      IF (PRESENT(max_ldaCoreRadius)) THEN
         max_ldaCoreRadius=0.0_dp
      END IF
      IF (PRESENT(max_ngrid_rad))  max_ngrid_rad = 0
      IF (PRESENT(max_sph_harm))   max_sph_harm = 0
      IF (PRESENT(lmax_rho0))      lmax_rho0 = 0

      nkind = SIZE(atomic_kind_set)
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             all_potential=all_potential,&
                             elp_potential=elp_potential,&
                             gth_potential=gth_potential,&
                             aux_basis_set=aux_basis_set,&
                             orb_basis_set=orb_basis_set,&
                             soft_basis_set=soft_basis_set,&
                             ldaCoreRadius=ldaCoreRadius,&
                             ngrid_rad=ngrid_rad,&
                             max_s_harm=max_s_harm, &
                             max_iso_not0=max_iso_not0,&
                             paw_atom=paw_atom,&
                             qm_atom=qm_atom,&
                             lmax_rho0=lmax_rho0_kind)
        IF (PRESENT(max_ldaCoreRadius)) THEN
           max_ldaCoreRadius=MAX(max_ldaCoreRadius,ldaCoreRadius)
        END IF
        IF (PRESENT(maxatom)) THEN
          maxatom = MAX(maxatom,atomic_kind%natom)
        END IF
        IF (PRESENT(maxlppl).AND.ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,nexp_ppl=n)
          maxlppl = MAX(maxlppl,2*(n - 1))
        END IF
        IF (PRESENT(maxlppl).AND.ASSOCIATED(elp_potential)) THEN
          CALL get_potential(potential=elp_potential,nexp_ppl=n_ppl)
          n=MAXVAL(n_ppl)
          maxlppl = MAX(maxlppl,2*(n - 1))
        END IF
        IF (PRESENT(maxlppnl).AND.ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,lprj_ppnl_max=imax)
          maxlppnl = MAX(maxlppnl,imax)
        END IF
        IF (ASSOCIATED(elp_potential)) THEN
          CALL get_potential(potential=elp_potential)
        END IF
        IF (PRESENT(maxcgf)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,ncgf=imax)
            maxcgf = MAX(maxcgf,imax)
          END IF
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,ncgf=imax)
            maxcgf = MAX(maxcgf,imax)
          END IF
        END IF
        IF (PRESENT(maxco)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            IF (PRESENT(maxder)) THEN
              CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                                     maxco=imax,maxder=maxder)
            ELSE
              CALL get_gto_basis_set(gto_basis_set=aux_basis_set,maxco=imax)
            END IF
            maxco = MAX(maxco,imax)
          END IF
          IF (ASSOCIATED(orb_basis_set)) THEN
            IF (PRESENT(maxder)) THEN
              CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                     maxco=imax,maxder=maxder)
            ELSE
              CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxco=imax)
            END IF
            maxco = MAX(maxco,imax)
          END IF
          IF (ASSOCIATED(gth_potential)) THEN
            CALL get_potential(potential=gth_potential,lprj_ppnl_max=imax)
            maxco = MAX(maxco,ncoset(imax))
          END IF
        END IF
        IF (PRESENT(maxco_proj)) THEN
          IF (ASSOCIATED(paw_proj_set)) THEN
            CALL get_paw_proj_set(paw_proj_set=paw_proj_set,ncgauprj=imax)
            maxco_proj = MAX(maxco_proj,imax)
          END IF
        ENDIF
        IF (PRESENT(maxgtop)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxco=imax,&
                              nset=n)
            maxgtop = MAX(maxgtop, n*imax)     
          END IF
        END IF
        IF (PRESENT(maxgtops)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxso=imax,&
                              nset=n)
            maxgtops = MAX(maxgtops, n*imax)
          END IF
        END IF
        IF (PRESENT(maxlgto)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,maxl=imax)
            maxlgto = MAX(maxlgto,imax)
          END IF
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxl=imax)
            maxlgto = MAX(maxlgto,imax)
          END IF
        END IF
        IF (PRESENT(maxnset)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                              nset=n)
            maxnset = MAX(maxnset, n)     
          END IF
        END IF
        IF (PRESENT(maxpgf)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,maxpgf=imax)
            maxpgf = MAX(maxpgf,imax)
          END IF
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxpgf=imax)
            maxpgf = MAX(maxpgf,imax)
          END IF
        END IF
        IF (PRESENT(maxppnl).AND.ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,nppnl=imax)
          maxppnl = MAX(maxppnl,imax)
        END IF
        IF (PRESENT(maxsgf)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,nsgf=imax)
            maxsgf = MAX(maxsgf,imax)
          END IF
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=imax)
            maxsgf = MAX(maxsgf,imax)
          END IF
        END IF
        IF (PRESENT(maxsgf_set)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,maxsgf_set=imax)
            maxsgf_set = MAX(maxsgf_set,imax)
          END IF
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxsgf_set=imax)
            maxsgf_set = MAX(maxsgf_set,imax)
          END IF
        END IF
        IF (PRESENT(maxshell)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,maxshell=imax)
            maxshell = MAX(maxshell,imax)
          END IF
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxshell=imax)
            maxshell = MAX(maxshell,imax)
          END IF
        END IF
        IF (PRESENT(natom)) THEN
          natom = natom + atomic_kind_set(ikind)%natom
        END IF
        IF (PRESENT(ncgf)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,ncgf=n)
            ncgf = ncgf + n*atomic_kind_set(ikind)%natom
          END IF
        END IF
        IF (PRESENT(ncgf_aux)) THEN
          IF (ASSOCIATED(aux_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=aux_basis_set,ncgf=n)
            ncgf_aux = ncgf_aux + n*atomic_kind_set(ikind)%natom
          END IF
        END IF
        IF (PRESENT(nelectron)) THEN
          IF (ASSOCIATED(atomic_kind%all_potential)) THEN
            CALL get_potential(potential=atomic_kind%all_potential,zeff=zeff)
          ELSE IF (ASSOCIATED(atomic_kind%gth_potential)) THEN
            CALL get_potential(potential=atomic_kind%gth_potential,zeff=zeff)
          ELSE
            zeff = 0
          END IF
          nelectron = nelectron + atomic_kind_set(ikind)%natom*zeff
        END IF
        IF (PRESENT(npgf)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,npgf_sum=n)
            npgf = npgf + n*atomic_kind_set(ikind)%natom
          END IF
        END IF
        IF (PRESENT(nset)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nset=n)
            nset = nset + n*atomic_kind_set(ikind)%natom
          END IF
        END IF
        IF (PRESENT(nsgf)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=n)
            nsgf = nsgf + n*atomic_kind_set(ikind)%natom
          END IF
        END IF
        IF (PRESENT(nshell)) THEN
          IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nshell_sum=n)
            nshell = nshell + n*atomic_kind_set(ikind)%natom
          END IF
        END IF
        IF (PRESENT(all_potential_present)) THEN
          IF (ASSOCIATED(all_potential)) THEN
            all_potential_present = .TRUE.
          END IF
        END IF
        IF (PRESENT(elp_potential_present)) THEN
          IF (ASSOCIATED(elp_potential)) THEN
            elp_potential_present = .TRUE.
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
        IF (PRESENT(qm_atom_present)) THEN
          IF (qm_atom) THEN
            qm_atom_present = .TRUE.
          END IF
        END IF
        IF (PRESENT(atom_of_kind)) THEN
          DO iatom=1,atomic_kind%natom
            atom_a = atomic_kind%atom_list(iatom)
            atom_of_kind(atom_a) = iatom
          END DO
        END IF
        IF (PRESENT(kind_of)) THEN
          DO iatom=1,atomic_kind%natom
            atom_a = atomic_kind%atom_list(iatom)
            kind_of(atom_a) = ikind
          END DO
        END IF
        IF (PRESENT(natom_of_kind)) THEN
          natom_of_kind(ikind) = atomic_kind_set(ikind)%natom
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
        IF( PRESENT(zetsoft_max) ) THEN
          CALL  get_gto_basis_set(gto_basis_set=soft_basis_set, &
                nset = n, npgf = npgf_kind, zet = zet_kind )
          DO iset = 1,n
            DO ipgf = 1,npgf_kind(iset)
               zetsoft_max= MAX(zetsoft_max,zet_kind(ipgf,iset))
            END DO
          END DO 
        END IF
      END DO
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind_set is not associated")

    END IF

  END SUBROUTINE get_atomic_kind_set

! *****************************************************************************

  SUBROUTINE init_atomic_kind(atomic_kind,globenv,method)

!   Purpose: Initialise an atomic kind data set.

!   History: - Creation (11.01.2002,MK)
!              20.09.2002 adapted for pol/kg use, gtb

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    CHARACTER(LEN=10),INTENT(IN), OPTIONAL   :: method


    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_atomic_kind"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind)) THEN

      IF (globenv%program_name == "QS") THEN
        IF (ASSOCIATED(atomic_kind%gth_potential))THEN
          CALL init_potential(atomic_kind%gth_potential)
        END IF
        IF (ASSOCIATED(atomic_kind%orb_basis_set))THEN
          IF(atomic_kind%orb_basis_set%norm_type < 0) &
                        atomic_kind%orb_basis_set%norm_type = 2
          CALL init_orb_basis_set(atomic_kind%orb_basis_set)
        END IF
      ELSE IF (globenv%program_name == "POL") THEN
        IF (ASSOCIATED(atomic_kind%orb_basis_set))THEN
          IF(atomic_kind%orb_basis_set%norm_type < 0) &
                        atomic_kind%orb_basis_set%norm_type = 1
          CALL init_orb_basis_set(atomic_kind%orb_basis_set)
        END IF
        IF (ASSOCIATED(atomic_kind%aux_basis_set))THEN
          IF(atomic_kind%aux_basis_set%norm_type < 0) &
                        atomic_kind%aux_basis_set%norm_type = 1
          CALL init_aux_basis_set(atomic_kind%aux_basis_set)
        END IF
      ELSE IF (globenv%program_name == "KG") THEN
        IF (ASSOCIATED(atomic_kind%gth_potential).AND. PRESENT(method))THEN
          IF(method == "KG_GPW")&  
          CALL init_potential(atomic_kind%gth_potential)
        END IF

        IF (ASSOCIATED(atomic_kind%orb_basis_set))THEN
          IF(atomic_kind%orb_basis_set%norm_type < 0) THEN
            IF(PRESENT(method)) THEN
              IF(method == "KG_GPW") THEN
                atomic_kind%orb_basis_set%norm_type = 2
              ELSE
                atomic_kind%orb_basis_set%norm_type = 1
              END IF
            ELSE
              atomic_kind%orb_basis_set%norm_type = 1
            END IF
          END IF
          CALL init_orb_basis_set(atomic_kind%orb_basis_set)
        END IF
        IF (ASSOCIATED(atomic_kind%aux_basis_set))THEN
          IF(atomic_kind%aux_basis_set%norm_type < 0) &
                        atomic_kind%aux_basis_set%norm_type = 1
          CALL init_aux_basis_set(atomic_kind%aux_basis_set)
        END IF
      ELSE IF (globenv%program_name == "FIST") THEN
        IF (ASSOCIATED(atomic_kind%elp_potential))THEN
          CALL allocate_potential(atomic_kind%elp_potential)
        END IF
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind is not associated",globenv)

    END IF

  END SUBROUTINE init_atomic_kind

! *****************************************************************************

  SUBROUTINE init_atomic_kind_set(atomic_kind_set,globenv,method)

!   Purpose: Initialise an atomic kind set data set.

!   History: - Creation (17.01.2002,MK)
!            - 20.09.2002 globenv passed (gt)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    CHARACTER(LEN=10),INTENT(IN), OPTIONAL   :: method

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_atomic_kind_set"

    INTEGER                                  :: ikind, nkind
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN

      nkind = SIZE(atomic_kind_set)

      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL init_atomic_kind(atomic_kind,globenv,method=method)
      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind_set is not associated",&
                        globenv)

    END IF

  END SUBROUTINE init_atomic_kind_set

! *****************************************************************************

  SUBROUTINE init_gapw_basis_set(atomic_kind_set,gapw,globenv)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(gapw_control_type),  POINTER        :: gapw
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_gapw_basis_set"

    INTEGER                                  :: ikind, nkind, output_unit
    LOGICAL                                  :: paw_atom
    REAL(dp)                                 :: covalent_radius,rc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis, soft_basis
    TYPE(gto_basis_set_type), POINTER        :: hard_basis
    TYPE(paw_proj_set_type), POINTER         :: paw_proj
    CHARACTER(LEN=2)                         :: current_symbol,kind_symbol

! *****************************************************************************

    IF (ASSOCIATED(atomic_kind_set)) THEN
      nkind = SIZE(atomic_kind_set)

      IF (globenv%print%basis_set) THEN
       output_unit = globenv%scr
       IF(globenv%ionode) WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
       "SOFT BASIS SET INFORMATION (Unnormalised Gaussian-type functions)"
      END IF

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL allocate_gto_basis_set(atomic_kind%soft_basis_set)
!        CALL allocate_gto_basis_set(atomic_kind%hard_basis_set)

        atomic_kind%hard_basis_set => atomic_kind%orb_basis_set

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis,&
                             hard_radius=rc)

!        bsname = orb_basis%name
!        atomic_kind%hard_basis_set%name = TRIM(bsname)//"_hard"

        paw_atom = .FALSE.

        CALL create_soft_basis(orb_basis,atomic_kind%soft_basis_set,gapw,rc,paw_atom)
        CALL set_atomic_kind(atomic_kind=atomic_kind,&
                             paw_atom=paw_atom)

        IF(paw_atom) THEN
          IF (globenv%print%basis_set) THEN
            IF(globenv%ionode) WRITE (UNIT=output_unit,FMT="(/,T2,I2,A)")&
              atomic_kind%kind_number,". Atomic kind: "//TRIM(atomic_kind%name)
          END IF

          CALL write_gto_basis_set(atomic_kind%soft_basis_set,globenv)
!         the hard_basis points to the orb_basis
          CALL allocate_paw_proj_set(atomic_kind%paw_proj_set)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               paw_proj_set=paw_proj)

          CALL projectors(paw_proj,orb_basis,rc,gapw,globenv)

         ENDIF

!       grid_atom and harmonics are allocated even if NOT PAW_ATOM
        CALL allocate_grid_atom(atomic_kind%grid_atom)
        CALL allocate_harmonics_atom(atomic_kind%harmonics)

      ENDDO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
               "The pointer atomic_kind_set is not associated",&
                globenv)

    ENDIF
  END SUBROUTINE init_gapw_basis_set

! *****************************************************************************

  SUBROUTINE read_atomic_kind(atomic_kind,globenv)

!   Purpose: Read an atomic kind data set from the input file.

!   History: - Creation (09.02.2002,MK)
!            - 20.09.2002,gt: adapted for POL/KG use (elp_potential)

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_atomic_kind", &
      routineP = module_name//':'//routine_name

    CHARACTER(LEN=2)                         :: current_symbol, element_symbol
    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=default_string_length) :: aux_basis_set_name, keyword, &
                                                norm_sym, orb_basis_set_name, &
                                                potential_name, potential_type
    INTEGER                                  :: i, ipos, &
                                                min_b_index, naux, &
                                                ngrid_ang, ngrid_rad, norb, z
    LOGICAL                                  :: failure, section_found
    TYPE(cp_error_type)                      :: error

!   ---------------------------------------------------------------------------

    failure = .FALSE.
    aux_basis_set_name = ""
    orb_basis_set_name = ""
    potential_name = ""
    potential_type = ""
    z = -1
    atomic_kind%mass=-1.0_dp

    CALL cp_error_init(error)

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="KIND",&
                      start_section_name=atomic_kind%name,&
                      section_found=section_found)
    IF (.NOT.section_found) THEN
       CALL finish_parser()
       IF (z == -1) THEN
          ipos = INDEX(atomic_kind%name,"_")
          IF ((ipos == 2).OR.(ipos == 3)) THEN
             element_symbol = atomic_kind%name(1:ipos-1)
          ELSE
             element_symbol = atomic_kind%name
          END IF
          CALL uppercase(element_symbol)
          DO i=LBOUND(ptable,1),UBOUND(ptable,1)
             current_symbol = ptable(i)%symbol
             CALL uppercase(current_symbol)
             IF (current_symbol == element_symbol) THEN
                atomic_kind%element_symbol = ptable(i)%symbol
                IF (atomic_kind%mass==-1.0_dp) THEN
                   atomic_kind%mass = ptable(i)%amass
                END IF
                z = i
                EXIT
             END IF
          END DO
       END IF
       IF (z /= -1) THEN
          CALL start_parser(file_name=globenv%input_file_name,&
                            globenv=globenv,&
                            start_section_label="KIND",&
                            start_section_name=atomic_kind%element_symbol,&
                            section_found=section_found)
       END IF
    END IF
    IF (.NOT.section_found) THEN
       CALL finish_parser() 
       CALL start_parser(file_name=globenv%input_file_name,&
                         globenv=globenv,&
                         start_section_label="KIND",&
                         start_section_name="DEFAULT",&
                         section_found=section_found)
    END IF

    IF (section_found) THEN
      DO
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine_name,module_name,__LINE__,"EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("AUXILIARY_BASIS_SET","AUX_BASIS_SET","AUX_BASIS")
            CALL read_object(aux_basis_set_name)
            test_result = test_object()
            naux = -1
            IF (test_result /= "EOL") THEN
              CALL read_object(norm_sym,lower_to_upper=.TRUE.)
              SELECT CASE (TRIM(norm_sym))
              CASE ("NONE","NO")
                naux = 0
              CASE ("WFN","WAVEFUNCTION","FUNCTION")
                naux = 1
              CASE ("DENS","DENSITY","SQUARE")
                naux = 2
              CASE DEFAULT
                CALL stop_parser(routine_name,module_name,__LINE__,&
                     "Invalid normalization method")
              END SELECT
            END IF
          CASE ("ORBITAL_BASIS_SET","ORB_BASIS_SET","ORB_BASIS","BASIS_SET","BAS")
            CALL read_object(orb_basis_set_name)
            test_result = test_object()
            norb = -1
            IF (test_result /= "EOL") THEN
              CALL read_object(norm_sym,lower_to_upper=.TRUE.)
              SELECT CASE (TRIM(norm_sym))
              CASE ("NONE","NO")
                norb = 0
              CASE ("WFN","WAVEFUNCTION","FUNCTION")
                norb = 1
              CASE ("DENS","DENSITY","SQUARE")
                norb = 2
              CASE DEFAULT
                CALL stop_parser(routine_name,module_name,__LINE__,&
                     "Invalid normalization method")
              END SELECT
            END IF
          CASE ("ELEMENT_SYMBOL","ELEMENT")
            CALL read_object(element_symbol,lower_to_upper=.TRUE.)
            DO i=LBOUND(ptable,1),UBOUND(ptable,1)
              current_symbol = ptable(i)%symbol
              CALL uppercase(current_symbol)
              IF (current_symbol == element_symbol) THEN
                atomic_kind%element_symbol = ptable(i)%symbol
                if (atomic_kind%mass==-1.0_dp) then
                   atomic_kind%mass = ptable(i)%amass
                end if
                z = i
                EXIT
              END IF
            END DO
          CASE ("ATOMIC_MASS","ATOMIC_WEIGHT","MASS","WEIGHT")
            CALL read_object(atomic_kind%mass)
          CASE ("POTENTIAL","POT")
            CALL read_object(potential_name)
          CASE ("LDA_CORE_RADIUS")
             DO i=1,SIZE(atomic_kind%ldaCoreRadius)
                CALL read_object(atomic_kind%ldaCoreRadius(i))
             END DO
          CASE ("HARD_EXP_RADIUS")
             CALL read_object(atomic_kind%hard_radius)
          CASE ("RHO0_EXP_RADIUS")
             CALL read_object(atomic_kind%hard0_radius)
             CALL cp_assert(atomic_kind%hard_radius>atomic_kind%hard0_radius,&
                       cp_failure_level,cp_assertion_failed,routineP,&
                      "rc0 should be <= rc",error,failure)
          CASE ("LEBEDEV_GRID")
             CALL read_object(atomic_kind%ngrid_ang)
             CALL cp_assert(atomic_kind%ngrid_ang > 0,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "# point lebedev grid < 0",error,failure)
          CASE ("RADIAL_GRID")
             CALL read_object(atomic_kind%ngrid_rad)
             CALL cp_assert(atomic_kind%ngrid_rad > 0,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "# point radial grid < 0",error,failure)
          CASE DEFAULT
            CALL stop_parser(routine_name,module_name,__LINE__,"INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    END IF

    CALL finish_parser()

!   *** The corresponding element of the current atomic kind is still   ***
!   *** undefined. Try the head first two characters of the kind name itself ***

    IF (z == -1) THEN
      ipos = INDEX(atomic_kind%name,"_")
      IF ((ipos == 2).OR.(ipos == 3)) THEN
        element_symbol = atomic_kind%name(1:ipos-1)
      ELSE
        element_symbol = atomic_kind%name
      END IF
      CALL uppercase(element_symbol)
      DO i=LBOUND(ptable,1),UBOUND(ptable,1)
        current_symbol = ptable(i)%symbol
        CALL uppercase(current_symbol)
        IF (current_symbol == element_symbol) THEN
          atomic_kind%element_symbol = ptable(i)%symbol
          if (atomic_kind%mass==-1.0_dp) then
             atomic_kind%mass = ptable(i)%amass
          end if
          z = i
          EXIT
        END IF
      END DO
    END IF

    atomic_kind%mass = atomic_kind%mass*massunit

    IF (z == -1) THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "No element could be assign to the atomic kind <"//&
                        TRIM(atomic_kind%name)//">. Specify a valid "//&
                        "element symbol in the input section <KIND "//&
                        TRIM(atomic_kind%name)//"> or select a valid "//&
                        "element symbol for the atomic kind name",globenv)
    END IF

!   *** Allocate and initialise the orbital basis set data set structure ***
    CALL init_orbital_pointers(5) ! debug the SUN optimizer

    IF (orb_basis_set_name /= '') THEN
      CALL allocate_gto_basis_set(atomic_kind%orb_basis_set)
      CALL read_gto_basis_set(atomic_kind%element_symbol,orb_basis_set_name,&
                              atomic_kind%orb_basis_set,globenv)
      CALL set_gto_basis_set(gto_basis_set=atomic_kind%orb_basis_set,&
                             norm_type=norb)
    END IF

!   *** Allocate and initialise the auxiliary basis set data set structure ***

    IF (aux_basis_set_name /= '') THEN
      CALL allocate_gto_basis_set(atomic_kind%aux_basis_set)
      CALL read_gto_basis_set(atomic_kind%element_symbol,aux_basis_set_name,&
                              atomic_kind%aux_basis_set,globenv)
      CALL set_gto_basis_set(gto_basis_set=atomic_kind%aux_basis_set,&
                             norm_type=naux)
    END IF

!   *** Allocate and initialise the potential data set structure ***

    IF (potential_name /= '') THEN

      ipos = INDEX(potential_name,"-")

      IF (ipos > 1) THEN
        potential_type = potential_name(:ipos-1)
      ELSE
        potential_type = potential_name
      END IF

      CALL uppercase(potential_type)

      SELECT CASE (TRIM(potential_type))
      CASE ("ALL")
        CALL allocate_potential(atomic_kind%all_potential)
        CALL read_potential(atomic_kind%element_symbol,potential_name,&
                            atomic_kind%all_potential,globenv)
        CALL set_potential(atomic_kind%all_potential,z=z)
      CASE ("GTH")
        CALL allocate_potential(atomic_kind%gth_potential)
        CALL read_potential(atomic_kind%element_symbol,potential_name,&
                            atomic_kind%gth_potential,globenv)
        CALL set_potential(atomic_kind%gth_potential,z=z)
      CASE ("GHOST")
        atomic_kind%mass = 0.0_dp
      CASE ("ELP")
        CALL allocate_potential(atomic_kind%elp_potential)
        CALL read_potential(atomic_kind%element_symbol,potential_name,&
                            atomic_kind%elp_potential,globenv)
      CASE DEFAULT
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "An invalid potential type <"//&
                          TRIM(potential_name)//"> was specified "//&
                          "for the atomic kind <"//&
                          TRIM(atomic_kind%name)//&
                          "> in the input file <"//&
                          TRIM(globenv%input_file_name)//">",globenv)
      END SELECT

    ELSE

      CALL allocate_potential(atomic_kind%all_potential)
      CALL set_default_all_potential(atomic_kind%all_potential,z)

    END IF

    CALL cp_error_dealloc_ref(error)

  END SUBROUTINE read_atomic_kind

! *****************************************************************************

  SUBROUTINE read_atomic_kind_set(atomic_kind_set,globenv)

!   Purpose: Read an atomic kind set data set from the input file.

!   History: - Creation (09.02.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_atomic_kind_set"

    INTEGER                                  :: ikind, nkind
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN

      nkind = SIZE(atomic_kind_set)

      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL read_atomic_kind(atomic_kind,globenv)
      END DO
    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind_set is not associated",&
                        globenv)

    END IF

  END SUBROUTINE read_atomic_kind_set

! *****************************************************************************
  SUBROUTINE check_atomic_kind(atomic_kind,dft_control,globenv)

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type),INTENT(IN)        :: dft_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "check_atomic_kind", &
      routineP = module_name//':'//routine_name

    REAL(dp)                             :: mass
    INTEGER                              :: z
    CHARACTER(LEN=default_string_length) :: name
    TYPE(gto_basis_set_type), POINTER    :: orb_basis_set
    TYPE(all_potential_type), POINTER    :: all_potential
    TYPE(elp_potential_type), POINTER    :: elp_potential
    TYPE(gth_potential_type), POINTER    :: gth_potential
    TYPE(semi_empirical_type), POINTER   :: se_parameter
    LOGICAL                              :: defined
!   ---------------------------------------------------------------------------

    NULLIFY ( orb_basis_set, all_potential, elp_potential, gth_potential )
    CALL get_atomic_kind(atomic_kind,name=name,z=z,mass=mass,&
                         orb_basis_set=orb_basis_set,&
                         all_potential=all_potential,&
                         elp_potential=elp_potential,&
                         gth_potential=gth_potential)

    IF (globenv%program_name == "QS") THEN
      IF (z < 1) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                     "No valid element is assign to the atomic kind ",globenv)
      END IF
      IF (.NOT. dft_control%qs_control%semi_empirical) THEN
        IF ( .NOT. ASSOCIATED(gth_potential) .AND. &
             .NOT. ASSOCIATED(all_potential) .AND. &
             .NOT. ASSOCIATED(elp_potential) .AND. &
             mass > 0._dp ) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                      "No potential type was defined for the "//&
                      "atomic kind <"//TRIM(atomic_kind%name),globenv)
        END IF
        IF ( .NOT. ASSOCIATED(orb_basis_set) ) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                      "No basis set type was defined for the "//&
                      "atomic kind <"//TRIM(atomic_kind%name),globenv)
        END IF
      ELSE
        CALL get_atomic_kind(atomic_kind,se_parameter=se_parameter)
        IF ( .NOT. ASSOCIATED(se_parameter) ) THEN
          CALL allocate_se_param(atomic_kind%se_parameter)
          CALL get_atomic_kind(atomic_kind,se_parameter=se_parameter)
        END IF
        CALL get_se_param(se_parameter,defined=defined)
        IF (.NOT.defined) THEN
          ! set default parameters 
          CALL se_param_set_default(se_parameter,z,&
                 dft_control%qs_control%method)
        END IF
        CALL write_se_param(se_parameter,globenv)
        CALL init_se_param(atomic_kind%se_parameter,&
                           atomic_kind%orb_basis_set)
        IF ( ASSOCIATED(all_potential) ) THEN
          CALL init_potential ( all_potential, itype="BARE", &
             zeff=se_parameter%zeff)
        END IF
      END IF
    ELSE IF (globenv%program_name == "KG") THEN
      !no checks programmed yet
    END IF

  END SUBROUTINE check_atomic_kind

! *****************************************************************************

  SUBROUTINE check_atomic_kind_set(atomic_kind_set,dft_control,globenv)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type),INTENT(IN)        :: dft_control
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "check_atomic_kind_set"

    INTEGER                                  :: ikind, nkind
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN
      nkind = SIZE(atomic_kind_set)
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL check_atomic_kind(atomic_kind,dft_control,globenv)
      END DO
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind_set is not associated",&
                        globenv)
    END IF

  END SUBROUTINE check_atomic_kind_set

! *****************************************************************************

  SUBROUTINE set_atomic_kind(atomic_kind,element_symbol,name,mass,kind_number,&
                             natom,atom_list,&
                             paw_atom,qm_atom,hard_radius,hard0_radius,&
                             elp_potential,soft_basis_set,hard_basis_set,&
                             lmax_rho0)

!   Purpose: Set the components of an atomic kind data set.

!   History: - Creation (11.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: element_symbol, name
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: mass
    INTEGER, INTENT(IN), OPTIONAL            :: kind_number, &
                                                lmax_rho0, natom
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: atom_list
    LOGICAL, INTENT(IN), OPTIONAL            :: paw_atom
    LOGICAL, INTENT(IN), OPTIONAL            :: qm_atom
    REAL(dp), INTENT(IN), OPTIONAL           :: hard_radius
    REAL(dp), INTENT(IN), OPTIONAL           :: hard0_radius
    TYPE(elp_potential_type), OPTIONAL, &
      POINTER                                :: elp_potential
    TYPE(gto_basis_set_type), OPTIONAL, &
      POINTER                                :: soft_basis_set
    TYPE(gto_basis_set_type), OPTIONAL, &
      POINTER                                :: hard_basis_set
    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_atomic_kind"

    INTEGER                                  :: istat, n

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind)) THEN

      IF (PRESENT(element_symbol)) atomic_kind%element_symbol = element_symbol
      IF (PRESENT(name)) atomic_kind%name = name
      IF (PRESENT(mass)) atomic_kind%mass = mass
      IF (PRESENT(kind_number)) atomic_kind%kind_number = kind_number
      IF (PRESENT(natom)) atomic_kind%natom = natom
      IF (PRESENT(atom_list)) THEN
        n = SIZE(atom_list)
        IF (n > 0) THEN
          IF (ASSOCIATED(atomic_kind%atom_list)) THEN
            DEALLOCATE (atomic_kind%atom_list,STAT=istat)
            IF (istat /= 0) THEN
              CALL stop_memory(routine_name,module_name,__LINE__,&
                               "atomic_kind%atom_list")
            END IF
          END IF
          ALLOCATE (atomic_kind%atom_list(n),STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "atomic_kind%atom_list",n*int_size)
          END IF
          atomic_kind%atom_list(:) = atom_list(:)
          atomic_kind%natom = n
        ELSE
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An invalid atom_list was supplied")
        END IF
      END IF
      IF (PRESENT(paw_atom)) atomic_kind%paw_atom = paw_atom
      IF (PRESENT(qm_atom)) atomic_kind%qm_atom = qm_atom
      IF (PRESENT(hard_radius)) atomic_kind%hard_radius = hard_radius
      IF (PRESENT(hard0_radius)) atomic_kind%hard0_radius = hard0_radius
      IF (PRESENT(elp_potential)) atomic_kind%elp_potential => elp_potential
      IF (PRESENT(soft_basis_set)) atomic_kind%soft_basis_set => soft_basis_set
      IF (PRESENT(hard_basis_set)) atomic_kind%hard_basis_set => hard_basis_set
      IF (PRESENT(lmax_rho0))  atomic_kind%lmax_rho0 = lmax_rho0

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind is not associated")

    END IF

  END SUBROUTINE set_atomic_kind

! *****************************************************************************

  SUBROUTINE write_atomic_kind(atomic_kind,globenv)

!   Purpose: Write an atomic kind data set to the output unit.

!   History: - Creation (09.02.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_atomic_kind"

    INTEGER                                  :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(atomic_kind)) THEN

      IF (globenv%print%kind_information) THEN
        output_unit = globenv%scr
        WRITE (UNIT=output_unit,FMT="(/,T2,I2,A,T45,A,T75,I6)")&
          atomic_kind%kind_number,". Atomic kind: "//TRIM(atomic_kind%name),&
          "Number of atoms of this kind: ",atomic_kind%natom
        IF (ANY(atomic_kind%ldaCoreRadius>0.0_dp)) THEN
           WRITE(UNIT=output_unit,fmt="(/,T16,A,T41,4(F10.3))")&
                "Lda core radius = ",atomic_kind%ldaCoreRadius
        END IF
        CALL write_orb_basis_set(atomic_kind%orb_basis_set,globenv)
        IF (atomic_kind%paw_atom) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,A)")&
          "***** GAPW  Atomic Kind: Soft Basis Set and",&
              "Atomic and Compensation Density are Required ***"
          WRITE (UNIT=output_unit,FMT="(/,T16,A,T71,F10.3)")&
          "GAPW|  Hard Gaussians Radius:",atomic_kind%hard_radius
          WRITE (UNIT=output_unit,FMT="(T16,A,T71,F10.3,/)")&
          "GAPW|  Rho0 Radius:",atomic_kind%hard0_radius
          WRITE (UNIT=output_unit,FMT="(T10,A))")&
          "***** Soft Basis Set *****"

          CALL write_orb_basis_set(atomic_kind%soft_basis_set,globenv)
        ENDIF
        CALL write_aux_basis_set(atomic_kind%aux_basis_set,globenv)
        CALL write_potential(atomic_kind%all_potential,globenv)
        CALL write_potential(atomic_kind%elp_potential,globenv)
        CALL write_potential(atomic_kind%gth_potential,globenv)
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind is not associated",globenv)

    END IF

  END SUBROUTINE write_atomic_kind

! *****************************************************************************

  SUBROUTINE write_atomic_kind_set(atomic_kind_set,globenv)

!   Purpose: Write an atomic kind set data set to the output unit.

!   History: - Creation (09.02.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_atomic_kind_set"

    INTEGER                                  :: ikind, nkind, output_unit
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%kind_information) THEN

      IF (ASSOCIATED(atomic_kind_set)) THEN

        output_unit = globenv%scr

        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "ATOMIC KIND INFORMATION"

        nkind = SIZE(atomic_kind_set)

        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL write_atomic_kind(atomic_kind,globenv)
        END DO

      ELSE

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The pointer atomic_kind_set is not associated",&
                          globenv)

      END IF

    END IF

  END SUBROUTINE write_atomic_kind_set

! *****************************************************************************

  SUBROUTINE write_gto_basis_sets(atomic_kind_set,globenv)

!   Purpose: Write all the GTO basis sets of an atomic kind set to the output
!            unit (for the printing of the unnormalized basis sets as read from
!            database).

!   History: - Creation (17.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_gto_basis_sets"

    INTEGER                                  :: ikind, nkind, output_unit
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%basis_set) THEN

      IF (ASSOCIATED(atomic_kind_set)) THEN

        output_unit = globenv%scr

        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
          "BASIS SET INFORMATION (Unnormalised Gaussian-type functions)"

        nkind = SIZE(atomic_kind_set)

        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          WRITE (UNIT=output_unit,FMT="(/,T2,I2,A)")&
            atomic_kind%kind_number,". Atomic kind: "//TRIM(atomic_kind%name)
          CALL write_gto_basis_set(atomic_kind_set(ikind)%orb_basis_set,&
                                   globenv)
          CALL write_gto_basis_set(atomic_kind_set(ikind)%aux_basis_set,&
                                   globenv)
        END DO

      ELSE

        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The pointer atomic_kind_set is not associated",&
                          globenv)

      END IF

    END IF

  END SUBROUTINE write_gto_basis_sets

! *****************************************************************************

END MODULE atomic_kind_types

! *****************************************************************************
