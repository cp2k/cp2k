!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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
                                             write_aux_basis_set,&
                                             write_gto_basis_set,&
                                             write_orb_basis_set
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_linked_list_int,              ONLY: cp_sll_int_create,&
                                             cp_sll_int_dealloc,&
                                             cp_sll_int_get_empty,&
                                             cp_sll_int_insert_ordered,&
                                             cp_sll_int_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
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
                                             write_potential
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             dp
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
  USE pao_bas_projection_types,        ONLY: cp_proj_dealloc,&
                                             cp_proj_create2,&
                                             cp_proj_get,&
                                             cp_proj_write,&
                                             pao_bas_projection_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: massunit
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE string_utilities,                ONLY: integer_to_string,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "atomic_kind_types"
  INTEGER, PARAMETER          :: max_name_length = 60

! *** Define the atomic kind types ***

  TYPE atomic_kind_type
    PRIVATE
    TYPE(all_potential_type), POINTER      :: all_potential
    TYPE(elp_potential_type), POINTER      :: elp_potential
    TYPE(gth_potential_type), POINTER      :: gth_potential
    TYPE(gto_basis_set_type), POINTER      :: aux_basis_set
    TYPE(gto_basis_set_type), POINTER      :: orb_basis_set
    TYPE(pao_bas_projection_type), POINTER :: pao_bas_proj
    CHARACTER(LEN=max_name_length)         :: name
    CHARACTER(LEN=2)                       :: element_symbol
    REAL(KIND = dp)                               :: mass
    INTEGER                                :: kind_number,&
                                              natom,&
                                              number_of_grid_points
    INTEGER, DIMENSION(:), POINTER         :: atom_list
  END TYPE atomic_kind_type

! *** Public subroutines ***

  PUBLIC :: allocate_atomic_kind_set,&
            deallocate_atomic_kind_set,&
            get_atomic_kind,&
            get_atomic_kind_set,&
            init_atomic_kind_set,&
            read_atomic_kind_set,&
            set_atomic_kind,&
            write_atomic_kind,&
            write_atomic_kind_set,&
            write_gto_basis_sets

! *** Public data types ***

  PUBLIC :: atomic_kind_type

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_atomic_kind_set(atomic_kind_set,nkind)

!   Purpose: Allocate an atomic kind set.

!   History: - Creation (08.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    INTEGER, INTENT(IN)                           :: nkind

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_atomic_kind_set"

!   *** Local variables ***

    INTEGER :: ikind,istat

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
      NULLIFY (atomic_kind_set(ikind)%aux_basis_set)
      NULLIFY (atomic_kind_set(ikind)%orb_basis_set)
      NULLIFY (atomic_kind_set(ikind)%atom_list)
      NULLIFY (atomic_kind_set(ikind)%pao_bas_proj)
      atomic_kind_set(ikind)%name = ""
      atomic_kind_set(ikind)%element_symbol = ""
      atomic_kind_set(ikind)%mass = 0.0_dp
      atomic_kind_set(ikind)%kind_number = 0
      atomic_kind_set(ikind)%natom = 0
      atomic_kind_set(ikind)%number_of_grid_points = 0
    END DO

  END SUBROUTINE allocate_atomic_kind_set

! *****************************************************************************

  SUBROUTINE deallocate_atomic_kind_set(atomic_kind_set)

!   Purpose: Deallocate an atomic kind set.

!   History: - Creation (09.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_atomic_kind_set"

!   *** Local variables ***

    INTEGER :: ikind,istat,nkind

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
        IF (ASSOCIATED(atomic_kind_set(ikind)%aux_basis_set)) THEN
          CALL deallocate_gto_basis_set(atomic_kind_set(ikind)%aux_basis_set)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%orb_basis_set)) THEN
          CALL deallocate_gto_basis_set(atomic_kind_set(ikind)%orb_basis_set)
        END IF
        IF (ASSOCIATED(atomic_kind_set(ikind)%pao_bas_proj)) THEN
          CALL cp_proj_dealloc(atomic_kind_set(ikind)%pao_bas_proj)
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
                             gth_potential,aux_basis_set,orb_basis_set,&
                             element_symbol,name,mass,kind_number,natom,&
                             atom_list,alpha_core_charge,ccore_charge,&
                             core_charge_radius,ncgf,nsgf,z,zeff,qeff,&
                             pao_bas_proj)

!   Purpose: Get informations about an atomic kind.

!   History: - Creation (10.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER                  :: atomic_kind
    TYPE(all_potential_type), OPTIONAL, POINTER      :: all_potential
    TYPE(elp_potential_type), OPTIONAL, POINTER      :: elp_potential
    TYPE(gth_potential_type), OPTIONAL, POINTER      :: gth_potential
    TYPE(gto_basis_set_type), OPTIONAL, POINTER      :: aux_basis_set,&
                                                        orb_basis_set
    CHARACTER(LEN=2), OPTIONAL, INTENT(OUT)          :: element_symbol
    CHARACTER(LEN=max_name_length), OPTIONAL,&
      INTENT(OUT)                                    :: name
    REAL(KIND = dp), OPTIONAL, INTENT(OUT)                  :: alpha_core_charge,&
                                                        ccore_charge,&
                                                        core_charge_radius,&
                                                        mass,qeff
    INTEGER, OPTIONAL, INTENT(OUT)                   :: kind_number,natom,&
                                                        ncgf,nsgf,z,zeff
    INTEGER, DIMENSION(:), OPTIONAL, POINTER         :: atom_list
    TYPE(pao_bas_projection_type), OPTIONAL, POINTER :: pao_bas_proj

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_atomic_kind"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind)) THEN

      IF (PRESENT(all_potential)) all_potential => atomic_kind%all_potential
      IF (PRESENT(elp_potential)) elp_potential => atomic_kind%elp_potential
      IF (PRESENT(gth_potential)) gth_potential => atomic_kind%gth_potential
      IF (PRESENT(aux_basis_set)) aux_basis_set => atomic_kind%aux_basis_set
      IF (PRESENT(orb_basis_set)) orb_basis_set => atomic_kind%orb_basis_set
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
      IF (PRESENT(pao_bas_proj)) pao_bas_proj => atomic_kind%pao_bas_proj

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind is not associated")

    END IF

  END SUBROUTINE get_atomic_kind

! *****************************************************************************

  SUBROUTINE get_atomic_kind_set(atomic_kind_set,maxatom,maxcgf,maxco,maxlgto,&
                                 maxlppl,maxlppnl,maxpgf,maxppnl,maxsgf,&
                                 maxsgf_set,maxshell,pao_max_min_bas,&
                                 pao_max_full_bas,pao_max_excl_bas,&
                                 pao_tot_min_bas,pao_tot_full_bas,&
                                 pao_tot_excl_bas,pao_max_l_angles,&
                                 natom,nelectron,ncgf,ncgf_aux,&
                                 npgf,nset,nsgf,nshell,&
                                 all_potential_present,gth_potential_present,&
                                 elp_potential_present,&
                                 atom_of_kind,kind_of,natom_of_kind,maxder)

!   Purpose: Get informations about an atomic kind set.

!   History: - Creation (14.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    INTEGER, OPTIONAL, INTENT(OUT)                :: maxatom,maxcgf,maxco,&
                                                     maxlgto,maxlppl,maxlppnl,&
                                                     maxpgf,maxppnl,maxsgf,&
                                                     maxsgf_set,maxshell,&
                                                     pao_max_min_bas,&
                                                     pao_max_full_bas,&
                                                     pao_max_excl_bas,&
                                                     pao_tot_min_bas,&
                                                     pao_tot_full_bas,&
                                                     pao_tot_excl_bas,&
                                                     pao_max_l_angles,&
                                                     natom,ncgf,ncgf_aux,&
                                                     nelectron,&
                                                     npgf,nset,nsgf,nshell
    LOGICAL, OPTIONAL, INTENT(OUT)                :: all_potential_present,&
                                                     elp_potential_present,&
                                                     gth_potential_present
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT)  :: atom_of_kind,kind_of,&
                                                     natom_of_kind
    INTEGER, OPTIONAL, INTENT(IN)                 :: maxder

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_atomic_kind_set"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(all_potential_type), POINTER :: all_potential
    TYPE(elp_potential_type), POINTER :: elp_potential
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: aux_basis_set,orb_basis_set
    INTEGER                           :: atom_a,iatom,ikind,imax,n,nkind,zeff,&
                                         p_min_bas,p_min_bas_max,p_excl_bas,&
                                         p_excl_bas_max,p_full_bas,&
                                         p_full_bas_max,p_min_bas_tot,&
                                         p_full_bas_tot,p_l_angles_max
    LOGICAL                           :: pao

    REAL(KIND = dp), DIMENSION(:,:), POINTER :: zet

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN

      IF (PRESENT(maxatom)) maxatom = 0
      IF (PRESENT(maxcgf)) maxcgf = 0
      IF (PRESENT(maxco)) maxco = 0
      IF (PRESENT(maxlgto)) maxlgto = -1
      IF (PRESENT(maxlppl)) maxlppl = -1
      IF (PRESENT(maxlppnl)) maxlppnl = -1
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
      IF (PRESENT(atom_of_kind)) atom_of_kind(:) = 0
      IF (PRESENT(kind_of)) kind_of(:) = 0
      IF (PRESENT(natom_of_kind)) natom_of_kind(:) = 0
      pao=.FALSE.
      IF (PRESENT(pao_max_excl_bas) .OR. PRESENT(pao_max_full_bas) .OR.&
           PRESENT(pao_max_min_bas).OR. PRESENT(pao_tot_min_bas).OR.&
           PRESENT(pao_tot_full_bas) .OR. PRESENT(pao_tot_excl_bas).OR.&
           PRESENT(pao_max_l_angles)) THEN
         IF (SIZE(atomic_kind_set)>0) THEN
            IF (ASSOCIATED(atomic_kind_set(1)%pao_bas_proj)) pao=.TRUE.
         END IF
         IF (.NOT. pao) THEN
            CALL cp_error_message(cp_warning_level,routine_name,&
                 "requested pao_bas_proj info, and pao_bas_proj is not allocated")
         END IF
         p_min_bas_max=0
         p_excl_bas_max=0
         p_full_bas_max=0
         p_min_bas_tot=0
         p_full_bas_tot=0
         p_l_angles_max=0
      END IF
      nkind = SIZE(atomic_kind_set)
      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             all_potential=all_potential,&
                             elp_potential=elp_potential,&
                             gth_potential=gth_potential,&
                             aux_basis_set=aux_basis_set,&
                             orb_basis_set=orb_basis_set)
        IF (PRESENT(maxatom)) THEN
          maxatom = MAX(maxatom,atomic_kind%natom)
        END IF
        IF (PRESENT(maxlppl).AND.ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,nexp_ppl=n)
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
        IF (pao) THEN
           CALL cp_proj_get(atomic_kind%pao_bas_proj,&
                size_min_bas=p_min_bas,&
                size_excl_bas=p_excl_bas, size_full_bas=p_full_bas)
           p_min_bas_max=MAX(p_min_bas_max,p_min_bas)
           p_excl_bas_max=MAX(p_excl_bas_max,p_excl_bas)
           p_full_bas_max=MAX(p_full_bas_max,p_full_bas)
           p_l_angles_max=MAX(p_l_angles_max,p_min_bas*p_excl_bas)
           p_min_bas_tot=p_min_bas_tot+p_min_bas*SIZE(atomic_kind%atom_list)
           p_full_bas_tot=p_full_bas_tot+p_full_bas*SIZE(atomic_kind%atom_list)
        END IF
      END DO
      IF (PRESENT(pao_max_excl_bas)) pao_max_excl_bas=p_excl_bas_max
      IF (PRESENT(pao_max_full_bas)) pao_max_full_bas=p_full_bas_max
      IF (PRESENT(pao_max_min_bas)) pao_max_min_bas=p_min_bas_max
      IF (PRESENT(pao_tot_excl_bas)) &
           pao_tot_excl_bas=p_full_bas_tot-p_min_bas_tot
      IF (PRESENT(pao_tot_full_bas)) pao_tot_full_bas=p_full_bas_tot
      IF (PRESENT(pao_tot_min_bas)) pao_tot_min_bas=p_min_bas_tot
      IF (PRESENT(pao_max_l_angles)) pao_max_l_angles=p_l_angles_max

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind_set is not associated")

    END IF

  END SUBROUTINE get_atomic_kind_set

! *****************************************************************************

  SUBROUTINE init_atomic_kind(atomic_kind,globenv)

!   Purpose: Initialise an atomic kind data set.

!   History: - Creation (11.01.2002,MK)
!              20.09.2002 adapted for pol/kg use, gtb

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_atomic_kind"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind)) THEN

      IF (globenv%program_name == "QS") THEN
        CALL init_potential(atomic_kind%gth_potential)
        CALL init_orb_basis_set(atomic_kind%orb_basis_set)
      ELSE IF (globenv%program_name == "POL") THEN
        CALL init_orb_basis_set(atomic_kind%orb_basis_set,globenv)
        CALL init_aux_basis_set(atomic_kind%aux_basis_set)
      ELSE IF (globenv%program_name == "KG") THEN
        CALL init_orb_basis_set(atomic_kind%orb_basis_set,globenv)
        IF (ASSOCIATED(atomic_kind%aux_basis_set))THEN
          CALL init_aux_basis_set(atomic_kind%aux_basis_set)
        END IF
      ELSE IF (globenv%program_name == "FIST") THEN
        CALL allocate_potential(atomic_kind%elp_potential)
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind is not associated",globenv)

    END IF

  END SUBROUTINE init_atomic_kind

! *****************************************************************************

  SUBROUTINE init_atomic_kind_set(atomic_kind_set,globenv)

!   Purpose: Initialise an atomic kind set data set.

!   History: - Creation (17.01.2002,MK)
!            - 20.09.2002 globenv passed (gt)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_atomic_kind_set"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER                         :: ikind,nkind

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(atomic_kind_set)) THEN

      nkind = SIZE(atomic_kind_set)

      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL init_atomic_kind(atomic_kind,globenv)
      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer atomic_kind_set is not associated",&
                        globenv)

    END IF

  END SUBROUTINE init_atomic_kind_set

! *****************************************************************************

  SUBROUTINE read_atomic_kind(atomic_kind,globenv)

!   Purpose: Read an atomic kind data set from the input file.

!   History: - Creation (09.02.2002,MK)
!            - 20.09.2002,gt: adapted for POL/KG use (elp_potential)

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER           :: atomic_kind
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_atomic_kind"

!   *** Local variables ***

    TYPE(cp_error_type)            :: error
    TYPE(cp_sll_int_type), POINTER :: min_bas_sll
    CHARACTER(LEN=max_name_length) :: aux_basis_set_name,keyword,&
                                      orb_basis_set_name,potential_name,&
                                      potential_type
    CHARACTER(LEN=3)               :: test_result
    CHARACTER(LEN=2)               :: current_symbol,element_symbol
    INTEGER                        :: i,ikind,ipos,istat,z, min_b_index
    LOGICAL                        :: found,section_found, failure

!   ---------------------------------------------------------------------------

    failure = .FALSE.
    aux_basis_set_name = ""
    orb_basis_set_name = ""
    potential_name = "GTH"
    potential_type = ""
    z = -1
    CALL cp_sll_int_create(min_bas_sll)
    CALL cp_error_init(error)

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="KIND",&
                      start_section_name=atomic_kind%name,&
                      section_found=section_found)
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
          CASE ("ORBITAL_BASIS_SET","ORB_BASIS_SET","ORB_BASIS","BASIS_SET","BAS")
            CALL read_object(orb_basis_set_name)
          CASE ("ELEMENT_SYMBOL","ELEMENT")
            CALL read_object(element_symbol,lower_to_upper=.TRUE.)
            DO i=LBOUND(ptable,1),UBOUND(ptable,1)
              current_symbol = ptable(i)%symbol
              CALL uppercase(current_symbol)
              IF (current_symbol == element_symbol) THEN
                atomic_kind%element_symbol = ptable(i)%symbol
                atomic_kind%mass = ptable(i)%amass
                z = i
                EXIT
              END IF
            END DO
          CASE ("ATOMIC_MASS","ATOMIC_WEIGHT","MASS","WEIGHT")
            CALL read_object(atomic_kind%mass)
          CASE ("POTENTIAL","POT")
            CALL read_object(potential_name)
          CASE ("PAO_MIN_BAS")
            DO
              test_result = test_object()
              IF (test_result == "INT") THEN
                CALL read_object(min_b_index)
                CALL cp_sll_int_insert_ordered(min_bas_sll,min_b_index,&
                                               error=error)
              ELSE
                EXIT
              END IF
            END DO
          CASE DEFAULT
            CALL stop_parser(routine_name,module_name,__LINE__,"INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    ELSE
      CALL stop_parser(routine_name,module_name,__LINE__,"SECTION_NOT_FOUND")
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
          atomic_kind%mass = ptable(i)%amass
          z = i
          EXIT
        END IF
      END DO
    END IF

    atomic_kind%mass = ptable(i)%amass*massunit

    IF (z == -1) THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "No element could be assign to the atomic kind <"//&
                        TRIM(atomic_kind%name)//">. Specify a valid "//&
                        "element symbol in the input section <KIND "//&
                        TRIM(atomic_kind%name)//"> or select a valid "//&
                        "element symbol for the atomic kind name",globenv)
    END IF

    ! this should be done in e.g. QS itself.  
    ! this routine should just read in what it finds in the input 
    IF (globenv%program_name == "QS") THEN 
      IF (orb_basis_set_name == '') THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "No orbital basis set was defined for the atomic "//&
                          "kind <"//TRIM(atomic_kind%name)//">. Specify a "//&
                          "valid orbital basis set name in the input "//&
                          "section <KIND "//TRIM(atomic_kind%name)//">",globenv)
      END IF
      IF (potential_name == '') THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "No atomic potential type was defined for the "//&
                          "atomic kind <"//TRIM(atomic_kind%name)//">. "//&
                          "Specify a valid potential type in the input "//&
                          "section <KIND "//TRIM(atomic_kind%name)//">",globenv)
      END IF
    END IF

!   *** Allocate and initialise the orbital basis set data set structure ***
    CALL init_orbital_pointers(5) ! debug the SUN optimizer

    IF (orb_basis_set_name /= '') THEN
      CALL allocate_gto_basis_set(atomic_kind%orb_basis_set)
      CALL read_gto_basis_set(atomic_kind%element_symbol,orb_basis_set_name,&
                              atomic_kind%orb_basis_set,globenv)
    END IF

!   *** Allocate and initialise the auxiliary basis set data set structure ***

    IF (aux_basis_set_name /= '') THEN
      CALL allocate_gto_basis_set(atomic_kind%aux_basis_set)
      CALL read_gto_basis_set(atomic_kind%element_symbol,aux_basis_set_name,&
                              atomic_kind%aux_basis_set,globenv)
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
      CASE ("GTH")
        CALL allocate_potential(atomic_kind%gth_potential)
        CALL read_potential(atomic_kind%element_symbol,potential_name,&
                            atomic_kind%gth_potential,globenv)
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

    END IF

    IF (.not.cp_sll_int_get_empty(min_bas_sll)) THEN
        CALL cp_proj_create2(atomic_kind%pao_bas_proj,&
                proj_indexes=min_bas_sll,&
                size_full_bas=atomic_kind%orb_basis_set%nsgf,&
                error=error)
     ELSE
        NULLIFY(atomic_kind%pao_bas_proj)
    END IF
    CALL cp_sll_int_dealloc(min_bas_sll,error=error)
    CALL cp_error_dealloc_ref(error)

  END SUBROUTINE read_atomic_kind

! *****************************************************************************

  SUBROUTINE read_atomic_kind_set(atomic_kind_set,globenv)

!   Purpose: Read an atomic kind set data set from the input file.

!   History: - Creation (09.02.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_atomic_kind_set"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER                         :: ikind,nkind

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

  SUBROUTINE set_atomic_kind(atomic_kind,element_symbol,name,mass,kind_number,&
                             natom,atom_list,pao_bas_proj)

!   Purpose: Set the components of an atomic kind data set.

!   History: - Creation (11.01.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), POINTER             :: atomic_kind
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)      :: element_symbol
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)      :: name
    REAL(KIND = dp), OPTIONAL, INTENT(IN)              :: mass
    INTEGER, OPTIONAL, INTENT(IN)               :: kind_number
    INTEGER, OPTIONAL, INTENT(IN)               :: natom
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: atom_list
    TYPE(pao_bas_projection_type), OPTIONAL,&
      INTENT(IN), TARGET                        :: pao_bas_proj

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_atomic_kind"

!   *** Local variables ***

    INTEGER :: istat,n

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
      IF (PRESENT(pao_bas_proj)) atomic_kind%pao_bas_proj => pao_bas_proj

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

    TYPE(atomic_kind_type), POINTER           :: atomic_kind
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_atomic_kind"

!   *** Local variables ***

    INTEGER :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (ASSOCIATED(atomic_kind)) THEN

      IF (globenv%print%kind_information) THEN
        output_unit = globenv%scr
        WRITE (UNIT=output_unit,FMT="(/,T2,I2,A,/,/,T6,A,I6)")&
          atomic_kind%kind_number,". Atomic kind: "//TRIM(atomic_kind%name),&
          "Number of atoms of this atomic kind in the cell: ",atomic_kind%natom
        CALL write_orb_basis_set(atomic_kind%orb_basis_set,globenv)
        CALL write_aux_basis_set(atomic_kind%aux_basis_set,globenv)
        CALL write_potential(atomic_kind%all_potential,globenv)
        CALL write_potential(atomic_kind%elp_potential,globenv)
        CALL write_potential(atomic_kind%gth_potential,globenv)
        IF (ASSOCIATED(atomic_kind%pao_bas_proj).AND.globenv%ionode) THEN
          CALL cp_proj_write(atomic_kind%pao_bas_proj, unit_nr=globenv%scr)
        END IF
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

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_atomic_kind_set"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER                         :: ikind,nkind,output_unit

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

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_gto_basis_sets"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    INTEGER                         :: ikind,nkind,output_unit

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
