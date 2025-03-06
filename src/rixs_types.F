!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> *************************************************************************************************
!> \brief Define Resonant Inelastic XRAY Scattering (RIXS) control type and associated create, 
!>        release, etc subroutines
!> \author BSG (02.2025)
!> *************************************************************************************************
MODULE rixs_types
    USE cp_fm_types,                           ONLY: cp_fm_type, cp_fm_release 
    USE kinds,                                 ONLY: dp, default_string_length
    USE qs_tddfpt2_types,                      ONLY: tddfpt_ground_state_mos
   ! USE qs_tddfpt2_utils,                      ONLY: tddfpt_release_ground_state_mos
    USE input_section_types,                   ONLY: section_vals_val_get, section_vals_type
    USE cp_dbcsr_api,                          ONLY: dbcsr_release_p, dbcsr_p_type
    USE cp_control_types,                      ONLY: dft_control_type!,&
                                                  !   tddfpt2_control_type
    USE input_constants,                       ONLY: xas_dip_len

! more uses

#include "./base/base_uses.f90"

    IMPLICIT NONE

    PRIVATE

! **************************************************************************************************
!> \brief
! **************************************************************************************************
    TYPE rixs_control_type
        ! types needed for xas_tdp part
        ! TODO: if dipole not lenght, cpabort
        INTEGER                                 :: xas_tdp_dipole_form = 0
        ! TODO: if not tda, cpabort
        LOGICAL                                 :: xas_tdp_tamm_dancoff = .FALSE.
        ! not sure if these are necessary yet
        LOGICAL                                 :: xas_tdp_do_loc = .FALSE.
        LOGICAL                                 :: is_periodic = .FALSE.
        INTEGER                                 :: batch_size = 0

        INTEGER, DIMENSION(:), POINTER          :: xas_tdp_list_ex_atoms => NULL()
        CHARACTER(len=default_string_length), &
            DIMENSION(:), POINTER               :: xas_tdp_list_ex_kinds => NULL()
        INTEGER, DIMENSION(:, :), POINTER       :: state_types => NULL()

        ! DUMMY
        REAL(dp)                                :: def_dum


        ! types needed from tddfpt2 part


    END TYPE rixs_control_type

! **************************************************************************************************
!> \brief
! **************************************************************************************************
    TYPE rixs_env_type
        TYPE(xas_tdp_core_type), POINTER               :: core_state => NULL()
        TYPE(tddfpt2_valence_type), POINTER            :: valence_state => NULL()

        ! TODO: TYPE of properties that will be calculated in rixs section

    END TYPE rixs_env_type

! **************************************************************************************************
!> \brief Core state coming from the xas_tdp routine
! reusing AB's code
! **************************************************************************************************
    TYPE xas_tdp_core_type
        CHARACTER(len=2), DIMENSION(3)          :: state_type_char = ""
        INTEGER                                 :: nex_atoms = 0
        INTEGER                                 :: nex_kinds = 0
        INTEGER, DIMENSION(:), POINTER          :: ex_atom_indices => NULL()
        INTEGER, DIMENSION(:), POINTER          :: ex_kind_indices => NULL()
        INTEGER, DIMENSION(:, :), POINTER       :: state_types => NULL()
        TYPE(rixs_donor_state_type), DIMENSION(:), &
            POINTER                             :: donor_states => NULL()        
        
        TYPE(dbcsr_p_type), DIMENSION(:), &
            POINTER                             :: dipmat => NULL() ! transition dipole matrix

    END TYPE xas_tdp_core_type

! **************************************************************************************************
!> \brief Donor state coming from the xas_tdp routine
! reusing AB's code
! **************************************************************************************************
    TYPE rixs_donor_state_type
        INTEGER                                 :: at_index = 0
        INTEGER                                 :: kind_index = 0
        INTEGER                                 :: ndo_mo = 0
        CHARACTER(LEN=default_string_length)    :: at_symbol = ""
        INTEGER                                 :: state_type = 0

        ! eigenvectors
        TYPE(cp_fm_type), POINTER               :: evects => NULL() ! *_coeffs
        ! TYPE(cp_fm_type), POINTER               :: sc_coeffs => NULL() ! spin-conserving
        ! TYPE(cp_fm_type), POINTER               :: sf_coeffs => NULL() ! spin-flip
        ! TYPE(cp_fm_type), POINTER               :: sg_coeffs => NULL() ! singlet
        ! TYPE(cp_fm_type), POINTER               :: tp_coeffs => NULL() ! triplet

        TYPE(cp_fm_type), POINTER               :: gs_coeffs => NULL()  ! gs mo coeffs
        INTEGER, DIMENSION(:, :), POINTER       :: mo_indices => NULL() 

        ! solutions to lr-tddfpt
        REAL(dp), DIMENSION(:), POINTER         :: evals => NULL() ! *_evals
        ! REAL(dp), DIMENSION(:), POINTER         :: sc_evals => NULL() 
        ! REAL(dp), DIMENSION(:), POINTER         :: sf_evals => NULL()
        ! REAL(dp), DIMENSION(:), POINTER         :: sg_evals => NULL()
        ! REAL(dp), DIMENSION(:), POINTER         :: tp_evals => NULL()
        ! REAL(dp), DIMENSION(:), POINTER         :: soc_evals => NULL()

        ! maybe energy_evals not necessary 
        REAL(dp), DIMENSION(:, :), POINTER      :: energy_evals => NULL() ! energy eigenvalue
        REAL(dp), DIMENSION(:, :), POINTER      :: osc_str => NULL()

    END TYPE rixs_donor_state_type


! **************************************************************************************************
!> \brief Valence state coming from the qs_tddfpt2 routine
! **************************************************************************************************
    TYPE tddfpt2_valence_type
    
        TYPE(cp_fm_type), POINTER                        :: evects => NULL() ! eigenvectors
        REAL(dp), DIMENSION(:), POINTER                  :: evals => NULL()  ! energies

        ! TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
        !    POINTER                                       :: gs_mos

        ! stuff below coming from tddfpt_ground_state_mos type
        TYPE(cp_fm_type), POINTER                        :: mos_occ => NULL()
        TYPE(cp_fm_type), POINTER                        :: mos_virt => NULL()
        TYPE(cp_fm_type), POINTER                        :: evals_occ_matrix => NULL()
        REAL(dp), DIMENSION(:), ALLOCATABLE              :: evals_occ
        REAL(dp), DIMENSION(:), ALLOCATABLE              :: evals_virt

    END TYPE tddfpt2_valence_type

    CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rixs_types'

    PUBLIC :: rixs_control_type, rixs_env_type, xas_tdp_core_type, tddfpt2_valence_type, &
              rixs_donor_state_type

    PUBLIC :: rixs_control_create, rixs_control_release, rixs_env_create, rixs_env_release, &
              read_rixs_control

CONTAINS 

! **************************************************************************************************
!> \brief Creates and initializes the rixs_control_type
!> \param rixs_control the type to initialize
! **************************************************************************************************
   SUBROUTINE rixs_control_create(rixs_control)
        TYPE(rixs_control_type), POINTER :: rixs_control

        CPASSERT(.NOT. ASSOCIATED(rixs_control))
        ALLOCATE(rixs_control)

        rixs_control%xas_tdp_dipole_form = xas_dip_len ! different from xas_tdp def value
        rixs_control%xas_tdp_tamm_dancoff = .FALSE. ! maybe not necessary, could be checked during the rixs call
        rixs_control%xas_tdp_do_loc = .FALSE. !
        rixs_control%is_periodic = .FALSE.
        rixs_control%batch_size = 64
        NULLIFY(rixs_control%xas_tdp_list_ex_atoms)
        NULLIFY(rixs_control%xas_tdp_list_ex_kinds)
        NULLIFY(rixs_control%state_types)
        
   END SUBROUTINE rixs_control_create

! **************************************************************************************************
!> \brief Releases the rixs_control_type
! **************************************************************************************************
   SUBROUTINE rixs_control_release(rixs_control)
        TYPE(rixs_control_type), POINTER :: rixs_control

        IF (ASSOCIATED(rixs_control)) THEN
            IF (ASSOCIATED(rixs_control%xas_tdp_list_ex_atoms)) THEN
                DEALLOCATE(rixs_control%xas_tdp_list_ex_atoms)
            END IF
            IF (ASSOCIATED(rixs_control%xas_tdp_list_ex_kinds)) THEN
                DEALLOCATE(rixs_control%xas_tdp_list_ex_kinds)
            END IF
            IF (ASSOCIATED(rixs_control%state_types)) THEN
                DEALLOCATE(rixs_control%state_types)
            END IF
        END IF

   END SUBROUTINE rixs_control_release

! **************************************************************************************************
!> \brief Reads the input and stores in the rixs_control_type
! **************************************************************************************************
    SUBROUTINE read_rixs_control(rixs_control, rixs_section)
        TYPE(rixs_control_type), POINTER             :: rixs_control
        TYPE(section_vals_type), POINTER             :: rixs_section

        REAL(dp)                                     :: def_dum

        LOGICAL                                      :: do_rixs

        CALL section_vals_val_get(rixs_section, "_SECTION_PARAMETERS_", l_val=do_rixs)
        CALL section_vals_val_get(rixs_section, "DUMMY", r_val=def_dum)

        rixs_control%def_dum = def_dum

        ! read input, empty for now
        ! maybe could read xas_tdp input?

    END SUBROUTINE read_rixs_control


! **************************************************************************************************
!> \brief Creates a rixs environment type
!> \param rixs_env the type to create
! **************************************************************************************************
    SUBROUTINE rixs_env_create(rixs_env)
        TYPE(rixs_env_type), POINTER                        :: rixs_env

        CPASSERT(.NOT. ASSOCIATED(rixs_env))

        ALLOCATE(rixs_env)

        NULLIFY(rixs_env%core_state)
        NULLIFY(rixs_env%valence_state)

    END SUBROUTINE rixs_env_create

! **************************************************************************************************
!> \brief Releases the rixs environment type
!> \param rixs_env the type to release
! **************************************************************************************************
    SUBROUTINE rixs_env_release(rixs_env)
        TYPE(rixs_env_type), POINTER                       :: rixs_env

        IF (ASSOCIATED(rixs_env)) THEN
            IF (ASSOCIATED(rixs_env%core_state)) THEN
                CALL xas_tdp_core_release(rixs_env%core_state)
            END IF
            IF (ASSOCIATED(rixs_env%valence_state)) THEN
                CALL tddfpt2_valence_state_release(rixs_env%valence_state)
            END IF
        END IF

    END SUBROUTINE rixs_env_release

! **************************************************************************************************
!> \brief Creates the core state type
!> \param core_state the type to create
! **************************************************************************************************
    SUBROUTINE xas_tdp_core_create(core_state)
        TYPE(xas_tdp_core_type), POINTER               :: core_state

        core_state%nex_atoms = 1
        core_state%nex_kinds = 1

        NULLIFY(core_state%ex_atom_indices)
        NULLIFY(core_state%ex_kind_indices)
        NULLIFY(core_state%state_types)
        NULLIFY(core_state%donor_states)
        NULLIFY(core_state%dipmat)

    END SUBROUTINE xas_tdp_core_create

! **************************************************************************************************
!> \brief Releases the core state type
!> \param core_state the type to release
! **************************************************************************************************
    SUBROUTINE xas_tdp_core_release(core_state)
        TYPE(xas_tdp_core_type), POINTER               :: core_state

        INTEGER                                        :: i

        IF (ASSOCIATED(core_state)) THEN
            IF (ASSOCIATED(core_state%ex_atom_indices)) THEN
                DEALLOCATE(core_state%ex_atom_indices)
            END IF
            IF (ASSOCIATED(core_state%ex_kind_indices)) THEN
                DEALLOCATE(core_state%ex_kind_indices)
            END IF
            IF (ASSOCIATED(core_state%state_types)) THEN
                DEALLOCATE(core_state%state_types)
            END IF
            IF (ASSOCIATED(core_state%donor_states)) THEN
                CALL rixs_donor_state_release(core_state%donor_states)
            END IF
            IF (ASSOCIATED(core_state%dipmat)) THEN
                DO i = 1, SIZE(core_state%dipmat)
                    CALL dbcsr_release_p(core_state%dipmat(i)%matrix)
                END DO
                DEALLOCATE(core_state%dipmat)
            END IF
        END IF

    END SUBROUTINE xas_tdp_core_release

! **************************************************************************************************
!> \brief Creates the rixs donor state type
!> \param rixs_env the type to create
! **************************************************************************************************
    SUBROUTINE rixs_donor_state_create(rixs_donor_state)
        TYPE(rixs_donor_state_type), POINTER               :: rixs_donor_state

        NULLIFY(rixs_donor_state%evects)
        NULLIFY(rixs_donor_state%gs_coeffs)
        NULLIFY(rixs_donor_state%mo_indices)
        NULLIFY(rixs_donor_state%evals)

        ! maybe needed
        NULLIFY(rixs_donor_state%energy_evals)
        NULLIFY(rixs_donor_state%osc_str)

    END SUBROUTINE rixs_donor_state_create

! **************************************************************************************************
!> \brief Releases the rixs donor state type
!> \param rixs_env the type to release
! **************************************************************************************************
    SUBROUTINE rixs_donor_state_release(donor_state_set)
        TYPE(rixs_donor_state_type), DIMENSION(:), &
            POINTER                                        :: donor_state_set
        INTEGER                                            :: i

        IF (ASSOCIATED(donor_state_set)) THEN
            DO i = 1, SIZE(donor_state_set)
                IF (ASSOCIATED(donor_state_set(i)%evects)) THEN
                    CALL cp_fm_release(donor_state_set(i)%evects)
                    DEALLOCATE(donor_state_set(i)%evects)
                END IF
                IF (ASSOCIATED(donor_state_set(i)%gs_coeffs)) THEN
                    CALL cp_fm_release(donor_state_set(i)%gs_coeffs)
                    DEALLOCATE(donor_state_set(i)%gs_coeffs)
                END IF
                IF (ASSOCIATED(donor_state_set(i)%evals)) THEN
                    DEALLOCATE(donor_state_set(i)%evals)
                END IF
                IF (ASSOCIATED(donor_state_set(i)%mo_indices)) THEN
                    DEALLOCATE(donor_state_set(i)%mo_indices)
                END IF

                ! maybe not necessary
                IF (ASSOCIATED(donor_state_set(i)%osc_str)) THEN
                    DEALLOCATE(donor_state_set(i)%osc_str)
                END IF
                IF (ASSOCIATED(donor_state_set(i)%energy_evals)) THEN
                    DEALLOCATE(donor_state_set(i)%energy_evals)
                END IF
            END DO
        END IF

    END SUBROUTINE rixs_donor_state_release

! **************************************************************************************************
!> \brief Creates the valence state type
! **************************************************************************************************
    SUBROUTINE tddfpt2_valence_state_create(valence_state)
        TYPE(tddfpt2_valence_type), POINTER               :: valence_state

        NULLIFY(valence_state%evects)
        NULLIFY(valence_state%evals)
        ! NULLIFY(valence_state%gs_mos)
        
        ! stuff below coming from tddfpt_ground_state_mos type
        NULLIFY(valence_state%mos_occ)
        NULLIFY(valence_state%mos_virt)
        NULLIFY(valence_state%evals_occ_matrix)
        IF (ALLOCATED(valence_state%evals_occ)) THEN
            DEALLOCATE(valence_state%evals_occ)
        END IF
        IF (ALLOCATED(valence_state%evals_virt)) THEN
            DEALLOCATE(valence_state%evals_virt)
        END IF

    END SUBROUTINE tddfpt2_valence_state_create

! **************************************************************************************************
!> \brief Releases the valence state type
! **************************************************************************************************
    SUBROUTINE tddfpt2_valence_state_release(valence_state)
        TYPE(tddfpt2_valence_type), POINTER               :: valence_state

        IF (ASSOCIATED(valence_state)) THEN
            IF (ASSOCIATED(valence_state%evects)) THEN
                CALL cp_fm_release(valence_state%evects)
                DEALLOCATE(valence_state%evects)
            END IF
            IF (ASSOCIATED(valence_state%evals)) THEN
                DEALLOCATE(valence_state%evals)
            END IF
            ! IF (ASSOCIATED(valence_state%gs_mos)) THEN
            !     DO i = 1, SIZE(valence_state%gs_mos)
            !         CALL tddfpt_release_ground_state_mos(valence_state%gs_mos(i))
            !     END DO
            !     DEALLOCATE(valence_state%gs_mos)
            ! END IF
            IF (ASSOCIATED(valence_state%mos_occ)) THEN
                DEALLOCATE(valence_state%mos_occ)
            END IF
            IF (ASSOCIATED(valence_state%mos_virt)) THEN
                DEALLOCATE(valence_state%mos_virt)
            END IF
            IF (ASSOCIATED(valence_state%evals_occ_matrix)) THEN
                DEALLOCATE(valence_state%evals_occ_matrix)
            END IF
            IF (ALLOCATED(valence_state%evals_occ)) THEN
                DEALLOCATE(valence_state%evals_occ)
            END IF
            IF (ALLOCATED(valence_state%evals_virt)) THEN
                DEALLOCATE(valence_state%evals_virt)
            END IF
        END IF

    END SUBROUTINE tddfpt2_valence_state_release


END MODULE rixs_types