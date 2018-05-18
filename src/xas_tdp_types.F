!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!> *************************************************************************************************
!> \brief Define XAS TDP control type and associated create, release, etc subroutines, as well as
!>        XAS TDP environment type and associated set, get, etc subroutines
!> \author AB (11.2017)
!> *************************************************************************************************
MODULE xas_tdp_types
   
   USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_release, cp_fm_p_type
   USE dbcsr_api,                       ONLY: dbcsr_type, dbcsr_release, dbcsr_p_type, &
                                              dbcsr_deallocate_matrix_set
   USE input_constants,                 ONLY: xas_dip_vel,&
                                              xas_tdp_by_index,&
                                              xas_tdp_by_kind, &
                                              tddfpt_singlet
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_get_subs_vals
   USE kinds,                           ONLY: dp ,default_string_length
   USE memory_utilities,                ONLY: reallocate
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type,&
                                              qs_loc_env_release
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief Type containing control information for TDP XAS calculations
!> \param define_excited whether excited atoms are chosen by kind or index
!> \param dipole_form whether the dipole moment is computed in the length or velocity representation
!> \param n_search # of lowest energy MOs to search for donor orbitals
!> \param n_gauss # of GTOs used in STO expansion for core donor state identification
!> \param check_donor_states whether one only checks MO to donor state assignment
!> \param tamm_dancoff whether the calculations should be done in the Tamm-Dancoff approximation
!> \param full_diag whether the LR matrix should be fully built and diagonalized
!> \param list_ex_atoms list of excited atom indices, kept empty if define_excited=by_kind
!> \param list_ex_kinds list of excited atom kinds, kept empty if define_excited=by_index
!> \param do_singlet whether the excitations to consider are singlet (versus triplet)
!> \param state_types columns correspond to the states to excite for each atom kind/index
!>                    the number of rows is the number of times the keyword is repeated
!> \param eri_mme_subsection the subsection dealing with the param for the eri_mme for the kernels
! **************************************************************************************************
   TYPE xas_tdp_control_type
      INTEGER                                 :: define_excited   
      INTEGER                                 :: dipole_form
      INTEGER                                 :: n_search
      INTEGER                                 :: n_gauss
      LOGICAL                                 :: check_donor_states
      LOGICAL                                 :: tamm_dancoff
      LOGICAL                                 :: full_diag
      LOGICAL                                 :: do_singlet
      INTEGER, DIMENSION(:), POINTER          :: list_ex_atoms
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER    :: list_ex_kinds
      INTEGER, DIMENSION(:,:), POINTER        :: state_types 
      TYPE(section_vals_type), POINTER        :: loc_subsection
      TYPE(section_vals_type), POINTER        :: print_loc_subsection
      TYPE(section_vals_type), POINTER        :: eri_mme_subsection

   END TYPE xas_tdp_control_type
   
!> *************************************************************************************************
!> \brief Type containing informations such as inputs and results for TDP XAS calculations
!> \param state_type_char an array containing the general donor state types as char (1s, 2s, 2p, ...)
!> \param nex_atoms number of excited atoms
!> \param nex_kinds number of excited kinds
!> \param ex_atom_indices array containing the indices of the excited atoms
!> \param ex_kind_indices array containing the indices of the excited kinds
!> \param state_types columns correspond to the different donor states of each excited atom
!> \param qs_loc_env the environemnt type dealing with the possible localization of donor orbitals
!> \param mos_of_ex_atoms links lowest energy MOs to excited atoms. Elements of value 1 mark the
!>        association between the MO irow and the excited atom icolumn
!> \param ri_inv_mats collection of kind specific electron repulsion inverse matrices for each
!>        excited kinds, based on basis functions of the RI auxiliary basis used for Coulomb and
!>        exchange kernels. TODO: adapt the type when switching to dbcsr_tensor
!> \param q_projector the projector on the unperturbed, unoccupied ground state as a dbcsr matrix
!> \param dipmat the dbcsr matrices containing the dipole in x,y,z directions evaluated on the 
!>        contracted spherical gaussians. It can either be in the length or the velocity 
!>        representation. For length representation, it has to be computed once with the origin on 
!>        each excited atom
!> \param ri_coul3 the 3D array (TODO: dbcsr tensor in the future) containing the 3-center coulomb
!>        integrals (ab|P) with the RI basis centered on P. Updated for each new excited atom
!> *************************************************************************************************
   TYPE xas_tdp_env_type
      CHARACTER(len=2), DIMENSION(3)          :: state_type_char
      INTEGER                                 :: nex_atoms
      INTEGER                                 :: nex_kinds
      INTEGER, DIMENSION(:), POINTER          :: ex_atom_indices
      INTEGER, DIMENSION(:), POINTEr          :: ex_kind_indices
      INTEGER, DIMENSION(:,:), POINTER        :: state_types
      REAL(dp), DIMENSION(:,:,:), POINTER     :: ri_coul3
      TYPE(donor_state_type), DIMENSION(:),&
                        POINTER               :: donor_states
      INTEGER, DIMENSION(:,:), POINTER        :: mos_of_ex_atoms
      TYPE(qs_loc_env_new_type), POINTER      :: qs_loc_env
      TYPE(cp_fm_p_type), DIMENSION(:), &
                                    POINTER   :: ri_inv_mats
      TYPE(dbcsr_type), POINTER               :: q_projector
      TYPE(dbcsr_p_type), DIMENSION(:), &
                                    POINTER   :: dipmat
   END TYPE xas_tdp_env_type

!> *************************************************************************************************
!> \brief Type containing informations about a single donor state
!> \param at_index the index of the atom to which the state belongs
!> \param kind_index the index of the atomic kind to which the state belongs
!> \param at_symbol the chemical symbol of the atom to which the state belongs
!> \param state_type whether this is a 1s, 2s, etc state
!> \param energy_evals the energy eigenvalue of the donor state (if canonical orbitals), or the set
!>        of epsilon_ij (if  non-canonical MOs)
!> \param mo_indices indices of associated MOs. Greater than 1 when not a s-type state.
!> \param lr_coeffs solutions of the linear-response TDDFT equation retricted to this state (vector)
!> \param lr_evals excitation energies => the eigenvalues of the linear response equation
!> \param osc_str oscillator strengths for the different excitation energies
!> \param matrix_tdp the dbcsr matrix to be diagonalized to solve the problem
!> \param metric the metric of the linear response problem M*c = omega*S*c
!> *************************************************************************************************
   TYPE donor_state_type
      INTEGER                                 :: at_index
      INTEGER                                 :: kind_index
      CHARACTER(len=2)                        :: at_symbol
      INTEGER                                 :: state_type
      REAL(dp), DIMENSION(:, :), POINTER       :: energy_evals
      INTEGER, DIMENSION(:), POINTER          :: mo_indices
      TYPE(cp_fm_type), POINTER               :: lr_coeffs
      REAL(dp), DIMENSION(:), POINTER         :: lr_evals
      REAL(dp), DIMENSION(:), POINTER         :: osc_str
      TYPE(dbcsr_type), POINTER               :: matrix_tdp
      TYPE(dbcsr_type), POINTER               :: metric
   END TYPE donor_state_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_types'

! *** Public data types ***
   PUBLIC :: xas_tdp_env_type, donor_state_type, xas_tdp_control_type

! *** Public subroutines ***
   PUBLIC :: set_donor_state, free_ds_memory, &
             xas_tdp_env_create, xas_tdp_env_release, set_xas_tdp_env, &
             xas_tdp_control_create, xas_tdp_control_release, read_xas_tdp_control

CONTAINS

! **************************************************************************************************
!> \brief Creates and initializes the xas_tdp_control_type
!> \param xas_tdp_control the type to initialize
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_create(xas_tdp_control)
      
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_control_create', &
                                                            routineP = moduleN//':'//routineN
      
      CPASSERT(.NOT. ASSOCIATED(xas_tdp_control))
      ALLOCATE (xas_tdp_control)
      
      xas_tdp_control%dipole_form = xas_dip_vel
      xas_tdp_control%define_excited = xas_tdp_by_index
      xas_tdp_control%n_search = -1
      xas_tdp_control%n_gauss = 3
      xas_tdp_control%check_donor_states = .FALSE.
      xas_tdp_control%tamm_dancoff = .FALSE.
      xas_tdp_control%full_diag = .FALSE.
      xas_tdp_control%do_singlet = .TRUE.
      NULLIFY (xas_tdp_control%state_types)
      NULLIFY (xas_tdp_control%list_ex_atoms)
      NULLIFY (xas_tdp_control%list_ex_kinds)
      NULLIFY (xas_tdp_control%loc_subsection)
      NULLIFY (xas_tdp_control%print_loc_subsection)
      NULLIFY (xas_tdp_control%eri_mme_subsection)

   END SUBROUTINE xas_tdp_control_create

! **************************************************************************************************
!> \brief Releases the xas_tdp_control_type
!> \param xas_tdp_control the type to release
! **************************************************************************************************
   SUBROUTINE xas_tdp_control_release(xas_tdp_control)
      
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xas_tdp_control_release', &
                                                            routineP = moduleN//':'//routineN
      
      IF (ASSOCIATED(xas_tdp_control)) THEN
         IF (ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_atoms)
         END IF
         IF (ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
            DEALLOCATE (xas_tdp_control%list_ex_kinds)
         END IF
         IF (ASSOCIATED(xas_tdp_control%state_types)) THEN
            DEALLOCATE(xas_tdp_control%state_types)
         END IF
         IF (ASSOCIATED(xas_tdp_control%loc_subsection)) THEN
            DEALLOCATE(xas_tdp_control%loc_subsection)
         END IF
         IF (ASSOCIATED(xas_tdp_control%print_loc_subsection)) THEN
            DEALLOCATE(xas_tdp_control%print_loc_subsection)
         END IF
         IF (ASSOCIATED(xas_tdp_control%eri_mme_subsection)) THEN
            DEALLOCATE(xas_tdp_control%eri_mme_subsection)
         END IF
         DEALLOCATE (xas_tdp_control)
      END IF

   END SUBROUTINE xas_tdp_control_release

! **************************************************************************************************
!> \brief Reads the inputs and stores in xas_tdp_control_type
!> \param xas_tdp_control the type where inputs are stored
!> \param xas_tdp_section the section from which input are read
! **************************************************************************************************
   SUBROUTINE read_xas_tdp_control(xas_tdp_control, xas_tdp_section)

      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(section_vals_type), POINTER                   :: xas_tdp_section
      
      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_xas_tdp_control', &
                                                            routineP = moduleN//':'//routineN

      INTEGER                                            :: irep, nrep, nexc, excitation 
      CHARACTER(len=default_string_length),&
                     DIMENSION(:), POINTER               :: k_list
      INTEGER, DIMENSION(:), POINTER                     :: a_list, t_list

      NULLIFY (k_list,a_list,t_list)
      
      CALL section_vals_val_get(xas_tdp_section, "DEFINE_EXCITED", &
                                i_val=xas_tdp_control%define_excited)

      CALL section_vals_val_get(xas_tdp_section, "DIPOLE_FORM", i_val=xas_tdp_control%dipole_form)

      CALL section_vals_val_get(xas_tdp_section, "ONLY_CHECK_DONOR_STATES", &
                                l_val=xas_tdp_control%check_donor_states)

      CALL section_vals_val_get(xas_tdp_section, "TAMM_DANCOFF", &
                                l_val=xas_tdp_control%tamm_dancoff)

      CALL section_vals_val_get(xas_tdp_section, "FULL_DIAG", &
                                l_val=xas_tdp_control%full_diag)

      CALL section_vals_val_get(xas_tdp_section, "EXCITATIONS", i_val=excitation)
      IF (excitation == tddfpt_singlet) THEN
         xas_tdp_control%do_singlet = .TRUE.
      ELSE
         xas_tdp_control%do_singlet = .FALSE.
      END IF

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_kinds)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            ALLOCATE (xas_tdp_control%list_ex_kinds(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN                   

            CALL section_vals_val_get(xas_tdp_section, "KIND_LIST", c_vals=k_list)
         
            IF (ASSOCIATED(k_list)) THEN
               nexc = SIZE(k_list)
               ALLOCATE (xas_tdp_control%list_ex_kinds(nexc))
               xas_tdp_control%list_ex_kinds = k_list
            END IF

         END IF
      END IF

      IF (.NOT. ASSOCIATED(xas_tdp_control%list_ex_atoms)) THEN
         IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN
            
            ALLOCATE (xas_tdp_control%list_ex_atoms(0))

         ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN

            CALL section_vals_val_get(xas_tdp_section, "ATOM_LIST", i_vals=a_list)

            IF (ASSOCIATED(a_list)) THEN
               nexc = SIZE(a_list)
               CALL reallocate(xas_tdp_control%list_ex_atoms, 1, nexc)
               xas_tdp_control%list_ex_atoms = a_list
            END IF

         END IF
      END IF

      CALL section_vals_val_get(xas_tdp_section, "STATE_TYPES", n_rep_val=nrep)

      IF(.NOT. ASSOCIATED(xas_tdp_control%state_types)) THEN
         ALLOCATE (xas_tdp_control%state_types(nrep,nexc))
         DO irep = 1,nrep
            CALL section_vals_val_get(xas_tdp_section, "STATE_TYPES", i_rep_val=irep, i_vals=t_list)
            IF (SIZE(t_list) .NE. nexc) THEN
               CPABORT("The STATE_TYPES keywords do not have the correct number of entries.")
            END IF
            xas_tdp_control%state_types(irep,:) = t_list
         END DO
      END IF
      
      CALL section_vals_val_get(xas_tdp_section, "STATE_SEARCH", i_val=xas_tdp_control%n_search)

      CALL section_vals_val_get(xas_tdp_section, "NGAUSS", i_val=xas_tdp_control%n_gauss)

      xas_tdp_control%loc_subsection => section_vals_get_subs_vals(xas_tdp_section, "LOCALIZE")       
      xas_tdp_control%print_loc_subsection => section_vals_get_subs_vals(&
                                                         xas_tdp_control%loc_subsection, "PRINT")                      
      
      xas_tdp_control%eri_mme_subsection => section_vals_get_subs_vals(xas_tdp_section, "ERI_MME")   

      IF (.NOT. ASSOCIATED(xas_tdp_control%eri_mme_subsection)) THEN
         CPABORT("The ERI_MME subsection is needed for the computation of the Coulomb kernel")
      END IF

   END SUBROUTINE read_xas_tdp_control

!> *************************************************************************************************
!> \brief Creates a TDP XAS environment type
!> \param xas_tdp_env the type to create
!> *************************************************************************************************
   SUBROUTINE xas_tdp_env_create(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_env_create', &
         routineP = moduleN//':'//routineN

      ALLOCATE (xas_tdp_env)

      xas_tdp_env%nex_atoms = 1
      xas_tdp_env%nex_kinds = 1

      NULLIFY (xas_tdp_env%ex_atom_indices)
      NULLIFY (xas_tdp_env%ex_kind_indices)
      NULLIFY (xas_tdp_env%state_types)
      NULLIFY (xas_tdp_env%donor_states)
      NULLIFY (xas_tdp_env%qs_loc_env)
      NULLIFY (xas_tdp_env%mos_of_ex_atoms)
      NULLIFY (xas_tdp_env%ri_inv_mats)
      NULLIFY (xas_tdp_env%q_projector)
      NULLIFY (xas_tdp_env%dipmat)
      NULLIFY (xas_tdp_env%ri_coul3)

!     Putting the state types as char manually
      xas_tdp_env%state_type_char(1) = "1s"
      xas_tdp_env%state_type_char(2) = "2s"
      xas_tdp_env%state_type_char(3) = "2p"

   END SUBROUTINE xas_tdp_env_create

!> *************************************************************************************************
!> \brief Releases the TDP XAS environment type
!> \param xas_tdp_env the type to release
!> *************************************************************************************************
   SUBROUTINE xas_tdp_env_release(xas_tdp_env)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(xas_tdp_env)) THEN
         IF (ASSOCIATED(xas_tdp_env%ex_atom_indices)) THEN
            DEALLOCATE (xas_tdp_env%ex_atom_indices)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ex_kind_indices)) THEN
            DEALLOCATE(xas_tdp_env%ex_kind_indices)
         END IF

         IF (ASSOCIATED(xas_tdp_env%state_types)) THEN
            DEALLOCATE (xas_tdp_env%state_types)
         END IF
         IF (ASSOCIATED(xas_tdp_env%donor_states)) THEN
            CALL deallocate_donor_state_set(xas_tdp_env%donor_states)
         END IF
         IF (ASSOCIATED(xas_tdp_env%qs_loc_env)) THEN
            CALL qs_loc_env_release(xas_tdp_env%qs_loc_env)
         END IF
         IF (ASSOCIATED(xas_tdp_env%mos_of_ex_atoms)) THEN
            DEALLOCATE(xas_tdp_env%mos_of_ex_atoms)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_inv_mats)) THEN
            DO i =1,SIZE(xas_tdp_env%ri_inv_mats)
               CALL cp_fm_release(xas_tdp_env%ri_inv_mats(i)%matrix)
            END DO
            DEALLOCATE(xas_tdp_env%ri_inv_mats)
         END IF
         IF (ASSOCIATED(xas_tdp_env%q_projector)) THEN
            CALL dbcsr_release(xas_tdp_env%q_projector)
            DEALLOCATE(xas_tdp_env%q_projector)
         END IF
         IF (ASSOCIATED(xas_tdp_env%dipmat)) THEN
            CALL dbcsr_deallocate_matrix_set(xas_tdp_env%dipmat)
         END IF
         IF (ASSOCIATED(xas_tdp_env%ri_coul3)) THEN
            DEALLOCATE(xas_tdp_env%ri_coul3)
         END IF
         DEALLOCATE (xas_tdp_env)
      END IF
   END SUBROUTINE xas_tdp_env_release

!> *************************************************************************************************
!> \brief Sets values of selected variables within the TDP XAS environment type
!> *************************************************************************************************
   SUBROUTINE set_xas_tdp_env(xas_tdp_env, nex_atoms, nex_kinds)

      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      INTEGER, INTENT(IN), OPTIONAL                      :: nex_atoms
      INTEGER, INTENT(IN), OPTIONAL                      :: nex_kinds

      CHARACTER(len=*), PARAMETER :: routineN = 'set_xas_tdp_env', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(xas_tdp_env))

      IF (PRESENT(nex_atoms)) xas_tdp_env%nex_atoms = nex_atoms
      IF (PRESENT(nex_kinds)) xas_tdp_env%nex_kinds = nex_kinds
      
   END SUBROUTINE set_xas_tdp_env

!> *************************************************************************************************
!> \brief sets specified values of the donor state type
!> \param donor_state the type which values should be set
!> \param at_index ...
!> \param at_symbol ...
!> \param state_type ...
!> *************************************************************************************************
   SUBROUTINE set_donor_state(donor_state, at_index, at_symbol, kind_index, state_type)

      TYPE(donor_state_type), POINTER                    :: donor_state
      INTEGER, INTENT(IN), OPTIONAL                      :: at_index
      CHARACTER(len=2), INTENT(IN), OPTIONAL             :: at_symbol
      INTEGER, INTENT(IN), OPTIONAL                      :: kind_index, state_type

      CHARACTER(len=*), PARAMETER :: routineN = 'set_donor_state', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(donor_state))

      IF (PRESENT(at_index)) donor_state%at_index = at_index
      IF (PRESENT(kind_index)) donor_state%kind_index = kind_index
      IF (PRESENT(state_type)) donor_state%state_type = state_type
      IF (PRESENT(at_symbol)) donor_state%at_symbol = at_symbol

   END SUBROUTINE set_donor_state

!> *************************************************************************************************
!> \brief Deallocate a set of donor states
!> \param donor_state_set the set of donor states to deallocate
!> *************************************************************************************************
   SUBROUTINE deallocate_donor_state_set(donor_state_set)

      TYPE(donor_state_type), DIMENSION(:), POINTER      :: donor_state_set

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_donor_state_set', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(donor_state_set)) THEN
         DO i = 1, SIZE(donor_state_set)

            IF (ASSOCIATED(donor_state_set(i)%lr_coeffs)) THEN
               CALL cp_fm_release(donor_state_set(i)%lr_coeffs)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%lr_evals)) THEN
               DEALLOCATE (donor_state_set(i)%lr_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%osc_str)) THEN
               DEALLOCATE (donor_state_set(i)%osc_str)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%energy_evals)) THEN
               DEALLOCATE (donor_state_set(i)%energy_evals)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%mo_indices)) THEN
               DEALLOCATE (donor_state_set(i)%mo_indices)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%matrix_tdp)) THEN
               CALL dbcsr_release(donor_state_set(i)%matrix_tdp)
               DEALLOCATE(donor_state_set(i)%matrix_tdp)
            END IF

            IF (ASSOCIATED(donor_state_set(i)%metric)) THEN  
               CALL dbcsr_release(donor_state_set(i)%metric)      
               DEALLOCATE(donor_state_set(i)%metric)           
            END IF
                        
         END DO
         DEALLOCATE (donor_state_set)
      END IF

   END SUBROUTINE deallocate_donor_state_set

!> *************************************************************************************************
!> \brief Deallocate a donor_state's heavy attributes
!> \param donor_state ...
!> *************************************************************************************************
   SUBROUTINE free_ds_memory(donor_state)
      
      TYPE(donor_state_type), POINTER                    :: donor_state

      CHARACTER(len=*), PARAMETER :: routineN = "free_ds_memory", routineP = moduleN//":"//routineN

      IF (ASSOCIATED(donor_state%lr_evals)) DEALLOCATE(donor_state%lr_evals)
      IF (ASSOCIATED(donor_state%osc_str)) DEALLOCATE(donor_state%osc_str)
      IF (ASSOCIATED(donor_state%lr_coeffs)) CALL cp_fm_release(donor_state%lr_coeffs)

      IF (ASSOCIATED(donor_state%matrix_tdp)) THEN
         CALL dbcsr_release(donor_state%matrix_tdp)
         DEALLOCATE(donor_state%matrix_tdp)
      END IF

      IF (ASSOCIATED(donor_state%metric)) THEN   
         CALL dbcsr_release(donor_state%metric) 
         DEALLOCATE(donor_state%metric)  
      END IF

   END SUBROUTINE free_ds_memory

END MODULE xas_tdp_types
