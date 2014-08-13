!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief driver for the xas calculation and xas_scf for the tp method
!> \par History
!>      created 05.2005
!>      replace overlap integral routine [07.2014,JGH]
!> \author MI (05.2005)
! *****************************************************************************
MODULE xas_methods

  USE ai_contraction,                  ONLY: block_add,&
                                             contraction
  USE ai_overlap,                      ONLY: overlap_ab
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: &
       allocate_sto_basis_set, create_gto_from_sto_basis, &
       deallocate_sto_basis_set, get_gto_basis_set, gto_basis_set_type, &
       init_orb_basis_set, set_sto_basis_set, srules, sto_basis_set_type
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_new, array_nullify, array_release, &
       convert_offsets_to_sizes, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_init, cp_dbcsr_p_type, cp_dbcsr_set, &
       dbcsr_distribution_obj, dbcsr_type_antisymmetric
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_external_control,             ONLY: external_control
  USE cp_fm_pool_types,                ONLY: fm_pool_create_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: &
       do_loc_none, state_loc_list, state_loc_range, xas_1s_type, &
       xas_2p_type, xas_2s_type, xas_dip_len2, xas_dip_vel, xas_dscf, &
       xas_scf_general, xas_tp_fh, xas_tp_hh, xas_tp_xfh, xas_tp_xhh, &
       xes_tp_val
  USE input_section_types,             ONLY: section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: evolt
  USE qs_density_mixing_types,         ONLY: direct_mixing_nr,&
                                             mixing_storage_create
  USE qs_diis,                         ONLY: qs_diis_b_clear,&
                                             qs_diis_b_create
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_loc_methods,                  ONLY: qs_loc_driver,&
                                             qs_print_cubes
  USE qs_loc_types,                    ONLY: localized_wfn_control_type,&
                                             qs_loc_env_create,&
                                             qs_loc_env_new_type,&
                                             qs_loc_env_release
  USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                             qs_loc_env_init,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_io,                        ONLY: write_mo_set
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_set
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                             rRc_xyz_ao
  USE qs_scf,                          ONLY: scf_env_cleanup
  USE qs_scf_initialization,           ONLY: qs_scf_env_initialize
  USE qs_scf_types,                    ONLY: general_diag_method_nr,&
                                             ot_method_nr,&
                                             qs_scf_env_type
  USE scf_control_types,               ONLY: init_smear,&
                                             read_smear_section,&
                                             scf_control_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             set_xas_env,&
                                             xas_env_create,&
                                             xas_env_release,&
                                             xas_environment_type
  USE xas_restart,                     ONLY: xas_read_restart
  USE xas_tp_scf,                      ONLY: xas_do_tp_scf,&
                                             xes_scf_once
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_methods'

! *** Public subroutines ***

  PUBLIC :: xas

CONTAINS

! *****************************************************************************
!> \brief Driver for xas calculations
!>      The initial mos are prepared
!>      A loop on the atoms to be excited is started
!>      For each atom the state to be excited is identified
!>      An scf optimization using the TP scheme or TD-DFT is used
!>      to evaluate the spectral energies and oscillator strengths
!> \param qs_env the qs_env, the xas_env lives in
!> \param dft_control ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2005 created [MI]
!> \author MI
!> \note
!>      the iteration counter is not finilized yet
!>      only the transition potential approach is active
!>      the localization can be switched off, otherwise
!>      it uses by default the berry phase approach
!>      The number of states to be localized is xas_control%nexc_search
!>      In general only the core states are needed
! *****************************************************************************
  SUBROUTINE xas(qs_env, dft_control, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i, iat, iatom, ispin, istat, method0, &
      my_homo(2), my_nelectron(2), nao, nelectron, nexc_atoms, nexc_search, &
      nmo, nspins, output_unit, state_to_be_excited
    INTEGER, DIMENSION(:), POINTER           :: state_of_atom
    LOGICAL :: ch_method_flags, converged, failure, my_uocc(2), should_stop, &
      skip_scf, transition_potential
    REAL(dp)                                 :: maxocc, occ_estate, &
                                                xas_nelectron
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, &
                                                occupation_numbers
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuffer
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, op_sm, ostrength_sm
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: groundstate_coeff
    TYPE(cp_fm_type), POINTER                :: all_vectors, excvec_coeff, &
                                                mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, loc_section, &
                                                print_loc_section, &
                                                scf_section, xas_section
    TYPE(xas_control_type), POINTER          :: xas_control
    TYPE(xas_environment_type), POINTER      :: xas_env

    CALL timeset(routineN,handle)

    failure = .FALSE.
    transition_potential = .FALSE.
    skip_scf = .FALSE.
    converged = .TRUE.
    should_stop = .FALSE.
    ch_method_flags = .FALSE.

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    NULLIFY(xas_env, groundstate_coeff, ostrength_sm, op_sm)
    NULLIFY(excvec_coeff, state_of_atom, qs_loc_env, cell)
    NULLIFY(occupation_numbers,matrix_ks)
    NULLIFY(all_vectors,state_of_atom,xas_control)
    NULLIFY(vecbuffer,op_sm)
    NULLIFY(dft_section, xas_section, loc_section, print_loc_section)
    dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)
    xas_section => section_vals_get_subs_vals(dft_section,"XAS",error=error)
    loc_section => section_vals_get_subs_vals(xas_section,"LOCALIZE",error=error)
    print_loc_section => section_vals_get_subs_vals(loc_section,"PRINT",error=error)

    xas_control => dft_control%xas_control
    output_unit = cp_print_key_unit_nr(logger,xas_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".Log",error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "START CORE LEVEL SPECTROSCOPY CALCULATION",&
            REPEAT("=",77)
    END IF

!   Create the xas environment
    CALL get_qs_env(qs_env,xas_env=xas_env,error=error)
    IF (.NOT.ASSOCIATED(xas_env)) THEN
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
           "Create and initialize the xas environment"
       END IF
      CALL xas_env_create(xas_env, error=error)
      CALL xas_env_init(xas_env, xas_control, qs_env, xas_section, logger, error=error)
      CALL set_qs_env(qs_env,xas_env=xas_env, error=error)
      CALL xas_env_release(xas_env,error=error)
      CALL get_qs_env(qs_env,xas_env=xas_env,error=error)
    END IF

!   Initialize the type of calculation
    NULLIFY(atomic_kind_set, scf_control, mos, para_env, particle_set)
    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
         cell = cell, scf_control=scf_control,&
         matrix_ks=matrix_ks,mos=mos, para_env=para_env, &
         particle_set=particle_set ,error=error)

!   The eigenstate of the KS Hamiltonian are nedeed
    NULLIFY(mo_coeff,eigenvalues)
    IF(scf_control%use_ot) THEN
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
               "Get eigenstates and eigenvalues from ground state MOs"
       END IF
      DO ispin = 1,dft_control%nspins
        CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
             eigenvalues=eigenvalues,homo=homo)
        CALL calculate_subspace_eigenvalues(mo_coeff,&
             matrix_ks(ispin)%matrix,eigenvalues, &
             do_rotation=.TRUE.,error=error)
      END DO
    END IF

!   Set initial occupation numbers, and store the original ones
    my_homo = 0
    my_nelectron = 0
    DO ispin = 1,dft_control%nspins
      CALL get_mo_set(mos(ispin)%mo_set,nelectron=my_nelectron(ispin),maxocc=maxocc,&
                      homo=my_homo(ispin),uniform_occupation=my_uocc(ispin))
    END DO

    nspins = dft_control%nspins
    transition_potential =  .TRUE. !(xas_control%xas_method==xas_tp_hh).OR.&
 !                           (xas_control%xas_method==xas_tp_fh).OR.&
  !                          (xas_control%xas_method==xas_tp_xhh).OR.&
   !                         (xas_control%xas_method==xas_tp_xfh).OR.&
    !                        (xas_control%xas_method==xas_dscf)
    IF(nspins==1 .AND. transition_potential) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "XAS with TP method requires LSD calculations")
    END IF

!   Set of states among which there is the state to be excited
    CALL get_mo_set(mos(1)%mo_set,nao=nao,homo=homo)
    IF(xas_control%nexc_search < 0) xas_control%nexc_search = homo
    nexc_search = xas_control%nexc_search

    CALL get_xas_env(xas_env=xas_env,&
         state_of_atom=state_of_atom,all_vectors=all_vectors,&
         groundstate_coeff=groundstate_coeff,excvec_coeff=excvec_coeff,&
         nexc_atoms=nexc_atoms,occ_estate=occ_estate,xas_nelectron=xas_nelectron,error=error)

    CALL set_xas_env(xas_env=xas_env,nexc_search=nexc_search,error=error)

!   SCF for only XES using occupied core and empty homo (only one SCF)
!   Probably better not to do the localization in this case, but only single out the
!   core orbital for the specific atom for which the spectrum is computed
    IF(xas_control%xas_method==xes_tp_val .AND. &
                   xas_control%xes_core_occupation==1.0_dp) THEN
      IF (output_unit>0) WRITE(UNIT=output_unit,FMT='(/,/,T10,A)') &
        "START Core Level Spectroscopy Calculation for the Emission Spectrum"
      IF(xas_control%xes_homo_occupation == 1) THEN
        IF (output_unit>0) WRITE(UNIT=output_unit,FMT='(T10,A,/,A)')&
         "The core state is fully occupied and XES from ground state calculation.",&
         " No SCF is needed, MOS already available"
      ELSE IF(xas_control%xes_homo_occupation == 0) THEN
        IF (output_unit>0) WRITE(UNIT=output_unit,FMT='(T10,A,/,A)')&
         "The core state is fully occupied and the homo is empty",&
         " (final state of the core hole decay). Only one SCF is needed (not one per atom)"
      END IF
      skip_scf = .TRUE.

      CALL set_xas_env(xas_env=xas_env,xas_estate=-1, homo_occ=xas_control%xes_homo_occupation, &
            error=error)
      CALL xes_scf_once(qs_env,xas_env,scf_section,converged,should_stop,error=error)

      IF(converged .AND. .NOT. should_stop .AND. xas_control%xes_homo_occupation == 0) THEN

        IF (output_unit>0) WRITE(UNIT=output_unit,FMT='(/,T10,A,I6)') &
              "SCF with empty homo converged "
      ELSE IF(.NOT. converged .OR. should_stop) THEN
        IF (output_unit>0) WRITE(UNIT=output_unit,FMT='(/,T10,A,I6)') &
              "SCF with empty homo NOT converged"
        GOTO 1000
      END IF

    END IF

   !Define the qs_loc_env : to find centers, spread and possibly localize them
    CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env,error=error)
    IF (qs_loc_env%do_localize) THEN
      IF(output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
              "Localize a sub-set of MOs with spin alpha,"//&
              " to better identify the core states"
          IF(&
              qs_loc_env%localized_wfn_control%set_of_states == state_loc_range ) THEN
              WRITE (UNIT=output_unit,FMT="( A , I7, A, I7)") " The sub-set contains states from ",&
                  qs_loc_env%localized_wfn_control%lu_bound_states(1,1), " to ",&
                  qs_loc_env%localized_wfn_control%lu_bound_states(2,1)
         ELSEIF (qs_loc_env%localized_wfn_control%set_of_states==state_loc_list) THEN
              WRITE (UNIT=output_unit,FMT="( A )") " The sub-set contains states given in the input list"
         END IF

      END IF
      CALL qs_loc_driver(qs_env,qs_loc_env,loc_section,print_loc_section,myspin=1,error=error)
    END IF

    ! Assign the character of the selected core states
    ! through the overlap with atomic-like states
    CALL cls_assign_core_states(xas_control,xas_env,qs_loc_env%localized_wfn_control,&
         qs_env,error=error)

    IF(skip_scf) THEN
      CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff)
      CALL cp_fm_to_fm(mo_coeff,all_vectors,ncol=nexc_search,&
             source_start=1,target_start=1)
    END IF

    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (op_sm(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  ! copy the coefficients of the mos in a temporary fm with the right structure
    IF(transition_potential) THEN
      CPPrecondition(ASSOCIATED(groundstate_coeff),cp_failure_level,routineP,error,failure)
      DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
       CALL cp_fm_to_fm(mo_coeff,groundstate_coeff(ispin)%matrix,nmo,1,1)
      END DO

      ! Calculate the operator
      CALL get_xas_env(xas_env=xas_env,ostrength_sm=ostrength_sm,error=error)
      DO i = 1,3
        NULLIFY(op_sm(i)%matrix)
        op_sm(i)%matrix => ostrength_sm(i)%matrix
      END DO
      IF(xas_control%dipole_form==xas_dip_vel) THEN
         CALL p_xyz_ao(op_sm,qs_env,error=error)
      END IF
    END IF
    ! for td-dft-xas the oscillator strength should be obtained from the linear response orbitals

    ! DO SCF if required
    NULLIFY(scf_env)
    CALL get_qs_env(qs_env,scf_env=scf_env,error=error)
     ! OT cannot be used with different occupation numbers in the core
     ! Therefore in case of XAS the general_diag_method_nr is imposed
    ch_method_flags = .FALSE.
    IF( occ_estate .LT. 1.0_dp .AND. scf_env%method == ot_method_nr) THEN
  CPPostcondition(xas_control%scf_method==xas_scf_general,cp_warning_level,routineP,error,failure)
       method0 = scf_env%method
       scf_env%method = general_diag_method_nr
       ch_method_flags = .TRUE.
    END IF

    !some further setup depending on the mixing method
    scf_env%mixing_method=xas_env%mixing_method
    scf_env%p_mix_alpha = xas_env%mixing_store%alpha
    scf_env%skip_diis = .FALSE.
    IF(xas_env%mixing_method>direct_mixing_nr) THEN
      scf_env%skip_diis = .TRUE.
      IF(xas_env%mixing_store%beta==0.0_dp) THEN
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
              routineP,"CLS: Mixing employing the Kerker damping factor needs BETA /= 0.0",error,failure)
      END IF
    END IF
    IF (xas_env%mixing_method==direct_mixing_nr) THEN
        IF(xas_env%eps_diis<qs_env%scf_control%eps_scf) THEN
          scf_env%skip_diis = .TRUE.
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                        "CLS: the DIIS scheme is disabled, since EPS_DIIS < EPS_SCF",&
                        error,failure)
        END IF
    END IF

    DO iat = 1,nexc_atoms
      iatom = xas_env%exc_atoms(iat)
      ! determine which state has to be excited in the global list
      state_to_be_excited = state_of_atom(iat)

      ! Take the state_to_be_excited vector from the full set and copy into excvec_coeff
      CALL get_mo_set(mos(1)%mo_set, &
           occupation_numbers=occupation_numbers,homo=homo,nmo=nmo,nelectron=nelectron)
      IF(xas_control%xas_method==xas_dscf .OR. xas_control%xas_method==xas_tp_xhh &
        .OR.  xas_control%xas_method==xas_tp_xfh) THEN
         CALL cp_assert(nmo>(nelectron+1),cp_failure_level,cp_assertion_failed,&
              routineP,"CLS: the required method needs added_mos to the ground state",error,failure)
      END IF
      occupation_numbers(1:homo)  = 1.0_dp
      ! If the restart file for this atom exists, the mos and the
      ! occupation numbers are overwritten
      ! It is necessary that the restart is for the same xas method
      ! otherwise the number of electrons and the occupation numbers
      ! may not  be consistent
      IF(xas_control%xas_restart) THEN
        CALL xas_read_restart(xas_env,xas_section, qs_env, xas_control%xas_method, iatom,&
             state_to_be_excited,error=error)
      END IF

      CALL set_xas_env(xas_env=xas_env,xas_estate=state_to_be_excited,error=error)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
      CPPrecondition(ASSOCIATED(excvec_coeff),cp_failure_level,routineP,error,failure)
      CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,state_to_be_excited,&
           nao,1,transpose=.TRUE.,error=error)
      CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
           nao,1,transpose=.TRUE.,error=error)

      IF (transition_potential) THEN
        IF(.NOT. skip_scf) THEN
          IF (output_unit>0) THEN
            IF(xas_control%xas_method==xas_dscf) THEN
              WRITE(UNIT=output_unit,FMT='(/,/,T10,A,I6)') &
              "START DeltaSCF for the first excited state from the core state of ATOM ", iatom
            ELSE
              WRITE(UNIT=output_unit,FMT='(/,/,T10,A,I6)') &
               "Start Core Level Spectroscopy Calculation with TP approach for ATOM ", iatom
              WRITE(UNIT=output_unit,FMT='(T10,A,f10.4)') "Occupation of the core orbital",&
                occ_estate
              IF(xas_nelectron==nelectron) THEN
                WRITE(UNIT=output_unit,FMT='(T10,A,f10.4)') &
                 "Occupation of lowest virtual orbitals to keep the same number of electrons"//&
                 " as in the groundstate"
              END IF
            END IF
          END IF

          CALL qs_scf_env_initialize(scf_env=scf_env,qs_env=qs_env,&
                                     do_xas_tp=.TRUE., error=error)
          ! if using mo_coeff_b then copy to fm
          DO ispin=1,SIZE(mos)!fm->dbcsr
            IF(mos(1)%mo_set%use_mo_coeff_b)THEN!fm->dbcsr
            !write(*,*) routinen//' copy_dbcsr_to_fm',__LINE__
              CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b,mos(ispin)%mo_set%mo_coeff,error=error)!fm->dbcsr
            ENDIF!fm->dbcsr
          ENDDO!fm->dbcsr

          IF(.NOT.scf_env%skip_diis) THEN
            IF (.NOT.ASSOCIATED(scf_env%scf_diis_buffer)) THEN
              CALL qs_diis_b_create(scf_env%scf_diis_buffer,nbuffer=scf_control%max_diis,error=error)
            END IF
            CALL qs_diis_b_clear(scf_env%scf_diis_buffer,error=error)
          END IF

          CALL xas_do_tp_scf(dft_control,xas_env,iatom,scf_env,qs_env,&
               xas_section,scf_section,converged,should_stop,error=error)

          CALL external_control(should_stop,"CLS",target_time=qs_env%target_time,&
               start_time=qs_env%start_time,error=error)
          IF(should_stop)THEN
             CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)
             GO TO 1000
          END IF

        END IF
       ! SCF DONE

!   *** Write last wavefunction to screen ***
        DO ispin=1,dft_control%nspins
           CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,4,&
                             dft_section,error=error)
        ENDDO

      ELSE
        ! Core level spectroscopy by TDDFT is not yet implemented
        ! the states defined by the rotation are the ground state orbitals
        ! the initial state from which I excite should be localized
        ! I take the excitations from lumo to nmo
      END IF

      IF(converged) THEN
 WRITE(*,*) 'qui si',  iat
!        CALL cls_prepare_states(xas_control,xas_env,qs_env,iatom,&
!             xas_section,output_unit,error=error)

        CALL cls_calculate_spectrum(xas_control,xas_env,qs_env,xas_section,&
             iatom,error=error)
      ELSE
        IF (output_unit>0) WRITE(UNIT=output_unit,FMT='(/,/,T10,A,I6)') &
              "SCF with core hole NOT converged for ATOM ", iatom
      END IF

      IF(.NOT.skip_scf) THEN
        ! Reset the initial core orbitals.
        ! The valence orbitals are taken from the last SCF,
        ! it should be a better initial guess
         DO ispin = 1,nspins
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
                 nao=nao, nmo=nmo)
            CALL cp_fm_to_fm(groundstate_coeff(ispin)%matrix,mo_coeff,nmo,1,1)
         END DO
         IF(iat == nexc_atoms) THEN
             CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)
         END IF
      END IF

    END DO

1000  CONTINUE  ! END of Calculation

    ! Set back the original method for the optimize
    IF(ch_method_flags) THEN
      scf_env%method = method0
    END IF

  ! Release what has to be released
    IF(ASSOCIATED(vecbuffer)) THEN
        DEALLOCATE(vecbuffer,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(op_sm,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF


    DO ispin = 1,dft_control%nspins
       CALL set_mo_set(mos(ispin)%mo_set, homo=my_homo(ispin),&
            uniform_occupation=my_uocc(ispin), error=error)
    END DO

    xas_env % xas_estate    = -1

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "END CORE LEVEL SPECTROSCOPY CALCULATION",&
            REPEAT("=",77)
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,xas_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE  xas

! *****************************************************************************
!> \brief allocate and initialize the structure needed for the xas calculation
!> \param xas_env the environment for XAS  calculations
!> \param xas_control parameters for the xas calculation
!> \param qs_env the qs_env, the xas_env lives in
!> \param xas_section ...
!> \param logger ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE xas_env_init(xas_env, xas_control, qs_env, xas_section, logger, error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(xas_control_type)                   :: xas_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xas_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_env_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name_sto
    INTEGER :: homo, i, iat, iatom, ik, ikind, ispin, istat, j, l, lfomo, &
      n_mo(2), n_rep, nao, natom, ncubes, nelectron, nexc_atoms, nexc_search, &
      nj, nk, nkind, nmo, nmoloc(2), norb(0:3), nsgf_gto, nsgf_sto, nspins, &
      nvirtual, nvirtual2, nwork
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, kind_type_tmp, &
                                                kind_z_tmp, last_sgf
    INTEGER, DIMENSION(4, 7)                 :: ne
    INTEGER, DIMENSION(:), POINTER           :: bounds, list, lq, nq, rbs
    LOGICAL                                  :: failure, ihavethis
    REAL(dp)                                 :: eps_diis, nele, occ_estate, &
                                                zatom
    REAL(dp), DIMENSION(:), POINTER          :: sto_zet
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: loc_section, mixing_section, &
                                                smear_section
    TYPE(sto_basis_set_type), POINTER        :: sto_basis_set

    failure=.FALSE.

    n_mo(1:2) = 0
    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

      NULLIFY( atomic_kind_set, dft_control, scf_control, matrix_s, mos, mpools )
      NULLIFY( para_env, particle_set, scf_env)
      NULLIFY(qs_loc_env)
      NULLIFY(sab_orb)
      CALL get_qs_env(qs_env=qs_env, &
                atomic_kind_set=atomic_kind_set, &
                dft_control = dft_control, &
                scf_control = scf_control, &
                mpools=mpools,&
                matrix_s=matrix_s, mos=mos, &
                para_env=para_env, particle_set=particle_set,&
                scf_env=scf_env, &
                sab_orb=sab_orb,&
                dbcsr_dist=dbcsr_dist,&
                error=error)

      nexc_search = xas_control%nexc_search
      IF(nexc_search < 0) THEN
        ! ground state occupation
        CALL get_mo_set(mos(1)%mo_set,nmo=nmo,lfomo=lfomo)
        nexc_search = lfomo -1
      END IF
      nexc_atoms = xas_control%nexc_atoms
      ALLOCATE(xas_env%exc_atoms(nexc_atoms),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      xas_env%exc_atoms = xas_control%exc_atoms
      eps_diis = xas_control%eps_diis
      IF(eps_diis .LT. 0.0_dp) eps_diis = scf_control%eps_diis
      CALL set_xas_env(xas_env=xas_env,eps_diis=eps_diis,nexc_search=nexc_search,&
           nexc_atoms=nexc_atoms,error=error)

      CALL mpools_get(mpools, ao_mo_fm_pools= xas_env%ao_mo_fm_pools,error=error)

      NULLIFY(mo_coeff)
      CALL get_mo_set(mos(1)%mo_set,nao=nao,homo=homo,nmo=nmo,mo_coeff=mo_coeff, nelectron=nelectron)

      nvirtual2 = 0
      IF(xas_control%added_mos .GT. 0) THEN
         nvirtual2 = MIN(xas_control%added_mos,nao-nmo)
         xas_env%unoccupied_eps = xas_control%eps_added
         xas_env%unoccupied_max_iter = xas_control%max_iter_added
      END IF
      nvirtual = nmo + nvirtual2

      n_mo(1:2) = nmo

      ALLOCATE(xas_env%centers_wfn(3,nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%atom_of_state(nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%type_of_state(nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%state_of_atom(nexc_atoms),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%mykind_of_atom(nexc_atoms),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      nkind = SIZE(atomic_kind_set,1)
      ALLOCATE(xas_env%mykind_of_kind(nkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      xas_env%mykind_of_kind = 0

      ! create a new matrix structure nao x 1
      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
           ncol_global=1,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
      CALL cp_fm_create (xas_env%excvec_coeff, tmp_fm_struct ,error=error)
      CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=1,&
           ncol_global=nexc_search,para_env=para_env,&
           context=mo_coeff%matrix_struct%context,error=error)
      CALL cp_fm_create (xas_env%excvec_overlap, tmp_fm_struct ,error=error)
      CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

      nspins = SIZE(mos,1)

     ! initialize operators for the calculation of the oscillator strengts
      IF (xas_control%xas_method==xas_tp_hh) THEN
       occ_estate = 0.5_dp
       nele = REAL(nelectron,dp) - 0.5_dp
      ELSEIF(xas_control%xas_method==xas_tp_xhh) THEN
       occ_estate = 0.5_dp
       nele = REAL(nelectron,dp)
      ELSEIF(xas_control%xas_method==xas_tp_fh) THEN
       occ_estate = 0.0_dp
       nele = REAL(nelectron,dp) - 1.0_dp
      ELSEIF(xas_control%xas_method==xas_tp_xfh) THEN
       occ_estate = 0.0_dp
       nele = REAL(nelectron,dp)
      ELSEIF(xas_control%xas_method==xes_tp_val) THEN
       occ_estate = xas_control%xes_core_occupation
       nele = REAL(nelectron,dp)-xas_control%xes_core_occupation
      ELSEIF(xas_control%xas_method==xas_dscf) THEN
       occ_estate = 0.0_dp
       nele = REAL(nelectron,dp)
      ENDIF
      CALL set_xas_env(xas_env=xas_env,occ_estate=occ_estate,xas_nelectron=nele,&
          nvirtual2=nvirtual2,nvirtual=nvirtual,error=error)

      ! Initialize the list of orbitals for cube files printing
      IF (BTEST(cp_print_key_should_output(logger%iter_info,xas_section,&
                    "PRINT%CLS_FUNCTION_CUBES",error=error),cp_p_file)) THEN
        NULLIFY(bounds,list)
        CALL section_vals_val_get(xas_section,&
             "PRINT%CLS_FUNCTION_CUBES%CUBES_LU_BOUNDS",&
             i_vals=bounds,error=error)
        ncubes = bounds(2) - bounds(1)  + 1
        IF(ncubes > 0 ) THEN
          ALLOCATE( xas_control%list_cubes(ncubes),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          DO ik = 1,ncubes
            xas_control%list_cubes(ik) = bounds(1) + (ik-1)
          END DO
        END IF

        IF(.NOT. ASSOCIATED(xas_control%list_cubes)) THEN
          CALL section_vals_val_get(xas_section, &
               "PRINT%CLS_FUNCTION_CUBES%CUBES_LIST",&
               n_rep_val=n_rep,error=error)
          ncubes = 0
          DO ik = 1,n_rep
            NULLIFY(list)
            CALL section_vals_val_get(xas_section,&
                 "PRINT%CLS_FUNCTION_CUBES%CUBES_LIST",&
                 i_rep_val=ik,i_vals=list,error=error)
            IF(ASSOCIATED(list)) THEN
              CALL reallocate(xas_control%list_cubes,1,ncubes+ SIZE(list))
              DO i = 1, SIZE(list)
                xas_control%list_cubes(i+ncubes) = list(i)
              END DO
              ncubes = ncubes + SIZE(list)
            END IF
          END DO  ! ik
        END IF

        IF(.NOT. ASSOCIATED(xas_control%list_cubes)) THEN
          ncubes = MAX(10,xas_control%added_mos/10)
          ncubes = MIN(ncubes,xas_control%added_mos)
          ALLOCATE( xas_control%list_cubes(ncubes),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ik = 1,ncubes
            xas_control%list_cubes(ik) = homo + ik
          END DO
        END IF
      ELSE
        NULLIFY(xas_control%list_cubes)
      ENDIF

      NULLIFY(tmp_fm_struct)
      nwork = nvirtual
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,ncol_global=nwork,&
           para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
      CALL cp_fm_create (xas_env%fm_work,tmp_fm_struct,error=error)
      CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

      ALLOCATE (xas_env%groundstate_coeff(nspins), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO ispin = 1,nspins
        NULLIFY(xas_env%groundstate_coeff(ispin)%matrix)
        CALL get_mo_set(mos(ispin)%mo_set,nao=nao,nmo=nmo)
        CALL fm_pool_create_fm(xas_env%ao_mo_fm_pools(ispin)%pool,&
             xas_env%groundstate_coeff(ispin)%matrix,&
           name="xas_env%mo0"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)
      END DO  ! ispin

      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=1,&
           ncol_global=nvirtual,para_env=para_env,&
           context=mo_coeff%matrix_struct%context,error=error)
      ALLOCATE (xas_env%dip_fm_set(2,3),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i = 1,3
        DO j = 1,2
           NULLIFY(xas_env%dip_fm_set(j,i)%matrix)
           CALL cp_fm_create (xas_env%dip_fm_set(j,i)%matrix, tmp_fm_struct ,error=error)
        END DO
      END DO
      CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

      !Array to store all the eigenstates: occupied and the required not occupied
      IF(nvirtual2 .GT. 0) THEN
        ALLOCATE(xas_env%unoccupied_evals(nvirtual2), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(tmp_fm_struct)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
             ncol_global=nvirtual2,&
             para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
        NULLIFY(xas_env%unoccupied_orbs)
        CALL cp_fm_create (xas_env%unoccupied_orbs,tmp_fm_struct,error=error)
        CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)
      END IF

      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao, &
           ncol_global=nvirtual,&
           para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
      NULLIFY(xas_env%all_vectors)
      CALL cp_fm_create (xas_env%all_vectors,tmp_fm_struct,error=error)
      CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

      ! Array to store all the energies needed  for the spectrum
      ALLOCATE(xas_env%all_evals(nvirtual), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      IF(xas_control%dipole_form==xas_dip_len2) THEN
        CALL cp_dbcsr_allocate_matrix_set(xas_env%ostrength_sm,3,error=error)
        DO i = 1,3
           ALLOCATE(xas_env%ostrength_sm(i)%matrix)
           CALL cp_dbcsr_init(xas_env%ostrength_sm(i)%matrix,error=error)
           CALL cp_dbcsr_copy(xas_env%ostrength_sm(i)%matrix,matrix_s(1)%matrix,&
                "xas_env%ostrength_sm"//&
                "-"//TRIM(ADJUSTL(cp_to_string(i))),error=error)
           CALL cp_dbcsr_set(xas_env%ostrength_sm(i)%matrix,0.0_dp,error=error)
        END DO
      ELSEIF(xas_control%dipole_form==xas_dip_vel) THEN
         !
         ! prepare for allocation
         natom = SIZE(particle_set,1)
         ALLOCATE (first_sgf(natom),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE (last_sgf(natom),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         CALL get_particle_set(particle_set=particle_set,&
              first_sgf=first_sgf,&
              last_sgf=last_sgf,error=error)
         ALLOCATE (rbs(natom), STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
         CALL array_nullify (row_blk_sizes)
         CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
         DEALLOCATE (first_sgf,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE (last_sgf,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         !
         !
         CALL cp_dbcsr_allocate_matrix_set(xas_env%ostrength_sm,3,error=error)
         ALLOCATE(xas_env%ostrength_sm(1)%matrix)
         CALL cp_dbcsr_init(xas_env%ostrength_sm(1)%matrix,error=error)
         CALL cp_dbcsr_create(matrix=xas_env%ostrength_sm(1)%matrix, &
              name="xas_env%ostrength_sm", &
              dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
              row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
              nblks=0, nze=0, mutable_work=.TRUE., &
              error=error)
         CALL cp_dbcsr_alloc_block_from_nbl(xas_env%ostrength_sm(1)%matrix,sab_orb,error=error)
         CALL cp_dbcsr_set(xas_env%ostrength_sm(1)%matrix,0.0_dp,error=error)
         DO i = 2,3
           ALLOCATE(xas_env%ostrength_sm(i)%matrix)
           CALL cp_dbcsr_init(xas_env%ostrength_sm(i)%matrix,error=error)
           CALL cp_dbcsr_copy(xas_env%ostrength_sm(i)%matrix,xas_env%ostrength_sm(1)%matrix,&
                "xas_env%ostrength_sm-"//TRIM(ADJUSTL(cp_to_string(i))),error=error)
           CALL cp_dbcsr_set(xas_env%ostrength_sm(i)%matrix,0.0_dp,error=error)
         END DO

         CALL array_release (row_blk_sizes)
      END IF

   !Define the qs_loc_env : to find centers, spread and possibly localize them
      IF(.NOT.(ASSOCIATED(xas_env%qs_loc_env))) THEN
        CALL qs_loc_env_create(qs_loc_env,error=error)
        CALL set_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env,error=error)
        CALL qs_loc_env_release(qs_loc_env,error=error)
        CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env,error=error)
        loc_section => section_vals_get_subs_vals(xas_section,"LOCALIZE",error=error)
        CALL qs_loc_control_init(qs_loc_env,loc_section,do_homo=.TRUE.,&
             do_xas=.TRUE.,nloc_xas=nexc_search,error=error)
        IF(.NOT. qs_loc_env%do_localize) THEN
          qs_loc_env%localized_wfn_control%localization_method = do_loc_none
        ELSE
          nmoloc = qs_loc_env%localized_wfn_control%nloc_states
          CALL set_loc_wfn_lists(qs_loc_env%localized_wfn_control,nmoloc,n_mo,nspins,error=error)
          CALL set_loc_centers(qs_loc_env%localized_wfn_control,nmoloc,nspins,error=error)
          CALL qs_loc_env_init(qs_loc_env,qs_loc_env%localized_wfn_control,&
               qs_env,myspin=1,do_localize=qs_loc_env%do_localize,error=error)
        END IF
      END IF

    ! Initialize Mixing scheme
      mixing_section => section_vals_get_subs_vals(xas_section,"MIXING",error=error)
      CALL section_vals_val_get(mixing_section,"METHOD",i_val=xas_env%mixing_method,error=error)
      CALL mixing_storage_create(xas_env%mixing_store, mixing_section, xas_env%mixing_method, &
              dft_control%qs_control%cutoff, error=error)

      CALL init_smear(xas_env%smear,error)
      smear_section => section_vals_get_subs_vals(xas_section,"SMEAR",error=error)
      CALL read_smear_section(xas_env%smear,smear_section,error=error)

      !Type of state
      norb = 0
      ALLOCATE(nq(1),lq(1),sto_zet(1),STAT=istat)
      IF( xas_control%state_type == xas_1s_type) THEN
        norb(0) = 1
        nq(1) = 1
        lq(1) = 0
      ELSEIF( xas_control%state_type == xas_2s_type ) THEN
        norb(0) = 2
        nq(1) = 2
        lq(1) = 0
      ELSEIF( xas_control%state_type == xas_2p_type ) THEN
        norb(0) = 2
        norb(1) = 1
        nq(1) = 2
        lq(1) = 1
      ELSE
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "XAS type of state not implemented")
      END IF

      ALLOCATE(kind_type_tmp(nkind),STAT=istat)
      ALLOCATE(kind_z_tmp(nkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      kind_type_tmp = 0
      kind_z_tmp = 0
      nk=0
      DO iat = 1,nexc_atoms
        iatom = xas_env%exc_atoms(iat)
        NULLIFY(atomic_kind)
        atomic_kind =>  particle_set(iatom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
             kind_number=ikind,zeff=zatom)
        ihavethis = .FALSE.
        DO ik = 1,nk
          IF(ikind==kind_type_tmp(ik)) THEN
            ihavethis = .TRUE.
            xas_env%mykind_of_atom(iat) = ik
            EXIT
          END IF
        END DO
        IF(.NOT. ihavethis) THEN
          nk = nk +1
          kind_type_tmp(nk) = ikind
          kind_z_tmp(nk) = INT(zatom)
          xas_env%mykind_of_atom(iat) = nk
          xas_env%mykind_of_kind(ikind) = nk
        END IF
      END DO  ! iat

      ALLOCATE(xas_env%my_gto_basis(nk),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%stogto_overlap(nk),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO ik = 1,nk
        NULLIFY(xas_env%my_gto_basis(ik)%gto_basis_set,sto_basis_set)
        ne = 0
        DO l=1,4 !lq(1)+1
          nj = 2*(l-1)+1
          DO i=l, 7! nq(1)
            ne(l,i) = ptable(kind_z_tmp(ik))%e_conv(l-1) - 2*nj*(i-l)
            ne(l,i) = MAX(ne(l,i),0)
            ne(l,i) = MIN(ne(l,i),2*nj)
          END DO
        END DO

        sto_zet(1) = srules(kind_z_tmp(ik),ne,nq(1),lq(1))
        CALL allocate_sto_basis_set(sto_basis_set,error)
        name_sto='xas_tmp_sto'
        CALL set_sto_basis_set(sto_basis_set,nshell=1,nq=nq,&
             lq=lq,zet=sto_zet,name=name_sto)
        CALL create_gto_from_sto_basis(sto_basis_set,&
             xas_env%my_gto_basis(ik)%gto_basis_set,xas_control%ngauss,error=error)
        CALL deallocate_sto_basis_set(sto_basis_set,error)
        xas_env%my_gto_basis(ik)%gto_basis_set%norm_type = 2
        CALL init_orb_basis_set(xas_env%my_gto_basis(ik)%gto_basis_set,error=error)

        atomic_kind => atomic_kind_set(kind_type_tmp(ik))
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
             orb_basis_set=orb_basis_set)

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nsgf=nsgf_gto)
        CALL get_gto_basis_set(gto_basis_set=xas_env%my_gto_basis(ik)%gto_basis_set,nsgf=nsgf_sto)
        ALLOCATE(xas_env%stogto_overlap(ik)%array(nsgf_sto,nsgf_gto),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        CALL calc_stogto_overlap(xas_env%my_gto_basis(ik)%gto_basis_set,orb_basis_set,&
             xas_env%stogto_overlap(ik)%array,error=error)
      END DO

      DEALLOCATE(nq,lq,sto_zet,STAT=istat)
      DEALLOCATE(kind_type_tmp,kind_z_tmp,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF  ! failure

  END SUBROUTINE xas_env_init

! *****************************************************************************
!> \brief Calculate and write the spectrum relative to the core level excitation
!>      of a specific atom. It works for TP approach, because of the definition
!>      of the oscillator strengths as  matrix elements of the dipole operator
!> \param xas_control ...
!> \param xas_env ...
!> \param qs_env ...
!> \param xas_section ...
!> \param iatom index of the excited atom
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param iter iteration of the xas calculation UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      03.2006 created [MI]
!> \author MI
!> \note
!>      for the tddft calculation should be re-thought
! *****************************************************************************
  SUBROUTINE cls_calculate_spectrum(xas_control,xas_env,qs_env,xas_section,&
             iatom,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xas_section
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cls_calculate_spectrum', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: homo, i, istat, lfomo, nabs, &
                                                nmo, nvirtual, output_unit, &
                                                xas_estate
    LOGICAL                                  :: append_cube, failure, length
    REAL(dp)                                 :: rc(3)
    REAL(dp), DIMENSION(:), POINTER          :: all_evals
    REAL(dp), DIMENSION(:, :), POINTER       :: sp_ab, sp_em
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op_sm, ostrength_sm
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dip_fm_set
    TYPE(cp_fm_type), POINTER                :: all_vectors, excvec_coeff, &
                                                fm_work
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure =.FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    NULLIFY(fm_work,ostrength_sm,op_sm, dip_fm_set)
    NULLIFY(all_evals,all_vectors,excvec_coeff)
    NULLIFY(mos,particle_set,sp_em,sp_ab)
    ALLOCATE (op_sm(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env, &
         mos=mos,particle_set=particle_set,error=error)

    CALL get_mo_set(mos(1)%mo_set, homo=homo,lfomo=lfomo,nmo=nmo)
    CALL get_xas_env(xas_env=xas_env,all_vectors=all_vectors,xas_estate=xas_estate,&
         all_evals=all_evals,dip_fm_set=dip_fm_set,excvec_coeff=excvec_coeff,&
         fm_work=fm_work, ostrength_sm=ostrength_sm,nvirtual=nvirtual,error=error)

    nabs = nvirtual -lfomo + 1
    ALLOCATE (sp_em(6,homo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (sp_ab(6,nabs),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(excvec_coeff),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

      IF(.NOT. xas_control%xas_method == xas_dscf) THEN
      ! Calculate the spectrum
        IF(xas_control%dipole_form==xas_dip_len2) THEN
           rc(1:3) = particle_set(iatom)%r(1:3)
           DO i = 1,3
             NULLIFY(op_sm(i)%matrix)
             op_sm(i)%matrix => ostrength_sm(i)%matrix
           END DO
           CALL rRc_xyz_ao(op_sm,qs_env,rc,order=1,minimum_image=.TRUE.,error=error)
           CALL spectrum_dip_vel(dip_fm_set,op_sm,mos,excvec_coeff,&
                all_vectors,all_evals,fm_work,&
                sp_em,sp_ab,xas_estate,nvirtual,error=error)
           DO i = 1,SIZE(ostrength_sm,1)
              CALL cp_dbcsr_set(ostrength_sm(i)%matrix,0.0_dp,error=error)
           END DO
        ELSE
           DO i = 1,3
             NULLIFY(op_sm(i)%matrix)
             op_sm(i)%matrix => ostrength_sm(i)%matrix
           END DO
           CALL spectrum_dip_vel(dip_fm_set,op_sm,mos,excvec_coeff,&
                all_vectors,all_evals,fm_work,&
                sp_em,sp_ab,xas_estate,nvirtual,error=error)
        END IF
      END IF

      CALL get_mo_set(mos(1)%mo_set, lfomo=lfomo)
      ! writw thw spectrum, if the file exists it is appended
      IF( .NOT. xas_control%xas_method == xas_dscf) THEN
        length = (.NOT. xas_control%dipole_form==xas_dip_vel)
        CALL xas_write(sp_em,sp_ab, xas_estate, &
             xas_section, iatom, lfomo, length=length, error=error)
      END IF

      DEALLOCATE(sp_em,STAT=istat)
      DEALLOCATE(sp_ab,STAT=istat)

      IF(BTEST(cp_print_key_should_output(logger%iter_info,xas_section,&
        "PRINT%CLS_FUNCTION_CUBES",error=error),cp_p_file)) THEN
        append_cube= section_get_lval(xas_section,"PRINT%CLS_FUNCTION_CUBES%APPEND",error=error)
        CALL xas_print_cubes(xas_control,qs_env,xas_section,mos,all_vectors,&
             iatom,append_cube,error=error)
      END IF

    END IF  ! failure

    DEALLOCATE (op_sm,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE cls_calculate_spectrum

! *****************************************************************************
!> \brief write the spectrum for each atom in a different output file
!> \param sp_em ...
!> \param sp_ab ...
!> \param estate ...
!> \param xas_section ...
!> \param iatom index of the excited atom
!> \param lfomo ...
!> \param length ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param ionode logical assigning the i/o node UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param spectrum temporary variable contaning the spectrum for the actual atom UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param xas_rootfname root of the filename where the spectrum is written UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      05.2005 created [MI]
!> \author MI
!> \note
!>      the iteration counter is not finilized yet
! *****************************************************************************
  SUBROUTINE xas_write(sp_em, sp_ab, estate, xas_section, iatom, &
             lfomo, length, error)

    REAL(dp), DIMENSION(:, :), POINTER       :: sp_em, sp_ab
    INTEGER, INTENT(IN)                      :: estate
    TYPE(section_vals_type), POINTER         :: xas_section
    INTEGER, INTENT(IN)                      :: iatom, lfomo
    LOGICAL, INTENT(IN)                      :: length
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: mittle_ab, mittle_em, my_act, &
                                                my_pos
    INTEGER                                  :: i, istate, out_sp_ab, &
                                                out_sp_em
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ene2
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    failure = .FALSE.

    my_pos = "APPEND"
    my_act = "WRITE"

    mittle_em = "xes_at"//TRIM(ADJUSTL(cp_to_string(iatom)))
    out_sp_em = cp_print_key_unit_nr(logger,xas_section,"PRINT%XES_SPECTRUM",&
               extension=".spectrum", file_position=my_pos, file_action=my_act,&
               file_form="FORMATTED", middle_name=TRIM(mittle_em), &
               error=error)

    IF(out_sp_em>0) THEN
      WRITE(out_sp_em,'(A,I6,A,I6,A,I6)') " Emission spectrum for atom ", iatom,&
             ", index of excited core MO is", estate, ", # of lines ", SIZE(sp_em,2)
      ene2=1.0_dp
      DO istate = estate,SIZE(sp_em,2)
        IF(length) ene2 = sp_em(1,istate)*sp_em(1,istate)
         WRITE(out_sp_em,'(I6,5F16.8,F10.5)') istate, sp_em(1,istate)*evolt, &
                 sp_em(2,istate)*ene2, sp_em(3,istate)*ene2,&
                 sp_em(4,istate)*ene2, sp_em(5,istate)*ene2,sp_em(6,istate)
      END DO
    END IF
    CALL cp_print_key_finished_output(out_sp_em,logger,xas_section,&
                "PRINT%XES_SPECTRUM", error=error)

    mittle_ab = "xas_at"//TRIM(ADJUSTL(cp_to_string(iatom)))
    out_sp_ab = cp_print_key_unit_nr(logger,xas_section,"PRINT%XAS_SPECTRUM",&
               extension=".spectrum", file_position=my_pos, file_action=my_act,&
               file_form="FORMATTED", middle_name=TRIM(mittle_ab), &
               error=error)

    IF(out_sp_ab>0) THEN
      WRITE(out_sp_ab,'(A,I6,A,I6,A,I6)') " Absorption spectrum for atom ", iatom,&
             ", index of excited core MO is", estate, ", # of lines ", SIZE(sp_ab,2)
      ene2=1.0_dp
      DO i = 1,SIZE(sp_ab,2)
        istate = lfomo -1 + i
        IF(length) ene2 = sp_ab(1,i)*sp_ab(1,i)
         WRITE(out_sp_ab,'(I6,5F16.8,F10.5)') istate, sp_ab(1,i)*evolt, &
                 sp_ab(2,i)*ene2, sp_ab(3,i)*ene2,&
                 sp_ab(4,i)*ene2, sp_ab(5,i)*ene2,sp_ab(6,i)
      END DO
    END IF

    CALL cp_print_key_finished_output(out_sp_ab,logger,xas_section,&
                "PRINT%XAS_SPECTRUM", error=error)

  END SUBROUTINE xas_write

! *****************************************************************************
!> \brief write the cube files for a set of selected states
!> \param xas_control provide number ant indexes of the states to be printed
!> \param qs_env ...
!> \param xas_section ...
!> \param mos mos from which the states to be printed are extracted
!> \param all_vectors ...
!> \param iatom index of the atom that has been excited
!> \param append_cube ...
!> \param error ...
!> \param ionod e UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param iter iteration step (e.g. along an md) UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \param qs_en v UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      08.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE  xas_print_cubes(xas_control,qs_env,xas_section,&
              mos,all_vectors,iatom,append_cube, error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xas_section
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_type), POINTER                :: all_vectors
    INTEGER, INTENT(IN)                      :: iatom
    LOGICAL, INTENT(IN)                      :: append_cube
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_print_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_mittle, my_pos
    INTEGER                                  :: homo, istat, istate0, nspins, &
                                                nstates
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: centers
    TYPE(section_vals_type), POINTER         :: print_key

    failure = .FALSE.
    nspins = SIZE(mos)

    print_key => section_vals_get_subs_vals(xas_section,"PRINT%CLS_FUNCTION_CUBES",error=error)
    my_mittle = 'at'//TRIM(ADJUSTL(cp_to_string(iatom)))
    nstates = SIZE(xas_control%list_cubes,1)

    IF(xas_control%do_centers) THEN
    ! one might like to calculate the centers of the xas orbital (without localizing them)
    ELSE
      ALLOCATE(centers(6,nstates),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      centers = 0.0_dp
    END IF

    CALL get_mo_set(mos(1)%mo_set, homo=homo)
    istate0 = 0

    my_pos = "REWIND"
    IF (append_cube) THEN
      my_pos = "APPEND"
    END IF

    CALL qs_print_cubes(qs_env,all_vectors,nstates,xas_control%list_cubes,&
         centers,print_key,my_mittle,state0=istate0,file_position=my_pos,error=error)

    DEALLOCATE(centers,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE xas_print_cubes

! *****************************************************************************
!> \brief Calculation of the spectrum when the dipole approximation
!>      in the length form is used. The position operator is calculated
!>      by the berry phase approach
!> \param fm_set components of the position operator in a full matrix form
!>                already multiplied by the coefficiets
!>                only the terms <C_i Op C_f> are calculated where
!>                C_i are the coefficients of the excited state
!> \param op_sm components of the position operator for the dipole
!>               in a sparse matrix form (cos and sin)
!>               calculated for the basis functions
!> \param mos wavefunctions coefficients
!> \param excvec coefficients of the excited orbital
!> \param all_vectors ...
!> \param all_evals ...
!> \param fm_work work space
!> \param cell parameters for the simulation cell
!> \param sp_em ...
!> \param sp_ab ...
!> \param estate index of the excited state
!> \param nstate ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param spectrum temporary variable contaning the spectrum for the actual atom UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      06.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE spectrum_dip_len(fm_set,op_sm,mos,excvec,&
             all_vectors,all_evals,fm_work,cell, sp_em, sp_ab,estate, nstate,error)

    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: fm_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: op_sm
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_type), POINTER                :: excvec, all_vectors
    REAL(dp), DIMENSION(:), POINTER          :: all_evals
    TYPE(cp_fm_type), POINTER                :: fm_work
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), POINTER       :: sp_em, sp_ab
    INTEGER, INTENT(IN)                      :: estate, nstate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'spectrum_dip_len', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: z
    INTEGER                                  :: homo, i, i_abs, istate, j, &
                                                lfomo, nao, nmo
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ene_f, ene_i, imagpart, &
                                                ra(3), realpart
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, &
                                                occupation_numbers

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(fm_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)

    NULLIFY(eigenvalues, occupation_numbers)
    IF(.NOT. failure) THEN
      CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues, occupation_numbers=occupation_numbers,&
           nao=nao, nmo=nmo, homo =homo, lfomo=lfomo)
      DO i=1,SIZE(fm_set,2)
        DO j = 1,SIZE(fm_set,1)
           CPPrecondition(ASSOCIATED(fm_set(j,i)%matrix),cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(fm_set(j,i)%matrix, 0.0_dp, error=error)
           CALL cp_fm_set_all(fm_work, 0.0_dp, error=error)
           CALL cp_dbcsr_sm_fm_multiply(op_sm(j,i)%matrix,all_vectors,fm_work,ncol=nstate,error=error)
           CALL cp_gemm("T","N",1,nstate,nao,1.0_dp,excvec,&
                  fm_work,0.0_dp, fm_set(j,i)%matrix,b_first_col=1,error=error)
        END DO
      END DO

      sp_em = 0.0_dp
      sp_ab = 0.0_dp
      ene_i = eigenvalues(estate)
      DO istate = 1,nstate
         ene_f = all_evals(istate)
         DO i = 1,3
           CALL cp_fm_get_element(fm_set(1,i)%matrix,1,istate,realpart)
           CALL cp_fm_get_element(fm_set(2,i)%matrix,1,istate,imagpart)
           z = CMPLX(realpart,imagpart,dp)
           ra(i) = ( cell % hmat ( i, i ) / twopi ) * AIMAG ( LOG ( z ) )
         END DO
         IF(istate<=homo) THEN
           sp_em(1,istate) = ene_f - ene_i
           sp_em(2,istate) = ra(1)*ra(1)
           sp_em(3,istate) = ra(2)*ra(2)
           sp_em(4,istate) = ra(3)*ra(3)
           sp_em(5,istate) = ra(1)*ra(1)+ra(2)*ra(2)+ra(3)*ra(3)
           sp_em(6,istate) = occupation_numbers(istate)
         END IF
         IF(istate>=lfomo ) THEN
           i_abs = istate - lfomo  + 1
           sp_ab(1,i_abs) = ene_f - ene_i
           sp_ab(2,i_abs) = ra(1)*ra(1)
           sp_ab(3,i_abs) = ra(2)*ra(2)
           sp_ab(4,i_abs) = ra(3)*ra(3)
           sp_ab(5,i_abs) = ra(1)*ra(1)+ra(2)*ra(2)+ra(3)*ra(3)
           IF(istate<=nmo) sp_ab(6,i_abs) = occupation_numbers(istate)
         END IF
      END DO
    END IF

  END SUBROUTINE spectrum_dip_len

! *****************************************************************************
!> \brief Calculation of the spectrum when the dipole approximation
!>      in the velocity form is used.
!> \param fm_set components of the position operator in a full matrix form
!>                already multiplied by the coefficiets
!>                only the terms <C_i Op C_f> are calculated where
!>                C_i are the coefficients of the excited state
!> \param op_sm components of the position operator for the dipole
!>               in a sparse matrix form (cos and sin)
!>               calculated for the basis functions
!> \param mos wavefunctions coefficients
!> \param excvec coefficients of the excited orbital
!> \param all_vectors ...
!> \param all_evals ...
!> \param fm_work work space
!> \param sp_em ...
!> \param sp_ab ...
!> \param estate index of the excited state
!> \param nstate ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param spectrum temporary variable contaning the spectrum for the actual atom UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      06.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE spectrum_dip_vel(fm_set,op_sm,mos,excvec, &
             all_vectors,all_evals,fm_work,sp_em, sp_ab,estate, nstate,error)

    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: fm_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op_sm
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_type), POINTER                :: excvec, all_vectors
    REAL(dp), DIMENSION(:), POINTER          :: all_evals
    TYPE(cp_fm_type), POINTER                :: fm_work
    REAL(dp), DIMENSION(:, :), POINTER       :: sp_em, sp_ab
    INTEGER, INTENT(IN)                      :: estate, nstate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'spectrum_dip_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: homo, i, i_abs, istate, &
                                                lfomo, nao, nmo
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dip(3), ene_f, ene_i
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, &
                                                occupation_numbers

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(fm_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    NULLIFY(eigenvalues, occupation_numbers)

    IF(.NOT. failure) THEN
      CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues,occupation_numbers=occupation_numbers,&
           nao=nao, nmo=nmo, homo =homo, lfomo=lfomo)

      DO i=1,SIZE(fm_set,2)
         CPPrecondition(ASSOCIATED(fm_set(1,i)%matrix),cp_failure_level,routineP,error,failure)
         CALL cp_fm_set_all(fm_set(1,i)%matrix, 0.0_dp, error=error)
         CALL cp_fm_set_all(fm_work, 0.0_dp, error=error)
         CALL cp_dbcsr_sm_fm_multiply(op_sm(i)%matrix,all_vectors,fm_work,ncol=nstate,error=error)
         CALL cp_gemm("T","N",1,nstate,nao,1.0_dp,excvec,&
              fm_work,0.0_dp, fm_set(1,i)%matrix,b_first_col=1,error=error)
      END DO

      sp_em = 0.0_dp
      sp_ab = 0.0_dp
      ene_i = eigenvalues(estate)

      DO istate = 1,nstate
         ene_f = all_evals(istate)

         DO i = 1,3
           CALL cp_fm_get_element(fm_set(1,i)%matrix,1,istate,dip(i))
         END DO
         IF(istate<=homo) THEN
           sp_em(1,istate) = ene_f - ene_i
           sp_em(2,istate) = dip(1)
           sp_em(3,istate) = dip(2)
           sp_em(4,istate) = dip(3)
           sp_em(5,istate) = dip(1)*dip(1)+dip(2)*dip(2)+dip(3)*dip(3)
           sp_em(6,istate) = occupation_numbers(istate)
         END IF
         IF(istate>=lfomo ) THEN
           i_abs = istate - lfomo + 1
           sp_ab(1,i_abs) = ene_f - ene_i
           sp_ab(2,i_abs) = dip(1)
           sp_ab(3,i_abs) = dip(2)
           sp_ab(4,i_abs) = dip(3)
           sp_ab(5,i_abs) = dip(1)*dip(1)+dip(2)*dip(2)+dip(3)*dip(3)
           IF(istate<=nmo) sp_ab(6,i_abs) = occupation_numbers(istate)
         END IF

      END DO
    END IF

  END SUBROUTINE spectrum_dip_vel

! *****************************************************************************
!> \brief ...
!> \param base_a ...
!> \param base_b ...
!> \param matrix ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calc_stogto_overlap(base_a, base_b, matrix,error)

    TYPE(gto_basis_set_type), POINTER        :: base_a, base_b
    REAL(dp), DIMENSION(:, :), POINTER       :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_stogto_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER :: iset, istat, jset, ldsab, maxcoa, maxcob, maxl, maxla, maxlb, &
      na, nb, nseta, nsetb, nsgfa, nsgfb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, &
                                                nsgfa_set, nsgfb_set
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rab(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: sab, work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, scon_a, scon_b, &
                                                sphi_a, sphi_b, zeta, zetb

    failure = .FALSE.

    NULLIFY(la_max,la_min,lb_max,lb_min)
    NULLIFY(npgfa,npgfb,nsgfa_set,nsgfb_set)
    NULLIFY(first_sgfa,first_sgfb)
    NULLIFY(rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb)

    CALL get_gto_basis_set(gto_basis_set=base_a,nsgf=nsgfa,nsgf_set=nsgfa_set,lmax=la_max,&
                          lmin=la_min,npgf=npgfa,pgf_radius=rpgfa, &
                          sphi=sphi_a,scon=scon_a,zet=zeta,first_sgf=first_sgfa, &
                          maxco=maxcoa,nset=nseta, maxl=maxla)

    CALL get_gto_basis_set(gto_basis_set=base_b,nsgf=nsgfb,nsgf_set=nsgfb_set, lmax=lb_max,&
                          lmin=lb_min,npgf=npgfb,pgf_radius=rpgfb, &
                          sphi=sphi_b,scon=scon_b,zet=zetb,first_sgf=first_sgfb, &
                          maxco=maxcob,nset=nsetb,maxl=maxlb)
    ! Initialize and allocate
    rab = 0.0_dp
    matrix = 0.0_dp

    ldsab = MAX(maxcoa,maxcob,nsgfa,nsgfb)
    maxl = MAX(maxla,maxlb)

    ALLOCATE(sab(ldsab,ldsab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work(ldsab,ldsab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO iset=1,nseta

       na = npgfa(iset)*(ncoset(la_max(iset))-ncoset(la_min(iset)-1))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          nb = npgfb(jset)*(ncoset(lb_max(jset))-ncoset(lb_min(iset)-1))
          sgfb = first_sgfb(1,jset)

          CALL overlap_ab(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                          rab,sab,error=error)
          CALL contraction(sab,work,ca=scon_a(:,sgfa:),na=na,ma=nsgfa_set(iset),&
                                    cb=scon_b(:,sgfb:),nb=nb,mb=nsgfb_set(jset),error=error)
          CALL block_add("IN",work,nsgfa_set(iset),nsgfb_set(jset),matrix,sgfa,sgfb,error=error)

       END DO  ! jset
    END DO  ! iset

    DEALLOCATE(sab,work,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE calc_stogto_overlap

! *****************************************************************************
!> \brief Starting from a set of mos, determine on which atom are centered
!>      and if they are of the right type (1s,2s ...)
!>      to be used in the specific core level spectrum calculation
!>      The set of states need to be from the core, otherwise the
!>      characterization of the type is not valid, since it assumes that
!>      the orbital is localizad on a specific atom
!>      It is probably reccomandable to run a localization cycle before
!>      proceeding to the assignment of the type
!>      The type is determined by computing the overalp with a
!>      type specific, minimal, STO bais set
!> \param xas_control ...
!> \param xas_env ...
!> \param localized_wfn_control ...
!> \param qs_env ...
!> \param error ...
!> \par History
!>      03.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE cls_assign_core_states(xas_control,xas_env,localized_wfn_control,qs_env,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cls_assign_core_states', &
      routineP = moduleN//':'//routineN

    INTEGER :: homo, i, iat, iatom, ikind, isgf, istat, istate, j, my_kind, &
      nao, natom, nexc_atoms, nexc_search, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
    INTEGER, DIMENSION(3)                    :: perd0
    INTEGER, DIMENSION(:), POINTER           :: atom_of_state, &
                                                mykind_of_kind, &
                                                state_of_atom, &
                                                state_of_mytype, type_of_state
    LOGICAL                                  :: failure
    REAL(dp)                                 :: component, dist, distmin, &
                                                maxocc, ra(3), rac(3), rc(3)
    REAL(dp), DIMENSION(:), POINTER          :: max_overlap, sto_state_overlap
    REAL(dp), DIMENSION(:, :), POINTER       :: centers_wfn
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: stogto_overlap
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    NULLIFY(cell,mos,particle_set)
    NULLIFY(atom_of_state,centers_wfn,mykind_of_kind,state_of_atom)
    NULLIFY(stogto_overlap,type_of_state,max_overlap)
    NULLIFY(state_of_mytype,type_of_state,sto_state_overlap)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CALL get_qs_env(qs_env=qs_env, cell=cell, mos=mos, particle_set=particle_set,error=error)

    ! The Berry operator can be used only for periodic systems
    ! If an isolated system is used the periodicity is overimposed
    perd0(1:3) = cell%perd(1:3)
    cell%perd(1:3) = 1

    CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, maxocc=maxocc, nao=nao, homo=homo)

    CALL get_xas_env(xas_env=xas_env,&
         centers_wfn=centers_wfn,atom_of_state=atom_of_state,&
         mykind_of_kind=mykind_of_kind,&
         type_of_state=type_of_state, state_of_atom=state_of_atom,&
         stogto_overlap=stogto_overlap,nexc_atoms=nexc_atoms, nexc_search=nexc_search,error=error)

    ! scratch array for the state
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    natom = SIZE(particle_set)

    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set=particle_set, first_sgf=first_sgf,error=error)
    ALLOCATE (sto_state_overlap(nexc_search),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (max_overlap(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    max_overlap = 0.0_dp
    ALLOCATE (state_of_mytype(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    state_of_mytype = 0
    atom_of_state = 0

    DO istate = 1,nexc_search
      centers_wfn(1,istate) =  localized_wfn_control%centers_set(1)%array(1,istate)
      centers_wfn(2,istate) =  localized_wfn_control%centers_set(1)%array(2,istate)
      centers_wfn(3,istate) =  localized_wfn_control%centers_set(1)%array(3,istate)

      ! Assign the state to the closest atom
      distmin = 100.0_dp

      DO iat = 1,nexc_atoms
        iatom = xas_control%exc_atoms(iat)
        ra(1:3) = particle_set(iatom)%r(1:3)
        rc(1:3) = centers_wfn(1:3,istate)
        rac = pbc(ra,rc,cell)
        dist = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)

        IF(dist < distmin) THEN
           atom_of_state(istate) = iatom
           distmin = dist
        END IF
      END DO
      IF(atom_of_state(istate) /= 0) THEN
      !Character of the state
         CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,istate,&
              nao,1,transpose=.TRUE.,error=error)

         iatom = atom_of_state(istate)

         NULLIFY(atomic_kind)
         atomic_kind =>  particle_set(iatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
                kind_number=ikind)

         my_kind = mykind_of_kind(ikind)

         sto_state_overlap(istate) = 0.0_dp
         DO i = 1,SIZE(stogto_overlap(my_kind)%array,1)
           component = 0.0_dp
           DO j = 1,SIZE(stogto_overlap(my_kind)%array,2)
             isgf = first_sgf(iatom) + j - 1
             component = component + stogto_overlap(my_kind)%array(i,j)*vecbuffer(1,isgf)
           END DO
           sto_state_overlap(istate) = sto_state_overlap(istate) + &
                   component * component
         END DO
         IF(sto_state_overlap(istate)>max_overlap(iatom)) THEN
           state_of_mytype(iatom) = istate
           max_overlap(iatom) = sto_state_overlap(istate)
         END IF
      END IF
    END DO  ! istate

  ! In the set of states, assign the index of the state to be excited for iatom
    IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
           "List the atoms to be excited and the relative of MOs index "
    END IF
    DO iat = 1,nexc_atoms
      iatom = xas_env%exc_atoms(iat)
      state_of_atom(iat) = state_of_mytype(iatom)
      IF (output_unit>0) THEN
         WRITE(UNIT=output_unit,FMT="(T10,A,I6,T32,A,I6)") 'Atom: ',iatom ,&
         "MO index", state_of_atom(iat)
      END IF
      IF(state_of_atom(iat)==0 .OR. state_of_atom(iat) .GT. homo) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
        "A wrong state has been selected for excitation, check the Wannier centers")
      END IF
    END DO

    ! Set back the correct periodicity
    cell%perd(1:3) =  perd0(1:3)

    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sto_state_overlap,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (max_overlap,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (state_of_mytype,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE cls_assign_core_states

END MODULE xas_methods
