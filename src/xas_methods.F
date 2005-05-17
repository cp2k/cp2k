!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/xas_methods   [1.0] *
!!
!!   NAME
!!      xas_methods
!!
!!   FUNCTION
!!      driver for the xas calculation and xas_scf for the tp method
!!
!!   AUTHOR
!!     MI (05.2005)
!!
!!   MODIFICATION HISTORY
!!     created 05.2005
!!
!!   SOURCE
!******************************************************************************


MODULE xas_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: xas_none, xas_tp_hh, xas_tp_fh, xas_tddft,&
                                             xas_1s_type, xas_2s_type,&
                                             do_loc_none, state_loc_upton
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_loc_control,                  ONLY: localized_wfn_control_type,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_retain
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_types,                    ONLY: qs_loc_env_release,&
                                             qs_loc_env_retain,&
                                             qs_loc_env_create,&
                                             qs_loc_env_new_type
  USE qs_loc_utils,                    ONLY: qs_loc_env_init
  USE qs_matrix_pools,                 ONLY: mpools_get, &
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,mo_set_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type, &
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE xas_control,                     ONLY: xas_control_type
  USE xas_env_types,                   ONLY: get_xas_env, set_xas_env, &  
                                             xas_env_create, xas_env_release,&
                                             xas_env_retain, &
                                             xas_environment_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: module_name='xas_methods'  

! *** Public subroutines ***

  PUBLIC :: xas 
 
!****************************************************************************

CONTAINS

!!****f* xas_methods/xas [1.0] *
!!
!!   NAME
!!     xas
!!
!!   FUNCTION
!!     Driver for xas calculations
!!     The initial mos are prepared
!!     A loop on the atoms to be excited is started
!!     For each atom the state to be excited is identified 
!!     An scf optimization using the TP scheme or TD-DFT is used
!!     to evaluate the spectral energies and oscillator strengths
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env, the xas_env lives in
!!     - globenv: global environment for i/o and error
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE xas(qs_env, xas_control, globenv, error)

    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(xas_control_type)                      :: xas_control
    TYPE(global_environment_type), &
      POINTER                                   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "xas", &
      routineP = module_name//':'//routineN

    INTEGER :: handle, homo, iat, iatom, ispin, istat, istate, iter, &
               nao, natoms, nmo, nspins, nexc_atoms, nexc_search, output_unit,&
               state_to_be_excited, xas_iter_count
    INTEGER, DIMENSION(:), POINTER              :: atom_of_state
    INTEGER, DIMENSION(:), POINTER              :: state_of_atom
    LOGICAL :: failure, ionode, transition_potential

    REAL(dp), DIMENSION(:,:), POINTER           :: spectrum
    REAL(dp), DIMENSION(:,:), POINTER           :: centers_wfn
    REAL(dp), DIMENSION(:,:), POINTER           :: vecbuffer
    TYPE( cp_fm_p_type ), DIMENSION(:),  &
      POINTER                                   :: groundstate_coeff
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(cp_fm_type), POINTER                   :: EXCVEC_COEFF
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(cp_para_env_type), POINTER             :: para_env
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(localized_wfn_control_type), POINTER   :: my_localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                   :: particle_set
    TYPE( qs_loc_env_new_type ), POINTER        :: qs_loc_env
    TYPE(qs_scf_env_type), POINTER              :: scf_env
    TYPE( real_matrix_p_type ), DIMENSION(:), &
      POINTER                                   :: ostrength_sm
    TYPE(xas_environment_type), POINTER         :: xas_env


    CALL timeset(routineN,"I",'',handle)

    failure = .FALSE.
    transition_potential = .FALSE.


    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)


    NULLIFY(xas_env, groundstate_coeff, ostrength_sm)
    NULLIFY(spectrum, state_of_atom, centers_wfn, qs_loc_env)

    CALL get_qs_env(qs_env,xas_env=xas_env,error=error)
    IF (.NOT.ASSOCIATED(xas_env)) THEN
      CALL xas_env_create(xas_env, error=error)
      CALL xas_env_init(xas_env, xas_control, qs_env, error=error)
      CALL set_qs_env(qs_env,xas_env=xas_env, error=error)
      CALL xas_env_release(xas_env,error=error)
      CALL get_qs_env(qs_env,xas_env=xas_env,error=error)
    END IF

    NULLIFY(dft_control, mos, para_env, particle_set)
    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
          mos=mos, para_env=para_env, particle_set=particle_set )

    nexc_search = xas_control%nexc_search
    nspins = dft_control%nspins
    transition_potential = (xas_control%xas_method==xas_tp_hh).OR.&
                           (xas_control%xas_method==xas_tp_fh)
    IF(nspins ==1 .AND. transition_potential) THEN
       CALL stop_program(routineP,"xas with tp method requires LSD calculations")
    END IF

    CALL get_mo_set(mos(ispin)%mo_set,nao=nao,homo=homo)
    IF(nexc_search < 0) nexc_search = homo

    CALL get_xas_env(xas_env=xas_env,iter_count=xas_iter_count,&
         centers_wfn=centers_wfn,atom_of_state=atom_of_state,&
         state_of_atom=state_of_atom,spectrum=spectrum,&
         groundstate_coeff=groundstate_coeff,&
         ostrength_sm=ostrength_sm)

    xas_iter_count =  xas_iter_count + 1
    CALL set_xas_env(xas_env=xas_env,iter_count=xas_iter_count)

    ! Take from xas_control what is needed
!     indexes_of_atoms and filenames for the output of spectra
!     exc_orb_type
!     excene and oscillators (size = nmo - homo = all virtual)
!     os_sm : sparse matrix of the oscillator strengths (for tp is needed for sure)
  
   !I know from which atoms I want to excite (O) and which type of states (1s) (input)

  !Consider the set of states among which I look for the one i want to excite (input or not)

   !Calculate os_sm 

   !Define the qs_loc_env : to find centers, spread and possibly localize them
    my_localized_wfn_control =>  xas_control%localized_wfn_control
    CALL localized_wfn_control_retain(my_localized_wfn_control,error=error)
    CPPrecondition(ASSOCIATED(my_localized_wfn_control),cp_failure_level,routineP,error,failure)
    IF(.NOT. xas_control%localized_wfn) THEN
      my_localized_wfn_control%localization_method = do_loc_none
    END IF
       ! The localization is needed only for the states up to nexc_search
    my_localized_wfn_control%set_of_states = state_loc_upton
    my_localized_wfn_control%nloc_states(1) = nexc_search
    my_localized_wfn_control%lu_bound_states(1,1) = 1
    my_localized_wfn_control%lu_bound_states(2,1) = nexc_search

    CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
    IF(.NOT.(ASSOCIATED(qs_loc_env))) THEN
      CALL qs_loc_env_create(qs_loc_env,error=error)
      CALL set_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
      CALL qs_loc_env_release(qs_loc_env,error=error)
      CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
      CALL qs_loc_env_init(qs_loc_env,my_localized_wfn_control,qs_env,error=error)
    END IF
    CALL qs_loc_driver(qs_env,qs_loc_env,ionode,onlyspin1=.TRUE.,error=error)

    DO istate = 1,nexc_search
      centers_wfn(1,istate) =  my_localized_wfn_control%centers_set(1)%array(1,istate)
      centers_wfn(2,istate) =  my_localized_wfn_control%centers_set(1)%array(2,istate)
      centers_wfn(3,istate) =  my_localized_wfn_control%centers_set(1)%array(3,istate)
    END DO

    CALL localized_wfn_control_release(my_localized_wfn_control,error=error)


  ! From the centers and spreads assign each state to one of the selected atoms

  ! If required check the character of the states of each atom
  ! by taking the overlap with atomic-like states

  !Now I know, in the actual set of states, which is the index of the state 
  ! to be excited, for each atom

  ! If loclization has been done copy the localized states in the used MOS array
  
  ! copy the coefficients of the mos in a temporary fm with the right structure
    IF(transition_potential) THEN
      CPPrecondition(ASSOCIATED(ostrength_sm),cp_failure_level,routineP,error,failure)
      CPPrecondition(ASSOCIATED(groundstate_coeff),cp_failure_level,routineP,error,failure)
      DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
            nao=nao, nmo=nmo)
       CALL cp_fm_to_fm(mo_coeff,groundstate_coeff(ispin)%matrix,nmo,1,1)
      END DO
      ! Calculate the operator for the oscillator strengths = sm
      ! In the case of td-dft it should come from the  solution  of tddft algorithm already
   !   CALL oscillator_strength_operator(qs_env,ostrength_sm,error=error)
    END IF

    ! Take the  state_to_be_excited vector from the full set and copy into excvec_coeff
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
  
    DO iat = 1,nexc_atoms
       iatom = xas_control%exc_atoms(iat)
       state_to_be_excited = state_of_atom(iat)

       dft_control % xas_estate  = state_to_be_excited
    !copy the stored initial mos in the  qs_env coefficient  array
    ! if localization has been performed they have to be the localized ones

       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
       CPPrecondition(ASSOCIATED(excvec_coeff),cp_failure_level,routineP,error,failure)
       CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,state_to_be_excited,&
            nao,1,transpose=.TRUE.,error=error)
       CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
            nao,1,transpose=.TRUE.,error=error)
    ! set to zero excene and oscillators

    !  copy the coefficients of state to be excited in a vector coefvec_fm(nao x 1)
      IF (transition_potential) THEN
         NULLIFY(scf_env)
         CALL get_qs_env(qs_env,scf_env=scf_env,error=error)
!         CALL xas_do_tp_scf(xas_env,scf_env,qs_env,error=error)
      ELSE
         ! the states defined by the rotation are the ground state orbitals
         ! the initial state from which I excite should be localized
         ! I take the excitations from lumo to nmo
!         CALL xas_tddft(error)
      END IF

      iter = xas_iter_count*xas_control%each_step
    ! append the spectrum of this iteration to the right output file
!     CALL xas_write(xas_control%excene,xas_control%oscillators, filename, iter, globenv, error)

    END DO

  ! Copy the stored initial MOS in the qs_env coeff array

  ! Release what has to be released 

    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    dft_control % xas_estate    = -1

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE  xas

!!***  *********************************************************************
!!****f* xas_methods/xas_env_init [1.0] *
!!
!!   NAME
!!     xas_env_init
!!
!!   FUNCTION
!!     allocate and initialize the structure needed for the xas calculation
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - xas_env: the environment for XAS  calculations 
!!     - xas_control: parameters for the xas calculation
!!     - qs_env: the qs_env, the xas_env lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [MI]
!!
!!***  *********************************************************************

  SUBROUTINE xas_env_init(xas_env, xas_control, qs_env, error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(xas_control_type)                   :: xas_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routine_name='xas_env_init',&
      routineP = module_name//':'//routine_name

    INTEGER :: homo, i, ispin, istat , nao, nexc_atoms, nexc_search, nmo, nspins
    LOGICAL :: failure
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                   :: ao_mo_fm_pools
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(cp_fm_struct_type), POINTER            :: tmp_fm_struct
    TYPE(cp_para_env_type), POINTER             :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(qs_matrix_pools_type), POINTER         :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s


    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

      NULLIFY(ao_mo_fm_pools, dft_control, matrix_s, mos, mpools, para_env)
      CALL get_qs_env(qs_env=qs_env, &
              dft_control = dft_control, &
              mpools=mpools,& 
              matrix_s=matrix_s, mos=mos, para_env=para_env)
      nexc_search = xas_control%nexc_search
      nexc_atoms = xas_control%nexc_atoms

      CALL mpools_get(mpools, ao_mo_fm_pools= ao_mo_fm_pools)

      NULLIFY(mo_coeff)
      CALL get_mo_set(mos(1)%mo_set,nao=nao,homo=homo,nmo=nmo,mo_coeff=mo_coeff)
      IF(nexc_search < 0) nexc_search = homo
      xas_env%nvirtual = nmo - homo

      ALLOCATE(xas_env%centers_wfn(3,nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%atom_of_state(nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%state_of_atom(nexc_atoms),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%spectrum(2,xas_env%nvirtual),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      ! create a new matrix structure nao x 1
      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
           ncol_global=1,para_env=para_env,context=mo_coeff%matrix_struct%context)
      CALL cp_fm_create (xas_env%excvec_coeff, tmp_fm_struct )
      CALL cp_fm_struct_release ( tmp_fm_struct )


      nspins = SIZE(mos,1)

      IF (xas_control%xas_method==xas_tp_hh .OR. xas_control%xas_method==xas_tp_fh) THEN
      ! initialize operators
        
        ALLOCATE (xas_env%groundstate_coeff(nspins), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DO ispin = 1,nspins
          CALL get_mo_set(mos(ispin)%mo_set,nao=nao,nmo=nmo)
          CALL fm_pool_create_fm(ao_mo_fm_pools(ispin)%pool,&
               xas_env%groundstate_coeff(ispin)%matrix,&
             name="xas_env%mo0"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)

        END DO  ! ispin

        ALLOCATE (xas_env%ostrength_sm(3),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DO i = 1,3
          NULLIFY(xas_env%ostrength_sm(i)%matrix)
          CALL replicate_matrix_structure(matrix_s(1)%matrix, &
                xas_env%ostrength_sm(i)%matrix,"xas_env%op_sm_dir"//&
                TRIM(ADJUSTL(cp_to_string(i))))
           CALL set_matrix(xas_env%ostrength_sm(i)%matrix,0.0_dp)
        END DO 
 
      ELSE IF (xas_control%xas_method==xas_tddft) THEN
      !initialize td-dft calculation
      END IF 

    END IF

  END SUBROUTINE xas_env_init

!***************************************************************************

END MODULE xas_methods
