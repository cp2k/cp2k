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

  USE ai_moments,                      ONLY: diffop
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
  USE cp_files,                        ONLY: open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_set_all,&
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
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: xas_none, xas_tp_hh, xas_tp_fh, xas_tddft,&
                                             xas_1s_type, xas_2s_type, xas_dip_len,&
                                             do_loc_none, op_loc_berry, state_loc_upton
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset
  USE particle_types,                  ONLY: get_particle_set,particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_loc_control,                  ONLY: localized_wfn_control_type,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_retain,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
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
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             block_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE qs_scf,                          ONLY: init_scf_run,&
                                             scf_env_cleanup
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE simulation_cell,                 ONLY: cell_type, pbc
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type, &
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
  USE xas_env_types,                   ONLY: get_xas_env, set_xas_env, &  
                                             xas_env_create, xas_env_release,&
                                             xas_env_retain, &
                                             xas_environment_type
  USE xas_tp_scf,                      ONLY: xas_do_tp_scf

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
!!     the iteration counter is not finilized yet
!!     only the transition potential approach is active
!!     the localization can be switched off, otherwise
!!     it uses by default the berry phase approach
!!     The number of states to be localized is xas_control%nexc_search
!!     In general only the core states are needed
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

    INTEGER :: handle, homo, i, iat, iatom, ikind, ikind_old,  iset, iset_old,&
               isgf, isgfa, ishell, ispin, istat, istate, iter, &
               j, la, la_old, maxnset, maxshell, maxl, maxlev, maxll, nao, natom, nelec, nmo,  &
               nset, nspins, nexc_search, nvirtual, output_unit,& 
               state_to_be_excited, xas_iter_count
    INTEGER, ALLOCATABLE, DIMENSION(:)          :: econf,first_sgf, nlev
    INTEGER, DIMENSION(:), POINTER              :: elec_conf, nshell
    INTEGER, DIMENSION(:, :), POINTER           :: first_sgfa, l, last_sgfa
    INTEGER, DIMENSION(:), POINTER              :: atom_of_state
    INTEGER, DIMENSION(:), POINTER              :: type_of_state
    INTEGER, DIMENSION(:), POINTER              :: state_of_atom
    LOGICAL :: failure, ionode, transition_potential
    REAL(dp) :: charmax, dist, distmin,  maxocc, &
                ra(3), rac(3), rc(3)
    REAL(dp), DIMENSION(:,:), ALLOCATABLE       :: orb_char
    REAL(dp), DIMENSION(:,:), POINTER           :: spectrum
    REAL(dp), DIMENSION(:,:), POINTER           :: centers_wfn
    REAL(dp), DIMENSION(:,:), POINTER           :: vecbuffer
    TYPE(all_potential_type), POINTER           :: all_potential
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                   :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER             :: atomic_kind    
    TYPE(cell_type), POINTER                    :: cell
    TYPE( cp_fm_p_type ), DIMENSION(:),  &
      POINTER                                   :: groundstate_coeff
    TYPE( cp_fm_p_type ), DIMENSION(:,:),  &
      POINTER                                   :: dip_fm_set
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(cp_fm_type), POINTER                   :: fm_work
    TYPE(cp_fm_type), POINTER                   :: excvec_coeff
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(cp_para_env_type), POINTER             :: para_env
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(gth_potential_type), POINTER           :: gth_potential
    TYPE(gto_basis_set_type), POINTER           :: orb_basis_set
    TYPE(localized_wfn_control_type), POINTER   :: my_localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                   :: particle_set
    TYPE( qs_loc_env_new_type ), POINTER        :: qs_loc_env
    TYPE(qs_scf_env_type), POINTER              :: scf_env
    TYPE( real_matrix_p_type ), DIMENSION(:,:), &
      POINTER                                   :: ostrength_sm
    TYPE(xas_environment_type), POINTER         :: xas_env

    INTEGER, DIMENSION(0:5,6) ::  typestate 
    INTEGER, DIMENSION(0:5,6) ::  typestate2 
    INTEGER, DIMENSION(0:5), PARAMETER ::  shift =(/0,1,2,3,4,5/)


    CALL timeset(routineN,"I",'',handle)

    failure = .FALSE.
    transition_potential = .FALSE.

    typestate(0:5,1) = (/ 1, 0, 0, 0, 0, 0/)
    typestate(0:5,2) = (/ 2, 3, 0, 0, 0, 0/)
    typestate(0:5,3) = (/ 4, 5, 6, 0, 0, 0/)
    typestate(0:5,4) = (/ 7, 8, 9,10, 0, 0/)
    typestate(0:5,5) = (/11,12,13,14,15, 0/)
    typestate(0:5,6) = (/16,17,18,19,20,21/)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)


    NULLIFY(xas_env, groundstate_coeff, ostrength_sm, dip_fm_set, fm_work)
    NULLIFY(excvec_coeff, spectrum, state_of_atom, centers_wfn, qs_loc_env, cell)

!   Create the xas environment
    CALL get_qs_env(qs_env,xas_env=xas_env,error=error)
    IF (.NOT.ASSOCIATED(xas_env)) THEN
      CALL xas_env_create(xas_env, error=error)
      CALL xas_env_init(xas_env, xas_control, qs_env, globenv, error=error)
      CALL set_qs_env(qs_env,xas_env=xas_env, error=error)
      CALL xas_env_release(xas_env,error=error)
      CALL get_qs_env(qs_env,xas_env=xas_env,error=error)
    END IF

!   Initialize the type of calculation
    NULLIFY(atomic_kind_set, cell, dft_control, mos, para_env, particle_set)
    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
         cell = cell, dft_control=dft_control, &
         mos=mos, para_env=para_env, particle_set=particle_set )
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxlgto=maxl,maxshell=maxshell,maxnset=maxnset)

    nexc_search = xas_control%nexc_search
    nspins = dft_control%nspins
    transition_potential = (xas_control%xas_method==xas_tp_hh).OR.&
                           (xas_control%xas_method==xas_tp_fh)
    IF(nspins ==1 .AND. transition_potential) THEN
       CALL stop_program(routineP,"xas with tp method requires LSD calculations")
    END IF

!   Set of states among which there is the state to be excited
    CALL get_mo_set(mos(1)%mo_set,nao=nao,homo=homo)
    IF(nexc_search < 0) nexc_search = homo

    CALL get_xas_env(xas_env=xas_env,iter_count=xas_iter_count,nvirtual=nvirtual,&
         centers_wfn=centers_wfn,atom_of_state=atom_of_state,&
         type_of_state = type_of_state, state_of_atom=state_of_atom,spectrum=spectrum,&
         groundstate_coeff=groundstate_coeff,excvec_coeff=excvec_coeff,&
         dip_fm_set=dip_fm_set,fm_work=fm_work)

    xas_iter_count =  xas_iter_count + 1
    CALL set_xas_env(xas_env=xas_env,iter_count=xas_iter_count)
    
   !I know from which atoms I want to excite (O) and which type of states (1s) (input)

   !Define the qs_loc_env : to find centers, spread and possibly localize them
    my_localized_wfn_control =>  xas_control%localized_wfn_control
    CALL localized_wfn_control_retain(my_localized_wfn_control,error=error)
    CPPrecondition(ASSOCIATED(my_localized_wfn_control),cp_failure_level,routineP,error,failure)

    CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)

    IF(.NOT.(ASSOCIATED(qs_loc_env))) THEN
      CALL qs_loc_env_create(qs_loc_env,error=error)
      CALL set_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
      CALL qs_loc_env_release(qs_loc_env,error=error)
      CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
      CALL qs_loc_env_init(qs_loc_env,my_localized_wfn_control,qs_env,onlyspin1=.TRUE.,error=error)

    END IF
    CALL qs_loc_driver(qs_env,qs_loc_env,ionode,onlyspin1=.TRUE.,error=error)

    ! scratch array for the state
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, maxocc=maxocc)
    natom = SIZE(particle_set)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set=particle_set, first_sgf=first_sgf)
    ALLOCATE(orb_char(0:maxl,maxnset),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nlev(0:maxl),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (econf(0:maxl),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ikind_old = 0
    DO istate = 1,nexc_search
      centers_wfn(1,istate) =  my_localized_wfn_control%centers_set(1)%array(1,istate)
      centers_wfn(2,istate) =  my_localized_wfn_control%centers_set(1)%array(2,istate)
      centers_wfn(3,istate) =  my_localized_wfn_control%centers_set(1)%array(3,istate)

      ! Assign the state to the closest atom
      distmin = 100.0_dp
      charmax = 0.0_dp

      DO iat = 1,xas_control%nexc_atoms
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

      !Character of the state
       CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,istate,&
            nao,1,transpose=.TRUE.,error=error)
  
       iatom = atom_of_state(istate)
      
       atomic_kind =>  particle_set(iatom)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
              all_potential=all_potential,&
              gth_potential=gth_potential,&
              orb_basis_set=orb_basis_set,&
              kind_number=ikind)

       IF(ikind /= ikind_old) THEN
         ikind_old = ikind
         IF (ASSOCIATED(all_potential)) THEN
            CALL get_potential(potential=all_potential,elec_conf=elec_conf)
         ELSE IF (ASSOCIATED(gth_potential)) THEN
            CALL get_potential(potential=gth_potential,elec_conf=elec_conf)
         END IF
         maxll = MIN(SIZE(elec_conf) - 1,maxl)
         econf(:) = 0.0_dp
         econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll),dp)

         typestate2 = typestate
         DO la = 0,maxll
            nelec = maxocc*REAL(2*la + 1,dp)
            maxlev = int(econf(la)/nelec) 
            IF(MOD(econf(la),nelec)/=0) maxlev = maxlev +1
            DO i = maxlev+shift(la)+1, 6
              typestate2(la,i ) = typestate(la,maxlev+shift(la))
            END DO 
         END DO  

         CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                nset=nset,&
                                nshell=nshell,&
                                l=l,&
                                first_sgf=first_sgfa,&
                                last_sgf=last_sgfa)
       END IF

       nlev(:) = 0
       orb_char(:,:) = 0.0_dp

       iset_old = 0
       DO iset=1,nset
         la_old = -1
         DO ishell=1,nshell(iset)
           la = l(ishell,iset)
           nelec = maxocc*REAL(2*la + 1,dp)
           IF (la/=la_old) THEN
             nlev(la) = nlev(la) + 1
             la_old = la
           END IF
           DO isgfa=first_sgfa(ishell,iset),last_sgfa(ishell,iset)
             isgf = first_sgf(iatom) + isgfa - 1
             orb_char(la , nlev(la)) = orb_char(la , nlev(la)) + &
                      vecbuffer(1,isgf)* vecbuffer(1,isgf)
           END DO  ! isgfa

           IF(orb_char(la , nlev(la)) >= charmax) THEN
              charmax = orb_char(la , nlev(la))
              type_of_state(istate) = typestate2(la , nlev(la)+shift(la))
           END IF

         END DO  ! ishell
       END DO  ! iset
    END DO  ! istate

    CALL localized_wfn_control_release(my_localized_wfn_control,error=error)

  ! If required check the character of the states of each atom
  ! by taking the overlap with atomic-like states

  ! In the set of states, which is the index of the one to be excited for iatom
    DO iat = 1,xas_control%nexc_atoms
      iatom = xas_control%exc_atoms(iat)

      DO istate = 1,nexc_search
         IF(atom_of_state(istate) == iatom .AND. & 
            type_of_state(istate) == xas_control%state_type) THEN
            state_of_atom(iat) = istate
            EXIT
         END IF
      END DO 
    END DO 


  ! copy the coefficients of the mos in a temporary fm with the right structure
    IF(transition_potential) THEN
      CPPrecondition(ASSOCIATED(groundstate_coeff),cp_failure_level,routineP,error,failure)
      DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
            nao=nao, nmo=nmo)
       CALL cp_fm_to_fm(mo_coeff,groundstate_coeff(ispin)%matrix,nmo,1,1)
      END DO

      IF(xas_control%dipole_form==xas_dip_len) THEN
        ostrength_sm =>  qs_loc_env%op_sm_set
      ELSE
      ! Calculate the operator 
        CALL get_xas_env(xas_env=xas_env,ostrength_sm=ostrength_sm)
        CALL get_dip_vel_operator(ostrength_sm,qs_env,error=error) 
      END IF
      ! In td-dft case the oscillator strength come from  linear response orbitals
    END IF

    DO iat = 1,xas_control%nexc_atoms
       iatom = xas_control%exc_atoms(iat)
       ! determine which state has to be excited in the global list
       state_to_be_excited = state_of_atom(iat)
       dft_control % xas_estate  = state_to_be_excited

      ! Take the state_to_be_excited vector from the full set and copy into excvec_coeff
       CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
       CPPrecondition(ASSOCIATED(excvec_coeff),cp_failure_level,routineP,error,failure)
       CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,state_to_be_excited,&
            nao,1,transpose=.TRUE.,error=error)
       CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
            nao,1,transpose=.TRUE.,error=error)

      IF (transition_potential) THEN
         IF(ionode) THEN
           WRITE(UNIT=output_unit,FMT='(/,/,T2,A,I6)') &
                  " START XAS CALCULATION WITH TP APPROACH FOR ATOM   ", iatom 
         END IF
         NULLIFY(scf_env)
         CALL get_qs_env(qs_env,scf_env=scf_env,error=error)
         CALL init_scf_run(scf_env=scf_env,qs_env=qs_env,&
                 globenv=globenv, error=error)

         CALL xas_do_tp_scf(xas_control,xas_env,scf_env,qs_env,globenv,error=error)

         IF(iat == xas_control%nexc_atoms) THEN 
            CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)
         END IF


         ! Calculate the spectrum
         IF(xas_control%dipole_form == xas_dip_len) THEN
           CALL spectrum_dip_len(dip_fm_set,ostrength_sm,mos,excvec_coeff,fm_work,&
                 cell,spectrum,dft_control%xas_estate,error=error)
         ELSE

         END IF
      ELSE
         ! the states defined by the rotation are the ground state orbitals
         ! the initial state from which I excite should be localized
         ! I take the excitations from lumo to nmo
!         CALL xas_tddft(error)
      END IF

      iter = xas_iter_count*xas_control%each_step
    ! append the spectrum of this iteration to the right output file
      CALL xas_write(spectrum, xas_control%xas_rootfname, &
           iatom, iter, ionode, error=error)
      DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
            nao=nao, nmo=nmo)
       CALL cp_fm_to_fm(groundstate_coeff(ispin)%matrix,mo_coeff,nmo,1,1)
      END DO

    END DO

  ! Release what has to be released 

    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(orb_char,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(nlev,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
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

  SUBROUTINE xas_env_init(xas_env, xas_control, qs_env, globenv, error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(xas_control_type)                   :: xas_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      POINTER                                :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routine_name='xas_env_init',&
      routineP = module_name//':'//routine_name

    INTEGER :: homo, i, ispin, istat , j, nao, nexc_atoms, nexc_search,&
               n_mo(2), nmo, nmoloc(2), nspins
    LOGICAL :: failure
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                   :: ao_mo_fm_pools
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(cp_fm_struct_type), POINTER            :: tmp_fm_struct
    TYPE(cp_para_env_type), POINTER             :: para_env
    TYPE(localized_wfn_control_type), POINTER   :: my_localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(qs_matrix_pools_type), POINTER         :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s

!------------------------------------------------------------------------------


    failure=.FALSE.

    n_mo(1:2) = 0

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

      n_mo(1:2) = nmo

      ALLOCATE(xas_env%centers_wfn(3,nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%atom_of_state(nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%type_of_state(nexc_search),STAT=istat)
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


      NULLIFY(tmp_fm_struct)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=1,&
           ncol_global=nexc_search,para_env=para_env,&
           context=mo_coeff%matrix_struct%context)
      CALL cp_fm_create (xas_env%excvec_overlap, tmp_fm_struct )
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

        CALL fm_pool_create_fm(ao_mo_fm_pools(1)%pool,&
               xas_env%fm_work,name="xas_env%fm_work",error=error)

        NULLIFY(tmp_fm_struct)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=1,&
             ncol_global=xas_env%nvirtual,para_env=para_env,&
             context=mo_coeff%matrix_struct%context)
        ALLOCATE (xas_env%dip_fm_set(2,3),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DO i = 1,3
          DO j = 1,2
             CALL cp_fm_create (xas_env%dip_fm_set(j,i)%matrix, tmp_fm_struct )
          END DO 
        END DO 
        CALL cp_fm_struct_release ( tmp_fm_struct )

        IF(xas_control%dipole_form /= xas_dip_len) THEN
          ALLOCATE (xas_env%ostrength_sm(1,3),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO i = 1,3
             NULLIFY(xas_env%ostrength_sm(1,i)%matrix)
             CALL replicate_matrix_structure(matrix_s(1)%matrix, &
                  xas_env%ostrength_sm(1,i)%matrix,"xas_env%op_sm_dir"//&
                  "-"//TRIM(ADJUSTL(cp_to_string(i))))
             CALL set_matrix(xas_env%ostrength_sm(1,i)%matrix,0.0_dp)
          END DO 
        END IF
 
      ELSE IF (xas_control%xas_method==xas_tddft) THEN
      !initialize td-dft calculation
      END IF 

   !Define the qs_loc_env : to find centers, spread and possibly localize them
    my_localized_wfn_control =>  xas_control%localized_wfn_control
    CALL localized_wfn_control_retain(my_localized_wfn_control,error=error)
    CPPrecondition(ASSOCIATED(my_localized_wfn_control),cp_failure_level,routineP,error,failure)
    IF(.NOT. xas_control%localized_wfn) THEN
      my_localized_wfn_control%localization_method = do_loc_none
    END IF

    ! The localization is needed only for the states up to nexc_search
    my_localized_wfn_control%set_of_states = state_loc_upton
    my_localized_wfn_control%operator_type = op_loc_berry
    my_localized_wfn_control%nloc_states(1) = nexc_search
    my_localized_wfn_control%nloc_states(2) = 0
    my_localized_wfn_control%lu_bound_states(1,1) = 1
    my_localized_wfn_control%lu_bound_states(2,1) = nexc_search

    nmoloc(1:2) = my_localized_wfn_control%nloc_states(1:2)
    CALL set_loc_wfn_lists(xas_control%localized_wfn_control,nmoloc,n_mo,nspins,error=error)
    CALL set_loc_centers(xas_control%localized_wfn_control,nmoloc,nspins,globenv,error=error)

    CALL localized_wfn_control_release(my_localized_wfn_control,error=error)
    END IF

  END SUBROUTINE xas_env_init

!***************************************************************************
!!****f* xas_methods/xas_write [1.0] *
!!
!!   NAME
!!     xas_write
!!
!!   FUNCTION
!!     write the spectrum for each atom in a different output file
!!
!!   NOTES
!!     the iteration counter is not finilized yet
!!
!!
!!   ARGUMENTS
!!     - spectrum: temporary variable contaning the spectrum for the actual atom
!!     - xas_rootfname: root of the filename where the spectrum is written
!!     - iatom: index of the excited atom
!!     - iter: iteration of the xas calculation
!!     - ionode: logical assigning the i/o node
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

  SUBROUTINE xas_write(spectrum, xas_rootfname, &
                       iatom, iter, ionode, error)

    REAL(dp), DIMENSION(:,:),POINTER            :: spectrum
    CHARACTER(len=50)                           :: xas_rootfname
    INTEGER, INTENT(IN)                         :: iatom, iter
    LOGICAL, INTENT(IN)                         :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'xas_write', &
      routineP = module_name//':'//routine_name

    CHARACTER(len = 50)                      :: filename
    INTEGER :: ia, ie, istate, out_sp
    LOGICAL :: failure, file_exists
    
    failure = .FALSE.
    out_sp = 23
    file_exists = .FALSE.

    IF(ionode) THEN
      ! check if file exists
      CALL xstring (xas_rootfname, ia, ie )
      filename = xas_rootfname(ia:ie)//'_at_'//TRIM(ADJUSTL(cp_to_string(iatom)))

      INQUIRE (FILE=filename,EXIST=file_exists)

      ! open file
      IF(file_exists) THEN
         CALL open_file(file_name=TRIM(filename),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="APPEND",&
                     file_status="OLD",&
                     unit_number=out_sp)
      ELSE
         CALL open_file(file_name=TRIM(filename),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="REWIND",&
                     file_status="REPLACE",&
                     unit_number=out_sp)
      END IF
      WRITE(out_sp,'(A,I10)') " Iteration Step ", iter

      DO istate = 1,SIZE(spectrum,2)
        WRITE(out_sp,'(I6,2F16.8)') istate, spectrum(1,istate), spectrum(2,istate)
      END DO 
        

      close(out_sp)
    END IF

  END SUBROUTINE xas_write

!***************************************************************************
!!****f* xas_methods/spectrum_dip_len [1.0] *
!!
!!   NAME
!!    spectrum_dip_len 
!!
!!   FUNCTION
!!     Calculation of the spectrum when the dipole approximation
!!     in the length form is used. The position operator is calculated 
!!     by the berry phase approach
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - fm_set: components of the position operator in a full matrix form
!!               already multiplied by the coefficiets
!!               only the terms <C_i Op C_f> are calculated where
!!               C_i are the coefficients of the excited state
!!     - op_sm: components of the position operator for the dipole 
!!              in a sparse matrix form (cos and sin)
!!              calculated for the basis functions
!!     - mos: wavefunctions coefficients
!!     - excvec: coefficients of the excited orbital
!!     - fm_work: work space
!!     - cell: parameters for the simulation cell
!!     - spectrum: temporary variable contaning the spectrum for the actual atom
!!     - estate: index of the excited state
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!!***  *********************************************************************

  SUBROUTINE spectrum_dip_len(fm_set,op_sm,mos,excvec,fm_work,cell,spectrum,estate,error)

    TYPE( cp_fm_p_type ), DIMENSION(:,:),  &
      POINTER                                   :: fm_set
    TYPE( real_matrix_p_type ), DIMENSION(:,:), &
      POINTER                                   :: op_sm
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(cp_fm_type), POINTER                   :: excvec
    TYPE(cp_fm_type), POINTER                   :: fm_work
    TYPE(cell_type), POINTER                    :: cell
    REAL(dp), DIMENSION(:,:), POINTER           :: spectrum
    INTEGER, INTENT(IN)                         :: estate

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "spectrum_dip_len", &
      routineP = module_name//':'//routineN

    COMPLEX(KIND=dp)                            :: z
    INTEGER :: homo, i, istate, j, nao, nvirtual, nmo
    LOGICAL :: failure
    REAL(dp) :: ene_i, ene_f, imagpart, ra(3), realpart
    REAL(dp), DIMENSION(:), POINTER             :: eigenvalues
    TYPE(cp_fm_type), POINTER                   :: mo_coeff

    failure = .FALSE.
    NULLIFY(eigenvalues,mo_coeff)

    CPPrecondition(ASSOCIATED(fm_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(spectrum),cp_failure_level,routineP,error,failure)
    nvirtual = SIZE(spectrum,2)
    
    IF(.NOT. failure) THEN
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo,&
           homo =homo, eigenvalues=eigenvalues)
      DO i=1,SIZE(fm_set,2)
        DO j = 1,SIZE(fm_set,1)
           CPPrecondition(ASSOCIATED(fm_set(j,i)%matrix),cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(fm_set(j,i)%matrix, 0.0_dp, error=error)
           CALL cp_sm_fm_multiply(op_sm(j,i)%matrix,mo_coeff,fm_work,ncol=nmo)
           CALL cp_fm_gemm("T","N",1,nvirtual,nao,1.0_dp,excvec,&
                  fm_work,0.0_dp, fm_set(j,i)%matrix,b_first_col=homo+1,error=error)
        END DO 
      END DO
      
      ene_i = eigenvalues(estate)
      DO istate = 1,nvirtual
         ene_f = eigenvalues(homo+istate)
         DO i = 1,3
           CALL cp_fm_get_element(fm_set(1,i)%matrix,1,istate,realpart)
           CALL cp_fm_get_element(fm_set(2,i)%matrix,1,istate,imagpart)
           z = CMPLX(realpart,imagpart,dp)
           ra(i) = ( cell % hmat ( i, i ) / twopi ) * AIMAG ( LOG ( z ) )
         END DO 
         spectrum(1,istate) = ene_f - ene_i
         spectrum(2,istate) = ra(1)*ra(1)+ra(2)*ra(2)+ra(3)*ra(3)
      END DO 
    END IF

  END SUBROUTINE spectrum_dip_len

!***************************************************************************
!!****f* xas_methods/get_dip_vel_operator [1.0] *
!!
!!   NAME
!!     get_dip_vel_operator
!!
!!   FUNCTION
!!     Calculation of the components of the dipole operator in the velocity form
!!     The elements of the  sparse matrices are the integrals in the
!!     basis functions
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - op_sm: components of the position operator for the dipole 
!!              in a sparse matrix form (cos and sin)
!!              calculated for the basis functions
!!     - qs_env: enviroment for the lists and the basis sets
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     06.2005 created [MI]
!!
!!***  *********************************************************************

  SUBROUTINE get_dip_vel_operator(op_sm, qs_env, error )

    TYPE( real_matrix_p_type ), DIMENSION(:,:), &
      POINTER                                   :: op_sm
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'get_berry_operator', &
      routineP = module_name//':'//routine_name

    INTEGER  :: handle, i, iab, iat, iatom, icol, ikind, ilist, inode, irow, &
                iset, istat, ithread, j, jat, jatom, jkind, jset, last_jatom, &
                ldab, ldsa, ldsb, ldwork, maxl, ncoa, ncob, nkind, nlist,& 
                nnode, nrow, nthread, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER             :: la_max, la_min, lb_max, lb_min
    INTEGER, DIMENSION(:), POINTER             :: npgfa, npgfb
    INTEGER, DIMENSION(:), POINTER             :: nsgfa, nsgfb
    INTEGER, DIMENSION(:,:), POINTER           :: first_sgfa,first_sgfb
    LOGICAL :: failure, new_atom_b
    REAL(KIND=dp) :: dab, rab2
    REAL(KIND=dp), DIMENSION(3)                :: ra, rb, rab
    REAL(KIND = dp), DIMENSION(:), POINTER     :: set_radius_a,set_radius_b
    REAL(KIND=dp), DIMENSION(:,:), POINTER     :: rpgfa, rpgfb, sphi_a,sphi_b, &
                                                  zeta, zetb 
    REAL(KIND=dp), DIMENSION(:,:), POINTER     :: work
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER   :: difab,workt 
    REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER :: difabt 
    TYPE(atomic_kind_type), DIMENSION(:),&
      POINTER                                  :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER            :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
          POINTER                              :: op_dip
    TYPE(block_p_type), DIMENSION(:,:), &
          POINTER                              :: op_dipt
    TYPE(cell_type), POINTER                   :: cell
    TYPE(gto_basis_set_type), POINTER          :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:),           POINTER          :: sab_orb
    TYPE(neighbor_list_type), POINTER          :: sab_orb_neighbor_list,&
                                                  sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER          :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                  :: particle_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routine_name,"I"," ",handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY (cell, particle_set, orb_basis_set)
    NULLIFY (sab_orb,sab_orb_neighbor_list,sab_orb_neighbor_list_local,&
             sab_orb_neighbor_node)
    NULLIFY (difab,difabt,op_dip,op_dipt,work,workt)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,sab_orb=sab_orb)

    nkind = SIZE(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=ldwork, maxlgto=maxl )
    ldab = ldwork
    nthread = 1
!$  nthread = omp_get_max_threads()
    CALL reallocate(difabt,1,ldab,1,ldab,1,3,0,nthread-1)
    CALL reallocate(workt,1,ldwork,1,ldwork,0,nthread-1)

    ALLOCATE (op_dipt(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, 3
      DO j=0,nthread-1
        NULLIFY (op_dipt(i,j)%block)
      END DO
    END DO 
 
    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             set_radius=set_radius_a,&
                             sphi=sphi_a,&
                             zet=zeta)
      ldsa = SIZE(sphi_a,1)
      DO jkind=1,nkind


        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgfb,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)


        iab = ikind + nkind*(jkind - 1)
        ldsb = SIZE(sphi_b,1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,difab,work,ilist,iatom,jatom)&
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb)&
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_dip,rab,rab2,dab,i)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,difabt,qs_env)&
!$OMP shared(op_dipt)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm_set)
!$OMP shared(atomic_kind_set,particle_set,cell)
        ithread = 0
!$      ithread = omp_get_thread_num()

        difab => difabt(:,:,:,ithread)
        work  => workt(:,:,ithread)
        op_dip => op_dipt(:,ithread)        

!$OMP do
       DO ilist = 1,nlist

!$OMP critical (qs_core_neighbor_list)
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
!$OMP end critical (qs_core_neighbor_list)

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)
          ra = pbc(particle_set(iatom)%r,cell)
          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom, r=rab)

!            rb = ra - rab
            rb = rab + ra

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF

            IF (new_atom_b) THEN
              IF (iatom <= jatom) THEN
                irow = iatom
                icol = jatom
              ELSE
                irow = jatom
                icol = iatom
              END IF

              DO i = 1,3
                 NULLIFY(op_dip(i)%block)
                 CALL get_block_node(matrix=op_sm(1,i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_dip(i)%block)
                END DO  
            END IF  ! new_atom_b
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            nrow = 0
            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

!            *** Calculate the primitive overlap integrals ***
                   CALL diffop(la_max(iset),npgfa(iset),zeta(:,iset),&
                        rpgfa(:,iset),la_min(iset), lb_max(jset),npgfb(jset),&
                        zetb(:,jset),rpgfb(:,jset),lb_min(jset),rab,difab)

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb,ncob,&
                        1.0_dp,difab(1,1,1),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa,nsgfb,ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(1)%block(sgfa,sgfb),&
                          SIZE(op_dip(1)%block,1))
                   ELSE
                     CALL dgemm("T","N",nsgfb,nsgfa,ncoa,&
                          1.0_dp,work(1,1),ldwork,&
                          sphi_a(1,sgfa),ldsa,&
                          1.0_dp,op_dip(1)%block(sgfb,sgfa),&
                          SIZE(op_dip(1)%block,1))
                   END IF

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb,ncob,&
                        1.0_dp,difab(1,1,2),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa,nsgfb,ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(2)%block(sgfa,sgfb),&
                          SIZE(op_dip(2)%block,1))
                   ELSE
                     CALL dgemm("T","N",nsgfb,nsgfa,ncoa,&
                          1.0_dp,work(1,1),ldwork,&
                          sphi_a(1,sgfa),ldsa,&
                          1.0_dp,op_dip(2)%block(sgfb,sgfa),&
                          SIZE(op_dip(2)%block,1))
                   END IF

!            *** Contraction ***
                   CALL dgemm("N","N",ncoa,nsgfb,ncob,&
                        1.0_dp,difab(1,1,3),ldsb,&
                        0.0_dp,work(1,1),ldwork)
                   IF(iatom<=jatom) THEN
                     CALL dgemm("T","N",nsgfa,nsgfb,ncoa,&
                          1.0_dp,sphi_a(1,sgfa),ldsa,&
                          work(1,1),ldwork,&
                          1.0_dp,op_dip(3)%block(sgfa,sgfb),&
                          SIZE(op_dip(3)%block,1))
                   ELSE
                     CALL dgemm("T","N",nsgfb,nsgfa,ncoa,&
                          1.0_dp,work(1,1),ldwork,&
                          sphi_a(1,sgfa),ldsa,&
                          1.0_dp,op_dip(3)%block(sgfb,sgfa),&
                          SIZE(op_dip(3)%block,1))
                   END IF

                END IF  !  >= dab

              END DO ! jset

              nrow = nrow + ncoa

            END DO ! iset

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

       END DO ! iatom => atom A

!$OMP end parallel

      END DO ! jkind

    END DO ! ikind

    DO i = 1,3
      DO j = 0,nthread-1
           NULLIFY(op_dipt(i,j)%block)
      END DO
    END DO 
    DEALLOCATE(op_dipt, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(difabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE get_dip_vel_operator

END MODULE xas_methods
