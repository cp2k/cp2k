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
  USE physcon,                         ONLY: evolt
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
                                             mo_set_p_type,mo_set_type,&
                                             write_mo_set
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
  USE qs_operators_ao,                 ONLY: p_xyz_ao
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
    INTEGER, DIMENSION(:), POINTER              :: elec_conf, nshell, npgf
    INTEGER, DIMENSION(:, :), POINTER           :: first_sgfa, l, last_sgfa
    INTEGER, DIMENSION(:), POINTER              :: atom_of_state
    INTEGER, DIMENSION(:), POINTER              :: type_of_state
    INTEGER, DIMENSION(:), POINTER              :: state_of_atom
    LOGICAL :: failure, ionode, transition_potential
    REAL(dp) :: charmax, dist, distmin,  maxocc, &
                ra(3), rac(3), rc(3)
    REAL(dp), DIMENSION(:), POINTER             :: occupation_numbers
    REAL(dp), DIMENSION(:,:), ALLOCATABLE       :: orb_char
    REAL(dp), DIMENSION(:,:), POINTER           :: spectrum, zet
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
    TYPE(real_matrix_p_type), DIMENSION(:),  &
      POINTER                                   :: op_sm
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

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)


    NULLIFY(xas_env, groundstate_coeff, ostrength_sm, op_sm, dip_fm_set, fm_work)
    NULLIFY(excvec_coeff, spectrum, state_of_atom, centers_wfn, qs_loc_env, cell)
    NULLIFY(occupation_numbers,my_localized_wfn_control)
    NULLIFY(elec_conf,l,nshell,zet,first_sgfa, last_sgfa)
    NULLIFY(atom_of_state,type_of_state,state_of_atom)

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
!IF(ionode) WRITE(*,*) 'homo ', homo, evolt

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
    ALLOCATE (op_sm(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ikind_old = 0
    DO istate = 1,nexc_search
      centers_wfn(1,istate) =  my_localized_wfn_control%centers_set(1)%array(1,istate)
      centers_wfn(2,istate) =  my_localized_wfn_control%centers_set(1)%array(2,istate)
      centers_wfn(3,istate) =  my_localized_wfn_control%centers_set(1)%array(3,istate)

 !if(ionode) WRITE(*,'(A,I4,3f10.5)') 'center ', istate, centers_wfn(1:3, istate)

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
! if(ionode) WRITE(*,'(I4,A,I4,6f10.5)') istate, 'atom ', atom_of_state(istate),  particle_set(atom_of_state(istate))%r(1:3), centers_wfn(1:3,istate)

      !Character of the state
       CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,istate,&
            nao,1,transpose=.TRUE.,error=error)
  
       iatom = atom_of_state(istate)
       NULLIFY(atomic_kind,all_potential,gth_potential,orb_basis_set)
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
                                l=l, zet=zet, npgf=npgf,&
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
           nelec = REAL(2*la + 1,dp)
           IF(la == 0) THEN
             IF(zet(npgf(iset),iset) .GT. 3._dp) THEN
                nlev(la) = 1
             ELSE
                nlev(la) = 2
             END IF
           ELSE
             EXIT
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
! if(ionode) WRITE(*,*) 'type ', type_of_state(istate)
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
! IF(ionode) WRITE(*,'(A,2I4)') 'state of iatom ',iatom , state_of_atom(iat)
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
        DO i = 1,3
          NULLIFY(op_sm(i)%matrix)
          op_sm(i)%matrix => ostrength_sm(1,i)%matrix
        END DO
        CALL p_xyz_ao(op_sm,qs_env,error=error)
      END IF

      ! In td-dft case the oscillator strength come from  linear response orbitals
    END IF

    DO iat = 1,xas_control%nexc_atoms
       iatom = xas_control%exc_atoms(iat)
       ! determine which state has to be excited in the global list
       state_to_be_excited = state_of_atom(iat)
       dft_control % xas_estate  = state_to_be_excited
!IF (ionode ) WRITE(*,*) 'estate ',  dft_control % xas_estate
      ! Take the state_to_be_excited vector from the full set and copy into excvec_coeff
       CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,&
            occupation_numbers=occupation_numbers,homo=homo)
       occupation_numbers(1:homo)  = 1.0_dp
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

         CALL xas_do_tp_scf(xas_control,xas_env,iatom,scf_env,qs_env,globenv,error=error)

         IF(iat == xas_control%nexc_atoms) THEN 
            CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)
         END IF

!   *** Write last wavefunction to screen ***
         DO ispin=1,dft_control%nspins
           CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
         ENDDO

         ! Calculate the spectrum
         IF(xas_control%dipole_form == xas_dip_len) THEN
           CALL spectrum_dip_len(xas_control,dip_fm_set,ostrength_sm,mos,excvec_coeff,fm_work,&
                 cell,spectrum,dft_control%xas_estate,error=error)
         ELSE
           CALL spectrum_dip_vel(xas_control,dip_fm_set,ostrength_sm,mos,excvec_coeff,fm_work,&
                spectrum,dft_control%xas_estate,error=error)
         END IF
      ELSE
         ! the states defined by the rotation are the ground state orbitals
         ! the initial state from which I excite should be localized
         ! I take the excitations from lumo to nmo
!         CALL xas_tddft(error)
      END IF

      iter = xas_iter_count*xas_control%each_step
    ! append the spectrum of this iteration to the right output file
      CALL xas_write(xas_control,spectrum, dft_control%xas_estate, &
           homo, iatom, iter, ionode, error=error)
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
    DEALLOCATE(op_sm,STAT=istat)
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

      NULLIFY( dft_control, matrix_s, mos, mpools, para_env)
      CALL get_qs_env(qs_env=qs_env, &
              dft_control = dft_control, &
              mpools=mpools,& 
              matrix_s=matrix_s, mos=mos, para_env=para_env)
      nexc_search = xas_control%nexc_search
      nexc_atoms = xas_control%nexc_atoms

      CALL mpools_get(mpools, ao_mo_fm_pools= xas_env%ao_mo_fm_pools)

      NULLIFY(mo_coeff)
      CALL get_mo_set(mos(1)%mo_set,nao=nao,homo=homo,nmo=nmo,mo_coeff=mo_coeff)
      IF(nexc_search < 0) nexc_search = homo
      IF(xas_control%emission) THEN
         xas_env%nvirtual = nmo
      ELSE
        xas_env%nvirtual = nmo - homo
      END IF

      n_mo(1:2) = nmo

      ALLOCATE(xas_env%centers_wfn(3,nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%atom_of_state(nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%type_of_state(nexc_search),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%state_of_atom(nexc_atoms),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(xas_env%spectrum(5,xas_env%nvirtual),STAT=istat)
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
          CALL fm_pool_create_fm(xas_env%ao_mo_fm_pools(ispin)%pool,&
               xas_env%groundstate_coeff(ispin)%matrix,&
             name="xas_env%mo0"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)

        END DO  ! ispin

        CALL fm_pool_create_fm(xas_env%ao_mo_fm_pools(1)%pool,&
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

  SUBROUTINE xas_write(xas_control, spectrum, estate, homo, iatom, iter, ionode, error)

    TYPE(xas_control_type)                      :: xas_control
    REAL(dp), DIMENSION(:,:),POINTER            :: spectrum
    INTEGER, INTENT(IN)                         :: estate, homo, iatom, iter
    LOGICAL, INTENT(IN)                         :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'xas_write', &
      routineP = module_name//':'//routine_name

    CHARACTER(len = 50)                      :: filename
    INTEGER :: ia, ie, istate, istate0, out_sp
    LOGICAL :: failure, file_exists
    
    failure = .FALSE.
    out_sp = 23
    file_exists = .FALSE.

    IF(ionode) THEN
      ! check if file exists
      CALL xstring (xas_control%xas_rootfname, ia, ie )
      filename = xas_control%xas_rootfname(ia:ie)//'_at_'//TRIM(ADJUSTL(cp_to_string(iatom)))

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

      IF(xas_control%emission) THEN
         istate0 = estate
         DO istate = estate,SIZE(spectrum,2)
            WRITE(out_sp,'(I6,5F16.8)') istate, spectrum(1,istate)*evolt, &
                  spectrum(2,istate), spectrum(3,istate),&    
                  spectrum(4,istate), spectrum(5,istate)
         END DO 
      ELSE
         istate0 = homo
         DO istate = 1,SIZE(spectrum,2)
            WRITE(out_sp,'(I6,5F16.8)') homo+istate, spectrum(1,istate)*evolt, &
                  spectrum(2,istate), spectrum(3,istate),&    
                  spectrum(4,istate), spectrum(5,istate)
         END DO 
      END IF

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

  SUBROUTINE spectrum_dip_len(xas_control,fm_set,op_sm,mos,excvec,fm_work,cell,spectrum,estate,error)

    TYPE(xas_control_type)                      :: xas_control
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
    INTEGER :: homo, i, istate, istate0, j, nao, nvirtual, nmo
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
      istate0 = homo
      IF(xas_control%emission) istate0 = 0
      DO i=1,SIZE(fm_set,2)
        DO j = 1,SIZE(fm_set,1)
           CPPrecondition(ASSOCIATED(fm_set(j,i)%matrix),cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(fm_set(j,i)%matrix, 0.0_dp, error=error)
           CALL cp_sm_fm_multiply(op_sm(j,i)%matrix,mo_coeff,fm_work,ncol=nmo)
           CALL cp_fm_gemm("T","N",1,nvirtual,nao,1.0_dp,excvec,&
                  fm_work,0.0_dp, fm_set(j,i)%matrix,b_first_col=istate0+1,error=error)
        END DO 
      END DO
      
      ene_i = eigenvalues(estate)
      DO istate = 1,nvirtual
         ene_f = eigenvalues(istate0+istate)
         DO i = 1,3
           CALL cp_fm_get_element(fm_set(1,i)%matrix,1,istate,realpart)
           CALL cp_fm_get_element(fm_set(2,i)%matrix,1,istate,imagpart)
           z = CMPLX(realpart,imagpart,dp)
           ra(i) = ( cell % hmat ( i, i ) / twopi ) * AIMAG ( LOG ( z ) )
         END DO 
         spectrum(1,istate) = ene_f - ene_i
         spectrum(2,istate) = ra(1)*ra(1)
         spectrum(3,istate) = ra(2)*ra(2)
         spectrum(4,istate) = ra(3)*ra(3)
         spectrum(5,istate) = ra(1)*ra(1)+ra(2)*ra(2)+ra(3)*ra(3)
      END DO 
    END IF

  END SUBROUTINE spectrum_dip_len


!***************************************************************************
!!****f* xas_methods/spectrum_dip_vel [1.0] *
!!
!!   NAME
!!    spectrum_dip_vel 
!!
!!   FUNCTION
!!     Calculation of the spectrum when the dipole approximation
!!     in the velocity form is used. 
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

  SUBROUTINE spectrum_dip_vel(xas_control, fm_set,op_sm,mos,excvec,fm_work,spectrum,estate,error)

    TYPE(xas_control_type)                      :: xas_control
    TYPE( cp_fm_p_type ), DIMENSION(:,:),  &
      POINTER                                   :: fm_set
    TYPE( real_matrix_p_type ), DIMENSION(:,:), &
      POINTER                                   :: op_sm
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(cp_fm_type), POINTER                   :: excvec
    TYPE(cp_fm_type), POINTER                   :: fm_work
    REAL(dp), DIMENSION(:,:), POINTER           :: spectrum
    INTEGER, INTENT(IN)                         :: estate

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "spectrum_dip_vel", &
      routineP = module_name//':'//routineN

    INTEGER :: homo, i, istate,istate0, j, nao, nvirtual, nmo
    LOGICAL :: failure
    REAL(dp) :: ene_i, ene_f,  dip(3)
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
      istate0 = homo
      IF(xas_control%emission) istate0 = 0
      DO i=1,SIZE(fm_set,2)
         CPPrecondition(ASSOCIATED(fm_set(1,i)%matrix),cp_failure_level,routineP,error,failure)
         CALL cp_fm_set_all(fm_set(1,i)%matrix, 0.0_dp, error=error)
         CALL cp_sm_fm_multiply(op_sm(1,i)%matrix,mo_coeff,fm_work,ncol=nmo)
         CALL cp_fm_gemm("T","N",1,nvirtual,nao,1.0_dp,excvec,&
              fm_work,0.0_dp, fm_set(1,i)%matrix,b_first_col=istate0+1,error=error)
      END DO
      
      ene_i = eigenvalues(estate)
      DO istate = 1,nvirtual
         ene_f = eigenvalues(istate0+istate)
         DO i = 1,3
           CALL cp_fm_get_element(fm_set(1,i)%matrix,1,istate,dip(i))
         END DO 
         spectrum(1,istate) = ene_f - ene_i
         spectrum(2,istate) = dip(1)*dip(1)
         spectrum(3,istate) = dip(2)*dip(2)
         spectrum(4,istate) = dip(3)*dip(3)
         spectrum(5,istate) = dip(1)*dip(1)+dip(2)*dip(2)+dip(3)*dip(3)
      END DO 
    END IF

  END SUBROUTINE spectrum_dip_vel

END MODULE xas_methods
