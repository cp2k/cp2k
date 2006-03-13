!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
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

  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: allocate_sto_basis_set,&
                                             create_gto_from_sto_basis,&
                                             deallocate_sto_basis_set,&
                                             get_gto_basis_set,&
                                             gto_basis_set_type,&
                                             init_orb_basis_set,&
                                             set_sto_basis_set,&
                                             sto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_files,                        ONLY: open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
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
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       do_loc_none, op_loc_berry, state_loc_upton, xas_1s_type, xas_2p_type, &
       xas_2s_type, xas_dip_len, xas_dscf, xas_scf_default, xas_scf_general, &
       xas_scf_ot, xas_tddft, xas_tp_fh, xas_tp_hh, xas_tp_xfh, xas_tp_xhh, &
       xes_tp_val
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_loc_control,                  ONLY: localized_wfn_control_release,&
                                             localized_wfn_control_retain,&
                                             localized_wfn_control_type,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_loc_methods,                  ONLY: qs_loc_driver,&
                                             qs_print_cubes
  USE qs_loc_types,                    ONLY: qs_loc_env_create,&
                                             qs_loc_env_new_type,&
                                             qs_loc_env_release
  USE qs_loc_utils,                    ONLY: qs_loc_env_init
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_set,&
                                             write_mo_set
  USE qs_operators_ao,                 ONLY: p_xyz_ao
  USE qs_scf,                          ONLY: init_scf_loop,&
                                             init_scf_run,&
                                             scf_env_cleanup
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE semi_empirical_parameters,       ONLY: srules
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE string_utilities,                ONLY: xstring
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             set_xas_env,&
                                             xas_env_create,&
                                             xas_env_release,&
                                             xas_environment_type
  USE xas_tp_scf,                      ONLY: xas_do_tp_scf
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_methods'  

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

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(xas_control_type)                   :: xas_control
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, i, iat, iatom, ikind, isgf, ispin, istat, &
      istate, iter, j, my_homo(2), my_kind, nao, natom, nelectron, &
      nexc_search, nmo, nspins, nvirtual2, output_unit, state_to_be_excited, &
      xas_iter_count
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf
    INTEGER, DIMENSION(:), POINTER           :: atom_of_state, &
                                                mykind_of_kind, &
                                                state_of_atom, &
                                                state_of_mytype, type_of_state
    LOGICAL                                  :: failure, ionode, my_uocc(2), &
                                                transition_potential, uocc
    REAL(dp)                                 :: component, dist, distmin, &
                                                maxocc, ra(3), rac(3), rc(3)
    REAL(dp), DIMENSION(:), POINTER          :: all_evals, eigenvalues, &
                                                max_overlap, &
                                                occupation_numbers, &
                                                sto_state_overlap
    REAL(dp), DIMENSION(:, :), POINTER       :: centers_wfn, spectrum, &
                                                vecbuffer
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: stogto_overlap
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: groundstate_coeff
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dip_fm_set
    TYPE(cp_fm_type), POINTER                :: all_vectors, excvec_coeff, &
                                                fm_work, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(localized_wfn_control_type), &
      POINTER                                :: my_localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, op_sm
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: ostrength_sm
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, scf_section, &
                                                xas_section
    TYPE(xas_environment_type), POINTER      :: xas_env

    CALL timeset(routineN,"I",'',handle)

    failure = .FALSE.
    transition_potential = .FALSE.


    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)


    NULLIFY(xas_env, groundstate_coeff, ostrength_sm, op_sm, dip_fm_set, fm_work)
    NULLIFY(excvec_coeff, spectrum, state_of_atom, centers_wfn, qs_loc_env, cell)
    NULLIFY(occupation_numbers,my_localized_wfn_control,matrix_ks)
    NULLIFY(all_evals,all_vectors)
    NULLIFY(atom_of_state,type_of_state,state_of_atom,stogto_overlap)
    NULLIFY(mykind_of_kind, state_of_mytype,max_overlap,sto_state_overlap)
    NULLIFY(dft_section, xas_section)
    dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)
    xas_section => section_vals_get_subs_vals(dft_section,"XAS",error=error)

    output_unit = cp_print_key_unit_nr(logger,xas_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".Log",error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "START XES/XAS CALCULATION",&
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
      CALL xas_env_init(xas_env, xas_control, qs_env, globenv, error=error)
      CALL set_qs_env(qs_env,xas_env=xas_env, error=error)
      CALL xas_env_release(xas_env,error=error)
      CALL get_qs_env(qs_env,xas_env=xas_env,error=error)
    END IF

!   Initialize the type of calculation
    NULLIFY(atomic_kind_set, cell, dft_control, scf_control, mos, para_env, particle_set)
    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
         cell = cell, dft_control=dft_control, scf_control=scf_control,&
         matrix_ks=matrix_ks,mos=mos, para_env=para_env, particle_set=particle_set )


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
             para_env=para_env, &
             do_rotation=.TRUE.)
      END DO
    END IF
    DO ispin = 1,dft_control%nspins
        CALL get_mo_set(mos(ispin)%mo_set,nelectron=nelectron,maxocc=maxocc,&
             homo=my_homo(ispin),uniform_occupation=my_uocc(ispin))
        my_homo(ispin) = homo
        my_uocc(ispin) = uocc
      IF(scf_control%smear /= 0.0_dp) THEN
        homo = NINT(nelectron/maxocc)
        uocc = .TRUE.
        CALL set_mo_set(mos(ispin)%mo_set, homo=homo,uniform_occupation=uocc)
      END IF
    END DO

    nexc_search = xas_control%nexc_search
    nspins = dft_control%nspins
    transition_potential = (xas_control%xas_method.LT.xas_tddft).OR.&
                           (xas_control%xas_method==xas_dscf)
    IF(nspins==1 .AND. transition_potential) THEN
       CALL stop_program(routineP,"xas with tp method requires LSD calculations")
    END IF
!   Set of states among which there is the state to be excited
    CALL get_mo_set(mos(1)%mo_set,nao=nao,homo=homo)
    IF(nexc_search < 0) nexc_search = homo

    CALL get_xas_env(xas_env=xas_env,iter_count=xas_iter_count,&
         centers_wfn=centers_wfn,atom_of_state=atom_of_state,&
         mykind_of_kind=mykind_of_kind,&
         type_of_state=type_of_state, state_of_atom=state_of_atom,&
         groundstate_coeff=groundstate_coeff,excvec_coeff=excvec_coeff,&
         fm_work=fm_work,stogto_overlap=stogto_overlap)

    xas_iter_count =  xas_iter_count + 1
    CALL set_xas_env(xas_env=xas_env,iter_count=xas_iter_count)
    
   !I know from which atoms I want to excite (O) and which type of states (1s) (input)

   !Define the qs_loc_env : to find centers, spread and possibly localize them
    my_localized_wfn_control =>  xas_control%localized_wfn_control
    CALL localized_wfn_control_retain(my_localized_wfn_control,error=error)
    CPPrecondition(ASSOCIATED(my_localized_wfn_control),cp_failure_level,routineP,error,failure)

    CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)

    IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,/,T10,A,I7,/)")&
           "Localize a set of MOs: those where to search the states to be excited",&
           "MOs to be localized from 1 to ",nexc_search
    END IF
    IF(.NOT.(ASSOCIATED(qs_loc_env))) THEN
      CALL qs_loc_env_create(qs_loc_env,error=error)
      CALL set_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
      CALL qs_loc_env_release(qs_loc_env,error=error)
      CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
      CALL qs_loc_env_init(qs_loc_env,my_localized_wfn_control,qs_env,onlyspin1=.TRUE.,error=error)
    END IF
    CALL qs_loc_driver(qs_env,qs_loc_env,onlyspin1=.TRUE.,error=error)

    ! scratch array for the state
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, maxocc=maxocc)
    natom = SIZE(particle_set)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL get_particle_set(particle_set=particle_set, first_sgf=first_sgf)
    ALLOCATE (op_sm(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
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
      centers_wfn(1,istate) =  my_localized_wfn_control%centers_set(1)%array(1,istate)
      centers_wfn(2,istate) =  my_localized_wfn_control%centers_set(1)%array(2,istate)
      centers_wfn(3,istate) =  my_localized_wfn_control%centers_set(1)%array(3,istate)

 !if(ionode) WRITE(*,'(A,I4,3f10.5)') 'center ', istate, centers_wfn(1:3, istate)

      ! Assign the state to the closest atom
      distmin = 100.0_dp

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
! WRITE(*,'(I4,A,I4,6f10.5)') istate, 'atom ', atom_of_state(istate), &
!         particle_set(atom_of_state(istate))%r(1:3), centers_wfn(1:3,istate)
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
!      write(*,*)  istate, iatom, sto_state_overlap(istate) ,  state_of_mytype(iatom)
      END IF
    END DO  ! istate

    CALL localized_wfn_control_release(my_localized_wfn_control,error=error)

  ! If required check the character of the states of each atom
  ! by taking the overlap with atomic-like states

  ! In the set of states, assign the index of the state to be excited for iatom
    IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,/)")&
           "List the atoms to be excited and the relative of MOs index "
    END IF
    DO iat = 1,xas_control%nexc_atoms
      iatom = xas_control%exc_atoms(iat)
     state_of_atom(iat) = state_of_mytype(iatom)
     IF (output_unit>0) THEN
       WRITE(UNIT=output_unit,FMT="(T10,A,I6,T32,A,I6)") 'Atom: ',iatom ,&
       "localized MO index", state_of_atom(iat)
     END IF
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
      ! Take the state_to_be_excited vector from the full set and copy into excvec_coeff
       CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,&
            occupation_numbers=occupation_numbers,homo=homo,nmo=nmo,nelectron=nelectron)
       IF(xas_control%xas_method==xas_dscf .OR. xas_control%xas_method==xas_tp_xhh &
         .OR.  xas_control%xas_method==xas_tp_xfh) THEN
          homo = nelectron + 1
          CPPrecondition(.NOT.(nmo<homo),cp_failure_level,routineP,error,failure)
          CALL set_mo_set(mos(1)%mo_set, homo = homo)
       END IF
       occupation_numbers(1:homo)  = 1.0_dp
       CPPrecondition(ASSOCIATED(excvec_coeff),cp_failure_level,routineP,error,failure)
       CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,state_to_be_excited,&
            nao,1,transpose=.TRUE.,error=error)
       CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
            nao,1,transpose=.TRUE.,error=error)

      IF (transition_potential) THEN

         IF (output_unit>0) THEN
           WRITE(UNIT=output_unit,FMT='(/,/,T10,A,I6)') &
                  "START XAS CALCULATION WITH TP APPROACH FOR ATOM   ", iatom 
         END IF
         NULLIFY(scf_env)
         CALL get_qs_env(qs_env,scf_env=scf_env,error=error)

         CALL init_scf_run(scf_env=scf_env,qs_env=qs_env,&
                 globenv=globenv, scf_section=scf_section, error=error)
         ! this is not standing optimally here, but for he time being retains functionality
         CALL init_scf_loop(scf_env=scf_env,qs_env=qs_env,&
                 globenv=globenv, scf_section=scf_section, error=error)

         CALL xas_do_tp_scf(xas_control,xas_env,iatom,scf_env,qs_env,globenv,error=error)



!         CALL qs_loc_driver(qs_env,qs_loc_env,onlyspin1=.TRUE.,error=error)
!dbg 
!  write(*,*) iat, state_to_be_excited
! stop 'check'
         IF(iat == xas_control%nexc_atoms) THEN 
            CALL scf_env_cleanup(scf_env,qs_env=qs_env,error=error)
         END IF

!   *** Write last wavefunction to screen ***
         DO ispin=1,dft_control%nspins
           CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,&
           4,6,dft_section,error=error)
         ENDDO

         CALL get_xas_env(xas_env=xas_env,all_vectors=all_vectors,nvirtual2=nvirtual2,&
              all_evals=all_evals,dip_fm_set=dip_fm_set,spectrum=spectrum)

         IF(.NOT. xas_control%xas_method == xas_dscf) THEN
         ! Calculate the spectrum
           IF(xas_control%dipole_form == xas_dip_len) THEN
             CALL spectrum_dip_len(xas_control,dip_fm_set,ostrength_sm,mos,excvec_coeff,&
                  all_vectors,all_evals,fm_work,&
                  cell,spectrum,dft_control%xas_estate,error=error)
           ELSE
             CALL spectrum_dip_vel(xas_control,dip_fm_set,ostrength_sm,mos,excvec_coeff,&
                  all_vectors,all_evals,fm_work,&
                  spectrum,dft_control%xas_estate,error=error)
           END IF
         END IF
      ELSE
         ! the states defined by the rotation are the ground state orbitals
         ! the initial state from which I excite should be localized
         ! I take the excitations from lumo to nmo
!         CALL xas_tddft(error)
      END IF

      !update iteration counter
      iter = xas_iter_count*xas_control%each_step

    ! append the spectrum of this iteration to the right output file
      IF( xas_control%xas_method == xas_dscf) THEN
!         homo = homo - 1
!         CALL set_mo_set(mos(1)%mo_set, homo = homo)
      ELSE
        CALL xas_write(xas_control,spectrum, dft_control%xas_estate, &
             homo, iatom, iter, ionode, error=error)
      END IF

      IF( xas_control%cubes) THEN
        CALL xas_print_cubes(xas_control,qs_env,mos,all_vectors,iter,iatom,ionode,error=error)
      END IF

      DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
            nao=nao, nmo=nmo)
!       CALL cp_fm_to_fm(groundstate_coeff(ispin)%matrix,mo_coeff,nmo,1,1)
! not sure that orthogonality is conserved
       CALL cp_fm_to_fm(groundstate_coeff(ispin)%matrix,mo_coeff,nexc_search,1,1)
      END DO

    END DO

    DO ispin = 1,dft_control%nspins
       CALL set_mo_set(mos(ispin)%mo_set, homo=my_homo(ispin),&
            uniform_occupation=my_uocc(ispin))
    END DO

  ! Release what has to be released 

    DEALLOCATE (first_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(op_sm,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (sto_state_overlap,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (max_overlap,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (state_of_mytype,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
   

    dft_control % xas_estate    = -1

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,T3,A,/,T25,A,/,T3,A,/)")&
            REPEAT("=",77),&
            "END XES/XAS CALCULATION",&
            REPEAT("=",77)
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,xas_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)
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
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: homo, i, iat, iatom, ik, ikind, ispin, istat, j, l, n_mo(2), &
      nao, nexc_atoms, nexc_search, nj, nk, nkind, nmo, nmoloc(2), norb(0:3), &
      nsgf_gto, nsgf_sto, nspins, nwork
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_type_tmp, kind_z_tmp
    INTEGER, DIMENSION(4, 7)                 :: ne
    INTEGER, DIMENSION(:), POINTER           :: lq, nq
    LOGICAL                                  :: failure, ihavethis
    REAL(dp)                                 :: zatom
    REAL(dp), DIMENSION(:), POINTER          :: sto_zet
    REAL(dp), DIMENSION(:, :), POINTER       :: sto_alpha
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(localized_wfn_control_type), &
      POINTER                                :: my_localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(sto_basis_set_type), POINTER        :: sto_basis_set

!------------------------------------------------------------------------------

    failure=.FALSE.

    n_mo(1:2) = 0

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

      NULLIFY( atomic_kind_set, dft_control, scf_control, matrix_s, mos, mpools )
      NULLIFY( para_env, particle_set)
      CALL get_qs_env(qs_env=qs_env, &
              atomic_kind_set=atomic_kind_set, &
              dft_control = dft_control, &
              scf_control = scf_control, &
              mpools=mpools,& 
              matrix_s=matrix_s, mos=mos, &
              para_env=para_env, particle_set=particle_set)
      nexc_search = xas_control%nexc_search
      nexc_atoms = xas_control%nexc_atoms

      xas_env%eps_diis = xas_control%eps_diis
      IF(xas_env%eps_diis .LT. 0.0_dp) xas_env%eps_diis = scf_control%eps_diis
      xas_env%p_mix = xas_control%p_mix
      IF(xas_env%p_mix .LT. 0.0_dp) xas_env%p_mix = scf_control%p_mix
      xas_env%scf_method = xas_control%scf_method
      IF(xas_env%scf_method==xas_scf_default) THEN
        xas_env%scf_method = xas_scf_general
        IF(scf_control%use_ot) xas_env%scf_method = xas_scf_ot
      END IF

      CALL mpools_get(mpools, ao_mo_fm_pools= xas_env%ao_mo_fm_pools)

      NULLIFY(mo_coeff)
      CALL get_mo_set(mos(1)%mo_set,nao=nao,homo=homo,nmo=nmo,mo_coeff=mo_coeff)
      IF(nexc_search < 0) nexc_search = homo
      
      IF(xas_control%emission) THEN
         xas_env%nvirtual1 = nmo
      ELSE
         xas_env%nvirtual1 = nmo - homo
      END IF
      IF(xas_control%added_mos .GT. 0) THEN
         xas_env%nvirtual2 = MIN(xas_control%added_mos,nao-nmo)
         xas_env%unoccupied_eps = xas_control%eps_added
         xas_env%unoccupied_max_iter = xas_control%max_iter_added
      END IF
      xas_env%nvirtual = xas_env%nvirtual1 +  xas_env%nvirtual2

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

      ! initialize operators for the calculation of the oscillator strengts
      IF (xas_control%xas_method.LT.xas_tddft &
          .OR. xas_control%xas_method==xas_dscf) THEN
        
        IF (xas_control%xas_method==xas_tp_hh) THEN
          xas_control%occ_estate = 0.5_dp
          xas_control%occ_homo = 1.0_dp
        ELSEIF(xas_control%xas_method==xas_tp_xhh) THEN
          xas_control%occ_estate = 0.5_dp
          xas_control%occ_homo = 0.5_dp
        ELSEIF(xas_control%xas_method==xas_tp_fh) THEN
          xas_control%occ_estate = 0.0_dp
          xas_control%occ_homo = 1.0_dp
        ELSEIF(xas_control%xas_method==xas_tp_xfh) THEN
          xas_control%occ_estate = 0.0_dp
          xas_control%occ_homo = 1.0_dp
        ELSEIF(xas_control%xas_method==xes_tp_val) THEN
          xas_control%occ_estate = xas_control%xes_core_occupation
          xas_control%occ_homo = 1.0_dp-xas_control%xes_core_occupation
        ELSEIF(xas_control%xas_method==xas_dscf) THEN
          xas_control%occ_estate = 0.0_dp
          xas_control%occ_homo = 1.0_dp
        ENDIF

        NULLIFY(tmp_fm_struct)
        nwork = MAX(xas_env%nvirtual,nmo)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,ncol_global=nwork,&
             para_env=para_env,context=mo_coeff%matrix_struct%context)
        CALL cp_fm_create (xas_env%fm_work,tmp_fm_struct)
        CALL cp_fm_struct_release ( tmp_fm_struct )


        ALLOCATE (xas_env%groundstate_coeff(nspins), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DO ispin = 1,nspins
          CALL get_mo_set(mos(ispin)%mo_set,nao=nao,nmo=nmo)
          CALL fm_pool_create_fm(xas_env%ao_mo_fm_pools(ispin)%pool,&
               xas_env%groundstate_coeff(ispin)%matrix,&
             name="xas_env%mo0"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)
        END DO  ! ispin

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

        !Array to store all the eigenstates: occupied and the required not occupied
        IF(xas_env%nvirtual2 .GT. 0) THEN
          ALLOCATE(xas_env%unoccupied_evals(xas_env%nvirtual2), STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          NULLIFY(tmp_fm_struct)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
               ncol_global=xas_env%nvirtual2,&
               para_env=para_env,context=mo_coeff%matrix_struct%context)
          CALL cp_fm_create (xas_env%unoccupied_orbs,tmp_fm_struct)
          CALL cp_fm_struct_release ( tmp_fm_struct )

        END IF

        NULLIFY(tmp_fm_struct)
        CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao, &
             ncol_global=xas_env%nvirtual,&
             para_env=para_env,context=mo_coeff%matrix_struct%context)
        CALL cp_fm_create (xas_env%all_vectors,tmp_fm_struct)
        CALL cp_fm_struct_release ( tmp_fm_struct )

        ! Array to store all the energies needed  for the spectrum
        ALLOCATE(xas_env%all_evals(xas_env%nvirtual), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

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

      !Type of state
      norb = 0
      ALLOCATE(nq(1),lq(1),sto_zet(1),STAT=istat)
      IF( xas_control%state_type == xas_1s_type) THEN
        norb(0) = 1
        ALLOCATE(sto_alpha(1,0:1),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        nq(1) = 1
        lq(1) = 0
      ELSEIF( xas_control%state_type == xas_2s_type ) THEN
        norb(0) = 2
        ALLOCATE(sto_alpha(2,0:1),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        nq(1) = 2
        lq(1) = 0
      ELSEIF( xas_control%state_type == xas_2p_type ) THEN
        norb(0) = 2
        norb(1) = 1
        ALLOCATE(sto_alpha(2,0:1),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        nq(1) = 2
        lq(1) = 1
      ELSE
        CALL stop_program(routineP,"xas type of state not implemented")
      END IF

      ALLOCATE(kind_type_tmp(nkind),STAT=istat)
      ALLOCATE(kind_z_tmp(nkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      kind_type_tmp = 0
      kind_z_tmp = 0
      nk=0
      DO iat = 1,xas_control%nexc_atoms
        iatom = xas_control%exc_atoms(iat)
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
        sto_alpha = 0.0_dp
        ne = 0
        DO l=1,4 !lq(1)+1
          nj = 2*(l-1)+1
          DO i=l, 7! nq(1)
            ne(l,i) = ptable(kind_z_tmp(ik))%e_conv(l-1) - 2*nj*(i-l)
            ne(l,i) = MAX(ne(l,i),0)
            ne(l,i) = MIN(ne(l,i),2*nj)
          END DO
        END DO

         sto_alpha(nq(1),lq(1))=srules(kind_z_tmp(ik),ne,nq(1),lq(1))
        sto_zet(1) = sto_alpha(nq(1),lq(1))
        CALL allocate_sto_basis_set(sto_basis_set)
        CALL set_sto_basis_set(sto_basis_set,nshell=1,nq=nq,&
             lq=lq,zet=sto_zet)
        CALL create_gto_from_sto_basis(sto_basis_set,&
             xas_env%my_gto_basis(ik)%gto_basis_set,xas_control%ngauss)
        CALL deallocate_sto_basis_set(sto_basis_set)
        xas_env%my_gto_basis(ik)%gto_basis_set%norm_type = 2
        CALL init_orb_basis_set(xas_env%my_gto_basis(ik)%gto_basis_set)

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
      DEALLOCATE(sto_alpha,STAT=istat)
      DEALLOCATE(kind_type_tmp,kind_z_tmp,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF  ! failure

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

    TYPE(xas_control_type)                   :: xas_control
    REAL(dp), DIMENSION(:, :), POINTER       :: spectrum
    INTEGER, INTENT(IN)                      :: estate, homo, iatom, iter
    LOGICAL, INTENT(IN)                      :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=50)                        :: filename
    INTEGER                                  :: ia, ie, istate, istate0, &
                                                out_sp
    LOGICAL                                  :: failure, file_exists

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

      CLOSE(out_sp)
    END IF

  END SUBROUTINE xas_write

!***************************************************************************
!!****f*  xas_methods/xas_print_cubes [1.0] *
!!
!!   NAME
!!     xas_print_cubes
!!
!!   FUNCTION
!!     write the cube files for a set of selected states
!!
!!   ARGUMENTS
!!    - xas_control : provide number ant indexes of the states to be printed
!!    - qs_env
!!    - mos     : mos from which the states to be printed are extracted
!!    - iter    : iteration step (e.g. along an md)
!!    - iatom   : index of the atom that has been excited
!!    - ionode
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE  xas_print_cubes(xas_control,qs_env,mos,all_vectors,iter,iatom,ionode,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_type), POINTER                :: all_vectors
    INTEGER, INTENT(IN)                      :: iter, iatom
    LOGICAL, INTENT(IN)                      :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_print_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=50)                        :: flag
    INTEGER                                  :: homo, istat, istate0, nspins, &
                                                nstates
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: centers

    failure = .FALSE.

    nspins = SIZE(mos)

    flag="at_"//TRIM(ADJUSTL(cp_to_string(iatom)))
    flag=TRIM(flag)

    nstates = SIZE(xas_control%list_cubes,1)

    IF(xas_control%do_centers) THEN
      ! one might like to calculate the centers of the xas orbital (without localizing them)
    ELSE
      ALLOCATE(centers(6,nstates),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      centers = 0.0_dp
    END IF

    CALL get_mo_set(mos(1)%mo_set, homo=homo)
    istate0 = homo
    IF(xas_control%emission) istate0 = 0

    CALL qs_print_cubes(qs_env,all_vectors,nstates,xas_control%list_cubes,&
               centers,iter,1,ionode,flag,stride=xas_control%stride,&
               state0=istate0,error=error)

    DEALLOCATE(centers,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE xas_print_cubes

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

  SUBROUTINE spectrum_dip_len(xas_control,fm_set,op_sm,mos,excvec,&
             all_vectors,all_evals,fm_work,cell,spectrum,estate,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: fm_set
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: op_sm
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_type), POINTER                :: excvec, all_vectors
    REAL(dp), DIMENSION(:), POINTER          :: all_evals
    TYPE(cp_fm_type), POINTER                :: fm_work
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:, :), POINTER       :: spectrum
    INTEGER, INTENT(IN)                      :: estate
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'spectrum_dip_len', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: z
    INTEGER                                  :: homo, i, istate, j, nao, nmo, &
                                                nvirtual
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ene_f, ene_i, imagpart, &
                                                ra(3), realpart
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues

!    TYPE(cp_fm_type), POINTER                   :: mo_coeff

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(fm_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(spectrum),cp_failure_level,routineP,error,failure)
    nvirtual = SIZE(spectrum,2)
    NULLIFY(eigenvalues) 
    IF(.NOT. failure) THEN
      CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues,nao=nao, nmo=nmo, homo =homo)
!      istate0 = homo
!      IF(xas_control%emission) istate0 = 0
      DO i=1,SIZE(fm_set,2)
        DO j = 1,SIZE(fm_set,1)
           CPPrecondition(ASSOCIATED(fm_set(j,i)%matrix),cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(fm_set(j,i)%matrix, 0.0_dp, error=error)
           CALL cp_fm_set_all(fm_work, 0.0_dp, error=error)
           CALL cp_sm_fm_multiply(op_sm(j,i)%matrix,all_vectors,fm_work,ncol=nvirtual)
           CALL cp_fm_gemm("T","N",1,nvirtual,nao,1.0_dp,excvec,&
                  fm_work,0.0_dp, fm_set(j,i)%matrix,b_first_col=1,error=error)
        END DO 
      END DO
      
      ene_i = eigenvalues(estate)
      DO istate = 1,nvirtual
         ene_f = all_evals(istate)
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

  SUBROUTINE spectrum_dip_vel(xas_control, fm_set,op_sm,mos,excvec, &
             all_vectors,all_evals,fm_work,spectrum,estate,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: fm_set
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: op_sm
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_fm_type), POINTER                :: excvec, all_vectors
    REAL(dp), DIMENSION(:), POINTER          :: all_evals
    TYPE(cp_fm_type), POINTER                :: fm_work
    REAL(dp), DIMENSION(:, :), POINTER       :: spectrum
    INTEGER, INTENT(IN)                      :: estate
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'spectrum_dip_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: homo, i, istate, nao, nmo, &
                                                nvirtual
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dip(3), ene_f, ene_i
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(fm_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(spectrum),cp_failure_level,routineP,error,failure)
    nvirtual = SIZE(spectrum,2)
    NULLIFY(eigenvalues) 

    IF(.NOT. failure) THEN
      CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues,nao=nao, nmo=nmo, homo =homo)
!      istate0 = homo
!      IF(xas_control%emission) istate0 = 0
  
!dbg
!    write(*,*) 'nv ',nvirtual 
!dbg
      DO i=1,SIZE(fm_set,2)
         CPPrecondition(ASSOCIATED(fm_set(1,i)%matrix),cp_failure_level,routineP,error,failure)
         CALL cp_fm_set_all(fm_set(1,i)%matrix, 0.0_dp, error=error)
         CALL cp_fm_set_all(fm_work, 0.0_dp, error=error)
         CALL cp_sm_fm_multiply(op_sm(1,i)%matrix,all_vectors,fm_work,ncol=nvirtual)
         CALL cp_fm_gemm("T","N",1,nvirtual,nao,1.0_dp,excvec,&
              fm_work,0.0_dp, fm_set(1,i)%matrix,b_first_col=1,error=error)
      END DO
      
      ene_i = eigenvalues(estate)
!dbg
!    write(*,*) 'ene i ',ene_i
!    write(*,*) associated(all_evals)
!    write(*,*) size(all_evals,1)
!dbg
      DO istate = 1,nvirtual
         ene_f = all_evals(istate)
!dbg
!    write(*,*) 'ene f ', istate , ene_f
!dbg
         DO i = 1,3
           CALL cp_fm_get_element(fm_set(1,i)%matrix,1,istate,dip(i))
         END DO 
         spectrum(1,istate) = ene_f - ene_i
         spectrum(2,istate) = dip(1)*dip(1)
!         spectrum(1,istate) = dip(1)
         spectrum(3,istate) = dip(2)*dip(2)
!         spectrum(3,istate) = dip(2)
         spectrum(4,istate) = dip(3)*dip(3)
!         spectrum(4,istate) = dip(3)
         spectrum(5,istate) = dip(1)*dip(1)+dip(2)*dip(2)+dip(3)*dip(3)
!dbg
!    write(*,*) 'sp ',  dip(1), dip(3), dip(3)
!dbg
      END DO 
    END IF
!dbg
!    STOP 'pippo2'
!dbg
  END SUBROUTINE spectrum_dip_vel

!***************************************************************************

  SUBROUTINE calc_stogto_overlap(base_a, base_b, matrix,error)

    TYPE(gto_basis_set_type), POINTER        :: base_a, base_b
    REAL(dp), DIMENSION(:, :), POINTER       :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_stogto_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER :: iset, istat, jset, ldai, ldsab, maxcoa, maxcob, maxl, maxla, &
      maxlb, ncoa, ncob, nseta, nsetb, nsgfa, nsgfb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, &
                                                nsgfa_set, nsgfb_set
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rab(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: sab, work
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb

    failure = .FALSE.

    NULLIFY(la_max,la_min,lb_max,lb_min)
    NULLIFY(npgfa,npgfb,nsgfa_set,nsgfb_set)
    NULLIFY(first_sgfa,first_sgfb)
    NULLIFY(rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb)

    CALL get_gto_basis_set(gto_basis_set=base_a,nsgf=nsgfa,nsgf_set=nsgfa_set,lmax=la_max,&
                          lmin=la_min,npgf=npgfa,pgf_radius=rpgfa, &
                          sphi=sphi_a,zet=zeta,first_sgf=first_sgfa, &
                          maxco=maxcoa,nset=nseta, maxl=maxla)

    CALL get_gto_basis_set(gto_basis_set=base_b,nsgf=nsgfb,nsgf_set=nsgfb_set, lmax=lb_max,&
                          lmin=lb_min,npgf=npgfb,pgf_radius=rpgfb, &
                          sphi=sphi_b,zet=zetb,first_sgf=first_sgfb, &
                          maxco=maxcob,nset=nsetb,maxl=maxlb)
    rab = 0.0_dp
    ! Initialize and allocate
    matrix = 0.0_dp

    ldsab = MAX(maxcoa,maxcob,nsgfa,nsgfb)
    maxl = MAX(maxla,maxlb)
    ldai = ncoset(maxl)

! write(*,*) 'size ', size(matrix,1),size(matrix,2)
! write(*,*) 'dim ', ldsab, maxl, ldai

    ALLOCATE(sab(ldsab,ldsab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work(ldsab,ldsab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ai_work(ldai,ldai,1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

! write(*,*) 'nset ', nseta , nsetb

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

! write(*,*) 'iset ', iset, ncoa, sgfa, la_max(iset)

       DO jset=1,nsetb

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)

! write(*,*) 'jset ', jset, ncob, sgfb, lb_max(jset)
          ai_work = 0.0_dp

          CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
               rpgfa(:,iset),zeta(:,iset),&
               lb_max(jset),lb_min(jset),npgfb(jset),&
               rpgfb(:,jset),zetb(:,jset),rab,0.0_dp,sab,0,.FALSE.,ai_work,ldai)
! write(*,*) 'sab '
!    do ic = 1,ncoa
!      write(*,'(I4,10F12.6)') ic, sab(ic,1:ncob)
!    end do
!
! write(*,*) 'sphi_a '
!  do ic = sgfa, sgfa+nsgfa_set(iset) - 1
! write(*,'(10F12.6)') sphi_a(1:ncoa,ic)
! END DO 
! write(*,*) 'sphi_b '
!  do ic = sgfb, sgfb+nsgfb_set(jset) - 1
! write(*,'(10F12.6)') sphi_b(1:ncob,ic)
! END DO 

          CALL dgemm("N","N",ncoa,nsgfb_set(jset),ncob,&
               1.0_dp,sab(1,1),SIZE(sab,1),&
               sphi_b(1,sgfb),SIZE(sphi_b,1),&
               0.0_dp,work(1,1),SIZE(work,1))
          CALL dgemm("T","N",nsgfa_set(iset),nsgfb_set(jset),ncoa,&
               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
               work(1,1),SIZE(work,1),&
               1.0_dp,matrix(sgfa,sgfb), SIZE(matrix,1))

! write(*,'(A,4I5)') 'matrix', sgfa,sgfa+nsgfa_set(iset)-1, sgfb,sgfb+nsgfb_set(jset)-1
!      DO ic = sgfb,sgfb+nsgfb_set(jset)-1
!         write(*,'(10F12.6)')   matrix(sgfa:sgfa+nsgfa_set(iset)-1,ic)
!      END DO 

       END DO  ! jset
    END DO  ! iset


    DEALLOCATE(sab,work,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE calc_stogto_overlap

END MODULE xas_methods
