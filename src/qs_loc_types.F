!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_loc_types [1.0] *
!!
!!   NAME
!!     qs_loc_types
!!
!!   FUNCTION
!!     New version of the module for the localization of the molecular orbitals
!!     This should be able to use different definition of the spread functional
!!     It should also calculate the integrals analytically so that it can be
!!     used irrespective of the pw_env and the collocation of wfn on the grids
!!     It should also work with a selected set of states, instead than all of them,
!!     in this case one should check that the selected states have the same occupation number
!!     The spread functional can be only estimated, or also optimized by minimization
!!     and in principle also maximization should be available.
!!     This operations can be required irrespective of the printing requirements
!!     It would be highly desirable to do all this along a MD run every N steps,
!!     and have a trajectory of the centeroids of the localized wfn
!!     In addition these functions can be used for properties calculations 
!!     like NMR and XAS. Therefore it is necessary that the rotated wfn are then copied
!!     in the mos fm matrix to be available for next use.
!!
!!   AUTHORS
!!     MI (05-2005)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_loc_types

! *****************************************************************************

  USE ai_moments,                      ONLY: contract_cossin,&
                                             cossin
  USE atomic_kind_types,               ONLY: atomic_kind_type, &
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE coefficient_types,               ONLY: coeff_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE external_potential_types,        ONLY: gth_potential_type, &
                                             get_potential
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE input_constants,                 ONLY: do_loc_min,&
                                             do_loc_max,&
                                             do_loc_jacobi,&
                                             do_loc_crazy,&
                                             do_loc_direct,&
                                             do_loc_none,&
                                             op_loc_berry,&
                                             op_loc_boys,&
                                             op_loc_pipek,&
                                             state_loc_all,&
                                             state_loc_upton,&
                                             state_loc_fromitoj,&
                                             state_loc_list
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_release,&
                                             localized_wfn_control_retain,&
                                             localized_wfn_control_type
  USE qs_localization_methods,         ONLY: initialize_weights
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             qlist_type,&
                                             block_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next, &
                                             reduced_3c_list_type
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE cp_files,                        ONLY: open_file, close_file
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             get_matrix_info
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: module_name = "qs_loc_types"


! *** Public ***
  PUBLIC :: qs_loc_env_create, qs_loc_env_destroy, qs_loc_env_init,&
            get_qs_loc_env, set_qs_loc_env
  PUBLIC :: qs_loc_env_new_type

!****t* qs_loc_types/qs_loc_env_new_type [1.0] *
!!
!!   NAME
!!     qs_loc_env_new_type
!!
!!   FUNCTION
!!     contains all the info needed by quickstep to calculate
!!     the spread of a selectedset of orbitals and if required
!!     to minimize or mazimize the spread by rotation of the orbitals
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     is_initialized : logical control variable
!!     para_env: info for the distribution of the calculations
!!     mo_coeff: full matrix containing only the selected subset of orbitals
!!     local_molecules: molecules distributed 
!!     cell : box that contains the system
!!     localized_wfn_control: variables and parameter that define the spread 
!!                            functional and the optimization algorithm
!!     particle_set: position, type, ao_indexes etc for each atom
!!     op_sm_set: set of sparse matrices used to define the spread operator
!!                 when the functional is defined by the use operator acting on the
!!                 basis functions, e.g. the Berry phase definition
!!                 The matrix element of the type <a|O|b> are computed in initialization
!!                 of qs_loc_env
!!     op_fm_set: set of full matrices used to define the spread operator
!!                when the functional has to be defined directly using the products of MOS
!!                as in the case of the Pipek-Mezek definition.
!!     weights: for a spread defined as extension of the orbitral in the box, these
!!              factors renormilize with respect to the box size 
!!get_qs_loc_env
!!   NOTE
!!     this type should replace the previous set up for the localization of the wfn
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04-05 created
!!
!!*** *********************************************************************


  TYPE qs_loc_env_new_type

       INTEGER :: is_initialized
       LOGICAL :: molecular_states
       TYPE ( cp_para_env_type ), POINTER       :: para_env
       TYPE ( cp_fm_p_type), DIMENSION(:),  POINTER :: moloc_coeff
       TYPE ( cp_fm_p_type), DIMENSION(:,:),&
         POINTER                                :: op_fm_set
       TYPE(distribution_1d_type), POINTER      :: local_molecules
       TYPE ( cell_type ), POINTER              :: cell
       TYPE ( localized_wfn_control_type ), &
         POINTER                                :: localized_wfn_control
       TYPE (particle_type), DIMENSION(:),&
         POINTER                                :: particle_set 
       TYPE (real_matrix_p_type), DIMENSION(:,:),&
         POINTER                                :: op_sm_set

       REAL (KIND = dp)                         :: weights ( 6 )
  END TYPE qs_loc_env_new_type
  

! *****************************************************************************

CONTAINS

!******************************************************************************
!!****f* qs_loc_types/qs_loc_env_create[1.0] *

  SUBROUTINE qs_loc_env_create(qs_loc_env,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'qs_loc_env_create', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: istat
    LOGICAL                                    :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(qs_loc_env),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN
       ALLOCATE(qs_loc_env,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       NULLIFY(qs_loc_env%para_env)
       NULLIFY(qs_loc_env%cell)
       NULLIFY(qs_loc_env%op_sm_set)
       NULLIFY(qs_loc_env%op_fm_set)
       NULLIFY(qs_loc_env%local_molecules)
       NULLIFY(qs_loc_env%moloc_coeff)
       NULLIFY(qs_loc_env%particle_set)
       NULLIFY(qs_loc_env%localized_wfn_control)
       qs_loc_env%weights = 0.0_dp

    END IF

  END SUBROUTINE qs_loc_env_create

!******************************************************************************
!!****f* qs_loc_types/qs_loc_env_destroy [1.0] *

  SUBROUTINE qs_loc_env_destroy(qs_loc_env,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'qs_loc_env_destroy', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: i, istat, j
    LOGICAL                                    :: failure

    failure =.FALSE.
    CPPrecondition(ASSOCIATED(qs_loc_env),cp_warning_level,routineP,error,failure)

    IF(.NOT. failure) THEN
       IF(ASSOCIATED(qs_loc_env%cell)) CALL cell_release(qs_loc_env%cell,error=error)
       IF(ASSOCIATED(qs_loc_env%local_molecules)) &
          CALL distribution_1d_release(qs_loc_env%local_molecules,error=error)
       IF (ASSOCIATED(qs_loc_env%localized_wfn_control)) THEN
         CALL localized_wfn_control_release(qs_loc_env%localized_wfn_control,&
                                            error=error)
       END IF
       IF(ASSOCIATED(qs_loc_env%para_env)) CALL cp_para_env_release(qs_loc_env%para_env)
       IF(ASSOCIATED(qs_loc_env%particle_set))  NULLIFY(qs_loc_env%particle_set)

       IF(ASSOCIATED(qs_loc_env%moloc_coeff)) THEN
         DO i=1,SIZE ( qs_loc_env % moloc_coeff,1) 
           CALL cp_fm_release(qs_loc_env%moloc_coeff(i)%matrix,error=error)
         END DO
         DEALLOCATE(qs_loc_env%moloc_coeff,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF(ASSOCIATED(qs_loc_env%op_fm_set)) THEN
         DO i=1,SIZE ( qs_loc_env % op_fm_set,2) 
           DO j=1,SIZE ( qs_loc_env % op_fm_set,1)
             CALL cp_fm_release(qs_loc_env%op_fm_set(j,i)%matrix,error=error)
           END DO
         END DO
         DEALLOCATE(qs_loc_env%op_fm_set,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF(ASSOCIATED(qs_loc_env%op_sm_set)) THEN
         DO i=1,SIZE ( qs_loc_env % op_sm_set, 2 )
           DO j=1,SIZE ( qs_loc_env % op_sm_set, 1 )
             CALL deallocate_matrix(qs_loc_env%op_sm_set(j,i)%matrix)
           ENDDO
         END DO
         DEALLOCATE(qs_loc_env%op_sm_set,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
       
    END IF 

  END SUBROUTINE qs_loc_env_destroy

!******************************************************************************
!!****f* qs_loc_types/get_qs_loc_env [1.0] *

  SUBROUTINE get_qs_loc_env(qs_loc_env,cell,local_molecules,&
             localized_wfn_control,moloc_coeff,&
             op_sm_set,op_fm_set,para_env,particle_set,weights,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cell_type), OPTIONAL, POINTER         :: cell
    TYPE(cp_fm_p_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: moloc_coeff
    TYPE(cp_fm_p_type), DIMENSION(:,:),&
      POINTER, OPTIONAL                        :: op_fm_set
   TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                  :: local_molecules
    TYPE(localized_wfn_control_type), &
      POINTER, OPTIONAL                        :: localized_wfn_control
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      OPTIONAL, POINTER                        :: op_sm_set
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                  :: para_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: particle_set
    REAL(dp), DIMENSION(6), OPTIONAL           :: weights
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'get_qs_loc_env', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: istat
    LOGICAL                                    :: failure

     failure =.FALSE.
     CPPrecondition(ASSOCIATED(qs_loc_env),cp_failure_level,routineP,error,failure)
 
     IF(.NOT.failure) THEN
       IF (PRESENT(cell)) cell =>  qs_loc_env%cell
       IF (PRESENT(moloc_coeff)) moloc_coeff => qs_loc_env%moloc_coeff
       IF (PRESENT(local_molecules)) local_molecules => qs_loc_env%local_molecules
       IF (PRESENT(localized_wfn_control)) &
          localized_wfn_control => qs_loc_env%localized_wfn_control
       IF (PRESENT(op_sm_set)) op_sm_set =>  qs_loc_env%op_sm_set
       IF (PRESENT(op_fm_set)) op_fm_set =>  qs_loc_env%op_fm_set
       IF (PRESENT(para_env)) para_env =>  qs_loc_env%para_env
       IF (PRESENT(particle_set)) particle_set =>  qs_loc_env%particle_set
       IF (PRESENT(weights)) weights(1:6)=  qs_loc_env%weights(1:6)
     END IF
 
  END SUBROUTINE get_qs_loc_env

!******************************************************************************
!!****f* qs_loc_types/set_qs_loc_env [1.0] *

  SUBROUTINE set_qs_loc_env(qs_loc_env,cell,local_molecules,&
             localized_wfn_control,moloc_coeff,&
             op_sm_set,op_fm_set,para_env,particle_set,weights,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cell_type), OPTIONAL, POINTER         :: cell
    TYPE(cp_fm_p_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: moloc_coeff
    TYPE(cp_fm_p_type), DIMENSION(:,:),&
      POINTER, OPTIONAL                        :: op_fm_set
   TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                  :: local_molecules
    TYPE(localized_wfn_control_type), &
      POINTER, OPTIONAL                        :: localized_wfn_control
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      OPTIONAL, POINTER                        :: op_sm_set
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                  :: para_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: particle_set
    REAL(dp), DIMENSION(6), OPTIONAL           :: weights
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error


    CHARACTER(len=*), PARAMETER :: routine_name = 'set_qs_loc_env', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: istat
    LOGICAL                                    :: failure

     failure =.FALSE.
     CPPrecondition(ASSOCIATED(qs_loc_env),cp_failure_level,routineP,error,failure)

     IF(.NOT.failure) THEN
     IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(qs_loc_env%cell,error=error)
       qs_loc_env%cell => cell
     END IF

     IF (PRESENT(local_molecules)) THEN
        CALL distribution_1d_retain(local_molecules,error=error)
        IF(ASSOCIATED(qs_loc_env%local_molecules)) &
             CALL distribution_1d_release(qs_loc_env%local_molecules,error=error)
        qs_loc_env%local_molecules => local_molecules
     END IF

     IF(PRESENT(localized_wfn_control)) THEN
       CALL localized_wfn_control_retain(localized_wfn_control,error=error)
       CALL localized_wfn_control_release(qs_loc_env%localized_wfn_control,error=error)
       qs_loc_env % localized_wfn_control => localized_wfn_control
     END IF
     IF(PRESENT(para_env)) THEN
      CALL cp_para_env_retain(para_env)
      CALL cp_para_env_release(qs_loc_env%para_env)
      qs_loc_env%para_env => para_env
     END IF
     IF (PRESENT(particle_set)) qs_loc_env%particle_set => particle_set
     IF(PRESENT(moloc_coeff)) THEN
       qs_loc_env%moloc_coeff => moloc_coeff
     END IF
     IF(PRESENT(op_sm_set)) THEN
       qs_loc_env%op_sm_set => op_sm_set
     END IF
     IF(PRESENT(op_fm_set)) THEN
       qs_loc_env%op_fm_set => op_fm_set
     END IF
     IF(PRESENT(weights)) THEN
       qs_loc_env%weights = weights
     END IF
     END IF
   
  END SUBROUTINE set_qs_loc_env
!******************************************************************************
!!****f* qs_loc_types/qs_loc_env_init [1.0] *
!!
!!   NAME
!!     qs_loc_env_init
!!
!!   FUNCTION
!!      allocates the data, and initializes the operators
!!
!!   NOTES
!!     similar to the old one, but not quite
!!
!!   ARGUMENTS
!!     - qs_loc_env: new environment for the localization calculations
!!     - localized_wfn_control: variables and directives for the localization 
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!*** **********************************************************************
  SUBROUTINE qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE ( localized_wfn_control_type ), &
      POINTER                                  :: localized_wfn_control
    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'qs_loc_env_init', &
      routineP = module_name//':'//routine_name

    INTEGER  :: handle, i, iatom, imo, imoloc, isos, ispin, istat, j, &
                lb, nao, naosub, natoms,  nmo, nmosub, nspins, ub
    LOGICAL  :: failure
    REAL(KIND=dp)                               :: my_occ, occ_imo
    REAL(KIND=dp), DIMENSION(:, :), POINTER     :: vecbuffer
    REAL(KIND=dp), DIMENSION(:), POINTER        :: occupations
    TYPE(cell_type), POINTER                    :: cell
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                   :: ao_mosub_fm_pools,&
                                                   mosub_mosub_fm_pools
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER   :: moloc_coeff
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(cp_para_env_type), POINTER             :: para_env
    TYPE(distribution_1d_type), POINTER         :: local_molecules
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                   :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s

!   ---------------------------------------------------------------------------

    CALL timeset("qs_loc_env_init","I",'',handle)

    failure = .FALSE.
    NULLIFY(mos,matrix_s,moloc_coeff,ao_mosub_fm_pools,mosub_mosub_fm_pools)

  
    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, cell=cell, &
                     local_molecules=local_molecules,&
                     particle_set=particle_set, &
                     para_env=para_env,&
                     mos=mos)

    CALL mpools_get(qs_env%mpools,ao_mosub_fm_pools=ao_mosub_fm_pools,&
                    mosub_mosub_fm_pools=mosub_mosub_fm_pools)

    nspins = SIZE(mos,1)
    ALLOCATE (moloc_coeff(nspins), STAT=istat)
    DO ispin = 1,nspins
      CALL get_mo_set(mos(ispin)%mo_set,nao=nao,nmo=nmo)
      CALL fm_pool_create_fm(ao_mosub_fm_pools(ispin)%pool,moloc_coeff(ispin)%matrix,&
           name="qs_loc_env%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)
      CALL cp_fm_get_info(moloc_coeff(ispin)%matrix, nrow_global=naosub, &
           ncol_global=nmosub,error=error)
      CPPostcondition(nao==naosub,cp_failure_level,routineP,error,failure)
      CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
      CALL cp_fm_set_all(moloc_coeff(ispin)%matrix,0.0_dp)

    END DO  ! ispin

    ! Copy the submatrix

    DO ispin = 1,nspins
      CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,&
           occupation_numbers=occupations, nao=nao, nmo=nmo)
      lb = localized_wfn_control%lu_bound_states(1,ispin)
      ub = localized_wfn_control%lu_bound_states(2,ispin)

      IF(localized_wfn_control%set_of_states == state_loc_list) THEN
        ALLOCATE(vecbuffer(1,nao),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        nmosub = SIZE(localized_wfn_control%loc_states,1)
        CPPostcondition(nmosub>0,cp_failure_level,routineP,error,failure)
        imoloc = 0
        my_occ =  occupations(localized_wfn_control%loc_states(1))
        DO i = lb,ub
           ! Get the index in the subset
           imoloc = imoloc + 1
           ! Get the index in the full set
           imo = localized_wfn_control%loc_states(i)
           occ_imo = occupations(imo)
           IF(ABS(occ_imo-my_occ) > EPSILON(1.0_dp)) THEN
             IF(localized_wfn_control%localization_method /= do_loc_none) &
               CALL stop_program(routineP,"States with different occupations"//&
               " cannot be rotated together")
           END IF
           ! Take the imo vector from the full set and copy in the imoloc vector of the subset
           CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,imo,&
                nao,1,transpose=.TRUE.,error=error)
           CALL cp_fm_set_submatrix(moloc_coeff(ispin)%matrix,vecbuffer,1,imoloc,&
                nao,1,transpose=.TRUE.,error=error)
         END DO
         DEALLOCATE(vecbuffer,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ELSE
         my_occ = occupations( lb )
         occ_imo = occupations( ub )
         IF(ABS(occ_imo-my_occ) >  EPSILON(1.0_dp)) THEN
            IF(localized_wfn_control%localization_method /= do_loc_none) &
             CALL stop_program(routineP,"States with different occupations"//&
             " cannot be rotated together")
         END IF
         nmosub = localized_wfn_control%nloc_states(ispin)
         CALL cp_fm_to_fm(mo_coeff,moloc_coeff(ispin)%matrix,nmosub,lb,1)
       END IF
    END DO


    CALL set_qs_loc_env(qs_loc_env=qs_loc_env, cell=cell, local_molecules=local_molecules,&
                        moloc_coeff=moloc_coeff, particle_set=particle_set,para_env=para_env,&
                        localized_wfn_control=localized_wfn_control)

    ! Prepare the operators

    IF ( localized_wfn_control%operator_type==op_loc_berry )  THEN
      ALLOCATE ( qs_loc_env % op_sm_set ( 2, 3 ), STAT = istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
      DO i=1, SIZE ( qs_loc_env % op_sm_set, 2 )
         DO j=1, SIZE ( qs_loc_env % op_sm_set, 1 )
           NULLIFY(qs_loc_env%op_sm_set(j,i)%matrix)
           CALL replicate_matrix_structure(matrix_s(1)%matrix, &
                qs_loc_env%op_sm_set(j,i)%matrix,"qs_loc_env%op_sm_"//&
                TRIM(ADJUSTL(cp_to_string(j)))//"-"//TRIM(ADJUSTL(cp_to_string(i))))
           CALL set_matrix(qs_loc_env%op_sm_set(j,i)%matrix,0.0_dp)
         ENDDO
      END DO

      ALLOCATE ( qs_loc_env % op_fm_set( 2, 3 ), STAT = istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i =   1, SIZE ( qs_loc_env % op_fm_set, 2 )
         DO j =   1, SIZE ( qs_loc_env % op_fm_set, 1 )
           NULLIFY(qs_loc_env % op_fm_set(j,i)%matrix)
           CALL fm_pool_create_fm(mosub_mosub_fm_pools(1)%pool,&
               qs_loc_env %op_fm_set(j,i)%matrix,&
               name="qs_loc_env%zij_fm_"//TRIM(ADJUSTL(cp_to_string(j)))//&
               "-"//TRIM(ADJUSTL(cp_to_string(i))),error=error)
           CALL cp_fm_get_info(qs_loc_env %op_fm_set(j,i)%matrix, nrow_global=nmosub, &
               error=error)
           CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(qs_loc_env %op_fm_set(j,i)%matrix,0.0_dp)
         END DO
      END DO
    ELSEIF ( localized_wfn_control%operator_type== op_loc_pipek) THEN
      natoms = SIZE ( qs_loc_env % particle_set, 1 )
!      ALLOCATE ( qs_loc_env % op_fm_set ( natoms , nspins ), STAT = istat )
      ALLOCATE ( qs_loc_env % op_fm_set ( natoms , 1 ), STAT = istat )
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
      DO ispin = 1, SIZE(qs_loc_env % op_fm_set ,2)
         CALL get_mo_set(mos(ispin)%mo_set,nmo=nmo)
         DO iatom = 1,natoms
           NULLIFY(qs_loc_env %op_fm_set(iatom,ispin)%matrix)
           CALL fm_pool_create_fm(mosub_mosub_fm_pools(1)%pool,&
                qs_loc_env %op_fm_set(iatom,ispin)%matrix,&
                name="qs_loc_env%op_fm_"//TRIM(ADJUSTL(cp_to_string(iatom)))//&
                "-"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)
           CALL cp_fm_get_info(qs_loc_env %op_fm_set(iatom,ispin)%matrix, nrow_global=nmosub, &
                error=error)
           CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
           CALL cp_fm_set_all(qs_loc_env %op_fm_set(iatom,ispin)%matrix,0.0_dp)
         END DO  ! iatom
      END DO  ! ispin
    ELSE
      CALL stop_program(routineP,"Type of operator not implemented")
    ENDIF


    IF (  localized_wfn_control%operator_type==op_loc_berry ) THEN 
 
      CALL initialize_weights ( qs_loc_env % cell, qs_loc_env % weights ) 
  
      CALL get_berry_operator ( qs_loc_env, qs_env , error=error)

    ELSEIF ( localized_wfn_control%operator_type== op_loc_pipek) THEN

!!    here we don't have to do anything
!!    CALL get_pipek_mezey_operator ( qs_loc_env, qs_env )

    ENDIF

    
    qs_loc_env%molecular_states = .FALSE.
  

! done
    qs_loc_env%is_initialized = 666
    CALL timestop(0.0_dp,handle) 

END SUBROUTINE qs_loc_env_init

!******************************************************************************
!!****f* qs_loc_types/get_berry_operator [1.0] *
!!
!!   NAME
!!     get_berry_operator
!!
!!   FUNCTION
!!      Computes the Berry operator for periodic systems
!!      used to define the spread of the MOS
!!      Here the matrix elements of the type <mu|cos(kr)|nu> and  <mu|sin(kr)|nu> 
!!      are computed, where mu and nu are the contracted basis functions.
!!      Namely the Berry operator is exp(ikr)
!!      k is defined somewhere
!!      the pair lists are exploited and sparse matrixes are constructed
!!
!!   NOTES
!!     The intgrals are computed analytically  using the primitives GTO
!!     The contraction is performed block-wise
!!
!!   ARGUMENTS
!!     - qs_loc_env: new environment for the localization calculations
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE get_berry_operator ( qs_loc_env, qs_env, error ) 

! globals
    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

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
    LOGICAL :: failure, new_atom_b, ortho
    REAL(KIND=dp)                              :: ooLx, ooLy, ooLz, rab2, dab
    REAL(KIND=dp)                              :: twopoLx, twopoLy, twopoLz
    REAL(KIND=dp), DIMENSION(3)                :: kvec, ra, rb, rab
    REAL(KIND=dp), DIMENSION(6)                :: weights
    REAL(KIND = dp), DIMENSION(:), POINTER     :: set_radius_a,set_radius_b
    REAL(KIND=dp), DIMENSION(:,:), POINTER     :: rpgfa, rpgfb, sphi_a,sphi_b, zeta, zetb 
    REAL(KIND=dp), DIMENSION(:,:), POINTER     :: cosab,sinab, work 
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER   :: cosabt,sinabt, workt
    TYPE(atomic_kind_type), DIMENSION(:),&
      POINTER                                  :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER            :: atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
          POINTER                              :: op_cos, op_sin
    TYPE(block_p_type), DIMENSION(:,:), &
      POINTER                                  :: op_cost, op_sint
    TYPE(cell_type), POINTER                   :: cell
    TYPE(gto_basis_set_type), POINTER          :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:),           POINTER          :: sab_orb
    TYPE(neighbor_list_type), POINTER          :: sab_orb_neighbor_list,&
                                                  sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER          :: sab_orb_neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                  :: particle_set
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      POINTER                                  :: op_sm_set

!$  INTEGER :: omp_get_max_threads,omp_get_thread_num

    CALL timeset(routine_name,"I"," ",handle)

    failure = .FALSE.
    NULLIFY (atomic_kind, atomic_kind_set)
    NULLIFY ( cell, sab_orb, op_sm_set , orb_basis_set, particle_set)
    NULLIFY (cosabt,sinabt, workt)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set, sab_orb=sab_orb)

    nkind = SIZE(atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=ldwork, maxlgto=maxl )
    ldab = ldwork
    nthread = 1
!$  nthread = omp_get_max_threads()
    CALL reallocate(cosabt,1,ldab,1,ldab,0,nthread-1) 
    CALL reallocate(sinabt,1,ldab,1,ldab,0,nthread-1) 
    CALL reallocate(workt,1,ldwork,1,ldwork,0,nthread-1)

    ALLOCATE (op_cost(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (op_sint(3,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, 3
      DO j=0,nthread-1
        NULLIFY (op_cost(i,j)%block)
        NULLIFY (op_sint(i,j)%block)
      END DO
    END DO 
 
    CALL get_qs_loc_env(qs_loc_env=qs_loc_env,op_sm_set=op_sm_set, cell=cell)
    ooLx = 1.0_dp/SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) ) 
    ooLy = 1.0_dp/SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )
    ooLz = 1.0_dp/SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )
    twopoLx = twopi*ooLx
    twopoLy = twopi*ooLy
    twopoLz = twopi*ooLz
    ortho = cell%orthorhombic

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
!$OMP private(ithread,cosab,sinab,work,ilist,iatom,jatom)&
!$OMP private(last_jatom,new_atom_b,istat,sgfa,sgfb,ra,rb)&
!$OMP private(nrow,iset,jset,ncoa,ncob,nnode,inode,irow,icol)&
!$OMP private(op_cos,op_sin,rab,rab2,dab,i,kvec)&
!$OMP private(sab_orb_neighbor_node,sab_orb_neighbor_list_local)&
!$OMP shared(ikind,jkind,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldab,ldwork)&
!$OMP shared(lb_max,npgfb,nsetb,iab,la_min,rpgfa,zeta)&
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ldsa,ldsb)&
!$OMP shared(nkind,sab_orb_neighbor_list,nlist,workt,cosabt,sinabt,qs_env)&
!$OMP shared(op_cost,op_sint,ortho,twopoLx,twopoLy,twopoLz)&
!$OMP shared(first_sgfa,first_sgfb,sab_orb,lb_min,la_max,op_sm_set)
!$OMP shared(atomic_kind_set,particle_set,cell)
        ithread = 0
!$      ithread = omp_get_thread_num()

        cosab => cosabt(:,:,ithread)
        sinab => sinabt(:,:,ithread)
        work  => workt(:,:,ithread)
        op_cos => op_cost(:,ithread)
        op_sin => op_sint(:,ithread)

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
                 NULLIFY(op_cos(i)%block)
                 CALL get_block_node(matrix=op_sm_set(1,i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_cos(i)%block)
                 NULLIFY(op_sin(i)%block)
                 CALL get_block_node(matrix=op_sm_set(2,i)%matrix,&
                                     block_row=irow,&
                                     block_col=icol,&
                                     block=op_sin(i)%block)
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

!           *** Calculate the primitive overlap integrals ***
                  IF(ortho)THEN

                    ! x component
                    kvec = 0.0_dp
                    kvec(1) = twopoLx
                    cosab = 0.0_dp  
                    sinab = 0.0_dp  

                    CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              ra,rb,kvec,cosab,sinab)

                    CALL contract_cossin(op_cos(1)%block,op_sin(1)%block,&
                               iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                               jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                               cosab,sinab,ldab,work,ldwork)

                    ! y component
                    kvec = 0.0_dp
                    kvec(2) = twopoLy
                    cosab = 0.0_dp  
                    sinab = 0.0_dp  


                    CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              ra,rb,kvec,cosab,sinab)

                    CALL contract_cossin(op_cos(2)%block,op_sin(2)%block,&
                               iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                               jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                               cosab,sinab,ldab,work,ldwork)

                    ! z component
                    kvec = 0.0_dp
                    kvec(3) = twopoLz
                    cosab = 0.0_dp  
                    sinab = 0.0_dp  


                    CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                              lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                              ra,rb,kvec,cosab,sinab)

                    CALL contract_cossin(op_cos(3)%block,op_sin(3)%block,&
                               iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                               jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                               cosab,sinab,ldab,work,ldwork)

                  ELSE
                    CALL stop_program(routineP,&
                      "spread calculation with berry phase in non orthorombic cell NYI")
                  END IF  ! ortho

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
         NULLIFY(op_cost(i,j)%block)
         NULLIFY(op_sint(i,j)%block)
      END DO
    END DO 
    DEALLOCATE(op_cost,op_sint, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(cosabt,sinabt,workt,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE get_berry_operator


!******************************************************************************
!
! Computes the Pipek-Mezey operator for periodic systems
!
!******************************************************************************
!!SUBROUTINE get_pipek_mezey_operator ( qs_loc_env, qs_env )
!!! globals
!!    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
!!    TYPE(qs_environment_type), POINTER       :: qs_env
!!
!!  CALL stop_program("qs_localization","PIPEK-MEZEY not implemented")
!!END SUBROUTINE get_pipek_mezey_operator

!******************************************************************************
!******************************************************************************
!
! maps wfc's to molecules 
!
!******************************************************************************
!!SUBROUTINE wfc_to_molecule ( qs_loc_env, center, molecule_set ) 
!!    TYPE(qs_loc_env_type), INTENT(IN)        :: qs_loc_env
!!    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
!!    TYPE(molecule_type), POINTER             :: molecule_set( : )
!!
!!    CHARACTER(LEN=default_string_length)     :: mol_name
!!    INTEGER :: counter, first_atom, i, iatom, ikind, imol, imol_now, isos, &
!!      istate, local_location, natom, natom_loc, natom_max, nkind, nmol, nstate
!!    INTEGER, POINTER                         :: wfc_to_atom_map( : )
!!    REAL(KIND=dp)                            :: dr(3), mydist(2), r_pbc( 3 ), dipole(3)
!!    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: dipole_set
!!    REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE :: charge_set
!!    REAL(KIND=dp), POINTER                   :: distance( : ), r( :, : )
!!    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
!!    TYPE(cp_para_env_type), POINTER          :: para_env
!!    TYPE(distribution_1d_type), POINTER      :: local_molecules
!!    TYPE(molecule_kind_type), POINTER        :: molecule_kind
!!    TYPE(particle_type), POINTER             :: particle_set( : )
!!    TYPE(atomic_kind_type), POINTER          :: atomic_kind
!!    TYPE(gth_potential_type), POINTER        :: gth_potential
!!    INTEGER                                  :: zeff,output_unit
!!
!!!---------------------------------------------------------------------------
!!!---------------------------------------------------------------------------
!!
!!  particle_set => qs_loc_env % particle_set
!!  para_env => qs_loc_env % para_env
!!  local_molecules => qs_loc_env % local_molecules
!!  nstate = SIZE ( center, 1 )
!!
!!
!!  ALLOCATE(wfc_to_atom_map(nstate),STAT=isos)
!!  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!!!---------------------------------------------------------------------------
!!!---------------------------------------------------------------------------
!!  nkind = SIZE (local_molecules%n_el)
!!  natom = 0
!!  natom_max = 0
!!  DO ikind=1,nkind
!!    nmol = SIZE ( local_molecules%list(ikind)%array )
!!    DO imol=1,nmol
!!      i = local_molecules%list(ikind)%array(imol)
!!!      WRITE(*,*) para_env%mepos,"ABC",ikind,imol,i
!!      molecule_kind => molecule_set(i)%molecule_kind
!!      CALL get_molecule_kind(molecule_kind=molecule_kind,&
!!                             natom=natom)
!!      natom_max = natom_max + natom
!!      molecule_set(i)%lmi%nstates=0
!!      IF(ASSOCIATED(molecule_set(i)%lmi%states)) THEN
!!        DEALLOCATE(molecule_set(i)%lmi%states,STAT=isos)
!!        IF(isos/=0) CALL stop_memory('wfc_to_molecule','states')
!!      END IF
!!    END DO
!!  END DO
!!  natom_loc = natom_max
!!  natom = natom_max
!!
!!  CALL mp_max(natom_max,para_env%group)
!!
!!!  WRITE(*,*) para_env%mepos,"nkind  = ",nkind
!!!  WRITE(*,*) para_env%mepos,"natom_loc  = ",natom_loc
!!!  WRITE(*,*) para_env%mepos,"natom_max  = ",natom_max
!!!  WRITE(*,*) para_env%mepos,"nstate = ",nstate
!!
!!  ALLOCATE(r(3,natom_max),STAT=isos)
!!  IF(isos/=0) CALL stop_memory('centers_and_spreads','r')
!!
!!  ALLOCATE(distance(natom_max),STAT=isos)
!!  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')
!!
!!  !Zero all the stuff
!!  r(:,:) = 0.0_dp
!!  distance(:) = 1.E10_dp
!!
!!!---------------------------------------------------------------------------
!!!---------------------------------------------------------------------------
!!  counter = 0
!!  nkind = SIZE (local_molecules%n_el)
!!  DO ikind=1,nkind
!!    nmol = SIZE ( local_molecules%list(ikind)%array )
!!    DO imol=1,nmol
!!      i = local_molecules%list(ikind)%array(imol)
!!      molecule_kind => molecule_set(i)%molecule_kind
!!      first_atom = molecule_set(i)%first_atom
!!      CALL get_molecule_kind(molecule_kind=molecule_kind,&
!!                             natom=natom,&
!!                             name=mol_name,&
!!                             atom_list=atom_list)
!!      DO iatom=1,natom
!!        counter = counter + 1
!!        r ( :, counter ) = particle_set(first_atom+iatom-1) % r ( : )
!!      END DO
!!    END DO
!!  END DO
!!
!!  !DO iatom=1,natom_loc
!!  !  WRITE(*,*) iatom,r(:,iatom)
!!  !END DO
!!
!!!---------------------------------------------------------------------------
!!!---------------------------------------------------------------------------
!!  DO istate=1,nstate
!!    distance(:) = 1.E10_dp
!!    DO iatom=1,natom_loc
!!      dr(1) = r(1,iatom) - center(istate,1)
!!      dr(2) = r(2,iatom) - center(istate,2)
!!      dr(3) = r(3,iatom) - center(istate,3)
!!      r_pbc = pbc ( dr, qs_loc_env%cell )
!!      distance(iatom) = SQRT ( DOT_PRODUCT ( r_pbc, r_pbc ) )
!!    END DO
!!
!!    !combine distance() from all procs
!!    local_location = MINLOC(distance,DIM=1)
!!    
!!    mydist(1) = distance(local_location)
!!    mydist(2) = para_env%mepos
!!
 !!   CALL mp_minloc(mydist,para_env%group)
 !!   
!!    IF(mydist(2)==para_env%mepos) THEN
!!      wfc_to_atom_map(istate)  =  local_location
!!    ELSE
!!      wfc_to_atom_map(istate)  =  0
!!    END IF
!!  END DO
!!!---------------------------------------------------------------------------
!!!---------------------------------------------------------------------------
!!  IF(natom_loc/=0) THEN
!!    DO istate = 1,nstate
!!      iatom = wfc_to_atom_map(istate)
!!      IF(iatom/=0) THEN
!!        counter = 0
!!        nkind = SIZE (local_molecules%n_el)
!!        DO ikind=1,nkind
!!          nmol = SIZE ( local_molecules%list(ikind)%array )
!!          DO imol=1,nmol
!!!!!!            imol_now = local_molecules%list(ikind)%array(imol)
!!            molecule_kind => molecule_set(imol_now)%molecule_kind
!!            CALL get_molecule_kind(molecule_kind=molecule_kind,&
!!                                   natom=natom)
!!            counter = counter + natom
!!            IF(counter>=iatom) EXIT
!!          END DO
!!          IF(counter>=iatom) EXIT
 !!       END DO
!! !!      
!!!        WRITE(*,*) "PP",para_env%mepos,istate,iatom,imol_now
!!     
!!        i = molecule_set(imol_now)%lmi%nstates
!!        i = i + 1
!!        molecule_set(imol_now)%lmi%nstates = i
!!        CALL reallocate(molecule_set(imol_now)%lmi%states,1,i)
!!        molecule_set(imol_now)%lmi%states(i) = istate
!!      END IF
!!    END DO
!!  END IF
!!
!!!---------------------------------------------------------------------------
!!! figure out dipole of the molecule. 
! Assume it is neutral (so origin doesn't matter, we take the first atom)
! if the molecule is larger than half the box there might be a problem folding stuff
! FIXME this should be moved out of this place
!---------------------------------------------------------------------------
!!  particle_set => qs_loc_env % particle_set
!!  para_env => qs_loc_env % para_env
!!  local_molecules => qs_loc_env % local_molecules
!!  nstate = SIZE ( center, 1 )
!!  ALLOCATE(dipole_set(3,SIZE(molecule_set)))
!!  ALLOCATE(charge_set(SIZE(molecule_set)))
!!  dipole_set=0.0_dp
!!  charge_set=0.0_dp
!!  nkind = SIZE (local_molecules%n_el)
!!  DO ikind=1,nkind ! loop over different molecules
!!     nmol = SIZE ( local_molecules%list(ikind)%array )
!!     DO imol=1,nmol ! all the molecules of the kind
!!        imol_now = local_molecules%list(ikind)%array(imol) ! index in the global array
!!        IF (.NOT.ASSOCIATED(molecule_set(imol_now)%lmi%states)) CYCLE
 !!       molecule_kind => molecule_set(imol_now)%molecule_kind
 !!!!       first_atom = molecule_set(imol_now)%first_atom
 !!       CALL get_molecule_kind(molecule_kind=molecule_kind,&
 !!                              natom=natom,atom_list=atom_list)
 !!       dipole=0.0_dp
 !!       ! charges of the atoms involved
 !!       DO iatom=1,natom
 !!           i=first_atom+iatom-1
 !!           atomic_kind=>particle_set(i)%atomic_kind
 !!           CALL get_atomic_kind(atomic_kind,gth_potential=gth_potential)
  !!          CALL get_potential(gth_potential,zeff=zeff)
  !!!!          r_pbc = particle_set(i) % r - particle_set(first_atom) % r
  !!          r_pbc = pbc ( r_pbc , qs_loc_env%cell )
  !!          dipole=dipole + zeff*r_pbc
  !!          charge_set(imol_now)=charge_set(imol_now)+zeff
!!        END DO
!!        ! charges of the wfc involved
!!        DO istate=1,SIZE(molecule_set(imol_now)%lmi%states) 
!!           i = molecule_set(imol_now)%lmi%states(istate)
!!           r_pbc = center(i,1:3) - particle_set(first_atom) % r 
!!           r_pbc = pbc ( r_pbc , qs_loc_env%cell )
!!           dipole=dipole - 2.0_dp * r_pbc ! XXXXXXXXXXXXXXXXXXXX double occupation FIXME XXXXXXXXX (where are the spins?)
!!           charge_set(imol_now)=charge_set(imol_now)-2.0_dp
!!        ENDDO
!!        dipole_set(:,imol_now)=dipole ! a.u.
 !!    ENDDO
 !! ENDDO
 !!!! CALL mp_sum(dipole_set,para_env%group)
 !! CALL mp_sum(charge_set,para_env%group)
 !! IF (para_env%mepos==para_env%source) THEN
 !!    CALL open_file(file_name="MOLECULAR_DIPOLES",&
 !!                       file_action="WRITE",&
 !!                       file_status="UNKNOWN",&
 !!                       file_position="APPEND",&
 !!                       unit_number=output_unit)
 !!    WRITE(UNIT=output_unit,FMT='(A40)') "# molecule nr, charge, dipole vector, dipole (Debye) "
 !!    dipole_set=dipole_set*2.541747_dp ! Debye
 !!    DO I=1,SIZE(dipole_set,2)
 !!       WRITE(UNIT=output_unit,FMT='(I6,5F12.6)') I,charge_set(I),dipole_set(1:3,I), &
 !!                                             SQRT(DOT_PRODUCT(dipole_set(1:3,I),dipole_set(1:3,I)))
 !!    ENDDO
  !!   CALL close_file(output_unit)
  !!ENDIF
 !! DEALLOCATE(dipole_set,charge_set)
!!!---------------------------------------------------------------------------
! end of molecular dipole calculation
!---------------------------------------------------------------------------

!  nkind = SIZE (local_molecules%n_el)
!  DO ikind=1,nkind
!    nmol = SIZE ( local_molecules%list(ikind)%array )
!    DO imol=1,nmol
!      i = local_molecules%list(ikind)%array(imol)
!      WRITE(*,*)"XX",para_env%mepos,i,molecule_set(i)%lmi%states
!    END DO
!  END DO
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!!  DEALLOCATE(distance,STAT=isos)
!!  IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')
!!  DEALLOCATE(r,STAT=isos)
!!  IF(isos/=0) CALL stop_memory('wfc_to_molecule','r')

!!  DEALLOCATE(wfc_to_atom_map,STAT=isos)
!!  IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  
  !CALL stop_program("wfc_to_atom","checkpoint")

!!END SUBROUTINE wfc_to_molecule
!------------------------------------------------------------------------------

END MODULE qs_loc_types

