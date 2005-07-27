
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_linres_methods *
!!
!!   NAME 
!!     qs_linres_methods
!!
!!   FUNCTION
!!     localize wavefunctions
!!     linear response scf
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     created 07-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_methods

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                             cp_fm_symm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release, &
                                             cp_fm_struct_p_type,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: global_environment_type 
  USE kinds,                           ONLY: dp
  USE input_constants,                 ONLY: do_loc_none,state_loc_all,&
                                             op_loc_berry, lr_diis, lr_pcg
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: linres_control_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_type,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_types,                    ONLY: qs_loc_env_release,&
                                             qs_loc_env_retain,&
                                             qs_loc_env_create,&
                                             qs_loc_env_destroy,&
                                             qs_loc_env_new_type
  USE qs_linres_nmr_op,                ONLY: nmr_operators
  USE qs_loc_utils,                    ONLY: qs_loc_env_init
  USE qs_matrix_pools,                 ONLY: mpools_get, qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_p_env_methods,                ONLY: p_postortho,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_tddfpt_eigensolver,           ONLY: apply_op
  USE qs_tddfpt_utils,                 ONLY: normalize
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE


  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: linres_localize, linres_scf

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_methods'

!!***
! *****************************************************************************

CONTAINS
                                      
! *****************************************************************************
!!****f* cp2k/qs_linres_methods/linres_localize
!!
!!   NAME
!!     linres_localize
!!
!!   FUNCTION
!!     Find the centers and spreads of the wfn, 
!!     if required apply a localization algorithm
!!
!!   NOTES
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE linres_localize(localized_wfn_control,qs_env,nspins,globenv,centers_only,error)

    TYPE(localized_wfn_control_type), POINTER   :: localized_wfn_control
    TYPE(qs_environment_type), POINTER          :: qs_env
    INTEGER, INTENT(IN)                         :: nspins
    LOGICAL, INTENT(IN), OPTIONAL               :: centers_only
    TYPE(global_environment_type),  &
      POINTER                                   :: globenv
    TYPE(cp_error_type), & 
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'boh', &
                                   routineP = moduleN//'/'//routineN

    INTEGER :: ispin, istat, n_mo(2), nelectron, nmoloc(2)
    LOGICAL :: failure, ionode, my_centers_only
    REAL(dp) :: maxocc
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(qs_loc_env_new_type), POINTER          :: qs_loc_env     

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%source==logger%para_env%mepos
    my_centers_only = .FALSE.
    IF(PRESENT(centers_only)) my_centers_only = centers_only
    nmoloc(1:2) = 0
    NULLIFY(mos,qs_loc_env)

    IF(my_centers_only) THEN
       localized_wfn_control%set_of_states = state_loc_all
       localized_wfn_control%localization_method = do_loc_none
       localized_wfn_control%operator_type = op_loc_berry
    END IF

    CALL get_qs_env(qs_env=qs_env, mos=mos)
    DO ispin = 1,nspins
       CALL get_mo_set(mos(ispin)%mo_set,nmo=n_mo(ispin),nelectron=nelectron,maxocc=maxocc)
       IF(localized_wfn_control%set_of_states == state_loc_all) THEN
         nmoloc(ispin) = NINT(nelectron/maxocc)
       ELSE
         nmoloc(ispin) = MIN(localized_wfn_control%nloc_states(1) ,n_mo(ispin))
       END IF
    END DO  ! ispin
    CALL set_loc_wfn_lists(localized_wfn_control,&
         nmoloc,n_mo,nspins,error=error)
    CALL set_loc_centers(localized_wfn_control,nmoloc,&
         nspins,globenv,error=error) 
    CALL qs_loc_env_create(qs_loc_env,error=error)
    CALL qs_loc_env_init(qs_loc_env,localized_wfn_control,qs_env,error=error)
    ! The localized orbitals are copied in the mos, 
    ! The orbital centers are stored in linres_control%localized_wfn_control
    CALL qs_loc_driver(qs_env,qs_loc_env,ionode,error=error)
    CALL qs_loc_env_destroy(qs_loc_env, error=error)
    DEALLOCATE(qs_loc_env,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE linres_localize
  
! *****************************************************************************
!!****f* cp2k/qs_linres_methods/linres_scf
!!
!!   NAME
!!     linres_scf
!!
!!   FUNCTION
!!     scf loop to optimize the first order wavefunctions (psi1)
!!     given a perturbation as an operator applied to the ground state orbitals (h1_psi0) 
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE linres_scf(p_env, qs_env, psi1, h1_psi0, error)
 
    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: psi1, h1_psi0
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_scf', &
                                   routineP = moduleN//'/'//routineN


    INTEGER ::  handle, i, ik, istat, ispin, iter, last_ik, nvectors, out_unit
    LOGICAL                                     :: ionode, failure
    TYPE(cp_blacs_env_type), POINTER            :: blacs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: R, X
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                   :: b, Ab
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                   :: ao_mo_fm_pools
    TYPE(cp_fm_struct_p_type), DIMENSION(:), &
      POINTER                                   :: kv_fm_struct
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(qs_matrix_pools_type), POINTER         :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s
 

    failure   = .FALSE.
    CALL timeset(routineN,"I","",handle)

     NULLIFY(ao_mo_fm_pools, blacs_env, dft_control, linres_control, matrix_s, mpools, para_env, R, X)
     NULLIFY(Ab, b, kv_fm_struct, logger)

    logger    => cp_error_get_logger(error)
    ionode    = (logger%para_env%mepos==logger%para_env%source)
    out_unit  = cp_logger_get_default_unit_nr(logger)


     CALL get_qs_env(qs_env=qs_env,  blacs_env=blacs_env, &
                     dft_control= dft_control, &
                     linres_control=linres_control,&
                     matrix_s=matrix_s,&
                     mpools=mpools, &
                     para_env=para_env)



     !----------------------!
     ! allocate the vectors !
     !----------------------!
     CALL mpools_get(mpools, ao_mo_fm_pools=ao_mo_fm_pools)
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools, X, name=routineP//":X")
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools, R, name=routineP//":R")
   
     ALLOCATE(kv_fm_struct(dft_control%nspins), STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     
     DO ispin=1, dft_control%nspins
        NULLIFY (kv_fm_struct(ispin)%struct)
        CALL cp_fm_struct_create(kv_fm_struct(ispin)%struct, para_env, &
             blacs_env, p_env%n_ao(ispin), p_env%n_mo(ispin))
     END DO

     ! The storage of older vectors and residues depends on the optimization algorithm (DIIS,PCG..)
     nvectors = 1
     ALLOCATE(b(nvectors,dft_control%nspins), Ab(nvectors,dft_control%nspins), STAT=istat)
!              Sb(nvectors,dft_control%nspins), stat=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     DO ispin=1,dft_control%nspins
       DO i=1, nvectors
          NULLIFY(b(i,ispin)%matrix, Ab(i,ispin)%matrix )
       END DO
     END DO
     ! Copy the initial guess
     DO ispin=1, dft_control%nspins
       CALL cp_fm_to_fm(psi1(ispin)%matrix, b(1,ispin)%matrix)
     END DO

     ik = 1
     last_ik = 1
     DO  iter = 1,linres_control%max_scf

     !Orthogonalize with respect to the psi0
     CALL p_preortho(p_env, qs_env, b(ik,:))
     CALL normalize(b(ik,:), R, matrix_s) ! R is temp

     DO ispin=1, dft_control%nspins
        CALL cp_fm_to_fm(b(ik,ispin)%matrix, X(ispin)%matrix)
     END DO
     ! Apply the operators that do not depend on the perturbationnn
     CALL apply_op(X, Ab(ik,:), p_env, qs_env, linres_control%do_kernel)
     ! Apply the perturbative operator
!     CALL apply_op_h1(Ab(ik,:),h1_psi0, p_env, qs_env)

     !Refine the orthogonalization constraint
     CALL p_postortho(p_env,qs_env, Ab(ik,:))

     !Optimizer
     IF(linres_control%opt_method == lr_diis) THEN

     ELSEIF(linres_control%opt_method == lr_pcg) THEN

     END IF

     !Check Convergence
     IF(iter /= 1) THEN

     END IF


     END DO  ! iter

    ! Copy the optimized wavefunctions
     DO ispin=1, dft_control%nspins
       CALL cp_fm_to_fm(b(1,ispin)%matrix,psi1(ispin)%matrix)
     END DO

    !----------!
    ! clean up !
    !----------!
     DO ispin=1, dft_control%nspins
        CALL cp_fm_struct_release(kv_fm_struct(ispin)%struct)
        DO i=1, nvectors
           IF (ASSOCIATED(b(i,ispin)%matrix)) &
                CALL cp_fm_release(b(i,ispin)%matrix)
           IF (ASSOCIATED(Ab(i,ispin)%matrix)) &
                CALL cp_fm_release(Ab(i,ispin)%matrix)
!           IF (ASSOCIATED(Sb(i,ispin)%matrix)) &
!                CALL cp_fm_release(Sb(i,ispin)%matrix)
        END DO
     END DO
     DEALLOCATE(b, Ab, kv_fm_struct, STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, X)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, R)
     NULLIFY(X, R, b, Ab, kv_fm_struct)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE linres_scf

END MODULE qs_linres_methods
