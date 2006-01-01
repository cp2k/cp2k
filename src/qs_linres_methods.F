
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_symm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release, &
                                             cp_fm_struct_p_type,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type,&
                                             cp_fm_to_fm
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: global_environment_type 
  USE kinds,                           ONLY: dp
  USE input_constants,                 ONLY: do_loc_none,state_loc_all,&
                                             op_loc_berry, lr_diis, lr_cg,&
                                             lr_sd,lr_precond_none,&
                                             lr_ls_none, lr_ls_2pnt, lr_ls_3pnt, lr_ls_gold
  USE preconditioner,                  ONLY: apply_preconditioner, &
                                             make_preconditioner
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
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

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
    LOGICAL                                     :: ionode, failure,energy_only
    REAL(dp) :: tmp
    REAL(dp), DIMENSION(:), POINTER             :: occupations
    TYPE(cp_blacs_env_type), POINTER            :: blacs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: R, X
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: b, Ab
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                   :: ao_mo_fm_pools
    TYPE(cp_fm_struct_p_type), DIMENSION(:), &
      POINTER                                   :: kv_fm_struct
    TYPE(cp_logger_type), POINTER               :: logger
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(dft_control_type), POINTER             :: dft_control
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER  :: mos
    TYPE(qs_matrix_pools_type), POINTER         :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_s
 

    failure   = .FALSE.
    energy_only = .FALSE.
    CALL timeset(routineN,"I","",handle)

     NULLIFY(ao_mo_fm_pools, blacs_env, dft_control, linres_control)
     NULLIFY(matrix_s, mos, mpools, para_env)
     NULLIFY(Ab, b, kv_fm_struct, logger,  R, X)

    logger    => cp_error_get_logger(error)
    ionode    = (logger%para_env%mepos==logger%para_env%source)
    out_unit  = cp_logger_get_default_unit_nr(logger)


     CALL get_qs_env(qs_env=qs_env,  blacs_env=blacs_env, &
                     dft_control= dft_control, &
                     linres_control=linres_control,&
                     matrix_s=matrix_s,&
                     mpools=mpools, mos=mos,&
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
     ALLOCATE(b(dft_control%nspins), Ab(dft_control%nspins), STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     DO ispin=1,dft_control%nspins
        NULLIFY(b(ispin)%matrix, Ab(ispin)%matrix )
     END DO
     ! Copy the initial guess
     DO ispin=1, dft_control%nspins
        CALL fm_pool_create_fm(ao_mo_fm_pools(ispin)%pool,&
               b(ispin)%matrix,&
               name="b"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)

        CALL cp_fm_to_fm(psi1(ispin)%matrix, b(ispin)%matrix)

        CALL fm_pool_create_fm(ao_mo_fm_pools(ispin)%pool,&
               Ab(ispin)%matrix,&
               name="Ab"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)
     END DO

     CALL check_p_env_init(p_env,linres_control, dft_control%nspins, error=error)

    !Orthogonalize with respect to the psi0
     CALL p_preortho(p_env, qs_env, b(:))
     CALL normalize(b(:), R, matrix_s) ! R is temp

     ik = 1
     last_ik = 1
     DO  iter = 1,linres_control%max_scf

       p_env%etotal = 0.0_dp
       p_env%gnorm = 0.0_dp
       p_env%gradient = 0.0_dp

       DO ispin=1, dft_control%nspins
          CALL cp_fm_to_fm(b(ispin)%matrix, X(ispin)%matrix)
       END DO
       ! Apply the operators that do not depend on the perturbation
       CALL apply_op(X, Ab(:), p_env, qs_env, linres_control%do_kernel)
       ! Apply the perturbative operator
       DO ispin=1, dft_control%nspins
          CALL cp_fm_scale_and_add(1.0_dp, Ab(ispin)%matrix, &
                                   1.0_dp, h1_psi0(ispin)%matrix) ! add h1_psi0 to Ab 
          CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                          occupation_numbers=occupations )
          CALL cp_fm_column_scale(Ab(ispin)%matrix, occupations)
          CALL cp_fm_trace(psi1(ispin)%matrix,Ab(ispin)%matrix,tmp,error)
          p_env%etotal = p_env%etotal + tmp 
       END DO

       !minimization
       CALL linres_scf_mini(p_env,qs_env,linres_control,b,Ab,&
            dft_control%nspins,iter,error=error)

       !Refine the orthogonalization constraint
  !     CALL p_postortho(p_env,qs_env, Ab(:))

       !Check Convergence
       IF(iter /= 1) THEN

       END IF


     END DO  ! iter

    ! Copy the optimized wavefunctions
     DO ispin=1, dft_control%nspins
       CALL cp_fm_to_fm(b(ispin)%matrix,psi1(ispin)%matrix)
     END DO

    !----------!
    ! clean up !
    !----------!
     DO ispin=1, dft_control%nspins
        CALL cp_fm_struct_release(kv_fm_struct(ispin)%struct)
        DO i=1, nvectors
           IF (ASSOCIATED(b(ispin)%matrix)) &
                CALL cp_fm_release(b(ispin)%matrix)
           IF (ASSOCIATED(Ab(ispin)%matrix)) &
                CALL cp_fm_release(Ab(ispin)%matrix)
        END DO
     END DO
     DEALLOCATE(b, Ab, kv_fm_struct, STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, X)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, R)
     NULLIFY(X, R, b, Ab, kv_fm_struct)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE linres_scf

  SUBROUTINE check_p_env_init(p_env,linres_control, nspins, error)

    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(linres_control_type), POINTER          :: linres_control
    INTEGER, INTENT(IN)                         :: nspins
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_p_env_init', &
                                   routineP = moduleN//'/'//routineN
    LOGICAL :: failure
    INTEGER :: ispin, nrow, ncol

    failure = .FALSE.

    p_env%iter = 0
    p_env%only_energy = .FALSE.
    p_env%new_preconditioner = .TRUE.
    p_env%os_valid = .FALSE.


    IF(linres_control%opt_method == lr_cg ) THEN

      p_env%ls_pos = 0.0_dp
      p_env%ls_energy = 0.0_dp
      p_env%ls_grad = 0.0_dp

      CPPrecondition(ASSOCIATED(p_env%g_psi1_old),cp_failure_level,routineP,error,failure)
      DO ispin = 1,nspins
        CALL cp_fm_get_info(p_env%g_psi1_old(ispin)%matrix,nrow_global=nrow,ncol_global=ncol)
        CPPrecondition(nrow==p_env%n_ao(ispin),cp_failure_level,routineP,error,failure)
        CPPrecondition(ncol==p_env%n_mo(ispin),cp_failure_level,routineP,error,failure)
        CALL cp_fm_set_all(p_env%g_psi1_old(ispin)%matrix,0.0_dp)
        CALL cp_fm_set_all(p_env%g_psi1(ispin)%matrix,0.0_dp)
        CALL cp_fm_set_all(p_env%d_psi1(ispin)%matrix,0.0_dp)
      END DO 

    ELSE IF (linres_control%opt_method == lr_sd) THEN

      p_env%ls_pos = 0.0_dp
      p_env%ls_energy = 0.0_dp
      p_env%ls_grad = 0.0_dp
      DO ispin = 1,nspins
        CALL cp_fm_get_info(p_env%g_psi1(ispin)%matrix,nrow_global=nrow,ncol_global=ncol)
        CPPrecondition(nrow==p_env%n_ao(ispin),cp_failure_level,routineP,error,failure)
        CPPrecondition(ncol==p_env%n_mo(ispin),cp_failure_level,routineP,error,failure)
        CALL cp_fm_set_all(p_env%g_psi1(ispin)%matrix,0.0_dp)
        CALL cp_fm_set_all(p_env%d_psi1(ispin)%matrix,0.0_dp)
      END DO

    ELSE IF (linres_control%opt_method == lr_diis) THEN
      DO ispin = 1,nspins
        CALL cp_fm_get_info( p_env%g_psi1(ispin)%matrix,nrow_global=nrow,ncol_global=ncol)
        CPPrecondition(nrow==p_env%n_ao(ispin),cp_failure_level,routineP,error,failure)
        CPPrecondition(ncol==p_env%n_mo(ispin),cp_failure_level,routineP,error,failure)
        CALL cp_fm_set_all(p_env%g_psi1(ispin)%matrix,0.0_dp)
        CALL cp_fm_set_all(p_env%d_psi1(ispin)%matrix,0.0_dp)
      END DO

    END IF

    IF(linres_control%preconditioner_type /= lr_precond_none) THEN

      CPPrecondition(ASSOCIATED(p_env%preconditioner),cp_failure_level,routineP,error,failure)
      Do ispin = 1,nspins
        CALL cp_fm_get_info(p_env%PS_psi0(ispin)%matrix,nrow_global=nrow,ncol_global=ncol)
        CPPrecondition(nrow==p_env%n_ao(ispin),cp_failure_level,routineP,error,failure)
        CPPrecondition(ncol==p_env%n_mo(ispin),cp_failure_level,routineP,error,failure)
      END DO
    END IF

  END SUBROUTINE check_p_env_init



! *****************************************************************************
  SUBROUTINE linres_scf_mini(p_env,qs_env,linres_control,b,Ab,nspins,iter,error)

    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(qs_environment_type), POINTER          :: qs_env
    TYPE(linres_control_type), POINTER          :: linres_control
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: b, Ab
    INTEGER, INTENT(IN)                         :: iter,nspins
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL                :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'linres_scf_mini', &
                                   routineP = moduleN//'/'//routineN

    INTEGER :: ip, ispin, istat, k, nao, nmo
    LOGICAL :: failure
    REAL(dp) :: polak_gamma, test_down, tmp
    TYPE(cp_fm_p_type), DIMENSION(2)            :: PS_psi0
    TYPE(cp_fm_type), POINTER                   :: mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                   :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                   :: matrix_ks, matrix_t, matrix_s



    failure = .FALSE.

    NULLIFY(PS_psi0(1)%matrix,PS_psi0(2)%matrix)
    NULLIFY(matrix_ks,matrix_s,matrix_t,mo_coeff,mos)
 
    CALL get_qs_env(qs_env=qs_env,&
         matrix_ks=matrix_ks,&
         matrix_s=matrix_s,&
         kinetic=matrix_t,&
         mos=mos,&
         error=error)

    IF(linres_control%preconditioner_type /= lr_precond_none) THEN

      IF(p_env%new_preconditioner) THEN
        p_env%os_valid = .FALSE.
        CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff)
        CALL make_preconditioner(p_env%preconditioner,linres_control%preconditioner_type,&
               matrix_ks(1)%matrix,matrix_s(1)%matrix,matrix_t(1)%matrix,&
               mo_coeff,linres_control%energy_gap)

        DO ispin = 1,nspins
           PS_psi0(ispin)%matrix => p_env%PS_psi0(ispin)%matrix
        END DO
        p_env%new_preconditioner = .FALSE.
      ENDIF
    ELSE
      DO ispin = 1,nspins
        PS_psi0(ispin)%matrix => p_env%S_psi0(ispin)%matrix
      END DO
    END IF


    IF(.NOT.p_env%os_valid) THEN
      DO ispin = 1,nspins
          IF(ASSOCIATED(p_env%preconditioner)) THEN
             CALL apply_preconditioner(p_env%preconditioner,p_env%S_psi0(ispin)%matrix,&
                   p_env%PS_psi0(ispin)%matrix)
          END IF
          CALL cp_fm_get_info(p_env%S_psi0(ispin)%matrix,nrow_global=nao,ncol_global=nmo)
          CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,p_env%S_psi0(ispin)%matrix,&
               PS_psi0(ispin)%matrix,0.0_dp,p_env%os_psi0(ispin)%matrix) 
          CALL cp_fm_cholesky_decompose(p_env%os_psi0(ispin)%matrix,error=error)
          CALL cp_fm_cholesky_invert(p_env%os_psi0(ispin)%matrix,error=error)
      END DO

      p_env%os_valid = .TRUE.
    END IF

    IF(p_env%ls_count==0) THEN
      CALL linres_derivative(Ab,PS_psi0,p_env,nspins,error=error)
    END IF
    
    SELECT CASE  (linres_control%opt_method)
      CASE (lr_sd)
         IF(p_env%ls_count==0) THEN
            p_env%gnorm = 0.0_dp
            DO ispin=1,nspins
               IF(linres_control%preconditioner_type /= lr_precond_none) THEN
                 CALL apply_preconditioner(p_env%preconditioner,&
                      p_env%g_psi1(ispin)%matrix,p_env%d_psi1(ispin)%matrix)
               ELSE
                 CALL cp_fm_to_fm(p_env%g_psi1(ispin)%matrix,&
                      p_env%d_psi1(ispin)%matrix,error=error) 
               END IF 
               CALL cp_fm_trace(p_env%g_psi1(ispin)%matrix,&
                    p_env%d_psi1(ispin)%matrix,tmp,error=error)
               p_env%gnorm = p_env%gnorm + tmp
               CALL cp_fm_scale_and_add(-1.0_dp,p_env%d_psi1(ispin)%matrix)
            ENDDO
            IF( p_env%gnorm .LT. 0.0) THEN
              WRITE(6,*) "WARNING Preconditioner not positive definite !"
            END IF
            
            k =  0 
            DO ispin = 1,nspins
              k = k + p_env%n_mo(ispin)
            END DO
            p_env%delta = SQRT(ABS(p_env%gnorm) / (p_env%n_ao(1)*k))
            p_env%gradient = -p_env%gnorm
         END IF
         CALL do_lr_line_search(p_env,b,linres_control,nspins,error=error)
      CASE (lr_cg)
         IF(p_env%ls_count==0) THEN
            p_env%gnorm = 0.0_dp
            p_env%gnorm_cross = 0.0_dp
            p_env%gradient = 0.0_dp
            DO ispin=1,nspins
               CALL cp_fm_trace(p_env%g_psi1(ispin)%matrix,&
                    p_env%d_psi1(ispin)%matrix,tmp,error=error)
               p_env%gradient = p_env%gradient + tmp
               CALL cp_fm_trace(p_env%g_psi1(ispin)%matrix,&
                    p_env%g_psi1_old(ispin)%matrix,tmp,error=error)
               p_env%gnorm_cross = p_env%gnorm_cross + tmp

               IF(linres_control%preconditioner_type /= lr_precond_none) THEN  
                  CALL apply_preconditioner(p_env%preconditioner,&
                       p_env%g_psi1(ispin)%matrix,p_env%g_psi1_old(ispin)%matrix)
                  CALL cp_fm_trace(p_env%g_psi1(ispin)%matrix,&
                       p_env%g_psi1_old(ispin)%matrix,tmp,error=error)
                  p_env%gnorm = p_env%gnorm + tmp 
                  IF(p_env%gnorm .LT. 0.0_dp) THEN
                      WRITE(6,*) "WARNING Preconditioner not positive definite !"
                  END IF
                  CALL cp_fm_to_fm(p_env%g_psi1_old(ispin)%matrix,&
                       p_env%g_psi1(ispin)%matrix,error=error)
               ELSE
                  CALL cp_fm_to_fm(p_env%g_psi1(ispin)%matrix,&
                       p_env%g_psi1_old(ispin)%matrix,error=error)
                  CALL cp_fm_trace(p_env%g_psi1(ispin)%matrix,&
                       p_env%g_psi1_old(ispin)%matrix,tmp,error=error)
                  p_env%gnorm = p_env%gnorm + tmp
               END IF
            END DO  ! ispin
            k =  0 
            DO ispin = 1,nspins
              k = k + p_env%n_mo(ispin)
            END DO
            p_env%delta = SQRT(ABS(p_env%gnorm) / (p_env%n_ao(1)*k))
            polak_gamma = (p_env%gnorm-p_env%gnorm_cross)/p_env%gnorm_old
            polak_gamma = MAX(polak_gamma,0.0_dp) ! reset to SD
            test_down = 0.0_dp
            DO ispin = 1,nspins
              CALL cp_fm_scale_and_add(polak_gamma,p_env%d_psi1(ispin)%matrix,&
                   -1.0_dp,p_env%g_psi1(ispin)%matrix,error=error)
              CALL cp_fm_trace(p_env%g_psi1(ispin)%matrix,&
                   p_env%d_psi1(ispin)%matrix,tmp,error=error)
              test_down = test_down + tmp
            END DO ! ispin
            IF (test_down.ge.0.0_dp) THEN ! reset to SD
              polak_gamma = 0.0_dp
              DO ispin = 1,nspins
                CALL cp_fm_scale_and_add(polak_gamma,p_env%d_psi1(ispin)%matrix,&
                     -1.0_dp,p_env%g_psi1(ispin)%matrix,error=error)
              END DO ! ispin
            END IF
            p_env%gradient = polak_gamma*p_env%gradient - p_env%gnorm
            p_env%gnorm_old = p_env%gnorm
          END IF
          CALL do_lr_line_search(p_env,b,linres_control,nspins,error=error) 
      CASE (lr_diis)
           CALL stop_program("linres_scf_mini","LR_DIIS NYI")
      CASE DEFAULT
           CALL stop_program("linres_scf_mini","OPT_METHOD NYI")
    END SELECT

  END SUBROUTINE linres_scf_mini


! *****************************************************************************
  SUBROUTINE linres_derivative(Ab,PS_psi0,p_env,nspins,error)

    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: Ab
    TYPE(cp_fm_p_type), DIMENSION(2)            :: PS_psi0
    TYPE(qs_p_env_type), POINTER                :: p_env
    INTEGER :: nspins
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    CHARACTER(len=*), PARAMETER :: routine_name='linres_derivative',&
      routineP = moduleN//':'//routine_name

    INTEGER :: fcol, ip, ispin, lcol, n_ao, n_mo, n_range
    LOGICAL :: failure

    failure = .FALSE.
    DO ispin = 1,nspins
       n_mo = p_env%n_mo(ispin)
       n_ao = p_env%n_ao(ispin)
       CALL cp_fm_to_fm(Ab(ispin)%matrix,p_env%g_psi1(ispin)%matrix,error=error) 
       CALL cp_fm_gemm('T','N',n_mo,n_mo,n_ao,1.0_dp,PS_psi0(ispin)%matrix,&
            Ab(ispin)%matrix,0.0_dp,p_env%buf1(ispin)%matrix) 
       CALL cp_fm_symm('L','U',n_mo,n_mo,1.0_dp,p_env%os_psi0(ispin)%matrix,&
            p_env%buf1(ispin)%matrix,0.0_dp,p_env%buf2(ispin)%matrix)
       CALL cp_fm_gemm('N','N',n_ao,n_mo,n_mo,-1.0_dp,p_env%S_psi0(ispin)%matrix,&
            p_env%buf2(ispin)%matrix,1.0_dp,p_env%g_psi1(ispin)%matrix)
    END DO ! ispin

  END SUBROUTINE linres_derivative

! *****************************************************************************
  SUBROUTINE do_lr_line_search(p_env,vectors,linres_control,nspins,error)

    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: vectors
    TYPE(linres_control_type), POINTER          :: linres_control
    INTEGER :: nspins
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                  :: error

    INTEGER  :: ispin
    REAL(dp) :: ds

    SELECT CASE(linres_control%ls_method)

      CASE(lr_ls_gold)
        ! to do
      CASE(lr_ls_3pnt)
         CALL do_lr_ls_3pnt(p_env,vectors,linres_control%ds_min,nspins,error=error)
      CASE(lr_ls_2pnt)
         CALL do_lr_ls_2pnt(p_env,vectors,linres_control%ds_min,nspins,error=error)
      CASE(lr_ls_none)
         ds=linres_control%ds_min
         DO ispin = 1,nspins
           CALL cp_fm_scale_and_add(1.0_dp,vectors(ispin)%matrix,ds,&
                p_env%d_psi1(ispin)%matrix,error=error)
         END DO 
      CASE DEFAULT
          CALL stop_program("do_lr_line_search","NYI")
    END SELECT

  END SUBROUTINE do_lr_line_search

! *****************************************************************************
  SUBROUTINE do_lr_ls_2pnt(p_env,vectors,ds_min,nspins,error)

    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: vectors
    REAL(dp)                                    :: ds_min
    INTEGER                                     :: nspins
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL                :: error

    INTEGER :: count, ispin
    REAL(dp) :: a, b, c, ds, pos, val, x0, x1
    REAL(DP), DIMENSION(:), POINTER             :: ls_pos 

    
    NULLIFY(ls_pos)
    p_env%ls_count = p_env%ls_count + 1
    count = p_env%ls_count 
    p_env%only_energy = .TRUE.
    p_env%ls_energy(count) = p_env%etotal
    ls_pos=>p_env%ls_pos
    SELECT CASE (count)

      CASE(1)
         ls_pos(count)   = 0.0_dp
         p_env%ls_grad(count)  = p_env%gnorm
         p_env%ls_pos(count+1) = ds_min
      CASE(2)

         x0 = 0.0_dp
         x1 = ls_pos(2)
         c = p_env%ls_energy(1)
         b = p_env%ls_grad(1)
         a = (p_env%ls_energy(2)-b*x1)/(x1*x1)
         IF(a .le. EPSILON(1.0_dp)) a = 1.E-15_dp
         pos = -b/(2.0_dp*a)
         val = a*pos*pos+b*pos+c 
         IF(val .LT. c .AND. val .LE. p_env%ls_energy(2)) THEN
           ls_pos(count+1) = MAX(MAXVAL(ls_pos(1:2))*0.01_dp,&
              MIN(pos,MAXVAL(ls_pos(1:2))*4.0_dp)) 
         ELSE
           ls_pos(count+1) = MAXVAL(ls_pos(1:2))*2.0_dp
         END IF
         p_env%only_energy = .FALSE.
         p_env%ls_count = 0 

      CASE DEFAULT
          CALL stop_program("lr_ls_2pnt","too many steps counted")
    END SELECT

    ds = ls_pos(count+1) - ls_pos(count)
    ds_min = ls_pos(count+1)
    DO ispin = 1,nspins
       CALL cp_fm_scale_and_add(1.0_dp,vectors(ispin)%matrix,ds,&
            p_env%d_psi1(ispin)%matrix,error=error)
    END DO 

  END SUBROUTINE do_lr_ls_2pnt


! *****************************************************************************
  SUBROUTINE do_lr_ls_3pnt(p_env,vectors,ds_min,nspins,error)

    TYPE(qs_p_env_type), POINTER                :: p_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                   :: vectors
    REAL(dp)                                    :: ds_min
    INTEGER                                     :: nspins
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL                :: error

    INTEGER :: count, ispin
    REAL(dp) :: a, b, c, den, ds, f1, f2, f3, nom, pos, val, x1, x2, x3
    REAL(DP), DIMENSION(:), POINTER             :: ls_pos 


    nullify(ls_pos)  
    p_env%ls_count = p_env%ls_count + 1
    count = p_env%ls_count 
    p_env%only_energy = .TRUE.
    p_env%ls_energy(count) = p_env%etotal

    SELECT CASE (count)
      CASE(1) 
        ls_pos(count) = 0.0_dp
        ls_pos(count+1) = ds_min*0.8_dp
      CASE(2) 
        IF(p_env%ls_energy(count) .GT. p_env%ls_energy(count-1)) THEN
          ls_pos(count+1) = ds_min*0.5_dp
        ELSE
          ls_pos(count+1) = ds_min*1.4_dp
        END IF
      CASE(3) 
         x1 = ls_pos(1)
         f1 = p_env%ls_energy(1)
         x2 = ls_pos(2)
         f2 = p_env%ls_energy(2)
         x3 = ls_pos(3)
         f3 = p_env%ls_energy(3)
         nom =  (x2-x1)**2*(f2-f3)-(x2-x3)**2*(f2-f1) 
         den =  (x2-x1)*(f2-f3)-(x2-x1)*(f2-f1)
         IF(ABS(den).LE.EPSILON(1.0_dp)*MAX(ABS(f2-f1),ABS(f2-f3))) THEN
            pos = x2
         ELSE
            CALL stop_program("b value not defined here...",&
                              "qs_linres_methods.F")
            pos = b-0.5_dp*nom/den
         END IF
         val = (pos-x1)*(pos-x2)*f3/((x3-x1)*(x3-x2))+ &
               (pos-x2)*(pos-x3)*f1/((x1-x2)*(x1-x3))+ &
               (pos-x3)*(pos-x1)*f2/((x2-x3)*(x2-x1))
         IF(val .LT. f1 .AND. val .LE. f2 .AND. val .LE. f3) THEN
             ls_pos(count+1) = MAX(MAXVAL(ls_pos(1:3))*0.01_dp,&
                   MIN(pos,MAXVAL(ls_pos(1:3))*4.0_dp))
         ELSE
             ls_pos(count+1) = MAXVAL(ls_pos(1:3))*2.0_dp
         END IF

         p_env%only_energy = .FALSE.
         p_env%ls_count = 0 
      CASE DEFAULT
         CALL stop_program("lr_ls_2pnt","too many steps counted")
    END SELECT

    ds = ls_pos(count+1) - ls_pos(count)
    ds_min = ls_pos(count+1)
    DO ispin = 1,nspins
       CALL cp_fm_scale_and_add(1.0_dp,vectors(ispin)%matrix,ds,&
            p_env%d_psi1(ispin)%matrix,error=error)
    END DO 

  END SUBROUTINE do_lr_ls_3pnt


! *****************************************************************************
END MODULE qs_linres_methods
