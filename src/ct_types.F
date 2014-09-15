!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types for all cayley transformation methods
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE ct_types
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_release,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: cg_polak_ribiere,&
                                             tensor_orthogonal,&
                                             tensor_up_down
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ct_types'

  ! Public types
  PUBLIC :: ct_step_env_type
  
  ! Public subroutines
  PUBLIC :: ct_step_env_init,ct_step_env_set,ct_step_env_get,ct_step_env_clean
  
  TYPE ct_step_env_type

    ! this type contains options for cayley transformation routines 

    ! use orbitals or projectors?
    LOGICAL :: use_occ_orbs, use_virt_orbs
    LOGICAL :: occ_orbs_orthogonal, virt_orbs_orthogonal
    ! tensor properties of matrix indeces:
    ! tensor_up_down, tensor_orthogonal
    INTEGER :: tensor_type
    ! neglect the quadratic term in riccati equations?
    LOGICAL :: neglect_quadratic_term
    ! what kind of output do we produce?
    LOGICAL :: update_p, update_q, calculate_energy_corr
    ! variety of conjugate gradient
    INTEGER :: conjugator

    ! type of preconditioner
    LOGICAL :: pp_preconditioner_full,&
               qq_preconditioner_full

    REAL(KIND=dp)         :: eps_convergence
    REAL(KIND=dp)         :: eps_filter
    INTEGER               :: max_iter
    !INTEGER               :: nspins
    LOGICAL               :: converged
    INTEGER               :: order_lanczos
    REAL(KIND=dp)         :: eps_lancsoz
    INTEGER               :: max_iter_lanczos

    REAL(KIND=dp)         :: energy_correction

!SPIN!!!    ! metric matrices for covariant to contravariant transformations 
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: p_index_up=>NULL()
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: p_index_down=>NULL()
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: q_index_up=>NULL()
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: q_index_down=>NULL()
!SPIN!!!
!SPIN!!!    ! kohn-sham, covariant-covariant representation
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: matrix_ks=>NULL()
!SPIN!!!    ! density, contravariant-contravariant representation
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: matrix_p=>NULL()
!SPIN!!!    ! occ orbitals, contravariant-covariant representation
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: matrix_t=>NULL()
!SPIN!!!    ! virt orbitals, contravariant-covariant representation
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: matrix_v=>NULL()
!SPIN!!!
!SPIN!!!    ! to avoid bulding Occ-by-N and Virt-vy-N matrices inside
!SPIN!!!    ! the ct routines get them from the external code
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: matrix_qp_template=>NULL()
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER  :: matrix_pq_template=>NULL()
!SPIN!!!
!SPIN!!!    ! single excitation amplitudes
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE  :: matrix_x
!SPIN!!!    ! residuals
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE  :: matrix_res

    ! metric matrices for covariant to contravariant transformations 
    TYPE(cp_dbcsr_type), POINTER  :: p_index_up=>NULL()
    TYPE(cp_dbcsr_type), POINTER  :: p_index_down=>NULL()
    TYPE(cp_dbcsr_type), POINTER  :: q_index_up=>NULL()
    TYPE(cp_dbcsr_type), POINTER  :: q_index_down=>NULL()

    ! kohn-sham, covariant-covariant representation
    TYPE(cp_dbcsr_type), POINTER  :: matrix_ks=>NULL()
    ! density, contravariant-contravariant representation
    TYPE(cp_dbcsr_type), POINTER  :: matrix_p=>NULL()
    ! occ orbitals, contravariant-covariant representation
    TYPE(cp_dbcsr_type), POINTER  :: matrix_t=>NULL()
    ! virt orbitals, contravariant-covariant representation
    TYPE(cp_dbcsr_type), POINTER  :: matrix_v=>NULL()

    ! to avoid bulding Occ-by-N and Virt-vy-N matrices inside
    ! the ct routines get them from the external code
    TYPE(cp_dbcsr_type), POINTER  :: matrix_qp_template=>NULL()
    TYPE(cp_dbcsr_type), POINTER  :: matrix_pq_template=>NULL()

    ! guess for single excitation amplitudes
    ! it is used exclusively as a guess, not modified
    ! it should be given in the up_down representation
    TYPE(cp_dbcsr_type), POINTER  :: matrix_x_guess=>NULL()

    ! single excitation amplitudes
    TYPE(cp_dbcsr_type)           :: matrix_x
    ! residuals
    TYPE(cp_dbcsr_type)           :: matrix_res

    TYPE(cp_para_env_type),  POINTER  :: para_env=>NULL()
    TYPE(cp_blacs_env_type), POINTER  :: blacs_env=>NULL()

  END TYPE

CONTAINS

  SUBROUTINE ct_step_env_init(env,error)
  
    TYPE(ct_step_env_type)                   :: env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ct_step_env_init', &
      routineP = moduleN//':'//routineN

    env%use_occ_orbs = .TRUE. 
    env%use_virt_orbs = .FALSE. 
    env%occ_orbs_orthogonal = .FALSE. 
    env%virt_orbs_orthogonal = .FALSE. 
    env%tensor_type = tensor_orthogonal 
    env%neglect_quadratic_term = .FALSE. 
    env%calculate_energy_corr = .TRUE. 
    env%update_p = .FALSE. 
    env%update_q = .FALSE. 
    env%pp_preconditioner_full=.TRUE.   
    env%qq_preconditioner_full=.FALSE.   
 
    env%eps_convergence = 1.0E-8_dp 
    env%eps_filter = 1.0E-8_dp
    env%max_iter = 400
    env%order_lanczos = 3
    env%eps_lancsoz = 1.0E-4_dp
    env%max_iter_lanczos = 40
    !env%nspins = -1
    env%converged = .FALSE.
    env%conjugator = cg_polak_ribiere

    NULLIFY(env%p_index_up)
    NULLIFY(env%p_index_down)
    NULLIFY(env%q_index_up)
    NULLIFY(env%q_index_down)
    
    NULLIFY(env%matrix_ks)
    NULLIFY(env%matrix_p)
    NULLIFY(env%matrix_t)
    NULLIFY(env%matrix_v)
    NULLIFY(env%matrix_x_guess)
    NULLIFY(env%matrix_qp_template)
    NULLIFY(env%matrix_pq_template)

    !RZK-warning read_parameters_from_input
    
  END SUBROUTINE ct_step_env_init

  SUBROUTINE ct_step_env_get(env,use_occ_orbs,use_virt_orbs,tensor_type,&
    occ_orbs_orthogonal,virt_orbs_orthogonal,neglect_quadratic_term,&
    update_p,update_q,eps_convergence,eps_filter,max_iter,&
    p_index_up,p_index_down,q_index_up,q_index_down,matrix_ks,matrix_p,&
    matrix_qp_template,matrix_pq_template,&
    matrix_t,matrix_v,matrix_x,copy_matrix_x,energy_correction,calculate_energy_corr,&
    converged,qq_preconditioner_full,pp_preconditioner_full,&
    error)

    TYPE(ct_step_env_type)                   :: env
    LOGICAL, OPTIONAL                        :: use_occ_orbs, use_virt_orbs
    INTEGER, OPTIONAL                        :: tensor_type
    LOGICAL, OPTIONAL :: occ_orbs_orthogonal, virt_orbs_orthogonal, &
      neglect_quadratic_term, update_p, update_q
    REAL(KIND=dp), OPTIONAL                  :: eps_convergence, eps_filter
    INTEGER, OPTIONAL                        :: max_iter
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER :: p_index_up, p_index_down, &
      q_index_up, q_index_down, matrix_ks, matrix_p, matrix_qp_template, &
      matrix_pq_template, matrix_t, matrix_v, matrix_x
    TYPE(cp_dbcsr_type), OPTIONAL            :: copy_matrix_x
    REAL(KIND=dp), OPTIONAL                  :: energy_correction
    LOGICAL, OPTIONAL                        :: calculate_energy_corr, &
                                                converged, &
                                                qq_preconditioner_full, &
                                                pp_preconditioner_full
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ct_step_env_get', &
      routineP = moduleN//':'//routineN

!INTEGER      , OPTIONAL                     :: nspins
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: p_index_up
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: p_index_down
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: q_index_up
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: q_index_down
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: matrix_ks
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: matrix_p
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: matrix_t
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: matrix_v
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: matrix_qp_template
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: matrix_pq_template
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), POINTER, OPTIONAL  :: matrix_x
!SPIN!!!    
!SPIN!!!    TYPE(cp_dbcsr_type), DIMENSION(:), OPTIONAL           :: copy_matrix_x
!INTEGER                                               :: ispin

    IF (PRESENT(use_occ_orbs)) use_occ_orbs = env%use_occ_orbs 
    IF (PRESENT(use_virt_orbs)) use_virt_orbs = env%use_virt_orbs 
    IF (PRESENT(occ_orbs_orthogonal)) occ_orbs_orthogonal = &
       env%occ_orbs_orthogonal
    IF (PRESENT(virt_orbs_orthogonal)) virt_orbs_orthogonal = &
       env%virt_orbs_orthogonal
    IF (PRESENT(tensor_type)) tensor_type = env%tensor_type 
    IF (PRESENT(neglect_quadratic_term)) neglect_quadratic_term = &
       env%neglect_quadratic_term
    IF (PRESENT(calculate_energy_corr)) calculate_energy_corr = &
       env%calculate_energy_corr 
    IF (PRESENT(update_p)) update_p = env%update_p 
    IF (PRESENT(update_q)) update_q = env%update_q
    IF (PRESENT(pp_preconditioner_full)) pp_preconditioner_full = &
       env%pp_preconditioner_full 
    IF (PRESENT(qq_preconditioner_full)) qq_preconditioner_full = &
       env%qq_preconditioner_full 
    IF (PRESENT(eps_convergence)) eps_convergence = env%eps_convergence
    IF (PRESENT(eps_filter)) eps_filter = env%eps_filter
    IF (PRESENT(max_iter)) max_iter = env%max_iter
    !IF (PRESENT(nspins)) nspins = env%nspins 
    IF (PRESENT(matrix_ks)) matrix_ks => env%matrix_ks 
    IF (PRESENT(matrix_p)) matrix_p => env%matrix_p 
    IF (PRESENT(matrix_t)) matrix_t => env%matrix_t 
    IF (PRESENT(matrix_v)) matrix_v => env%matrix_v 
    IF (PRESENT(matrix_qp_template)) matrix_qp_template => &
       env%matrix_qp_template
    IF (PRESENT(matrix_pq_template)) matrix_pq_template => &
       env%matrix_pq_template
    IF (PRESENT(p_index_up)) p_index_up => env%p_index_up 
    IF (PRESENT(q_index_up)) q_index_up => env%q_index_up 
    IF (PRESENT(p_index_down)) p_index_down => env%p_index_down 
    IF (PRESENT(q_index_down)) q_index_down => env%q_index_down
    IF (PRESENT(copy_matrix_x)) THEN
       !DO ispin=1,env%nspins
          !CALL cp_dbcsr_copy(copy_matrix_x(ispin),env%matrix_x(ispin),error=error)
          CALL cp_dbcsr_copy(copy_matrix_x,env%matrix_x,error=error)
       !ENDDO
    ENDIF
    !IF (PRESENT(matrix_x)) matrix_x => env%matrix_x
    IF (PRESENT(energy_correction)) energy_correction = env%energy_correction
    IF (PRESENT(converged)) converged = env%converged

  END SUBROUTINE ct_step_env_get

  SUBROUTINE ct_step_env_set(env,para_env,blacs_env,use_occ_orbs,&
    use_virt_orbs,tensor_type,&
    occ_orbs_orthogonal,virt_orbs_orthogonal,neglect_quadratic_term,&
    update_p,update_q,eps_convergence,eps_filter,max_iter,&
    p_index_up,p_index_down,q_index_up,q_index_down,matrix_ks,matrix_p,&
    matrix_qp_template,matrix_pq_template,&
    matrix_t,matrix_v,matrix_x_guess,calculate_energy_corr,conjugator,&
    qq_preconditioner_full,pp_preconditioner_full,&
    error)

    TYPE(ct_step_env_type)                   :: env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    LOGICAL, OPTIONAL                        :: use_occ_orbs, use_virt_orbs
    INTEGER, OPTIONAL                        :: tensor_type
    LOGICAL, OPTIONAL :: occ_orbs_orthogonal, virt_orbs_orthogonal, &
      neglect_quadratic_term, update_p, update_q
    REAL(KIND=dp), OPTIONAL                  :: eps_convergence, eps_filter
    INTEGER, OPTIONAL                        :: max_iter
    TYPE(cp_dbcsr_type), OPTIONAL, TARGET :: p_index_up, p_index_down, &
      q_index_up, q_index_down, matrix_ks, matrix_p, matrix_qp_template, &
      matrix_pq_template, matrix_t, matrix_v, matrix_x_guess
    LOGICAL, OPTIONAL                        :: calculate_energy_corr
    INTEGER, OPTIONAL                        :: conjugator
    LOGICAL, OPTIONAL                        :: qq_preconditioner_full, &
                                                pp_preconditioner_full
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ct_step_env_set', &
      routineP = moduleN//':'//routineN

!INTEGER      , OPTIONAL                     :: nspins
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: p_index_up
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: p_index_down
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: q_index_up
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: q_index_down
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: matrix_ks
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: matrix_p
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: matrix_t
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: matrix_v
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: matrix_qp_template
!SPIN!!!    TYPE(cp_dbcsr_type), TARGET, DIMENSION(:), OPTIONAL  :: matrix_pq_template
! set para_env and blacs_env which are needed to operate with full matrices
! it would be nice to have everything with cp_dbcsr matrices, well maybe later

    env%para_env => para_env 
    env%blacs_env => blacs_env 

    IF (PRESENT(use_occ_orbs)) env%use_occ_orbs = use_occ_orbs 
    IF (PRESENT(use_virt_orbs)) env%use_virt_orbs = use_virt_orbs 
    IF (PRESENT(occ_orbs_orthogonal)) env%occ_orbs_orthogonal = &
       occ_orbs_orthogonal 
    IF (PRESENT(virt_orbs_orthogonal)) env%virt_orbs_orthogonal = &
       virt_orbs_orthogonal 
    IF (PRESENT(tensor_type)) env%tensor_type = tensor_type 
    IF (PRESENT(neglect_quadratic_term)) env%neglect_quadratic_term = &
       neglect_quadratic_term
    IF (PRESENT(calculate_energy_corr)) env%calculate_energy_corr = &
       calculate_energy_corr
    IF (PRESENT(update_p)) env%update_p = update_p 
    IF (PRESENT(update_q)) env%update_q = update_q
    IF (PRESENT(pp_preconditioner_full)) env%pp_preconditioner_full = &
       pp_preconditioner_full 
    IF (PRESENT(qq_preconditioner_full)) env%qq_preconditioner_full = &
       qq_preconditioner_full 
    IF (PRESENT(eps_convergence)) env%eps_convergence = eps_convergence
    IF (PRESENT(eps_filter)) env%eps_filter = eps_filter
    IF (PRESENT(max_iter)) env%max_iter = max_iter
    !IF (PRESENT(nspins)) env%nspins = nspins 
    IF (PRESENT(conjugator)) env%conjugator = conjugator 
    IF (PRESENT(matrix_ks)) env%matrix_ks => matrix_ks 
    IF (PRESENT(matrix_p)) env%matrix_p => matrix_p 
    IF (PRESENT(matrix_t)) env%matrix_t => matrix_t 
    IF (PRESENT(matrix_v)) env%matrix_v => matrix_v 
    IF (PRESENT(matrix_x_guess)) env%matrix_x_guess => matrix_x_guess
    IF (PRESENT(matrix_qp_template)) env%matrix_qp_template => &
       matrix_qp_template 
    IF (PRESENT(matrix_pq_template)) env%matrix_pq_template => &
       matrix_pq_template 
    IF (PRESENT(p_index_up)) env%p_index_up => p_index_up 
    IF (PRESENT(q_index_up)) env%q_index_up => q_index_up 
    IF (PRESENT(p_index_down)) env%p_index_down => p_index_down 
    IF (PRESENT(q_index_down)) env%q_index_down => q_index_down 

  END SUBROUTINE ct_step_env_set

  SUBROUTINE ct_step_env_clean(env,error) 

    TYPE(ct_step_env_type)                   :: env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ct_step_env_clean', &
      routineP = moduleN//':'//routineN

!INTEGER                                     :: ispin

    NULLIFY(env%para_env)
    NULLIFY(env%blacs_env)

    !DO ispin=1,env%nspins
       CALL cp_dbcsr_release(env%matrix_x,error=error)
       CALL cp_dbcsr_release(env%matrix_res,error=error)
       !CALL cp_dbcsr_release(env%matrix_x(ispin),error=error)
       !CALL cp_dbcsr_release(env%matrix_res(ispin),error=error)
    !ENDDO
    !DEALLOCATE(env%matrix_x,env%matrix_res)
    
    NULLIFY(env%p_index_up)
    NULLIFY(env%p_index_down)
    NULLIFY(env%q_index_up)
    NULLIFY(env%q_index_down)
    
    NULLIFY(env%matrix_ks)
    NULLIFY(env%matrix_p)
    NULLIFY(env%matrix_t)
    NULLIFY(env%matrix_v)
    NULLIFY(env%matrix_x_guess)
    NULLIFY(env%matrix_qp_template)
    NULLIFY(env%matrix_pq_template)

  END SUBROUTINE ct_step_env_clean

END MODULE ct_types

