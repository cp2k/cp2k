!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****m* cp2k/qs_tddfpt_utils *
!!
!!   NAME
!!     qs_tddfpt_utils
!!
!!   FUNCTION
!!
!!   SYNOPSIS
!!     
!!   MODIFICATION HISTORY
!!     09-JUL-2002, TCH, development started
!!
!!   SEE ALSO
!!     
!!****

MODULE qs_tddfpt_utils

  USE cp_control_types,                ONLY: gapw_control_type,&
                                             tddfpt_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale, &
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_pool_types,                ONLY: fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type, &
                                             cp_fm_to_fm
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE hartree_local_methods,           ONLY: init_coulomb_local
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: qs_environment_type, &
                                             get_qs_env
  USE qs_kpp1_env_methods,             ONLY: kpp1_create
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_rho_atom_methods,             ONLY: init_rho_atom
  USE qs_rho0_methods,                 ONLY: init_rho0
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_allocate,&
                                             tddfpt_env_deallocate,&
                                             tddfpt_env_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: moduleN="qs_tddfpt_utils"
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

  TYPE simple_solution_sorter
     INTEGER                               :: orbit
     INTEGER                               :: lumo
     REAL(KIND=DP)                        :: value
     TYPE(simple_solution_sorter), POINTER :: next
  END TYPE simple_solution_sorter

  PRIVATE
  
  ! METHODS
  PUBLIC :: tddfpt_cleanup, &
            tddfpt_init, &
            tddfpt_write_banner, &
            tddfpt_env_init, &
            co_initial_guess, &
            calc_norm, &
            find_contributions,&
            normalize, &
            reorthogonalize, &
            dot_prod

CONTAINS

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_utils/tddfpt_init *
!!
!!   NAME
!!     tddfpt_init
!!
!!   SYNOPSIS
!!     Subroutine tddfpt_init(p_env, t_env, qs_env, n_ev, glob_env)
!!       Implicit None
!!       Type(qs_p_env_type), Pointer:: p_env
!!       Type(tddfpt_env_type), Intent (OUT):: t_env
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (IN):: glob_env
!!     End Subroutine tddfpt_init
!!
!!   FUNCTION
!!     Initialize some necessary structures for a tddfpt calculation.
!!
!!   ARGUMENTS
!!     - p_env    : perturbation environment to be initialized
!!     - t_env    : tddfpt environment to be initialized
!!     - qs_env   : Quickstep environment with the results of a
!!                  ground state calcualtion
!!     - glob_env : global environment
!!
!!   NOTES
!!     -
!!
!!****
  SUBROUTINE tddfpt_init(p_env, t_env, qs_env)


    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(tddfpt_env_type), INTENT(out)       :: t_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = "tddfpt_init", &
                                   routineP = moduleN//"/"//routineN

! --------------------------------------------------------------------

    !------------------!
    ! create the p_env !
    !------------------!
    CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE.)
    CALL p_env_psi0_changed(p_env, qs_env) ! update the m_epsilon matrix

    !------------------!
    ! create the t_env !
    !------------------!
    CALL tddfpt_env_allocate(t_env, p_env, qs_env)
    CALL tddfpt_env_init(t_env, p_env, qs_env)

  END SUBROUTINE tddfpt_init

! *****************************************************************************

!!****f* cp2k/qs_tddfpt_types/tddfpt_env_init *
!!
!!   FUNCTION
!!     Initialize t_env with meaningfull values.
!!
!!***
  SUBROUTINE tddfpt_env_init(t_env, p_env, qs_env)

  
    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_init', &
      routineP = moduleN//'/'//routineN

    INTEGER                                  :: n_spins, spin

! -------------------------------------------------------------------

    n_spins = qs_env%dft_control%nspins
    
    IF (qs_env%dft_control%tddfpt_control%invert_S) THEN
       DO spin=1, n_spins
          CALL copy_sm_to_fm(qs_env%matrix_s(spin)%matrix, t_env%invS(spin)%matrix)
          CALL cp_fm_cholesky_decompose(t_env%invS(spin)%matrix)
          CALL cp_fm_cholesky_invert(t_env%invS(spin)%matrix)
       END DO
    END IF

  END SUBROUTINE tddfpt_env_init

! *****************************************************************************

  SUBROUTINE tddfpt_write_banner(globenv)

    
    TYPE(global_environment_type), &
      INTENT(in)                             :: globenv

    IF (.NOT.globenv%ionode) RETURN
    IF (globenv%print%program_banner) THEN
       
       WRITE (unit=globenv%scr, fmt="(2X,A)")&
       "",&
       "-----------------------------------------------------------------------------",&
       "-                    Time-dependent calculation started                     -",&
       "-----------------------------------------------------------------------------",&
       ""    

    END IF

  END SUBROUTINE tddfpt_write_banner

! *****************************************************************************

  SUBROUTINE tddfpt_cleanup(t_env, p_env)


    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env

    CHARACTER(LEN=*), PARAMETER :: routineN = "cleanup", &
                                   routineP = moduleN//'/'//routineN

! -------------------------------------------------------------------

    CALL tddfpt_env_deallocate(t_env)
    CALL p_env_release(p_env)
    
  END SUBROUTINE tddfpt_cleanup

! *****************************************************************************

  SUBROUTINE co_initial_guess(matrices, energies, n_v, qs_env, p_env)


    TYPE(cp_fm_p_type), DIMENSION(:,:),pointer :: matrices
    INTEGER, INTENT(IN)                      :: n_v ! number of vectors to initialize
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    REAL(kind=DP), DIMENSION(:), &
      INTENT(OUT)                            :: energies

    INTEGER                                  :: i, j, n_cols, n_lumos, &
                                                n_orbits, n_rows, n_spins, &
                                                spin, stat, oo, vo
    REAL(KIND=DP)                            :: scale_factor, evd
    REAL(KIND=DP), ALLOCATABLE, &
      DIMENSION(:, :)                        :: lumos, guess
    REAL(KIND=DP), DIMENSION(:), POINTER     :: orbital_eigenvalues
    TYPE(simple_solution_sorter), POINTER    :: sorter_iterator, &
                                                sorter_pointer, sorter_start
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

    CHARACTER(len=*), PARAMETER :: routineN = "co_initial_guess", &
                                   routineP = moduleN//"/"//routineN

    ! *************************

    tddfpt_control => qs_env%dft_control%tddfpt_control
    n_spins = qs_env%dft_control%nspins
    energies = 0.0_dp

    IF (.NOT.ASSOCIATED(tddfpt_control%lumos)) THEN
       CALL stop_program(routineP,"LUMOS missing")
    END IF

    DO spin=1, n_spins

       n_cols = matrices(1,spin)%matrix%matrix_struct%ncol_global
       n_rows = matrices(1,spin)%matrix%matrix_struct%nrow_global
    
       DO i=1, n_v
          CALL cp_fm_set_all(matrices(i,spin)%matrix, 0.0_dp)
       END DO

       CALL get_mo_set(qs_env%mos(spin)%mo_set, eigenvalues=orbital_eigenvalues)
    
       n_lumos = tddfpt_control%lumos(spin)%matrix%matrix_struct%ncol_global

       n_orbits = SIZE(orbital_eigenvalues)

       !-----------------------------------------!
       ! create a SORTED list of initial guesses !
       !-----------------------------------------!
       ! first element
       evd = tddfpt_control%lumos_eigenvalues(1,spin)-orbital_eigenvalues(n_orbits)
       ALLOCATE(sorter_start)
       sorter_start%orbit=n_orbits
       sorter_start%lumo=1
       sorter_start%value=evd
       NULLIFY(sorter_start%next)
       ! rest of the elements
       DO oo=n_orbits, 1, -1
          DO vo=1, n_lumos

             IF (oo==n_orbits .AND. vo==1) CYCLE ! already in list

             evd = tddfpt_control%lumos_eigenvalues(vo,spin)-orbital_eigenvalues(oo)

             sorter_iterator => sorter_start
             NULLIFY(sorter_pointer)
             DO WHILE (ASSOCIATED(sorter_iterator%next))
                IF (sorter_iterator%next%value > evd) THEN
                   sorter_pointer => sorter_iterator%next
                   EXIT
                END IF
                sorter_iterator => sorter_iterator%next
             END DO

             ALLOCATE(sorter_iterator%next)
             sorter_iterator%next%orbit = oo
             sorter_iterator%next%lumo  = vo
             sorter_iterator%next%value = evd
             sorter_iterator%next%next  => sorter_pointer

          END DO
       END DO

       ALLOCATE(lumos(n_rows, n_lumos), guess(n_rows, n_orbits), stat=stat)
       CALL cp_fm_get_submatrix(tddfpt_control%lumos(spin)%matrix, lumos, &
                                start_col=1, n_cols=n_lumos)

       !-------------------!
       ! fill the matrices !
       !-------------------!
       sorter_iterator=>sorter_start
       DO i=1, MIN(n_v,n_orbits*n_lumos)
          CALL dcopy(n_rows*n_orbits, 0.0_dp, 0, guess, 1)
          CALL dcopy(n_rows, lumos(:,sorter_iterator%lumo), 1, &
                             guess(:,sorter_iterator%orbit),1)
          CALL cp_fm_set_submatrix(matrices(i,spin)%matrix, &
                                   guess(:,1:n_orbits))
          energies(i) = sorter_iterator%value
          sorter_iterator=>sorter_iterator%next
       END DO
       IF (n_v > n_orbits*n_lumos) THEN
          DO i=n_orbits*n_lumos+1, n_v
             CALL cp_fm_init_random(matrices(i,spin)%matrix, n_orbits)
             energies(i)=1.0E38_dp
          END DO
       END IF


       !--------------!
       ! some cleanup !
       !--------------!
       DEALLOCATE(lumos, guess)
       sorter_iterator=>sorter_start
       DO WHILE(ASSOCIATED(sorter_iterator))
          sorter_pointer=>sorter_iterator
          sorter_iterator=>sorter_iterator%next
          DEALLOCATE(sorter_pointer)
       END DO
          

    END DO

  END SUBROUTINE co_initial_guess

! *****************************************************************************
  
  SUBROUTINE find_contributions(qs_env, t_env, error)
    
    TYPE(tddfpt_env_type), INTENT(IN)  :: t_env
    TYPE(qs_environment_type), POINTER :: qs_env
    TYPE(cp_error_type), &
         INTENT(INOUT), OPTIONAL       :: error

    CHARACTER(LEN=*), PARAMETER        :: routineN = 'find_contributions', &
                                          routineP = moduleN//'/'//routineN

    INTEGER                            :: i, j, n_ev, spin, n_spins, stat
    INTEGER                            :: occ, virt
    INTEGER, DIMENSION(2)              :: nrows, nhomos, nlumos
    LOGICAL                            :: failure, ionode 
    REAL(KIND=dp)                      :: contribution, summed_contributions
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                  :: homo_coeff_col, lumo_coeff_col
    TYPE(tddfpt_control_type)          :: t_control
    TYPE(cp_fm_p_type), DIMENSION(:), &
         POINTER                       :: S_lumos
    TYPE(cp_logger_type), POINTER      :: logger

    NULLIFY(S_lumos, logger)
    failure = .FALSE.
    logger    => cp_error_get_logger(error)
    ionode    = (logger%para_env%mepos==logger%para_env%source)
    IF (ionode) WRITE (*,*)
    IF (ionode) WRITE (*,*)

    t_control = qs_env%dft_control%tddfpt_control
    n_ev      = t_control%n_ev
    n_spins   = qs_env%dft_control%nspins

    ALLOCATE(S_lumos(n_spins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO spin=1, n_spins
       nrows(spin) = t_control%lumos(spin)%matrix%matrix_struct%nrow_global
       nhomos(spin) = t_env%evecs(1,spin)%matrix%matrix_struct%ncol_global
       nlumos(spin) = t_control%lumos(spin)%matrix%matrix_struct%ncol_global
       CALL cp_fm_create(S_lumos(spin)%matrix, t_control%lumos(spin)%matrix%matrix_struct, &
                         "S times lumos", error=error)
       CALL cp_sm_fm_multiply(qs_env%matrix_s(spin)%matrix, t_control%lumos(spin)%matrix, &
                              S_lumos(spin)%matrix, nlumos(spin), 1.0_dp, 0.0_dp)
    END DO
       
    DO spin=1, n_spins
       ALLOCATE(homo_coeff_col(nrows(spin),1),lumo_coeff_col(nrows(spin),1), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO

    DO i=1, n_ev
       IF (ionode) THEN
          WRITE (*,'(A,I3,5X,F15.6)') "  excited state : ", i, t_env%evals(i)*evolt
          WRITE (*,*)
       END IF
       summed_contributions = 0.0_dp
       DO spin=1, n_spins
          searchloop: DO occ=nhomos(spin), 1, -1
             CALL cp_fm_get_submatrix(t_env%evecs(i,spin)%matrix, homo_coeff_col, &
                                      1, occ, nrows(spin), 1, error=error)
             DO virt=1, nlumos(spin)
                CALL cp_fm_get_submatrix(S_lumos(spin)%matrix, lumo_coeff_col, &
                                         1, virt, nrows(spin), 1, error=error)
                contribution = 0.0_dp
                DO j=1, nrows(spin)
                   contribution = contribution + homo_coeff_col(j,spin)*lumo_coeff_col(j,spin)
                END DO
                summed_contributions = summed_contributions + (contribution)**2
                IF (ABS(contribution) > 1.0e-1_dp) THEN
                   IF (ionode) WRITE (*,'(18X,I3,A,I3,10X,F8.3,5X,F8.3)') &
                        occ, " ->", nhomos(spin)+virt, ABS(contribution), summed_contributions
                END IF
                IF (ABS(summed_contributions - 1.0_dp) < 1.0e-2_dp) CYCLE searchloop
             END DO
          END DO searchloop
       END DO
       IF (ionode) WRITE (*,*)
    END DO

    DO spin=1, n_spins
       DEALLOCATE(homo_coeff_col, lumo_coeff_col, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_fm_release(S_lumos(spin)%matrix, error=error)
    END DO

    DEALLOCATE(S_lumos, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE find_contributions

! *****************************************************************************
  
  SUBROUTINE normalize(X, tmp_vec, metric)

    ! ARGUMENTS
    TYPE(cp_fm_p_type), DIMENSION(:)         :: x, tmp_vec
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: metric

    ! LOCALS
    INTEGER                                  :: n_spins, spin
    REAL(KIND=dp)                           :: norm, tmp


    n_spins = SIZE(x)
    norm = 0.0_dp

    DO spin=1, n_spins
       tmp = 0.0_dp
       CALL cp_sm_fm_multiply(metric(spin)%matrix, X(spin)%matrix, &
            tmp_vec(spin)%matrix, &
            X(spin)%matrix%matrix_struct%ncol_global, &
            1.0_dp, 0.0_dp)
       CALL cp_fm_trace(X(spin)%matrix, tmp_vec(spin)%matrix, tmp)
       norm = norm + tmp
    END DO

    norm = SQRT(norm)
    DO spin=1, n_spins
       CALL cp_fm_scale((1.0_dp/norm), X(spin)%matrix)
    END DO

  END SUBROUTINE normalize

! *****************************************************************************
  
  !---------------------------------------!
  ! x must not be changed in this routine !
  ! tmp_vec may be changed                !
  !---------------------------------------!
  SUBROUTINE calc_norm(norm, X, metric, tmp_vec)


    TYPE(cp_fm_p_type), DIMENSION(:)         :: x, tmp_vec
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: metric
    REAL(KIND=dp), INTENT(OUT)              :: norm

    INTEGER :: spin

    CALL dot_prod(norm, X, X, metric, tmp_vec)
    norm = SQRT(norm)

  END SUBROUTINE calc_norm

! *****************************************************************************

  SUBROUTINE reorthogonalize(X, V_set, SV_set, work, n)

    TYPE(cp_fm_p_type), DIMENSION(:)         :: X
    TYPE(cp_fm_p_type), DIMENSION(:, :)      :: V_set
    TYPE(cp_fm_p_type), DIMENSION(:, :)      :: SV_set
    TYPE(cp_fm_p_type), DIMENSION(:)         :: work
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(LEN=*), PARAMETER :: routineN = "reorthogonalize", &
                                   routineP = moduleN//"/"//routineN

    INTEGER                                  :: i, n_spins, spin, time_handle
    REAL(DP)                                 :: dot_product, tmp

    CALL timeset(routineN,"I","",time_handle)

    IF (n>0) THEN

       n_spins = SIZE(X)
       DO spin=1, n_spins
          CALL cp_fm_to_fm(X(spin)%matrix, work(spin)%matrix)
       END DO
       
       DO i=1, n
          dot_product = 0.0_dp
          DO spin=1, n_spins
             CALL cp_fm_trace(SV_set(i,spin)%matrix, work(spin)%matrix, tmp)
             dot_product = dot_product + tmp
          END DO
          DO spin=1, n_spins
             CALL cp_fm_scale_and_add(1.0_dp, X(spin)%matrix, &
                                      -1.0_dp*dot_product, V_set(i,spin)%matrix)
          END DO
       END DO
       
    END IF

    CALL timestop(0.0_dp, time_handle)

  END SUBROUTINE reorthogonalize

! *****************************************************************************

  SUBROUTINE dot_prod(res, X, Y, metric, tmp_vec)

    REAL(KIND=dp), INTENT(OUT)      :: res
    TYPE(cp_fm_p_type), DIMENSION(:) :: X, Y, tmp_vec
    TYPE(real_matrix_p_type), DIMENSION(:), &
         POINTER                     :: metric

    CHARACTER(LEN=*), PARAMETER :: routineN = "dot_prod", &
                                   routineP = moduleN//"/"//routineN

    INTEGER                          :: spin, n_spins
    REAL(KIND=dp)                   :: tmp

    n_spins = SIZE(X)
    res = 0.0_dp
    DO spin=1, n_spins
       tmp = 0.0_dp
       CALL cp_sm_fm_multiply(metric(spin)%matrix, Y(spin)%matrix, &
                              tmp_vec(spin)%matrix, &
                              Y(spin)%matrix%matrix_struct%ncol_global, &
                              1.0_dp, 0.0_dp)
       CALL cp_fm_trace(X(spin)%matrix, tmp_vec(spin)%matrix, tmp)
       res = res + tmp
    END DO

  END SUBROUTINE dot_prod

! *****************************************************************************

END MODULE qs_tddfpt_utils
