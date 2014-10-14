!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      09-JUL-2002, TCH, development started
! *****************************************************************************
MODULE qs_tddfpt_utils
  USE cp_control_types,                ONLY: dft_control_type,&
                                             tddfpt_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: evolt
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_tddfpt_types,                 ONLY: tddfpt_env_allocate,&
                                             tddfpt_env_deallocate,&
                                             tddfpt_env_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_utils'
  LOGICAL, PARAMETER          :: DEBUG_THIS_MODULE = .TRUE.

! *****************************************************************************
  TYPE simple_solution_sorter
     INTEGER                               :: orbit
     INTEGER                               :: lumo
     REAL(KIND=DP)                        :: value
     TYPE(simple_solution_sorter), POINTER :: next
  END TYPE simple_solution_sorter

  PRIVATE

  ! METHODS
  PUBLIC :: tddfpt_cleanup,&
            tddfpt_init,&
            co_initial_guess,&
            find_contributions,&
            normalize,&
            reorthogonalize

CONTAINS

! *****************************************************************************
!> \brief Initialize some necessary structures for a tddfpt calculation.
!> \param p_env perturbation environment to be initialized
!> \param t_env tddfpt environment to be initialized
!> \param qs_env Quickstep environment with the results of a
!>                   ground state calcualtion
!> \param error ...
! *****************************************************************************
  SUBROUTINE tddfpt_init(p_env, t_env, qs_env, error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(tddfpt_env_type), INTENT(out)       :: t_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init', &
      routineP = moduleN//':'//routineN

!------------------!
! create the p_env !
!------------------!

    CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE.,error=error)
    CALL p_env_psi0_changed(p_env, qs_env,error=error) ! update the m_epsilon matrix

    !------------------!
    ! create the t_env !
    !------------------!
    CALL tddfpt_env_allocate(t_env, p_env, qs_env,error=error)
    CALL tddfpt_env_init(t_env, qs_env,error=error)

  END SUBROUTINE tddfpt_init

! *****************************************************************************
!> \brief Initialize t_env with meaningfull values.
!> \param t_env ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE tddfpt_env_init(t_env, qs_env,error)

    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_spins, spin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control

    NULLIFY(matrix_s, dft_control)

    CALL get_qs_env(qs_env, matrix_s=matrix_s, dft_control=dft_control, error=error)
    n_spins = dft_control%nspins
    IF (dft_control%tddfpt_control%invert_S) THEN
       DO spin=1, n_spins
          CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, t_env%invS(spin)%matrix,error=error)
          CALL cp_fm_cholesky_decompose(t_env%invS(spin)%matrix,error=error)
          CALL cp_fm_cholesky_invert(t_env%invS(spin)%matrix,error=error)
       END DO
    END IF

  END SUBROUTINE tddfpt_env_init

! *****************************************************************************
!> \brief ...
!> \param t_env ...
!> \param p_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE tddfpt_cleanup(t_env, p_env, error)

    TYPE(tddfpt_env_type), INTENT(inout)     :: t_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_cleanup', &
      routineP = moduleN//':'//routineN

    CALL tddfpt_env_deallocate(t_env,error=error)
    CALL p_env_release(p_env,error=error)

  END SUBROUTINE tddfpt_cleanup

! *****************************************************************************
!> \brief ...
!> \param matrices ...
!> \param energies ...
!> \param n_v ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE co_initial_guess(matrices, energies, n_v, qs_env, error)

    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: matrices
    REAL(kind=DP), DIMENSION(:), INTENT(OUT) :: energies
    INTEGER, INTENT(IN)                      :: n_v
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'co_initial_guess', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n_cols, n_lumos, n_orbits, &
                                                n_rows, n_spins, oo, spin, &
                                                stat, vo
    LOGICAL                                  :: failure
    REAL(KIND=DP)                            :: evd
    REAL(KIND=DP), ALLOCATABLE, &
      DIMENSION(:, :)                        :: guess, lumos
    REAL(KIND=DP), DIMENSION(:), POINTER     :: orbital_eigenvalues
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(simple_solution_sorter), POINTER    :: sorter_iterator, &
                                                sorter_pointer, sorter_start
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

! number of vectors to initialize

    NULLIFY(dft_control)
    failure=.FALSE.

    CALL get_qs_env(qs_env, dft_control=dft_control, error=error)
    tddfpt_control => dft_control%tddfpt_control
    n_spins = dft_control%nspins
    energies = 0.0_dp

    IF (.NOT.ASSOCIATED(tddfpt_control%lumos)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"LUMOS missing")
    END IF

    DO spin=1, n_spins

       n_cols = matrices(1,spin)%matrix%matrix_struct%ncol_global
       n_rows = matrices(1,spin)%matrix%matrix_struct%nrow_global

       DO i=1, n_v
          CALL cp_fm_set_all(matrices(i,spin)%matrix, 0.0_dp,error=error)
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
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_fm_get_submatrix(tddfpt_control%lumos(spin)%matrix, lumos, &
                                start_col=1, n_cols=n_lumos,error=error)

       !-------------------!
       ! fill the matrices !
       !-------------------!
       sorter_iterator=>sorter_start
       DO i=1, MIN(n_v,n_orbits*n_lumos)
          CALL dcopy(n_rows*n_orbits, 0.0_dp, 0, guess, 1)
          CALL dcopy(n_rows, lumos(:,sorter_iterator%lumo), 1, &
                             guess(:,sorter_iterator%orbit),1)
          CALL cp_fm_set_submatrix(matrices(i,spin)%matrix, &
                                   guess(:,1:n_orbits),error=error)
          energies(i) = energies(i) + sorter_iterator%value/REAL(n_spins,dp)
          sorter_iterator=>sorter_iterator%next
       END DO
       IF (n_v > n_orbits*n_lumos) THEN
          DO i=n_orbits*n_lumos+1, n_v
             CALL cp_fm_init_random(matrices(i,spin)%matrix, n_orbits,error=error)
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
!> \brief ...
!> \param qs_env ...
!> \param t_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE find_contributions(qs_env, t_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(tddfpt_env_type), INTENT(IN)        :: t_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_contributions', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, n_ev, n_spins, occ, &
                                                output_unit, spin, stat, virt
    INTEGER, DIMENSION(2)                    :: nhomos, nlumos, nrows
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: contribution, &
                                                summed_contributions
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: homo_coeff_col, lumo_coeff_col
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: S_lumos
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(tddfpt_control_type)                :: t_control

    NULLIFY(S_lumos, logger, matrix_s, dft_control)
    failure = .FALSE.
    logger    => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    CALL get_qs_env(qs_env, matrix_s=matrix_s, dft_control=dft_control, error=error)

    IF (output_unit>0) WRITE (output_unit,*)
    IF (output_unit>0) WRITE (output_unit,*)

    t_control = dft_control%tddfpt_control
    n_ev      = t_control%n_ev
    n_spins   = dft_control%nspins

    ALLOCATE(S_lumos(n_spins), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO spin=1, n_spins
       nrows(spin) = t_control%lumos(spin)%matrix%matrix_struct%nrow_global
       nhomos(spin) = t_env%evecs(1,spin)%matrix%matrix_struct%ncol_global
       nlumos(spin) = t_control%lumos(spin)%matrix%matrix_struct%ncol_global
       CALL cp_fm_create(S_lumos(spin)%matrix, t_control%lumos(spin)%matrix%matrix_struct, &
                         "S times lumos", error=error)
       CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, t_control%lumos(spin)%matrix, &
                              S_lumos(spin)%matrix, nlumos(spin), 1.0_dp, 0.0_dp,error=error)
    END DO

    ALLOCATE(homo_coeff_col(MAXVAL(nrows(1:n_spins)),1), &
             lumo_coeff_col(MAXVAL(nrows(1:n_spins)),1), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1, n_ev
       IF (output_unit>0) THEN
          WRITE (output_unit,'(A,I3,5X,F15.6)') "  excited state : ", i, t_env%evals(i)*evolt
          WRITE (output_unit,*)
       END IF
       summed_contributions = 0.0_dp
       DO spin=1, n_spins
          IF (n_spins == 2) THEN
             IF (spin == 1) THEN
                IF (output_unit>0) WRITE (output_unit,*) 'alpha:'
             ELSE
                IF (output_unit>0) WRITE (output_unit,*) 'beta:'
             END IF
          END IF
          searchloop: DO occ=nhomos(spin), 1, -1
             CALL cp_fm_get_submatrix(t_env%evecs(i,spin)%matrix, homo_coeff_col, &
                                      1, occ, nrows(spin), 1, error=error)
             DO virt=1, nlumos(spin)
                CALL cp_fm_get_submatrix(S_lumos(spin)%matrix, lumo_coeff_col, &
                                         1, virt, nrows(spin), 1, error=error)
                contribution = 0.0_dp
                DO j=1, nrows(spin)
                   contribution = contribution + homo_coeff_col(j,1)*lumo_coeff_col(j,1)
                END DO
                summed_contributions = summed_contributions + (contribution)**2
                IF (ABS(contribution) > 5.0e-2_dp) THEN
                   IF (output_unit>0) WRITE (output_unit,'(14X,I5,A,I5,10X,F8.3,5X,F8.3)') &
                        occ, " ->", nhomos(spin)+virt, ABS(contribution), summed_contributions
                END IF
                IF (ABS(summed_contributions - 1.0_dp) < 1.0e-3_dp) CYCLE searchloop
             END DO
          END DO searchloop
       END DO
       IF (output_unit>0) WRITE (output_unit,*)
    END DO

    !
    ! punch a checksum for the regs
    IF(output_unit > 0) THEN
       WRITE(output_unit,'(T2,A,E14.6)') ' TDDFPT : CheckSum  =',SQRT(SUM(t_env%evals**2))
    ENDIF

    DO spin=1, n_spins
       CALL cp_fm_release(S_lumos(spin)%matrix, error=error)
    END DO
    DEALLOCATE(homo_coeff_col, lumo_coeff_col, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(S_lumos, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE find_contributions

! *****************************************************************************
!> \brief ...
!> \param X ...
!> \param tmp_vec ...
!> \param metric ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE normalize(X, tmp_vec, metric, error)

    TYPE(cp_fm_p_type), DIMENSION(:)         :: x, tmp_vec
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: metric
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n_spins, spin
    REAL(KIND=dp)                            :: norm, tmp

    n_spins = SIZE(x)
    norm = 0.0_dp

    DO spin=1, n_spins
       tmp = 0.0_dp
       CALL cp_dbcsr_sm_fm_multiply(metric(1)%matrix, X(spin)%matrix, &
            tmp_vec(spin)%matrix, &
            X(spin)%matrix%matrix_struct%ncol_global, &
            1.0_dp, 0.0_dp,error=error)
       CALL cp_fm_trace(X(spin)%matrix, tmp_vec(spin)%matrix, tmp,error=error)
       norm = norm + tmp
    END DO

    norm = SQRT(norm)
    DO spin=1, n_spins
       CALL cp_fm_scale((1.0_dp/norm), X(spin)%matrix,error=error)
    END DO

  END SUBROUTINE normalize

  !---------------------------------------!
  ! x must not be changed in this routine !
  ! tmp_vec may be changed                !
  !---------------------------------------!
! *****************************************************************************
!> \brief ...
!> \param X ...
!> \param V_set ...
!> \param SV_set ...
!> \param work ...
!> \param n ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE reorthogonalize(X, V_set, SV_set, work, n, error)

    TYPE(cp_fm_p_type), DIMENSION(:)         :: X
    TYPE(cp_fm_p_type), DIMENSION(:, :)      :: V_set, SV_set
    TYPE(cp_fm_p_type), DIMENSION(:)         :: work
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reorthogonalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, n_spins, spin
    REAL(DP)                                 :: dot_product, tmp

    CALL timeset(routineN,handle)

    IF (n>0) THEN

       n_spins = SIZE(X)
       DO spin=1, n_spins
          CALL cp_fm_to_fm(X(spin)%matrix, work(spin)%matrix,error=error)
       END DO

       DO i=1, n
          dot_product = 0.0_dp
          DO spin=1, n_spins
             CALL cp_fm_trace(SV_set(i,spin)%matrix, work(spin)%matrix, tmp,error=error)
             dot_product = dot_product + tmp
          END DO
          DO spin=1, n_spins
             CALL cp_fm_scale_and_add(1.0_dp, X(spin)%matrix, &
                                      -1.0_dp*dot_product, V_set(i,spin)%matrix,error=error)
          END DO
       END DO

    END IF

    CALL timestop(handle)

  END SUBROUTINE reorthogonalize

! *****************************************************************************

END MODULE qs_tddfpt_utils
