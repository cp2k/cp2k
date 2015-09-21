!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains ADMM methods which only require the density matrix
!> \par History
!>      11.2014 created [Ole Schuett]
!> \author Ole Schuett
! *****************************************************************************
MODULE admm_dm_methods
  USE admm_dm_types,                   ONLY: admm_dm_type,&
                                             mcweeny_history_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_deallocate_matrix_set, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_get_block_p, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_scale, cp_dbcsr_set, &
       cp_dbcsr_type
  USE input_constants,                 ONLY: do_admm_basis_projection,&
                                             do_admm_blocked_projection
  USE iterate_matrix,                  ONLY: invert_Hotelling
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_env_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_set,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: admm_dm_calc_rho_aux, admm_dm_merge_ks_matrix

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_dm_methods'

  CONTAINS

! *****************************************************************************
!> \brief Entry methods: Calculates auxilary density matrix from primary one.
!> \param ks_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE admm_dm_calc_rho_aux(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CHARACTER(len=*), PARAMETER :: routineN = 'admm_dm_calc_rho_aux', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(admm_dm_type), POINTER              :: admm_dm

    NULLIFY(admm_dm)
    CALL timeset(routineN,handle)
    CALL get_ks_env(ks_env, admm_dm=admm_dm)

    SELECT CASE(admm_dm%method)
      CASE(do_admm_basis_projection)
        CALL map_dm_projection(ks_env)

      CASE(do_admm_blocked_projection)
        CALL map_dm_blocked(ks_env)

      CASE DEFAULT
         CPABORT("admm_dm_calc_rho_aux: unknown method")
    END SELECT

    IF(admm_dm%purify)&
       CALL purify_mcweeny(ks_env)

    CALL update_rho_aux(ks_env)

    CALL timestop(handle)
 END SUBROUTINE admm_dm_calc_rho_aux

! *****************************************************************************
!> \brief Entry methods: Merges auxilary Kohn-Sham matrix into primary one.
!> \param ks_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE admm_dm_merge_ks_matrix(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_dm_merge_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_merge

    CALL timeset(routineN,handle)
    NULLIFY(admm_dm, matrix_ks_merge)

    CALL get_ks_env(ks_env, admm_dm=admm_dm)

    IF(admm_dm%purify) THEN
        CALL revert_purify_mcweeny(ks_env, matrix_ks_merge)
    ELSE
        CALL get_ks_env(ks_env, matrix_ks_aux_fit=matrix_ks_merge)
    ENDIF

    SELECT CASE(admm_dm%method)
      CASE(do_admm_basis_projection)
        CALL merge_dm_projection(ks_env, matrix_ks_merge)

      CASE(do_admm_blocked_projection)
        CALL merge_dm_blocked(ks_env, matrix_ks_merge)

      CASE DEFAULT
         CPABORT("admm_dm_merge_ks_matrix: unknown method")
    END SELECT

    IF(admm_dm%purify)&
       CALL cp_dbcsr_deallocate_matrix_set(matrix_ks_merge)

    CALL timestop(handle)

  END SUBROUTINE admm_dm_merge_ks_matrix

! *****************************************************************************
!> \brief Calculates auxilary density matrix via basis projection.
!> \param ks_env ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE map_dm_projection(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CHARACTER(len=*), PARAMETER :: routineN = 'map_dm_projection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: s_mstruct_changed
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux, matrix_s_mixed, &
                                                rho_ao, rho_ao_aux
    TYPE(cp_dbcsr_type)                      :: matrix_s_aux_inv, matrix_tmp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux

    NULLIFY(dft_control, admm_dm, matrix_s_aux, matrix_s_mixed, rho, rho_aux)
    NULLIFY(rho_ao, rho_ao_aux)

    CALL get_ks_env(ks_env,&
                    admm_dm=admm_dm,&
                    dft_control=dft_control,&
                    matrix_s_aux_fit=matrix_s_aux,&
                    matrix_s_aux_fit_vs_orb=matrix_s_mixed,&
                    s_mstruct_changed=s_mstruct_changed,&
                    rho=rho,&
                    rho_aux_fit=rho_aux)

    CALL qs_rho_get(rho, rho_ao=rho_ao)
    CALL qs_rho_get(rho_aux, rho_ao=rho_ao_aux)

    IF(s_mstruct_changed) THEN
        ! Calculate A = S_aux^(-1) * S_mixed
        CALL cp_dbcsr_init(matrix_s_aux_inv)
        CALL cp_dbcsr_create(matrix_s_aux_inv, template=matrix_s_aux(1)%matrix, matrix_type="N")
        CALL invert_Hotelling(matrix_s_aux_inv, matrix_s_aux(1)%matrix, admm_dm%eps_filter)

        IF(.NOT. ASSOCIATED(admm_dm%matrix_A)) THEN
           ALLOCATE(admm_dm%matrix_A)
           CALL cp_dbcsr_init(admm_dm%matrix_A)
           CALL cp_dbcsr_create(admm_dm%matrix_A, template=matrix_s_mixed(1)%matrix, matrix_type="N")
        ENDIF
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s_aux_inv, matrix_s_mixed(1)%matrix, &
                               0.0_dp, admm_dm%matrix_A)
        CALL cp_dbcsr_release(matrix_s_aux_inv)
    ENDIF

    ! Calculate P_aux = A * P * A^T
    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp, template=admm_dm%matrix_A)
    DO ispin=1, dft_control%nspins
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, admm_dm%matrix_A, rho_ao(ispin)%matrix,&
                           0.0_dp, matrix_tmp)
        CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, admm_dm%matrix_A, &
                           0.0_dp, rho_ao_aux(ispin)%matrix)
    END DO
    CALL cp_dbcsr_release(matrix_tmp)

 END SUBROUTINE map_dm_projection

! *****************************************************************************
!> \brief Calculates auxilary density matrix via blocking.
!> \param ks_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE map_dm_blocked(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    INTEGER                                  :: blk, iatom, ispin, jatom
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block, sparse_block_aux
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao, rho_ao_aux
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux

    NULLIFY(dft_control, admm_dm, rho, rho_aux, rho_ao, rho_ao_aux)

    CALL get_ks_env(ks_env,&
                    admm_dm=admm_dm,&
                    dft_control=dft_control,&
                    rho=rho,&
                    rho_aux_fit=rho_aux)

    CALL qs_rho_get(rho, rho_ao=rho_ao)
    CALL qs_rho_get(rho_aux, rho_ao=rho_ao_aux)

    ! ** set blocked density matrix to 0
    DO ispin=1, dft_control%nspins
       CALL cp_dbcsr_set(rho_ao_aux(ispin)%matrix, 0.0_dp)
       ! ** now loop through the list and copy corresponding blocks
       CALL cp_dbcsr_iterator_start(iter, rho_ao(ispin)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
         IF(admm_dm%block_map(iatom,jatom) == 1) THEN
           CALL cp_dbcsr_get_block_p(rho_ao_aux(ispin)%matrix,&
                       row=iatom,col=jatom,BLOCK=sparse_block_aux,found=found)
           IF(found)&
              sparse_block_aux = sparse_block
         END IF
       END DO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDDO

  END SUBROUTINE map_dm_blocked

! *****************************************************************************
!> \brief Call calculate_rho_elec() for auxilary density
!> \param ks_env ...
! *****************************************************************************
  SUBROUTINE update_rho_aux(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    INTEGER                                  :: ispin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r_aux
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao_aux
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g_aux, rho_r_aux
    TYPE(qs_rho_type), POINTER               :: rho_aux

    NULLIFY(dft_control,admm_dm,rho_aux,rho_ao_aux,rho_r_aux,rho_g_aux,tot_rho_r_aux)

    CALL get_ks_env(ks_env,&
                    admm_dm=admm_dm,&
                    dft_control=dft_control,&
                    rho_aux_fit=rho_aux)

    CALL qs_rho_get(rho_aux,&
                    rho_ao=rho_ao_aux,&
                    rho_r=rho_r_aux,&
                    rho_g=rho_g_aux,&
                    tot_rho_r=tot_rho_r_aux)

    DO ispin=1, dft_control%nspins
       CALL calculate_rho_elec(ks_env=ks_env,&
                               matrix_p=rho_ao_aux(ispin)%matrix,&
                               rho=rho_r_aux(ispin),&
                               rho_gspace=rho_g_aux(ispin),&
                               total_rho=tot_rho_r_aux(ispin),&
                               soft_valid=.FALSE.,&
                               basis_type="AUX_FIT")
    END DO

    CALL qs_rho_set(rho_aux, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

 END SUBROUTINE update_rho_aux

! *****************************************************************************
!> \brief Merges auxilary Kohn-Sham matrix via basis projection.
!> \param ks_env ...
!> \param matrix_ks_merge Input: The KS matrix to be merged
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE merge_dm_projection(ks_env, matrix_ks_merge)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_merge

    INTEGER                                  :: ispin
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type)                      :: matrix_tmp
    TYPE(dft_control_type), POINTER          :: dft_control

    NULLIFY(admm_dm, dft_control, matrix_ks)

    CALL get_ks_env(ks_env,&
                    admm_dm=admm_dm,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks)

    ! Calculate K += A^T * K_aux * A
    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp, template=admm_dm%matrix_A, matrix_type="N")

    DO ispin=1, dft_control%nspins
        CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ks_merge(ispin)%matrix, admm_dm%matrix_A, &
                           0.0_dp, matrix_tmp)
        CALL cp_dbcsr_multiply("T", "N", 1.0_dp, admm_dm%matrix_A, matrix_tmp, &
                           1.0_dp, matrix_ks(ispin)%matrix)
    END DO

    CALL cp_dbcsr_release(matrix_tmp)

  END SUBROUTINE merge_dm_projection

! *****************************************************************************
!> \brief Merges auxilary Kohn-Sham matrix via blocking.
!> \param ks_env ...
!> \param matrix_ks_merge Input: The KS matrix to be merged
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE merge_dm_blocked(ks_env, matrix_ks_merge)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_merge

    INTEGER                                  :: blk, iatom, ispin, jatom
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(dft_control_type), POINTER          :: dft_control

    NULLIFY(admm_dm, dft_control, matrix_ks)

    CALL get_ks_env(ks_env,&
                    admm_dm=admm_dm,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks)

    DO ispin=1, dft_control%nspins
       CALL cp_dbcsr_iterator_start(iter, matrix_ks_merge(ispin)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
          IF(admm_dm%block_map(iatom,jatom) == 0) &
            sparse_block = 0.0_dp
       END DO
       CALL cp_dbcsr_iterator_stop(iter)
       CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_ks_merge(ispin)%matrix, 1.0_dp, 1.0_dp)
    ENDDO

  END SUBROUTINE merge_dm_blocked

! *****************************************************************************
!> \brief Apply McWeeny purification to auxilary density matrix
!> \param ks_env ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE purify_mcweeny(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, istep, nspins, &
                                                unit_nr
    REAL(KIND=dp)                            :: frob_norm
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, rho_ao_aux
    TYPE(cp_dbcsr_type)                      :: matrix_ps, matrix_psp, &
                                                matrix_test
    TYPE(cp_dbcsr_type), POINTER             :: matrix_p, matrix_s
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mcweeny_history_type), POINTER      :: history, new_hist_entry
    TYPE(qs_rho_type), POINTER               :: rho_aux_fit

    CALL timeset(routineN,handle)
    NULLIFY(dft_control, admm_dm, matrix_s_aux_fit, rho_aux_fit, new_hist_entry,&
            matrix_p, matrix_s,rho_ao_aux)

    unit_nr = cp_logger_get_default_unit_nr()
    CALL get_ks_env(ks_env,&
                    dft_control=dft_control,&
                    admm_dm=admm_dm,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    rho_aux_fit=rho_aux_fit)

    CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux)

    matrix_p => rho_ao_aux(1)%matrix
    CALL cp_dbcsr_init(matrix_PS)
    CALL cp_dbcsr_create(matrix_PS,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_PSP)
    CALL cp_dbcsr_create(matrix_PSP,template=matrix_p,matrix_type="S")
    CALL cp_dbcsr_init(matrix_test)
    CALL cp_dbcsr_create(matrix_test,template=matrix_p,matrix_type="S")

    nspins = dft_control%nspins
    DO ispin=1, nspins
       matrix_p => rho_ao_aux(ispin)%matrix
       matrix_s => matrix_s_aux_fit(1)%matrix
       history => admm_dm%mcweeny_history(ispin)%p
       IF(ASSOCIATED(history)) STOP "purify_dm_mcweeny: history already associated"
       IF(nspins==1) CALL cp_dbcsr_scale(matrix_p,0.5_dp)

       DO istep=1, admm_dm%mcweeny_max_steps
          ! allocate new element in linked list
           ALLOCATE(new_hist_entry)
           new_hist_entry%next => history
           history => new_hist_entry
           history%count = istep
           NULLIFY(new_hist_entry)
           CALL cp_dbcsr_init(history%m)
           CALL cp_dbcsr_create(history%m, template=matrix_p ,matrix_type="N")
           CALL cp_dbcsr_copy(history%m, matrix_p, name="P from McWeeny")

           ! calc PS and PSP
           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p, matrix_s,&
                0.0_dp, matrix_ps)

           CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_ps, matrix_p,&
                0.0_dp, matrix_psp)

           !test convergence
           CALL cp_dbcsr_copy(matrix_test, matrix_psp)
           CALL cp_dbcsr_add(matrix_test, matrix_p, 1.0_dp,-1.0_dp)
           frob_norm=cp_dbcsr_frobenius_norm(matrix_test)
           IF(unit_nr>0) WRITE(unit_nr,'(t3,a,i5,a,f16.8)')"McWeeny-Step",istep,&
                               ": Deviation of idempotency", frob_norm
           IF(frob_norm < 1000_dp*admm_dm%eps_filter .AND. istep>1) EXIT

           ! build next P matrix
           CALL cp_dbcsr_copy(matrix_p, matrix_PSP, name="P from McWeeny")
           CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_PS, matrix_PSP,&
                3.0_dp, matrix_p)
       END DO
       admm_dm%mcweeny_history(ispin)%p => history
       IF(nspins==1) CALL cp_dbcsr_scale(matrix_p,2.0_dp)
    END DO

    ! clean up
    CALL cp_dbcsr_release(matrix_PS)
    CALL cp_dbcsr_release(matrix_PSP)
    CALL cp_dbcsr_release(matrix_test)
    CALL timestop(handle)
  END SUBROUTINE purify_mcweeny

! *****************************************************************************
!> \brief Prepare auxilary KS-matrix for merge using reverse McWeeny
!> \param ks_env ...
!> \param matrix_ks_merge Output: The KS matrix for the merge
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE revert_purify_mcweeny(ks_env, matrix_ks_merge)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_merge

    CHARACTER(LEN=*), PARAMETER :: routineN = 'revert_purify_mcweeny', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins, unit_nr
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mcweeny_history_type), POINTER      :: history_curr, history_next

    CALL timeset(routineN,handle)
    unit_nr = cp_logger_get_default_unit_nr()
    NULLIFY(admm_dm, dft_control, matrix_ks, matrix_ks_aux_fit,&
            matrix_s_aux_fit, matrix_s_aux_fit_vs_orb,&
            history_next, history_curr,  matrix_k)

    CALL get_ks_env(ks_env,&
                    admm_dm=admm_dm,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb)

    nspins = dft_control%nspins
    ALLOCATE(matrix_ks_merge(nspins))

    DO ispin=1, nspins
       ALLOCATE(matrix_ks_merge(ispin)%matrix)
       matrix_k => matrix_ks_merge(ispin)%matrix
       CALL cp_dbcsr_init(matrix_k)
       CALL cp_dbcsr_copy(matrix_k, matrix_ks_aux_fit(ispin)%matrix, name="K")
       history_curr => admm_dm%mcweeny_history(ispin)%p
       NULLIFY(admm_dm%mcweeny_history(ispin)%p)

       ! reverse McWeeny iteration
       DO WHILE (ASSOCIATED(history_curr))
          IF(unit_nr>1) WRITE(unit_nr,'(t3,a,i5)')"Reverse McWeeny-Step ", history_curr%count
          CALL reverse_mcweeny_step(matrix_k=matrix_k,&
                                    matrix_s=matrix_s_aux_fit(1)%matrix,&
                                    matrix_p=history_curr%m)
          CALL cp_dbcsr_release(history_curr%m)
          history_next => history_curr%next
          DEALLOCATE(history_curr)
          history_curr => history_next
          NULLIFY(history_next)
       END DO

    END DO

    ! clean up
    CALL timestop(handle)

  END SUBROUTINE revert_purify_mcweeny

! *****************************************************************************
!> \brief Multiply matrix_k with partial derivative of McWeeny by reversing it.
!> \param matrix_k ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE reverse_mcweeny_step(matrix_k, matrix_s, matrix_p)
    TYPE(cp_dbcsr_type)                      :: matrix_k, matrix_s, matrix_p

    CHARACTER(LEN=*), PARAMETER :: routineN = 'reverse_mcweeny_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: matrix_ps, matrix_sp, &
                                                matrix_sum, matrix_tmp

    CALL timeset(routineN,handle)
    CALL cp_dbcsr_init(matrix_ps)
    CALL cp_dbcsr_create(matrix_ps,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_sp)
    CALL cp_dbcsr_create(matrix_sp,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_tmp)
    CALL cp_dbcsr_create(matrix_tmp,template=matrix_p,matrix_type="N")
    CALL cp_dbcsr_init(matrix_sum)
    CALL cp_dbcsr_create(matrix_sum,template=matrix_p,matrix_type="N")

    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_p, matrix_s,&
             0.0_dp, matrix_ps)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s, matrix_p,&
             0.0_dp, matrix_sp)

    !TODO: can we exploid more symmetry?
    CALL cp_dbcsr_multiply("N", "N", 3.0_dp, matrix_k, matrix_ps,&
             0.0_dp, matrix_sum)
    CALL cp_dbcsr_multiply("N", "N", 3.0_dp, matrix_sp, matrix_k,&
             1.0_dp, matrix_sum)

    !matrix_tmp = KPS
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_k, matrix_ps,&
             0.0_dp, matrix_tmp)
    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_tmp, matrix_ps,&
            1.0_dp, matrix_sum)
    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_sp, matrix_tmp,&
            1.0_dp, matrix_sum)

    !matrix_tmp = SPK
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_sp, matrix_k,&
             0.0_dp, matrix_tmp)
    CALL cp_dbcsr_multiply("N", "N", -2.0_dp, matrix_sp, matrix_tmp,&
            1.0_dp, matrix_sum)

    ! overwrite matrix_k
    CALL cp_dbcsr_copy(matrix_k, matrix_sum, name="K from reverse McWeeny")

    ! clean up
    CALL cp_dbcsr_release(matrix_sum)
    CALL cp_dbcsr_release(matrix_tmp)
    CALL cp_dbcsr_release(matrix_ps)
    CALL cp_dbcsr_release(matrix_sp)
    CALL timestop(handle)
  END SUBROUTINE reverse_mcweeny_step

END MODULE admm_dm_methods
