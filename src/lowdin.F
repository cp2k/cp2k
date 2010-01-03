!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief compute lowdin charges
!>      sum(mu in C, (Shalf P Shalf)mumu)
!> \author Joost VandeVondele, Feb 2006
! *****************************************************************************
MODULE lowdin
  USE atomic_charges,                  ONLY: print_atomic_charges
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: cp_fm_power
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_get_info,&
                                             dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop
  USE dbcsr_types,                     ONLY: dbcsr_iterator,&
                                             dbcsr_obj
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lowdin'

! *** Public subroutines ***

  PUBLIC :: lowdin_charges

CONTAINS

! *****************************************************************************
!> \brief compute the lowdin charges
!> \param p_matrix , s_matrix, para_env
!> \param charges previously allocated with the right size (natom,nspin)
!> \note
!>      charges are computed per spin in the LSD case
!>      no effort to do this in a computational very efficient way.
!> \par History
!>      02.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE compute_charges(s_matrix,mos,charges,error)
    TYPE(dbcsr_obj), POINTER                 :: s_matrix
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: blk, i, icol, irow, ispin, &
                                                n_dependent, ncol_global, &
                                                ncol_local, nrow_global, &
                                                nrow_local, nspin
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: DATA, local_data
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_shalf, matrix_shalfc, &
                                                matrix_work, mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_iterator)                     :: iter

    nspin=SIZE(mos)
    charges=0.0_dp

    ! compute S**(1./2.)
    CALL get_mo_set(mo_set=mos(1)%mo_set,mo_coeff=mo_coeff)
    CALL cp_fm_get_info(mo_coeff,para_env=para_env,context=context,nrow_global=nrow_global,error=error)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                                   nrow_global=nrow_global,ncol_global=nrow_global,error=error)
    CALL cp_fm_create(matrix_shalf, fm_struct_tmp,name="matrix_shalf",error=error)
    CALL cp_fm_create(matrix_work, fm_struct_tmp,name="matrix_work",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL copy_dbcsr_to_fm(s_matrix,matrix_shalf,error=error)
    CALL cp_fm_power(matrix_shalf,matrix_work,0.5_dp,EPSILON(0.0_dp),n_dependent,error=error)
    IF (n_dependent/=0) THEN
       ! ... bad luck, I'd say
    ENDIF

    DO ispin=1,nspin
      ! compute S**(1./2.)*C
      CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL cp_fm_get_info(mo_coeff,para_env=para_env, &
                          context=context,nrow_global=nrow_global,ncol_global=ncol_global,error=error)
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                                   nrow_global=nrow_global,ncol_global=ncol_global,error=error)
      CALL cp_fm_create(matrix_shalfc,fm_struct_tmp,name="matrix_shalfc",error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)
      CALL cp_fm_gemm("N","N",nrow_global,ncol_global,nrow_global, &
                           1.0_dp,matrix_shalf,mo_coeff,0.0_dp,matrix_shalfc,error=error)
      ! and now the matrix (shalf*C)(shalf*C)^T only the trace is really needed, but well.
      CALL cp_fm_gemm("N","T",nrow_global,nrow_global,ncol_global, &
                           1.0_dp,matrix_shalfc,matrix_shalfc,0.0_dp,matrix_work,error=error)
      CALL cp_fm_release(matrix_shalfc,error=error)
      ! Now condense on atomic charges
      CALL cp_fm_get_info(matrix_work,row_indices=row_indices,col_indices=col_indices, &
                                      local_data=local_data,nrow_local=nrow_local,ncol_local=ncol_local,error=error)

      CALL dbcsr_iterator_start(iter, s_matrix)
      DO WHILE (dbcsr_iterator_blocks_left (iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, DATA, tr, blk)
         IF(irow.EQ.icol) THEN
            !DO iblock=1,nblock
            DO i = 1, SIZE(DATA,1)
               charges(irow,ispin)=charges(irow,ispin)+DATA(i,i)
            ENDDO
         ENDIF
      ENDDO
      CALL dbcsr_iterator_stop(iter)

    ENDDO

    CALL mp_sum(charges,para_env%group)

    IF (nspin.EQ.1) charges=2.0_dp*charges

    CALL cp_fm_release(matrix_work,error=error)
    CALL cp_fm_release(matrix_shalf,error=error)

  END SUBROUTINE compute_charges

! *****************************************************************************
!> \brief print the lowdin charges to scr on ionode
!> \param p_matrix , s_matrix, para_env
!> \param particle_set (needed for Z)
!> \param scr unit for output
!> \param ionode true if scr is connected, and output is needed
!> \par History
!>      02.2006 Joost VandeVondele
! *****************************************************************************
  SUBROUTINE lowdin_charges(s_matrix,mos,particle_set, &
                              scr,title,error)
    TYPE(dbcsr_obj), POINTER                 :: s_matrix
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lowdin_charges', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nblock, nspin
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges

    CALL timeset(routineN,handle)

    ! here we get the numbers for charges
    nspin=SIZE(mos)

    CALL dbcsr_get_info(s_matrix,nblkrows_total=nblock)

    ALLOCATE(charges(nblock,nspin))
    CALL compute_charges(s_matrix,mos,charges,error=error)

    CALL print_atomic_charges(particle_set,scr,title,electronic_charges=charges)

    DEALLOCATE(charges)
    CALL timestop(handle)

  END SUBROUTINE

END MODULE lowdin
