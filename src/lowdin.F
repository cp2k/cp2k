!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K group                                       !
!-----------------------------------------------------------------------------!


!!****** cp2k/lowdin [1.0] *
!!
!!   NAME
!!     lowdin
!!
!!   FUNCTION
!!     compute lowdin charges
!!     sum(mu in C, (Shalf P Shalf)mumu)
!!
!!   AUTHOR
!!     Joost VandeVondele, Feb 2006
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE lowdin
  USE atomic_charges,                  ONLY: print_atomic_charges
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
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE sparse_matrix_types,             ONLY: get_matrix_info,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lowdin'

! *** Public subroutines ***

  PUBLIC :: lowdin_charges

CONTAINS
!!****f* lowdin/compute_charges *
!!
!!   NAME
!!     compute_charges
!!
!!   FUNCTION
!!     compute the lowdin charges
!!
!!   NOTES
!!     charges are computed per spin in the LSD case
!!     no effort to do this in a computational very efficient way.
!!
!!   INPUTS
!!    - p_matrix, s_matrix, para_env
!!    - charges : previously allocated with the right size (natom,nspin)
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE compute_charges(s_matrix,mos,charges)
    TYPE(real_matrix_type), POINTER          :: s_matrix
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges

    INTEGER :: iblock, icol, irow, ispin, n_dependent, nblock, ncol_global, &
      ncol_local, nrow_global, nrow_local, nspin
    INTEGER, DIMENSION(:), POINTER           :: col_indices, first_row, &
                                                last_row, row_indices
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_shalf, matrix_shalfc, &
                                                matrix_work, mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env

    nspin=SIZE(mos)
    charges=0.0_dp

    ! compute S**(1./2.) 
    CALL get_mo_set(mo_set=mos(1)%mo_set,mo_coeff=mo_coeff)
    CALL cp_fm_get_info(mo_coeff,para_env=para_env,context=context,nrow_global=nrow_global)
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                                   nrow_global=nrow_global,ncol_global=nrow_global)
    CALL cp_fm_create(matrix_shalf, fm_struct_tmp,name="matrix_shalf")
    CALL cp_fm_create(matrix_work, fm_struct_tmp,name="matrix_work")
    CALL cp_fm_struct_release(fm_struct_tmp)
    CALL copy_sm_to_fm(s_matrix,matrix_shalf)
    CALL cp_fm_power(matrix_shalf,matrix_work,0.5_dp,EPSILON(0.0_dp),n_dependent)
    IF (n_dependent/=0) THEN 
       ! ... bad luck, I'd say
    ENDIF 

    DO ispin=1,nspin
      ! compute S**(1./2.)*C
      CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL cp_fm_get_info(mo_coeff,para_env=para_env, &
                          context=context,nrow_global=nrow_global,ncol_global=ncol_global)
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                                   nrow_global=nrow_global,ncol_global=ncol_global)
      CALL cp_fm_create(matrix_shalfc,fm_struct_tmp,name="matrix_shalfc")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_gemm("N","N",nrow_global,ncol_global,nrow_global, &
                           1.0_dp,matrix_shalf,mo_coeff,0.0_dp,matrix_shalfc)
      ! and now the matrix (shalf*C)(shalf*C)^T only the trace is really needed, but well.
      CALL cp_fm_gemm("N","T",nrow_global,nrow_global,ncol_global, &
                           1.0_dp,matrix_shalfc,matrix_shalfc,0.0_dp,matrix_work)
      CALL cp_fm_release(matrix_shalfc) 
      ! Now condense on atomic charges
      CALL cp_fm_get_info(matrix_work,row_indices=row_indices,col_indices=col_indices, &
                                      local_data=local_data,nrow_local=nrow_local,ncol_local=ncol_local)
      CALL get_matrix_info(s_matrix,nblock_row=nblock,first_row=first_row,last_row=last_row)
      DO icol=1,ncol_local
       DO irow=1,nrow_local
        IF (row_indices(irow).EQ.col_indices(icol)) THEN
         DO iblock=1,nblock
            IF (first_row(iblock).LE.row_indices(irow) .AND. &
                last_row(iblock).GE.row_indices(irow)) THEN
               charges(iblock,ispin)=charges(iblock,ispin)+local_data(irow,icol)
            ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDDO
    ENDDO

    CALL mp_sum(charges,para_env%group)

    IF (nspin.EQ.1) charges=2.0_dp*charges

    CALL cp_fm_release(matrix_work)
    CALL cp_fm_release(matrix_shalf)

  END SUBROUTINE compute_charges

!!****f* lowdin/lowdin_charges *
!!
!!   NAME
!!     lowdin_charges
!!
!!   FUNCTION
!!     print the lowdin charges to scr on ionode
!!
!!   NOTES
!!
!!   INPUTS
!!    - p_matrix, s_matrix, para_env
!!    - particle_set (needed for Z)
!!    - scr : unit for output
!!    - ionode : true if scr is connected, and output is needed
!!
!!   MODIFICATION HISTORY
!!     02.2006 Joost VandeVondele
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE lowdin_charges(s_matrix,mos,particle_set, &
                              scr,title)
    TYPE(real_matrix_type), POINTER          :: s_matrix
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title

    INTEGER                                  :: handle, i, nblock, nspin
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges

! GO

    CALL timeset("lowdin_charges","I","",handle)

    ! here we get the numbers for charges
    nspin=SIZE(mos)
    CALL get_matrix_info(s_matrix,nblock_row=nblock)

    ALLOCATE(charges(nblock,nspin))
    CALL compute_charges(s_matrix,mos,charges)

    CALL print_atomic_charges(particle_set,scr,title,electronic_charges=charges)

    DEALLOCATE(charges)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE

END MODULE lowdin
