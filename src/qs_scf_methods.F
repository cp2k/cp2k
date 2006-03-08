!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2006  CP2K developers group
!-----------------------------------------------------------------------------!

!!****** cp2k/qs_scf_methods [1.0] *
!!
!!   NAME
!!     qs_scf_methods
!!
!!   FUNCTION
!!     groups fairly general SCF methods, so that modules other than qs_scf can use them too
!!     split off from qs_scf to reduce dependencies
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     - Joost VandeVondele (03.2006)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf_methods
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_symm
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_block_jacobi_classic,&
                                             cp_fm_power,&
                                             cp_fm_syevd,&
                                             cp_fm_syevx
  USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_mo_types,                     ONLY: correct_mo_eigenvalues,&
                                             get_mo_set,&
                                             mo_set_type,&
                                             set_mo_occupation,&
                                             set_mo_set
  USE qs_scf_types,                    ONLY: ao_mixing_nr,&
                                             no_mixing_nr
  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_program,&
                                             stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_methods'
  REAL(KIND=dp), PARAMETER             :: ratio=0.25_dp

  PUBLIC :: eigensolver, simple_eigensolver, scf_env_density_mixing

! *****************************************************************************

CONTAINS

!!****f* qs_scf/scf_env_density_mixing [1.0] *
!!
!!   NAME
!!     scf_env_density_mixing
!!
!!   FUNCTION
!!     perform (if requested) a density mixing
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE scf_env_density_mixing(p_mix_new,mixing_method,p_mix,iter_delta,&
                                    qs_env,invert,error)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_mix_new
    INTEGER, INTENT(IN)                      :: mixing_method
    REAL(KIND=dp), INTENT(IN)                :: p_mix
    REAL(KIND=dp), INTENT(INOUT)             :: iter_delta
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: invert
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_density_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure, my_invert
    REAL(KIND=dp)                            :: tmp

    CALL timeset(routineN,handle)

    failure=.FALSE.
    my_invert=.FALSE.
    IF (PRESENT(invert)) my_invert=invert

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       SELECT CASE (mixing_method)
       CASE(no_mixing_nr)
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"don't call mixing in no_mixing_nr:",error,failure)
       CASE(ao_mixing_nr)
          CPPrecondition(ASSOCIATED(p_mix_new),cp_failure_level,routineP,error,failure)
!FM       iter_delta=0.0_dp
          DO ispin=1,SIZE(p_mix_new)
             IF (my_invert) THEN
                CPPrecondition(p_mix/=0.0_dp,cp_failure_level,routineP,error,failure)
                IF (p_mix/=1.0_dp) THEN
                   CALL cp_sm_scale_and_add(matrix_a=p_mix_new(ispin)%matrix,&
                        alpha=1.0_dp/p_mix,&
                        matrix_b=qs_env%rho%rho_ao(ispin)%matrix,&
                        beta=(p_mix-1.0_dp)/p_mix,&
                        error=error)
                END IF
             ELSE
                CALL cp_sm_mix(m1=p_mix_new(ispin)%matrix,&
                     m2=qs_env%rho%rho_ao(ispin)%matrix,&
                     p_mix=p_mix,&
                     delta=tmp,&
                     para_env=qs_env%para_env,&
                     error=error)
                iter_delta=MAX(iter_delta,tmp)
             END IF
          END DO
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"unknown density mixing method:"//&
               cp_to_string(mixing_method),error,failure)
       END SELECT
    END IF

    CALL timestop(handle)
  END SUBROUTINE scf_env_density_mixing


! *****************************************************************************

  SUBROUTINE eigensolver(matrix_ks,mo_set,ortho,work,do_level_shift,level_shift,&
       use_cholesky, work_syevx, use_jacobi,jacobi_threshold,smear,&
       error)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.
!   ks will be modified

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: matrix_ks
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: ortho, work
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift
    LOGICAL, INTENT(INOUT)                   :: use_cholesky
    REAL(KIND=dp), INTENT(IN)                :: work_syevx
    LOGICAL, INTENT(IN)                      :: use_jacobi
    REAL(KIND=dp), INTENT(IN)                :: jacobi_threshold, smear
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, homo, imo, lfomo, &
                                                nao, nelectron, nmo
    LOGICAL                                  :: mo_uocc
    REAL(KIND=dp)                            :: alpha, maxocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, mo_occupation
    TYPE(cp_fm_type), POINTER                :: mo_coeff

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I"," ",handle)

!   *** Diagonalise the Kohn-Sham matrix ***

    CALL get_mo_set(mo_set=mo_set,&
         nao=nao,&
         nmo=nmo,&
         homo=homo, &
         nelectron=nelectron,&
         maxocc=maxocc,&
         eigenvalues=mo_eigenvalues,&
         occupation_numbers=mo_occupation,&
         mo_coeff=mo_coeff)

    IF (do_level_shift) THEN
       IF (use_cholesky) THEN
          CALL stop_program("eigensolver","level shift not implemented")
       ENDIF
    ENDIF

    IF (use_cholesky) THEN
       CALL cp_fm_cholesky_reduce(matrix_ks,ortho)

       IF (use_jacobi) THEN
          CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues)
          use_cholesky = .FALSE.
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
          CALL cp_fm_to_fm(mo_coeff,ortho)
       ELSE IF (nmo==nao) THEN
          CALL cp_fm_syevd(matrix_ks,work,mo_eigenvalues)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
       ELSE
          CALL cp_fm_syevx(matrix_ks,work,mo_eigenvalues,nmo,work_syevx)
          CALL cp_fm_cholesky_restore(work,nmo,ortho,mo_coeff,"SOLVE")
       END IF

    ELSE IF (use_jacobi) THEN

       CALL cp_fm_symm("L","U",nao,homo,1.0_dp,matrix_ks,mo_coeff,0.0_dp,work)
       CALL cp_fm_gemm("T","N",homo,nao-homo,nao,1.0_dp,work,mo_coeff,&
            0.0_dp,matrix_ks,b_first_col=homo+1,c_first_col=homo+1)

!     epsi =1e-11
!     Anz_max=4
!     N=matrix_ks%matrix_struct%nrow_global

!     *** Scale the elements with the inverse ***
!     *** difference of the old eigenvalues   ***

!     CALL blacs_block_jacobi_scaled(matrix_ks,work,mo_eigenvalues,&
!                                    jacobi_threshold,homo+1)

!     *** Klassisches Jacobi-Block-Verfahren ***
!     *** Rotiere Elemente nur, wenn sie groesser als thresh sind ***

       CALL cp_fm_block_jacobi_classic(matrix_ks,mo_coeff,mo_eigenvalues,&
            jacobi_threshold,homo+1)

    ELSE ! full S^-1 has been computed

       CALL cp_fm_symm("L","U",nao,nao,1.0_dp,matrix_ks,ortho,0.0_dp,work)
       CALL cp_fm_gemm("T","N",nao,nao,nao,1.0_dp,ortho,work,0.0_dp,matrix_ks)

       IF (do_level_shift) THEN
          DO imo=homo+1,nmo
             ! now unnecessary broadcast of element
             CALL cp_fm_get_element(matrix_ks,imo,imo,alpha)
             alpha=alpha+level_shift
             CALL cp_fm_set_element(matrix_ks,imo,imo,alpha)
          END DO
       END IF

       CALL cp_fm_syevx(matrix_ks,work,mo_eigenvalues,nmo,work_syevx)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,ortho,work,0.0_dp,&
            mo_coeff)

       IF (do_level_shift) THEN

          CALL correct_mo_eigenvalues(mo_set,level_shift)

!     *** Use last MO set as orthogonalization matrix ***

          CALL cp_fm_to_fm(mo_coeff,ortho)

       END IF

    END IF

    mo_uocc = .TRUE.
    IF (smear /= 0.0_dp) THEN
      CALL set_mo_occupation(mo_occupation,mo_eigenvalues,&
                             homo,lfomo,maxocc,nelectron,&
                             mo_uocc,&
                             smear=smear)
      CALL set_mo_set(mo_set=mo_set,&
                      homo=homo,lfomo=lfomo,&
                      uniform_occupation=mo_uocc)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE eigensolver

! *****************************************************************************
  SUBROUTINE simple_eigensolver(matrix_ks,mo_set,work,&
             do_level_shift,level_shift,work_syevx,smear,error)

    TYPE(cp_fm_type), POINTER                :: matrix_ks
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: work
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift, work_syevx, smear
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, homo, imo, lfomo, &
                                                nao, nelectron, nmo, istat
    LOGICAL                                  :: mo_uocc
    REAL(KIND=dp)                            :: alpha, maxocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, mo_occupation
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: eigvals

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I"," ",handle)

    NULLIFY(mo_eigenvalues,mo_occupation,mo_coeff)

    CALL get_mo_set(mo_set=mo_set,&
         nao=nao,&
         nmo=nmo,&
         homo=homo, &
         nelectron=nelectron,&
         maxocc=maxocc,&
         eigenvalues=mo_eigenvalues,&
         occupation_numbers=mo_occupation,&
         mo_coeff=mo_coeff)

    IF (do_level_shift) THEN
      DO imo=homo+1,nmo
         ! now unnecessary broadcast of element
         CALL cp_fm_get_element(matrix_ks,imo,imo,alpha)
         alpha=alpha+level_shift
         CALL cp_fm_set_element(matrix_ks,imo,imo,alpha)
      END DO
    END IF
    
    ALLOCATE(eigvals(nao),STAT=istat)
    IF (istat.NE.0) CALL stop_memory("simple_eigensolver","eigvals",nao)

    IF ( REAL(nmo,dp)/REAL(nao,dp) < ratio ) THEN
      CALL cp_fm_syevx(matrix_ks,work,eigvals,nmo,work_syevx,error)
    ELSE
      CALL cp_fm_syevd(matrix_ks,work,eigvals,error)
    END IF
    mo_eigenvalues(1:nmo) = eigvals(1:nmo)
    CALL cp_fm_to_fm(work,mo_coeff,nmo,1,1)

    DEALLOCATE(eigvals,STAT=istat)
    IF (istat.NE.0) CALL stop_memory("simple_eigensolver","eigvals")

    IF (do_level_shift) THEN
      CALL correct_mo_eigenvalues(mo_set,level_shift)
    END IF

    mo_uocc = .TRUE.
    IF (smear /= 0.0_dp) THEN
      CALL set_mo_occupation(mo_occupation,mo_eigenvalues,&
                             homo,lfomo,maxocc,nelectron,&
                             mo_uocc, smear=smear)
      CALL set_mo_set(mo_set=mo_set,&
                      homo=homo,lfomo=lfomo,&
                      uniform_occupation=mo_uocc)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE simple_eigensolver

! *****************************************************************************

!!****f* qs_scf/cp_sm_mix [1.0] *
!!
!!   NAME
!!     cp_sm_mix
!!
!!   FUNCTION
!!     Perform a mixing of the given matrixes into the first matrix
!!     m1 = m2 + p_mix (m1-m2)
!!
!!   NOTES
!!     if you what to store the result in m2 swap m1 and m2 an use 
!!     (1-pmix) as pmix
!!     para_env should be removed (embedded in matrix)
!!
!!   ARGUMENTS
!!     - m1: first (new) matrix, is modified
!!     - m2: the second (old) matrix
!!     - p_mix:how much m1 is conserved (0: none, 1: all)
!!     - delta: maximum norm of m1-m2
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 rewamped [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_sm_mix(m1,m2,p_mix,delta,para_env,error)

    TYPE(real_matrix_type), POINTER          :: m1, m2
    REAL(KIND=dp), INTENT(IN)                :: p_mix
    REAL(KIND=dp), INTENT(OUT)               :: delta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_mix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iblock_col_m1, &
                                                iblock_col_m2, iblock_row, j, &
                                                nblock_row
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: r_mix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block, p_old_block
    TYPE(real_block_node_type), POINTER      :: block_node_m1, block_node_m2

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)
    failure = .FALSE.
    delta = 0.0_dp

    r_mix = 1.0_dp - p_mix

    CALL get_matrix_info(matrix=m2,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

       block_node_m1 => first_block_node(matrix=m1,&
            block_row=iblock_row)
       block_node_m2 => first_block_node(matrix=m2,&
            block_row=iblock_row)

       DO WHILE (ASSOCIATED(block_node_m1))

          CALL get_block_node(block_node=block_node_m1,&
               block_col=iblock_col_m1,&
               BLOCK=p_new_block)

          CALL get_block_node(block_node=block_node_m2,&
               block_col=iblock_col_m2,&
               BLOCK=p_old_block)
          CPPostcondition(iblock_col_m1==iblock_col_m2,cp_failure_level,routineP,error,failure)

          DO j=1,SIZE(p_new_block,2)
             DO i=1,SIZE(p_new_block,1)
                p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
                delta = MAX(delta,ABS(p_new_block(i,j)))
                p_new_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
             END DO
          END DO

          block_node_m1 => next_block_node(block_node_m1)
          block_node_m2 => next_block_node(block_node_m2)

       END DO

    END DO

    CALL mp_max(delta,para_env%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE cp_sm_mix

END MODULE qs_scf_methods
