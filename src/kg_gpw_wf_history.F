!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/kg_gpw_wf_history [1.0] *
!!
!!   NAME
!!     kg_gpw_wf_history
!!
!!   FUNCTION
!!     Routines to store the history of optimized wfn when the KG_GPW method is used.
!!     The method requires that the wfn are stored in molecular blocks and not as full matrix
!!     
!!   AUTHOR
!!     MI (01.04.2005)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_wf_history

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kg_gpw_fm_mol_methods,           ONLY: calculate_mol_density_matrix,&
                                             fm_mol_eigensolver,&
                                             fm_mol_orthonormality,&
                                             mol_make_basis,&
                                             multiply_sparse_mol_mo
  USE kg_gpw_fm_mol_types,             ONLY: get_fm_mol_block,&
                                             get_kg_fm_mol_set,&
                                             get_mol_mo_set,&
                                             fm_mol_blocks_type,&
                                             kg_fm_mol_set_type,&
                                             kg_fm_p_type,&
                                             mol_mo_set_p_type,&
                                             mol_mo_set_type
  USE kg_gpw_fm_mol_utils,             ONLY: copy_sparse2mol_block
  USE kinds,                           ONLY: dp,dp_size,int_size
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_wf_history"

! *** Public subroutines ***

  PUBLIC :: kg_gpw_fm_mol_linear, kg_gpw_fm_mol_ps, kg_gpw_prev_wf

!***************************************************************************

CONTAINS

!***************************************************************************
!!****f* kg_gpw_wf_history/kg_gpw_fm_mol_linear
!!
!!   NAME
!!    kg_gpw_fm_mol_linear
!!
!!   FUNCTION
!!     The linear extrapolation of the mos is calculated by using the two
!!     previous steps. 
!!      set_a%...%mos =  alpha*set_a%...%mos + beta*set_b%...%mos
!!     If required the density is also calculated
!!     
!!   NOTE
!!     matrix_p and matrix_s are given as input only if the density has to be calculated
!!     i.e. when all the terms have been added to the final mos
!!
!!   INPUTS
!!     - alpha : first multiplicative constant
!!     - set_a : contains the initial mos and the output mos (inout)
!!     - beta : second multiplicative constant
!!     - set_b : contains the mos of a previous step that will be added (in)
!!     - matrix_p: density matrix in real sparse matrix format, optional 
!!     - matrix_s: overlapmatrix in real sparse matrix format, optional
!!     - error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE kg_gpw_fm_mol_linear(alpha,set_a,beta,set_b,matrix_p,matrix_s,error)

    REAL(dp), INTENT(IN), OPTIONAL            :: alpha, beta
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                 :: set_a
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                 :: set_b
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER, OPTIONAL                       :: matrix_p,&
                                                 matrix_s
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name ='kg_gpw_fm_mol_linear', &
      routineP = module_name//':'//routine_name

    INTEGER  :: imol, imolecule_kind, ispin, istat, n_ao, n_mo(2), natom, &
                nmo_a, nmo_b, nmo_max, nmol, nmolecule_kind, nspins
    LOGICAL  :: calculate_density, failure
    REAL(dp)                                  :: my_alpha, my_beta
    REAL(dp), DIMENSION(:,:), POINTER         :: mo_coeff_a, mo_coeff_b
    REAL(dp), DIMENSION(:,:), POINTER         :: s_mo
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                 :: blocks_a, blocks_b
    TYPE(fm_mol_blocks_type), POINTER         :: mol_block_a, mol_block_b
    TYPE(kg_fm_mol_set_type), POINTER         :: fm_mol_a, fm_mol_b
    TYPE(mol_mo_set_p_type), DIMENSION(:),&
      POINTER                                 :: mos_a, mos_b
    TYPE(mol_mo_set_type), POINTER            :: mo_set_a, mo_set_b

!
!   ---------------------------------------------------------------------------
!

    failure = .FALSE. 
    calculate_density = .FALSE.
    IF(PRESENT(matrix_p) .AND. PRESENT(matrix_s)) calculate_density = .TRUE.
    nmolecule_kind=SIZE(set_a,1)
    CPPostcondition(SIZE(set_b,1)==nmolecule_kind,cp_failure_level,routineP,error,failure)

    my_alpha = 1.0_dp
    my_beta  = 1.0_dp
    IF (PRESENT(alpha)) my_alpha=alpha
    IF (PRESENT(beta)) my_beta=beta

    DO imolecule_kind = 1, nmolecule_kind
      NULLIFY(fm_mol_a,fm_mol_b, blocks_a, blocks_b)
      fm_mol_a => set_a(imolecule_kind)
      fm_mol_b => set_b(imolecule_kind)
      
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_a, nmolecule_local=nmol, natom=natom,&
           n_ao = n_ao, n_mo = n_mo, fm_mol_blocks=blocks_a)
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_b,fm_mol_blocks=blocks_b)
      CPPostcondition(SIZE(blocks_b,1)==nmol,cp_failure_level,routineP,error,failure)
      nspins = 1
      IF(n_mo(2) /= 0) nspins = 2

      IF(calculate_density) THEN
        nmo_max = MAX(n_mo(1),n_mo(2))
        ALLOCATE(s_mo(n_ao,nmo_max), STAT = istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      DO imol = 1, nmol
        NULLIFY(mol_block_a, mos_a, mol_block_b, mos_b)
        mol_block_a => blocks_a(imol)
        mol_block_b => blocks_b(imol)
        CALL get_fm_mol_block(fm_mol_block = mol_block_a,&
             mos = mos_a )
        CALL get_fm_mol_block(fm_mol_block = mol_block_b,&
             mos = mos_b)

        DO ispin=1,nspins
          NULLIFY(mo_set_a, mo_set_b)
          mo_set_a => mos_a(ispin)%mo_set
          mo_set_b => mos_b(ispin)%mo_set

          CALL get_mol_mo_set(mo_set=mo_set_a, mo = mo_coeff_a,&
               nmo=nmo_a)
          CALL get_mol_mo_set(mo_set=mo_set_b, mo = mo_coeff_b,&
               nmo=nmo_b)

          CPPostcondition(nmo_a==nmo_b,cp_failure_level,routineP,error,failure)

          IF(my_alpha == 1.0_dp) THEN
            CALL daxpy(SIZE(mo_coeff_a,1)*SIZE(mo_coeff_a,2),&
                       my_beta, mo_coeff_b(1,1),1,&
                       mo_coeff_a(1,1),1)
          ELSE
            CALL dscal(SIZE(mo_coeff_a,1)*SIZE(mo_coeff_a,2),&
                       my_alpha,mo_coeff_a(1,1),1)
            CALL daxpy(SIZE(mo_coeff_a,1)*SIZE(mo_coeff_a,2),&
                       my_beta, mo_coeff_b(1,1),1,&
                       mo_coeff_a(1,1),1)
          END IF

          IF(calculate_density) THEN

            ! Multiply Overlap matric and Coefficient matrix
            CALL  multiply_sparse_mol_mo(matrix_s(1)%matrix,mol_block_a,&
                                         natom,mo_coeff_a,nmo_a,s_mo)
            ! orthogonalize the molecular orbitals 
            CALL mol_make_basis(mo_coeff_a, nmo_a, matrix_ortho=s_mo, otype="SV")


            CALL calculate_mol_density_matrix(mo_set_a, mol_block_a, &
                 natom, matrix_p(ispin)%matrix, error=error)
          END IF
        END DO  ! ispin

      END DO  ! imol

      IF(calculate_density) THEN
        DEALLOCATE(s_mo, STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

    END DO  ! imolecule_kind


  END SUBROUTINE  kg_gpw_fm_mol_linear

!***************************************************************************
!!****f* kg_gpw_wf_history/kg_gpw_fm_mol_ps
!!
!!   NAME
!!    kg_gpw_fm_mol_ps
!!
!!   FUNCTION
!!     The ps extrapolation of the mos is calculated by using a
!!     number of previous steps given by the order of the extrapolation (nvec) 
!!     At each i-th step in the loop to compose the final mos, the very last
!!     wfn are multiplied by the wfn of the (nlast - i)th step
!!     C(next) = C(next) - beta *  C(i)[C(i)S(i)C(last)]
!!     If required the density is also calculated
!!     
!!   NOTE
!!     matrix_p and matrix_s are given as input only if the density has to be calculated
!!     i.e. when all the terms have been added to the final mos
!!
!!   INPUTS
!!     - set_a : contains the output mos that wiill be used as guess in the next step
!!     - beta : second multiplicative constant
!!     - set_b : contains the mos of the (nlast - i)th step, to be multiplied by the last mos
!!     - ov_b    : contains the overlap matrix of the (nlast - i)th
!!     - set_last : contains the mos of the step
!!     - matrix_p: density matrix in real sparse matrix format, optional
!!     - matrix_s: overlapmatrix in real sparse matrix format, optional
!!     - error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE kg_gpw_fm_mol_ps(set_a, beta, set_b, ov_b, set_last, &
                              matrix_p, matrix_s, error )

    REAL(dp), INTENT(IN), OPTIONAL            :: beta
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                 :: set_a
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                 :: set_b
    TYPE(real_matrix_type), POINTER           :: ov_b
    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                 :: set_last
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER, OPTIONAL                       :: matrix_p,&
                                                 matrix_s
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name ='kg_gpw_fm_mol_linear', &
      routineP = module_name//':'//routine_name

    INTEGER :: icol, imol, imolecule_kind,&
               ispin, istat, irow, ldsc, ldcsc, &
               n_ao, n_mo(2), natom, nmo_a, nmo_b, nmo_l, nmo_max, nmol, &
               nmolecule_kind, nspins
    LOGICAL                                   :: calculate_density, failure
    REAL(dp)                                  :: my_beta
    REAL(dp), DIMENSION(:,:), POINTER         :: mo_coeff_a, mo_coeff_b, mo_coeff_l
    REAL(dp), DIMENSION(:,:),  POINTER        :: ov_b_mol
    REAL(dp), DIMENSION(:,:), POINTER         :: mol_csc
    REAL(dp), DIMENSION(:,:), POINTER         :: sc_tmp
    REAL(dp), DIMENSION(:,:), POINTER         :: s_mo
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                 :: blocks_a, blocks_b, blocks_l
    TYPE(fm_mol_blocks_type), POINTER         :: mol_block_a, mol_block_b, mol_block_l 
    TYPE(kg_fm_mol_set_type), POINTER         :: fm_mol_a, fm_mol_b, fm_mol_l 
    TYPE(mol_mo_set_p_type), DIMENSION(:),&
      POINTER                                 :: mos_a, mos_b, mos_l
    TYPE(mol_mo_set_type), POINTER            :: mo_set_a, mo_set_b, mo_set_l
!
!   ---------------------------------------------------------------------------
!


    failure = .FALSE.
    calculate_density = .FALSE.
    IF(PRESENT(matrix_p) .AND. PRESENT(matrix_s)) calculate_density = .TRUE.
    nmolecule_kind=SIZE(set_a,1)
    CPPostcondition(SIZE(set_b,1)==nmolecule_kind,cp_failure_level,routineP,error,failure)
    CPPostcondition(SIZE(set_last,1)==nmolecule_kind,cp_failure_level,routineP,error,failure)

    my_beta  = 1.0_dp
    IF (PRESENT(beta)) my_beta=beta

    DO imolecule_kind = 1, nmolecule_kind
        NULLIFY(fm_mol_a,fm_mol_b, fm_mol_l, blocks_a, blocks_b, blocks_l)
      fm_mol_a => set_a(imolecule_kind)
      fm_mol_b => set_b(imolecule_kind)
      fm_mol_l => set_last(imolecule_kind)
      
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_a, nmolecule_local=nmol, natom=natom,&
           n_ao = n_ao, n_mo = n_mo, fm_mol_blocks=blocks_a)
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_b,fm_mol_blocks=blocks_b)
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_l,fm_mol_blocks=blocks_l)
      CPPostcondition(SIZE(blocks_b,1)==nmol,cp_failure_level,routineP,error,failure)
      CPPostcondition(SIZE(blocks_l,1)==nmol,cp_failure_level,routineP,error,failure)
      nspins = 1
      IF(n_mo(2) /= 0) nspins = 2

      ! Allocate the array to store the overlap matrix of step b
      ALLOCATE(ov_b_mol(n_ao,n_ao),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ! Allocate the array CSC 
      ldcsc = MAX(n_mo(1),n_mo(2))
      ALLOCATE(mol_csc(ldcsc,ldcsc), STAT= istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      IF(calculate_density) THEN
        nmo_max = MAX(n_mo(1),n_mo(2))
        ALLOCATE(s_mo(n_ao,nmo_max), STAT = istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      DO imol = 1, nmol
        NULLIFY(mol_block_a, mos_a, mol_block_b, mos_b, mol_block_l, mos_l)
        mol_block_a => blocks_a(imol)
        mol_block_b => blocks_b(imol)
        mol_block_l => blocks_l(imol)
        CALL get_fm_mol_block(fm_mol_block = mol_block_a,&
             mos = mos_a, work = sc_tmp )
        CPPrecondition(ASSOCIATED(sc_tmp),cp_failure_level,routineP,error,failure)
        ldsc = SIZE(sc_tmp,1)
        CPPrecondition(ldsc>=n_ao,cp_failure_level,routineP,error,failure)
        CALL get_fm_mol_block(fm_mol_block = mol_block_b,&
             mos = mos_b)
        CALL get_fm_mol_block(fm_mol_block = mol_block_l,&
             mos = mos_l)

        ! Copy the right block of ov_b into the temporary space  
        CALL copy_sparse2mol_block(matrix=ov_b, mol_block=mol_block_b,&
             matrix_mol_block=ov_b_mol, natom=natom, nrow=n_ao, ncol=n_ao, error=error)


        DO ispin=1,nspins
          NULLIFY(mo_set_a, mo_set_b, mo_set_l)
          NULLIFY(mo_coeff_a, mo_coeff_b, mo_coeff_l)

          mo_set_a => mos_a(ispin)%mo_set
          mo_set_b => mos_b(ispin)%mo_set
          mo_set_l => mos_l(ispin)%mo_set

          CALL get_mol_mo_set(mo_set=mo_set_a, mo = mo_coeff_a,&
               nmo=nmo_a)
          CALL get_mol_mo_set(mo_set=mo_set_b, mo = mo_coeff_b,&
               nmo=nmo_b)
          CALL get_mol_mo_set(mo_set=mo_set_l, mo = mo_coeff_l,&
               nmo=nmo_l)

          CPPostcondition(nmo_a==nmo_b,cp_failure_level,routineP,error,failure)
          CPPostcondition(nmo_a==nmo_l,cp_failure_level,routineP,error,failure)

          CALL DGEMM('N','N',n_ao,nmo_l,n_ao,1.0_dp,ov_b_mol(1,1),n_ao,&
               mo_coeff_l(1,1),SIZE(mo_coeff_l,1),0.0_dp,sc_tmp(1,1),ldsc)
          CALL DGEMM('T','N',nmo_b,nmo_b,n_ao,1.0_dp,mo_coeff_b(1,1),n_ao,&
               sc_tmp(1,1),ldsc,0.0_dp,mol_csc(1,1),ldcsc)
          CALL DGEMM('N','N',n_ao,nmo_b,nmo_b,1.0_dp,mo_coeff_b(1,1),n_ao,&
               mol_csc(1,1),ldcsc,0.0_dp,sc_tmp(1,1),ldsc)
          DO icol = 1,nmo_a
            DO irow = 1,n_ao
              mo_coeff_a(irow,icol) = mo_coeff_a(irow,icol) + beta*sc_tmp(irow,icol)
            END DO
          END DO

          ! Calculate the density matrix with the completed new guess for the mos
          IF(calculate_density) THEN
            ! Multiply Overlap matric and Coefficient matrix
            CALL  multiply_sparse_mol_mo(matrix_s(1)%matrix,mol_block_a,&
                                         natom,mo_coeff_a,nmo_a,s_mo)
            ! orthogonalize the molecular orbitals 
            CALL mol_make_basis(mo_coeff_a, nmo_a, matrix_ortho=s_mo, otype="SV")

            CALL calculate_mol_density_matrix(mo_set_a, mol_block_a, &
                 natom, matrix_p(ispin)%matrix, error=error)
          END IF

        END DO  ! ispin

      END DO ! imol
      IF(calculate_density) THEN
        DEALLOCATE(s_mo, STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      DEALLOCATE(mol_csc,STAT = istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      DEALLOCATE(ov_b_mol,STAT = istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END DO  ! imolecule_kind


  END SUBROUTINE  kg_gpw_fm_mol_ps

!***************************************************************************
!!****f* kg_gpw_wf_history/kg_gpw_prev_wf
!!
!!   NAME
!!    kg_gpw_prev_wf 
!!
!!   FUNCTION
!!     Calculate the density from the previous wfn given in molecular blocks 
!!     
!!   NOTE
!!
!!   INPUTS
!!     - kg_fm_mol : containes the blocks of mos
!!     - matrix_p: density matrix in real sparse matrix format
!!     - matrix_s: overlapmatrix in real sparse matrix format
!!     - error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE kg_gpw_prev_wf(kg_fm_mol_set,matrix_p,matrix_s,error)

    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                    :: kg_fm_mol_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                    :: matrix_p,&
                                                    matrix_s
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error

   CHARACTER(LEN=*), PARAMETER :: routine_name = "kg_gpw_prev_wf", &
      routineP = module_name//":"//routine_name

    INTEGER :: imol, imolecule_kind, ispin, istat, n_ao, nat_mol, natom, nmo, nmo_kind(2),&
               nmo_max, nmol_local, nmolecule_kind, nspins 
    LOGICAL :: failure
    REAL(dp), DIMENSION(:,:), POINTER            :: s_mo
    REAL(dp), DIMENSION(:,:), POINTER            :: mo_coeff
    TYPE(fm_mol_blocks_type), POINTER            :: mol_block
    TYPE(fm_mol_blocks_type), DIMENSION(:),&
      POINTER                                    :: fm_mol_blocks
    TYPE(kg_fm_mol_set_type), POINTER            :: fm_mol_set
    TYPE(mol_mo_set_p_type), DIMENSION(:),&
      POINTER                                    :: mos
    TYPE(mol_mo_set_type), POINTER               :: mo_set

!
!   ---------------------------------------------------------------------------
!
    failure = .FALSE.
    nmolecule_kind = SIZE(kg_fm_mol_set)

      DO imolecule_kind = 1,nmolecule_kind
        NULLIFY(fm_mol_set,fm_mol_blocks)
        fm_mol_set => kg_fm_mol_set(imolecule_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set, &
                              natom = nat_mol, &
                              n_ao = n_ao , n_mo = nmo_kind, &
                              fm_mol_blocks = fm_mol_blocks,&
                              nmolecule_local = nmol_local)
        nspins = 1
        IF(nmo_kind(2)/=0)  nspins = 2

        nmo_max = MAX(nmo_kind(1),nmo_kind(2))
        ALLOCATE(s_mo(n_ao,nmo_max), STAT = istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        DO imol = 1, nmol_local
          NULLIFY(mol_block,mos)
          mol_block => fm_mol_blocks(imol)
          CALL get_fm_mol_block(fm_mol_block = mol_block,&
                                mos = mos)
          DO ispin=1,nspins
            NULLIFY(mo_set,mo_coeff)
            mo_set => mos(ispin)%mo_set
            CALL get_mol_mo_set(mo_set=mo_set, mo = mo_coeff,&
                                nmo=nmo)

            ! Multiply Overlap matric and Coefficient matrix
            CALL  multiply_sparse_mol_mo(matrix_s(1)%matrix,mol_block,&
                                         nat_mol,mo_coeff,nmo,s_mo)
            ! orthogonalize the molecular orbitals 
            CALL mol_make_basis(mo_coeff, nmo, matrix_ortho=s_mo, otype="SV")

            CALL calculate_mol_density_matrix(mo_set,mol_block,&
                                       nat_mol,matrix_p(ispin)%matrix,error=error)
          ENDDO  ! ispin
        END DO  ! imol

        DEALLOCATE(s_mo,STAT = istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      END DO  ! imolecule_kind

  END SUBROUTINE  kg_gpw_prev_wf

!!*** **********************************************************************

END MODULE kg_gpw_wf_history
