!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****** cp2k/kg_gpw_fm_mol_methods [1.0] *
!!
!!   NAME
!!     kg_gpw_fm_mol_methods
!!
!!   FUNCTION
!!     Some utilities to handle the molecular blocks of the full matrix
!!     of the molecular orbitals
!!     
!!   AUTHOR
!!     MI (29.11.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_fm_mol_utils

  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_val_get
  USE kg_gpw_fm_mol_types,             ONLY: fm_mol_blocks_type,&
                                             get_fm_mol_block,&
                                             get_kg_fm_mol_set,&
                                             get_mol_mo_set,&
                                             kg_fm_mol_set_type,&
                                             kg_fm_p_type,&
                                             mol_mo_set_p_type,&
                                             mol_mo_set_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size, default_path_length
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_type,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_gpw_fm_mol_utils'

! *** Public subroutines ***

  PUBLIC :: calculate_w_matrix_per_molecule, &
            copy_sparse2mol_block, fm_mol_cholesky_restore, &
            fm_mol_power, fm_mol_syevd, fm_mol_syevx, mol_mo_random,&
            read_mo_mol_restart, write_mo_mol_restart

! *****************************************************************************

CONTAINS
!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/calculate_w_matrix_per_molecule
!!
!!   NAME
!!     calculate_w_matrix_per_molecule
!!
!!   FUNCTION
!!     Create the sparse matrix W =  fE*MO by MO^T, where fE is the
!!     diagonal matrix of the eigenvalues multiplied by the occupation numbers
!!     Here each molecular block of the MO is treated independently and the 
!!     result is copied inthe right part of the sparse W matrix.
!!     Therefore everything is done in a loop over the molecules
!!
!!   INPUTS
!!     - fm_mol_set: molecular MO blocks and other info about the molecules
!!     - matrix_w: sparse matrix, the structure must be already there
!!     - ispin
!!     - error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE calculate_w_matrix_per_molecule(fm_mol_set,w_matrix,ispin,error)

    TYPE(kg_fm_mol_set_type), DIMENSION(:), &
      POINTER                                :: fm_mol_set
    TYPE(real_matrix_type), POINTER          :: w_matrix
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'calculate_w_matrix_per_molecule', &
      routineP = moduleN//':'//routineN

    INTEGER :: homo, iao, iat_mol, iatom, ifs, ils, imo, imol, &
      imolecule_kind, inao, istat, jao, jat_mol, jatom, jfs, jls, jnao, ldao, &
      ldmo, nao_mol, nat_mol, nmo_mol(2), nmol_local, nmolecule_kind
    INTEGER, DIMENSION(:), POINTER           :: i_atom, ifirst_ao, ilast_ao
    LOGICAL                                  :: failure
    REAL(dp)                                 :: eigocc
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, occupation
    REAL(dp), DIMENSION(:, :), POINTER       :: mo, w_block, wmo, wmo_x_mot
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos_mol
    TYPE(mol_mo_set_type), POINTER           :: mo_set

!   ---------------------------------------------------------------------------

    failure = .FALSE.
    nmolecule_kind = SIZE(fm_mol_set,1)
    NULLIFY(wmo_x_mot,mo,wmo)
    ldmo = 1
    ldao = 1
    ALLOCATE(wmo(ldao,ldmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    wmo = 0.0_dp
    CALL set_matrix(w_matrix,0.0_dp)

    DO imolecule_kind = 1, nmolecule_kind

      NULLIFY (fm_mol, fm_mol_blocks)
      fm_mol => fm_mol_set(imolecule_kind)
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol,&
                             nmolecule_local=nmol_local,&
                             natom=nat_mol, n_ao=nao_mol, n_mo=nmo_mol,&
                             fm_mol_blocks=fm_mol_blocks)

      IF((ldao < nao_mol) .OR. (ldmo < nmo_mol(ispin))) THEN
        ldao = nao_mol
        ldmo = nmo_mol(ispin)
        CALL reallocate(wmo,1,ldao,1,ldmo)
      END IF

      DO imol = 1, nmol_local
        NULLIFY(i_atom,ifirst_ao,ilast_ao,mos_mol, mol_block)
        mol_block => fm_mol_blocks(imol)
        CALL get_fm_mol_block(fm_mol_block = mol_block,&
                              mos=mos_mol,&
                              work = wmo_x_mot, &
                              index_atom = i_atom,&
                              ifirst_ao = ifirst_ao,&
                              ilast_ao = ilast_ao )

        NULLIFY(eigenvalues, mo, mo_set, occupation)
        mo_set => mos_mol(ispin)%mo_set
        CALL get_mol_mo_set(mo_set, homo = homo, mo = mo ,&
                            eigenvalues = eigenvalues,&
                            occupation_numbers = occupation ) 

        DO imo = 1,homo
          eigocc = eigenvalues(imo)*occupation(imo)
          CALL DCOPY(nao_mol,mo(1,imo),1,wmo(1,imo),1)
          CALL DSCAL(nao_mol,eigocc,wmo(1,imo),1)
        END DO 
        CALL DGEMM('N','T',nao_mol,nao_mol,homo,1.0_dp,wmo(1,1),ldao,&
                   mo(1,1),ldao,0.0_dp,wmo_x_mot,SIZE(wmo_x_mot,1))

        DO iat_mol = 1,nat_mol
          iatom = i_atom(iat_mol)
          ifs = ifirst_ao(iat_mol)
          ils = ilast_ao(iat_mol)
          inao = ils - ifs + 1

          DO jat_mol = iat_mol,nat_mol
            jatom = i_atom(jat_mol)
            jfs = ifirst_ao(jat_mol)
            jls = ilast_ao(jat_mol)
            jnao = jls - jfs + 1
 

            NULLIFY(w_block)
            IF(iatom <= jatom) THEN

              CALL get_block_node(matrix=w_matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=w_block)
              CPPrecondition(ASSOCIATED(w_block),cp_failure_level,routineP,error,failure)
              CPPrecondition(inao==SIZE(w_block,1),cp_failure_level,routineP,error,failure)
              CPPrecondition(jnao==SIZE(w_block,2),cp_failure_level,routineP,error,failure)

              DO jao = 1,jnao
                DO iao = 1,inao
                  w_block(iao,jao) = wmo_x_mot(ifs+iao-1,jfs+jao-1)
                END DO  ! iao
              END DO  ! jao

            ELSE

              CALL get_block_node(matrix=w_matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=w_block)
              CPPrecondition(ASSOCIATED(w_block),cp_failure_level,routineP,error,failure)
              CPPrecondition(jnao==SIZE(w_block,1),cp_failure_level,routineP,error,failure)
              CPPrecondition(inao==SIZE(w_block,2),cp_failure_level,routineP,error,failure)

              DO iao = 1,inao
                DO jao = 1,jnao
                  w_block(jao,iao) = wmo_x_mot(jfs+jao-1,ifs+iao-1)
                END DO  ! jao
              END DO  ! iao

            END IF

          END DO  ! jat_mol
        END DO  ! iat_mol
      END DO  ! imol

    END DO  ! imolecule_kind

    DEALLOCATE(wmo,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE calculate_w_matrix_per_molecule

!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/fm_mol_cholesky_restore
!!
!!   NAME
!!     fm_mol_cholesky_restore
!!
!!   FUNCTION
!!     Multiply one block (molecular block for KG_GPW) of the MO full matrix by the 
!!     corresponding cholesky decomposition of the overlap matrix, U
!!
!!   INPUTS
!!     -vectors_in: set of MO as input
!!     -nmo:  number of MO we need
!!     -upper_dec: cholesky decomposition of the overlap matrix, U
!!     -vectors_out:  set of MO as output after restore
!!     - op : can be "SOLVE" (out = U^-1 * in ) or "MULTIPLY"   (out = U * in )
!!     - pos: can be "LEFT" or "RIGHT" (U at the left or at the right)
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE fm_mol_cholesky_restore(vectors_in,nmo,upper_dec,vectors_out,&
                                     op,pos,error)

    REAL(dp), DIMENSION(:, :)                :: vectors_in
    INTEGER, INTENT(IN)                      :: nmo
    REAL(dp), DIMENSION(:, :), POINTER       :: upper_dec, vectors_out
    CHARACTER(LEN=*), INTENT(IN)             :: op
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: pos
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_mol_cholesky_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: my_pos
    INTEGER                                  :: nao
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha

!  --------------------------------------------------------------------------------

    failure = .FALSE.

    CALL cp_assert(op == "SOLVE" .OR. op == "MULTIPLY",&
         cp_failure_level,cp_assertion_failed,routineP,&
         "wrong argument op", error,failure)

    IF (PRESENT(pos)) THEN
       SELECT CASE(pos)
       CASE("LEFT")
         my_pos='L'
       CASE("RIGHT")
         my_pos='R'
       CASE DEFAULT
          CALL cp_assert(.FALSE.,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "wrong argument pos", error,failure)
       END SELECT
    ELSE
       my_pos='L'
    ENDIF

    nao = SIZE(upper_dec,1)
    CPPostcondition(SIZE(vectors_out,1)==nao,cp_failure_level,routineP,error,failure)
    CPPostcondition(SIZE(vectors_out,2)>=nmo,cp_failure_level,routineP,error,failure)
    alpha = 1.0_dp
    CALL dcopy(nmo*nao,vectors_in(1,1),1,vectors_out(1,1),1)
    IF (op.EQ."SOLVE") THEN
      CALL dtrsm(my_pos,'U','N','N',nao,nmo,alpha,upper_dec(1,1),nao,&
                 vectors_out(1,1),nao)
    ELSE
      CALL dtrmm(my_pos,'U','N','N',nao,nmo,alpha,upper_dec(1,1),nao,&
                 vectors_out(1,1),nao)
    ENDIF

  END SUBROUTINE fm_mol_cholesky_restore

!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/fm_mol_power
!!
!!   NAME
!!     fm_mol_power
!!
!!   FUNCTION
!!     Raise the real symmetric n by n matrix to the power given by
!!     exponent. All eigenvectors with a corresponding eigenvalue lower
!!     than threshold are quenched. Here the n by n matrix is the block of
!!     the overlap matrix corresponding to one molecule 
!!
!!   INPUTS
!!     -matrix: block of the overlap matrix
!!     -eigenvectors: 
!!     -ldv: leading dimension of the matrix
!!     -ncv: number of vectors
!!     -exponent
!!     -threshold: to quench the eigenvalues
!!     -n_dependent
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE fm_mol_power(matrix,eigenvectors,ldv,ncv,exponent,threshold,n_dependent)

    REAL(dp), DIMENSION(:, :), POINTER       :: matrix
    REAL(dp), DIMENSION(:, :)                :: eigenvectors
    INTEGER, INTENT(IN)                      :: ldv, ncv
    REAL(dp), INTENT(IN)                     :: exponent, threshold
    INTEGER, INTENT(OUT)                     :: n_dependent

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_mol_power', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icol, istat, ncol, &
                                                nrow
    REAL(dp)                                 :: f, p
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eigenvalues

!   --------------------------------------------------------------------------

    CALL timeset("fm_mol_power","I","",handle)

    nrow = SIZE(matrix,1)
    ncol = SIZE(matrix,2)

    n_dependent = 0

    ALLOCATE(eigenvalues(ncol),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"eigenvalues",ncol*dp_size)
    eigenvalues(:) = 0.0_dp

    CALL fm_mol_syevx(matrix,eigenvectors,ldv,ncv,eigenvalues,ncol)

    p = 0.5_dp*exponent

    DO icol = 1,ncol
       IF(eigenvalues(icol) < threshold) THEN 
         n_dependent = n_dependent + 1
         eigenvectors(1:nrow,icol) = 0.0_dp
       ELSE

         f = eigenvalues(icol)**p
         eigenvectors(1:nrow,icol) =&
             f*eigenvectors(1:nrow,icol)
       END IF
    END DO

    CALL dsyrk("U","N",nrow,ncol,1.0_dp,eigenvectors(1,1),ldv,&
                                 0.0_dp,matrix(1,1),nrow)

    DEALLOCATE(eigenvalues, STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineP,"eigenvalues")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE fm_mol_power

!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/fm_mol_syevd
!!
!!   NAME
!!     fm_mol_syevd
!!
!!   FUNCTION
!!     computes all eigenvalues and vectors of a real symmetric matrix    
!!
!!   INPUTS
!!     -matrix: block of the full matrix
!!     -nrow: leading dimension of the matrix
!!     -eigenvalues
!!     -matrix_out
!!     -error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************


  SUBROUTINE fm_mol_syevd(matrix,nrow,eigenvalues,matrix_out,error)

    REAL(dp), DIMENSION(:, :), POINTER       :: matrix
    INTEGER, INTENT(IN)                      :: nrow
    REAL(dp), DIMENSION(:)                   :: eigenvalues
    REAL(dp), DIMENSION(:, :)                :: matrix_out
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_mol_syevd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info, istat, lda, liwork, &
                                                lwork
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work

!   --------------------------------------------------------------------------

    failure = .FALSE.

    IF(ASSOCIATED(matrix)) THEN
      lda = SIZE(matrix,1)
      CPPostcondition(nrow<=lda,cp_failure_level,routineP,error,failure)

      lwork=1+6*nrow+2*nrow**2+5000
      liwork=5*nrow+3
      ALLOCATE(work(lwork),STAT=istat)
      IF (istat.NE.0) CALL stop_memory("eigensolver","work",lwork*dp_size)
      ALLOCATE(iwork(liwork),STAT=istat)
      IF (istat.NE.0) CALL stop_memory("eigensolver","iwork",liwork*int_size)
      CALL DSYEVD('V','U', nrow, matrix(1,1), nrow, eigenvalues(1),&
                  work(1), lwork, iwork(1), liwork, info)

      !   *** Error handling ***
      IF (info /= 0) CALL stop_program(routineP,"Error in dsyevd")

      CALL DCOPY (nrow*nrow,matrix(1,1),1,matrix_out(1,1),1)

    ELSE
      !stop error
      CALL stop_program(routineP,"Error in dsyevd")
    END IF

  END SUBROUTINE fm_mol_syevd

!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/fm_mol_syevx
!!
!!   NAME
!!     fm_mol_syevx
!!
!!   FUNCTION
!!     Diagonalise the symmetric n by n matrix using the LAPACK library 
!!     Compute ncol egenvalues and eigenvectors
!!     
!!   INPUTS
!!     -matrix: block of the full matrix
!!     -eigenvetors
!!     -ldv: leading dimension of the matrix
!!     -ncv: number of columns available
!!     -eigenvalues
!!     -ncol: number of vector to be calculated
!!     -error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE fm_mol_syevx(matrix,eigenvetors,ldv,ncv,eigenvalues,ncol,error)

    REAL(dp), DIMENSION(:, :), POINTER       :: matrix
    REAL(dp), DIMENSION(:, :)                :: eigenvetors
    INTEGER, INTENT(IN)                      :: ldv, ncv
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: eigenvalues
    INTEGER, INTENT(IN)                      :: ncol
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_mol_syevx', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: vl = 0.0_dp, vu = 0.0_dp 

    INTEGER                                  :: info, istat, ldb, liwork, &
                                                lwork, m, nb, nrow
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ifail, iwork
    INTEGER, EXTERNAL                        :: ilaenv
    LOGICAL                                  :: failure
    REAL(dp)                                 :: abstol
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: w, work
    REAL(KIND=dp), EXTERNAL                  :: dlamch

!   --------------------------------------------------------------------------

    failure =  .FALSE.
    abstol = 2.0_dp*dlamch("S")
    nrow = SIZE(matrix,1)
    ldb = ldv  ! SIZE(eigenvetors,1)
!    CPPostcondition(ncol<=SIZE(eigenvetors,2),cp_failure_level,routineP,error,failure)
    CPPostcondition(ncol<=ncv,cp_failure_level,routineP,error,failure)
    CPPostcondition(ncol<=SIZE(eigenvalues,1),cp_failure_level,routineP,error,failure)

    ALLOCATE(w(nrow), STAT =istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w",nrow*dp_size)
    w(:) = 0.0_dp

    eigenvalues(:) = 0.0_dp

    nb = MAX(ilaenv(1,"DSYTRD","U",nrow,-1,-1,-1),&
             ilaenv(1,"DORMTR","U",nrow,-1,-1,-1))

    lwork = MAX((nb + 3)*nrow,8*nrow)+nrow ! sun bug fix
    liwork = 5*nrow

    ALLOCATE (ifail(nrow),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"ifail",nrow*int_size)
    ALLOCATE (iwork(liwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"iwork",liwork*int_size)
    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"work",lwork*dp_size)

    CALL dsyevx("V","I","U",nrow,matrix(1,1),nrow,vl,vu,1,ncol,abstol,&
                m,w,eigenvetors(1,1),ldb,work(1),lwork,&
                iwork(1),ifail(1),info)

!   *** Error handling ***
    IF (info /= 0) CALL stop_program(routineP,"Error in dsyevx")

    eigenvalues(1:ncol) = w(1:ncol)

!   *** Release work storage ***
    DEALLOCATE(w,ifail,iwork,work,STAT =istat)
    IF (istat /= 0) CALL stop_memory(routineP,"w,ifail,iwork,work")

  END SUBROUTINE fm_mol_syevx 

!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/mol_mo_random
!!
!!   NAME
!!     mol_mo_random
!!
!!   FUNCTION
!!     Assign a random value to the MO block matrix
!!     
!!   INPUTS
!!     -mo: block of the full matrix
!!     -mepos
!!     -num_pe
!!     -ncol: number of vector to be calculated
!!     -start_col: from which column to start
!!     -error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE mol_mo_random(mo,mepos,num_pe,ncol,start_col,error)

    REAL(dp), DIMENSION(:, :), POINTER       :: mo
    INTEGER, INTENT(IN)                      :: mepos, num_pe
    INTEGER, INTENT(IN), OPTIONAL            :: ncol, start_col
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER              :: routineN = 'mol_mo_random'

    INTEGER                                  :: icol, my_ncol, my_start, nrow
    INTEGER, DIMENSION(4), SAVE              :: ISEED
    LOGICAL                                  :: failure
    LOGICAL, SAVE                            :: FIRST = .TRUE.

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    IF (FIRST) THEN
      ISEED(1)= mepos
      ISEED(2)= mepos + num_pe
      ISEED(3)= mepos + 2*num_pe
      ISEED(4)= mepos + 3*num_pe ! last one has to be odd
      IF(MOD(ISEED(4),2) .EQ. 0) ISEED(4)=ISEED(4)+1
      FIRST=.FALSE.
    ENDIF

    CPPostcondition(ASSOCIATED(mo),cp_failure_level,routineN,error,failure)
    my_start = 1
    IF(PRESENT(start_col)) my_start = start_col
    my_ncol = SIZE(mo,2)
    IF(PRESENT(ncol)) THEN
      CPPostcondition(my_ncol>=ncol,cp_failure_level,routineN,error,failure)
      my_ncol = ncol
    END IF

    nrow = SIZE(mo,1)
    DO icol = my_start ,my_start + my_ncol - 1
      CALL dlarnv(1,iseed,nrow,mo(1,icol))
    END DO 

  END SUBROUTINE mol_mo_random

!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/copy_sparse2mol_block
!!
!!   NAME
!!     copy_sparse2mol_block
!!
!!   FUNCTION
!!     Copy a block of a sparse matrix into the right part of the
!!     molecular block of a full matrix (full for one molecule not the whole system)
!!     
!!   INPUTS
!!     -matrix
!!     -mol_block
!!     -matrix_mol_block
!!     -natom
!!     -nrow
!!     -ncol
!!     -error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE copy_sparse2mol_block(matrix, mol_block, matrix_mol_block, &
                                   natom, nrow, ncol, error)

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    REAL(dp), DIMENSION(:, :), POINTER       :: matrix_mol_block
    INTEGER, INTENT(IN)                      :: natom, nrow, ncol
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_sparse2mol_block', &
      routineP = moduleN//':'//routineN

    INTEGER :: col_first, col_last, iat, iatom, icol, icol_mol_s, irow, &
      irow_mol_s, jat, jatom, row_first, row_last
    INTEGER, DIMENSION(:), POINTER           :: i_atom, ifirst_ao, ilast_ao
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: s_block

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    NULLIFY(i_atom,ifirst_ao,ilast_ao)

    CALL get_fm_mol_block(fm_mol_block = mol_block,&
                          index_atom = i_atom,&
                          ifirst_ao = ifirst_ao,&
                          ilast_ao = ilast_ao )

    CPPrecondition(ASSOCIATED(matrix_mol_block),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(matrix_mol_block,1)>=nrow,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(matrix_mol_block,2)>=ncol,cp_failure_level,routineP,error,failure)

 ! Transfer the sparse matrix in the molecular block
    DO iat = 1,natom
      iatom = i_atom(iat)
      row_first = ifirst_ao(iat)
      row_last  = ilast_ao(iat)
      DO jat = 1 , natom
        jatom = i_atom(jat)
        col_first = ifirst_ao(jat)
        col_last  = ilast_ao(jat)

        NULLIFY(s_block)
        IF (iatom <= jatom) THEN
          CALL get_block_node(matrix=matrix,&
                              block_row=iatom,&
                              block_col=jatom,&
                              BLOCK=s_block)
          icol =  1 
          DO icol_mol_s = col_first, col_last
            irow = 1
            DO irow_mol_s = row_first, row_last
              matrix_mol_block(irow_mol_s,icol_mol_s) = s_block(irow,icol)
              irow = irow + 1
            END DO  ! irow_mol_s 
            icol = icol + 1
          END DO   ! icol_mol_s
        ELSE
          CALL get_block_node(matrix=matrix,&
                              block_row=jatom,&
                              block_col=iatom,&
                              BLOCK=s_block)
          icol =  1 
          DO icol_mol_s = col_first, col_last
            irow = 1
            DO irow_mol_s = row_first, row_last
              matrix_mol_block(irow_mol_s,icol_mol_s) = s_block(icol,irow)
              irow = irow + 1
            END DO  ! irow_mol_s 
            icol = icol + 1
          END DO  ! icol_mol_s
        END IF      
      END DO  ! jat
    END DO  ! iat

  END SUBROUTINE copy_sparse2mol_block


!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/write_mo_mol_restart
!!
!!   NAME
!!     write_mo_mol_restart
!!
!!   FUNCTION
!!     Write a RESTART file for the KG_GPW calculation
!!     The MO are written in a list of blocks, each one corresponding
!!     to one molecule
!!     Since the molecule are distributed over the projectors, some commuication
!!     is required to write the list in the correct order (order of the molecule)
!!     
!!   NOTE
!!     The KG_GPW RESTART file has a special form. It should be read and written only
!!     in a KG_GPW calculation. However the defoult name is still RESTART.
!!     An initial flag identifies this as a KG_GPW RESTART.
!!
!!   INPUTS
!!     - kg_fm : containes all the info abut the molecules and the MOof each molecule
!!     - globenv: parallel distribution and i/o information
!!     - id_nr: give an index to the restart if =/ 1
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE write_mo_mol_restart(kg_fm,globenv,nspins,id_nr, error)

    TYPE(kg_fm_p_type)                       :: kg_fm
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: nspins, id_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_mol_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)                       :: file_name
    CHARACTER(LEN=6)                         :: method
    INTEGER :: group, handle, imo, imol_found, imol_global, imol_local, &
      imol_proc, imolecule_kind, ispin, istat, n_ao, n_mo(2), nao_max, &
      nelectron_global, nelectron_spin(2), nmo_max, nmolecule_global, &
      nmolecule_kind, nmolecule_local, restart_unit, source
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nao_molkind, natom_molkind, &
                                                nmol_molkind
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: nel_molkind
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: eig_buffer, occ_buffer
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: mo_buffer
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, &
                                                occupation_numbers
    REAL(dp), DIMENSION(:, :), POINTER       :: mo
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                :: mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(mol_mo_set_type), POINTER           :: mo_set

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.

    group = globenv%group
    source = globenv%source

    method ="KG_GPW"
    nmolecule_kind = SIZE(kg_fm%kg_fm_mol_set)
    nmolecule_global = kg_fm%nmolecule_global
    nelectron_global = kg_fm%nelectron_global
    nao_max = kg_fm%nao_max
    nmo_max = kg_fm%nmo_max

    IF (globenv%ionode) THEN

       CALL section_vals_val_get(globenv%input_file,"FORCE_EVAL%DFT%RESTART_FILE_NAME",&
            c_val=file_name,error=error)
      IF (id_nr/=1) THEN
        file_name(LEN_TRIM(file_name):LEN(file_name)) = "-"//&
                    ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="WRITE",&
                     file_form="UNFORMATTED",&
                     file_status="REPLACE",&
                     unit_number=restart_unit)

!     *** create some info about the molecules ***


      ALLOCATE(nmol_molkind(nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(nao_molkind(nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(natom_molkind(nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(nel_molkind(2,nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      DO imolecule_kind = 1,nmolecule_kind
        NULLIFY(fm_mol_set)
        fm_mol_set => kg_fm%kg_fm_mol_set(imolecule_kind)
        CALL  get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,&
                                nmol_kind_global=nmol_molkind(imolecule_kind),&
                                n_ao = nao_molkind(imolecule_kind),&
                                natom = natom_molkind(imolecule_kind), &
                                nelectron_spin=nel_molkind(1:2,imolecule_kind) )
      END DO 

      WRITE (restart_unit) method, nmolecule_kind, nmolecule_global, nelectron_global, nspins
      WRITE (restart_unit) nmol_molkind
      WRITE (restart_unit) nao_molkind
      WRITE (restart_unit) natom_molkind
      WRITE (restart_unit) nel_molkind(1,1:nmolecule_kind)
      WRITE (restart_unit) nel_molkind(2,1:nmolecule_kind)

      DEALLOCATE(nmol_molkind,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(nao_molkind,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(natom_molkind,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(nel_molkind,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF  ! ionode

    ! The following is done by all the processor
    ! The MOs etc are communicated to all the processor, one molecule at a time 

    ALLOCATE(mo_buffer(nao_max,nmo_max,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(eig_buffer(nmo_max,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(occ_buffer(nmo_max,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO imol_global = 1,nmolecule_global
      imol_found = 0
      imol_proc = -1
      CALL dcopy(nao_max*nmo_max*nspins,0.0_dp,0,mo_buffer(1,1,1),1)
      CALL dcopy(nmo_max*nspins,0.0_dp,0,eig_buffer(1,1),1)
      CALL dcopy(nmo_max*nspins,0.0_dp,0,occ_buffer(1,1),1)
      DO imolecule_kind = 1,nmolecule_kind
        NULLIFY(fm_mol_set,mol_blocks,mos, mo_set)
        NULLIFY(mo,eigenvalues,occupation_numbers)
        fm_mol_set => kg_fm%kg_fm_mol_set(imolecule_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,  &
                               nmolecule_local = nmolecule_local,&
                               n_ao=n_ao, n_mo=n_mo,nelectron_spin=nelectron_spin,&
                               fm_mol_blocks=mol_blocks) 
        DO imol_local = 1,nmolecule_local
           mol_block => mol_blocks(imol_local)
           IF(imol_global == mol_block%imol_global) THEN
             ! Molecule found, now retrieve the MO etc
             ! fill mo_buffer, eig_buffer, occ_buffer one spin at a time
             ! write them down in the restart
             ! and ciaociao
             imol_found = 1
             imol_proc = globenv%mepos
             CALL get_fm_mol_block(fm_mol_block=mol_block, mos = mos)
             DO ispin = 1,nspins
                mo_set => mos(ispin)%mo_set
                CALL get_mol_mo_set(mo_set=mo_set,mo=mo,eigenvalues=eigenvalues,&
                                    occupation_numbers=occupation_numbers)

!dbg
!  WRITE(*,*) imol_local, mo(1:n_ao,1)
!dbg
                DO imo = 1 , n_mo(ispin)
                  CALL DCOPY(n_ao,mo(1,imo),1,mo_buffer(1,imo,ispin),1)
                  eig_buffer(imo,ispin) = eigenvalues(imo)
                  occ_buffer(imo,ispin) = occupation_numbers(imo)
                END DO
             END DO
           END IF
                      
        END DO  ! imol_local

      END DO  ! imolecule_kind

      CALL mp_sum(imol_found,group)

      IF(imol_found==0 ) THEN
        IF(globenv%ionode) &
           WRITE(globenv%scr,*) " WRITE RESTART : ", imol_global, " MOLECULE NOT FOUND ON ANY PROC. "
        CALL stop_program (routineP,"Error ")
      ELSEIF(imol_found >=2 ) THEN
        IF(globenv%ionode) &
             WRITE(globenv%scr,*) " WRITE RESTART : ", imol_global, " ON MORE THAN ONE PROC. "
        CALL stop_program (routineP,"Error ")
      ELSE 

        CALL mp_sum(mo_buffer,group)
        CALL mp_sum(eig_buffer,group)
        CALL mp_sum(occ_buffer,group)

        DO ispin = 1,nspins
         IF(globenv%ionode.AND.(n_mo(ispin)>0)) THEN
           WRITE (restart_unit) imol_global , ispin
           WRITE (restart_unit) n_mo(ispin),n_ao,&
                                nelectron_spin(ispin)
           WRITE (restart_unit) eig_buffer(1:n_mo(ispin),ispin)
           WRITE (restart_unit) occ_buffer(1:n_mo(ispin),ispin)
           DO imo = 1,n_mo(ispin)
             WRITE (restart_unit) mo_buffer(1:n_ao,imo,ispin)
          END DO

         END IF
        END DO

      END IF

    END DO ! imol_global
 
    DEALLOCATE(mo_buffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(eig_buffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(occ_buffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (globenv%ionode) CALL close_file(unit_number=restart_unit)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE write_mo_mol_restart

!***************************************************************************
!!****f* kg_gpw_fm_mol_utils/read_mo_mol_restart
!!
!!   NAME
!!     read_mo_mol_restart
!!
!!   FUNCTION
!!     Read  the RESTART file for the KG_GPW calculation
!!     The MO are written in a list of blocks, each one corresponding
!!     to one molecule
!!     Since the molecule are distributed over the projectors, some commuication
!!     is required to write the list in the correct order (order of the molecule)
!!     
!!   NOTE
!!     The KG_GPW RESTART file has a special form. It should be read and written only
!!     in a KG_GPW calculation. However the defoult name is still RESTART.
!!     An initial flag identifies this as a KG_GPW RESTART.
!!
!!   INPUTS
!!     - kg_fm : containes all the info abut the molecules and the MOof each molecule
!!     - globenv: parallel distribution and i/o information
!!     - nspins
!!     - id_nr: give an index to the restart if =/ 1
!!     - error
!!
!!   AURHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************

  SUBROUTINE read_mo_mol_restart(kg_fm,globenv,nspins,id_nr, error)

    TYPE(kg_fm_p_type)                       :: kg_fm
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: nspins, id_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mo_mol_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=200)                       :: file_name
    CHARACTER(LEN=6)                         :: method
    INTEGER :: group, handle, imo, imol_found, imol_global, imol_global_read, &
      imol_local, imol_proc, imolecule_kind, ispin, ispin_read, istat, n_ao, &
      n_ao_read, n_mo(2), n_mo_read(2), nao_max, natom, nelectron_global, &
      nelectron_global_read, nelectron_spin(2), nmo_max, nmo_min, &
      nmol_kind_global, nmolecule_global, nmolecule_global_read, &
      nmolecule_kind, nmolecule_kind_read, nmolecule_local, nspins_read, &
      restart_unit, source
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nao_molkind_read, &
                                                natom_molkind_read, &
                                                nmol_molkind_read
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: nel_molkind_read
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: eig_buffer, occ_buffer
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: mo_buffer
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, &
                                                occupation_numbers
    REAL(dp), DIMENSION(:, :), POINTER       :: mo
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                :: mol_blocks
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(mol_mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(mol_mo_set_type), POINTER           :: mo_set

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.

    group = globenv%group
    source = globenv%source

    method ="KG_GPW"
    nmolecule_kind = SIZE(kg_fm%kg_fm_mol_set)
    nmolecule_global = kg_fm%nmolecule_global
    nelectron_global = kg_fm%nelectron_global
    nao_max = kg_fm%nao_max
    nmo_max = kg_fm%nmo_max

    IF (globenv%ionode) THEN

      CALL section_vals_val_get(globenv%input_file,"FORCE_EVAL%DFT%RESTART_FILE_NAME",&
            c_val=file_name,error=error)
      IF (id_nr/=1) THEN
        file_name(LEN_TRIM(file_name):LEN(file_name)) = "-"//&
                    ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

!     *** read some info about the molecules ***

      ALLOCATE(nmol_molkind_read(nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(nao_molkind_read(nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(natom_molkind_read(nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(nel_molkind_read(2,nmolecule_kind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      READ (restart_unit) method, nmolecule_kind_read, nmolecule_global_read, nelectron_global_read, nspins_read
 
      IF (.NOT. (method=="KG_GPW")) THEN
        CALL stop_program(routineP,"READ RESTART: this is not a KG_GPW restart file")
      END IF
      IF (nmolecule_kind_read  /= nmolecule_kind) THEN
        CALL stop_program(routineP,"READ RESTART: different number of molecule kinds")
      END IF
      IF (nmolecule_global_read /= nmolecule_global_read) THEN
        CALL stop_program(routineP,"READ RESTART: different number of molecules ")
      END IF
      IF (nelectron_global_read /= nelectron_global_read) THEN
        CALL stop_program(routineP,"READ RESTART: different number of electrons ")
      END IF
      IF (nspins_read /= nspins) THEN
        CALL stop_program(routineP,"READ RESTART: nspin is not equal (program this !!)")
      END IF
      
      READ (restart_unit) nmol_molkind_read
      READ (restart_unit) nao_molkind_read
      READ (restart_unit) natom_molkind_read
      READ (restart_unit) nel_molkind_read(1,1:nmolecule_kind)
      READ (restart_unit) nel_molkind_read(2,1:nmolecule_kind)

      DO imolecule_kind = 1,nmolecule_kind
        NULLIFY(fm_mol_set)
        fm_mol_set => kg_fm%kg_fm_mol_set(imolecule_kind)
        CALL  get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,&
                                nmol_kind_global=nmol_kind_global,&
                                n_ao = n_ao,&
                                natom = natom, &
                                nelectron_spin=nelectron_spin) 
        CPPostcondition(nmol_kind_global==nmol_molkind_read(imolecule_kind),cp_warning_level,routineP,error,failure)
        CPPostcondition(n_ao==nao_molkind_read(imolecule_kind),cp_warning_level,routineP,error,failure)
        CPPostcondition(natom==natom_molkind_read(imolecule_kind),cp_warning_level,routineP,error,failure)
        IF(failure) THEN
         WRITE(globenv%scr,*) " READ RESTART : molecule kind ",imolecule_kind , " is not consistent "
         CALL stop_program(routineP,"error")
        END IF
      END DO 

      DEALLOCATE(nmol_molkind_read,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(nao_molkind_read,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(natom_molkind_read,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(nel_molkind_read,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF  ! ionode

    ALLOCATE(mo_buffer(nao_max,nmo_max,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(eig_buffer(nmo_max,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(occ_buffer(nmo_max,nspins),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO imol_global = 1,nmolecule_global
      imol_found = 0
      imol_proc = 0
      CALL dcopy(nao_max*nmo_max*nspins,0.0_dp,0,mo_buffer(1,1,1),1)
      CALL dcopy(nmo_max*nspins,0.0_dp,0,eig_buffer(1,1),1)
      CALL dcopy(nmo_max*nspins,0.0_dp,0,occ_buffer(1,1),1)
 
      IF(globenv%ionode) THEN
        DO ispin = 1,nspins
          READ (restart_unit) imol_global_read , ispin_read
          READ (restart_unit) n_mo_read(ispin),n_ao_read,&
                                nelectron_spin(ispin)
          READ (restart_unit) eig_buffer(1:n_mo_read(ispin),ispin)
          READ (restart_unit) occ_buffer(1:n_mo_read(ispin),ispin)
          DO imo = 1,n_mo_read(ispin)
           READ (restart_unit) mo_buffer(1:n_ao_read,imo,ispin)
          END DO
        END DO
      END IF
      CALL mp_bcast(n_mo_read,source,group)
      CALL mp_bcast(n_ao_read,source,group)

      CALL mp_bcast(mo_buffer,source,group)
      CALL mp_bcast(occ_buffer,source,group)
      CALL mp_bcast(eig_buffer,source,group)

      DO imolecule_kind = 1,nmolecule_kind
        NULLIFY(fm_mol_set,mol_blocks,mos, mo_set)
        NULLIFY(mo,eigenvalues,occupation_numbers)
        fm_mol_set => kg_fm%kg_fm_mol_set(imolecule_kind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,  &
                               nmolecule_local = nmolecule_local,&
                               n_ao=n_ao, n_mo=n_mo,nelectron_spin=nelectron_spin,&
                               fm_mol_blocks=mol_blocks) 
        DO imol_local = 1,nmolecule_local
           mol_block => mol_blocks(imol_local)
           IF(imol_global == mol_block%imol_global) THEN
             ! Molecule found, now take the read  MO etc
             ! check consistency
             ! fill the correct mo_set one spin at a time
             ! and ciaociao
             imol_found = 1
             imol_proc = globenv%mepos

             IF(n_ao /= n_ao_read) CALL stop_program(routineP,"RESTART READ: different number of AO")

             CALL get_fm_mol_block(fm_mol_block=mol_block, mos = mos)

             DO ispin = 1,nspins
                mo_set => mos(ispin)%mo_set
                CALL get_mol_mo_set(mo_set=mo_set,mo=mo,eigenvalues=eigenvalues,&
                                    occupation_numbers=occupation_numbers)
                nmo_min = MIN( n_mo(ispin),n_mo_read(ispin)) 

                DO imo = 1 , nmo_min
                  CALL DCOPY(n_ao,mo_buffer(1,imo,ispin),1,mo(1,imo),1)
                  eigenvalues(imo) = eig_buffer(imo,ispin)
                  occupation_numbers(imo) = occ_buffer(imo,ispin)
                END DO
                DO imo =  nmo_min+1 , n_mo(ispin)
                  CALL DCOPY(n_ao,0.0_dp,0,mo(1,imo),1)
                  eigenvalues(imo) = 0.0_dp
                  occupation_numbers(imo) = 0.0_dp
                END DO
             END DO
           END IF
                      
        END DO  ! imol_local

      END DO  ! imolecule_kind

    END DO  ! imol_global

    DEALLOCATE(mo_buffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(eig_buffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(occ_buffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (globenv%ionode) CALL close_file(unit_number=restart_unit)

    CALL timestop(0.0_dp,handle)


  END SUBROUTINE read_mo_mol_restart


!***************************************************************************
END MODULE  kg_gpw_fm_mol_utils




