!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE kg_gpw_fm_mol_types,             ONLY: get_fm_mol_block,&
                                             get_kg_fm_mol_set,&
                                             get_mol_mo_set,&
                                             fm_mol_blocks_type,&
                                             kg_fm_mol_set_type,&
                                             mol_mo_set_p_type,&
                                             mol_mo_set_type
  USE kinds,                           ONLY: dp,dp_size,int_size
  USE memory_utilities,                ONLY: reallocate
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

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_fm_mol_utils"

! *** Public subroutines ***

  PUBLIC :: calculate_w_matrix_per_molecule, &
            copy_sparse2mol_block, fm_mol_cholesky_restore, &
            fm_mol_power, fm_mol_syevd, fm_mol_syevx, mol_mo_random

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

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_w_matrix_per_molecule", &
      routineP = module_name//':'//routine_name

    INTEGER :: homo, iao, iat_mol, iatom, ifs, ils, imo, imol, &
              imolecule_kind, inao, istat, jao, jatom, jat_mol,&
              jfs, jls, jnao, ldao, ldmo, nao_mol, nat_mol, &
               nmo_mol(2), nmol_local, nmolecule_kind, nocc
    INTEGER, DIMENSION(:), POINTER            :: i_atom
    INTEGER, DIMENSION(:), POINTER            :: ifirst_ao
    INTEGER, DIMENSION(:), POINTER            :: ilast_ao
    LOGICAL                                   :: failure
    REAL(dp)                                  :: eigocc
    REAL(dp), DIMENSION(:), POINTER           :: eigenvalues,&
                                                 occupation
    REAL(dp), DIMENSION(:,:), POINTER         :: wmo_x_mot
    REAL(dp), DIMENSION(:,:), POINTER         :: mo, wmo
    REAL(dp), DIMENSION(:,:), POINTER         :: w_block

    TYPE(fm_mol_blocks_type), DIMENSION(:),&
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

    REAL(dp), DIMENSION(:,:)                       :: vectors_in
    INTEGER, INTENT(IN)                            :: nmo
    REAL(dp), DIMENSION(:,:), POINTER              :: upper_dec
    REAL(dp), DIMENSION(:,:), POINTER              :: vectors_out
    CHARACTER ( LEN = * ), INTENT ( IN )           :: op
    CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: pos
    TYPE(cp_error_type), INTENT(inout), OPTIONAL   :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "fm_mol_cholesky_restore",&
      routineP = module_name//":"//routineN

    CHARACTER                                      :: my_pos
    INTEGER                                        :: nao
    LOGICAL                                        :: failure
    REAL(dp)                                       :: alpha
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

    REAL(dp), DIMENSION(:,:), POINTER            :: matrix
    REAL(dp), DIMENSION(:,:)                     :: eigenvectors
    INTEGER, INTENT(IN)                          :: ldv, ncv
    REAL(dp), INTENT(IN)                         :: exponent, threshold
    INTEGER , INTENT(OUT)                        :: n_dependent

    CHARACTER(LEN=*), PARAMETER :: routineN ="fm_mol_power",&
         routineP = module_name//':'//routineN

    INTEGER                                      :: handle, icol, istat,&
                                                    ncol,nrow
    REAL(dp)                                     :: f,p
    REAL(dp), DIMENSION(:), ALLOCATABLE          :: eigenvalues

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

    REAL(dp), DIMENSION(:,:), POINTER        :: matrix
    INTEGER, INTENT(IN)                      :: nrow
    REAL(dp), DIMENSION(:)                   :: eigenvalues
    REAL(dp), DIMENSION(:,:)                 :: matrix_out
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

!   *** Local parameters ***
    CHARACTER(LEN=*), PARAMETER :: routineN = "fm_mol_syevd", &
      routineP = module_name//":"//routineN

!   *** Local variables ***
    INTEGER                                  :: istat, info, lda, liwork, lwork
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

    REAL(dp), DIMENSION(:,:), POINTER            :: matrix
    REAL(dp), DIMENSION(:,:)                     :: eigenvetors
    INTEGER,INTENT(IN)                           :: ldv,ncv
    REAL(dp), DIMENSION(:), INTENT(OUT)          :: eigenvalues
    INTEGER, INTENT(IN)                          :: ncol
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN ="fm_mol_syevx",&
         routineP = module_name//':'//routineN

    REAL(KIND = dp), PARAMETER                   :: vl = 0.0_dp,&
                                                    vu = 0.0_dp 

!   *** Local variables ***
    INTEGER  :: info, istat, ldb, liwork, lwork, m, nb, nrow
    INTEGER,  DIMENSION(:), ALLOCATABLE          :: ifail,iwork
    REAL(dp)                                     :: abstol
    REAL(dp), DIMENSION(:), ALLOCATABLE          :: w, work
    REAL(dp), DIMENSION(:,:), POINTER            :: a,z

    REAL(KIND = dp), EXTERNAL                    :: dlamch
    INTEGER, EXTERNAL                            :: ilaenv
    LOGICAL :: failure

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

    REAL(dp), DIMENSION(:,:), POINTER        :: mo
    INTEGER, INTENT(IN)                      :: mepos, num_pe
    INTEGER, INTENT(IN), OPTIONAL            :: start_col, ncol
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN ="mol_mo_random"

    INTEGER                                  :: icol, my_start, my_ncol, nrow
    INTEGER, DIMENSION(4), SAVE              :: ISEED
    LOGICAL, SAVE                            :: FIRST = .TRUE.
    LOGICAL                                  :: failure

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
    REAL(dp), DIMENSION(:,:), POINTER        :: matrix_mol_block
    INTEGER, INTENT(IN)                      :: natom, nrow, ncol
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "copy_sparse2mol_block",&
      routineP = module_name//":"//routineN

    INTEGER                                  :: col_first, col_last, &
                                                iat, iatom,&
                                                icol, icol_mol_s, irow, &
                                                irow_mol_s, & 
                                                jat, jatom, &
                                                row_first, row_last , i
    INTEGER, DIMENSION(:), POINTER           :: i_atom,ifirst_ao,ilast_ao
    REAL(dp), DIMENSION(:,:),  POINTER       :: s_block
    LOGICAL                                  :: failure

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
END MODULE  kg_gpw_fm_mol_utils




