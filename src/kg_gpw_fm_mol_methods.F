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
!!     Construct, initialize and use the full matrix blocks for the KG_GPW method 
!!     Using this method all the operations on the full matrixes can be
!!     splitted inoperations on the single molecular blocks, because
!!     each molecule is treated as independent, i.e. it does not interact
!!     directly with the other molecules 
!!
!!     It should contain also the routines for the :
!!                           initialization of mos (guess or restart)
!!                           orthogonalization
!!                           diagonalization
!!                           transfer of data to and from the sparse matrix
!!     
!!   AUTHOR
!!     MI (20.11.2004)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE kg_gpw_fm_mol_methods

  USE atomic_kind_types,               ONLY: get_atomic_kind
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
  USE distribution_1d_types,           ONLY: distribution_1d_type 
  USE kg_gpw_fm_mol_types,             ONLY: allocate_fm_mol_blocks,&
                                             allocate_mol_mo_set,&
                                             get_fm_mol_block,&
                                             get_mol_mo_set,&
                                             fm_mol_blocks_type,&
                                             init_mol_mo_set,&
                                             kg_fm_mol_set_type, &
                                             mol_mo_set_type,&
                                             set_kg_fm_mol_set
  USE kinds,                           ONLY: dp,dp_size,int_size
  USE molecule_kind_types,             ONLY: get_molecule_kind,molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "kg_gpw_fm_mol_methods"

! *** Define the derived structure types ***

! *** Public subroutines ***

  PUBLIC ::  build_local_fm_mol, calculate_mol_density_matrix,&
             mol_make_basis, multiply_sparse_mol_mo

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE  build_local_fm_mol(fm_mol_set,molecule_kind, imolkind,&
                                 local_molecules, molecule_set,&
                                 particle_set, nspins, nmo_eq_nao, &
                                 added_mos,use_cholesky,error)

!   Purpose: 

!   History: - Creation (20.11.2004,MI)

!   ***************************************************************************

    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind 
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: imolkind, nspins
    LOGICAL, INTENT(IN)                      :: nmo_eq_nao
    INTEGER, DIMENSION(2), INTENT(IN)        :: added_mos
    LOGICAL, INTENT(IN)                      :: use_cholesky
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_local_fm_mol',&
      routineP = module_name//':'//routineN

    TYPE(molecule_type), POINTER              :: molecule
    TYPE(fm_mol_blocks_type), POINTER         :: fm_mol_block

    INTEGER  :: iao, iat, iatom, iglobal, ikind, imol, ispin, istat, &
               multiplicity, n_ao, n_mo(2), natom, nelectron,&
               nelectron_spin(2), nmol_global, nmol_local, nsgf
    REAL(dp) :: maxocc
    INTEGER, DIMENSION(:), POINTER            :: index_mol_local
    INTEGER, DIMENSION(:), POINTER            :: molecule_list
    LOGICAL                                   :: failure

    failure = .FALSE.
 
    NULLIFY(molecule_list, index_mol_local, fm_mol_block, molecule)

    CALL get_molecule_kind(molecule_kind=molecule_kind,&
                            molecule_list=molecule_list,&
                            natom=natom,&
                            nelectron=nelectron,nsgf=n_ao)
 
!   If this molecule kind is charged this has to appear here.
!   In this case the number of elctrons is modified and we can calculate
!   the right number of molecular orbitals and te occupation numbers

    IF ((MODULO(nelectron,2) /= 0).AND.(nspins == 1)) THEN
      CALL stop_program(routineN,module_name,__LINE__,&
                   "Use the LSD option for an odd number of electrons.")
    END IF

    nmol_local = local_molecules%n_el(imolkind)
    index_mol_local => local_molecules%list(imolkind)%array

    nmol_global = SIZE(molecule_list,1)

    IF (MODULO(nelectron,2) == 0) THEN
      multiplicity = 1
    ELSE
      multiplicity = 2
    END IF


    IF (nspins == 1) THEN

      maxocc = 2.0_dp
      nelectron_spin(1) = nelectron
      nelectron_spin(2) = 0
      n_mo(1) = nelectron/2
      n_mo(2) = 0

    ELSE

      maxocc=1.0_dp

!     *** the simplist spin distribution is written here. Special cases will
!     *** need additional user input

      IF (MODULO(nelectron + multiplicity - 1,2) /= 0) THEN
        CALL stop_program(routineN,module_name,__LINE__,&
                          "LSD: try to use a different multiplicity.")
      END IF

      nelectron_spin(1) = (nelectron + multiplicity - 1)/2
      nelectron_spin(2) = (nelectron - multiplicity + 1)/2

      IF (nelectron_spin(2) < 0) THEN
        CALL stop_program(routineN,module_name,__LINE__,&
                         "LSD: too few electrons for this multiplicity.")
      END IF

      n_mo(1) = nelectron_spin(1)
      n_mo(2) = nelectron_spin(2)

    END IF

!    *** if additional mos need to be taken along, we do it here ***
     n_mo(1) = n_mo(1) + added_mos(1)
     IF (n_mo(2) /= 0) n_mo(2) = n_mo(2) + added_mos(2)

     IF(nmo_eq_nao) THEN
        n_mo(1) = n_ao
        IF (n_mo(2) /= 0) n_mo(2) = n_ao
     END IF
 
     CALL set_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,&
                           nmolecule_global=nmol_global,&
                           nmolecule_local=nmol_local,&
                           natom = natom, n_ao=n_ao,&
                           n_mo =n_mo , nelectron_spin=nelectron_spin,& 
                           maxocc=maxocc, multiplicity=multiplicity)

     CALL allocate_fm_mol_blocks(fm_mol_set,nmol_local,error=error)

     DO imol = 1,nmol_local
        fm_mol_block => fm_mol_set%fm_mol_blocks(imol)
        iglobal = index_mol_local(imol)
        molecule => molecule_set(iglobal)
        ALLOCATE(fm_mol_block%mos(nspins),STAT=istat)
        DO ispin = 1,nspins
          NULLIFY(fm_mol_block%mos(ispin)%mo_set)
          CALL allocate_mol_mo_set(fm_mol_block%mos(ispin)%mo_set,n_ao,n_mo(ispin),&
                                   maxocc,nelectron_spin(ispin),error=error)
          CALL init_mol_mo_set(fm_mol_block%mos(ispin)%mo_set,n_ao,&
                               n_mo(ispin),error=error)
        END DO
        ALLOCATE(fm_mol_block%ortho(n_ao,n_ao),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "ortho",n_ao*n_ao*dp_size)

        ALLOCATE(fm_mol_block%index_atom(natom),STAT=istat)
        ALLOCATE(fm_mol_block%index_kind(natom),STAT=istat)
        ALLOCATE(fm_mol_block%ifirst_ao(natom),STAT=istat)
        ALLOCATE(fm_mol_block%ilast_ao(natom),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
             "index_atom,index_kind,ifirst_ao,ilast_ao",natom*4*int_size)

        ALLOCATE(fm_mol_block%work(n_ao,n_ao),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                          "work",n_ao*n_ao*dp_size)

        iao = 0
        DO iat = 1,natom
          iatom = molecule%first_atom +iat -1
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               kind_number=ikind,nsgf=nsgf)
          fm_mol_block%index_atom(iat) = iatom
          fm_mol_block%index_kind(iat) = ikind
          fm_mol_block%ifirst_ao(iat)  = iao + 1
          iao = iao + nsgf
          fm_mol_block%ilast_ao(iat)  = iao
        END DO

     END DO

  END SUBROUTINE build_local_fm_mol

!***************************************************************************

  SUBROUTINE mol_make_basis(matrix_mo, ncol, matrix_ortho, otype,&
                            matrix_p, work, error)

    REAL(dp), DIMENSION(:,:), POINTER            :: matrix_mo
    INTEGER, INTENT(IN)                          :: ncol
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL  :: matrix_ortho
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL       :: otype
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL  :: matrix_p,work
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                  :: error

    CHARACTER(LEN=*), PARAMETER :: routineN ="mol_make_basis",&
         routineP = module_name//':'//routineN

    INTEGER :: i, info, istat, lda, ldb, ldc, nrow 
    LOGICAL                                      :: failure, found_type, &
                                                    use_cholesky, use_pmatrix, &
                                                    use_simple, use_sv 
    REAL(dp)                                     :: rone, rzero
    REAL(dp), DIMENSION(:,:), ALLOCATABLE        :: matrix_cc

!   ----------------------------------------------------------------------

    failure = .FALSE. 

    nrow = SIZE(matrix_mo,1)
    
    rone = 1.0_dp
    rzero = 0.0_dp

    IF (PRESENT(matrix_ortho)) THEN
       IF (.NOT. PRESENT(otype)) THEN
           CALL stop_program("mo_make_basis","ortho needs type !?")
       ENDIF
       found_type=.FALSE.
       use_sv=.FALSE.
       use_cholesky=.FALSE.
       IF (otype .EQ. "CHOLESKY") THEN
          found_type=.TRUE.
          use_cholesky=.TRUE.
       ENDIF
       IF (otype .EQ. "SV") THEN
          found_type=.TRUE.
          use_sv=.TRUE.
       ENDIF
       IF (.NOT. found_type) CALL stop_program("make_basis","sorry wrong type")
    ELSE
       use_simple=.TRUE.
    ENDIF

    ALLOCATE(matrix_cc(ncol,ncol) , STAT = istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_cc",ncol*ncol*dp_size)

    lda = SIZE(matrix_mo,1)
    ldb = SIZE(matrix_ortho,1)
    ldc = SIZE(matrix_cc,1)

    IF (use_pmatrix) THEN

    ELSE

    END IF

    IF (use_sv) THEN

      CALL DGEMM('T','N',ncol,ncol,nrow,rone,matrix_mo(1,1),lda,&
                 matrix_ortho(1,1),ldb,rzero,matrix_cc(1,1),ldc)

      CALL DPOTRF('U',ncol,matrix_cc(1,1),ldc,info)
      CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    ELSE
      write(*,*) 'make_basis, nosv option not checked' 
      stop 'how am I suppose to know matrix_cc'
      CALL DSYRK('U','T',ncol,nrow,rone,matrix_cc(1,1),lda,rzero,matrix_cc,ldc)
      CALL DPOTRF('U',ncol,matrix_cc(1,1),ldc,info)
      CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

    END IF

    CALL DTRSM('R','U','N','N',nrow,ncol,rone,matrix_cc(1,1),ncol,&
                 matrix_mo(1,1),nrow)


    DEALLOCATE(matrix_cc, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"matrix_cc")

  END SUBROUTINE mol_make_basis

!***************************************************************************

  SUBROUTINE calculate_mol_density_matrix(mo_set,mol_block,natom,p_ispin,error)

    TYPE(mol_mo_set_type), POINTER           :: mo_set 
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    TYPE(real_matrix_type), POINTER          ::  p_ispin
    INTEGER                                  ::  natom
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error


    CHARACTER(LEN=*), PARAMETER :: routineN = "calculate_mol_density_matrix",&
      routineP = module_name//':'//routineN

    INTEGER                                  :: col_first, col_last, &
                                                handle, homo, i, iat, iatom,&
                                                icol, icol_mol_p, irow, &
                                                irow_mol_p, ispin, istat, & 
                                                jat, jatom, lfomo, ncol, &
                                                ncol_max, nrow, nspins, &
                                                row_first, row_last

    INTEGER, DIMENSION(:), POINTER           :: i_atom,ifirst_ao,ilast_ao
    REAL(dp) :: alpha, maxocc
    REAL(dp), DIMENSION(:),  POINTER         :: occupation
    REAL(dp), DIMENSION(:,:),  POINTER       :: p_block
    REAL(dp), DIMENSION(:,:), POINTER        :: mol_p
    REAL(dp), DIMENSION(:,:), POINTER        :: mo_coeff
    REAL(dp), DIMENSION(:,:), POINTER        :: fmo_coeff
    LOGICAL :: failure, smear

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.
    smear = .FALSE.

    NULLIFY(i_atom,ifirst_ao,ilast_ao,mol_p)

    CALL get_fm_mol_block(fm_mol_block = mol_block,&
                          index_atom = i_atom,&
                          work = mol_p,&
                          ifirst_ao = ifirst_ao,&
                          ilast_ao = ilast_ao )

    CPPrecondition(ASSOCIATED(mol_p),cp_failure_level,routineP,error,failure)
    
    NULLIFY(mo_coeff,occupation)
    CALL get_mol_mo_set(mo_set=mo_set, mo = mo_coeff,&
                        occupation_numbers = occupation,&
                        maxocc = maxocc, &
                        homo=homo, lfomo=lfomo, nmo=ncol)

    CPPrecondition(ASSOCIATED(mo_coeff),cp_failure_level,routineP,error,failure)
    nrow = SIZE(mo_coeff,1)

    IF(lfomo .LE. homo) THEN 
      smear = .TRUE.
      NULLIFY(fmo_coeff)
      ALLOCATE(fmo_coeff(nrow,ncol),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,"fmo_coeff",nrow*ncol*dp_size)

      CALL DCOPY(nrow*ncol,mo_coeff,1,fmo_coeff,1) 
      DO icol = 1,ncol
        CALL DSCAL(nrow,occupation(icol),fmo_coeff(1,icol),1)
      END DO
      alpha = 1.0_dp
      ncol_max = ncol
    ELSE
      alpha = maxocc
      fmo_coeff => mo_coeff
      ncol_max = homo
    END IF

    CALL DGEMM('N','T',nrow,nrow,ncol_max,alpha,mo_coeff,nrow,&
               fmo_coeff,nrow,0.0_dp,mol_p,SIZE(mol_p,1))

    IF(smear) THEN
      DEALLOCATE(fmo_coeff,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,"fmo_coeff")
    ELSE
      NULLIFY(fmo_coeff)
    END IF

! Transfer the molecular block of the density matrix in the sparse matrix
    DO iat = 1,natom
      iatom = i_atom(iat)
      row_first = ifirst_ao(iat)
      row_last  = ilast_ao(iat)
      DO jat = iat , natom
        jatom = i_atom(jat)
        col_first = ifirst_ao(jat)
        col_last  = ilast_ao(jat)

        NULLIFY(p_block)
        IF (iatom <= jatom) THEN
          CALL get_block_node(matrix=p_ispin,&
                              block_row=iatom,&
                              block_col=jatom,&
                              BLOCK=p_block)
          icol =  1 
          DO icol_mol_p = col_first, col_last
            irow = 1
            DO irow_mol_p = row_first, row_last
              p_block(irow,icol) = mol_p(irow_mol_p,icol_mol_p)
              irow = irow + 1
            END DO  ! irow_mol_p 
            icol = icol + 1
          END DO   ! icol_mol_p
        ELSE
          CALL get_block_node(matrix=p_ispin,&
                              block_row=jatom,&
                              block_col=iatom,&
                              BLOCK=p_block)
          icol =  1 
          DO icol_mol_p = col_first, col_last
          irow = 1
            DO irow_mol_p = row_first, row_last
              p_block(icol,irow) = mol_p(irow_mol_p,icol_mol_p)
              irow = irow + 1
            END DO  
            icol = icol + 1
          END DO  ! icol_mol_p
        END IF      
      END DO  ! jat
    END DO  ! iat

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_mol_density_matrix

!***************************************************************************

  SUBROUTINE multiply_sparse_mol_mo(matrix,mol_block,natom,vectors,&
                                    nmo,mat_x_vec,error)

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(fm_mol_blocks_type), POINTER        :: mol_block
    INTEGER, INTENT(IN)                      :: natom,nmo
    REAL(dp), DIMENSION(:,:), POINTER        :: vectors
    REAL(dp), DIMENSION(:,:), POINTER        :: mat_x_vec
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "multiply_sparse_mol_mo",&
      routineP = module_name//":"//routineN

    INTEGER                                  :: col_first, col_last, &
                                                iat, iatom,&
                                                icol, icol_mol_s, irow, &
                                                irow_mol_s, & 
                                                jat, jatom, ldb, ldc, nrow, &
                                                row_first, row_last , i

    INTEGER, DIMENSION(:), POINTER           :: i_atom,ifirst_ao,ilast_ao
    REAL(dp), DIMENSION(:,:),  POINTER       :: s_block
    REAL(dp), DIMENSION(:,:), POINTER        :: mol_s
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    NULLIFY(i_atom,ifirst_ao,ilast_ao,mol_s)

    CALL get_fm_mol_block(fm_mol_block = mol_block,&
                          index_atom = i_atom,&
                          work = mol_s,&
                          ifirst_ao = ifirst_ao,&
                          ilast_ao = ilast_ao )

    CPPrecondition(ASSOCIATED(mol_s),cp_failure_level,routineP,error,failure)

    mol_s = 0.0_dp 
    nrow = SIZE(mol_s,1)
    ldb  =  SIZE(vectors,1)
    ldc  =  SIZE(mat_x_vec,1)
    
    CPPrecondition(ldc>=nrow,cp_failure_level,routineP,error,failure)
    CPPrecondition(ldb>=nrow,cp_failure_level,routineP,error,failure)

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
              mol_s(irow_mol_s,icol_mol_s) = s_block(irow,icol)
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
              mol_s(irow_mol_s,icol_mol_s) = s_block(icol,irow)
              irow = irow + 1
            END DO
            icol = icol + 1
          END DO
        END IF      
      END DO  ! jat
    END DO  ! iat
 
    CALL DGEMM('N','N',nrow,nmo,nrow,1.0_dp,mol_s,nrow,vectors,ldb,&
                0.0_dp,mat_x_vec,ldc) 

  END SUBROUTINE multiply_sparse_mol_mo

!***************************************************************************
END MODULE  kg_gpw_fm_mol_methods
