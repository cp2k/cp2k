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
                                             get_kg_fm_mol_set,&
                                             get_mol_mo_set,&
                                             fm_mol_blocks_type,&
                                             init_mol_mo_set,&
                                             kg_fm_mol_set_type, &
                                             mol_mo_set_p_type,&
                                             mol_mo_set_type,&
                                             set_kg_fm_mol_set,&
                                             set_mol_mo_set
  USE kg_gpw_fm_mol_utils,             ONLY: copy_sparse2mol_block,&
                                             fm_mol_cholesky_restore,&
                                             fm_mol_syevd,&
                                             fm_mol_syevx
  USE kinds,                           ONLY: dp,dp_size,int_size
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: get_molecule_kind,molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qs_mo_types,                     ONLY: set_mo_occupation
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

! *** Public subroutines ***

  PUBLIC ::  build_local_fm_mol, calculate_mol_density_matrix,&
             fm_mol_eigensolver, fm_mol_orthonormality, &
             mol_make_basis, multiply_sparse_mol_mo

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE  build_local_fm_mol(fm_mol_set,molecule_kind, imolkind,&
                                 local_molecules, molecule_set,&
                                 particle_set, nspins, nmo_eq_nao, &
                                 added_mos,use_cholesky,&
                                 nelectron_global,nao_global,&
                                 nao_max, nmo_max, &
                                 maxocc_global,error)

!   Purpose: 

!   History: - Creation (20.11.2004,MI)

!   ***************************************************************************

    TYPE(kg_fm_mol_set_type), TARGET           :: fm_mol_set
    TYPE(molecule_kind_type), POINTER          :: molecule_kind 
    TYPE(distribution_1d_type), POINTER        :: local_molecules
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                  :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                  :: particle_set
    INTEGER, INTENT(IN)                        :: imolkind, nspins
    LOGICAL, INTENT(IN)                        :: nmo_eq_nao
    INTEGER, DIMENSION(2), INTENT(IN)          :: added_mos
    LOGICAL, INTENT(IN)                        :: use_cholesky
    INTEGER , INTENT(INOUT)                    :: nelectron_global
    INTEGER , INTENT(INOUT)                    :: nao_global
    INTEGER , INTENT(INOUT)                    :: nao_max, nmo_max
    REAL(dp), DIMENSION(2), INTENT(INOUT)      :: maxocc_global
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_local_fm_mol',&
      routineP = module_name//':'//routineN

    TYPE(molecule_type), POINTER              :: molecule

    INTEGER  :: iao, iat, iatom, iglobal, ikind, imol, ispin, istat, &
               multiplicity, n_ao, n_mo(2), natom, nelectron,&
               nelectron_spin(2), nmol_global, nmol_local, nsgf
    REAL(dp) :: maxocc
    INTEGER, DIMENSION(:), POINTER            :: index_mol_local
    INTEGER, DIMENSION(:), POINTER            :: molecule_list
    LOGICAL                                   :: failure
    TYPE(kg_fm_mol_set_type), POINTER         :: fm_mol_set_p

    failure = .FALSE.
 
    NULLIFY(molecule_list, index_mol_local, molecule)

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
      maxocc_global(1) = max( maxocc_global(1),maxocc)
      maxocc_global(2) = 0.0_dp
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
      maxocc_global(1) = max( maxocc_global(1),maxocc)
      maxocc_global(2) = max( maxocc_global(2),maxocc)

    END IF

!    *** if additional mos need to be taken along, we do it here ***
     n_mo(1) = n_mo(1) + added_mos(1)
     IF (n_mo(2) /= 0) n_mo(2) = n_mo(2) + added_mos(2)

     IF(nmo_eq_nao) THEN
        n_mo(1) = n_ao
        IF (n_mo(2) /= 0) n_mo(2) = n_ao
     END IF
     fm_mol_set_p => fm_mol_set
     CALL set_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set_p,&
                           nmol_kind_global=nmol_global,&
                           nmolecule_local=nmol_local,&
                           natom = natom, n_ao=n_ao,&
                           n_mo =n_mo , nelectron_spin=nelectron_spin,& 
                           maxocc=maxocc, multiplicity=multiplicity)

     CALL allocate_fm_mol_blocks(fm_mol_set,nmol_local,error=error)

     nelectron_global = nelectron_global + nmol_local * nelectron
     nao_global = nao_global + nmol_local * n_ao
     nao_max = MAX(nao_max,n_ao)
     nmo_max = MAX(nmo_max,n_mo(1),n_mo(2))
     DO imol = 1,nmol_local
        iglobal = index_mol_local(imol)
        fm_mol_set%fm_mol_blocks(imol)%imol_global = iglobal
        molecule => molecule_set(iglobal)
        ALLOCATE(fm_mol_set%fm_mol_blocks(imol)%mos(nspins),STAT=istat)
        DO ispin = 1,nspins
          NULLIFY(fm_mol_set%fm_mol_blocks(imol)%mos(ispin)%mo_set)
          CALL allocate_mol_mo_set(fm_mol_set%fm_mol_blocks(imol)%mos(ispin)%mo_set,&
                                   n_ao,n_mo(ispin),&
                                   maxocc,nelectron_spin(ispin),error=error)
          CALL init_mol_mo_set(fm_mol_set%fm_mol_blocks(imol)%mos(ispin)%mo_set,n_ao,&
                               n_mo(ispin),error=error)
        END DO
        ALLOCATE(fm_mol_set%fm_mol_blocks(imol)%ortho(n_ao,n_ao),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "ortho",n_ao*n_ao*dp_size)

        ALLOCATE(fm_mol_set%fm_mol_blocks(imol)%index_atom(natom),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "index_atom",natom*int_size)
        ALLOCATE(fm_mol_set%fm_mol_blocks(imol)%index_kind(natom),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "index_kind",natom*int_size)
        ALLOCATE(fm_mol_set%fm_mol_blocks(imol)%ifirst_ao(natom),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "ifirst_ao",natom*int_size)
        ALLOCATE(fm_mol_set%fm_mol_blocks(imol)%ilast_ao(natom),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "ilast_ao",natom*int_size)

        ALLOCATE(fm_mol_set%fm_mol_blocks(imol)%work(n_ao,n_ao),STAT=istat)
        IF(istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
                                        "work",n_ao*n_ao*dp_size)

        iao = 0
        DO iat = 1,natom
          iatom = molecule%first_atom +iat -1
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               kind_number=ikind,nsgf=nsgf)
          fm_mol_set%fm_mol_blocks(imol)%index_atom(iat) = iatom
          fm_mol_set%fm_mol_blocks(imol)%index_kind(iat) = ikind
          fm_mol_set%fm_mol_blocks(imol)%ifirst_ao(iat)  = iao + 1
          iao = iao + nsgf
          fm_mol_set%fm_mol_blocks(imol)%ilast_ao(iat)  = iao
        END DO

     END DO

  END SUBROUTINE build_local_fm_mol

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

   SUBROUTINE fm_mol_eigensolver(ks_mol_block, mo_set, ortho, aux_nxn,ldaux,&
                                 do_level_shift, level_shift, use_cholesky, &
                                 use_jacobi, jacobi_threshold ,&
                                 smear, error)
!
!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.
!   ks will be modified

!   History: - Creation (01.12.2004, Matthias Krack)
!
!   ***************************************************************************

    REAL(dp), DIMENSION(:,:), POINTER        :: ks_mol_block
    TYPE(mol_mo_set_type), POINTER           :: mo_set
    REAL(dp), DIMENSION(:,:), POINTER        :: ortho
    REAL(dp), DIMENSION(:,:)                 :: aux_nxn
    INTEGER, INTENT(IN)                      :: ldaux
    LOGICAL, INTENT(IN)                      :: do_level_shift
    REAL(KIND=dp), INTENT(IN)                :: level_shift
    LOGICAL, INTENT(INOUT)                   :: use_cholesky
    LOGICAL, INTENT(IN)                      :: use_jacobi
    REAL(KIND=dp), INTENT(IN)                :: jacobi_threshold, smear
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "fm_mol_eigensolver", &
      routineP = module_name//":"//routineN

    INTEGER                                  :: handle, homo, info, imo, &
                                                itype, istat, lfomo, &
                                                mo_ncol, nao, nelectron, nmo, p, q
    LOGICAL                                  :: failure,  mo_uocc
    REAL(KIND=dp)                            :: alpha, c, maxocc, tan_theta, tau, s
    REAL(dp), DIMENSION(:), POINTER          :: mo_eigenvalues,&
                                                mo_occupation
    REAL(dp), DIMENSION(:), ALLOCATABLE      :: vec_tmp
    REAL(dp), DIMENSION(:,:), POINTER        :: mo_coeff

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

!   *** Diagonalise the Kohn-Sham matrix ***
    NULLIFY (mo_coeff,mo_eigenvalues,mo_occupation)
    CALL get_mol_mo_set(mo_set=mo_set, mo = mo_coeff,&
                        eigenvalues = mo_eigenvalues,&
                        occupation_numbers = mo_occupation,&
                        nmo=mo_ncol)

    nao = SIZE(ks_mol_block,1)
!    ldaux = SIZE(aux_nxn,1)
    CPPostcondition(ldaux>=nao,cp_failure_level,routineP,error,failure)
    CPPostcondition(ldaux>=mo_ncol,cp_failure_level,routineP,error,failure)

    IF (do_level_shift) THEN
       IF (use_cholesky) THEN
          CALL stop_program("eigensolver","level shift not implemented")
       ENDIF
    ENDIF

    IF (use_cholesky) THEN
       ! Transform the generalized eigenvalue problem into a simple one
       itype = 1
       CALL DSYGST(itype,'U',nao,ks_mol_block(1,1),nao,ortho(1,1),nao,info)
       CPPostcondition(info==0,cp_failure_level,routineP,error,failure)

       IF (use_jacobi .OR. smear /= 0.0_dp) THEN
         ! Calculate all eigenvalues and eigenvectors
         CALL fm_mol_syevd(ks_mol_block,nao,mo_eigenvalues,aux_nxn,error=error)

         ! Transform the UC coefficients into the C coefficients 
         ! (U=ortho cholesky decomposition of the overlap matrix)
         CALL fm_mol_cholesky_restore(aux_nxn,mo_ncol,ortho,mo_coeff,"SOLVE")

         IF (use_jacobi ) THEN
           use_cholesky = .FALSE.
           CALL dcopy(nao*nao,mo_coeff(1,1),1,ortho(1,1),1)
         END IF
       ELSE
         ! Calculate selected eigenvalues and eigenvectors
         CALL fm_mol_syevx(ks_mol_block,aux_nxn,ldaux,ldaux,mo_eigenvalues,mo_ncol,error=error)
         ! Transform the UC coefficients into the C coefficients 
         ! (U=ortho cholesky decomposition of the overlap matrix)
         CALL fm_mol_cholesky_restore(aux_nxn,mo_ncol,ortho,mo_coeff,"SOLVE")
       END IF

    ELSE IF (use_jacobi) THEN
      CALL get_mol_mo_set(mo_set=mo_set,homo = homo) 

      CALL DSYMM("L","U",nao,homo,1.0_dp,ks_mol_block(1,1),nao,&
                 mo_coeff(1,1),nao,0.0_dp,aux_nxn(1,1),ldaux)
      CALL DGEMM("T","N",homo,nao-homo,nao,1.0_dp,aux_nxn(1,1),ldaux,&
                 mo_coeff(1,homo+1),nao,0.0_dp,ks_mol_block(1,homo+1),nao)

!     *** Klassisches Jacobi-Block-Verfrahren ***
!     *** Rotiere Elemente nur, wenn sie groesser als thresh sind ***

      ALLOCATE(vec_tmp(nao), STAT =istat)
      IF(istat/= 0.0_dp) CALL stop_memory(routineP,"allocate vec_tmp",nao*dp_size)
      DO q = homo+1,nao
        DO p =1,homo
          IF(ABS(ks_mol_block(p,q))>jacobi_threshold) THEN
            tau = (mo_eigenvalues(q)-mo_eigenvalues(p))/(2.0_dp*ks_mol_block(p,q))
            tan_theta = SIGN(1.0_dp,tau)/(ABS(tau)+SQRT(1.0_dp+tau*tau)) 
            ! cos theta
            c = 1.0_dp/SQRT(1.0_dp+tan_theta*tan_theta)
            s = tan_theta*c
           ! Und jetzt noch die Eigenvektoren produzieren:
           ! Q * J
           !  Verstaendliche Version (bevor die BLAS-Aufrufe sie ersetzt haben)
           !  vec_tmp(p) = c*mo_coeff(:,p) - s*mo_coeff(:,q)
           !  vec_tmp(q) = s*mo_coeff(:,p) + c*mo_coeff(:,q)

           !  mo_coeff(:,p)=vec_tmp(p)
           !  mo_coeff(:,q)=vec_tmp(q)

           CALL DCOPY(nao,mo_coeff(1,p),1,vec_tmp(1),1)
           CALL DSCAL(nao,c,mo_coeff(1,p),1)
           CALL DAXPY(nao,-s,mo_coeff(1,q),1,mo_coeff(1,p),1)
           CALL DSCAL(nao,c,mo_coeff(1,q),1)
           CALL DAXPY(nao,s,vec_tmp(1),1,mo_coeff(1,q),1)

          END IF
        END DO  ! p
      END DO  ! q
      DEALLOCATE(vec_tmp, STAT =istat)
      IF(istat/= 0.0_dp) CALL stop_memory(routineP,"deallocate vec_tmp")
!
    ELSE

      CALL DSYMM("L","U",nao,nao,1.0_dp,ks_mol_block(1,1),nao,&
                 ortho(1,1),nao,0.0_dp,aux_nxn(1,1),ldaux)
      CALL DGEMM("T","N",nao,nao,nao,1.0_dp,ortho(1,1),nao,aux_nxn(1,1),& 
                 ldaux,0.0_dp,ks_mol_block(1,1),nao)

       IF (do_level_shift) THEN
          CALL get_mol_mo_set(mo_set=mo_set,homo = homo) 
          DO imo=homo+1,mo_ncol
            ks_mol_block(imo,imo)=ks_mol_block(imo,imo)+level_shift
          END DO
       END IF

       ! Calculate selected eigenvalues and eigenvectors (U^-1C)
       CALL fm_mol_syevx(ks_mol_block,aux_nxn,ldaux,ldaux,mo_eigenvalues,mo_ncol,error=error)
       ! Transform the eigenvectors U(U^-1C)=C
       CALL DGEMM("N","N",nao,mo_ncol,nao,1.0_dp,ortho(1,1),nao,aux_nxn(1,1),nao,&
                  0.0_dp,mo_coeff(1,1),nao)

       IF (do_level_shift) THEN

!     *** Use last MO set as orthogonalization matrix ***
          DO imo=homo+1,mo_ncol
            mo_eigenvalues(imo) = mo_eigenvalues(imo) - level_shift
          END DO

          CALL DCOPY(nao*nao,mo_coeff(1,1),1,ortho(1,1),1)

       END IF

    END IF

    mo_uocc = .TRUE.
    IF (smear /= 0.0_dp) THEN
     CALL get_mol_mo_set(mo_set=mo_set,lfomo=lfomo, homo=homo,&
                         maxocc=maxocc, nelectron=nelectron)
     CALL set_mo_occupation(mo_occupation,mo_eigenvalues,&
                          homo,lfomo,maxocc,nelectron,mo_uocc,&
                          halfhole=-1,smear=smear)
     CALL  set_mol_mo_set(mo_set=mo_set,homo=homo,lfomo=lfomo)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE fm_mol_eigensolver

!***************************************************************************

  SUBROUTINE fm_mol_orthonormality(orthonormality,kg_fm_mol_set,overlap,nspins,error)

    REAL(dp), DIMENSION(:)                     :: orthonormality
    TYPE(kg_fm_mol_set_type), DIMENSION(:),&
      POINTER                                  :: kg_fm_mol_set
    TYPE(real_matrix_type), POINTER, OPTIONAL  :: overlap
    INTEGER, INTENT(IN)                        :: nspins
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = 'fm_mol_orthonormality', &
      routineP = module_name//':'//routine_name

    INTEGER :: icol, imol, imolecule_kind, irow, istat, ispin, ldaux2, ldc, &
               ldcsc, mo_ncol, mo_ncol_old, n_ao_kind, nat_mol, nmol_local, nmolecule_kind
    LOGICAL :: failure
    REAL(dp)                                   :: csc_val, max_csc_offd
    REAL(dp), DIMENSION(:,:), POINTER          :: aux1_matrix,&
                                                  aux2_matrix,&
                                                  csc_matrix,&
                                                  mo_coeff
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                  :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER          :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER          :: fm_mol_set
    TYPE(mol_mo_set_p_type), DIMENSION(:),&
      POINTER                                  :: mos
    TYPE(mol_mo_set_type), POINTER             :: mo_set

! -----------------------------------------------------------------------------------

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(kg_fm_mol_set),cp_failure_level,routineP,error,failure)

    nmolecule_kind = SIZE(kg_fm_mol_set,1)

    DO imolecule_kind = 1,nmolecule_kind
      NULLIFY(fm_mol_set,fm_mol_blocks)
      fm_mol_set => kg_fm_mol_set(imolecule_kind)
      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set, &
                             natom = nat_mol, &
                             n_ao = n_ao_kind , &
                             fm_mol_blocks = fm_mol_blocks,&
                             nmolecule_local = nmol_local)

      ALLOCATE(aux2_matrix(n_ao_kind,n_ao_kind), STAT = istat)
      IF (istat /= 0) CALL stop_memory(routineP,"aux2_matrix",&
                                       n_ao_kind*n_ao_kind*dp_size)
      ldaux2 = n_ao_kind
      mo_ncol_old = 0
      NULLIFY(aux1_matrix,csc_matrix)
      DO imol = 1, nmol_local
        NULLIFY(mol_block,mos,mos)
        mol_block => fm_mol_blocks(imol)
        CALL get_fm_mol_block(fm_mol_block = mol_block,mos=mos)

        ! Copy the overlap matrix of this molecule in the aux2_matrix 
         CALL copy_sparse2mol_block(overlap, mol_block,&
                                    aux2_matrix, nat_mol, n_ao_kind, &
                                    n_ao_kind, error=error)
        max_csc_offd = 0.0_dp
        DO ispin=1,nspins
           NULLIFY(mo_set,mo_coeff)
           mo_set => mos(ispin)%mo_set
           CALL get_mol_mo_set(mo_set=mo_set, mo = mo_coeff,&
                               nmo=mo_ncol)
           ldc = SIZE(mo_coeff,1)
           CPPrecondition(ldc>=n_ao_kind,cp_failure_level,routineP,error,failure) 

           IF(mo_ncol > mo_ncol_old) THEN
             CALL reallocate(aux1_matrix,1,n_ao_kind,1,mo_ncol)
             CALL reallocate(csc_matrix,1,mo_ncol,1,mo_ncol)
             ldcsc = mo_ncol
             mo_ncol_old = mo_ncol
           END IF
 
           CALL DGEMM('N','N',n_ao_kind,mo_ncol,n_ao_kind,1.0_dp,aux2_matrix,ldaux2,&
                      mo_coeff,ldc,0.0_dp,aux1_matrix,n_ao_kind)
           CALL DGEMM('T','N',mo_ncol,mo_ncol,n_ao_kind,1.0_dp,mo_coeff,ldc,&
                      aux1_matrix,n_ao_kind,0.0_dp,csc_matrix,ldcsc)

           DO icol = 1,mo_ncol
             DO irow =  1,mo_ncol
               csc_val = csc_matrix(irow,icol)
               IF(irow==icol) csc_val  = csc_val - 1.0_dp
               max_csc_offd = MAX(max_csc_offd,ABS(csc_val))
             END DO
           END DO
        END DO  ! ispin

        orthonormality(fm_mol_blocks(imol)%imol_global) = max_csc_offd
      END DO  ! imol

      IF(ASSOCIATED(aux1_matrix)) THEN
        DEALLOCATE(aux1_matrix, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineP,"deall. aux1_matrix")
      END IF
      IF(ASSOCIATED(csc_matrix)) THEN
        DEALLOCATE(csc_matrix, STAT = istat)
        IF (istat /= 0) CALL stop_memory(routineP,"deall. csc_matrix")
      END IF
      DEALLOCATE(aux2_matrix, STAT = istat)
      IF (istat /= 0) CALL stop_memory(routineP,"deall. aux2_matrix")

    END DO  ! imolecule_kind

  END SUBROUTINE fm_mol_orthonormality

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

! *****************************************************************************

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
