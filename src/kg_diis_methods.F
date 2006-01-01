!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004   CP2K developers group                         ! 
!-----------------------------------------------------------------------------!


!!****** cp2k/kg_diis_methods [1.0] *
!!
!!   NAME
!!     kg_diis_methods
!!
!!   FUNCTION
!!     Perform a direct inversion in the iterative subspace (DIIS).
!!
!!   AUTHOR
!!    MI 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_diis_methods
  USE kg_diis_types,                   ONLY: allocate_diis_block,&
                                             init_mol_diis_set,&
                                             kg_diis_buffer_type,&
                                             mol_diis_block_type,&
                                             mol_diis_buffer_set_type,&
                                             get_mol_diis_set
  USE kg_gpw_fm_mol_methods,           ONLY: calculate_mol_density_matrix,&
                                             fm_mol_eigensolver
  USE kg_gpw_fm_mol_types,             ONLY: fm_mol_blocks_type,&
                                             get_kg_fm_mol_set,&
                                             get_fm_mol_block,&
                                             get_mol_mo_set,&
                                             kg_fm_p_type,&
                                             kg_fm_mol_set_type,&
                                             mol_mo_set_p_type,&
                                             mol_mo_set_type
  USE kg_gpw_fm_mol_utils,             ONLY: copy_sparse2mol_block
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"


  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: module_name="kg_diis_methods"

! *** Global variables ***

  INTEGER, SAVE, PRIVATE :: last_diis_b_id=0

  ! *** Public subroutines ***

  PUBLIC :: alloc_mol_diis_iblock, kg_diis_step_A, kg_diis_step_B, &
            kg_diis_b_create, kg_diis_b_clear

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE alloc_mol_diis_iblock(scf_diis_buffer,mol_diis_iblock,nao,error)

    TYPE(kg_diis_buffer_type), POINTER         :: scf_diis_buffer
    TYPE(mol_diis_block_type) , POINTER        :: mol_diis_iblock
    INTEGER, INTENT(IN)                        :: nao
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error
    CHARACTER(len=*), PARAMETER :: routine_name="alloc_mol_diis_iblock",&
        routineP = module_name//':'//routine_name

    INTEGER                                    :: ibuffer, ispin, istat,&
                                                  nbuffer, nspins
    LOGICAL                                    :: failure

    failure = .FALSE.
    nbuffer = scf_diis_buffer%nbuffer
    nspins = scf_diis_buffer%nspins

    CPPrecondition(ASSOCIATED(mol_diis_iblock),cp_failure_level,routineP,error,failure)

    IF (.NOT.ASSOCIATED(mol_diis_iblock%ks_buffer)) THEN
      ALLOCATE(mol_diis_iblock%ks_buffer(nbuffer,nspins), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"mol_diis_iblock%ks_buffer(",0)

      DO ispin = 1,nspins
        DO ibuffer = 1,nbuffer
          NULLIFY(mol_diis_iblock%ks_buffer(ibuffer,ispin)%array)
          ALLOCATE(mol_diis_iblock%ks_buffer(ibuffer,ispin)%array(nao,nao),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineP,"ks_buffer(ibuffer,ispin)%array",&
                                           nao*nao*dp_size)
          mol_diis_iblock%ks_buffer(ibuffer,ispin)%array = 0.0_dp
        END DO 
      END DO

    END IF

    IF (.NOT.ASSOCIATED(mol_diis_iblock%er_buffer)) THEN
      ALLOCATE(mol_diis_iblock%er_buffer(nbuffer,nspins), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"mol_diis_iblock%er_buffer(",0)

      DO ispin = 1,nspins
        DO ibuffer = 1,nbuffer
          NULLIFY(mol_diis_iblock%er_buffer(ibuffer,ispin)%array)
          ALLOCATE(mol_diis_iblock%er_buffer(ibuffer,ispin)%array(nao,nao),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineP,"er_buffer(ibuffer,ispin)%array",&
                                           nao*nao*dp_size)
          mol_diis_iblock%er_buffer(ibuffer,ispin)%array = 0.0_dp
        END DO 
      END DO

    END IF

    IF (.NOT.ASSOCIATED(mol_diis_iblock%b_buffer)) THEN
      ALLOCATE (mol_diis_iblock%b_buffer(nbuffer+1,nbuffer+1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineP,"mol_diis_iblock%b_buffer",&
                                      (nbuffer+1)*(nbuffer+1)*dp_size)
      mol_diis_iblock%b_buffer = 0.0_dp
    END IF

  END SUBROUTINE alloc_mol_diis_iblock

! *****************************************************************************

  SUBROUTINE kg_diis_b_clear(scf_diis_buffer,error)

    TYPE(kg_diis_buffer_type), POINTER         :: scf_diis_buffer
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name="kg_diis_b_clear",&
        routineP = module_name//':'//routine_name

    LOGICAL                                    :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_diis_buffer),cp_failure_level,routine_name,error,failure)
    CPPrecondition(scf_diis_buffer%ref_count>=1,cp_failure_level,routine_name,error,failure)
    IF (.NOT. failure) THEN
       scf_diis_buffer%ncall=0
    END IF
    
  END SUBROUTINE kg_diis_b_clear

! *****************************************************************************

  SUBROUTINE kg_diis_b_create(scf_diis_buffer, kg_fm_mol_set,&
                              nbuffer,nspins, error)

!    FUNCTION
!       Allocates the kg_scf buffer for each molecule
! 
! *****************************************************************************

    TYPE(kg_diis_buffer_type), POINTER         :: scf_diis_buffer
    TYPE(kg_fm_mol_set_type), DIMENSION(:),&
      POINTER                                  :: kg_fm_mol_set
    INTEGER, INTENT(IN)                        :: nbuffer, nspins
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'kg_diis_b_create',&
     routineP = module_name//':'//routine_name

    INTEGER :: ikind, imol, istat, natom, nmol_global, nmol_kind, nmol_local
    LOGICAL :: failure
    TYPE(kg_fm_mol_set_type), POINTER          :: fm_mol_set
    TYPE(mol_diis_buffer_set_type), POINTER    :: diis_set

    failure=.FALSE.

    ALLOCATE (scf_diis_buffer, STAT = istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      scf_diis_buffer%nbuffer = nbuffer
      scf_diis_buffer%ncall = 0
      last_diis_b_id = last_diis_b_id+1
      scf_diis_buffer%id_nr=last_diis_b_id
      scf_diis_buffer%ref_count=1
      scf_diis_buffer%nspins = nspins
      nmol_kind =SIZE(kg_fm_mol_set,1)
      NULLIFY(scf_diis_buffer%mol_diis_buffer_set)

      ALLOCATE(scf_diis_buffer%mol_diis_buffer_set(nmol_kind), STAT = istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL init_mol_diis_set(scf_diis_buffer%mol_diis_buffer_set,nmol_kind)

      DO ikind = 1,nmol_kind
        fm_mol_set => kg_fm_mol_set(ikind)
        CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,natom=natom, &
                               nmol_kind_global=nmol_global,&
                               nmolecule_local=nmol_local)

        NULLIFY(diis_set)
        diis_set => scf_diis_buffer%mol_diis_buffer_set(ikind)
        CALL allocate_diis_block(diis_set,nmol_global,nmol_local)

      END DO   ! ikind
    END IF

  END SUBROUTINE kg_diis_b_create

! *****************************************************************************

  SUBROUTINE kg_diis_step_A(diis_buffer,ibuffer,nb_active,&
                            kg_fm_mol_set,diis_mat,&
                            nspins,diis_error,matrix_ks,overlap,error)

!   First part of the Diis update : create the diis matrix for each molecule

    TYPE(kg_diis_buffer_type), POINTER         :: diis_buffer
    INTEGER, INTENT(IN)                        :: ibuffer, nb_active
    TYPE(kg_fm_mol_set_type), DIMENSION(:),&
      POINTER                                  :: kg_fm_mol_set
    REAL(dp), DIMENSION(:,:), POINTER          :: diis_mat
    INTEGER, INTENT(IN)                        :: nspins
    REAL(dp), INTENT(INOUT)                    :: diis_error
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                  :: matrix_ks
    TYPE(real_matrix_type), POINTER, OPTIONAL  :: overlap
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = 'kg_diis_step_A', &
      routineP = module_name//':'//routine_name

    INTEGER :: handle, homo, icol, ib, irow, imol, imolecule_kind, ispin, istat, jb, &
               ldaux1, ldaux2, ldc, ldks, lder, mo_ncol, n_ao_kind, nat_mol,&
               nmol_local, nmolecule_kind, i
    LOGICAL                                    :: failure
    REAL(dp)                                   :: tmp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: aux2_matrix
    REAL(dp), DIMENSION(:,:), POINTER          :: aux1_matrix,&
                                                  mo_coeff
    REAL(dp), DIMENSION(:,:), POINTER          :: b_buffer
    REAL(dp), DIMENSION(:,:), POINTER          :: ks_ib
    REAL(dp), DIMENSION(:,:), POINTER          :: er_ib, er_jb
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                  :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER          :: mol_block
    TYPE(kg_fm_mol_set_type), POINTER          :: fm_mol_set
    TYPE(mol_mo_set_p_type), DIMENSION(:),&
      POINTER                                  :: mos
    TYPE(mol_mo_set_type), POINTER             :: mo_set
    TYPE(mol_diis_buffer_set_type), POINTER    :: mol_diis_set
    TYPE(mol_diis_block_type) , POINTER        :: mol_diis_iblock


!   ---------------------------------------------------------------------------

    CALL timeset("kg_diis_step_A","I","",handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(kg_fm_mol_set),cp_failure_level,routineP,error,failure)

    nmolecule_kind = SIZE(kg_fm_mol_set,1)

    !   *** Quick return, if no DIIS is requested ***
    IF (diis_buffer%nbuffer < 1) THEN
       CALL timestop(0.0_dp,handle)
       RETURN
    END IF

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

      NULLIFY(mol_diis_set)
      mol_diis_set => diis_buffer%mol_diis_buffer_set(imolecule_kind)

      DO imol = 1, nmol_local
        NULLIFY(mol_block,mos,aux1_matrix,mos)
        mol_block => fm_mol_blocks(imol)
        CALL get_fm_mol_block(fm_mol_block = mol_block,&
                              work = aux1_matrix,&
                              mos = mos)
        ldaux1 = SIZE(aux1_matrix,1)
        CPPrecondition(ldaux1>=n_ao_kind,cp_failure_level,routineP,error,failure)

        NULLIFY(mol_diis_iblock)
        CALL get_mol_diis_set(mol_diis_buffer_set=mol_diis_set,imol=imol,&
                              mol_diis_iblock=mol_diis_iblock,error=error)

        CALL alloc_mol_diis_iblock(diis_buffer,mol_diis_iblock,n_ao_kind,error=error)

        NULLIFY(b_buffer)
        mol_diis_iblock%mol_error_max = 0.0_dp 
        b_buffer => mol_diis_iblock%b_buffer

        DO ispin=1,nspins

          ! Copy the ks matrix of this molecule in the aux1_matrix 
           CALL copy_sparse2mol_block(matrix_ks(ispin)%matrix, mol_block,&
                                      aux1_matrix, nat_mol, n_ao_kind, &
                                      n_ao_kind, error=error)
           NULLIFY(mo_set,mo_coeff)
           mo_set => mos(ispin)%mo_set
           CALL get_mol_mo_set(mo_set=mo_set, mo = mo_coeff,&
                               nmo=mo_ncol, homo = homo)
           ldc = SIZE(mo_coeff,1)
           CPPrecondition(ldc>=n_ao_kind,cp_failure_level,routineP,error,failure) 

           ks_ib => mol_diis_iblock%ks_buffer(ibuffer,ispin)%array
           ldks = SIZE(ks_ib,1)
           CPPrecondition(ldks>=n_ao_kind,cp_failure_level,routineP,error,failure)

           CALL DCOPY(n_ao_kind*n_ao_kind,aux1_matrix,1,ks_ib,1)
           CALL DSYMM("L","U",n_ao_kind,homo,2.0_dp,ks_ib(1,1),ldks,mo_coeff(1,1),ldc,&
                      0.0_dp,aux1_matrix(1,1),ldaux1)

           er_ib => mol_diis_iblock%er_buffer(ibuffer,ispin)%array
           lder = SIZE(er_ib,1)
           CPPrecondition(lder>=n_ao_kind,cp_failure_level,routineP,error,failure)

           IF(PRESENT(overlap)) THEN
             ! Copy the overlap matrix of this molecule in the er_ib
             CALL copy_sparse2mol_block(overlap,mol_block,er_ib,nat_mol,n_ao_kind,&
                                        n_ao_kind, error=error)
             CALL DSYMM("L","U",n_ao_kind,homo,2.0_dp,er_ib(1,1),lder,&
                        mo_coeff(1,1),ldc,0.0_dp,aux2_matrix(1,1),ldaux2) 

             CALL DGEMM("N","T",n_ao_kind,n_ao_kind,homo,1.0_dp,aux2_matrix(1,1),ldaux2,&
                        aux1_matrix(1,1),ldaux1,0.0_dp,er_ib(1,1),lder)
             CALL DGEMM("N","T",n_ao_kind,n_ao_kind,homo,1.0_dp,aux1_matrix(1,1),ldaux1,&
                        aux2_matrix(1,1),ldaux2,-1.0_dp,er_ib(1,1),lder)

           ELSE

             CALL DGEMM("N","T",n_ao_kind,n_ao_kind,homo,1.0_dp,mo_coeff(1,1),ldc,&
                        aux1_matrix(1,1),ldaux1,0.0_dp,er_ib,lder)
             CALL DGEMM("N","T",n_ao_kind,n_ao_kind,homo,1.0_dp,aux1_matrix(1,1),ldaux1,&
                        mo_coeff(1,1),ldc,0.0_dp,er_ib,lder)
           END IF

           tmp = MAXVAL(ABS(er_ib(1:n_ao_kind,1:n_ao_kind)))
           mol_diis_iblock%mol_error_max = MAX(mol_diis_iblock%mol_error_max,tmp)

           DO jb = 1,nb_active
             IF(ispin==1) b_buffer(jb,ibuffer) = 0.0_dp
             er_ib => mol_diis_iblock%er_buffer(ibuffer,ispin)%array
             er_jb => mol_diis_iblock%er_buffer(jb,ispin)%array
             tmp = 0.0_dp
             DO icol = 1,n_ao_kind
                DO irow = 1,n_ao_kind
                   tmp = tmp + er_jb(irow,icol)*er_ib(irow,icol)
                END DO
              END DO
              b_buffer(jb,ibuffer) =  b_buffer(jb,ibuffer) + tmp
              IF(ispin==nspins) b_buffer(ibuffer,jb) = b_buffer(jb,ibuffer) 
            END DO

          END DO  ! ispin

         ! Update the global diis matrix with the contribution of this molecule
         DO ib = 1,nb_active
           DO jb = 1,nb_active
             diis_mat(jb,ib) = diis_mat(jb,ib) +  b_buffer(jb,ib)
           END DO
         END DO 
         ! Update the diis error
         diis_error = MAX(diis_error,mol_diis_iblock%mol_error_max) 

       END DO  ! imol

       DEALLOCATE(aux2_matrix,STAT = istat)
       IF (istat /= 0) CALL stop_memory(routineP,"aux2_matrix")

     END DO  ! imolecule_kind

     CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_diis_step_A

! *****************************************************************************

  SUBROUTINE kg_diis_step_B(diis_buffer,ibuffer,nb_active,&
                            kg_fm_set,scf_control,diis_ev,&
                            nspins,matrix_ks,matrix_p,&
                            do_iter,diis_step,do_level_shift,use_jacobi,error)

!   Second part of the DIIS update : modify the mos

    TYPE(kg_diis_buffer_type), POINTER         :: diis_buffer
    INTEGER, INTENT(IN)                        :: ibuffer, nb_active
    TYPE(kg_fm_p_type), POINTER                :: kg_fm_set
    TYPE(scf_control_type), POINTER            :: scf_control
    REAL(dp), DIMENSION(:), POINTER            :: diis_ev
    INTEGER, INTENT(IN)                        :: nspins
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                  :: matrix_ks
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                  :: matrix_p
    LOGICAL, INTENT(IN)                        :: do_iter,diis_step,&
                                                  do_level_shift,use_jacobi
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = 'kg_diis_step_B', &
      routineP = module_name//':'//routine_name

    INTEGER :: handle, imol, imolecule_kind, ispin, istat, jb, ldaux1, ldaux2, &
               ldks, n_ao_kind, nat_mol, nmol_local, nmolecule_kind, i
    LOGICAL                                    :: ionode, failure
    REAL(dp)                                   :: tmp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: aux2_matrix
    REAL(dp), DIMENSION(:,:), POINTER          :: aux1_matrix,&
                                                  ortho
    REAL(dp), DIMENSION(:,:), POINTER          :: ks_ib
    TYPE(cp_logger_type), POINTER              :: logger
    TYPE(fm_mol_blocks_type), DIMENSION(:), &
      POINTER                                  :: fm_mol_blocks
    TYPE(fm_mol_blocks_type), POINTER          :: mol_block
    TYPE(kg_fm_mol_set_type), DIMENSION(:),&
      POINTER                                  :: kg_fm_mol_set
    TYPE(kg_fm_mol_set_type), POINTER          :: fm_mol_set
    TYPE(mol_mo_set_p_type), DIMENSION(:),&
      POINTER                                  :: mos
    TYPE(mol_mo_set_type), POINTER             :: mo_set
    TYPE(mol_diis_buffer_set_type), POINTER    :: mol_diis_set
    TYPE(mol_diis_block_type) , POINTER        :: mol_diis_iblock

!   ---------------------------------------------------------------------------

    CALL timeset("kg_diis_step_B","I","",handle)

    NULLIFY(logger)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%source==logger%para_env%mepos
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(kg_fm_set),cp_failure_level,routineP,error,failure)

    kg_fm_mol_set => kg_fm_set%kg_fm_mol_set
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

      NULLIFY(mol_diis_set)
      mol_diis_set => diis_buffer%mol_diis_buffer_set(imolecule_kind)

      DO imol = 1, nmol_local
        NULLIFY(mol_block,mos,ortho,aux1_matrix)
        mol_block => fm_mol_blocks(imol)
        CALL get_fm_mol_block(fm_mol_block = mol_block,&
                               ortho = ortho,&
                               work = aux1_matrix,&
                               mos = mos)
        ldaux1 = SIZE(aux1_matrix,1)
        CPPrecondition(ldaux1>=n_ao_kind,cp_failure_level,routineP,error,failure)

        CALL get_mol_diis_set(mol_diis_buffer_set=mol_diis_set,imol=imol,&
                              mol_diis_iblock=mol_diis_iblock,error=error)

        DO ispin=1,nspins

          NULLIFY(mo_set)
          mo_set => mos(ispin)%mo_set

          IF(do_iter) THEN
             IF(diis_step) THEN

               aux1_matrix = 0.0_dp
               DO jb = 1,nb_active
                 ks_ib => mol_diis_iblock%ks_buffer(jb,ispin)%array
                 ldks = SIZE(ks_ib,1)
                 CPPrecondition(ldks>=n_ao_kind,cp_failure_level,routineP,error,failure)

                 CALL DAXPY(n_ao_kind*n_ao_kind,-diis_ev(jb),&
                            ks_ib(1,1),1,aux1_matrix(1,1),1)
               END DO 

             ELSE
               ks_ib => mol_diis_iblock%ks_buffer(ibuffer,ispin)%array
               ldks = SIZE(ks_ib,1)
               CPPrecondition(ldks>=n_ao_kind,cp_failure_level,routineP,error,failure)
               CALL DCOPY(n_ao_kind*n_ao_kind,ks_ib(1,1),1,aux1_matrix(1,1),1)
             END IF

          ELSE
             ! Copy the core hamiltonian of this molecule in the aux1_matrix 
              CALL copy_sparse2mol_block(matrix_ks(ispin)%matrix, mol_block,&
                                         aux1_matrix, nat_mol, n_ao_kind, &
                                         n_ao_kind, error=error)
          END IF

          ! Diagonalize the core hamiltonian and copy the eigenvectors in mo_set 
          CALL fm_mol_eigensolver(aux1_matrix,mo_set,ortho,aux2_matrix,n_ao_kind,&
                                   do_level_shift,scf_control%level_shift,&
                                   use_cholesky=scf_control%use_cholesky, &
                                   use_jacobi=use_jacobi,&
                                   jacobi_threshold=scf_control%jacobi_threshold,&
                                   smear=scf_control%smear, error=error)

          IF (logger%print_keys%each_scf_step) THEN
             ! give data to ionode
!             CALL give_mo_set(kg_fm_set,mo_set,ispin,mol_diis_iblock%imol_global,4,6,globenv)
          END IF

          CALL calculate_mol_density_matrix(mo_set,mol_block,&
                            nat_mol,matrix_p(ispin)%matrix,error=error)
        END DO  ! ispin

      END DO  ! imol

      DEALLOCATE(aux2_matrix,STAT = istat)
      IF (istat /= 0) CALL stop_memory(routineP,"aux2_matrix")

    END DO  ! imolecule_kind

    IF (logger%print_keys%each_scf_step .AND. ionode) THEN
!      CALL write_mo_set(kg_fm_set,4,6,globenv)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE kg_diis_step_B

! *****************************************************************************
END MODULE kg_diis_methods
!******************************************************************************

