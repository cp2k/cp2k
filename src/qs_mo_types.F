!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qs_mo_types [1.0] *
!!
!!   NAME
!!     qs_mo_types
!!
!!   FUNCTION
!!     Definition and initialisation of the mo data type.
!!
!!   AUTHOR
!!     Matthias Krack (09.05.2001,MK)
!!
!!   MODIFICATION HISTORY
!!     - adapted to the new QS environment data structure (02.04.2002,MK)
!!     - set_mo_occupation added (17.04.02,MK)
!!     - correct_mo_eigenvalues added (18.04.02,MK)
!!     - calculate_density_matrix moved from qs_scf to here (22.04.02,MK)
!!     - mo_set_p_type added (23.04.02,MK)
!!     - PRIVATE attribute set for TYPE mo_set_type (23.04.02,MK)
!!     - calculate_w_matrix added for the calculation of the Pulay forces
!!       (23.05.02,MK)
!!     - started conversion to LSD (1.2003, Joost VandeVondele)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_mo_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_sm_fm_interactions,           ONLY: cp_sm_plus_fm_fm_t
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm, &
                                             cp_fm_column_scale, &
                                             cp_fm_syrk
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_types,                     ONLY: cp_fm_set_all, &
                                             cp_fm_get_element, &
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_type,&
                                             cp_fm_create, &
                                             cp_fm_to_fm, &
                                             cp_fm_release, &
                                             cp_fm_get_submatrix, &
                                             cp_fm_set_submatrix
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_sum, &
                                             mp_max
  USE orbital_pointers,                ONLY: nso
  USE output_utilities,                ONLY: print_warning
  USE particle_types,                  ONLY: particle_type
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             allocate_matrix,&
                                             get_block_node,&
                                             real_matrix_type,&
                                             scale_matrix,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: module_name = "qs_mo_types"

  TYPE mo_set_type
    PRIVATE
    REAL(wp)                           :: maxocc
    INTEGER                            :: homo,lfomo,nao,nelectron,nmo
    REAL(wp), DIMENSION(:), POINTER    :: eigenvalues,occupation_numbers
    TYPE(cp_fm_type), POINTER :: eigenvectors
  END TYPE mo_set_type

  TYPE mo_set_p_type
    TYPE(mo_set_type), POINTER :: mo_set
  END TYPE mo_set_p_type

! *** Public data types ***

  PUBLIC :: mo_set_p_type,&
            mo_set_type

! *** Public subroutines ***

  PUBLIC :: allocate_mo_set,&
            calculate_density_matrix,&
            calculate_w_matrix,&
            correct_mo_eigenvalues,&
            deallocate_mo_set,&
            get_mo_set,&
            init_mo_set,&
            read_mo_set,&
            set_mo_occupation,&
            write_mo_set, &
            calculate_subspace_eigenvalues, &
            calculate_orthonormality

! *****************************************************************************

  INTERFACE calculate_density_matrix
    MODULE PROCEDURE calculate_density_matrix_sparse,&
                     calculate_density_matrix_blacs
  END INTERFACE

  INTERFACE calculate_subspace_eigenvalues
    MODULE PROCEDURE subspace_eigenvalues_ks ! outdated ,subspace_eigenvalues_vec
  END INTERFACE

  INTERFACE calculate_w_matrix
    MODULE PROCEDURE calculate_w_matrix_1,&
                     calculate_w_matrix_2
  END INTERFACE

  INTERFACE read_mo_set
    MODULE PROCEDURE read_mo_set_from_restart
  END INTERFACE

  INTERFACE write_mo_set
    MODULE PROCEDURE write_mo_set_to_output_unit,write_mo_set_to_restart
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

!!****f* qs_mo_types/allocate_mo_set [1.0] *
!!
!!   NAME
!!     allocate_mo_set
!!
!!   SYNOPSIS
!!     Subroutine allocate_mo_set(mo_set, nao, nmo, nelectron, maxocc,&
!!         error)
!!       Type(mo_set_type), Pointer:: mo_set
!!       Real(Kind=wp), Intent (IN):: maxocc
!!       Integer, Intent (IN):: nao, nelectron, nmo
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine allocate_mo_set
!!
!!   FUNCTION
!!     allocate a mo set and partially initializes it (nao,nmo,nelectron
!!     are vaild).
!!     For the full initialization you need to call init_mo_set
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - mo_set: the mo_set to allocate
!!     - nao: number of atom orbitals
!!     - nmo: number of molecular orbitals
!!     - nelectrons: number of electrons
!!     - maxocc: maximum occupation of an orbital (LDA: 2, LSD:1)
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     11.2002 splitted initialization in two phases [fawzi]
!!
!!*** *********************************************************************

  SUBROUTINE allocate_mo_set(mo_set,nao,nmo,nelectron,maxocc, error)

!   Purpose: Allocate a wavefunction data structure.

!   History: - Creation (15.05.2001, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER                :: mo_set
    REAL(wp), INTENT(IN)                      :: maxocc
    INTEGER, INTENT(IN)                       :: nao,nelectron,nmo
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_mo_set (MODULE qs_mo_types)"

!   *** Local variables ***

    INTEGER :: istat
    LOGICAL :: failure

!   ---------------------------------------------------------------------------

    failure=.FALSE.

    ALLOCATE (mo_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set",0)

    mo_set%maxocc = maxocc
    mo_set%nelectron = nelectron
    mo_set%nao = nao
    mo_set%nmo = nmo
    mo_set%homo = 0
    mo_set%lfomo = 0

    NULLIFY (mo_set%eigenvalues)
    NULLIFY (mo_set%occupation_numbers)
    NULLIFY (mo_set%eigenvectors)

  END SUBROUTINE allocate_mo_set

!****************************************************************************

!!****f* qs_mo_types/init_mo_set [1.0] *
!!
!!   NAME
!!     init_mo_set
!!
!!   SYNOPSIS
!!     Subroutine init_mo_set(mo_set, fm_pool, name, error)
!!       Type(mo_set_type), Pointer:: mo_set
!!       Type(cp_fm_pool_type), Pointer:: fm_pool
!!       Character(Len=*), Intent (IN):: name
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine init_mo_set
!!
!!   FUNCTION
!!     initializes an allocated mo_set.
!!     eigenvalues, eigenvectors, occupation_numbers are valid only
!!     after this call.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - mo_set: the mo_set to initialize
!!     - fm_pool: a pool out which you initialize the mo_set
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 rewamped [fawzi]
!!
!!*** *********************************************************************

  SUBROUTINE init_mo_set(mo_set, fm_pool, name, error)

    TYPE(mo_set_type), POINTER                   :: mo_set
    TYPE(cp_fm_pool_type), POINTER               :: fm_pool
    CHARACTER(len=*), INTENT(in)                 :: name
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

    LOGICAL :: failure
    INTEGER :: istat, nmo,nao
    CHARACTER(len=*), PARAMETER :: routineN='init_mo_set',&
         routineP=module_name//':'//routineN

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(mo_set%eigenvectors),cp_failure_level,routineP,error,failure)

    ALLOCATE (mo_set%eigenvalues(mo_set%nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"mo_set%eigenvalues",mo_set%nmo*wp_size)
    END IF
    mo_set%eigenvalues(:) = 0.0_wp

    ALLOCATE (mo_set%occupation_numbers(mo_set%nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"mo_set%occupation_numbers",mo_set%nmo*wp_size)
    END IF
    CALL set_mo_occupation(mo_set=mo_set,smear=0.0_wp)

    CALL fm_pool_create_fm(fm_pool,mo_set%eigenvectors,&
         name=name,error=error)
    CALL cp_fm_get_info(mo_set%eigenvectors, nrow_global=nao, &
         ncol_global=nmo,error=error)
    CPPostcondition(nao>=mo_set%nao,cp_failure_level,routineP,error,failure)
    CPPostcondition(nmo>=mo_set%nmo,cp_failure_level,routineP,error,failure)

  END SUBROUTINE init_mo_set

! *****************************************************************************

  SUBROUTINE calculate_density_matrix_blacs(mo_set,density_matrix,error)

!   Purpose: Calculate the density matrix from the MO eigenvectors and the
!            MO occupation numbers.

!   History: - Creation (01.05.2001,MK)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER        :: density_matrix
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,imo

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix_blacs","I","",handle)

    CALL cp_fm_syrk("U","N",mo_set%lfomo-1,mo_set%maxocc,mo_set%eigenvectors,&
                    1,1,0.0_wp,density_matrix)

!   *** Add the density contribution from the fractional occupied MOs ***

    DO imo=mo_set%lfomo,mo_set%homo
      CALL cp_fm_syrk("U","N",1,mo_set%occupation_numbers(imo),&
                      mo_set%eigenvectors,1,imo,1.0_wp,density_matrix)
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_density_matrix_blacs

! *****************************************************************************

  SUBROUTINE calculate_density_matrix_sparse(mo_set,density_matrix,error)

!   Purpose: Calculate the density matrix

!   History: - Creation (06.2002, Joost VandeVondele)
!            - Fractional occupied orbitals (MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER           :: density_matrix
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local variables ***

    INTEGER :: handle,imo

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix_sparse","I","",handle)

    CALL set_matrix(density_matrix,0.0_wp)

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=mo_set%eigenvectors,&
                            first_col=1,&
                            last_col=mo_set%lfomo-1,&
                            alpha=mo_set%maxocc)

    DO imo=mo_set%lfomo,mo_set%homo
      CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=mo_set%eigenvectors,&
                              first_col=imo,&
                              last_col=imo,&
                              alpha=mo_set%occupation_numbers(imo))
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_density_matrix_sparse

! *****************************************************************************

  SUBROUTINE calculate_w_matrix_1(mo_set,w_matrix)

!   Purpose: Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!            and the MO occupation numbers.

!   History: - Creation (03.03.03,MK)
!            - Modification that computes it as a full block, several times (e.g. 20) 
!              faster at the cost of some additional memory

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER           :: w_matrix
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(cp_fm_type), POINTER        :: weighted_vectors

!   *** Local variables ***

    REAL(wp), DIMENSION(:), POINTER :: eigocc
    INTEGER  :: handle,imo,istat

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_w_matrix_1"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   ---------------------------------------------------------------------------
    CALL timeset(routine_name,"I","",handle)
    NULLIFY(weighted_vectors)

    CALL set_matrix(w_matrix,0.0_wp)
    CALL cp_fm_create(weighted_vectors,mo_set%eigenvectors%matrix_struct,"weighted_vectors")
    CALL cp_fm_to_fm(mo_set%eigenvectors,weighted_vectors)

    ! scale every column with the occupation
    ALLOCATE(eigocc(mo_set%homo),stat=istat)
    IF (istat /= 0) CALL stop_memory(routine,"eigocc")
    DO imo=1,mo_set%homo
       eigocc(imo) = mo_set%eigenvalues(imo)*mo_set%occupation_numbers(imo)
    ENDDO
    CALL cp_fm_column_scale(weighted_vectors,eigocc)
    DEALLOCATE(eigocc)

    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%eigenvectors,&
                            matrix_g=weighted_vectors,&
                            first_col=1,&
                            last_col=mo_set%homo)

    CALL cp_fm_release(weighted_vectors)
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_w_matrix_1

! *****************************************************************************

  SUBROUTINE calculate_w_matrix_2(mo_set,ks_matrix,w_matrix)

!   Purpose: Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!            and the MO occupation numbers.

!   History: - Creation (17.07.2002,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER           :: w_matrix,ks_matrix
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(cp_fm_type), POINTER        :: weighted_vectors,h_block
    TYPE(cp_fm_struct_type), POINTER          :: fm_struct_tmp

!   *** Local variables ***

    INTEGER :: handle,imo,ncol_global,nrow_global,nrow_block,ncol_block

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_w_matrix_2","I","",handle)

    NULLIFY(weighted_vectors,h_block,fm_struct_tmp)

    CALL cp_fm_get_info(matrix=mo_set%eigenvectors,&
                        ncol_global=ncol_global,&
                        nrow_global=nrow_global,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block)

    CALL cp_fm_create(weighted_vectors,mo_set%eigenvectors%matrix_struct,"weighted_vectors")
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                            para_env=mo_set%eigenvectors%matrix_struct%para_env, &
                                            context=mo_set%eigenvectors%matrix_struct%context)
    CALL cp_fm_create(h_block,fm_struct_tmp, name="h block")
    CALL cp_fm_struct_release(fm_struct_tmp)


    IF (mo_set%homo /= ncol_global) THEN
      CALL stop_program("calculate_w_matrix_2","NYI")
    END IF

    CALL cp_sm_fm_multiply(ks_matrix,mo_set%eigenvectors,weighted_vectors, &
                           mo_set%homo,mo_set%eigenvectors%matrix_struct%para_env)
    CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_wp, &
                    mo_set%eigenvectors,weighted_vectors,0.0_wp,h_block)
    CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_wp, &
                    mo_set%eigenvectors,h_block,0.0_wp,weighted_vectors)

    ! scale every column with the occupation
    CALL cp_fm_column_scale(weighted_vectors,mo_set%occupation_numbers)

    CALL set_matrix(w_matrix,0.0_wp)
    CALL cp_sm_plus_fm_fm_t(sparse_matrix=w_matrix,&
                            matrix_v=mo_set%eigenvectors,&
                            matrix_g=weighted_vectors,&
                            first_col=1,&
                            last_col=mo_set%homo)

    CALL cp_fm_release(weighted_vectors)
    CALL cp_fm_release(h_block)
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_w_matrix_2

! *****************************************************************************

  SUBROUTINE subspace_eigenvalues_ks(orbitals,ks_matrix,evals_arg,para_env,ctxt,ionode,scr)

!   Purpose: Calculates the eigenvalues (ritz values) for a set of orbitals
!            and rotates them to make them eigenvectors
!            if evals_arg is present its filled
!            if scr and ionode are present they go to the screen as well.

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER           :: ks_matrix
    TYPE(cp_fm_type), POINTER          :: orbitals
    TYPE(cp_para_env_type),  POINTER   :: para_env
    TYPE(cp_blacs_env_type), POINTER   :: ctxt
    REAL(wp), DIMENSION(:),  POINTER, OPTIONAL :: evals_arg
    LOGICAL, OPTIONAL :: ionode
    INTEGER, OPTIONAL :: scr

!   *** Local variables ***

    TYPE(cp_fm_type), POINTER          :: weighted_vectors,h_block,e_vectors
    REAL(wp), DIMENSION(:), ALLOCATABLE       :: evals
    INTEGER  :: handle,imo,ncol_global,nrow_global,nrow_block,ncol_block,i,j
    TYPE(cp_fm_struct_type), POINTER :: fm_struct_tmp

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_subspace_eigenvalues","I","",handle)

    NULLIFY(weighted_vectors,h_block,e_vectors,fm_struct_tmp)
    CALL cp_fm_get_info(matrix=orbitals, &
                                 ncol_global=ncol_global, &
                                 nrow_global=nrow_global, &
                                 nrow_block=nrow_block, &
                                 ncol_block=ncol_block)

    ALLOCATE(evals(ncol_global))

    CALL cp_fm_create(weighted_vectors,orbitals%matrix_struct,"weighted_vectors")
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global,  ncol_global=ncol_global, &
                                            para_env=orbitals%matrix_struct%para_env, &
                                            context=orbitals%matrix_struct%context)
    CALL cp_fm_create(h_block,fm_struct_tmp, name="h block")
    CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e vectors")
    CALL cp_fm_struct_release(fm_struct_tmp)

    ! h subblock and diag
    CALL cp_sm_fm_multiply(ks_matrix,orbitals,weighted_vectors, &
                                       ncol_global,para_env=para_env)


    CALL cp_fm_gemm('T','N',ncol_global,ncol_global,nrow_global,1.0_wp, &
                    orbitals,weighted_vectors,0.0_wp,h_block)

    CALL cp_fm_syevd(h_block,e_vectors,evals)

    ! rotate the orbitals
    CALL cp_fm_gemm('N','N',nrow_global,ncol_global,ncol_global,1.0_wp, &
                    orbitals,e_vectors,0.0_wp,weighted_vectors)
    CALL cp_fm_to_fm(weighted_vectors,orbitals)

    ! give output
    IF (PRESENT(evals_arg)) THEN
      evals_arg(:)=evals(:)
    ENDIF

    IF (PRESENT(ionode) .OR. PRESENT(scr)) THEN
      IF (.NOT. PRESENT(ionode)) CALL stop_program("subspace evals","IONODE?")
      IF (.NOT. PRESENT(scr)) CALL stop_program("subspace evals","SCR?")
      IF (ionode) THEN
         DO i=1,ncol_global,4
            j=MIN(3,ncol_global-i)
            SELECT CASE (j)
            CASE(3)
             WRITE(scr,'(1X,4F12.8)') evals(i:i+j)
            CASE(2)
             WRITE(scr,'(1X,3F12.8)') evals(i:i+j)
            CASE(1)
             WRITE(scr,'(1X,2F12.8)') evals(i:i+j)
            CASE(0)
             WRITE(scr,'(1X,1F12.8)') evals(i:i+j)
            END SELECT
         ENDDO
      ENDIF
    ENDIF

    CALL cp_fm_release(weighted_vectors)
    CALL cp_fm_release(h_block)
    CALL cp_fm_release(e_vectors)

    DEALLOCATE(evals)
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE subspace_eigenvalues_ks

! *****************************************************************************
! computes the effective orthonormality of a set of mos given an s-matrix
! orthonormality is the max deviation from unity of the C^T S C
! *****************************************************************************
  SUBROUTINE calculate_orthonormality(mo_array,matrix_s,orthonormality)
   TYPE(real_matrix_type), POINTER                 :: matrix_s
   TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: mo_array
   REAL(wp)                                        :: orthonormality
   
   INTEGER :: ispin,nspin,i,j,n,k,handle
   REAL(wp) :: alpha,max_alpha
   TYPE(cp_fm_struct_type), POINTER :: tmp_fm_struct
   TYPE(cp_fm_type), POINTER :: svec,overlap
   INTEGER :: nrow_local,ncol_local
   INTEGER, DIMENSION(:), POINTER :: row_indices,col_indices
   NULLIFY(tmp_fm_struct,svec,overlap)

   CALL timeset("calculate_orthonormality","I","",handle)

   nspin=size(mo_array)
   max_alpha=0.0_wp

   DO ispin=1,nspin
          ! get S*C
          CALL cp_fm_create(svec,mo_array(ispin)%mo_set%eigenvectors%matrix_struct)       
          CALL cp_fm_get_info(mo_array(ispin)%mo_set%eigenvectors,nrow_global=n,ncol_global=k)
          CALL cp_sm_fm_multiply(matrix_s,mo_array(ispin)%mo_set%eigenvectors, svec, &
                                      k, para_env=mo_array(ispin)%mo_set%eigenvectors%matrix_struct%para_env)
          ! get C^T (S*C)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
                                   para_env=mo_array(ispin)%mo_set%eigenvectors%matrix_struct%para_env, &
                                   context=mo_array(ispin)%mo_set%eigenvectors%matrix_struct%context)
          CALL cp_fm_create(overlap,tmp_fm_struct)       
          CALL cp_fm_struct_release(tmp_fm_struct)
          CALL cp_fm_gemm('T','N',k,k,n,1.0_wp, mo_array(ispin)%mo_set%eigenvectors,svec,0.0_wp,overlap)
          CALL cp_fm_get_info(overlap,nrow_local=nrow_local,ncol_local=ncol_local, &
                              row_indices=row_indices,col_indices=col_indices)
          do i=1,nrow_local
           do j=1,ncol_local
              alpha=overlap%local_data(i,j)
              if (row_indices(i).eq.col_indices(j)) alpha=alpha-1.0_wp
              max_alpha=MAX(max_alpha,abs(alpha))
           enddo
          enddo
          CALL cp_fm_release(svec)
          CALL cp_fm_release(overlap)
    ENDDO
    CALL mp_max(max_alpha,mo_array(1)%mo_set%eigenvectors%matrix_struct%para_env%group)
    orthonormality=max_alpha
    ! write(6,*) "max deviation from orthonormalization ",orthonormality

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_orthonormality

! *****************************************************************************


  SUBROUTINE correct_mo_eigenvalues(mo_set,level_shift)

!   Purpose: Correct MO eigenvalues after MO level shifting.

!   History: - Creation (19.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER :: mo_set
    REAL(wp), INTENT(IN)       :: level_shift

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "correct_mo_eigenvalues"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: imo

!   ---------------------------------------------------------------------------

    IF (level_shift == 0.0_wp) RETURN

    DO imo=mo_set%homo+1,mo_set%nmo
      mo_set%eigenvalues(imo) = mo_set%eigenvalues(imo) - level_shift
    END DO

  END SUBROUTINE correct_mo_eigenvalues

! *****************************************************************************

  SUBROUTINE deallocate_mo_set(mo_set)

!   Purpose: Deallocate a wavefunction data structure.

!   History: - Creation (15.05.2001,MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER :: mo_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_mo_set"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(mo_set)) THEN
       IF (ASSOCIATED(mo_set%eigenvalues)) THEN
          DEALLOCATE (mo_set%eigenvalues,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"mo_set%eigenvalues")
       END IF

       IF (ASSOCIATED(mo_set%occupation_numbers)) THEN
          DEALLOCATE (mo_set%occupation_numbers,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"mo_set%occupation_numbers")
       END IF

       CALL cp_fm_release(mo_set%eigenvectors)

       DEALLOCATE (mo_set,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routine,"mo_set")
    END IF

  END SUBROUTINE deallocate_mo_set

! *****************************************************************************

  SUBROUTINE get_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,nmo,&
                        eigenvalues,occupation_numbers,eigenvectors)

!   Purpose: Get the components of a MO set data structure.

!   History: - Creation (22.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER                   :: mo_set
    REAL(wp), OPTIONAL, INTENT(OUT)              :: maxocc
    INTEGER, OPTIONAL, INTENT(OUT)               :: homo,lfomo,nao,nelectron,&
                                                    nmo
    REAL(wp), OPTIONAL, DIMENSION(:), POINTER    :: eigenvalues,&
                                                    occupation_numbers
    TYPE(cp_fm_type), OPTIONAL, POINTER :: eigenvectors

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

    IF (PRESENT(maxocc)) maxocc = mo_set%maxocc
    IF (PRESENT(homo)) homo = mo_set%homo
    IF (PRESENT(lfomo)) lfomo = mo_set%lfomo
    IF (PRESENT(nao)) nao = mo_set%nao
    IF (PRESENT(nelectron)) nelectron = mo_set%nelectron
    IF (PRESENT(nmo)) nmo = mo_set%nmo
    IF (PRESENT(eigenvalues)) eigenvalues => mo_set%eigenvalues
    IF (PRESENT(occupation_numbers)) THEN
      occupation_numbers => mo_set%occupation_numbers
    END IF
    IF (PRESENT(eigenvectors)) eigenvectors => mo_set%eigenvectors

  END SUBROUTINE get_mo_set

! *****************************************************************************

  SUBROUTINE set_mo_occupation(mo_set,smear)

!   Purpose: Smearing of the MO occupation due to the given smear value.

!   History: - Creation (17.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER :: mo_set
    REAL(wp), INTENT(IN)       :: smear

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_mo_occupation"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    REAL(wp) :: e1,e2,edelta,edist,nelec
    INTEGER  :: imo,nomo

!   ---------------------------------------------------------------------------

    nomo = NINT(mo_set%nelectron/mo_set%maxocc)

    mo_set%occupation_numbers(:) = 0.0_wp
    mo_set%occupation_numbers(1:nomo) = mo_set%maxocc

    mo_set%homo = nomo
    mo_set%lfomo = nomo + 1

    IF (smear /= 0.0_wp) THEN

!     *** Define the energy window for the eigenvalues ***

      e1 = mo_set%eigenvalues(nomo) - 0.5_wp*smear
      e2 = mo_set%eigenvalues(nomo) + 0.5_wp*smear

!     *** Find the lowest fractional occupied MO ***

      DO imo=1,nomo
        IF (mo_set%eigenvalues(imo) > e1) THEN
          mo_set%lfomo = imo
          EXIT
        END IF
      END DO

!     *** Find the highest fractional occupied MO which will be the HOMO ***

      DO imo=mo_set%nmo,nomo,-1
        IF (mo_set%eigenvalues(imo) < e2) THEN
          mo_set%homo = imo
          EXIT
        END IF
      END DO

!     *** Get the number of electrons to be smeared ***

      edist = 0.0_wp
      nelec = 0.0_wp

      DO imo=mo_set%lfomo,mo_set%homo
        nelec = nelec + mo_set%occupation_numbers(imo)
        edist = edist + ABS(e2 - mo_set%eigenvalues(imo))
      END DO

!     *** Smear electrons inside the window ***

      DO imo=mo_set%lfomo,mo_set%homo
        edelta = ABS(e2 - mo_set%eigenvalues(imo))
        mo_set%occupation_numbers(imo) = MIN(mo_set%maxocc,nelec*edelta/edist)
        nelec = nelec - mo_set%occupation_numbers(imo)
        edist = edist - edelta
      END DO

!     *** Check, if the smearing involves more than one MO ***

      IF (mo_set%lfomo == mo_set%homo) THEN
        mo_set%homo = nomo
        mo_set%lfomo = nomo + 1
      END IF

    END IF

  END SUBROUTINE set_mo_occupation

! *****************************************************************************

  SUBROUTINE write_mo_set_to_restart(mo_array,atomic_kind_set,particle_set,&
                                     globenv,id_nr)

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, INTENT(in)                           :: id_nr
    TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: mo_array

!   *** Local paramters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_mo_set_to_restart"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit
    INTEGER            :: i,iatom,icol,irow,iset,isgf,iso,ishell,ispin,istat,&
                          jcol,left,lshell,nao,natom,ncol,nmo,nset,nset_max,&
                          nshell_max,nsgf,nspin,output_unit,right

    INTEGER, DIMENSION(:), POINTER    :: nshell
    REAL(wp), DIMENSION(:,:), POINTER :: matrix
    INTEGER, DIMENSION(:,:), POINTER  :: l

    INTEGER, DIMENSION(:), POINTER     :: nset_info
    INTEGER, DIMENSION(:,:), POINTER   :: nshell_info
    INTEGER, DIMENSION(:,:,:), POINTER :: nso_info

    REAL(wp) , DIMENSION(:,:), POINTER :: vecbuffer

!   ---------------------------------------------------------------------------

    nspin = SIZE(mo_array)

    nao = mo_array(1)%mo_set%nao

    IF (globenv%ionode) THEN

      IF (id_nr==1) THEN
        file_name = globenv%restart_file_name
      ELSE
        file_name = TRIM(globenv%restart_file_name)//"-"//&
                    ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="WRITE",&
                     file_form="UNFORMATTED",&
                     file_status="REPLACE",&
                     unit_number=restart_unit)

!     *** create some info about the basis set first ***

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom=natom,&
                               nsgf=nsgf)

      nset_max = 0
      nshell_max = 0

      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             orb_basis_set=orb_basis_set)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=l)
        nset_max = MAX(nset_max,nset)
        DO iset=1,nset
          nshell_max = MAX(nshell_max,nshell(iset))
        END DO
      END DO

      ALLOCATE (nso_info(nshell_max,nset_max,natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info",&
                                       nshell_max*nset_max*natom*int_size)
      nso_info(:,:,:) = 0

      ALLOCATE (nshell_info(nset_max,natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info",&
                                       nset_max*natom*int_size)
      nshell_info(:,:) = 0

      ALLOCATE (nset_info(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info",natom*int_size)
      nset_info(:) = 0

      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             orb_basis_set=orb_basis_set)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=l)
        nset_info(iatom) = nset
        DO iset=1,nset
          nshell_info(iset,iatom)=nshell(iset)
          DO ishell=1,nshell(iset)
            lshell = l(ishell,iset)
            nso_info(ishell,iset,iatom) = nso(lshell)
          END DO
        END DO
      END DO

      WRITE (restart_unit) natom,nspin,nao,nset_max,nshell_max
      WRITE (restart_unit) nset_info
      WRITE (restart_unit) nshell_info
      WRITE (restart_unit) nso_info

      DEALLOCATE (nset_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info")

      DEALLOCATE (nshell_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info")

      DEALLOCATE (nso_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info")

    END IF

    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer",nao*wp_size)

    DO ispin=1,nspin

      nmo=mo_array(ispin)%mo_set%nmo

      IF (globenv%ionode) THEN
        WRITE (restart_unit) nmo,&
                             mo_array(ispin)%mo_set%homo,&
                             mo_array(ispin)%mo_set%lfomo,&
                             mo_array(ispin)%mo_set%nelectron
        WRITE (restart_unit) mo_array(ispin)%mo_set%eigenvalues(1:nmo),&
                             mo_array(ispin)%mo_set%occupation_numbers(1:nmo)
      END IF
      ! this and the corresponding in read should be blocked for efficiency
      DO i=1,nmo
        CALL cp_fm_get_submatrix(mo_array(ispin)%mo_set%eigenvectors,vecbuffer,1,i,nao,1,transpose=.TRUE.)
        IF (globenv%ionode) WRITE (restart_unit) vecbuffer
      END DO

    END DO

    DEALLOCATE (vecbuffer,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer")

    IF (globenv%ionode) CALL close_file(unit_number=restart_unit)

  END SUBROUTINE write_mo_set_to_restart

! *****************************************************************************

  SUBROUTINE read_mo_set_from_restart(mo_array,atomic_kind_set,particle_set,&
                                      globenv,id_nr)

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, INTENT(in)                           :: id_nr
    TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: mo_array

!   *** Local paramters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_mo_set_from_restart"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: mype,npe,restart_unit
    INTEGER            :: iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                          left,lshell,natom,ncol,nset,nsgf,output_unit,right,&
                          nao,nmo,i,nset_max,nshell_max,natom_read,nao_read,&
                          nmo_read,nspin_read,nspin,ispin,group,source

    INTEGER, DIMENSION(:), POINTER    :: nshell
    REAL(wp), DIMENSION(:,:), POINTER :: matrix
    INTEGER, DIMENSION(:,:), POINTER  :: l

    INTEGER, DIMENSION(:), POINTER     :: nset_info
    INTEGER, DIMENSION(:,:), POINTER   :: nshell_info
    INTEGER, DIMENSION(:,:,:), POINTER :: nso_info,offset_info

    REAL(wp) , DIMENSION(:,:), POINTER :: vecbuffer,vecbuffer_read

    INTEGER :: irow_read,nset_read,ishell_read,iset_read,offset_read
    LOGICAL :: use_this,did_change

!   ---------------------------------------------------------------------------

    group = globenv%group
    source = globenv%source

    nspin = SIZE(mo_array)
    nao = mo_array(1)%mo_set%nao

    IF (globenv%ionode) THEN

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom=natom,&
                               nsgf=nsgf)
      IF (id_nr==1) THEN
        file_name=globenv%restart_file_name
      ELSE
        file_name = TRIM(globenv%restart_file_name)//"-"//&
                         ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

      READ (restart_unit) natom_read,nspin_read,nao_read,nset_max,nshell_max

      ! we should allow for restarting with different spin settings
      IF (nspin_read /= nspin) THEN
        CALL stop_program(routine,"nspin is not equal (program this !!)")
      END IF

      IF (natom_read /= natom) THEN
        CALL stop_program(routine,"incorrect natom")
      END IF

      ALLOCATE (nso_info(nshell_max,nset_max,natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info",&
                                       nshell_max*nset_max*natom_read*int_size)

      ALLOCATE (nshell_info(nset_max,natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info",&
                                       nset_max*natom_read*int_size)

      ALLOCATE (nset_info(natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info",natom_read*int_size)


      ALLOCATE (offset_info(nshell_max,nset_max,natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info",&
                                       nshell_max*nset_max*natom_read*int_size)

      IF (nao_read /= nao) THEN
        WRITE(globenv%scr,*) " READ RESTART : WARNING : DIFFERENT # AOs ",nao,nao_read
        did_change=.TRUE.
      END IF

      READ (restart_unit) nset_info
      READ (restart_unit) nshell_info
      READ (restart_unit) nso_info

      i=1
      DO iatom=1,natom
        DO iset=1,nset_info(iatom)
          DO ishell=1,nshell_info(iset,iatom)
            offset_info(ishell,iset,iatom) = i
            i=i+nso_info(ishell,iset,iatom)
          END DO
        END DO
      END DO

      ALLOCATE(vecbuffer_read(1,nao_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vecbuffer_read",&
                                       nao_read*wp_size)

    END IF

    ALLOCATE (vecbuffer(1,nao),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer",nao*wp_size)

    DO ispin=1,nspin

      nmo=mo_array(ispin)%mo_set%nmo

      mo_array(ispin)%mo_set%eigenvalues(:) = 0.0_wp
      mo_array(ispin)%mo_set%occupation_numbers(:) = 0.0_wp
      CALL cp_fm_set_all(mo_array(ispin)%mo_set%eigenvectors,0.0_wp)

      IF (globenv%ionode) THEN
        READ (restart_unit) nmo_read,&
                            mo_array(ispin)%mo_set%homo,&
                            mo_array(ispin)%mo_set%lfomo,&
                            mo_array(ispin)%mo_set%nelectron
        IF (mo_array(ispin)%mo_set%homo > nmo) THEN
          CALL stop_program(routine,"The number of the occupied MOs on the "//&
                                    "restart unit is greater than the "//&
                                    "number of the allocated MOs.")
        END IF
        IF (nmo_read < nmo) THEN
          nmo = nmo_read
          CALL print_warning(routine,"The number of MOs on the restart "//&
                                     "unit is smaller than the number of "//&
                                     "the allocated MOs. The MO set will "//&
                                     "be padded with zeros!")
        ELSE IF (nmo_read > nmo) THEN
          CALL print_warning(routine,"The number of MOs on the restart "//&
                                     "unit is greater than the number of "//&
                                     "the allocated MOs. The read MO set "//&
                                     "will be truncated!")
        END IF
        READ (restart_unit) mo_array(ispin)%mo_set%eigenvalues(1:nmo),&
                            mo_array(ispin)%mo_set%occupation_numbers(1:nmo)
      END IF

      CALL mp_bcast(nmo,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%homo,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%lfomo,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%nelectron,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%eigenvalues,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%occupation_numbers,source,group)

      DO i=1,nmo
        IF (globenv%ionode) THEN
          READ (restart_unit) vecbuffer_read
          ! now, try to assign the read to the real vector
          ! in case the basis set changed this involves some guessing
          irow=1
          DO iatom=1,natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                                 orb_basis_set=orb_basis_set)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   nshell=nshell,&
                                   l=l)

            use_this = .TRUE.
            iset_read = 1
            DO iset=1,nset
              ishell_read = 1
              DO ishell=1,nshell(iset)
                lshell = l(ishell,iset)
                IF (iset_read > nset_info(iatom)) use_this = .FALSE.
                IF (use_this) THEN ! avoids out of bound access of the lower line if false
                   IF (nso(lshell) == nso_info(ishell_read,iset_read,iatom)) THEN
                      offset_read=offset_info(ishell_read,iset_read,iatom)
                      ishell_read=ishell_read+1
                      IF (ishell_read > nshell_info(iset,iatom)) THEN
                        ishell_read = 1
                        iset_read = iset_read+1
                      END IF
                    ELSE
                      use_this = .FALSE.
                    END IF
                END IF
                DO iso=1,nso(lshell)
                  IF (use_this) THEN
                    IF (offset_read-1+iso.LT.1 .OR. offset_read-1+iso.GT.nao_read) THEN
                      ! this indicates a bug
                      ! IF (did_change) WRITE(globenv%scr,*) "-?-",irow,offset_read-1+iso
                      vecbuffer(1,irow)=0.0_wp
                    ELSE
                      vecbuffer(1,irow)=vecbuffer_read(1,offset_read-1+iso)
                      ! IF (did_change) WRITE(globenv%scr,*) "   ",irow,offset_read-1+iso
                    END IF
                  ELSE
                    vecbuffer(1,irow) = 0.0_wp
                  END IF
                  irow = irow + 1
                END DO
                use_this = .TRUE.
              END DO
            END DO
          END DO

        ELSE

          vecbuffer(1,:) = 0.0_wp

        END IF

        CALL mp_bcast(vecbuffer,source,group)
        CALL cp_fm_set_submatrix(mo_array(ispin)%mo_set%eigenvectors,vecbuffer,1,i,nao,1,transpose=.TRUE.)
        did_change=.FALSE.

      END DO

      IF (globenv%print%level > LOW) THEN
        CALL write_mo_set(mo_array(ispin)%mo_set,atomic_kind_set,particle_set,&
                          4,6,globenv)
      END IF

    END DO

    DEALLOCATE(vecbuffer,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer")

    IF (globenv%ionode) THEN
      DEALLOCATE(vecbuffer_read,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vecbuffer_read")
      DEALLOCATE(offset_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"offset_info")
      DEALLOCATE(nso_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info")
      DEALLOCATE(nshell_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info")
      DEALLOCATE(nset_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info")
      CALL close_file(unit_number=restart_unit)
    END IF

  END SUBROUTINE read_mo_set_from_restart

! *****************************************************************************

  SUBROUTINE write_mo_set_to_output_unit(mo_set,atomic_kind_set,particle_set,&
                                         before,after,globenv)

!   Purpose: Write the MO eigenvalues, MO occupation numbers and
!            MO eigenvectors.

!   History: - Creation (15.05.2001,MK)

!   ***************************************************************************

!   List of variables:

!   after : Number of digits after point.
!   before: Number of digits before point.

!   ---------------------------------------------------------------------------

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(mo_set_type), POINTER                    :: mo_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    INTEGER, INTENT(IN)                           :: after,before

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_mo_set_to_output_unit (MODULE qs_mo_types)"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    CHARACTER(LEN=60) :: name
    CHARACTER(LEN=25) :: fmtstr1
    CHARACTER(LEN=22) :: fmtstr2
    CHARACTER(LEN=35) :: fmtstr3
    CHARACTER(LEN=2)  :: element_symbol
    INTEGER           :: from,iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                         left,lshell,natom,ncol,nset,nsgf,output_unit,right,&
                         to,width

    CHARACTER(LEN=6), DIMENSION(:), POINTER :: sgf_symbol
    INTEGER, DIMENSION(:), POINTER          :: nshell
    REAL(wp), DIMENSION(:,:), POINTER       :: matrix
    INTEGER, DIMENSION(:,:), POINTER        :: l
    INTEGER                                 :: nrow_global,ncol_global

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%mo_eigenvalues) RETURN

    output_unit = globenv%scr

    IF (globenv%print%mo_eigenvectors) THEN
      CALL cp_fm_get_info(mo_set%eigenvectors,nrow_global=nrow_global,ncol_global=ncol_global)
      ALLOCATE(matrix(nrow_global,ncol_global))
      CALL cp_fm_get_submatrix(mo_set%eigenvectors,matrix)
      IF (.NOT.globenv%ionode) DEALLOCATE(matrix)
    END IF

    IF (.NOT.globenv%ionode) RETURN

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,&
                             nsgf=nsgf)

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,23X,  (1X,F  .  ))"
    fmtstr3 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    IF (globenv%print%mo_eigenvectors) THEN
      name = "MO EIGENVALUES, MO OCCUPATION NUMBERS AND MO EIGENVECTORS"
    ELSE IF (globenv%print%mo_occupation_numbers) THEN
      name = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
    ELSE
      name = "MO EIGENVALUES"
    END IF

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
    WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

    WRITE (UNIT=fmtstr2(9:10),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr2(16:17),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr2(19:20),FMT="(I2)") after

    WRITE (UNIT=fmtstr3(22:23),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr3(29:30),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr3(32:33),FMT="(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol

      from = icol
      to = MIN((from+ncol-1),nsgf)

      WRITE (UNIT=output_unit,FMT=fmtstr1)&
        (jcol,jcol=from,to)
      WRITE (UNIT=output_unit,FMT=fmtstr2)&
        (mo_set%eigenvalues(jcol),jcol=from,to)
      WRITE (UNIT=output_unit,FMT="(A)") ""

      IF (globenv%print%mo_occupation_numbers) THEN
        WRITE (UNIT=output_unit,FMT=fmtstr2)&
          (mo_set%occupation_numbers(jcol),jcol=from,to)
        WRITE (UNIT=output_unit,FMT="(A)") ""
      END IF

      IF (globenv%print%mo_eigenvectors) THEN

        irow = 1

        DO iatom=1,natom

          IF (iatom /= 1) WRITE (UNIT=output_unit,FMT="(A)") ""

          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol,&
                               orb_basis_set=orb_basis_set)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 nset=nset,&
                                 nshell=nshell,&
                                 l=l,&
                                 sgf_symbol=sgf_symbol)

          isgf = 1

          DO iset=1,nset
            DO ishell=1,nshell(iset)
              lshell = l(ishell,iset)
              DO iso=1,nso(lshell)
                WRITE (UNIT=output_unit,FMT=fmtstr3)&
                  irow,iatom,element_symbol,sgf_symbol(isgf),&
                  (matrix(irow,jcol),jcol=from,to)
                isgf = isgf + 1
                irow = irow + 1
              END DO
            END DO
          END DO

        END DO

      END IF

    END DO

    WRITE (UNIT=output_unit,FMT="(/)")

!   *** Release work storage ***

    IF (globenv%print%mo_eigenvectors) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"matrix")
    END IF

  END SUBROUTINE write_mo_set_to_output_unit

! *****************************************************************************

END MODULE qs_mo_types
