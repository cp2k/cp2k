!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a linear scaling quickstep SCF run based on the density
!>        matrix, with a focus on the interface between dm_ls_scf and qs
!> \par History
!>       2011.04 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf_qs
  USE array_types,                     ONLY: array_data,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_release
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_complete_redistribute, cp_dbcsr_copy, &
       cp_dbcsr_copy_into_existing, cp_dbcsr_create, cp_dbcsr_desymmetrize, &
       cp_dbcsr_distribution, cp_dbcsr_finalize, cp_dbcsr_init, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_alloc_block_from_nbl,&
                                             cp_dbcsr_allocate_matrix_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_col_dist,&
                                             dbcsr_distribution_init,&
                                             dbcsr_distribution_mp,&
                                             dbcsr_distribution_new,&
                                             dbcsr_distribution_release,&
                                             dbcsr_distribution_row_dist,&
                                             dbcsr_name
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8,&
                                             dbcsr_type_symmetric
  USE dm_ls_scf_types,                 ONLY: ls_cluster_atomic,&
                                             ls_cluster_molecular,&
                                             ls_mstruct_type,&
                                             ls_scf_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_initial_guess,                ONLY: calculate_atomic_block_dm,&
                                             calculate_mopac_dm
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_qs'

  PUBLIC :: matrix_ls_create, matrix_qs_to_ls, matrix_ls_to_qs, ls_scf_init_qs, ls_scf_dm_to_ks

CONTAINS

! *****************************************************************************
!> \brief create a matrix for use (and as a template) in ls based on a qs template
!> \par History
!>       2011.03 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE matrix_ls_create(matrix_ls,matrix_qs,ls_mstruct,error)
    TYPE(cp_dbcsr_type)                      :: matrix_ls, matrix_qs
    TYPE(ls_mstruct_type), INTENT(IN)        :: ls_mstruct
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_ls_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_type, handle, iatom, &
                                                imol, jatom, natom, nmol
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_to_cluster, &
                                                atom_to_cluster_primus, &
                                                cl_dist, cluster_sizes, &
                                                primus_of_mol
    INTEGER, DIMENSION(:), POINTER           :: blk_dist, blk_sizes
    TYPE(array_i1d_obj)                      :: clustered_blk_sizes, &
                                                clustered_col_dist, &
                                                clustered_row_dist
    TYPE(dbcsr_distribution_obj)             :: block_distribution, &
                                                clustered_distribution

    CALL timeset(routineN,handle)

    IF (ls_mstruct%single_precision) THEN
       data_type=dbcsr_type_real_4
    ELSE
       data_type=dbcsr_type_real_8
    ENDIF

    ! here we allocate our basic matrices, right now, use QS S as a template
    ! later, we might want to use our own format (e.g. clustered)
    SELECT CASE(ls_mstruct%cluster_type)
    CASE(ls_cluster_atomic)
       CALL cp_dbcsr_init(matrix_ls,error=error)
       CALL cp_dbcsr_create(matrix_ls,template=matrix_qs, &
                            data_type=data_type,&
                            error=error)
    CASE(ls_cluster_molecular)
       ! create format of the clustered matrix
       natom=cp_dbcsr_nblkrows_total(matrix_qs)
       nmol=MAXVAL(ls_mstruct%atom_to_molecule) 
       ALLOCATE(atom_to_cluster_primus(natom))
       ALLOCATE(atom_to_cluster(natom))
       ALLOCATE(primus_of_mol(nmol))
       DO iatom=1,natom
          atom_to_cluster(iatom)=ls_mstruct%atom_to_molecule(iatom)
          ! the first atom of the molecule is the primus
          ! if the number of atoms per molecule is independent of system size, this is not a quadratic loop
          ! it assumes that all atoms of the molecule are consecutive.
          DO jatom=iatom,1,-1
             IF (ls_mstruct%atom_to_molecule(jatom)==atom_to_cluster(iatom)) THEN
                atom_to_cluster_primus(iatom)=jatom
             ELSE
                EXIT
             ENDIF
          ENDDO
          primus_of_mol(atom_to_cluster(iatom))=atom_to_cluster_primus(iatom)
       ENDDO

       block_distribution = cp_dbcsr_distribution (matrix_qs)
       CALL dbcsr_distribution_init (clustered_distribution)
       ALLOCATE(cl_dist(nmol))
       ! row
       blk_dist => array_data(dbcsr_distribution_row_dist(block_distribution))
       DO imol = 1, nmol
         cl_dist(imol) = blk_dist(primus_of_mol(imol))
       ENDDO
       CALL array_new (clustered_row_dist, cl_dist, lb=1)
       ! col
       blk_dist => array_data(dbcsr_distribution_col_dist(block_distribution))
       DO imol = 1, nmol
         cl_dist(imol) = blk_dist(primus_of_mol(imol))
       ENDDO
       CALL array_new (clustered_col_dist, cl_dist, lb=1)

       blk_sizes => array_data(cp_dbcsr_row_block_sizes(matrix_qs))
       ALLOCATE(cluster_sizes(nmol))
       cluster_sizes=0
       DO iatom = 1, natom
          cluster_sizes(atom_to_cluster(iatom))=cluster_sizes(atom_to_cluster(iatom))+&
             blk_sizes(iatom)
       ENDDO

       ! the parallel distribution
       CALL array_new(clustered_blk_sizes, cluster_sizes, lb=1)
       CALL dbcsr_distribution_new(clustered_distribution,&
            dbcsr_distribution_mp(block_distribution),&
            clustered_row_dist, clustered_col_dist)
       CALL array_release (clustered_row_dist)
       CALL array_release (clustered_col_dist)

       ! the matrix
       CALL cp_dbcsr_init (matrix_ls, error=error)
       CALL cp_dbcsr_create (matrix_ls, dbcsr_name (matrix_qs%matrix),&
            clustered_distribution,&
            dbcsr_type_symmetric,&
            clustered_blk_sizes, clustered_blk_sizes,&
            data_type=data_type,&
            error=error)
       CALL array_release (clustered_blk_sizes)
       CALL dbcsr_distribution_release (clustered_distribution)
       CALL cp_dbcsr_finalize (matrix_ls, error=error)
    CASE DEFAULT
       STOP " BUG cluster option"
    END SELECT


    CALL timestop (handle)

  END SUBROUTINE matrix_ls_create

! *****************************************************************************
!> \brief first link to QS, copy a QS matrix to LS matrix
!>        used to isolate QS style matrices from LS style
!>        will be useful for future features (e.g. precision, symmetry, blocking, ...)
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE matrix_qs_to_ls(matrix_ls,matrix_qs,ls_mstruct,error)
    TYPE(cp_dbcsr_type)                      :: matrix_ls, matrix_qs
    TYPE(ls_mstruct_type), INTENT(IN)        :: ls_mstruct
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_qs_to_ls', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: matrix_qs_nosym, matrix_tmp

    CALL timeset(routineN,handle)

    SELECT CASE(ls_mstruct%cluster_type)
    CASE(ls_cluster_atomic)
       ! takes care of an eventual data_type conversion
       CALL cp_dbcsr_copy(matrix_ls,matrix_qs,error=error)
    CASE(ls_cluster_molecular)
       ! desymmetrize the qs matrix
       CALL cp_dbcsr_init (matrix_qs_nosym, error=error)
       CALL cp_dbcsr_create (matrix_qs_nosym, template=matrix_qs,&
            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL cp_dbcsr_desymmetrize (matrix_qs, matrix_qs_nosym,&
            error=error)

       ! perform the magic complete redistribute copy
       IF (ls_mstruct%single_precision) THEN
          CALL cp_dbcsr_init(matrix_tmp, error=error)
          CALL cp_dbcsr_create(matrix_tmp,template=matrix_qs_nosym, &
                               data_type=dbcsr_type_real_4, error=error)
          CALL cp_dbcsr_copy(matrix_tmp,matrix_qs_nosym,error=error)
          CALL cp_dbcsr_complete_redistribute(matrix_tmp, matrix_ls, error=error);
          CALL cp_dbcsr_release(matrix_tmp, error=error)
       ELSE
          CALL cp_dbcsr_complete_redistribute(matrix_qs_nosym, matrix_ls, error=error);
       ENDIF
       CALL cp_dbcsr_release (matrix_qs_nosym, error=error)

    CASE DEFAULT 
       STOP
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief second link to QS, copy a LS matrix to QS matrix
!>        used to isolate QS style matrices from LS style
!>        will be useful for future features (e.g. precision, symmetry, blocking, ...)
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE matrix_ls_to_qs(matrix_qs,matrix_ls,ls_mstruct,error)
    TYPE(cp_dbcsr_type)                      :: matrix_qs, matrix_ls
    TYPE(ls_mstruct_type), INTENT(IN)        :: ls_mstruct
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_ls_to_qs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(cp_dbcsr_type)                      :: matrix_tmp

    CALL timeset(routineN,handle)

    IF (ls_mstruct%single_precision) THEN
       CALL cp_dbcsr_init (matrix_tmp, error=error)
       CALL cp_dbcsr_create (matrix_tmp, template=matrix_ls,&
            data_type=dbcsr_type_real_8, error=error)
       CALL cp_dbcsr_copy (matrix_tmp, matrix_ls, error=error)
    ENDIF


    SELECT CASE(ls_mstruct%cluster_type)
    CASE(ls_cluster_atomic)
       IF (ls_mstruct%single_precision) THEN
          CALL cp_dbcsr_copy_into_existing (matrix_qs, matrix_tmp, error=error)
       ELSE
          CALL cp_dbcsr_copy_into_existing (matrix_qs, matrix_ls, error=error)
       ENDIF
    CASE(ls_cluster_molecular)
       CALL cp_dbcsr_set(matrix_qs,0.0_dp,error=error)
       IF (ls_mstruct%single_precision) THEN
          CALL cp_dbcsr_complete_redistribute(matrix_tmp, matrix_qs, keep_sparsity=.TRUE., error=error) 
       ELSE
          CALL cp_dbcsr_complete_redistribute(matrix_ls, matrix_qs, keep_sparsity=.TRUE., error=error) 
       ENDIF
    CASE DEFAULT 
       STOP "BUG"
    END SELECT

    IF (ls_mstruct%single_precision) THEN
       CALL cp_dbcsr_release(matrix_tmp,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE matrix_ls_to_qs

! *****************************************************************************
!> \brief further required initialization of QS.
!>        Might be factored-out since this seems common code with the other SCF.
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_init_qs(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_init_qs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, unit_nr
    INTEGER, DIMENSION(2)                    :: nelectron_spin
    LOGICAL                                  :: has_unit_metric
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! get basic quantities from the qs_env
    CALL get_qs_env(qs_env,dft_control=dft_control,&
                           matrix_s=matrix_s,&
                           matrix_ks=matrix_ks,&
                           ks_env=ks_env,&
                           input=input,&
                           energy=energy,&
                           atomic_kind_set=atomic_kind_set,&
                           particle_set=particle_set,&
                           has_unit_metric=has_unit_metric,&
                           para_env=para_env,&
                           nelectron_spin=nelectron_spin,error=error)

    nspin=dft_control%nspins


    ! we might have to create matrix_ks
    IF (.NOT.ASSOCIATED(matrix_ks)) THEN
       CALL cp_dbcsr_allocate_matrix_set(matrix_ks,nspin,error)
       DO ispin=1,nspin
          ALLOCATE(matrix_ks(ispin)%matrix)
          CALL cp_dbcsr_init(matrix_ks(ispin)%matrix,error=error)
          CALL cp_dbcsr_create(matrix_ks(ispin)%matrix,template=matrix_s(1)%matrix,error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(ispin)%matrix,qs_env%sab_orb,error=error)
          CALL cp_dbcsr_set(matrix_ks(ispin)%matrix,0.0_dp,error=error)
       ENDDO
       CALL set_qs_env(qs_env=qs_env,matrix_ks=matrix_ks,error=error)
    ENDIF
    ! and ks_env
    IF (.NOT.ASSOCIATED(ks_env)) THEN
       CALL qs_ks_create(ks_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env, ks_env=ks_env,error=error)
       CALL qs_ks_release(ks_env,error=error)
    END IF

    ! create an initial atomic guess
    DO ispin=1,nspin
       CALL cp_dbcsr_set(qs_env%rho%rho_ao(ispin)%matrix,0.0_dp,error=error)
       IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
          CALL calculate_mopac_dm(qs_env%rho%rho_ao(ispin)%matrix,matrix_s(1)%matrix, has_unit_metric, &
                                  dft_control,particle_set, atomic_kind_set,&
                                  ispin, nspin, nelectron_spin(ispin),&
                                  para_env, error)
       ELSE
          CALL calculate_atomic_block_dm(qs_env%rho%rho_ao(ispin)%matrix,matrix_s(1)%matrix, &
                                         particle_set, atomic_kind_set, &
                                         ispin, nspin, nelectron_spin(ispin), &
                                         unit_nr, error)
       ENDIF
    ENDDO

    CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
                             just_energy=.FALSE.,error=error)

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'(T2,A,F20.9)') "Energy with the initial guess:",energy%total
       WRITE(unit_nr,'()') 
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init_qs

! *****************************************************************************
!> \brief use the density matrix in ls_scf_env to compute the new energy and KS matrix
!> \par History
!>       2011.04 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_dm_to_ks(qs_env,ls_scf_env,energy_new,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    REAL(KIND=dp)                            :: energy_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_dm_to_ks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    TYPE(qs_energy_type), POINTER            :: energy

    CALL timeset(routineN,handle)

    nspin=ls_scf_env%nspins
    CALL get_qs_env(qs_env, energy=energy, error=error)

    ! set the new density matrix 
    DO ispin=1,nspin
         CALL matrix_ls_to_qs(qs_env%rho%rho_ao(ispin)%matrix, ls_scf_env%matrix_p(ispin), ls_scf_env%ls_mstruct, error=error)
    END DO

    ! compute the corresponding KS matrix and new energy
    CALL qs_rho_update_rho(qs_env%rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
                             just_energy=.FALSE., print_active=.TRUE., error=error)
    energy_new=energy%total

    CALL timestop(handle)

  END SUBROUTINE ls_scf_dm_to_ks


END MODULE dm_ls_scf_qs
