!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a linear scaling quickstep SCF run based on the density
!>        matrix, with a focus on the interface between dm_ls_scf and qs
!> \par History
!>       2011.04 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf_qs
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_data, array_i1d_obj, array_new, array_release, &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_complete_redistribute, &
       cp_dbcsr_copy, cp_dbcsr_copy_into_existing, cp_dbcsr_create, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_finalize, &
       cp_dbcsr_init, cp_dbcsr_nblkrows_total, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_row_block_sizes, cp_dbcsr_set, &
       cp_dbcsr_type, dbcsr_distribution_col_dist, dbcsr_distribution_init, &
       dbcsr_distribution_mp, dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_distribution_row_dist, dbcsr_name, &
       dbcsr_type_no_symmetry, dbcsr_type_real_4, dbcsr_type_real_8, &
       dbcsr_type_symmetric
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE dm_ls_scf_types,                 ONLY: ls_cluster_atomic,&
                                             ls_cluster_molecular,&
                                             ls_mstruct_type,&
                                             ls_scf_env_type
  USE kinds,                           ONLY: dp
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_initial_guess,                ONLY: calculate_atomic_block_dm,&
                                             calculate_mopac_dm
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                             qs_ks_env_type,&
                                             set_ks_env
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_qs'

  PUBLIC :: matrix_ls_create, matrix_qs_to_ls, matrix_ls_to_qs, ls_scf_init_qs,&
           ls_scf_dm_to_ks, ls_scf_qs_atomic_guess, write_matrix_to_cube

CONTAINS

! *****************************************************************************
!> \brief create a matrix for use (and as a template) in ls based on a qs template
!> \param matrix_ls ...
!> \param matrix_qs ...
!> \param ls_mstruct ...
!> \param error ...
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
!> \param matrix_ls ...
!> \param matrix_qs ...
!> \param ls_mstruct ...
!> \param error ...
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

       ! complete_redistribute does not zero blocks that might be present in the
       ! target but not in the source
       CALL cp_dbcsr_set(matrix_ls,0.0_dp,error=error)

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
!> \param matrix_qs ...
!> \param matrix_ls ...
!> \param ls_mstruct ...
!> \param error ...
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
!> \param qs_env ...
!> \param error ...
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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    NULLIFY(sab_orb)
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
                           sab_orb=sab_orb,&
                           error=error)

    nspin=dft_control%nspins

    ! we might have to create matrix_ks
    IF (.NOT.ASSOCIATED(matrix_ks)) THEN
       CALL cp_dbcsr_allocate_matrix_set(matrix_ks,nspin,error)
       DO ispin=1,nspin
          ALLOCATE(matrix_ks(ispin)%matrix)
          CALL cp_dbcsr_init(matrix_ks(ispin)%matrix,error=error)
          CALL cp_dbcsr_create(matrix_ks(ispin)%matrix,template=matrix_s(1)%matrix,error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(ispin)%matrix,sab_orb,error=error)
          CALL cp_dbcsr_set(matrix_ks(ispin)%matrix,0.0_dp,error=error)
       ENDDO
       CALL set_ks_env(ks_env,matrix_ks=matrix_ks,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE ls_scf_init_qs

! *****************************************************************************
!> \brief get an atomic initial guess
!> \param qs_env ...
!> \param energy ...
!> \param error ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_qs_atomic_guess(qs_env,energy,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp)                            :: energy
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_qs_atomic_guess', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin, unit_nr
    INTEGER, DIMENSION(2)                    :: nelectron_spin
    LOGICAL                                  :: has_unit_metric
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)
    NULLIFY(rho,rho_ao)

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
                           energy=qs_energy,&
                           atomic_kind_set=atomic_kind_set,&
                           qs_kind_set=qs_kind_set,&
                           particle_set=particle_set,&
                           has_unit_metric=has_unit_metric,&
                           para_env=para_env,&
                           nelectron_spin=nelectron_spin,&
                           rho=rho,&
                           error=error)

    CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)

    nspin=dft_control%nspins

    ! create an initial atomic guess
    DO ispin=1,nspin
       CALL cp_dbcsr_set(rho_ao(ispin)%matrix,0.0_dp,error=error)
       IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
          CALL calculate_mopac_dm(rho_ao(ispin)%matrix,matrix_s(1)%matrix, has_unit_metric, &
                                  dft_control,particle_set, atomic_kind_set, qs_kind_set,&
                                  ispin, nspin, nelectron_spin(ispin),&
                                  para_env, error)
       ELSE
          CALL calculate_atomic_block_dm(rho_ao(ispin)%matrix,matrix_s(1)%matrix, &
                                         particle_set, atomic_kind_set, qs_kind_set, &
                                         ispin, nspin, nelectron_spin(ispin), &
                                         unit_nr, error)
       ENDIF
    ENDDO

    CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., error=error)

    energy=qs_energy%total

    CALL timestop(handle)

  END SUBROUTINE ls_scf_qs_atomic_guess

! *****************************************************************************
!> \brief use the density matrix in ls_scf_env to compute the new energy and KS matrix
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param energy_new ...
!> \param error ...
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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho

    NULLIFY(energy, rho, rho_ao)
    CALL timeset(routineN,handle)

    nspin=ls_scf_env%nspins
    CALL get_qs_env(qs_env, energy=energy, rho=rho, error=error)
    CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)

    ! set the new density matrix
    DO ispin=1,nspin
         CALL matrix_ls_to_qs(rho_ao(ispin)%matrix, ls_scf_env%matrix_p(ispin), ls_scf_env%ls_mstruct, error=error)
    END DO

    ! compute the corresponding KS matrix and new energy
    CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)
    CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE., error=error)
    CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE.,&
                             just_energy=.FALSE., print_active=.TRUE., error=error)
    energy_new=energy%total

    CALL timestop(handle)

  END SUBROUTINE ls_scf_dm_to_ks


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_p_ls ...
!> \param unit_nr ...
!> \param title ...
!> \param stride ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_matrix_to_cube(qs_env,ls_scf_env,matrix_p_ls,unit_nr,title,stride,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_p_ls
    INTEGER, INTENT(IN)                      :: unit_nr
    CHARACTER(LEN=*), INTENT(IN)             :: title
    INTEGER, DIMENSION(:), POINTER           :: stride
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_matrix_to_cube', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: tot_rho
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type), POINTER             :: matrix_p_qs
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_subsys_type), POINTER            :: subsys

    CALL timeset(routineN,handle)

    NULLIFY(ks_env, pw_env,auxbas_pw_pool, pw_pools, particles, subsys, matrix_ks)

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    subsys=subsys,&
                    pw_env=pw_env,&
                    matrix_ks=matrix_ks,&
                    error=error)

    CALL qs_subsys_get(subsys,particles=particles,error=error)

    ! convert the density matrix (ls style) to QS style
    ALLOCATE(matrix_p_qs)
    CALL cp_dbcsr_init(matrix_p_qs,error=error) !initialization
    CALL cp_dbcsr_copy(matrix_p_qs,matrix_ks(1)%matrix,error=error)
    CALL cp_dbcsr_set(matrix_p_qs,0.0_dp,error=error) !zero matrix creation   
    CALL matrix_ls_to_qs(matrix_p_qs, matrix_p_ls,ls_scf_env%ls_mstruct,error=error) 

    ! Print total electronic density
    CALL pw_env_get(pw_env=pw_env,&
                    auxbas_pw_pool=auxbas_pw_pool,&
                    pw_pools=pw_pools,&
                    error=error)
    CALL pw_pool_create_pw(pool=auxbas_pw_pool,&
                          pw=wf_r%pw,&
                          use_data=REALDATA3D,&
                          in_space=REALSPACE,&
                          error=error)
    CALL pw_zero(wf_r%pw,error=error)
    CALL pw_pool_create_pw(pool=auxbas_pw_pool,&
                          pw=wf_g%pw,&
                          use_data=COMPLEXDATA1D,&
                          in_space=RECIPROCALSPACE,&
                          error=error)
    CALL pw_zero(wf_g%pw,error=error)
    CALL calculate_rho_elec(matrix_p=matrix_p_qs,&
                            rho=wf_r,&
                            rho_gspace=wf_g,&
                            total_rho=tot_rho,&
                            ks_env=ks_env, &
                            error=error)

    ! write this to a cube
    CALL cp_pw_to_cube(wf_r%pw, unit_nr=unit_nr, title=title,&
         particles=particles, stride=stride, error = error)

    !free memory
    CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw,error=error) 
    CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_g%pw,error=error)   
    CALL cp_dbcsr_release(matrix_p_qs,error=error)
    DEALLOCATE(matrix_p_qs)

    CALL timestop(handle)

  END SUBROUTINE write_matrix_to_cube

END MODULE dm_ls_scf_qs
