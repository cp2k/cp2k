!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines needed for kpoint calculation
!> \par History
!>       2014.07 created [JGH]
!>       2014.11 unified k-point and gamma-point code [Ole Schuett]
!> \author JGH
! *****************************************************************************
MODULE kpoint_methods
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_get_block_p, &
       cp_dbcsr_init, cp_dbcsr_p_type, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_type_antisymmetric, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                             cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_env,                     ONLY: cp_cart_create,&
                                             cp_para_env_create,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_cart_type,&
                                             cp_para_env_type
  USE fermi_utils,                     ONLY: fermikp
  USE k290,                            ONLY: csym_type,&
                                             kp_sym_gen,&
                                             release_csym_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_env_create,&
                                             kpoint_env_p_type,&
                                             kpoint_env_type,&
                                             kpoint_type
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_max,&
                                             mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_matrix_pools,                 ONLY: mpools_create,&
                                             mpools_get,&
                                             mpools_rebuild_fm_pools,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             get_neighbor_list_set_p,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE scf_control_types,               ONLY: smear_type
  USE termination,                     ONLY: print_warning
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_methods'

  PUBLIC :: kpoint_initialize, kpoint_env_initialize, kpoint_initialize_mos
  PUBLIC :: kpoint_init_cell_index, kpoint_set_mo_occupation
  PUBLIC :: kpoint_density_matrices, kpoint_density_transform
  PUBLIC :: rskp_transform

! *****************************************************************************  

CONTAINS

! *****************************************************************************
!> \brief Generate the kpoints and initialize the kpoint environment
!> \param kpoint       The kpoint environment
!> \param particle_set Particle types and coordinates
!> \param cell         Computational cell information
!> \param error        CP2K error handling
! *****************************************************************************
  SUBROUTINE kpoint_initialize(kpoint, particle_set, cell, error)

    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_initialize', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2), ALLOCATABLE, &
      DIMENSION(:)                           :: element
    INTEGER                                  :: i, ierr, ik, iounit, natom
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atype
    LOGICAL                                  :: failure, inv, spez, symmo
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: weight, zeff
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: coord
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(csym_type), POINTER                 :: crys_sym

    failure=.FALSE.

    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)

    SELECT CASE (kpoint%kp_scheme)
       CASE ("NONE")
          ! do nothing
       CASE ("GAMMA")
          kpoint%nkp = 1
          ALLOCATE(kpoint%xkp(3,1),kpoint%wkp(1),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          kpoint%xkp(1:3,1) = 0.0_dp
          kpoint%wkp(1) = 1.0_dp
       CASE ("MONKHORST-PACK","MACDONALD")
          NULLIFY(crys_sym)

          IF(kpoint%full_grid) THEN
             ! we set up a random molecule to avoid any possible symmetry
             natom = 10
             ALLOCATE(coord(3,natom),zeff(natom),weight(natom),atype(natom),element(natom),STAT=ierr)
             CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
             DO i=1,natom
                atype(i) = i
                coord(1,i) = SIN(i*0.12345_dp)
                coord(2,i) = COS(i*0.23456_dp)
                coord(3,i) = SIN(i*0.34567_dp)
             END DO
          ELSE
             natom = SIZE(particle_set)
             ALLOCATE(coord(3,natom),zeff(natom),weight(natom),atype(natom),element(natom),STAT=ierr)
             CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
             DO i=1,natom
                CALL get_atomic_kind(atomic_kind=particle_set(i)%atomic_kind,kind_number=atype(i))
                coord(1:3,i) = particle_set(i)%r(1:3)
             END DO
          END IF
          IF(kpoint%verbose) THEN
             NULLIFY(logger)
             logger => cp_error_get_logger(error)
             iounit = cp_logger_get_default_io_unit(logger)
          ELSE
             iounit = -1
          END IF

          CALL kp_sym_gen(crys_sym,coord,atype,cell%hmat,&
                          nk=kpoint%nkp_grid,&
                          shift=kpoint%kp_shift,&
                          symm=kpoint%symmetry,&
                          delta=kpoint%eps_geo,&
                          iounit=iounit,&
                          xkpoint=kpoint%xkp,&
                          wkpoint=kpoint%wkp, error=error)
          kpoint%nkp = SIZE(kpoint%wkp)

          ! inversion symmetry
          inv = (crys_sym%li > 0)
          ! symmorphic group without special translation
          symmo = (crys_sym%isy ==1)

          CALL release_csym_type(crys_sym,error)
          DEALLOCATE(coord,zeff,weight,atype,element,STAT=ierr)
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

       CASE ("GENERAL")
          ! do nothing
       CASE DEFAULT
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT

    ! check for consistency of options
    SELECT CASE (kpoint%kp_scheme)
       CASE ("NONE")
          ! don't use k-point code
       CASE ("GAMMA")
          CPPostcondition(kpoint%nkp==1, cp_failure_level, routineP, error, failure)
          CPPostcondition(SUM(ABS(kpoint%xkp))<=1.e-12_dp, cp_failure_level, routineP, error, failure)
          CPPostcondition(kpoint%wkp(1)==1.0_dp, cp_failure_level, routineP, error, failure)
          CPPostcondition(.NOT.kpoint%symmetry, cp_failure_level, routineP, error, failure)
       CASE ("GENERAL")
          CPPostcondition(.NOT.kpoint%symmetry, cp_failure_level, routineP, error, failure)
          CPPostcondition(kpoint%nkp>=1, cp_failure_level, routineP, error, failure)
       CASE ("MONKHORST-PACK","MACDONALD")
          CPPostcondition(kpoint%nkp>=1, cp_failure_level, routineP, error, failure)
    END SELECT
    IF(kpoint%use_real_wfn) THEN
       ! what about inversion symmetry?
       ikloop : DO ik=1,kpoint%nkp
          DO i=1,3
             spez = (kpoint%xkp(i,ik) == 0.0_dp .OR. kpoint%xkp(i,ik) == 0.5_dp)
             IF (.NOT. spez) EXIT ikloop
          END DO
       END DO ikloop
       IF (.NOT.spez .AND. .NOT.(inv .AND. symmo)) THEN
          ! Warning: real wfn might be wrong for this system
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "A calculation using real wavefunctions is requested. "//&
               "We could not determine if the symmetry of the system allows real wavefunctions. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)

       END IF
    END IF

  END SUBROUTINE kpoint_initialize

! *****************************************************************************
!> \brief Initialize the kpoint environment
!> \param kpoint       Kpoint environment
!> \param error        CP2K error handling
! *****************************************************************************
  SUBROUTINE kpoint_env_initialize(kpoint, error)

    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_env_initialize', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3*default_string_length)   :: message
    INTEGER :: comm_cart, comm_inter_kp, comm_kp, ierr, igr, ik, ikk, ngr, &
      niogrp, nkp, nkp_grp, nkp_loc, npe, unit_nr
    INTEGER, DIMENSION(2)                    :: dims, pos
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(2)                    :: rdim
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_cart_type), POINTER         :: cart
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_full, &
                                                para_env_inter_kp, para_env_kp
    TYPE(kpoint_env_p_type), DIMENSION(:), &
      POINTER                                :: kp_env
    TYPE(kpoint_env_type), POINTER           :: kp

    failure = .FALSE.
    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)
    para_env => kpoint%para_env
    CPAssert(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPAssert(.NOT.ASSOCIATED(kpoint%kp_env),cp_failure_level,routineP,error,failure)

    NULLIFY(kp_env)
    nkp = kpoint%nkp
    npe = para_env%num_pe
    IF(npe == 1) THEN
       ! only one process availabe -> owns all kpoints
       ALLOCATE(kp_env(nkp),STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
       DO ik=1,nkp
          NULLIFY(kp_env(ik)%kpoint_env)
          CALL kpoint_env_create(kp_env(ik)%kpoint_env, error)
          kp => kp_env(ik)%kpoint_env
          kp%nkpoint = ik
          kp%wkp = kpoint%wkp(ik)
          kp%xkp(1:3) = kpoint%xkp(1:3,ik)
          kp%is_local = .TRUE.
       END DO
       ALLOCATE(kpoint%kp_dist(2,1),STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
       kpoint%kp_dist(1,1) = 1
       kpoint%kp_dist(2,1) = nkp
       kpoint%kp_range(1) = 1
       kpoint%kp_range(2) = nkp
       kpoint%kp_env => kp_env
       ! parallel environments
       kpoint%para_env_full => para_env
       CALL cp_para_env_retain(para_env, error)
       kpoint%para_env_kp => para_env
       CALL cp_para_env_retain(para_env, error)
       kpoint%para_env_inter_kp => para_env
       CALL cp_para_env_retain(para_env, error)
       kpoint%iogrp = .TRUE.
       kpoint%nkp_groups = 1
    ELSE
       IF(kpoint%parallel_group_size == -1) THEN
          ! maximum parallelization over kpoints
          ngr = MAX(npe/nkp,1)
       ELSE IF(kpoint%parallel_group_size == 0) THEN
          ! no parallelization over kpoints
          ngr = npe
       ELSE IF(kpoint%parallel_group_size > 0) THEN
          ngr = MIN(kpoint%parallel_group_size,npe)
       ELSE
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       nkp_grp = npe/ngr
       ! processor dimensions
       dims(1)=ngr
       dims(2)=nkp_grp
       CPPostcondition(MOD(nkp,nkp_grp)==0,cp_failure_level,routineP,error,failure)
       nkp_loc=nkp/nkp_grp

       IF ((dims(1)*dims(2)/=npe)) THEN
          message = "Number of processors is not divisible by the kpoint group size."
          CALL print_warning(routineN,moduleN,__LINE__,TRIM(message),para_env)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF

       NULLIFY(cart,para_env_full,para_env_kp,para_env_inter_kp)
       CALL mp_cart_create(comm_old=para_env%group, ndims=2, dims=dims, pos=pos, &
                           comm_cart=comm_cart)

       CALL cp_cart_create(cart,comm_cart,ndims=2,owns_group=.TRUE.,error=error)
       CALL cp_para_env_create(para_env_full,comm_cart,owns_group=.FALSE.,error=error)
       rdim(2)=.FALSE.
       rdim(1)=.TRUE.
       CALL mp_cart_sub(comm=comm_cart,rdim=rdim,sub_comm=comm_kp)
       CALL cp_para_env_create(para_env_kp,comm_kp,owns_group=.TRUE.,error=error)
       rdim(2)=.TRUE.
       rdim(1)=.FALSE.
       CALL mp_cart_sub(comm=comm_cart,rdim=rdim,sub_comm=comm_inter_kp)
       CALL cp_para_env_create(para_env_inter_kp,comm_inter_kp,&
            owns_group=.TRUE.,error=error)

       niogrp = 0
       IF(para_env%ionode) niogrp=1
       CALL mp_sum(niogrp,para_env_kp%group)
       kpoint%iogrp = (niogrp==1)

       ! parallel groups
       kpoint%cart => cart
       kpoint%para_env_full => para_env_full
       kpoint%para_env_kp => para_env_kp
       kpoint%para_env_inter_kp => para_env_inter_kp

       ! distribution of kpoints
       ALLOCATE(kpoint%kp_dist(2,nkp_grp),STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
       DO igr = 1,nkp_grp
          kpoint%kp_dist(1:2,igr) = get_limit(nkp,nkp_grp,igr-1)
       END DO
       ! local kpoints
       kpoint%kp_range(1:2) = kpoint%kp_dist(1:2,para_env_inter_kp%mepos+1)
       ALLOCATE(kp_env(nkp_loc),STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
       DO ik=1,nkp_loc
          NULLIFY(kp_env(ik)%kpoint_env)
          ikk = kpoint%kp_range(1) + ik -1
          CALL kpoint_env_create(kp_env(ik)%kpoint_env, error)
          kp => kp_env(ik)%kpoint_env
          kp%nkpoint = ikk
          kp%wkp = kpoint%wkp(ikk)
          kp%xkp(1:3) = kpoint%xkp(1:3,ikk)
          kp%is_local = (ngr == 1)
       END DO

       kpoint%kp_env => kp_env

       logger => cp_error_get_logger(error)
       unit_nr=cp_logger_get_default_io_unit(logger)
       IF (unit_nr>0) THEN
          WRITE(unit_nr,FMT="(T2,A,T71,I10)") "KPOINTS| Number of kpoint groups ",nkp_grp
          WRITE(unit_nr,FMT="(T2,A,T71,I10)") "KPOINTS| Size of each kpoint group",ngr
          WRITE(unit_nr,FMT="(T2,A,T71,I10)") "KPOINTS| Number of kpoints per group",nkp_loc
       END IF
       kpoint%nkp_groups = nkp_grp

    END IF

  END SUBROUTINE kpoint_env_initialize

! *****************************************************************************
!> \brief Initialize a set of MOs and density matrix for each kpoint (kpoint group)
!> \param kpoint  Kpoint environment
!> \param mos     Reference MOs (global)
!> \param error   CP2K error handling
! *****************************************************************************
  SUBROUTINE kpoint_initialize_mos(kpoint, mos, error)

    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_initialize_mos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ic, ierr, ik, is, nao, nc, &
                                                nelectron, nkp_loc, nmo, nspin
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: flexible_electron_count, &
                                                maxocc, n_el_f
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_ao_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_fm_type), POINTER                :: fmlocal
    TYPE(kpoint_env_type), POINTER           :: kp
    TYPE(qs_matrix_pools_type), POINTER      :: mpools

    failure = .FALSE.
    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)

    IF (ASSOCIATED(mos)) THEN

       IF(kpoint%use_real_wfn) THEN
          nc = 1
       ELSE
          nc = 2
       END IF
       nspin = SIZE(mos,1)
       nkp_loc = kpoint%kp_range(2) - kpoint%kp_range(1) + 1
       IF (nkp_loc > 0) THEN
          CPAssert(SIZE(kpoint%kp_env)==nkp_loc,cp_failure_level,routineP,error,failure)
          ! allocate the mo sets, correct number of kpoints (local), real/complex, spin
          DO ik=1,nkp_loc
             kp => kpoint%kp_env(ik)%kpoint_env
             ALLOCATE(kp%mos(nc,nspin),STAT=ierr)
             CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
             DO is=1,nspin
                CALL get_mo_set(mos(is)%mo_set,nao=nao,nmo=nmo,nelectron=nelectron,&
                     n_el_f=n_el_f,maxocc=maxocc,flexible_electron_count=flexible_electron_count)
                DO ic=1,nc
                   NULLIFY(kp%mos(ic,is)%mo_set)
                   CALL allocate_mo_set(kp%mos(ic,is)%mo_set,nao,nmo,nelectron,n_el_f,maxocc,&
                                flexible_electron_count,error)
                END DO
             END DO
          END DO

          ! generate the blacs environment for the kpoint group      
          ! we generate a blacs env for each kpoint group in parallel
          ! we assume here that the group para_env_inter_kp will connect
          ! equivalent parts of fm matrices, i.e. no reshuffeling of processors
          NULLIFY(blacs_env)
          CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=kpoint%para_env_kp, error=error)
          kpoint%blacs_env => blacs_env
   
          ! matrix pools for the kpoint group, information on MOs is transfered using 
          ! generic mos structure
          NULLIFY(mpools)
          CALL mpools_create(mpools=mpools,error=error)
          CALL mpools_rebuild_fm_pools(mpools=mpools,mos=mos,&
               blacs_env=blacs_env,para_env=kpoint%para_env_kp,error=error)
          kpoint%mpools => mpools

          ! allocate density matrices
          CALL mpools_get(mpools,ao_ao_fm_pools=ao_ao_fm_pools,error=error)
          CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool,fmlocal,error=error)
          CALL cp_fm_get_info(fmlocal,matrix_struct=matrix_struct,error=error)
          DO ik=1,nkp_loc
             kp => kpoint%kp_env(ik)%kpoint_env
             ! density matrix
             IF(ASSOCIATED(kp%pmat)) THEN
                DO is=1,SIZE(kp%pmat,2)
                   DO ic=1,SIZE(kp%pmat,1)
                      CALL cp_fm_release(kp%pmat(ic,is)%matrix,error)
                   END DO
                END DO
                DEALLOCATE(kp%pmat,STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
             END IF
             ALLOCATE(kp%pmat(nc,nspin),STAT=ierr)
             CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
             DO is=1,nspin
                DO ic=1,nc
                   CALL cp_fm_create(kp%pmat(ic,is)%matrix,matrix_struct,error=error)
                END DO
             END DO
             ! energy weighted density matrix
             IF(ASSOCIATED(kp%wmat)) THEN
                DO is=1,SIZE(kp%wmat,2)
                   DO ic=1,SIZE(kp%wmat,1)
                      CALL cp_fm_release(kp%wmat(ic,is)%matrix,error)
                   END DO
                END DO
                DEALLOCATE(kp%wmat,STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
             END IF
             ALLOCATE(kp%wmat(nc,nspin),STAT=ierr)
             CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
             DO is=1,nspin
                DO ic=1,nc
                   CALL cp_fm_create(kp%wmat(ic,is)%matrix,matrix_struct,error=error)
                END DO
             END DO
          END DO
          CALL fm_pool_give_back_fm(ao_ao_fm_pools(1)%pool, fmlocal, error=error)

       END IF

    END IF

  END SUBROUTINE kpoint_initialize_mos

! *****************************************************************************
!> \brief Generates the mapping of cell indices and linear RS index
!>        CELL (0,0,0) is always mapped to index 1
!> \param kpoint    Kpoint environment
!> \param sab_nl    Defining neighbour list
!> \param para_env  Parallel environment
!> \param dft_control ...
!> \param error     CP2K error handling
! *****************************************************************************
  SUBROUTINE kpoint_init_cell_index(kpoint, sab_nl, para_env, dft_control, error)

    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_init_cell_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, i3, ic, icount, it, &
                                                ncount
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: index_to_cell
    INTEGER, DIMENSION(3)                    :: cell, itm
    INTEGER, DIMENSION(:, :), POINTER        :: list
    INTEGER, DIMENSION(:, :, :), POINTER     :: cti
    LOGICAL                                  :: failure, new
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    failure = .FALSE.
    CPAssert(ASSOCIATED(kpoint),cp_failure_level,routineP,error,failure)

    ALLOCATE(list(3,125))
    list = 0
    icount = 1

    CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,cell=cell)

       new = .TRUE.
       DO ic=1,icount
          IF (cell(1) == list(1,ic) .AND. cell(2) == list(2,ic) .AND. &
             cell(3) == list(3,ic)) THEN
             new = .FALSE.
             EXIT
          END IF
       END DO
       IF (new) THEN
          icount = icount + 1
          IF(icount > SIZE(list,2)) THEN
             CALL reallocate(list,1,3,1,2*SIZE(list,2))
          END IF
          list(1:3,icount) = cell(1:3)
       END IF

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    itm(1) = MAXVAL(ABS(list(1,1:icount)))
    itm(2) = MAXVAL(ABS(list(2,1:icount)))
    itm(3) = MAXVAL(ABS(list(3,1:icount)))
    CALL mp_max(itm,para_env%group)
    it = MAXVAL(itm(1:3))
    IF(ASSOCIATED(kpoint%cell_to_index)) THEN
       DEALLOCATE(kpoint%cell_to_index)
    END IF
    ALLOCATE(kpoint%cell_to_index(-itm(1):itm(1),-itm(2):itm(2),-itm(3):itm(3)))
    cti => kpoint%cell_to_index
    cti(:,:,:) = 0
    DO ic=1,icount
       i1 = list(1,ic)
       i2 = list(2,ic)
       i3 = list(3,ic)
       cti(i1,i2,i3) = ic
    END DO
    CALL mp_sum(cti,para_env%group)
    ncount = 0
    DO i1=-itm(1),itm(1)
      DO i2=-itm(2),itm(2)
        DO i3=-itm(3),itm(3)
          IF(cti(i1,i2,i3)==0) THEN
            cti(i1,i2,i3) = 1000000
          ELSE
            ncount = ncount + 1
            cti(i1,i2,i3) = (ABS(i1)+ABS(i2)+ABS(i3))*1000+ABS(i3)*100+ABS(i2)*10+ABS(i1)
            cti(i1,i2,i3) = cti(i1,i2,i3) + (i1+i2+i3)
          END IF
        END DO
      END DO
    END DO

    ALLOCATE(index_to_cell(3,ncount))
    DO ic=1,ncount
       cell = MINLOC(cti)
       i1 = cell(1) - 1 - itm(1)
       i2 = cell(2) - 1 - itm(2)
       i3 = cell(3) - 1 - itm(3)
       cti(i1,i2,i3) = 1000000
       index_to_cell(1,ic) = i1
       index_to_cell(2,ic) = i2
       index_to_cell(3,ic) = i3
    END DO
    cti(:,:,:) = 0
    DO ic=1,ncount
       i1 = index_to_cell(1,ic)
       i2 = index_to_cell(2,ic)
       i3 = index_to_cell(3,ic)
       cti(i1,i2,i3) = ic
    END DO

    ! keep pointer to this neighborlist
    kpoint%sab_nl => sab_nl

    ! set number of images
    dft_control%nimages = SIZE(index_to_cell, 2)

    DEALLOCATE(list,index_to_cell)
  END SUBROUTINE kpoint_init_cell_index

! *****************************************************************************
!> \brief Transformation of real space matrices to a kpoint
!> \param rmatrix  Real part of kpoint matrix
!> \param cmatrix  Complex part of kpoint matrix (optional)
!> \param rsmat    Real space matrices
!> \param ispin    Spin index
!> \param xkp      Kpoint coordinates
!> \param cell_to_index   mapping of cell indices to RS index 
!> \param sab_nl   Defining neighbor list
!> \param error    CP2K error handling
! *****************************************************************************
  SUBROUTINE rskp_transform(rmatrix,cmatrix,rsmat,ispin,&
                            xkp,cell_to_index,sab_nl,error)

    TYPE(cp_dbcsr_type), POINTER             :: rmatrix
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: cmatrix
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rsmat
    INTEGER, INTENT(IN)                      :: ispin
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: xkp
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rskp_transform', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ic, icol, irow, jatom
    INTEGER, DIMENSION(3)                    :: cell
    LOGICAL                                  :: do_symmetric, failure, found, &
                                                wfn_real_only
    REAL(KIND=dp)                            :: arg, coskl, fsym, sinkl
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cblock, rblock, rsblock
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    wfn_real_only = .TRUE.
    IF(PRESENT(cmatrix)) wfn_real_only = .FALSE.
    CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl,symmetric=do_symmetric)

    IF(wfn_real_only) THEN
       CALL cp_dbcsr_set(rmatrix,0.0_dp,error=error)
    ELSE
       CALL cp_dbcsr_set(rmatrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(cmatrix,0.0_dp,error=error)
    END IF

    CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,iatom=iatom,jatom=jatom,cell=cell)

       fsym = 1.0_dp
       irow = iatom
       icol = jatom
       IF( do_symmetric .AND. (iatom > jatom) ) THEN
          irow = jatom
          icol = iatom
          fsym = -1.0_dp
       END IF

       arg = REAL(cell(1),dp)*xkp(1)+REAL(cell(2),dp)*xkp(2)+REAL(cell(3),dp)*xkp(3) 
       coskl = COS(twopi*arg)
       sinkl = fsym*SIN(twopi*arg)

       ic = cell_to_index(cell(1),cell(2),cell(3))
       CALL cp_dbcsr_get_block_p(matrix=rsmat(ispin,ic)%matrix,row=irow,col=icol,&
                                 block=rsblock,found=found)
       CPPostcondition(found,cp_failure_level,routineP,error,failure)

       IF(wfn_real_only) THEN
          CALL cp_dbcsr_get_block_p(matrix=rmatrix,row=irow,col=icol,&
                                    block=rblock,found=found)
          CPPostcondition(found,cp_failure_level,routineP,error,failure)
          rblock = rblock + coskl*rsblock
       ELSE
          CALL cp_dbcsr_get_block_p(matrix=rmatrix,row=irow,col=icol,&
                                    block=rblock,found=found)
          CPPostcondition(found,cp_failure_level,routineP,error,failure)
          CALL cp_dbcsr_get_block_p(matrix=cmatrix,row=irow,col=icol,&
                                    block=cblock,found=found)
          CPPostcondition(found,cp_failure_level,routineP,error,failure)
          rblock = rblock + coskl*rsblock
          cblock = cblock + sinkl*rsblock
       END IF

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

  END SUBROUTINE rskp_transform

! *****************************************************************************
!> \brief Given the eigenvalues of all kpoints, calculates the occupation numbers
!> \param kpoint  Kpoint environment
!> \param smear   Smearing information
!> \param error   CP2K error handling
! *****************************************************************************
  SUBROUTINE kpoint_set_mo_occupation(kpoint,smear,error)

    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(smear_type), POINTER                :: smear
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_set_mo_occupation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, ikpgr, ispin, kplocal, &
                                                nb, nelectron, nkp, nmo, &
                                                nspin, stat
    INTEGER, DIMENSION(2)                    :: kp_range
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: kTS, mu, nel
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: weig, wocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues, occupation, wkp
    TYPE(cp_para_env_type), POINTER          :: para_env_inter_kp
    TYPE(kpoint_env_type), POINTER           :: kp
    TYPE(mo_set_type), POINTER               :: mo_set

    failure = .FALSE.
    ! first collect all the eigenvalues
    CALL get_kpoint_info(kpoint,nkp=nkp,error=error)
    kp => kpoint%kp_env(1)%kpoint_env
    nspin = SIZE(kp%mos,2)
    mo_set => kp%mos(1,1)%mo_set
    CALL get_mo_set(mo_set,nmo=nmo,nelectron=nelectron)
    IF(nspin==2) THEN
       CALL get_mo_set(kp%mos(1,2)%mo_set,nmo=nb)
       CPAssert(nmo==nb,cp_failure_level,routineP,error,failure)
    END IF
    ALLOCATE(weig(nmo,nkp,nspin),wocc(nmo,nkp,nspin),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    weig = 0.0_dp
    wocc = 0.0_dp
    CALL get_kpoint_info(kpoint,kp_range=kp_range,error=error)
    kplocal = kp_range(2) - kp_range(1) + 1
    DO ikpgr = 1,kplocal
       ik = kp_range(1)+ikpgr-1
       kp => kpoint%kp_env(ikpgr)%kpoint_env
       DO ispin=1,nspin
          mo_set => kp%mos(1,ispin)%mo_set
          CALL get_mo_set(mo_set,eigenvalues=eigenvalues)
          weig(1:nmo,ik,ispin) = eigenvalues(1:nmo)
       END DO
    END DO
    CALL get_kpoint_info(kpoint,para_env_inter_kp=para_env_inter_kp,error=error)
    CALL mp_sum(weig,para_env_inter_kp%group)

    nel = REAL(nelectron,KIND=dp)
    CALL get_kpoint_info(kpoint,wkp=wkp,error=error)
    IF(smear%do_smear) THEN
       ! finite electronic temperature
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    ELSE
       ! fixed occupations (2/1)
       IF(nspin==1) THEN
          CALL Fermikp(wocc(:,:,1),mu,kTS,weig(:,:,1),nel,wkp,0.0_dp,2.0_dp)
       ELSE
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    DO ikpgr = 1,kplocal
       ik = kp_range(1)+ikpgr-1
       kp => kpoint%kp_env(ikpgr)%kpoint_env
       DO ispin=1,nspin
          mo_set => kp%mos(1,ispin)%mo_set
          CALL get_mo_set(mo_set,eigenvalues=eigenvalues,occupation_numbers=occupation)
          eigenvalues(1:nmo) = weig(1:nmo,ik,ispin)
          occupation(1:nmo) = wocc(1:nmo,ik,ispin)
       END DO
    END DO

    DEALLOCATE(weig,wocc,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE kpoint_set_mo_occupation

! *****************************************************************************
!> \brief Calculate kpoint density matrices (rho(k), owned by kpoint groups)
!> \param kpoint    kpoint environment
!> \param energy_weighted  calculate energy weighted density matrix
!> \param error     CP2K error handling
! *****************************************************************************
  SUBROUTINE kpoint_density_matrices(kpoint,energy_weighted,error)

    TYPE(kpoint_type), POINTER               :: kpoint
    LOGICAL, OPTIONAL                        :: energy_weighted
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_density_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikpgr, ispin, kplocal, nao, &
                                                nmo, nspin
    INTEGER, DIMENSION(2)                    :: kp_range
    LOGICAL                                  :: failure, wtype
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues, occupation
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_fm_type), POINTER                :: cpmat, fwork, pmat, rpmat
    TYPE(kpoint_env_type), POINTER           :: kp
    TYPE(mo_set_type), POINTER               :: mo_set

    failure = .FALSE.

    ! which type of density matrix
    IF(PRESENT(energy_weighted)) THEN
       wtype = energy_weighted
    ELSE
       ! default is normal density matrix
       wtype = .FALSE.
    END IF

    ! work matrix
    mo_set => kpoint%kp_env(1)%kpoint_env%mos(1,1)%mo_set
    CALL get_mo_set(mo_set,nao=nao,nmo=nmo)
    CALL cp_fm_get_info(mo_set%mo_coeff,matrix_struct=matrix_struct,error=error)
    CALL cp_fm_create(fwork,matrix_struct,error=error)

    CALL get_kpoint_info(kpoint,kp_range=kp_range,error=error)
    kplocal = kp_range(2) - kp_range(1) + 1
    DO ikpgr = 1,kplocal
       kp => kpoint%kp_env(ikpgr)%kpoint_env
       nspin = SIZE(kp%mos,2)
       DO ispin=1,nspin
          mo_set => kp%mos(1,ispin)%mo_set
          IF(wtype) THEN
             CALL get_mo_set(mo_set,eigenvalues=eigenvalues)
          END IF
          IF(kpoint%use_real_wfn) THEN
             IF(wtype) THEN
                pmat => kp%wmat(1,ispin)%matrix
             ELSE
                pmat => kp%pmat(1,ispin)%matrix
             END IF
             CALL get_mo_set(mo_set,occupation_numbers=occupation)
             CALL cp_fm_to_fm(mo_set%mo_coeff,fwork,error)
             CALL cp_fm_column_scale(fwork,occupation)
             IF(wtype) THEN
                CALL cp_fm_column_scale(fwork,eigenvalues)
             END IF
             CALL cp_gemm("N","T",nao,nao,nmo,1.0_dp,mo_set%mo_coeff,fwork,0.0_dp,pmat,error)
          ELSE
             IF(wtype) THEN
                rpmat => kp%wmat(1,ispin)%matrix
                cpmat => kp%wmat(2,ispin)%matrix
             ELSE
                rpmat => kp%pmat(1,ispin)%matrix
                cpmat => kp%pmat(2,ispin)%matrix
             END IF
             CALL get_mo_set(mo_set,occupation_numbers=occupation)
             CALL cp_fm_to_fm(mo_set%mo_coeff,fwork,error)
             CALL cp_fm_column_scale(fwork,occupation)
             IF(wtype) THEN
                CALL cp_fm_column_scale(fwork,eigenvalues)
             END IF
             ! Re(c)*Re(c)
             CALL cp_gemm("N","T",nao,nao,nmo,1.0_dp,mo_set%mo_coeff,fwork,0.0_dp,rpmat,error)
             mo_set => kp%mos(2,ispin)%mo_set
             ! Im(c)*Re(c)
             CALL cp_gemm("N","T",nao,nao,nmo,-1.0_dp,mo_set%mo_coeff,fwork,0.0_dp,cpmat,error)
             ! Re(c)*Im(c)
             CALL cp_gemm("N","T",nao,nao,nmo,1.0_dp,fwork,mo_set%mo_coeff,1.0_dp,cpmat,error)
             CALL cp_fm_to_fm(mo_set%mo_coeff,fwork,error)
             CALL cp_fm_column_scale(fwork,occupation)
             IF(wtype) THEN
                CALL cp_fm_column_scale(fwork,eigenvalues)
             END IF
             ! Im(c)*Im(c)
             CALL cp_gemm("N","T",nao,nao,nmo,1.0_dp,mo_set%mo_coeff,fwork,1.0_dp,rpmat,error)
          END IF
       END DO
    END DO

    CALL cp_fm_release(fwork,error=error)

  END SUBROUTINE kpoint_density_matrices

! *****************************************************************************
!> \brief generate real space desnity matrices in DBCSR format
!> \param kpoint  Kpoint environment
!> \param denmat  Real space (DBCSR) density matrices
!> \param wtype   True = energy weighted density matrix
!>                False = normal density matrix
!> \param tempmat DBCSR matrix to be used as template
!> \param sab_nl ...
!> \param fmwork  FM work matrices (kpoint group)
!> \param error   CP2K error handling
! *****************************************************************************
  SUBROUTINE kpoint_density_transform(kpoint,denmat,wtype,tempmat,sab_nl,fmwork,error)

    TYPE(kpoint_type), POINTER               :: kpoint
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: denmat
    LOGICAL, INTENT(IN)                      :: wtype
    TYPE(cp_dbcsr_type), POINTER             :: tempmat
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fmwork
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_density_transform', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bcntxt, iatom, ic, icell, &
                                                icol, ik, ikk, irow, ispin, &
                                                jatom, nc, nimg, nkp, nspin, &
                                                stat
    INTEGER, DIMENSION(3)                    :: cell
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: do_symmetric, failure, found, &
                                                my_kpgrp, real_only
    REAL(KIND=dp)                            :: arg, coskl, fc, sinkl
    REAL(KIND=dp), DIMENSION(:), POINTER     :: wkp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cblock, dblock, rblock, xkp
    TYPE(cp_dbcsr_type), POINTER             :: cpmat, rpmat
    TYPE(cp_fm_type), POINTER                :: fmdummy
    TYPE(cp_para_env_type), POINTER          :: para_env_inter_kp
    TYPE(kpoint_env_type), POINTER           :: kp
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator

    failure = .FALSE.

    ! work storage
    ALLOCATE(rpmat,cpmat,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(rpmat,error=error)
    CALL cp_dbcsr_init(cpmat,error=error)
    CALL cp_dbcsr_create(rpmat,template=tempmat,matrix_type=dbcsr_type_symmetric,error=error)
    CALL cp_dbcsr_create(cpmat,template=tempmat,matrix_type=dbcsr_type_antisymmetric,error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(rpmat,sab_nl,error)
    CALL cp_dbcsr_alloc_block_from_nbl(cpmat,sab_nl,error)

    CALL get_kpoint_info(kpoint,nkp=nkp,xkp=xkp,wkp=wkp,&
                         para_env_inter_kp=para_env_inter_kp,&
                         sab_nl=sab_nl,cell_to_index=cell_to_index,&
                         error=error)
    ! initialize real space density matrices
    kp => kpoint%kp_env(1)%kpoint_env
    nspin = SIZE(kp%mos,2)
    nc = SIZE(kp%mos,1)
    nimg = SIZE(denmat, 2)
    real_only = (nc==1)
    bcntxt = kpoint%blacs_env_all%group
    NULLIFY(fmdummy)
    DO ispin=1,nspin
       DO ic=1,nimg
          CALL cp_dbcsr_set(denmat(ispin,ic)%matrix,0.0_dp,error=error)
       END DO
       !
       DO ik=1,nkp
          my_kpgrp = (ik >= kpoint%kp_range(1) .AND. ik <= kpoint%kp_range(2))
          IF(my_kpgrp) THEN
            ikk = ik - kpoint%kp_range(1) + 1
            kp => kpoint%kp_env(ikk)%kpoint_env
          ELSE
            NULLIFY(kp)
          END IF
          ! collect this density matrix on all processors
          CPPrecondition(SIZE(fmwork)>=nc,cp_failure_level,routineP,error,failure)
          DO ic=1,nc
             IF(my_kpgrp) THEN
                ! choose density matrix
                IF(wtype) THEN
                   CALL cp_fm_copy_general(kp%wmat(ic,ispin)%matrix,fmwork(ic)%matrix,bcntxt,error)
                ELSE
                   CALL cp_fm_copy_general(kp%pmat(ic,ispin)%matrix,fmwork(ic)%matrix,bcntxt,error)
                END IF
             ELSE
                CALL cp_fm_copy_general(fmdummy,fmwork(ic)%matrix,bcntxt,error)
             END IF
          END DO
          ! reduce to dbcsr storage
          IF(real_only) THEN
             CALL copy_fm_to_dbcsr(fmwork(1)%matrix,rpmat,keep_sparsity=.TRUE.,error=error)
          ELSE
             CALL copy_fm_to_dbcsr(fmwork(1)%matrix,rpmat,keep_sparsity=.TRUE.,error=error)
             ! it seems this copy to a antisymmetric dbcsr matrix changes sign!
             CALL copy_fm_to_dbcsr(fmwork(2)%matrix,cpmat,keep_sparsity=.TRUE.,error=error)
          END IF

          ! transformation
          CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl,symmetric=do_symmetric)
          CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
          DO WHILE (neighbor_list_iterate(nl_iterator)==0)
             CALL get_iterator_info(nl_iterator,iatom=iatom,jatom=jatom,cell=cell)

             ! fc = +/-1 sign change from above, is used to get correct sign due to
             !           a<->b ==> cell = -cell
             irow = iatom
             icol = jatom
             fc = -1.0_dp
             IF( do_symmetric .AND. iatom > jatom) THEN
                irow = jatom
                icol = iatom
                fc = 1.0_dp
             END IF
             icell = cell_to_index(cell(1),cell(2),cell(3))
             arg = REAL(cell(1),dp)*xkp(1,ik)+REAL(cell(2),dp)*xkp(2,ik)+REAL(cell(3),dp)*xkp(3,ik) 
             coskl = wkp(ik)*COS(twopi*arg)
             sinkl = wkp(ik)*fc*SIN(twopi*arg)

             CALL cp_dbcsr_get_block_p(matrix=denmat(ispin,icell)%matrix,row=irow,col=icol,&
                                       block=dblock,found=found)
             CPPostcondition(found,cp_failure_level,routineP,error,failure)

             IF(real_only) THEN
                CALL cp_dbcsr_get_block_p(matrix=rpmat,row=irow,col=icol,&
                                          block=rblock,found=found)
                CPPostcondition(found,cp_failure_level,routineP,error,failure)
                dblock = dblock + coskl*rblock
             ELSE
                CALL cp_dbcsr_get_block_p(matrix=rpmat,row=irow,col=icol,&
                                          block=rblock,found=found)
                CPPostcondition(found,cp_failure_level,routineP,error,failure)
                CALL cp_dbcsr_get_block_p(matrix=cpmat,row=irow,col=icol,&
                                          block=cblock,found=found)
                CPPostcondition(found,cp_failure_level,routineP,error,failure)
                dblock = dblock + coskl*rblock
                dblock = dblock + sinkl*cblock
             END IF
          END DO
          CALL neighbor_list_iterator_release(nl_iterator)

       END DO
    END DO

    CALL cp_dbcsr_deallocate_matrix(rpmat,error)
    CALL cp_dbcsr_deallocate_matrix(cpmat,error)

  END SUBROUTINE kpoint_density_transform

! *****************************************************************************

END MODULE kpoint_methods
