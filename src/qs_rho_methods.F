!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods of the rho structure (defined in qs_rho_types)
!> \par History
!>      08.2002 created [fawzi]
!>      08.2014 kpoints [JGH]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_rho_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_distribution, cp_dbcsr_init, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_type_symmetric
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE kpoint_types,                    ONLY: kpoint_type
  USE lri_environment_methods,         ONLY: calculate_lri_densities
  USE lri_environment_types,           ONLY: lri_density_type,&
                                             lri_environment_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_collocate_density,            ONLY: calculate_drho_elec,&
                                             calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_local_rho_types,              ONLY: local_rho_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_types,                    ONLY: qs_rho_clear,&
                                             qs_rho_get,&
                                             qs_rho_set,&
                                             qs_rho_type
  USE task_list_types,                 ONLY: task_list_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_methods'

  PUBLIC :: qs_rho_update_rho, qs_rho_rebuild, duplicate_rho_type

CONTAINS

! *****************************************************************************
!> \brief rebuilds rho (if necessary allocating and initializing it)
!> \param rho the rho type to rebuild (defaults to qs_env%rho)
!> \param qs_env the environment to which rho belongs
!> \param rebuild_ao if it is necessary to rebuild rho_ao. Defaults to true.
!> \param rebuild_grids if it in necessary to rebuild rho_r and rho_g.
!>        Defaults to false.
!> \param basis_set_id ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created replacing qs_rho_create and qs_env_rebuild_rho[fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs updated  pw pools, s, s_mstruct and h in qs_env.
!>      The use of p to keep the structure of h (needed for the forces)
!>      is ugly and should be removed.
!>      Change so that it does not allocate a subcomponent if it is not
!>      associated and not requested?
! *****************************************************************************
  SUBROUTINE qs_rho_rebuild(rho, qs_env, rebuild_ao, rebuild_grids, basis_set_id, error)
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: rebuild_ao, rebuild_grids
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_rho_rebuild', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: headline
    INTEGER                                  :: handle, i, ic, &
                                                my_basis_set_id, nimg, nspins
    LOGICAL                                  :: do_kpoints, failure, &
                                                my_rebuild_ao, &
                                                my_rebuild_grids
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_s_kp, rho_ao_kp
    TYPE(cp_dbcsr_type), POINTER             :: refmatrix, tmatrix
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: drho_g, drho_r, rho_g, rho_r, &
                                                tau_g, tau_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    NULLIFY (pw_env,auxbas_pw_pool,matrix_s,matrix_s_kp,dft_control)
    NULLIFY (tot_rho_r,rho_ao_kp,rho_r,rho_g,drho_r,drho_g,tau_r,tau_g)
    NULLIFY (sab_orb)
    CALL timeset(routineN,handle)
    failure=.FALSE.
    my_rebuild_ao=.TRUE.
    my_rebuild_grids=.TRUE.
    IF (PRESENT(rebuild_ao)) my_rebuild_ao=rebuild_ao
    IF (PRESENT(rebuild_grids)) my_rebuild_grids=rebuild_grids

    IF (PRESENT(basis_set_id)) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    CALL get_qs_env(qs_env,&
                    kpoints=kpoints,&
                    do_kpoints=do_kpoints,&
                    pw_env=pw_env,&
                    dft_control=dft_control,&
                    error=error)

    nimg = dft_control%nimages

    SELECT CASE (my_basis_set_id)
      CASE (use_orb_basis_set)
        CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, sab_orb=sab_orb, error=error)
        refmatrix => matrix_s_kp(1,1)%matrix
      CASE (use_aux_fit_basis_set)
        CPPostcondition(.NOT.do_kpoints, cp_failure_level,routineP,error,failure)
        CALL get_qs_env(qs_env, matrix_s_aux_fit=matrix_s, sab_aux_fit=sab_orb, error=error)
        refmatrix => matrix_s(1)%matrix
    END SELECT


    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
    nspins=dft_control%nspins

    IF (.NOT.ASSOCIATED(rho)) STOP routineP//" rho not associated"
    CALL qs_rho_get(rho,&
                    tot_rho_r=tot_rho_r,&
                    rho_ao_kp=rho_ao_kp,&
                    rho_r=rho_r,&
                    rho_g=rho_g,&
                    drho_r=drho_r,&
                    drho_g=drho_g,&
                    tau_r=tau_r,&
                    tau_g=tau_g,&
                    error=error)

    IF (.NOT.ASSOCIATED(tot_rho_r)) THEN
       ALLOCATE(tot_rho_r(nspins))
       tot_rho_r=0.0_dp
       CALL qs_rho_set(rho, tot_rho_r=tot_rho_r, error=error)
    ENDIF

    ! rho_ao
    IF (my_rebuild_ao.OR.(.NOT.ASSOCIATED(rho_ao_kp))) THEN
       IF (ASSOCIATED(rho_ao_kp)) &
          CALL cp_dbcsr_deallocate_matrix_set(rho_ao_kp,error=error)
       ! Create a new density matrix set
       CALL cp_dbcsr_allocate_matrix_set(rho_ao_kp,nspins,nimg,error=error)
       CALL qs_rho_set(rho, rho_ao_kp=rho_ao_kp, error=error)
       DO i=1,nspins
          DO ic=1,nimg
             IF (nspins > 1) THEN
                IF (i == 1) THEN
                   headline = "DENSITY MATRIX FOR ALPHA SPIN"
                ELSE
                   headline = "DENSITY MATRIX FOR BETA SPIN"
                END IF
             ELSE
                headline = "DENSITY MATRIX"
             END IF
             ALLOCATE(rho_ao_kp(i,ic)%matrix)
             tmatrix => rho_ao_kp(i,ic)%matrix
             CALL cp_dbcsr_init(tmatrix, error=error)
             CALL cp_dbcsr_create(matrix=tmatrix, name=TRIM(headline),&
                  dist=cp_dbcsr_distribution(refmatrix), matrix_type=dbcsr_type_symmetric,&
                  row_blk_size=cp_dbcsr_row_block_sizes(refmatrix),&
                  col_blk_size=cp_dbcsr_col_block_sizes(refmatrix),&
                  nblks=0, nze=0, error=error)
             CALL cp_dbcsr_alloc_block_from_nbl(tmatrix,sab_orb,error=error)
             CALL cp_dbcsr_set(tmatrix,0.0_dp,error=error)
          END DO
       END DO
    END IF

    ! rho_r
    IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho_r)) THEN
       IF (.NOT.failure) THEN
          IF (ASSOCIATED(rho_r)) THEN
             DO i=1,SIZE(rho_r)
                CALL pw_release(rho_r(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_r)
          END IF
          ALLOCATE(rho_r(nspins))
          CALL qs_rho_set(rho, rho_r=rho_r, error=error)
       END IF
       IF (.NOT.failure) THEN
          DO i=1,nspins
             CALL pw_pool_create_pw(auxbas_pw_pool,rho_r(i)%pw,&
                  use_data=REALDATA3D,in_space=REALSPACE,error=error)
          END DO
       END IF
    END IF

    ! rho_g
    IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho_g)) THEN
       IF (.NOT.failure) THEN
          IF (ASSOCIATED(rho_g)) THEN
             DO i=1,SIZE(rho_g)
                CALL pw_release(rho_g(i)%pw,error=error)
             END DO
             DEALLOCATE(rho_g)
          END IF
          ALLOCATE(rho_g(nspins))
          CALL qs_rho_set(rho, rho_g=rho_g, error=error)
       END IF
       IF (.NOT.failure) THEN
          DO i=1,nspins
             CALL pw_pool_create_pw(auxbas_pw_pool,rho_g(i)%pw,&
                  use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
          END DO
       END IF
    END IF

    ! allocate drho_r and drho_g if xc_deriv_collocate
    IF (dft_control%drho_by_collocation) THEN
       ! drho_r
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(drho_r)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(drho_r)) THEN
                DO i=1,SIZE(drho_r)
                   CALL pw_release(drho_r(i)%pw,error=error)
                END DO
                DEALLOCATE(drho_r)
             END IF
             ALLOCATE(drho_r(3*nspins))
             CALL qs_rho_set(rho, drho_r=drho_r, error=error)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,3*nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,drho_r(i)%pw,&
                     use_data=REALDATA3D,in_space=REALSPACE,error=error)
             END DO
          END IF
       END IF
       ! drho_g
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(drho_g)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(drho_g)) THEN
                DO i=1,SIZE(drho_g)
                   CALL pw_release(drho_g(i)%pw,error=error)
                END DO
                DEALLOCATE(drho_g)
             END IF
             ALLOCATE(drho_g(3*nspins))
             CALL qs_rho_set(rho, drho_g=drho_g, error=error)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,3*nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,drho_g(i)%pw,&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
             END DO
          END IF
       END IF
    END IF

    ! allocate tau_r and tau_g if use_kinetic_energy_density
    IF (dft_control%use_kinetic_energy_density) THEN
       ! tau_r
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(tau_r)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(tau_r)) THEN
                DO i=1,SIZE(tau_r)
                   CALL pw_release(tau_r(i)%pw,error=error)
                END DO
                DEALLOCATE(tau_r)
             END IF
             ALLOCATE(tau_r(nspins))
             CALL qs_rho_set(rho, tau_r=tau_r, error=error)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,tau_r(i)%pw,&
                     use_data=REALDATA3D,in_space=REALSPACE,error=error)
             END DO
          END IF
       END IF

       ! tau_g
       IF (my_rebuild_grids.OR..NOT.ASSOCIATED(tau_g)) THEN
          IF (.NOT.failure) THEN
             IF (ASSOCIATED(tau_g)) THEN
                DO i=1,SIZE(tau_g)
                   CALL pw_release(tau_g(i)%pw,error=error)
                END DO
                DEALLOCATE(tau_g)
             END IF
             ALLOCATE(tau_g(nspins))
             CALL qs_rho_set(rho, tau_g=tau_g, error=error)
          END IF
          IF (.NOT.failure) THEN
             DO i=1,nspins
                CALL pw_pool_create_pw(auxbas_pw_pool,tau_g(i)%pw,&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
             END DO
          END IF
       END IF
    END IF  ! use_kinetic_energy_density

    CALL timestop(handle)
END SUBROUTINE qs_rho_rebuild

! *****************************************************************************
!> \brief updates rho_r and rho_g to the rho%rho_ao.
!>      if use_kinetic_energy_density also computes tau_r and tau_g
!> \param rho_struct the rho structure that should be updated
!> \param qs_env the qs_env rho_struct refers to
!>        the integrated charge in r space
!> \param local_rho_set ...
!> \param task_list_external ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_rho_update_rho(rho_struct, qs_env, local_rho_set, task_list_external, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(local_rho_type), OPTIONAL, POINTER  :: local_rho_set
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list_external
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_update_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins
    LOGICAL                                  :: failure, gapw, gapw_xc
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r, tot_rho_r_xc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao, rho_xc_ao
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(lri_density_type), POINTER          :: lri_density
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: drho_g, drho_r, drho_xc_g, &
                                                rho_g, rho_r, rho_xc_g, &
                                                rho_xc_r, tau_g, tau_r, &
                                                tau_xc_g, tau_xc_r
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho_xc
    TYPE(task_list_type), POINTER            :: task_list

    CALL timeset(routineN,handle)

    failure=.FALSE.

    NULLIFY(dft_control,rho_xc, ks_env, rho_ao, rho_r, rho_g, drho_r, drho_g, tau_r, tau_g)
    NULLIFY(rho_xc_ao,rho_xc_g,rho_xc_r,drho_xc_g, tau_xc_r, tau_xc_g, tot_rho_r, tot_rho_r_xc)
    NULLIFY(lri_env, lri_density, para_env, atomic_kind_set)

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

      CALL get_qs_env(qs_env,&
                      ks_env=ks_env,&
                      dft_control=dft_control,&
                      task_list=task_list,&
                      lri_density=lri_density,&
                      lri_env=lri_env,&
                      atomic_kind_set=atomic_kind_set,&
                      para_env=para_env,&
                      error=error)

      CALL qs_rho_get(rho_struct,&
                      rho_r=rho_r,&
                      rho_g=rho_g,&
                      tot_rho_r=tot_rho_r,&
                      drho_r=drho_r,&
                      drho_g=drho_g,&
                      tau_r=tau_r,&
                      tau_g=tau_g,&
                      error=error)

      IF (PRESENT(task_list_external)) task_list => task_list_external 

      nspins=dft_control%nspins
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      IF(dft_control%qs_control%semi_empirical) THEN
        !
        CALL qs_rho_set(rho_struct, rho_r_valid=.FALSE., rho_g_valid=.FALSE., error=error)
      ELSEIF(dft_control%qs_control%dftb) THEN
        !
        CALL qs_rho_set(rho_struct, rho_r_valid=.FALSE., rho_g_valid=.FALSE., error=error)
      ELSEIF(dft_control%qs_control%scptb) THEN
        !
        CALL qs_rho_set(rho_struct, rho_r_valid=.FALSE., rho_g_valid=.FALSE., error=error)
      ELSEIF(dft_control%qs_control%lrigpw) THEN
        CALL qs_rho_get(rho_struct, rho_ao=rho_ao, error=error)
        CALL calculate_lri_densities(lri_env,lri_density,qs_env,rho_ao,&
                                     lri_rho_struct=rho_struct,&
                                     atomic_kind_set=atomic_kind_set,&
                                     para_env=para_env,error=error)
        CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE., error=error)
      ELSE
        CALL qs_rho_get(rho_struct, rho_ao=rho_ao, error=error)
        DO ispin=1,nspins
           CALL calculate_rho_elec(matrix_p=rho_ao(ispin)%matrix,&
                rho=rho_r(ispin),&
                rho_gspace=rho_g(ispin),&
                total_rho=tot_rho_r(ispin),&
                ks_env=ks_env,soft_valid=gapw,&
                task_list_external=task_list, error=error)
        END DO
        CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE., error=error)
      END IF

      ! if needed compute also the gradient of the density
      IF (dft_control%drho_by_collocation) THEN
         CALL qs_rho_get(rho_struct, rho_ao=rho_ao, error=error)
         DO ispin=1,nspins
           CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
           CALL calculate_drho_elec(matrix_p=rho_ao(ispin)%matrix,&
             drho=drho_r(3*(ispin-1)+1:3*ispin),&
             drho_gspace=drho_g(3*(ispin-1)+1:3*ispin),&
             qs_env=qs_env,soft_valid=gapw,error=error)
         END DO
         CALL qs_rho_set(rho_struct, drho_r_valid=.TRUE., drho_g_valid=.TRUE., error=error)
      ENDIF

      ! if needed compute also the kinetic energy density
      IF (dft_control%use_kinetic_energy_density) THEN
         IF(dft_control%qs_control%lrigpw) THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
             "LRIGPW not implemented for meta-GGAs")
         ENDIF
         CALL qs_rho_get(rho_struct, rho_ao=rho_ao, error=error)
         DO ispin=1,nspins
            CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
            CALL calculate_rho_elec(matrix_p=rho_ao(ispin)%matrix,&
                 rho=tau_r(ispin),&
                 rho_gspace=tau_g(ispin),&
                 total_rho=dum, &  ! presumably not meaningful
                 ks_env=ks_env, soft_valid=gapw, &
                 compute_tau=.TRUE.,&
                 error=error)
         END DO
         CALL qs_rho_set(rho_struct, tau_r_valid=.TRUE., tau_g_valid=.TRUE., error=error)
      ENDIF

     ! GAPW o GAPW_XC require the calculation of hard and soft local densities
      IF(gapw ) THEN
         CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
         CALL qs_rho_get(rho_struct, rho_ao=rho_ao, error=error)
         IF(PRESENT(local_rho_set)) THEN
            CALL calculate_rho_atom_coeff(qs_env,rho_ao,local_rho_set%rho_atom_set,error=error)
         ELSE
            CALL calculate_rho_atom_coeff(qs_env,rho_ao,error=error)
         ENDIF
      ENDIF
      IF(gapw_xc) THEN
         CPPrecondition(.NOT.PRESENT(task_list_external),cp_failure_level,routineP,error,failure)
         CALL get_qs_env(qs_env=qs_env,rho_xc=rho_xc,error=error)
         CALL qs_rho_get(rho_xc,&
                         rho_ao=rho_xc_ao,&
                         rho_r=rho_xc_r,&
                         rho_g=rho_xc_g,&
                         tot_rho_r=tot_rho_r_xc,&
                         drho_g=drho_xc_g,&
                         tau_r=tau_xc_r,&
                         tau_g=tau_xc_g,&
                         error=error)
         CALL calculate_rho_atom_coeff(qs_env,rho_ao,error=error)
         ! copy rho_ao into rho_xc_ao
         DO ispin=1,nspins
           CALL cp_dbcsr_copy(rho_xc_ao(ispin)%matrix, rho_ao(ispin)%matrix, error=error)
         END DO
         DO ispin=1,nspins
           CALL calculate_rho_elec(matrix_p=rho_xc_ao(ispin)%matrix,&
                rho=rho_xc_r(ispin),&
                rho_gspace=rho_xc_g(ispin),&
                total_rho=tot_rho_r_xc(ispin),&
                ks_env=ks_env,soft_valid=gapw_xc,error=error)
         END DO
         CALL qs_rho_set(rho_xc, rho_r_valid=.TRUE., rho_g_valid=.TRUE., error=error)
         ! if needed compute also the gradient of the density
         IF (dft_control%drho_by_collocation) THEN
            DO ispin=1,nspins
              CALL calculate_drho_elec(matrix_p=rho_xc_ao(ispin)%matrix,&
                drho=rho_xc_r(3*(ispin-1)+1:3*ispin),&
                drho_gspace=drho_xc_g(3*(ispin-1)+1:3*ispin),&
                qs_env=qs_env,soft_valid=gapw_xc,error=error)
            END DO
            CALL qs_rho_set(rho_xc, drho_r_valid=.TRUE., drho_g_valid=.TRUE., error=error)
         ENDIF
         ! if needed compute also the kinetic energy density
         IF (dft_control%use_kinetic_energy_density) THEN
           DO ispin=1,nspins
             CALL calculate_rho_elec(matrix_p=rho_xc_ao(ispin)%matrix,&
                 rho=tau_xc_r(ispin),&
                 rho_gspace=tau_xc_g(ispin),&
                 total_rho=dum, &  ! presumably not meaningful
                 ks_env=ks_env, soft_valid=gapw_xc, &
                 compute_tau=.TRUE., error=error)
           END DO
           CALL qs_rho_set(rho_xc, tau_r_valid=.TRUE., tau_g_valid=.TRUE., error=error)
         ENDIF
      ENDIF
   END IF

    CALL timestop(handle)

END SUBROUTINE qs_rho_update_rho

! *****************************************************************************
!> \brief Duplicates a pointer physically
!> \param rho_input The rho structure to be duplicated
!> \param rho_output The duplicate rho structure
!> \param qs_env The QS environment from which the auxilliary PW basis-set
!>                pool is taken
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      07.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!> \note
!>      Associated pointers are deallocated, nullified pointers are NOT accepted!
! *****************************************************************************
  SUBROUTINE duplicate_rho_type(rho_input, rho_output, qs_env, error)

    TYPE(qs_rho_type), POINTER               :: rho_input, rho_output
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'duplicate_rho_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins, &
                                                rebuild_each_in
    LOGICAL :: drho_g_valid_in, drho_r_valid_in, failure, rho_g_valid_in, &
      rho_r_valid_in, soft_valid_in, tau_g_valid_in, tau_r_valid_in
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_g_in, tot_rho_g_out, &
                                                tot_rho_r_in, tot_rho_r_out
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao_in, rho_ao_out
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER :: drho_g_in, drho_g_out, &
      drho_r_in, drho_r_out, rho_g_in, rho_g_out, rho_r_in, rho_r_out, &
      tau_g_in, tau_g_out, tau_r_in, tau_r_out
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(dft_control, pw_env, auxbas_pw_pool)
    NULLIFY(rho_ao_in, rho_ao_out)
    NULLIFY(rho_r_in, rho_r_out, rho_g_in, rho_g_out, drho_r_in, drho_r_out)
    NULLIFY(drho_g_in, drho_g_out, tau_r_in,  tau_r_out, tau_g_in, tau_g_out)
    NULLIFY(tot_rho_r_in, tot_rho_r_out, tot_rho_g_in, tot_rho_g_out)

    CPPrecondition(ASSOCIATED(rho_input), cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(rho_output), cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, &
                    error=error)
    CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
    nspins = dft_control%nspins

    CALL qs_rho_clear(rho_output, error=error)

    CALL qs_rho_get(rho_input,&
                    rho_ao       = rho_ao_in,&
                    rho_r        = rho_r_in,&
                    rho_g        = rho_g_in,&
                    drho_r       = drho_r_in,&
                    drho_g       = drho_g_in,&
                    tau_r        = tau_r_in,&
                    tau_g        = tau_g_in,&
                    tot_rho_r    = tot_rho_r_in,&
                    tot_rho_g    = tot_rho_g_in,&
                    rho_g_valid  = rho_g_valid_in,&
                    rho_r_valid  = rho_r_valid_in,&
                    drho_g_valid = drho_g_valid_in,&
                    drho_r_valid = drho_r_valid_in,&
                    tau_r_valid  = tau_r_valid_in,&
                    tau_g_valid  = tau_g_valid_in,&
                    soft_valid   = soft_valid_in,&
                    rebuild_each = rebuild_each_in,&
                    error=error)

    ! rho_ao
    IF (ASSOCIATED(rho_ao_in)) THEN
      CALL cp_dbcsr_allocate_matrix_set(rho_ao_out, nspins,error=error)
      CALL qs_rho_set(rho_output, rho_ao=rho_ao_out, error=error)
      DO i = 1, nspins
         ALLOCATE(rho_ao_out(i)%matrix)
         CALL cp_dbcsr_init(rho_ao_out(i)%matrix, error=error)
         CALL cp_dbcsr_copy(rho_ao_out(i)%matrix,rho_ao_in(i)%matrix, &
              name="myDensityMatrix_for_Spin_"//TRIM(ADJUSTL(cp_to_string(i))), &
              error=error)
         CALL cp_dbcsr_set(rho_ao_out(i)%matrix,0.0_dp,error=error)
      END DO
    END IF

    ! rho_r
    IF (ASSOCIATED(rho_r_in)) THEN
      ALLOCATE(rho_r_out(nspins))
      CALL qs_rho_set(rho_output, rho_r=rho_r_out, error=error)
      DO i = 1, nspins
        CALL pw_pool_create_pw(auxbas_pw_pool, rho_r_out(i)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE, &
                                error=error)
        rho_r_out(i)%pw%cr3d(:,:,:) = rho_r_in(i)%pw%cr3d(:,:,:)
      END DO
    END IF

    ! rho_g
    IF (ASSOCIATED(rho_g_in)) THEN
      ALLOCATE(rho_g_out(nspins))
      CALL qs_rho_set(rho_output, rho_g=rho_g_out, error=error)
      DO i = 1, nspins
        CALL pw_pool_create_pw(auxbas_pw_pool, rho_g_out(i)%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE, error=error)
        rho_g_out(i)%pw%cc(:)      = rho_g_in(i)%pw%cc(:)
      END DO
    END IF

    ! drho_r and drho_g are only needed if calculated by collocation
    IF (dft_control%drho_by_collocation) THEN
      ! drho_r
      IF (ASSOCIATED(drho_r_in)) THEN
        ALLOCATE(drho_r_out(3*nspins))
        CALL qs_rho_set(rho_output, drho_r=drho_r_out, error=error)
        DO i = 1, 3*nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, drho_r_out(i)%pw, &
                                  use_data=REALDATA3D, in_space=REALSPACE, &
                                  error=error)
          drho_r_out(i)%pw%cr3d(:,:,:) = drho_r_in(i)%pw%cr3d(:,:,:)
        END DO
      END IF

      ! drho_g
      IF (ASSOCIATED(drho_g_in)) THEN
        ALLOCATE(drho_g_out(3*nspins))
        CALL qs_rho_set(rho_output, drho_g=drho_g_out, error=error)
        DO i = 1, 3*nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, drho_g_out(i)%pw, &
                                  use_data=COMPLEXDATA1D, &
                                  in_space=RECIPROCALSPACE, error=error)
          drho_g_out(i)%pw%cc(:) = drho_g_in(i)%pw%cc(:)
        END DO
      END IF
    END IF

    ! tau_r and tau_g are only needed in the case of Meta-GGA XC-functionals
    ! are used. Therefore they are only allocated if
    ! dft_control%use_kinetic_energy_density is true
    IF (dft_control%use_kinetic_energy_density) THEN
      ! tau_r
      IF (ASSOCIATED(tau_r_in)) THEN
        ALLOCATE(tau_r_out(nspins))
        CALL qs_rho_set(rho_output, tau_r=tau_r_out, error=error)
        DO i = 1, nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, tau_r_out(i)%pw, &
                                  use_data=REALDATA3D, in_space=REALSPACE, &
                                  error=error)
          tau_r_out(i)%pw%cr3d(:,:,:) = tau_r_in(i)%pw%cr3d(:,:,:)
        END DO
      END IF

      ! tau_g
      IF (ASSOCIATED(tau_g_in)) THEN
        ALLOCATE(tau_g_out(nspins))
        CALL qs_rho_set(rho_output, tau_g=tau_g_out, error=error)
        DO i = 1, nspins
          CALL pw_pool_create_pw(auxbas_pw_pool, tau_g_out(i)%pw, &
                                  use_data=COMPLEXDATA1D, &
                                  in_space=RECIPROCALSPACE, error=error)
          tau_g_out(i)%pw%cc(:) = tau_g_in(i)%pw%cc(:)
        END DO
      END IF
    END IF

    CALL qs_rho_set(rho_output,&
                    rho_g_valid  = rho_g_valid_in,&
                    rho_r_valid  = rho_r_valid_in,&
                    drho_g_valid = drho_g_valid_in,&
                    drho_r_valid = drho_r_valid_in,&
                    tau_r_valid  = tau_r_valid_in,&
                    tau_g_valid  = tau_g_valid_in,&
                    soft_valid   = soft_valid_in,&
                    rebuild_each = rebuild_each_in,&
                    error=error)

    ! tot_rho_r
    IF (ASSOCIATED(tot_rho_r_in)) THEN
      ALLOCATE(tot_rho_r_out(nspins))
      CALL qs_rho_set(rho_output, tot_rho_r=tot_rho_r_out, error=error)
      DO i = 1, nspins
        tot_rho_r_out(i) = tot_rho_r_in(i)
      END DO
    END IF

    ! tot_rho_g
    IF (ASSOCIATED(tot_rho_g_in)) THEN
      ALLOCATE(tot_rho_g_out(nspins))
      CALL qs_rho_set(rho_output, tot_rho_g=tot_rho_g_out, error=error)
      DO i = 1, nspins
         tot_rho_g_out(i) = tot_rho_g_in(i)
      END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE duplicate_rho_type

END MODULE qs_rho_methods

