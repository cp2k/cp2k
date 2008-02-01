!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief methods of the rho structure (defined in qs_rho_types)
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_rho_methods
  USE cp_control_types,                ONLY: dft_control_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_release
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             allocate_matrix_set,&
                                             copy_matrix,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix,&
                                             replicate_matrix_structure
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_methods'
  INTEGER, PRIVATE, SAVE :: last_rho_id_nr=0

  PUBLIC :: qs_rho_update_rho, qs_rho_rebuild, duplicate_rho_type, &
            diff_rho_type, copy_rho_type

!***
CONTAINS

! *****************************************************************************
!> \brief rebuilds rho (if necessary allocating and initializing it)
!> \param rho the rho type to rebuild (defaults to qs_env%rho)
!> \param qs_env the environment to which rho belongs
!> \param rebuild_ao if it is necessary to rebuild rho_ao. Defaults to true.
!> \param rebuild_grids if it in necessary to rebuild rho_r and rho_g.
!>        Defaults to false.
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      needs updated  pw pools, s, s_mstruct and h in qs_env.
!>      The use of p to keep the structure of h (needed for the forces)
!>      is ugly and should be removed.
!>      Change so that it does not allocate a subcomponent if it is not
!>      associated and not requested?
!> \par History
!>      11.2002 created replacing qs_rho_create and qs_env_rebuild_rho[fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_rho_rebuild(rho, qs_env, rebuild_ao, rebuild_grids, gapw_xc, error)
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: rebuild_ao, rebuild_grids, &
                                                gapw_xc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_rebuild', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins, stat
    LOGICAL                                  :: failure, my_gapw_xc, &
                                                my_rebuild_ao, &
                                                my_rebuild_grids
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

  NULLIFY(pw_env,auxbas_pw_pool,matrix_s,dft_control)
  CALL timeset(routineN,'I',' ',handle)
  failure=.FALSE.
  my_gapw_xc=.FALSE.
  my_rebuild_ao=.TRUE.
  my_rebuild_grids=.TRUE.
  IF (PRESENT(rebuild_ao)) my_rebuild_ao=rebuild_ao
  IF (PRESENT(rebuild_grids)) my_rebuild_grids=rebuild_grids
  IF (PRESENT(gapw_xc)) my_gapw_xc=gapw_xc

  CALL get_qs_env(qs_env, pw_env=pw_env,dft_control=dft_control,&
       matrix_s=matrix_s,error=error)
  CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
  nspins=dft_control%nspins

  IF (.NOT.ASSOCIATED(rho)) THEN
     my_rebuild_grids=.TRUE.
     my_rebuild_ao=.TRUE.
     ALLOCATE(rho,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.NOT.failure) THEN
        rho%ref_count=1
        rho%rho_g_valid=.FALSE.
        rho%rho_r_valid=.FALSE.
        rho%tau_g_valid=.FALSE.
        rho%tau_r_valid=.FALSE.
        rho%soft_valid = .FALSE.
        last_rho_id_nr=last_rho_id_nr+1
        rho%id_nr=last_rho_id_nr
        rho%rebuild_each=5
        NULLIFY(rho%rho_r, rho%rho_g, rho%rho_ao, rho%tot_rho_r, rho%tot_rho_g, rho%tau_r, rho%tau_g)
        ALLOCATE(rho%tot_rho_r(nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (.not.failure) rho%tot_rho_r=0.0_dp
     END IF
  END IF

  ! rho_ao
  IF(.NOT. my_gapw_xc) THEN
    IF (my_rebuild_ao.OR..NOT.ASSOCIATED(rho%rho_ao)) THEN
       IF (.NOT.failure) THEN
          IF (ASSOCIATED(rho%rho_ao)) THEN
             CALL deallocate_matrix_set(rho%rho_ao,error=error)
          END IF
          CALL allocate_matrix_set(rho%rho_ao,nspins,error=error)
          DO i=1,nspins
             NULLIFY (rho%rho_ao(i)%matrix)
             CALL replicate_matrix_structure(matrix_s(1)%matrix,&
                rho%rho_ao(i)%matrix,&
                  "DENSITY_MATRIX_FOR_SPIN_"//TRIM(ADJUSTL(cp_to_string(i))),error=error)
          END DO
       END IF
     END IF
  END IF

  ! rho_r
  IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%rho_r)) THEN
     IF (.NOT.failure) THEN
        IF (ASSOCIATED(rho%rho_r)) THEN
           DO i=1,SIZE(rho%rho_r)
              CALL pw_release(rho%rho_r(i)%pw,error=error)
           END DO
           DEALLOCATE(rho%rho_r,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        ALLOCATE(rho%rho_r(nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (.NOT.failure) THEN
        DO i=1,nspins
           CALL pw_pool_create_pw(auxbas_pw_pool,rho%rho_r(i)%pw,&
                use_data=REALDATA3D,in_space=REALSPACE,error=error)
        END DO
     END IF
  END IF

  ! rho_g
  IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%rho_g)) THEN
     IF (.NOT.failure) THEN
        IF (ASSOCIATED(rho%rho_g)) THEN
           DO i=1,SIZE(rho%rho_g)
              CALL pw_release(rho%rho_g(i)%pw,error=error)
           END DO
           DEALLOCATE(rho%rho_g,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        ALLOCATE(rho%rho_g(nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF (.NOT.failure) THEN
        DO i=1,nspins
           CALL pw_pool_create_pw(auxbas_pw_pool,rho%rho_g(i)%pw,&
                use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
        END DO
     END IF
  END IF

  ! allocate tau_r and tau_g if use_kinetic_energy_density
  IF (dft_control%use_kinetic_energy_density) THEN
     ! tau_r
     IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%tau_r)) THEN
        IF (.NOT.failure) THEN
           IF (ASSOCIATED(rho%tau_r)) THEN
              DO i=1,SIZE(rho%tau_r)
                 CALL pw_release(rho%tau_r(i)%pw,error=error)
              END DO
              DEALLOCATE(rho%tau_r,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           ALLOCATE(rho%tau_r(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           DO i=1,nspins
              CALL pw_pool_create_pw(auxbas_pw_pool,rho%tau_r(i)%pw,&
                   use_data=REALDATA3D,in_space=REALSPACE,error=error)
           END DO
        END IF
     END IF

     ! tau_g
     IF (my_rebuild_grids.OR..NOT.ASSOCIATED(rho%tau_g)) THEN
        IF (.NOT.failure) THEN
           IF (ASSOCIATED(rho%tau_g)) THEN
              DO i=1,SIZE(rho%tau_g)
                 CALL pw_release(rho%tau_g(i)%pw,error=error)
              END DO
              DEALLOCATE(rho%tau_g,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           ALLOCATE(rho%tau_g(nspins),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF (.NOT.failure) THEN
           DO i=1,nspins
              CALL pw_pool_create_pw(auxbas_pw_pool,rho%tau_g(i)%pw,&
                   use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
           END DO
        END IF
     END IF
  END IF  ! use_kinetic_energy_density

  CALL timestop(0.0_dp,handle)
END SUBROUTINE qs_rho_rebuild

! *****************************************************************************
!> \brief updates rho_r and rho_g to the rho%rho_ao.
!>      if use_kinetic_energy_density also computes tau_r and tau_g
!> \param rho_struct the rho structure that should be updated
!> \param qs_env the qs_env rho_struct refers to
!>        the integrated charge in r space
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_rho_update_rho(rho_struct, qs_env, error)
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_update_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nspins
    LOGICAL                                  :: failure, gapw, gapw_xc
    REAL(KIND=dp)                            :: dum
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho_xc

    failure=.FALSE.

    NULLIFY(dft_control,rho_xc)

    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

      CALL get_qs_env(qs_env, dft_control=dft_control,error=error)
      nspins=dft_control%nspins
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
    END IF
    IF (.NOT.failure) THEN
      ! get tmp grids
      IF(dft_control%qs_control%semi_empirical) THEN
        !
      ELSEIF(dft_control%qs_control%dftb) THEN
        !
      ELSE
        DO ispin=1,nspins
           CALL calculate_rho_elec(matrix_p=rho_struct%rho_ao(ispin)%matrix,&
                rho=rho_struct%rho_r(ispin),&
                rho_gspace=rho_struct%rho_g(ispin),&
                total_rho=rho_struct%tot_rho_r(ispin),&
                qs_env=qs_env,soft_valid=gapw,error=error)
        END DO
      END IF
      rho_struct%rho_r_valid=.TRUE.
      rho_struct%rho_g_valid=.TRUE.

      ! if needed compute also the kinetic energy density
      IF (dft_control%use_kinetic_energy_density) THEN
         DO ispin=1,nspins
            CALL calculate_rho_elec(matrix_p=rho_struct%rho_ao(ispin)%matrix,&
                 rho=rho_struct%tau_r(ispin),&
                 rho_gspace=rho_struct%tau_g(ispin),&
                 total_rho=dum, &  ! presumably not meaningful
                 qs_env=qs_env, soft_valid=gapw, &
                 compute_tau=.TRUE., error=error)
         END DO
         rho_struct%tau_r_valid=.TRUE.
         rho_struct%tau_g_valid=.TRUE.
      ENDIF
!FM      call cp_sm_write(matrix=rho_struct%rho_ao(1)%matrix, unit_nr=122, &
!FM           long_description=.TRUE., local=.true., error=error)
!FM      WRITE (123,"(5F30.20)") rho_struct%rho_r(1)%pw%cr3d
!FM      stop "pippo"

     ! GAPW o GAPW_XC require the calculation of hard and soft local densities
      IF(gapw ) THEN
         CALL calculate_rho_atom_coeff(qs_env,rho_struct%rho_ao,error=error)
      ENDIF
      IF(gapw_xc ) THEN
         CALL get_qs_env(qs_env=qs_env,rho_xc=rho_xc,error=error)
         CALL calculate_rho_atom_coeff(qs_env,rho_struct%rho_ao,error=error)
         DO ispin=1,nspins
           CALL calculate_rho_elec(matrix_p=rho_xc%rho_ao(ispin)%matrix,&
                rho=rho_xc%rho_r(ispin),&
                rho_gspace=rho_xc%rho_g(ispin),&
                total_rho=rho_xc%tot_rho_r(ispin),&
                qs_env=qs_env,soft_valid=gapw_xc,error=error)
         END DO
         rho_xc%rho_r_valid=.TRUE.
         rho_xc%rho_g_valid=.TRUE.
         ! if needed compute also the kinetic energy density
         IF (dft_control%use_kinetic_energy_density) THEN
           DO ispin=1,nspins
             CALL calculate_rho_elec(matrix_p=rho_xc%rho_ao(ispin)%matrix,&
                 rho=rho_xc%tau_r(ispin),&
                 rho_gspace=rho_xc%tau_g(ispin),&
                 total_rho=dum, &  ! presumably not meaningful
                 qs_env=qs_env, soft_valid=gapw_xc, &
                 compute_tau=.TRUE., error=error)
           END DO
           rho_xc%tau_r_valid=.TRUE.
           rho_xc%tau_g_valid=.TRUE.
         ENDIF
      ENDIF
    END IF

END SUBROUTINE qs_rho_update_rho

! *****************************************************************************
!> \brief Duplicates a pointer physically
!> \param rho_input The rho structure to be duplicated
!> \param rho_output The duplicate rho structure
!> \param qs_env The QS environment from which the auxilliary PW basis-set
!>                pool is taken
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \note
!>      Associated pointers are deallocated, nullified pointers are NOT accepted!
!> \par History
!>      07.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE duplicate_rho_type(rho_input, rho_output, qs_env, error)

    TYPE(qs_rho_type), POINTER               :: rho_input, rho_output
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'duplicate_rho_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins, stat
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(dft_control, pw_env, auxbas_pw_pool)

    CPPrecondition(ASSOCIATED(rho_input), cp_failure_level, routineP, error, failure)
    CPPrecondition(rho_input%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, &
                      error=error)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      nspins = dft_control%nspins

      IF (ASSOCIATED(rho_output)) THEN
        CALL qs_rho_release(rho_output, error)
      END IF
      ALLOCATE(rho_output, stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      NULLIFY(rho_output%rho_ao, rho_output%rho_g, rho_output%rho_r, &
              rho_output%tau_g, rho_output%tau_r, rho_output%tot_rho_r, &
              rho_output%tot_rho_g)

      ! rho_ao
      IF (ASSOCIATED(rho_input%rho_ao)) THEN
        !CALL allocate_matrix_set(rho_output%rho_ao, SIZE(rho_input%rho_ao))
        CALL allocate_matrix_set(rho_output%rho_ao, nspins,error=error)
        !DO i = 1, SIZE(rho_input%rho_ao)
        DO i = 1, nspins
          CALL replicate_matrix(source=rho_input%rho_ao(i)%matrix, &
                                TARGET=rho_output%rho_ao(i)%matrix, &
                                target_name="myDensityMatrix_for_Spin_"//&
                                TRIM(ADJUSTL(cp_to_string(i))), &
                                allocate_blocks=.TRUE.,error=error)
        END DO
      END IF
    END IF

    ! rho_r
    IF (.NOT. failure) THEN
      IF (ASSOCIATED(rho_input%rho_r)) THEN
        ALLOCATE(rho_output%rho_r(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_r(i)%pw, &
                                    use_data=REALDATA3D, in_space=REALSPACE, &
                                    error=error)
            rho_output%rho_r(i)%pw%cr3d(:,:,:) = rho_input%rho_r(i)%pw%cr3d(:,:,:)
          END DO
        END IF
      END IF
    END IF

    ! rho_g
    IF (.NOT. failure) THEN
      IF (ASSOCIATED(rho_input%rho_g)) THEN
        ALLOCATE(rho_output%rho_g(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_g(i)%pw, &
                                    use_data=COMPLEXDATA1D, &
                                    in_space=RECIPROCALSPACE, error=error)
            rho_output%rho_g(i)%pw%cc(:)      = rho_input%rho_g(i)%pw%cc(:)
          END DO
        END IF
      END IF
    END IF

    ! tau_r and tau_g are only needed in the case of Meta-GGA XC-functionals
    ! are used. Therefore they are only allocated if
    ! dft_control%use_kinetic_energy_density is true
    IF (.NOT.failure) THEN
      IF (dft_control%use_kinetic_energy_density) THEN
        ! tau_r
        IF (ASSOCIATED(rho_input%tau_r)) THEN
          ALLOCATE(rho_output%tau_r(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_r(i)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE, &
                                      error=error)
              rho_output%tau_r(i)%pw%cr3d(:,:,:) = rho_input%tau_r(i)%pw%cr3d(:,:,:)
            END DO
          END IF
        END IF

        ! tau_g
        IF (ASSOCIATED(rho_input%tau_g)) THEN
          ALLOCATE(rho_output%tau_g(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_g(i)%pw, &
                                      use_data=COMPLEXDATA1D, &
                                      in_space=RECIPROCALSPACE, error=error)
              rho_output%tau_g(i)%pw%cc(:) = rho_input%tau_g(i)%pw%cc(:)
            END DO
          END IF
        END IF
      END IF

      rho_output%rho_g_valid = rho_input%rho_g_valid
      rho_output%rho_r_valid = rho_input%rho_r_valid
      rho_output%tau_r_valid = rho_input%tau_r_valid
      rho_output%tau_g_valid = rho_input%tau_g_valid
      rho_output%soft_valid  = rho_input%soft_valid

      rho_output%ref_count    = 1
      last_rho_id_nr          = last_rho_id_nr+1
      rho_output%id_nr        = last_rho_id_nr
      rho_output%rebuild_each = rho_input%rebuild_each

      ! tot_rho_r
      IF (ASSOCIATED(rho_input%tot_rho_r)) THEN
        ALLOCATE(rho_output%tot_rho_r(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
            rho_output%tot_rho_r(i) = rho_input%tot_rho_r(i)
          END DO
        END IF
      END IF

      ! tot_rho_g
      IF (ASSOCIATED(rho_input%tot_rho_g)) THEN
        ALLOCATE(rho_output%tot_rho_g(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
             rho_output%tot_rho_g(i) = rho_input%tot_rho_g(i)
          END DO
        END IF
      END IF
    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE duplicate_rho_type

! *****************************************************************************
!> \brief Creates a qs_rho_type, which is equivalent to
!>      alpha*rho_input1 - beta*rho_input2
!> \param alpha & beta: Multiplier for the rho input structures
!> \param rho_input 1 & rho_input2: The rho input structure
!> \param rho_output The difference rho structure
!> \param qs_env The QS environment from which the auxilliary PW basis-set
!>                pool is taken
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \note
!>      Associated pointers are deallocated, nullified pointers are NOT accepted!
!>      Inputs and output must not be the same!
!> \par History
!>      09.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE diff_rho_type(alpha, rho_input1, beta, rho_input2, rho_output, &
                           qs_env, error)

    REAL(KIND=dp), OPTIONAL                  :: alpha
    TYPE(qs_rho_type), POINTER               :: rho_input1
    REAL(KIND=dp), OPTIONAL                  :: beta
    TYPE(qs_rho_type), POINTER               :: rho_input2, rho_output
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'diff_rho_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(dft_control, pw_env, auxbas_pw_pool)

    CPPrecondition(ASSOCIATED(rho_input1), cp_failure_level, routineP, error, failure)
    CPPrecondition(rho_input1%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(rho_input2), cp_failure_level, routineP, error, failure)
    CPPrecondition(rho_input2%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (PRESENT(alpha)) THEN
      a = alpha
    ELSE
      a = 1.0_dp
    ENDIF

    IF (PRESENT(beta)) THEN
      b= beta
    ELSE
      b = -1.0_dp
    END IF

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, &
                      error=error)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      nspins = dft_control%nspins

      IF (ASSOCIATED(rho_output)) THEN
        CALL qs_rho_release(rho_output, error)
      END IF
      ALLOCATE(rho_output, stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    END IF

    IF (.NOT. failure) THEN
      NULLIFY(rho_output%rho_ao, rho_output%rho_g, rho_output%rho_r, &
              rho_output%tau_g, rho_output%tau_r, rho_output%tot_rho_r, &
              rho_output%tot_rho_g)

      ! rho_ao
      IF (ASSOCIATED(rho_input1%rho_ao) .AND. ASSOCIATED(rho_input2%rho_ao)) THEN
        !CALL allocate_matrix_set(rho_output%rho_ao, SIZE(rho_input%rho_ao))
        CALL allocate_matrix_set(rho_output%rho_ao, nspins,error=error)
        !DO i = 1, SIZE(rho_input%rho_ao)
        DO i = 1, nspins
          CALL replicate_matrix(source=rho_input1%rho_ao(i)%matrix, &
                                TARGET=rho_output%rho_ao(i)%matrix, &
                                target_name="myDensityMatrix_for_Spin_"//&
                                TRIM(ADJUSTL(cp_to_string(i))), &
                                allocate_blocks=.TRUE.,error=error)
          CALL add_matrices(c=rho_output%rho_ao(i)%matrix, alpha=a, &
                            a=rho_input1%rho_ao(i)%matrix, beta=b, &
                            b=rho_input2%rho_ao(i)%matrix)
        END DO
      END IF
    END IF

    ! rho_r
    IF (.NOT. failure) THEN
      IF (ASSOCIATED(rho_input1%rho_r) .AND. ASSOCIATED(rho_input2%rho_r)) THEN
        ALLOCATE(rho_output%rho_r(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_r(i)%pw, &
                                    use_data=REALDATA3D, in_space=REALSPACE, &
                                    error=error)
            rho_output%rho_r(i)%pw%cr3d(:,:,:) = a*rho_input1%rho_r(i)%pw%cr3d(:,:,:) - &
                                                 b*rho_input2%rho_r(i)%pw%cr3d(:,:,:)
          END DO
        END IF
      END IF
    END IF

    ! rho_g
    IF (.NOT. failure) THEN
      IF (ASSOCIATED(rho_input1%rho_g) .AND. ASSOCIATED(rho_input2%rho_g)) THEN
        ALLOCATE(rho_output%rho_g(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_g(i)%pw, &
                                    use_data=COMPLEXDATA1D, &
                                    in_space=RECIPROCALSPACE, error=error)
            rho_output%rho_g(i)%pw%cc(:)      = a*rho_input1%rho_g(i)%pw%cc(:) - &
                                                b*rho_input2%rho_g(i)%pw%cc(:)
          END DO
        END IF
      END IF
    END IF

    ! tau_r and tau_g are only needed in the case of Meta-GGA XC-functionals
    ! are used. Therefore they are only allocated if
    ! dft_control%use_kinetic_energy_density is true
    IF (.NOT.failure) THEN
      IF (dft_control%use_kinetic_energy_density) THEN
        ! tau_r
        IF (ASSOCIATED(rho_input1%tau_r) .AND. ASSOCIATED(rho_input2%tau_r)) THEN
          ALLOCATE(rho_output%tau_r(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_r(i)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE, &
                                      error=error)
              rho_output%tau_r(i)%pw%cr3d(:,:,:) = a*rho_input1%tau_r(i)%pw%cr3d(:,:,:) - &
                                                   b*rho_input2%tau_r(i)%pw%cr3d(:,:,:)
            END DO
          END IF
        END IF

        ! tau_g
        IF (ASSOCIATED(rho_input1%tau_g) .AND. ASSOCIATED(rho_input2%tau_g)) THEN
          ALLOCATE(rho_output%tau_g(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_g(i)%pw, &
                                      use_data=COMPLEXDATA1D, &
                                      in_space=RECIPROCALSPACE, error=error)
              rho_output%tau_g(i)%pw%cc(:) = a*rho_input1%tau_g(i)%pw%cc(:) - &
                                             b*rho_input2%tau_g(i)%pw%cc(:)
            END DO
          END IF
        END IF
      END IF

      rho_output%rho_g_valid = rho_input1%rho_g_valid
      rho_output%rho_r_valid = rho_input1%rho_r_valid
      rho_output%tau_r_valid = rho_input1%tau_r_valid
      rho_output%tau_g_valid = rho_input1%tau_g_valid
      rho_output%soft_valid  = rho_input1%soft_valid

      rho_output%ref_count    = 1
      last_rho_id_nr          = last_rho_id_nr+1
      rho_output%id_nr        = last_rho_id_nr

      rho_output%rebuild_each = rho_input1%rebuild_each

      ! tot_rho_r
      IF (ASSOCIATED(rho_input1%tot_rho_r) .AND. ASSOCIATED(rho_input2%tot_rho_r)) THEN
        ALLOCATE(rho_output%tot_rho_r(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
            rho_output%tot_rho_r(i) = a*rho_input1%tot_rho_r(i) - &
                                      b*rho_input2%tot_rho_r(i)
          END DO
        END IF
      END IF

      ! tot_rho_g
      IF (ASSOCIATED(rho_input1%tot_rho_g) .AND. ASSOCIATED(rho_input2%tot_rho_g)) THEN
        ALLOCATE(rho_output%tot_rho_g(nspins), stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
        IF (.NOT. failure) THEN
          DO i = 1, nspins
             rho_output%tot_rho_g(i) = a*rho_input1%tot_rho_g(i) - &
                                       b*rho_input2%tot_rho_g(i)
          END DO
        END IF
      END IF
    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE diff_rho_type

! *****************************************************************************
!> \brief Copies the whole rho structure
!> \param rho_input The rho structure to be copied
!> \param rho_output The copy rho structure
!> \param qs_env The QS environment from which the auxilliary PW basis-set
!>                pool is taken
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      09.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE copy_rho_type(rho_input, rho_output, qs_env, error)

    TYPE(qs_rho_type), POINTER               :: rho_input, rho_output
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_rho_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nspins, stat
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(dft_control, pw_env, auxbas_pw_pool)

    CPPrecondition(ASSOCIATED(rho_input), cp_failure_level, routineP, error, failure)
    CPPrecondition(rho_input%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control, &
                      error=error)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      nspins = dft_control%nspins

      IF (.NOT.  ASSOCIATED(rho_output)) THEN
        ALLOCATE(rho_output, stat=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      END IF
    END IF

    IF (.NOT. failure) THEN

      ! rho_ao
      IF (ASSOCIATED(rho_input%rho_ao)) THEN
        IF (.NOT. ASSOCIATED(rho_output%rho_ao)) THEN
          CALL allocate_matrix_set(rho_output%rho_ao, nspins,error=error)
          DO i = 1, nspins
            CALL replicate_matrix(source=rho_input%rho_ao(i)%matrix, &
                                  TARGET=rho_output%rho_ao(i)%matrix, &
                                  target_name="myDensityMatrix_for_Spin_"//&
                                  TRIM(ADJUSTL(cp_to_string(i))), &
                                  allocate_blocks=.TRUE.,error=error)
          END DO
        ELSE
          DO i = 1,nspins
            CALL copy_matrix(source=rho_input%rho_ao(i)%matrix, &
                             TARGET=rho_output%rho_ao(i)%matrix)
            !CALL transfer_matrix(source_matrix=rho_input%rho_ao(i)%matrix, &
            !                     target_matrix=rho_output%rho_ao(i)%matrix, &
            !                     error=error)
          END DO
        END IF
      END IF
    END IF

    ! rho_r
    IF (.NOT. failure) THEN
      IF (ASSOCIATED(rho_input%rho_r)) THEN
        IF (.NOT. ASSOCIATED(rho_output%rho_r)) THEN
          ALLOCATE(rho_output%rho_r(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_r(i)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE, &
                                      error=error)
              rho_output%rho_r(i)%pw%cr3d(:,:,:) = rho_input%rho_r(i)%pw%cr3d(:,:,:)
            END DO
          END IF
        ELSE
          DO i = 1,nspins
            rho_output%rho_r(i)%pw%cr3d(:,:,:) = rho_input%rho_r(i)%pw%cr3d(:,:,:)
          END DO
        END IF
      END IF
    END IF

    ! rho_g
    IF (.NOT. failure) THEN
      IF (ASSOCIATED(rho_input%rho_g)) THEN
        IF (.NOT. ASSOCIATED(rho_output%rho_g)) THEN
          ALLOCATE(rho_output%rho_g(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%rho_g(i)%pw, &
                                      use_data=COMPLEXDATA1D, &
                                      in_space=RECIPROCALSPACE, error=error)
              rho_output%rho_g(i)%pw%cc(:)      = rho_input%rho_g(i)%pw%cc(:)
            END DO
          ELSE
            DO i = 1,nspins
              rho_output%rho_g(i)%pw%cc(:)      = rho_input%rho_g(i)%pw%cc(:)
            END DO
          END IF
        END IF
      END IF
    END IF

    ! tau_r and tau_g are only needed in the case of Meta-GGA XC-functionals
    ! are used. Therefore they are only allocated if
    ! dft_control%use_kinetic_energy_density is true
    IF (.NOT.failure) THEN
      IF (dft_control%use_kinetic_energy_density) THEN
        ! tau_r
        IF (ASSOCIATED(rho_input%tau_r)) THEN
          IF (.NOT. ASSOCIATED(rho_output%tau_r)) THEN
            ALLOCATE(rho_output%tau_r(nspins), stat=stat)
            CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
            IF (.NOT. failure) THEN
              DO i = 1, nspins
                CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_r(i)%pw, &
                                        use_data=REALDATA3D, in_space=REALSPACE, &
                                        error=error)
                rho_output%tau_r(i)%pw%cr3d(:,:,:) = rho_input%tau_r(i)%pw%cr3d(:,:,:)
              END DO
            ELSE
              DO i = 1,nspins
                rho_output%rho_g(i)%pw%cc(:)      = rho_input%rho_g(i)%pw%cc(:)
              END DO
            END IF
          END IF
        END IF

        ! tau_g
        IF (ASSOCIATED(rho_input%tau_g)) THEN
          IF (.NOT. ASSOCIATED(rho_output%tau_g)) THEN
            ALLOCATE(rho_output%tau_g(nspins), stat=stat)
            CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
            IF (.NOT. failure) THEN
              DO i = 1, nspins
                CALL pw_pool_create_pw(auxbas_pw_pool, rho_output%tau_g(i)%pw, &
                                        use_data=COMPLEXDATA1D, &
                                        in_space=RECIPROCALSPACE, error=error)
                rho_output%tau_g(i)%pw%cc(:) = rho_input%tau_g(i)%pw%cc(:)
              END DO
            ELSE
              DO i = 1,nspins
                rho_output%rho_g(i)%pw%cc(:)      = rho_input%rho_g(i)%pw%cc(:)
              END DO
            END IF
          END IF
        END IF
      END IF

      rho_output%rho_g_valid = rho_input%rho_g_valid
      rho_output%rho_r_valid = rho_input%rho_r_valid
      rho_output%tau_r_valid = rho_input%tau_r_valid
      rho_output%tau_g_valid = rho_input%tau_g_valid
      rho_output%soft_valid  = rho_input%soft_valid

      rho_output%ref_count    = 1
      last_rho_id_nr          = last_rho_id_nr+1
      rho_output%id_nr        = last_rho_id_nr
      rho_output%rebuild_each = rho_input%rebuild_each

      ! tot_rho_r
      IF (ASSOCIATED(rho_input%tot_rho_r)) THEN
        IF (.NOT. ASSOCIATED(rho_output%tot_rho_r)) THEN
          ALLOCATE(rho_output%tot_rho_r(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              rho_output%tot_rho_r(i) = rho_input%tot_rho_r(i)
            END DO
          END IF
        ELSE
          DO i = 1,nspins
            rho_output%tot_rho_r(i) = rho_input%tot_rho_r(i)
          END DO
        END IF
      END IF

      ! tot_rho_g
      IF (ASSOCIATED(rho_input%tot_rho_g)) THEN
        IF (.NOT. ASSOCIATED(rho_output%tot_rho_g)) THEN
          ALLOCATE(rho_output%tot_rho_g(nspins), stat=stat)
          CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
          IF (.NOT. failure) THEN
            DO i = 1, nspins
              rho_output%tot_rho_g(i) = rho_input%tot_rho_g(i)
            END DO
          END IF
        ELSE
          DO i = 1,nspins
            rho_output%tot_rho_g(i) = rho_input%tot_rho_g(i)
          END DO
        END IF
      END IF
    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE copy_rho_type

END MODULE qs_rho_methods

