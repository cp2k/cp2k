!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Driver for the localization that should be general
!>      for all the methods available and all the definition of the
!>      spread functional
!>      Write centers, spread and cubes only if required and for the
!>      selected states
!>      The localized functions are copied in the standard mos array
!>      for the next use
!> \par History
!>      01.2008 Teodoro Laino [tlaino] - University of Zurich
!>        - Merging the two localization codes and updating to new structures
!> \author MI (04.2005)
! *****************************************************************************
MODULE qs_loc_methods
  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             deallocate_atomic_kind_set,&
                                             set_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_get_submatrix, &
       cp_fm_p_type, cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, &
       cp_fm_to_fm, cp_fm_type
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_loc_crazy, do_loc_direct, do_loc_jacobi, do_loc_l1_norm_sd, &
       do_loc_none, dump_dcd, dump_dcd_aligned_cell, dump_xmol, op_loc_berry, &
       op_loc_boys, op_loc_pipek, state_loc_list
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_get_ivals,&
                                             section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             sp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: twopi
  USE memory_utilities,                ONLY: reallocate
  USE motion_utils,                    ONLY: get_output_format
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             deallocate_particle_set,&
                                             get_particle_set,&
                                             particle_type,&
                                             write_particle_coordinates
  USE physcon,                         ONLY: angstrom
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                             localized_wfn_control_type,&
                                             qs_loc_env_new_type
  USE qs_loc_utils,                    ONLY: jacobi_rotation_pipek
  USE qs_localization_methods,         ONLY: approx_l1_norm_sd,&
                                             crazy_rotations,&
                                             direct_mini,&
                                             jacobi_rotations
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE string_utilities,                ONLY: xstring
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_methods'

! *** Public ***
  PUBLIC ::  qs_loc_driver, qs_print_cubes

CONTAINS

! *****************************************************************************
!> \brief Calculate and optimize the spread functional as calculated from
!>       the selected mos  and by the definition using the berry phase
!>       as given by silvestrelli et al
!>       If required the centers and the spreads for each mos selected
!>       are calculated from z_ij and printed to file.
!>       The centers files is appended if already exists
!> \param method indicates localization algorithm
!> \param qs_loc_env new environment for the localization calculations
!> \param vectors selected mos to be localized
!> \param op_sm_set sparse matrices containing the integrals of the kind <mi e{iGr} nu>
!> \param zij_fm_set set of full matrix of size nmoloc x nmoloc, will contain the z_ij numbers
!>                    as defined by Silvestrelli et al
!> \param ionode logical that indicates the node responsable for the i/o
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>       This definition need the use of complex numbers, therefore the
!>       optimization routines are specific for this case
!>       The file for the centers and the spreads have a xyz format
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE optimize_loc_berry( method, qs_loc_env, vectors, op_sm_set, &
       zij_fm_set, para_env, cell, weights,  ispin, loc_section, print_loc_section,&
       error)

    INTEGER, INTENT(IN)                      :: method
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(cp_fm_type), POINTER                :: vectors
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: op_sm_set
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: zij_fm_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:)                   :: weights
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(section_vals_type), POINTER         :: loc_section, print_loc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'optimize_loc_berry', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, max_iter, nao, &
                                                nmoloc, out_each, &
                                                output_unit, sweeps
    LOGICAL                                  :: converged, crazy_use_diag
    REAL(dp)                                 :: crazy_scale, &
                                                eps_localization, &
                                                max_crazy_angle, start_time, &
                                                target_time
    TYPE(cp_fm_type), POINTER                :: opvec
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,print_loc_section,"PROGRAM_RUN_INFO",&
         extension=".locInfo",error=error)

    NULLIFY(opvec)
    ! get rows and cols of the input
    CALL cp_fm_get_info(vectors,nrow_global=nao,ncol_global=nmoloc,error=error)
    ! replicate the input kind of matrix
    CALL cp_fm_create(opvec,vectors%matrix_struct,error=error)

    ! Compute zij here
    DO i=1,SIZE ( zij_fm_set, 2 )
       DO j=1,SIZE ( zij_fm_set, 1 )
          CALL cp_fm_set_all(zij_fm_set(j,i)%matrix, 0.0_dp, error=error)
          CALL cp_dbcsr_sm_fm_multiply(op_sm_set(j,i)%matrix, vectors,opvec,ncol=nmoloc ,error=error)
          CALL cp_fm_gemm("T","N",nmoloc,nmoloc,nao,1.0_dp,vectors,opvec,0.0_dp,&
               zij_fm_set(j,i)%matrix,error=error)
       ENDDO
    ENDDO

    max_iter = qs_loc_env%localized_wfn_control%max_iter
    max_crazy_angle = qs_loc_env%localized_wfn_control%max_crazy_angle
    crazy_use_diag = qs_loc_env%localized_wfn_control%crazy_use_diag
    crazy_scale = qs_loc_env%localized_wfn_control%crazy_scale
    eps_localization = qs_loc_env%localized_wfn_control%eps_localization
    out_each = qs_loc_env%localized_wfn_control%out_each
    target_time = qs_loc_env%target_time
    start_time = qs_loc_env%start_time
    CALL centers_spreads_berry( qs_loc_env, zij_fm_set, nmoloc, cell, weights,&
         ispin, loc_section, print_loc_section, only_initial_out=.TRUE., error=error )
    SELECT CASE ( method )
    CASE ( do_loc_jacobi )
        CALL jacobi_rotations ( weights, zij_fm_set, vectors, para_env, max_iter=max_iter,&
            eps_localization=eps_localization, sweeps=sweeps, &
            out_each=out_each, target_time=target_time, start_time=start_time, error=error )
    CASE ( do_loc_crazy )
       CALL crazy_rotations( weights, zij_fm_set, vectors, max_iter=max_iter, max_crazy_angle = max_crazy_angle, &
            crazy_scale=crazy_scale, crazy_use_diag=crazy_use_diag, &
            eps_localization=eps_localization, iterations= sweeps ,converged=converged,&
            error=error)
       ! Possibly fallback to jacobi if the crazy rotation fails
       IF ( .NOT.converged ) THEN
          IF (qs_loc_env%localized_wfn_control%jacobi_fallback) THEN
             IF (output_unit>0) WRITE (output_unit,"(A,I6,A)")&
                  " Crazy Wannier localization not converged after ",sweeps,&
                  " iterations, switching to jacobi rotations"
             CALL jacobi_rotations ( weights, zij_fm_set, vectors, para_env, max_iter=max_iter,&
                 eps_localization=eps_localization, sweeps=sweeps, &
                 out_each=out_each, target_time=target_time, start_time=start_time, error=error )
          ELSE
             IF (output_unit>0) WRITE (output_unit,"(A,I6,A)")&
                  " Crazy Wannier localization not converged after ",sweeps,&
                  " iterations, and jacobi_fallback switched off "
          ENDIF
       END IF
    CASE ( do_loc_direct )
       CALL direct_mini( weights, zij_fm_set, vectors, max_iter=max_iter,&
            eps_localization=eps_localization, iterations=sweeps,  out_each=out_each, error=error  )
    CASE( do_loc_l1_norm_sd )
       IF(.NOT. cell%orthorhombic ) THEN
           CALL stop_program(routineN,moduleN,__LINE__,&
                             "Non-orthorhombic cell with the selected method NYI")
       ELSE
       CALL approx_l1_norm_sd(vectors,max_iter,eps_localization,converged,sweeps,error=error)
       ! here we need to set zij for the computation of the centers and spreads
       DO i=1,SIZE ( zij_fm_set, 2 )
          DO j=1,SIZE ( zij_fm_set, 1 )
             CALL cp_fm_set_all(zij_fm_set(j,i)%matrix, 0.0_dp, error=error)
             CALL cp_dbcsr_sm_fm_multiply(op_sm_set(j,i)%matrix, vectors,opvec,ncol=nmoloc ,error=error)
             CALL cp_fm_gemm("T","N",nmoloc,nmoloc,nao,1.0_dp,vectors,opvec,0.0_dp,&
                  zij_fm_set(j,i)%matrix,error=error)
          ENDDO
       ENDDO
       END IF
    CASE ( do_loc_none )
       IF (output_unit>0) THEN
          WRITE(output_unit,'(A,I6,A)') " No MOS localization applied "
       ENDIF
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown localization method")
    END SELECT
    IF (output_unit>0) THEN
       IF(sweeps<=max_iter) WRITE(output_unit,'(T15,A,I3,A,I6,A)') " Localization  for spin ",ispin ,&
            " converged in ",sweeps," iterations"
    ENDIF

    CALL centers_spreads_berry( qs_loc_env, zij_fm_set, nmoloc, cell, weights,&
         ispin, loc_section, print_loc_section, error=error )
    CALL cp_fm_release(opvec,error=error)
    CALL cp_print_key_finished_output(output_unit,logger,print_loc_section,"PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE optimize_loc_berry

! *****************************************************************************
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE  optimize_loc_pipek(qs_env, method, qs_loc_env, vectors, zij_fm_set,&
       ispin, loc_section, print_loc_section, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: method
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(cp_fm_type), POINTER                :: vectors
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: zij_fm_set
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(section_vals_type), POINTER         :: loc_section, print_loc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'optimize_loc_pipek', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, isgf, istat, &
                                                ldz, nao, natom, ncol, &
                                                nmoloc, output_unit, sweeps
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf, nsgf
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_ao_fm_pools
    TYPE(cp_fm_type), POINTER                :: opvec, ov_fm
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,print_loc_section,"PROGRAM_RUN_INFO",&
         extension=".locInfo",error=error)

    NULLIFY(particle_set)
    ! get rows and cols of the input
    CALL cp_fm_get_info(vectors,nrow_global=nao,ncol_global=nmoloc,error=error)
    ! replicate the input kind of matrix
    CALL cp_fm_create(opvec,vectors%matrix_struct,error=error)
    CALL cp_fm_set_all(opvec, 0.0_dp, error=error)

    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, particle_set=particle_set,error=error)

    natom = SIZE(particle_set,1)
    ALLOCATE (first_sgf(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (nsgf(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    !   construction of
    CALL get_particle_set(particle_set=particle_set,&
         first_sgf=first_sgf,last_sgf=last_sgf,nsgf=nsgf,error=error)

    !   Copy the overlap sparse matrix in a full matrix
    CALL mpools_get(qs_env%mpools,ao_ao_fm_pools=ao_ao_fm_pools,error=error)
    CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, ov_fm, name=" ", error=error)
    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,ov_fm,error=error)

    !   Compute zij here
    DO iatom = 1, natom
       CALL cp_fm_set_all(zij_fm_set(iatom,1)%matrix, 0.0_dp, error=error)
       CALL cp_fm_get_info(zij_fm_set(iatom,1)%matrix,ncol_global=ldz,error=error)
       isgf = first_sgf(iatom)
       ncol = nsgf(iatom)
       ! multiply fmxfm, using only part of the ao : Ct x S
       CALL cp_fm_gemm('N','N',nao,nmoloc,nao,1.0_dp,ov_fm,vectors,0.0_dp,opvec,&
            b_first_col=1,a_first_row=1,b_first_row=1,error=error)
       CALL cp_fm_gemm('T','N',nmoloc,nmoloc,ncol,0.5_dp,vectors,opvec,&
            0.0_dp,zij_fm_set(iatom,1)%matrix,&
            b_first_col=1,a_first_row=isgf,b_first_row=isgf,error=error)

       CALL cp_fm_gemm('N','N',nao,nmoloc,ncol,1.0_dp,ov_fm,vectors,0.0_dp,opvec,&
            b_first_col=1,a_first_row=isgf,b_first_row=isgf,error=error)
       CALL cp_fm_gemm('T','N',nmoloc,nmoloc,nao,0.5_dp,vectors,opvec,&
            1.0_dp,zij_fm_set(iatom,1)%matrix,&
            b_first_col=1,a_first_row=1,b_first_row=1,error=error)
    END DO  ! iatom

    !   And now perform the optimization and rotate the orbitals
    SELECT CASE ( method )
    CASE ( do_loc_jacobi )
       CALL jacobi_rotation_pipek(zij_fm_set, vectors, sweeps ,error=error)
    CASE ( do_loc_crazy )
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Crazy and Pipek not implemented.", &
            error=error, error_level=cp_failure_level)
    CASE( do_loc_l1_norm_sd )
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="L1 norm and Pipek not implemented.", &
            error=error, error_level=cp_failure_level)
    CASE ( do_loc_direct )
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Direct and Pipek not implemented.", &
            error=error, error_level=cp_failure_level)
    CASE ( do_loc_none )
       IF (output_unit>0) WRITE(output_unit,'(A,I6,A)') " No MOS localization applied "
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown localization method")
    END SELECT

    IF (output_unit>0) WRITE(output_unit,'(A,I3,A,I6,A)') " Localization  for spin ",ispin ,&
         " converged in ",sweeps," iterations"

    CALL centers_spreads_pipek( qs_loc_env, zij_fm_set, particle_set, ispin,&
          print_loc_section, error=error )

    DEALLOCATE (first_sgf, last_sgf, nsgf, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL cp_fm_release(opvec,error=error)
    CALL cp_print_key_finished_output(output_unit,logger,print_loc_section,"PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE optimize_loc_pipek

! *****************************************************************************
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE centers_spreads_berry( qs_loc_env, zij, nmoloc, cell, weights, ispin,&
       loc_section, print_loc_section, only_initial_out, error)
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: zij( :, : )
    INTEGER, INTENT(IN)                      :: nmoloc
    TYPE(cell_type), POINTER                 :: cell
    REAL(dp), DIMENSION(:)                   :: weights
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(section_vals_type), POINTER         :: loc_section, print_loc_section
    LOGICAL, INTENT(IN), OPTIONAL            :: only_initial_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'centers_spreads_berry', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: file_tmp, iter
    COMPLEX(KIND=dp)                         :: z
    INTEGER                                  :: idir, istate, jdir, nstates, &
                                                output_unit, unit_out_s
    LOGICAL                                  :: failure, my_only_init
    REAL(dp)                                 :: spread_i, spread_ii, &
                                                sum_spread_i, sum_spread_ii
    REAL(dp), DIMENSION(3)                   :: c, c2, cpbc
    REAL(dp), DIMENSION(:, :), POINTER       :: centers
    REAL(KIND=dp)                            :: imagpart, realpart
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    NULLIFY(centers, logger, print_key)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    my_only_init=.FALSE.
    IF(PRESENT(only_initial_out)) my_only_init=only_initial_out

    file_tmp = TRIM(qs_loc_env%tag_mo)//"_spreads_s"//TRIM(ADJUSTL(cp_to_string(ispin)))
    output_unit = cp_print_key_unit_nr(logger,print_loc_section,"PROGRAM_RUN_INFO",&
         extension=".locInfo",error=error)
    unit_out_s = cp_print_key_unit_nr(logger,print_loc_section,"WANNIER_SPREADS",&
         middle_name=file_tmp,extension=".data",error=error)
    iter = cp_iter_string(logger%iter_info,error=error)
    IF(unit_out_s>0 .AND. .NOT. my_only_init) WRITE(unit_out_s,'(i6,/,A)') nmoloc, TRIM(iter)

    CALL  cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstates ,error=error)
    CPPrecondition(nstates>=nmoloc,cp_failure_level,routineP,error,failure)

    centers => qs_loc_env%localized_wfn_control%centers_set(ispin)%array
    CPPrecondition(ASSOCIATED(centers),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(centers,2)==nmoloc,cp_failure_level,routineP,error,failure)
    sum_spread_i  = 0.0_dp
    sum_spread_ii  = 0.0_dp
    DO istate=1,nmoloc
       c = 0.0_dp
       c2 = 0.0_dp
       spread_i = 0.0_dp
       spread_ii = 0.0_dp
       DO jdir = 1,SIZE(zij,2)
         CALL cp_fm_get_element(zij(1,jdir)%matrix,istate,istate,realpart)
         CALL cp_fm_get_element(zij(2,jdir)%matrix,istate,istate,imagpart)
         z=CMPLX(realpart,imagpart,dp)
         spread_i =  spread_i -weights ( jdir ) *&
                     LOG (realpart*realpart+imagpart*imagpart ) / twopi/twopi
         spread_ii =  spread_ii + weights ( jdir ) *&
                     (1.0_dp - (realpart*realpart+imagpart*imagpart )) / twopi/twopi
         IF(jdir < 4) THEN
           DO idir=1,3
             c(idir) =  c(idir) + &
                    ( cell % hmat ( idir, jdir ) / twopi ) * AIMAG ( LOG ( z ) )
           ENDDO
         END IF
       END DO
       cpbc = pbc(c,cell)
       centers (1:3, istate ) = cpbc(1:3)
       centers (4, istate ) = spread_i
       centers (5, istate ) = spread_ii
       sum_spread_i = sum_spread_i +  centers (4, istate )
       sum_spread_ii = sum_spread_ii +  centers (5, istate )
       IF (unit_out_s>0 .AND. .NOT. my_only_init) WRITE(unit_out_s,'(I6,2F16.8)') istate, centers(4:5,istate)
    ENDDO

    ! Print of wannier centers
    print_key => section_vals_get_subs_vals(print_loc_section,"WANNIER_CENTERS",error=error)
    IF (.NOT. my_only_init)  CALL print_wannier_centers(qs_loc_env, print_key, centers, logger, ispin, error)

    IF (output_unit>0) THEN
        WRITE(output_unit,'(T10, A, 2x, A26, A26)') " Spread Functional ","sum_in -w_i ln(|z_in|^2)",&
                             "sum_in w_i(1-|z_in|^2)"
        IF(my_only_init) THEN
          WRITE(output_unit,'(T4,A,10x,2F20.10)') " Initial Spread (Berry) : ", sum_spread_i, sum_spread_ii
        ELSE
          WRITE(output_unit,'(T4,A,10x,2F20.10)') " Total Spread (Berry) : ", sum_spread_i, sum_spread_ii
        END IF
    END IF

    IF (unit_out_s>0 .AND. .NOT. my_only_init) WRITE (unit_out_s,'(A,2F16.10)') " Total ",  sum_spread_i, sum_spread_ii

    CALL cp_print_key_finished_output(unit_out_s,logger,print_loc_section,"WANNIER_SPREADS",error=error)
    CALL cp_print_key_finished_output(output_unit,logger,print_loc_section,"PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE centers_spreads_berry

! *****************************************************************************
!> \brief define and print the centers and spread
!>       when the pipek operator is used
!> \param zij_fm_set matrix elements that define the populations on atoms
!> \param localized_wfn_control control parameter for a localization calculation
!> \param ispin spin 1 or 2
!> \param ionode i/o proc. should write the output if required
!>       error:
!> \par History
!>      05.2005 created [MI]
! *****************************************************************************
  SUBROUTINE centers_spreads_pipek(qs_loc_env,zij_fm_set,particle_set,ispin,&
       print_loc_section,error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: zij_fm_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(section_vals_type), POINTER         :: print_loc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'centers_spreads_pipek', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: file_tmp, iter
    INTEGER                                  :: iatom, istat, istate, natom, &
                                                nstate, unit_out_s
    INTEGER, DIMENSION(:), POINTER           :: atom_of_state
    LOGICAL                                  :: failure
    REAL(dp)                                 :: r(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: Qii, ziimax
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: diag
    REAL(dp), DIMENSION(:, :), POINTER       :: centers
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key

    NULLIFY(centers,logger, print_key)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CALL cp_fm_get_info ( zij_fm_set(1,1)%matrix, nrow_global = nstate,error=error)
    natom = SIZE(zij_fm_set,1)

    centers => qs_loc_env%localized_wfn_control%centers_set(ispin)%array
    CPPrecondition(ASSOCIATED(centers),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(centers,2)==nstate,cp_failure_level,routineP,error,failure)

    file_tmp = TRIM(qs_loc_env%tag_mo)//"_spreads_s"//TRIM(ADJUSTL(cp_to_string(ispin)))
    unit_out_s = cp_print_key_unit_nr(logger,print_loc_section,"WANNIER_SPREADS",&
         middle_name=file_tmp,extension=".data",log_filename=.FALSE.,error=error)
    iter = cp_iter_string(logger%iter_info,error=error)
    IF(unit_out_s>0) WRITE(unit_out_s,'(i6,/,A)') TRIM(iter)

    ALLOCATE (atom_of_state(nstate),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    atom_of_state = 0

    ALLOCATE(diag(nstate,natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    diag = 0.0_dp

    DO iatom = 1,natom
       DO istate = 1,nstate
          CALL cp_fm_get_element(zij_fm_set(iatom,1)%matrix,istate,istate,diag(istate,iatom))
       END DO
    END DO

    ALLOCATE (Qii(nstate),ziimax(nstate),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ziimax = 0.0_dp
    Qii    = 0.0_dp

    DO iatom = 1,natom
       DO istate = 1,nstate
          Qii(istate) = Qii(istate) + diag(istate,iatom)*diag(istate,iatom)
          IF(ABS(diag(istate,iatom)) > ziimax(istate)) THEN
             ziimax(istate) = ABS(diag(istate,iatom))
             atom_of_state(istate) = iatom
          END IF
       END DO
    END DO

    DO istate = 1,nstate
       iatom = atom_of_state(istate)
       r(1:3) = particle_set(iatom)%r(1:3)
       centers(1:3,istate)=r(1:3)
       centers(4,istate)=1.0_dp/Qii(istate)
       IF (unit_out_s>0) WRITE(unit_out_s,'(I6,F16.8)') istate, angstrom*centers(4,istate)
    END DO

    ! Print the wannier centers
    print_key => section_vals_get_subs_vals(print_loc_section,"WANNIER_CENTERS",error=error)
    CALL print_wannier_centers(qs_loc_env, print_key, centers, logger, ispin, error)

    CALL cp_print_key_finished_output(unit_out_s,logger,print_loc_section,"WANNIER_SPREADS",error=error)

    DEALLOCATE(Qii,ziimax,atom_of_state,diag,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE centers_spreads_pipek

! *****************************************************************************
!> \brief set up the calculation of localized orbitals
!> \par History
!>      04.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE qs_loc_driver(qs_env,qs_loc_env,loc_section,print_loc_section,myspin,&
       ext_mo_coeff,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_section, print_loc_section
    INTEGER, INTENT(IN), OPTIONAL            :: myspin
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: ext_mo_coeff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_loc_driver', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_pos
    INTEGER :: dim_op, handle, i, imo, imoloc, ir, ispin, istat, istate, j, &
      jstate, l_spin, lb, loc_method, n_rep, nao, ncubes, nmo, nmosub, &
      s_spin, ub
    INTEGER, DIMENSION(:), POINTER           :: bounds, list, list_cubes
    LOGICAL                                  :: append_cube, failure, &
                                                list_cubes_setup
    LOGICAL, SAVE                            :: first_time = .TRUE.
    REAL(dp), DIMENSION(6)                   :: weights
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: centers, vecbuffer
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: op_sm_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: moloc_coeff
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: op_fm_set
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)
    NULLIFY(para_env, mos, dft_control, list_cubes)
    NULLIFY(cell, localized_wfn_control, moloc_coeff, op_sm_set, op_fm_set)
    failure = .FALSE.
    qs_loc_env%first_time = first_time
    qs_loc_env%target_time = qs_env%target_time
    qs_loc_env%start_time = qs_env%start_time

    CALL get_qs_loc_env(qs_loc_env=qs_loc_env, &
         localized_wfn_control=localized_wfn_control,&
         moloc_coeff=moloc_coeff,op_sm_set=op_sm_set, op_fm_set=op_fm_set,cell=cell,&
         weights=weights, dim_op=dim_op, error=error)

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control,&
         para_env=para_env, mos=mos, error=error)

    s_spin = 1
    l_spin = dft_control%nspins
    IF(PRESENT(myspin))  THEN
       s_spin = myspin
       l_spin = myspin
    END IF

    DO ispin = s_spin,l_spin
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo)
       loc_method = localized_wfn_control%localization_method
       SELECT CASE (localized_wfn_control%operator_type)
       CASE(op_loc_berry)
          ! Here we allocate op_fm_set with the RIGHT size for uks
          NULLIFY(tmp_fm_struct,mo_coeff)
          nmosub = localized_wfn_control%nloc_states(ispin)
          IF (PRESENT(ext_mo_coeff)) THEN
             CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmosub,&
                  ncol_global=nmosub,para_env=para_env,context=ext_mo_coeff%matrix_struct%context,&
                  error=error)
          ELSE
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmosub,&
                  ncol_global=nmosub,para_env=para_env,context=mo_coeff%matrix_struct%context,&
                  error=error)
          END IF
          !
          ALLOCATE ( op_fm_set( 2, dim_op), STAT = istat )
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, dim_op
             DO j = 1, SIZE ( op_fm_set, 1 )
                NULLIFY(op_fm_set(j,i)%matrix)
                CALL cp_fm_create (op_fm_set(j,i)%matrix , tmp_fm_struct ,error=error)
                CALL cp_fm_get_info(op_fm_set(j,i)%matrix, nrow_global=nmosub, &
                     error=error)
                CPPostcondition(nmo>=nmosub,cp_failure_level,routineP,error,failure)
                CALL cp_fm_set_all(op_fm_set(j,i)%matrix,0.0_dp,error=error)
             END DO
          END DO
          CALL cp_fm_struct_release(tmp_fm_struct,error=error)

          CALL optimize_loc_berry(loc_method, qs_loc_env, moloc_coeff(ispin)%matrix, op_sm_set,&
               op_fm_set, para_env, cell, weights, ispin, loc_section, print_loc_section,  error=error)

          ! Here we dealloctate op_fm_set
          IF(ASSOCIATED(op_fm_set)) THEN
             DO i=1,dim_op
                DO j=1,SIZE ( op_fm_set,1)
                   CALL cp_fm_release(op_fm_set(j,i)%matrix,error=error)
                END DO
             END DO
             DEALLOCATE(op_fm_set,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF

       CASE(op_loc_boys)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Boys localization not implemented", &
               error=error, error_level=cp_failure_level)

       CASE(op_loc_pipek)
          CALL  optimize_loc_pipek(qs_env, loc_method, qs_loc_env, moloc_coeff(ispin)%matrix,&
               op_fm_set,  ispin, loc_section, print_loc_section, error=error)

       END SELECT

       ! give back the localized orbitals
       IF (.NOT.PRESENT(ext_mo_coeff)) THEN
          CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
       END IF

       lb = localized_wfn_control%lu_bound_states(1,ispin)
       ub = localized_wfn_control%lu_bound_states(2,ispin)

       IF(localized_wfn_control%set_of_states == state_loc_list) THEN
          ALLOCATE(vecbuffer(1,nao),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          nmosub = SIZE(localized_wfn_control%loc_states,1)
          imoloc = 0
          DO i = lb,ub
             ! Get the index in the subset
             imoloc = imoloc + 1
             ! Get the index in the full set
             imo = localized_wfn_control%loc_states(i,ispin)

             CALL cp_fm_get_submatrix(moloc_coeff(ispin)%matrix,vecbuffer,1,imoloc,&
                  nao,1,transpose=.TRUE.,error=error)
             IF (PRESENT(ext_mo_coeff)) THEN
                CALL cp_fm_set_submatrix(ext_mo_coeff,vecbuffer,1,imo,&
                     nao,1,transpose=.TRUE.,error=error)
             ELSE
                CALL cp_fm_set_submatrix(mo_coeff,vecbuffer,1,imo,&
                     nao,1,transpose=.TRUE.,error=error)
             END IF
          END DO
          DEALLOCATE(vecbuffer,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       ELSE
          nmosub = localized_wfn_control%nloc_states(ispin)
          IF (PRESENT(ext_mo_coeff)) THEN
             CALL cp_fm_to_fm(moloc_coeff(ispin)%matrix,ext_mo_coeff,nmosub,1,lb)
          ELSE
             CALL cp_fm_to_fm(moloc_coeff(ispin)%matrix,mo_coeff,nmosub,1,lb)
          END IF
       END IF

       ! Write cube files if required
       IF(localized_wfn_control%print_cubes) THEN
          list_cubes_setup = .FALSE.
          NULLIFY(bounds, list, list_cubes)

         ! Provides boundaries of MOs
          CALL section_vals_val_get(print_loc_section,"WANNIER_CUBES%CUBES_LU_BOUNDS",&
               i_vals=bounds,error=error)
          ncubes = bounds(2) - bounds(1)  + 1
          IF(ncubes > 0) THEN
             list_cubes_setup =.TRUE.
             ALLOCATE( list_cubes(ncubes),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ir = 1,ncubes
                list_cubes(ir) = bounds(1) + (ir-1)
             END DO
          END IF

          ! Provides the list of MOs
          CALL section_vals_val_get(print_loc_section,"WANNIER_CUBES%CUBES_LIST",&
               n_rep_val=n_rep,error=error)
          IF (.NOT.list_cubes_setup) THEN
             ncubes = 0
             DO ir = 1,n_rep
                CALL section_vals_val_get(print_loc_section,"WANNIER_CUBES%CUBES_LIST",&
                     i_rep_val=ir,i_vals=list,error=error)
                IF(ASSOCIATED(list)) THEN
                   CALL reallocate(list_cubes,1,ncubes+ SIZE(list))
                   DO i = 1, SIZE(list)
                      list_cubes(i+ncubes) = list(i)
                   END DO
                   ncubes = ncubes + SIZE(list)
                END IF
             END DO
             IF (ncubes>0) list_cubes_setup =.TRUE.
          END IF

          ! Full list of Mos
          IF (.NOT.list_cubes_setup) THEN
             list_cubes_setup =.TRUE.
             ncubes = localized_wfn_control%nloc_states(1)
             IF(ncubes>0) THEN
                ALLOCATE(list_cubes(ncubes),STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             END IF
             DO i = 1, ncubes
                list_cubes(i) = i
             END DO
          END IF

          ncubes = SIZE(list_cubes)
          ncubes = MIN(ncubes,nmo)
          ALLOCATE(centers(6,ncubes),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, ncubes
             istate = list_cubes(i)
             DO j = 1,localized_wfn_control%nloc_states(ispin)
                jstate = localized_wfn_control%loc_states(j,ispin)
                IF(istate==jstate) THEN
                   centers(1:6,i) = localized_wfn_control%centers_set(ispin)%array(1:6,j)
                   EXIT
                ENDIF
             END DO
          END DO ! ncubes

          ! Real call for dumping the cube files
          print_key => section_vals_get_subs_vals(print_loc_section,"WANNIER_CUBES",error=error)
          append_cube= section_get_lval(print_loc_section,"WANNIER_CUBES%APPEND",error=error)
          my_pos = "REWIND"
          IF (append_cube) THEN
            my_pos = "APPEND"
          END IF

          CALL qs_print_cubes(qs_env,moloc_coeff(ispin)%matrix,ncubes,list_cubes,centers,&
               print_key,"loc"//TRIM(ADJUSTL(qs_loc_env%tag_mo)),&
               ispin=ispin,file_position=my_pos,error=error)

          DEALLOCATE(centers,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(list_cubes,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO  ! ispin
    first_time = .FALSE.
    CALL timestop(handle)
  END SUBROUTINE qs_loc_driver

! *****************************************************************************
!> \brief write the cube files for a set of selected states
!> \param qs_en v
!> \param mo_coeff set mos from which the states to be printed are extracted
!> \param nstates number of states to be printed
!> \param state_list list of the indexes of the states to be printed
!> \param centers centers and spread, all=0 if they hva not been calculated
!> \param ispi n
!> \param ionod e
!> \param root initial part of the cube file names
!> \note
!>      This routine shoul not be in this module
!> \par History
!>      08.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE qs_print_cubes(qs_env,mo_coeff,nstates,state_list,centers,&
       print_key,root,ispin,idir,state0,file_position,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    INTEGER, INTENT(IN)                      :: nstates
    INTEGER, DIMENSION(:), POINTER           :: state_list
    REAL(dp), DIMENSION(:, :), POINTER       :: centers
    TYPE(section_vals_type), POINTER         :: print_key
    CHARACTER(LEN=*)                         :: root
    INTEGER, INTENT(IN), OPTIONAL            :: ispin, idir
    INTEGER, OPTIONAL                        :: state0
    CHARACTER(LEN=default_string_length), &
      OPTIONAL                               :: file_position
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_print_cubes', &
      routineP = moduleN//':'//routineN
    CHARACTER, DIMENSION(3), PARAMETER       :: labels = (/'x','y','z'/)

    CHARACTER                                :: label
    CHARACTER(LEN=default_path_length)       :: file_tmp, filename, my_pos
    CHARACTER(LEN=default_string_length)     :: title
    INTEGER                                  :: handle, ia, ie, istate, &
                                                ivector, my_ispin, my_state0, &
                                                unit_out_c
    LOGICAL                                  :: add_idir, add_spin, failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(auxbas_pw_pool, pw_env, logger)
    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,subsys=subsys,error=error)
    CALL cp_subsys_get(subsys,particles=particles,error=error)



    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)

    CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                           use_data = REALDATA3D,&
                           in_space = REALSPACE, error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                           use_data = COMPLEXDATA1D,&
                           in_space = RECIPROCALSPACE, error=error)

    my_state0 = 0
    IF(PRESENT(state0)) my_state0 = state0

    add_spin = .FALSE.
    my_ispin = 1
    IF(PRESENT(ispin)) THEN
       add_spin=.TRUE.
       my_ispin=ispin
    END IF
    add_idir = .FALSE.
    IF(PRESENT(idir)) THEN
       add_idir = .TRUE.
       label=labels(idir)
    END IF

    my_pos = "REWIND"
    IF(PRESENT(file_position)) THEN
      my_pos = file_position
    END IF

    DO istate = 1, nstates
       ivector = state_list(istate) - my_state0
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,cell=cell,&
            dft_control=dft_control,particle_set=particle_set,pw_env=pw_env,error=error)

       CALL calculate_wavefunction(mo_coeff, ivector, wf_r, wf_g, atomic_kind_set,&
            cell, dft_control, particle_set, pw_env, error=error)

       ! Formatting the middle part of the name
       ivector = state_list(istate)
       CALL xstring (root,ia,ie)
       IF(add_idir) THEN
          filename=root(ia:ie)//"_"//label//"_w"//TRIM(ADJUSTL(cp_to_string(ivector)))
       ELSE
          filename=root(ia:ie)//"_w"//TRIM(ADJUSTL(cp_to_string(ivector)))
       END IF
       IF(add_spin) THEN
          file_tmp=filename
          CALL xstring (file_tmp,ia,ie)
          filename=file_tmp(ia:ie)//"_s"//TRIM(ADJUSTL(cp_to_string(ispin)))
       END IF

       ! Using the print_key tools to open the file with the proper name
       unit_out_c = cp_print_key_unit_nr(logger,print_key,"",middle_name=filename,&
            extension=".cube",file_position=my_pos,log_filename=.FALSE.,error=error)
       IF(SIZE(centers,1)==6) THEN
          WRITE(title,'(A7,I5.5,A1,I1.1,A1,6F10.4)') "WFN ",ivector,"_s",my_ispin," ",&
               centers(1:3,istate)*angstrom, centers(4:6,istate)*angstrom
       ELSE
          WRITE(title,'(A7,I5.5,A1,I1.1,A1,3F10.4)') "WFN ",ivector,"_s",my_ispin," ",&
               centers(1:3,istate)*angstrom
       END IF
       CALL cp_pw_to_cube(wf_r%pw,unit_out_c,title,&
            particles=particles,&
            stride=section_get_ivals(print_key,"STRIDE",error=error),error=error)
       CALL cp_print_key_finished_output(unit_out_c,logger,print_key,"",error=error)
    END DO

    CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_g%pw, error=error)
    CALL timestop(handle)
  END SUBROUTINE qs_print_cubes

! *****************************************************************************
!> \brief Prints wannier centers
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_wannier_centers(qs_loc_env, print_key, center, logger, ispin, error)
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: print_key
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_wannier_centers', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: iter, unit_str
    CHARACTER(LEN=default_string_length)     :: my_ext, my_form
    INTEGER                                  :: iunit, l, nstates
    LOGICAL                                  :: first_time, init_traj
    REAL(KIND=dp)                            :: unit_conv

    nstates = SIZE(center,2)
    my_form = "formatted"
    my_ext  = ".data"
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,first_time=first_time,&
         error=error),cp_p_file)) THEN
       ! Find out if we want to print IONS+CENTERS or ONLY CENTERS..
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,"/IONS+CENTERS",&
            error=error),cp_p_file)) THEN
          CALL get_output_format(print_key,my_form=my_form,my_ext=my_ext,error=error)
       END IF
       IF (first_time.OR.(.NOT.qs_loc_env%first_time)) THEN
          iunit=cp_print_key_unit_nr(logger,print_key,"",extension=my_ext,file_form=my_form,&
               middle_name=TRIM(qs_loc_env%tag_mo)//"_centers_s"//TRIM(ADJUSTL(cp_to_string(ispin))),&
               log_filename=.FALSE.,on_file=.TRUE.,is_new_file=init_traj,error=error)
          IF (iunit>0) THEN
             ! Gather units of measure for output (if available)
             CALL section_vals_val_get(print_key,"UNIT",c_val=unit_str,error=error)
             unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

             IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,"/IONS+CENTERS",error=error),cp_p_file)) THEN
                ! Different possible formats
                CALL print_wannier_traj(qs_loc_env, print_key, center, iunit, init_traj, unit_conv, error)
             ELSE
                ! Default print format
                iter=cp_iter_string(logger%iter_info,error=error)
                WRITE(iunit,'(i6,/,a)') nstates,TRIM(iter)
                DO l=1,nstates
                   WRITE(iunit,'(A,4F16.8)') "X", unit_conv*center(1:4,l)
                END DO
             END IF
          END IF
          CALL cp_print_key_finished_output(iunit,logger,print_key,on_file=.TRUE.,error=error)
       END IF
    END IF
  END SUBROUTINE print_wannier_centers

! *****************************************************************************
!> \brief computes spread and centers
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE print_wannier_traj(qs_loc_env, print_key, center, iunit, init_traj, unit_conv, error)
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: print_key
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    INTEGER, INTENT(IN)                      :: iunit
    LOGICAL, INTENT(IN)                      :: init_traj
    REAL(KIND=dp), INTENT(IN)                :: unit_conv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_wannier_traj', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: iter, remark1, remark2, title
    INTEGER                                  :: i, iskip, natom, ntot, &
                                                outformat
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(particle_set, atomic_kind_set, atomic_kind, logger)
    logger => cp_error_get_logger(error)
    natom = SIZE(qs_loc_env%particle_set)
    ntot =  natom + SIZE(center,2)
    CALL allocate_particle_set(particle_set, ntot,error)
    CALL allocate_atomic_kind_set(atomic_kind_set,1,error)
    atomic_kind => atomic_kind_set(1)
    CALL set_atomic_kind(atomic_kind=atomic_kind,kind_number=0,&
         name="X",element_symbol="X",mass=0.0_dp)
    ! Particles
    DO i = 1, natom
       particle_set(i)%atomic_kind => qs_loc_env%particle_set(i)%atomic_kind
       particle_set(i)%r = pbc(qs_loc_env%particle_set(i)%r,qs_loc_env%cell)
    END DO
    ! Wannier Centers
    DO i = natom+1, ntot
       particle_set(i)%atomic_kind => atomic_kind
       particle_set(i)%r = pbc(center(1:3,i-natom),qs_loc_env%cell)
    END DO
    ! Dump the structure
    CALL section_vals_val_get(print_key,"FORMAT",i_val=outformat,error=error)

    ! Header file
    SELECT CASE (outformat)
    CASE (dump_dcd,dump_dcd_aligned_cell)
       IF ( init_traj ) THEN
          !Lets write the header for the coordinate dcd
          ! Note (TL) : even the new DCD format is unfortunately too poor
          !             for our capabilities.. for example here the printing
          !             of the geometry could be nested inside several iteration
          !             levels.. this cannot be exactly reproduce with DCD.
          !             Just as a compromise let's pick-up the value of the MD iteration
          !             level. In any case this is not any sensible information for the standard..
          iskip = section_get_ival(print_key,"EACH%MD",error=error)
          WRITE(iunit) "CORD",0,-1,iskip,&
               0,0,0,0,0,0,REAL(0,KIND=sp),1,0,0,0,0,0,0,0,0,24
          remark1= "REMARK FILETYPE CORD DCD GENERATED BY CP2K"
          remark2= "REMARK Support new DCD format with cell information"
          WRITE(iunit) 2,remark1,remark2
          WRITE(iunit) ntot
          CALL m_flush(iunit)
       END IF
    CASE (dump_xmol)
       iter = cp_iter_string(logger%iter_info,error=error)
       WRITE (UNIT=title,FMT="(A)")" Particles+Wannier centers. Iteration:"//TRIM(iter)
    CASE DEFAULT
       title = ""
    END SELECT
    CALL write_particle_coordinates(particle_set,iunit,outformat,"POS",title,qs_loc_env%cell,&
         unit_conv=unit_conv, error=error)
    CALL m_flush(iunit)
    CALL deallocate_particle_set(particle_set,error)
    CALL deallocate_atomic_kind_set(atomic_kind_set,error)
  END SUBROUTINE print_wannier_traj

END MODULE qs_loc_methods
