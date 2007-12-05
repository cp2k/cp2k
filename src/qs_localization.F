!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_localization [1.0] *
!!
!!   NAME
!!     qs_localization
!!
!!   FUNCTION
!!
!!     should be able to find properties related to the position operator in PBC
!!     e.g. location and spread of periodic orbitals
!!          Wannier functions
!!
!!   AUTHORS
!!     Joost VandeVondele with some help of Leonardo Bernasconi (01.2003)
!!     CJM (04.2003)
!!     CJM & IKUO (08.2003)
!!     Teo (07.2007) - New IO output for wannier centers
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_localization
  ! *****************************************************************************
  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             deallocate_atomic_kind_set,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE input_constants,                 ONLY: do_loc_crazy,&
                                             do_loc_direct,&
                                             do_loc_jacobi,&
                                             do_loc_l1_norm_sd,&
                                             dump_dcd,&
                                             dump_xmol
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             sp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: twopi
  USE md_energies,                     ONLY: get_output_format
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_max,&
                                             mp_minloc,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             deallocate_particle_set,&
                                             particle_type,&
                                             write_particle_coordinates
  USE physcon,                         ONLY: angstrom,&
                                             debye
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_localization_methods,         ONLY: approx_l1_norm_sd,&
                                             crazy_rotations,&
                                             direct_mini,&
                                             initialize_weights,&
                                             jacobi_rotations
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE f77_blas
#include "cp_common_uses.h"

  IMPLICIT NONE

  ! *** Public ***
  PUBLIC :: qs_loc_init
  PUBLIC :: qs_loc_env_type, PIPEK_MEZEY, POSITION_OPERATOR
  PUBLIC :: qs_loc_compute_properties
  PUBLIC :: qs_loc_destroy

  ! *** Localization Data Type ***
  TYPE qs_loc_env_type
     LOGICAL :: wannier
     LOGICAL :: molecular_states
     LOGICAL :: jacobi_fallback
     INTEGER :: max_iter
     INTEGER :: is_initialized
     INTEGER :: operator_type
     INTEGER :: localization_method
     REAL(KIND=dp) :: eps_localization
     REAL (KIND = dp) :: weights ( 3 )
     CHARACTER(LEN=default_path_length) :: filename
     TYPE ( real_matrix_p_type ),POINTER :: OPERATOR(:,:)
     TYPE ( cp_para_env_type ), POINTER :: para_env
     TYPE ( cell_type ), POINTER :: cell
     TYPE ( particle_type ), POINTER :: particle_set ( : )
     TYPE(distribution_1d_type), POINTER         :: local_molecules
  END TYPE qs_loc_env_type

  INTEGER, PARAMETER :: PIPEK_MEZEY=301, POSITION_OPERATOR=302
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_localization'
  ! *****************************************************************************

CONTAINS

!!****f* qs_localization/qs_loc_init [1.0] *
!!
!!   NAME
!!     qs_loc_init
!!
!!   FUNCTION
!!     allocates the data, and initializes the operators
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE qs_loc_init(qs_loc_env,qs_env,error)

    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_loc_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, isos, j, natoms
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

    CALL timeset("qs_loc_init","I",'',handle)
    NULLIFY(qs_loc_env%para_env)
    NULLIFY(qs_loc_env%cell)
    NULLIFY(qs_loc_env%particle_set)
    !
    ! start replicating the s structure, we do have
    ! 6 operators (x,y,z)*(real,imag) for POSITION_OPERATOR
    ! and we have NATOMS operators for PIPEK_MEZEY
    ! notice that the operator should have the same structure as S
    ! borrow para env from S
    !
    CALL get_qs_env(qs_env=qs_env,matrix_s=matrix_s,cell=qs_loc_env%cell, &
         local_molecules=qs_loc_env%local_molecules,&
         particle_set = qs_loc_env % particle_set, &
         para_env=qs_loc_env%para_env,error=error)
    CALL cp_para_env_retain(qs_loc_env%para_env,error=error)

    qs_loc_env % operator_type = POSITION_OPERATOR
    IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  THEN
       ALLOCATE ( qs_loc_env % operator ( 2, 3 ), stat = isos )
       IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
    ELSEIF ( qs_loc_env % operator_type == PIPEK_MEZEY ) THEN
       natoms = SIZE ( qs_loc_env % particle_set )
       ALLOCATE ( qs_loc_env % operator ( 1, natoms ), stat = isos )
       IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'operator' )
    ENDIF

    DO i=1, SIZE ( qs_loc_env % operator, 2 )
       DO j=1, SIZE ( qs_loc_env % operator, 1 )
          NULLIFY(qs_loc_env%operator(j,i)%matrix)
          CALL replicate_matrix_structure(matrix_s(1)%matrix, &
               qs_loc_env%operator(j,i)%matrix,"operator",error=error)
          CALL set_matrix(qs_loc_env%operator(j,i)%matrix,0.0_dp)
       ENDDO
    END DO


    IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN
       CALL initialize_weights ( qs_loc_env % cell, qs_loc_env % weights )
       CALL get_boys_foster_operator ( qs_loc_env, qs_env ,error=error)
    ENDIF

    IF ( qs_loc_env % operator_type == PIPEK_MEZEY ) &
         CALL get_pipek_mezey_operator ( qs_loc_env, qs_env )

    qs_loc_env%molecular_states = .FALSE.
    qs_loc_env%is_initialized = 777
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qs_loc_init

!!****f* qs_localization/get_boys_foster_operator [1.0] *
!!
!!   NAME
!!     get_boys_foster_operator
!!
!!   FUNCTION
!!     Computes the Boys-Foster operator for periodic systems
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE get_boys_foster_operator ( qs_loc_env, qs_env, error )
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_boys_foster_operator', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ix, iy, iz, j
    REAL(KIND=dp)                            :: thearg, theval, x, y, z
    TYPE(pw_p_type)                          :: op_rspace
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    NULLIFY ( pw_env, auxbas_pw_pool )
    !
    ! OK, now we construct these operators the easy way, put them on the pw grid and integrate
    ! analytical construction should be faster.
    !
    ! get a pw grid
    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
    CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,op_rspace%pw,&
         use_data=REALDATA3D,in_space=REALSPACE,error=error)

    DO i=1,3
       DO j=1,2
          ! CALL coeff_zero(op_rspace)
          ! the definition of the position on the PW grid has to be
          ! compatible with the one used in the integration procedure / collocation procedure
          ! integrate_v_rspace calculate_rho_elec
          DO iz=op_rspace%pw%pw_grid%bounds_local(1,3),op_rspace%pw%pw_grid%bounds_local(2,3)
             DO iy=op_rspace%pw%pw_grid%bounds_local(1,2),op_rspace%pw%pw_grid%bounds_local(2,2)
                DO ix=op_rspace%pw%pw_grid%bounds_local(1,1),op_rspace%pw%pw_grid%bounds_local(2,1)
                   x=REAL(ix-op_rspace%pw%pw_grid%bounds(1,1),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(1),dp)
                   y=REAL(iy-op_rspace%pw%pw_grid%bounds(1,2),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(2),dp)
                   z=REAL(iz-op_rspace%pw%pw_grid%bounds(1,3),dp)*twopi/REAL(op_rspace%pw%pw_grid%npts(3),dp)
                   SELECT CASE (i)
                   CASE (1)
                      thearg=x
                   CASE (2)
                      thearg=y
                   CASE (3)
                      thearg=z
                   END SELECT
                   SELECT CASE (j)
                   CASE (1)
                      theval=COS(thearg)
                   CASE (2)
                      theval=SIN(thearg)
                   END SELECT
                   ! notice the volume element, slightly ugly I think
                   op_rspace%pw%cr3d(ix,iy,iz)=theval*op_rspace%pw%pw_grid%dvol
                ENDDO
             ENDDO
          ENDDO
          ! do the actual integration
          CALL set_matrix(qs_loc_env%operator(j,i)%matrix,0.0_dp)
          CALL integrate_v_rspace(op_rspace,h=qs_loc_env%operator(j,i),qs_env=qs_env,calculate_forces=.FALSE.,error=error)
       ENDDO
    ENDDO

    CALL pw_pool_give_back_pw ( auxbas_pw_pool, op_rspace%pw ,error=error)

  END SUBROUTINE get_boys_foster_operator

!!****f* qs_localization/get_pipek_mezey_operator [1.0] *
!!
!!   NAME
!!     get_pipek_mezey_operator
!!
!!   FUNCTION
!!     Computes the Pipek-Mezey operator for periodic systems
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE get_pipek_mezey_operator ( qs_loc_env, qs_env )
    ! globals
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'get_pipek_mezey_operator', &
      routineP = moduleN//':'//routineN

    CALL stop_program("qs_localization","PIPEK-MEZEY not implemented")
  END SUBROUTINE get_pipek_mezey_operator

!!****f* qs_localization/qs_loc_compute_properties [1.0] *
!!
!!   NAME
!!     qs_loc_compute_properties
!!
!!   FUNCTION
!!     applies the operators to obtain the necessary data
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE qs_loc_compute_properties(qs_loc_env, molecule_set, vectors, center, ionode,&
       loc_print_key, error)

    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(cp_fm_type), POINTER                :: vectors
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: center
    LOGICAL                                  :: ionode
    TYPE(section_vals_type), POINTER         :: loc_print_key
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_loc_compute_properties', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, isos, j, k, n, &
                                                natoms, scr, sweeps
    LOGICAL                                  :: converged
    TYPE(cp_fm_p_type), ALLOCATABLE          :: ZIJ( :, : )
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: opvec
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I",'',handle)
    NULLIFY(tmp_fm_struct,opvec)
    logger => cp_error_get_logger(error)

    ! get rows and cols of the input
    CALL cp_fm_get_info(vectors,nrow_global=n,ncol_global=k,error=error)
    ! replicate the input kind of matrix
    CALL cp_fm_create(opvec,vectors%matrix_struct,error=error)
    ! create a new matrix structure k x k in size, borrowing most of info
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k, &
         ncol_global=k,para_env=vectors%matrix_struct%para_env, &
         context=vectors%matrix_struct%context,error=error)

    IF (ANY(SHAPE(center).NE.(/k,6/))) CALL stop_program(routineP,'center')

    ! Allocating ZIJ
    IF ( qs_loc_env % operator_type == POSITION_OPERATOR ) THEN
       ALLOCATE ( ZIJ ( 2, 3 ), stat = isos )
       IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
    ELSEIF ( qs_loc_env % operator_type == PIPEK_MEZEY ) THEN
       natoms = SIZE ( qs_loc_env % operator, 2 )
       ALLOCATE ( ZIJ ( 1, natoms ), stat = isos )
       IF ( isos /= 0 ) CALL stop_memory ( 'qs_localization', 'ZIJ' )
    ENDIF
    ! make the actual matrices
    DO i=1,SIZE ( ZIJ, 2 )
       DO j=1, SIZE ( ZIJ, 1 )
          NULLIFY ( ZIJ ( j, i ) % matrix )
          CALL cp_fm_create ( ZIJ ( j, i ) % matrix, tmp_fm_struct ,error=error)
       ENDDO
    ENDDO

    ! release one instance of the fm_struct
    CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

    ! Compute zij here
    DO i=1,SIZE ( ZIJ, 2 )
       DO j=1,SIZE ( ZIJ, 1 )
          CALL cp_sm_fm_multiply(qs_loc_env%operator(j,i)%matrix, &
               vectors,opvec,ncol=k ,error=error)
          CALL cp_fm_gemm("T","N",k,k,n,1.0_dp,vectors,opvec,0.0_dp,ZIJ(j,i)%matrix,error=error)
       ENDDO
    ENDDO

    IF ( qs_loc_env % wannier ) THEN
       SELECT CASE (qs_loc_env % localization_method)
       CASE ( do_loc_jacobi )
          CALL jacobi_rotations ( qs_loc_env % weights, zij, vectors, &
               max_iter=qs_loc_env%max_iter,eps_localization=qs_loc_env%eps_localization, &
               sweeps=sweeps )
       CASE ( do_loc_crazy )
          CALL crazy_rotations  ( qs_loc_env % weights, zij, vectors, &
               max_iter=qs_loc_env%max_iter,eps_localization=qs_loc_env%eps_localization, &
               iterations=sweeps, converged=converged, error=error )
          IF ( .NOT.converged ) THEN
             IF (qs_loc_env % jacobi_fallback) THEN
                WRITE (cp_logger_get_default_unit_nr(logger,local=.FALSE.),"(A,I6,A)")&
                     " Crazy Wannier localization not converged after ",sweeps,&
                     " iterations, switching to jacobi rotations"
                CALL jacobi_rotations ( qs_loc_env % weights, zij, vectors,&
                     max_iter=qs_loc_env%max_iter,eps_localization=qs_loc_env%eps_localization, sweeps = sweeps )
             ELSE
                WRITE (cp_logger_get_default_unit_nr(logger,local=.FALSE.),"(A,I6,A)")&
                     " Crazy Wannier localization not converged after ",sweeps,&
                     " iterations, and jacobi_fallback switched off "
             ENDIF
          END IF
       CASE ( do_loc_direct )
          CALL direct_mini      ( qs_loc_env % weights, zij, vectors, sweeps )
       CASE( do_loc_l1_norm_sd )
          CALL approx_l1_norm_sd(vectors,qs_loc_env%max_iter,qs_loc_env%eps_localization,&
               &                 converged,sweeps,error=error)
          ! here we need to set zij for the computation of the centers and spreads
          DO i=1,SIZE ( ZIJ, 2 )
             DO j=1,SIZE ( ZIJ, 1 )
                CALL cp_sm_fm_multiply(qs_loc_env%operator(j,i)%matrix, &
                     vectors,opvec,ncol=k ,error=error)
                CALL cp_fm_gemm("T","N",k,k,n,1.0_dp,vectors,opvec,&
                     &          0.0_dp,ZIJ(j,i)%matrix,error=error)
             ENDDO
          ENDDO
       CASE DEFAULT
          CALL stop_program("qs_localization","unknown localization method")
       END SELECT
       IF (ionode) THEN
          scr=cp_logger_get_default_unit_nr(logger,.FALSE.)
          WRITE(scr,'(A,I6,A)') " Wannier localization finished in ",sweeps," iterations"
       ENDIF
    ELSE
       IF (ionode) THEN
          scr=cp_logger_get_default_unit_nr(logger,.FALSE.)
          WRITE(scr,'(A,I6,A)') " No Wannier localization applied "
       ENDIF
    ENDIF

    IF ( qs_loc_env % operator_type == POSITION_OPERATOR )  THEN
       CALL centers_and_spreads ( ionode, qs_loc_env, loc_print_key, zij, center, error )
       IF (qs_loc_env % molecular_states ) &
            CALL wfc_to_molecule ( qs_loc_env, loc_print_key,center, molecule_set,error )
    ENDIF

    ! release results
    DO i=1,SIZE ( ZIJ, 2 )
       DO j=1, SIZE ( ZIJ, 1 )
          CALL cp_fm_release(ZIJ(j,i)%matrix,error=error)
       ENDDO
    ENDDO

    DEALLOCATE(zij, STAT = isos)
    IF(isos/=0) CALL stop_memory('qs_loc_compute_properties','zij')
    CALL cp_fm_release(opvec,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_loc_compute_properties

!!****f* qs_localization/centers_and_spreads [1.0] *
!!
!!   NAME
!!     centers_and_spreads
!!
!!   FUNCTION
!!     computes spread and centers
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE centers_and_spreads ( ionode, qs_loc_env, loc_print_key, zij, center, error )
    LOGICAL, INTENT(IN)                      :: ionode
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_print_key
    TYPE(cp_fm_p_type), INTENT(INOUT)        :: zij( :, : )
    REAL(KIND=dp), INTENT(OUT)               :: center( :, : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'centers_and_spreads', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: almost_determinant, z
    INTEGER                                  :: i, l, nstates, scr
    REAL(KIND=dp)                            :: imagpart, realpart
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_type), POINTER               :: complexresult
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ! get nstates
    CALL cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstates ,error=error)

    ! initialize complexresult for boys-foster
    NULLIFY ( complexresult )
    CALL cp_cfm_create ( complexresult, zij ( 1, 1 ) % matrix % matrix_struct ,error=error)

    ! OK print out the centers ( center ( 1:3) )  and spreads ( center ( 4:6 ) ) of these friends
    cell => qs_loc_env % cell
    DO l=1,nstates
       DO i=1,3
          CALL cp_fm_get_element(zij(1,i)%matrix,l,l,realpart)
          CALL cp_fm_get_element(zij(2,i)%matrix,l,l,imagpart)
          z=CMPLX(realpart,imagpart,dp)
          center (l ,i )=( cell % hmat ( i, i ) / twopi ) * AIMAG ( LOG ( z ) )
          center (l, i+3 )=SQRT( -qs_loc_env % weights ( i ) * 2.0_dp * LOG ( ABS ( Z ) ) ) / twopi
       ENDDO
    ENDDO
    
    ! Print wannier centers
    CALL print_wannier_centers(qs_loc_env, loc_print_key, center, logger, error)

    ! now, find the spread of the orbitals
    DO i=1,3
       complexresult%local_data ( :, : ) = CMPLX (ZIJ(1,i)%matrix%local_data(:,:), &
            ZIJ(2,i)%matrix%local_data(:,:),dp)
       CALL cp_cfm_lu_decompose ( complexresult, almost_determinant )
       IF (ionode) THEN
          scr=cp_logger_get_default_unit_nr()
          WRITE(scr,'(A,I6,F16.10)') " D Eq. 11 ",i,SQRT(-2.0*nstates*LOG(ABS(almost_determinant)))
       ENDIF
    ENDDO
    CALL cp_cfm_release(complexresult,error=error)

  END SUBROUTINE centers_and_spreads

!!****f* qs_localization/print_wannier_centers [1.0] *
!!
!!   NAME
!!     print_wannier_centers
!!
!!   FUNCTION
!!     Prints wannier centers
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE print_wannier_centers(qs_loc_env, loc_print_key, center, logger, error)
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_print_key
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_wannier_centers', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: iter
    CHARACTER(LEN=default_string_length)     :: my_ext, my_form
    INTEGER                                  :: iunit, l, nstates
    LOGICAL                                  :: init_traj

    nstates = SIZE(center,1)
    my_form = "formatted"
    my_ext  = ".data"
    IF (BTEST(cp_print_key_should_output(logger%iter_info,loc_print_key,"/IONS+CENTERS",error=error),cp_p_file)) THEN
       CALL get_output_format(loc_print_key,my_form=my_form,my_ext=my_ext,error=error)
    END IF
    iunit=cp_print_key_unit_nr(logger,loc_print_key,"",extension=my_ext,file_form=my_form,&
         middle_name=TRIM(qs_loc_env%filename),log_filename=.FALSE.,on_file=.TRUE.,&
         is_new_file=init_traj,error=error)
    IF (iunit>0) THEN
       IF (BTEST(cp_print_key_should_output(logger%iter_info,loc_print_key,"/IONS+CENTERS",error=error),cp_p_file)) THEN
          ! Different possible formats
          CALL print_wannier_traj(qs_loc_env, loc_print_key, center, iunit, init_traj, error)
       ELSE
          ! Default print format
          iter=cp_iter_string(logger%iter_info,error=error)
          WRITE(iunit,'(i6,/,a)') nstates,TRIM(iter)
          DO l=1,nstates
             WRITE(iunit,'(A,6F16.8)') "X", angstrom*center(l, 1 : 6 )
          END DO
       END IF
    END IF
    CALL cp_print_key_finished_output(iunit,logger,loc_print_key,on_file=.TRUE.,error=error)

  END SUBROUTINE print_wannier_centers

!!****f* qs_localization/print_wannier_traj [1.0] *
!!
!!   NAME
!!     print_wannier_traj
!!
!!   FUNCTION
!!     computes spread and centers
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE print_wannier_traj(qs_loc_env, loc_print_key, center, iunit, init_traj, error)
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_print_key
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    INTEGER, INTENT(IN)                      :: iunit
    LOGICAL, INTENT(IN)                      :: init_traj
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
    ntot =  natom + SIZE(center,1)
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
       particle_set(i)%r = pbc(center(i-natom,1:3),qs_loc_env%cell)
    END DO
    ! Dump the structure
    CALL section_vals_val_get(loc_print_key,"FORMAT",i_val=outformat,error=error)
    IF(outformat==dump_dcd) THEN
       IF ( init_traj ) THEN
          !Lets write the header for the coordinate dcd
          iskip = section_get_ival(loc_print_key,"EACH",error=error)
          WRITE(iunit) "CORD",0,-1,iskip,&
               0,0,0,0,0,0,REAL(0,KIND=sp),1,0,0,0,0,0,0,0,0,24
          remark1= "REMARK FILETYPE CORD DCD GENERATED BY CP2K"
          remark2= "REMARK Support new DCD format with cell information"
          WRITE(iunit) 2,remark1,remark2
          WRITE(iunit) ntot
          CALL m_flush(iunit)
       ENDIF
    ELSEIF(outformat==dump_xmol) THEN
       iter = cp_iter_string(logger%iter_info,error=error) 
       WRITE (UNIT=title,FMT="(A)")" Particles+Wannier centers. Iteration:"//TRIM(iter)
    ELSE
       title = ""
    END IF
    CALL write_particle_coordinates(particle_set,iunit,outformat,"POS",title,qs_loc_env%cell)
    CALL m_flush(iunit)
    CALL deallocate_particle_set(particle_set,error)
    CALL deallocate_atomic_kind_set(atomic_kind_set,error)
  END SUBROUTINE print_wannier_traj

!!****f* qs_localization/qs_loc_destroy [1.0] *
!!
!!   NAME
!!     qs_loc_destroy
!!
!!   FUNCTION
!!     Deallocates all data
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE qs_loc_destroy ( qs_loc_env, error )

    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_loc_destroy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

    IF (qs_loc_env%is_initialized.ne.777) CALL stop_program(routineP,"problems")

    DO i=1,SIZE ( qs_loc_env % operator, 2 )
       DO j=1,SIZE ( qs_loc_env % operator, 1 )
          CALL deallocate_matrix(qs_loc_env%operator(j,i)%matrix,error=error)
       ENDDO
    END DO
    DEALLOCATE(qs_loc_env%operator)

    CALL cp_para_env_release(qs_loc_env%para_env,error=error)

  END SUBROUTINE qs_loc_destroy

!!****f* qs_localization/wfc_to_molecule [1.0] *
!!
!!   NAME
!!     wfc_to_molecule
!!
!!   FUNCTION
!!     maps wfc's to molecules
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - qs_env: the qs_env in which the qs_env lives
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!*** **********************************************************************
  SUBROUTINE wfc_to_molecule ( qs_loc_env, loc_print_key, center, molecule_set, error )
    TYPE(qs_loc_env_type), INTENT(IN)        :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_print_key
    REAL(KIND=dp), INTENT(IN)                :: center( :, : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'wfc_to_molecule', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: mol_name
    INTEGER :: counter, first_atom, i, iatom, ikind, imol, imol_now, isos, &
      istate, local_location, natom, natom_loc, natom_max, nkind, nmol, &
      nstate, output_unit
    INTEGER, POINTER                         :: wfc_to_atom_map( : )
    REAL(KIND=dp)                            :: dipole(3), dr(3), mydist(2), &
                                                r_pbc( 3 ), zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: charge_set
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dipole_set
    REAL(KIND=dp), POINTER                   :: distance( : ), r( :, : )
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(particle_type), POINTER             :: particle_set( : )

    logger => cp_error_get_logger(error)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    particle_set => qs_loc_env % particle_set
    para_env => qs_loc_env % para_env
    local_molecules => qs_loc_env % local_molecules
    nstate = SIZE ( center, 1 )
    ALLOCATE(wfc_to_atom_map(nstate),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    nkind = SIZE (local_molecules%n_el)
    natom = 0
    natom_max = 0
    DO ikind=1,nkind
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol
          i = local_molecules%list(ikind)%array(imol)
          molecule_kind => molecule_set(i)%molecule_kind
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               natom=natom)
          natom_max = natom_max + natom
          molecule_set(i)%lmi%nstates=0
          IF(ASSOCIATED(molecule_set(i)%lmi%states)) THEN
             DEALLOCATE(molecule_set(i)%lmi%states,STAT=isos)
             IF(isos/=0) CALL stop_memory('wfc_to_molecule','states')
          END IF
       END DO
    END DO
    natom_loc = natom_max
    natom = natom_max

    CALL mp_max(natom_max,para_env%group)

    ALLOCATE(r(3,natom_max),STAT=isos)
    IF(isos/=0) CALL stop_memory('centers_and_spreads','r')

    ALLOCATE(distance(natom_max),STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')

    !Zero all the stuff
    r(:,:) = 0.0_dp
    distance(:) = 1.E10_dp

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    counter = 0
    nkind = SIZE (local_molecules%n_el)
    DO ikind=1,nkind
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol
          i = local_molecules%list(ikind)%array(imol)
          molecule_kind => molecule_set(i)%molecule_kind
          first_atom = molecule_set(i)%first_atom
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               natom=natom,&
               name=mol_name,&
               atom_list=atom_list)
          DO iatom=1,natom
             counter = counter + 1
             r ( :, counter ) = particle_set(first_atom+iatom-1) % r ( : )
          END DO
       END DO
    END DO

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    DO istate=1,nstate
       distance(:) = 1.E10_dp
       DO iatom=1,natom_loc
          dr(1) = r(1,iatom) - center(istate,1)
          dr(2) = r(2,iatom) - center(istate,2)
          dr(3) = r(3,iatom) - center(istate,3)
          r_pbc = pbc ( dr, qs_loc_env%cell )
          distance(iatom) = SQRT ( DOT_PRODUCT ( r_pbc, r_pbc ) )
       END DO

       !combine distance() from all procs
       local_location = MINLOC(distance,DIM=1)

       mydist(1) = distance(local_location)
       mydist(2) = para_env%mepos

       CALL mp_minloc(mydist,para_env%group)

       IF(mydist(2)==para_env%mepos) THEN
          wfc_to_atom_map(istate)  =  local_location
       ELSE
          wfc_to_atom_map(istate)  =  0
       END IF
    END DO
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    IF(natom_loc/=0) THEN
       DO istate = 1,nstate
          iatom = wfc_to_atom_map(istate)
          IF(iatom/=0) THEN
             counter = 0
             nkind = SIZE (local_molecules%n_el)
             DO ikind=1,nkind
                nmol = SIZE ( local_molecules%list(ikind)%array )
                DO imol=1,nmol
                   imol_now = local_molecules%list(ikind)%array(imol)
                   molecule_kind => molecule_set(imol_now)%molecule_kind
                   CALL get_molecule_kind(molecule_kind=molecule_kind,&
                        natom=natom)
                   counter = counter + natom
                   IF(counter>=iatom) EXIT
                END DO
                IF(counter>=iatom) EXIT
             END DO
             i = molecule_set(imol_now)%lmi%nstates
             i = i + 1
             molecule_set(imol_now)%lmi%nstates = i
             CALL reallocate(molecule_set(imol_now)%lmi%states,1,i)
             molecule_set(imol_now)%lmi%states(i) = istate
          END IF
       END DO
    END IF

    !---------------------------------------------------------------------------
    ! figure out dipole of the molecule.
    ! Assume it is neutral (so origin doesn't matter, we take the first atom)
    ! if the molecule is larger than half the box there might be a problem folding stuff
    ! FIXME this should be moved out of this place
    !---------------------------------------------------------------------------
    particle_set => qs_loc_env % particle_set
    para_env => qs_loc_env % para_env
    local_molecules => qs_loc_env % local_molecules
    nstate = SIZE ( center, 1 )
    ALLOCATE(dipole_set(3,SIZE(molecule_set)))
    ALLOCATE(charge_set(SIZE(molecule_set)))
    dipole_set=0.0_dp
    charge_set=0.0_dp
    nkind = SIZE (local_molecules%n_el)
    DO ikind=1,nkind ! loop over different molecules
       nmol = SIZE ( local_molecules%list(ikind)%array )
       DO imol=1,nmol ! all the molecules of the kind
          imol_now = local_molecules%list(ikind)%array(imol) ! index in the global array
          IF (.NOT.ASSOCIATED(molecule_set(imol_now)%lmi%states)) CYCLE
          molecule_kind => molecule_set(imol_now)%molecule_kind
          first_atom = molecule_set(imol_now)%first_atom
          CALL get_molecule_kind(molecule_kind=molecule_kind,&
               natom=natom,atom_list=atom_list)
          dipole=0.0_dp
          ! charges of the atoms involved
          DO iatom=1,natom
             i=first_atom+iatom-1
             atomic_kind=>particle_set(i)%atomic_kind
             CALL get_atomic_kind(atomic_kind,gth_potential=gth_potential)
             CALL get_potential(gth_potential,zeff=zeff)
             r_pbc = particle_set(i) % r - particle_set(first_atom) % r
             r_pbc = pbc ( r_pbc , qs_loc_env%cell )
             dipole=dipole + zeff*r_pbc
             charge_set(imol_now)=charge_set(imol_now)+zeff
          END DO
          ! charges of the wfc involved
          DO istate=1,SIZE(molecule_set(imol_now)%lmi%states)
             i = molecule_set(imol_now)%lmi%states(istate)
             r_pbc = center(i,1:3) - particle_set(first_atom) % r
             r_pbc = pbc ( r_pbc , qs_loc_env%cell )
             ! XXXXXXXXXXXXXXXXXXXX double occupation FIXME XXXXXXXXX (where are the spins?)
             dipole=dipole - 2.0_dp * r_pbc 
             charge_set(imol_now)=charge_set(imol_now)-2.0_dp
          ENDDO
          dipole_set(:,imol_now)=dipole ! a.u.
       ENDDO
    ENDDO
    CALL mp_sum(dipole_set,para_env%group)
    CALL mp_sum(charge_set,para_env%group)
    
    output_unit=cp_print_key_unit_nr(logger,loc_print_key,"/MOLECULAR_DIPOLES",&
         extension=".data",middle_name="MOLECULAR_DIPOLES",error=error)
    IF (output_unit>0) THEN
       WRITE(UNIT=output_unit,FMT='(A80)') "# molecule nr, charge, dipole vector, dipole (Debye) "
       dipole_set=dipole_set*debye ! Debye
       DO I=1,SIZE(dipole_set,2)
          WRITE(UNIT=output_unit,FMT='(I6,5F12.6)') I,charge_set(I),dipole_set(1:3,I), &
               SQRT(DOT_PRODUCT(dipole_set(1:3,I),dipole_set(1:3,I)))
       ENDDO
    ENDIF
    CALL cp_print_key_finished_output(output_unit,logger,loc_print_key,&
         "/MOLECULAR_DIPOLES",error=error)
    DEALLOCATE(dipole_set,charge_set)
    !---------------------------------------------------------------------------
    ! end of molecular dipole calculation
    !---------------------------------------------------------------------------

    DEALLOCATE(distance,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','distance')
    DEALLOCATE(r,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','r')

    DEALLOCATE(wfc_to_atom_map,STAT=isos)
    IF(isos/=0) CALL stop_memory('wfc_to_molecule','wfc_to_atom_map')

  END SUBROUTINE wfc_to_molecule
  !------------------------------------------------------------------------------

END MODULE qs_localization

