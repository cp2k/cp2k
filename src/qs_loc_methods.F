!   CP2K: A general program to perform moleocular dynamics simulations         !
!   Copyright (C) 2004 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_loc_methods_new [1.0] *
!!
!!   NAME
!!     qs_loc_methods_new
!!
!!   FUNCTION
!!     Driver for the localization that should be general
!!     for all the methods available and all the definition of the
!!     spread functional
!!     Write centers, spread and cubes only if required and for the
!!     selected states
!!     The localized functions are copied in the standard mos array
!!     for the next use
!!
!!   AUTHORS
!!     MI (04.2005)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
MODULE qs_loc_methods
! *****************************************************************************

  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             localized_wfn_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_files,                        ONLY: open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace,&
                                             cp_fm_transpose
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                             cp_fm_pool_p_type,&
                                             fm_pool_create_fm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_maxabsval,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr
  USE input_constants,                 ONLY: do_loc_min,&
                                             do_loc_max,&
                                             do_loc_jacobi,&
                                             do_loc_crazy,&
                                             do_loc_direct,&
                                             do_loc_none,&
                                             op_loc_berry,&
                                             op_loc_boys,&
                                             op_loc_pipek,&
                                             state_loc_all,&
                                             state_loc_upton,&
                                             state_loc_fromitoj,&
                                             state_loc_list,&
                                             state_loc_none
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: pi, twopi
  USE message_passing,                 ONLY: mp_max
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                             qs_loc_env_new_type
  USE qs_localization_methods,         ONLY: crazy_rotations,&
                                             direct_mini,&
                                             jacobi_rotations
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type

  USE qs_loc_types,                    ONLY: get_qs_loc_env,& 
                                             qs_loc_env_new_type, &
                                             set_qs_loc_env
  USE physcon,                         ONLY: angstrom
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: module_name = "qs_loc_methods"


! *** Public ***
  PUBLIC ::  qs_loc_driver
!!  PUBLIC :: centers_spreads_berry, optimize_loc_berry, qs_loc_driver

! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****f* qs_loc_methods/centers_spreads_berry

  SUBROUTINE centers_spreads_berry( localized_wfn_control, zij, nmoloc, &
                                    cell, weights, ispin, ionode, error)

     TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
     TYPE(cp_fm_p_type), INTENT(INOUT)        :: zij( :, : )
     TYPE(cell_type), POINTER                 :: cell
     REAL(dp), DIMENSION(:)                   :: weights
     INTEGER, INTENT(IN)                      :: ispin, nmoloc
     LOGICAL, INTENT(IN)                      :: ionode
     TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'centers_spreads_berry', &
      routineP = module_name//':'//routine_name


     CHARACTER(len = 50)                      :: filename1, filename2
     COMPLEX(KIND=dp)                         :: almost_determinant, z
     INTEGER :: idir, istate, iter, nstates, out_c, out_s, scr
     LOGICAL :: failure, file_exists
     REAL(KIND=dp)                            :: imagpart, realpart
     REAL(dp), DIMENSION(3)                   :: spreadtot
     REAL(dp), DIMENSION(:,:), POINTER        :: centers
     TYPE(cp_cfm_type), POINTER               :: complexresult

     failure = .FALSE.
     out_c = 21
     out_s = 22
!  It shoud be the step at which the property localization is calculated
!    iter = localized_wfn_control%iteration_step
     iter = 1
     file_exists = .FALSE.
 

     CALL  cp_fm_get_info ( zij ( 1, 1 ) % matrix, nrow_global = nstates )
     CPPrecondition(nstates>=nmoloc,cp_failure_level,routineP,error,failure)

     centers => localized_wfn_control%centers_set(ispin)%array
     CPPrecondition(ASSOCIATED(centers),cp_failure_level,routineP,error,failure)
     CPPrecondition(SIZE(centers,1)==nmoloc,cp_failure_level,routineP,error,failure)

     IF(ionode) THEN
       ! check if file exists
       filename1 = localized_wfn_control%centers_file_name(ispin)
       filename2 = localized_wfn_control%spreads_file_name(ispin)
       INQUIRE (FILE=filename1,EXIST=file_exists)


       ! open file
       IF(file_exists) THEN
         CALL open_file(file_name=TRIM(filename1),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="APPEND",&
                     file_status="OLD",&
                     unit_number=out_c)  
         WRITE(out_c,'(I8)') nmoloc
         WRITE(out_c,'(A,I10)') " Iteration Step ", iter
         
         CALL open_file(file_name=TRIM(filename2),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="APPEND",&
                     file_status="OLD",&
                     unit_number=out_s)
         WRITE(out_s,'(A,I10)') " Iteration Step ", iter
       ELSE
         CALL open_file(file_name=TRIM(filename1),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="REWIND",&
                     file_status="REPLACE",&
                     unit_number=out_c)
         WRITE(out_c,'(I8)') nmoloc
         WRITE(out_c,'(A,I10)') " Iteration Step ", iter
         CALL open_file(file_name=TRIM(filename2),&
                     file_action="READWRITE",&
                     file_form="FORMATTED",&
                     file_position="REWIND",&
                     file_status="REPLACE",&
                     unit_number=out_s)
         WRITE(out_s,'(A,I10)') " Iteration Step ", iter

       END IF 
     END IF

     DO istate=1,nmoloc
       DO idir=1,3
         CALL cp_fm_get_element(zij(1,idir)%matrix,istate,istate,realpart)
         CALL cp_fm_get_element(zij(2,idir)%matrix,istate,istate,imagpart)
         z=CMPLX(realpart,imagpart,dp)
         centers (istate ,idir )=( cell % hmat ( idir, idir ) / twopi ) * AIMAG ( LOG ( z ) )
         centers (istate, idir+3 )=SQRT( -weights ( idir ) * 2.0_dp * LOG ( ABS ( Z ) ) ) / twopi
      ENDDO
      IF (ionode) THEN
          WRITE(out_c,'(A,3F16.8)') "He ", angstrom*centers(istate, 1 : 3 )
          WRITE(out_s,'(I6,3F16.8)') istate, angstrom*centers(istate, 4 : 6 )
      END IF
     ENDDO

! now, find the spread of the orbitals
! one should consider whether this is really useful
     CALL cp_cfm_create ( complexresult, zij ( 1, 1 ) % matrix % matrix_struct )

     DO idir=1,3
        complexresult%local_data ( :, : ) = CMPLX (zij(1,idir)%matrix%local_data(:,:), &
                                         zij(2,idir)%matrix%local_data(:,:),dp)
        CALL cp_cfm_lu_decompose ( complexresult, almost_determinant )
        spreadtot(idir) =  SQRT(-2.0*nmoloc*LOG(ABS(almost_determinant)))
     ENDDO

    CALL cp_cfm_release(complexresult)

     IF (ionode) THEN
        scr=cp_logger_get_default_unit_nr()
          
       WRITE(scr,'(A,3F16.10)') " Total Spread of the selected mos (Berry) x,y,z ",&
             spreadtot(1:3)
       WRITE(out_s,'(A,3F16.10)')  " Total ",  spreadtot(1:3) 

       close(out_c)
       close(out_s)
     ENDIF

  END SUBROUTINE centers_spreads_berry

!******************************************************************************
!!****f* qs_loc_methods/optimize_loc_berry [1.0] *
!!
!!   NAME
!!     optimize_loc_berry
!!
!!   FUNCTION
!!      Calculate and optimize the spread functional as calculated from
!!      the selected mos  and by the definition using the berry phase
!!      as given by silvestrelli et al
!!      If required the centers and the spreads for each mos selected
!!      are calculated from z_ij and printed in the proper files.
!!      The files should be appended if already exist, in this way 
!!      it would be possible to follow the evlution of the localization
!!      e.g. along a MD trajectory
!!
!!   NOTES
!!      This definition need the use of complex numbers, therefore the 
!!      optimization routines are specific for this case
!!      The file for the centers and the spreads have a xyz format
!!
!!   ARGUMENTS
!!     - method: indicates localization algorithm
!!     - qs_loc_env: new environment for the localization calculations
!!     - vectors: selected mos to be localized
!!     - op_sm_set: sparse matrices containing the integrals of the kind <mi e{iGr} nu>
!!     - zij_fm_set: set of full matrix of size nmoloc x nmoloc, will contain the z_ij numbers 
!!                   as defined by Silvestrelli et al
!!     - ionode: logical that indicates the node responsable for the i/o
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE optimize_loc_berry( method, localized_wfn_control, vectors, op_sm_set, &
                                 zij_fm_set, cell, weights,  ispin, ionode, error)

    INTEGER, INTENT(IN)                         :: method 
    TYPE(localized_wfn_control_type), POINTER   :: localized_wfn_control
    TYPE(cp_fm_type), POINTER                   :: vectors
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      POINTER                                   :: op_sm_set
    TYPE(cp_fm_p_type), DIMENSION(:,:),POINTER  :: zij_fm_set
    TYPE(cell_type), POINTER                    :: cell
    REAL(dp), DIMENSION(:)                      :: weights
    INTEGER , INTENT(IN)                        :: ispin
    LOGICAL, INTENT(IN)                         :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'optimize_loc_berry', &
      routineP = module_name//':'//routine_name


    INTEGER                                   :: handle, i, j, k, l, nao, &
                                                 nmoloc, sweeps, scr
    TYPE(cp_fm_type), POINTER                 :: opvec

    CALL timeset(routineP,"I",'',handle)

    NULLIFY(opvec)

    ! get rows and cols of the input
    CALL cp_fm_get_info(vectors,nrow_global=nao,ncol_global=nmoloc)
    ! replicate the input kind of matrix
    CALL cp_fm_create(opvec,vectors%matrix_struct)


! Compute zij here
    DO i=1,SIZE ( zij_fm_set, 2 )
      DO j=1,SIZE ( zij_fm_set, 1 )
         CALL cp_fm_set_all(zij_fm_set(j,i)%matrix, 0.0_dp, error=error)
         CALL cp_sm_fm_multiply(op_sm_set(j,i)%matrix, vectors,opvec,ncol=nmoloc )
         CALL cp_fm_gemm("T","N",nmoloc,nmoloc,nao,1.0_dp,vectors,opvec,0.0_dp,&
              zij_fm_set(j,i)%matrix)
      ENDDO
    ENDDO

!      CALL centers_spreads_berry( localized_wfn_control, zij_fm_set, nmoloc, &
!           cell, weights, ispin, ionode, error=error ) 

    SELECT CASE ( method ) 
    CASE ( do_loc_jacobi )
       CALL jacobi_rotations ( weights, zij_fm_set, vectors, sweeps ) 
    CASE ( do_loc_crazy )
       CALL crazy_rotations  ( weights, zij_fm_set, vectors, sweeps ) 
    CASE ( do_loc_direct )
       CALL direct_mini      ( weights, zij_fm_set, vectors, sweeps ) 
    CASE ( do_loc_none ) 
      IF (ionode) THEN
          scr=cp_logger_get_default_unit_nr()
          WRITE(scr,'(A,I6,A)') " No MOS localization applied "
      ENDIF
    CASE DEFAULT
       CALL stop_program("qs_localization","unknown localization method")
    END SELECT
    IF (ionode) THEN
       scr=cp_logger_get_default_unit_nr()
       WRITE(scr,'(A,I3,A,I6,A)') " Localization  for spin ",ispin ,&
            " converged in ",sweeps," iterations"
    ENDIF

    IF ( localized_wfn_control%print_centers)  THEN
      CALL centers_spreads_berry( localized_wfn_control, zij_fm_set, nmoloc, &
           cell, weights, ispin, ionode, error=error ) 
    ENDIF 

    CALL cp_fm_release(opvec)
 
    CALL timestop(0.0_dp,handle) 

  END SUBROUTINE optimize_loc_berry

!******************************************************************************

  SUBROUTINE  optimize_loc_pipek(qs_env, method, localized_wfn_control,&
              vectors, zij_fm_set,  ispin, ionode, error)

    TYPE(qs_environment_type), POINTER           :: qs_env
    INTEGER, INTENT(IN)                          :: method 
    TYPE(localized_wfn_control_type), POINTER    :: localized_wfn_control
    TYPE(cp_fm_type), POINTER                    :: vectors
    TYPE(cp_fm_p_type), DIMENSION(:,:),POINTER   :: zij_fm_set
    INTEGER, INTENT(IN)                          :: ispin
    LOGICAL, INTENT(IN)                          :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'optimize_loc_pipek', &
      routineP = module_name//':'//routine_name


    INTEGER                                      :: handle, iatom, isgf, istat, ldz, nao, natom, ncol, nmoloc, scr, sweeps
    INTEGER, DIMENSION(:), ALLOCATABLE           :: first_sgf, last_sgf, nsgf
    LOGICAL :: failure
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
       POINTER                                   :: ao_ao_fm_pools
    TYPE(cp_fm_type), POINTER                    :: ov_fm
    TYPE(cp_fm_type), POINTER                    :: opvec
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                    :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                    :: matrix_s



    CALL timeset(routineP,"I",'',handle)

    failure = .FALSE. 
    NULLIFY(particle_set)
    ! get rows and cols of the input
    CALL cp_fm_get_info(vectors,nrow_global=nao,ncol_global=nmoloc)
    ! replicate the input kind of matrix
    CALL cp_fm_create(opvec,vectors%matrix_struct)
    CALL cp_fm_set_all(opvec, 0.0_dp, error=error)

    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, particle_set=particle_set) 
    natom = SIZE(particle_set,1)
    ALLOCATE (first_sgf(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (nsgf(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

!   construction of 
    CALL get_particle_set(particle_set=particle_set,&
         first_sgf=first_sgf,last_sgf=last_sgf,nsgf=nsgf)

!   Copy the overlap sparse matrix in a full matrix
    CALL mpools_get(qs_env%mpools,ao_ao_fm_pools=ao_ao_fm_pools)
    CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, ov_fm, name=" ", error=error)
    CALL copy_sm_to_fm(matrix_s(1)%matrix,ov_fm)

!   Compute zij here
    DO iatom = 1, natom
      CALL cp_fm_set_all(zij_fm_set(iatom,1)%matrix, 0.0_dp, error=error)
      CALL cp_fm_get_info(zij_fm_set(iatom,1)%matrix,ncol_global=ldz)
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
!       CALL jacobi_rotation_pipek (  zij_fm_set, vectors, sweeps ) 
    CASE ( do_loc_crazy )
!       CALL crazy_rotations_pipek ( zij_fm_set, vectors, sweeps ) 
    CASE ( do_loc_direct )
!       CALL direct_mini_pipek     ( zij_fm_set, vectors, sweeps ) 
    CASE ( do_loc_none ) 
      IF (ionode) THEN
          scr=cp_logger_get_default_unit_nr()
          WRITE(scr,'(A,I6,A)') " No MOS localization applied "
      ENDIF
    CASE DEFAULT
       CALL stop_program("qs_localization","unknown localization method")
    END SELECT
    IF (ionode) THEN
       scr=cp_logger_get_default_unit_nr()
       WRITE(scr,'(A,I3,A,I6,A)') " Localization  for spin ",ispin ,&
            " converged in ",sweeps," iterations"
    ENDIF


    IF ( localized_wfn_control%print_centers)  THEN
!      CALL centers_spreads_pipek( localized_wfn_control, zij_fm_set, nmoloc, ispin, &
!           ionode, error=error ) 
    ENDIF 


    DEALLOCATE (first_sgf, last_sgf, nsgf, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL cp_fm_release(opvec)

    CALL timestop(0.0_dp,handle) 

  END SUBROUTINE optimize_loc_pipek

!******************************************************************************
!!****f* qs_loc_methods/qs_loc_driver()


!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE qs_loc_driver(qs_env,qs_loc_env,ionode,error)

    TYPE(qs_environment_type), POINTER           :: qs_env
    TYPE(qs_loc_env_new_type), POINTER           :: qs_loc_env
    LOGICAL, INTENT(IN)                          :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                   :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'qs_loc_driver', &
      routineP = module_name//':'//routine_name

    INTEGER :: i, ispin, imo, imoloc, istat, lb, loc_method, nao, nmosub, ub
    LOGICAL :: failure
    TYPE(cell_type), POINTER                     :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: moloc_coeff
    TYPE(cp_fm_type),  POINTER                   :: mo_coeff
    TYPE(dft_control_type), POINTER              :: dft_control
    TYPE(localized_wfn_control_type), POINTER    :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                    :: mos
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      POINTER                                    :: op_sm_set
    TYPE(cp_fm_p_type), DIMENSION(:,:), POINTER  :: op_fm_set
    REAL(dp), DIMENSION(6)                       :: weights
    REAL(KIND=dp), DIMENSION(:, :), POINTER      :: vecbuffer



    failure = .FALSE.
    NULLIFY(cell, localized_wfn_control, moloc_coeff, op_sm_set, op_fm_set)
    NULLIFY(mos, dft_control)

    CALL get_qs_loc_env(qs_loc_env=qs_loc_env,&
                        localized_wfn_control=localized_wfn_control,&
                        moloc_coeff=moloc_coeff,&
                        op_sm_set=op_sm_set,op_fm_set=op_fm_set,&
                        cell=cell, weights=weights)

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, mos=mos)
 
    DO ispin = 1,dft_control%nspins

       loc_method = localized_wfn_control%localization_method
       SELECT CASE (localized_wfn_control%operator_type)
       CASE(op_loc_berry)
         CALL optimize_loc_berry(loc_method, localized_wfn_control, &
              moloc_coeff(ispin)%matrix, op_sm_set, &
              op_fm_set, cell, weights, ispin, ionode, error=error)
 
       CASE(op_loc_boys) 

       CASE(op_loc_pipek)
         CALL  optimize_loc_pipek(qs_env, loc_method, localized_wfn_control,&
               moloc_coeff(ispin)%matrix, op_fm_set,  ispin, ionode, error=error)
       END SELECT 

       ! give back the localized orbitals
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, nao=nao, mo_coeff=mo_coeff)

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
            imo = localized_wfn_control%loc_states(i)

            CALL cp_fm_get_submatrix(moloc_coeff(ispin)%matrix,vecbuffer,1,imoloc,&
                 nao,1,transpose=.TRUE.,error=error)
            CALL cp_fm_set_submatrix(mo_coeff,vecbuffer,1,imo,&
                 nao,1,transpose=.TRUE.,error=error)
          END DO
          DEALLOCATE(vecbuffer,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       ELSE
         nmosub = localized_wfn_control%nloc_states(ispin)

         CALL cp_fm_to_fm(moloc_coeff(ispin)%matrix,mo_coeff,1,lb)

       END IF
    END DO  ! ispin

  END SUBROUTINE qs_loc_driver

!******************************************************************************

END MODULE qs_loc_methods
