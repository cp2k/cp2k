!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Apply the direct inversion in the iterative subspace (DIIS) of Pulay
!>      in the framework of an SCF iteration for convergence acceleration
!> \par Literature
!>      - P. Pulay, Chem. Phys. Lett. 73, 393 (1980)
!>      - P. Pulay, J. Comput. Chem. 3, 556 (1982)
!> \par History
!>      - Changed to BLACS matrix usage (08.06.2001,MK)
!>      - rewritten to include LSD (1st attempt) (01.2003, Joost VandeVondele)
!>      - DIIS for ROKS (05.04.06,MK)
!> \author Matthias Krack (28.06.2000)
! *****************************************************************************
MODULE qs_diis
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_init, &
       cp_dbcsr_maxabs, cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, &
       cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_transposed, cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_trace
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_maxabsval,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size
  USE mathlib,                         ONLY: diamat_all
  USE qs_diis_types,                   ONLY: qs_diis_buffer_type,&
                                             qs_diis_buffer_type_sparse
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_diis'
  INTEGER, SAVE               :: last_diis_b_id = 0

  ! Public subroutines

  PUBLIC :: qs_diis_b_clear,&
            qs_diis_b_create,&
            qs_diis_b_step
  PUBLIC :: qs_diis_b_clear_sparse, &
            qs_diis_b_create_sparse,&
            qs_diis_b_step_4lscf

CONTAINS

! *****************************************************************************
!> \brief Allocates an SCF DIIS buffer
!> \param diis_buffer the buffer to create
!> \param nbuffer ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE qs_diis_b_create(diis_buffer,nbuffer,error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    INTEGER, INTENT(in)                      :: nbuffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    ALLOCATE(diis_buffer,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
      NULLIFY (diis_buffer%b_matrix)
      NULLIFY (diis_buffer%error)
      NULLIFY (diis_buffer%parameter)
      diis_buffer%nbuffer = nbuffer
      diis_buffer%ncall = 0
      last_diis_b_id = last_diis_b_id + 1
      diis_buffer%id_nr = last_diis_b_id
      diis_buffer%ref_count = 1
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_diis_b_create

! *****************************************************************************
!> \brief Allocate and initialize a DIIS buffer for nao*nao parameter
!>      variables and with a buffer size of nbuffer.
!> \param diis_buffer the buffer to initialize
!> \param matrix_struct the structure for the matrix of the buffer
!> \param nspin ...
!> \param scf_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      - Creation (07.05.2001, Matthias Krack)
!>      - Changed to BLACS matrix usage (08.06.2001,MK)
!>      - DIIS for ROKS (05.04.06,MK)
!> \author Matthias Krack
!> \note
!>      check to allocate matrixes only when needed, using a linked list?
! *****************************************************************************
  SUBROUTINE qs_diis_b_check_i_alloc(diis_buffer,matrix_struct,nspin,&
                                     scf_section,error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    INTEGER, INTENT(IN)                      :: nspin
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_check_i_alloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ibuffer, ispin, &
                                                istat, nbuffer, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)

    nbuffer = diis_buffer%nbuffer

    IF (.NOT.ASSOCIATED(diis_buffer%error)) THEN
      ALLOCATE (diis_buffer%error(nbuffer,nspin),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "diis_buffer%error",0)

      DO ispin=1,nspin
        DO ibuffer=1,nbuffer
          NULLIFY (diis_buffer%error(ibuffer,ispin)%matrix)
          CALL cp_fm_create(diis_buffer%error(ibuffer,ispin)%matrix,&
                            name="qs_diis_b"//&
                                 TRIM(ADJUSTL(cp_to_string(diis_buffer%id_nr)))//&
                                 "%error("//&
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//","//&
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//")",&
                            matrix_struct=matrix_struct,error=error)
        END DO
      END DO
    END IF

    IF (.NOT.ASSOCIATED(diis_buffer%parameter)) THEN
      ALLOCATE (diis_buffer%parameter(nbuffer,nspin),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "diis_buffer%parameter",0)

      DO ispin=1,nspin
        DO ibuffer=1,nbuffer
          NULLIFY (diis_buffer%parameter(ibuffer,ispin)%matrix)
          CALL cp_fm_create(diis_buffer%parameter(ibuffer,ispin)%matrix,&
                            name="qs_diis_b"//&
                                 TRIM(ADJUSTL(cp_to_string(diis_buffer%id_nr)))//&
                                 "%parameter("//&
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//","//&
                                 TRIM(ADJUSTL(cp_to_string(ibuffer)))//")",&
                            matrix_struct=matrix_struct,error=error)
        END DO
      END DO
    END IF

    IF (.NOT.ASSOCIATED(diis_buffer%b_matrix)) THEN
      ALLOCATE (diis_buffer%b_matrix(nbuffer+1,nbuffer+1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "diis_buffer%b_matrix",&
                                       (nbuffer + 1)**2*dp_size)
      diis_buffer%b_matrix = 0.0_dp
      output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%DIIS_INFO",&
                                         extension=".scfLog",error=error)
      IF (output_unit > 0) THEN
        WRITE (UNIT=output_unit,FMT="(/,T9,A)")&
          "DIIS | The SCF DIIS buffer was allocated and initialized"
      END IF
      CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
                                        "PRINT%DIIS_INFO",error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_diis_b_check_i_alloc

! *****************************************************************************
!> \brief Update the SCF DIIS buffer, and if appropriate does a diis step.
!> \param diis_buffer ...
!> \param mo_array ...
!> \param kc ...
!> \param sc ...
!> \param delta ...
!> \param error_max ...
!> \param diis_step ...
!> \param eps_diis ...
!> \param nmixing ...
!> \param s_matrix ...
!> \param scf_section ...
!> \param roks ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      - Creation (07.05.2001, Matthias Krack)
!>      - Changed to BLACS matrix usage (08.06.2001, MK)
!>      - 03.2003 rewamped [fawzi]
!>      - Adapted for high-spin ROKS (08.04.06,MK)
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE qs_diis_b_step(diis_buffer,mo_array,kc,sc,delta,error_max,&
                            diis_step,eps_diis,nmixing,s_matrix,scf_section,roks,&
                            error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: kc
    TYPE(cp_fm_type), POINTER                :: sc
    REAL(KIND=dp), INTENT(IN)                :: delta
    REAL(KIND=dp), INTENT(OUT)               :: error_max
    LOGICAL, INTENT(OUT)                     :: diis_step
    REAL(KIND=dp), INTENT(IN)                :: eps_diis
    INTEGER, INTENT(IN), OPTIONAL            :: nmixing
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: s_matrix
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(IN), OPTIONAL            :: roks
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_step', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: eigenvalue_threshold = 1.0E-12_dp

    CHARACTER(LEN=2*default_string_length)   :: message
    INTEGER                                  :: handle, homo, ib, imo, ispin, &
                                                istat, jb, my_nmixing, nao, &
                                                nb, nb1, nmo, nspin, &
                                                output_unit
    LOGICAL                                  :: eigenvectors_discarded, &
                                                failure, mo_uocc, my_roks
    REAL(KIND=dp)                            :: maxocc, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ev, occ
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occa, occb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_fm_type), POINTER                :: c, new_errors, old_errors, &
                                                parameters
    TYPE(cp_logger_type), POINTER            :: logger

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    nspin = SIZE(mo_array)
    diis_step = .FALSE.
    failure = .FALSE.

    IF (PRESENT(roks)) THEN
      my_roks = .TRUE.
      nspin = 1
    ELSE
      my_roks = .FALSE.
    END IF

    my_nmixing = 2
    IF (PRESENT(nmixing)) my_nmixing = nmixing

    NULLIFY (c,new_errors,old_errors,parameters,matrix_struct,a,b,occa,occb)
    logger => cp_error_get_logger(error)

    ! Quick return, if no DIIS is requested

    IF (diis_buffer%nbuffer < 1) THEN
      CALL timestop(handle)
      RETURN
    END IF

    CALL cp_fm_get_info(kc(1)%matrix,&
                        matrix_struct=matrix_struct,&
                        error=error)
    CALL qs_diis_b_check_i_alloc(diis_buffer,&
                                 matrix_struct=matrix_struct,&
                                 nspin=nspin,&
                                 scf_section=scf_section,&
                                 error=error)

    error_max = 0.0_dp

    ib = MODULO(diis_buffer%ncall,diis_buffer%nbuffer) + 1
    diis_buffer%ncall = diis_buffer%ncall + 1
    nb = MIN(diis_buffer%ncall,diis_buffer%nbuffer)

    DO ispin=1,nspin

      CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                      nao=nao,&
                      nmo=nmo,&
                      homo=homo,&
                      mo_coeff=c,&
                      occupation_numbers=occa,&
                      uniform_occupation=mo_uocc,&
                      maxocc=maxocc)

      new_errors => diis_buffer%error(ib,ispin)%matrix
      parameters => diis_buffer%parameter(ib,ispin)%matrix

      ! Copy the Kohn-Sham matrix K to the DIIS buffer

      CALL cp_fm_to_fm(kc(ispin)%matrix,parameters,error=error)

      IF (my_roks) THEN

        ALLOCATE (occ(nmo),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "occ",nmo*dp_size)

        CALL get_mo_set(mo_set=mo_array(2)%mo_set,&
                        occupation_numbers=occb)

        DO imo=1,nmo
          occ(imo) = SQRT(occa(imo) + occb(imo))
        END DO

        CALL cp_fm_to_fm(c,sc,error=error)
        CALL cp_fm_column_scale(sc,occ(1:homo))

        ! KC <- K*C
        CALL cp_fm_symm("L","U",nao,homo,1.0_dp,parameters,sc,0.0_dp,kc(ispin)%matrix,error=error)

        IF (PRESENT(s_matrix)) THEN
          CALL copy_dbcsr_to_fm(s_matrix(1)%matrix,new_errors,error=error)
          ! SC <- S*C
          CALL cp_fm_symm("L","U",nao,homo,1.0_dp,new_errors,c,0.0_dp,sc,error=error)
          CALL cp_fm_column_scale(sc,occ(1:homo))
        END IF

        ! new_errors <- KC*(SC)^T - (SC)*(KC)^T = K*P*S - S*P*K
        ! or for an orthogonal basis
        ! new_errors <- KC*C^T - C*(KC)^T = K*P - P*K with S = I
        CALL cp_gemm("N","T",nao,nao,homo,1.0_dp,sc,kc(ispin)%matrix, 0.0_dp,new_errors,error=error)
        CALL cp_gemm("N","T",nao,nao,homo,1.0_dp,kc(ispin)%matrix,sc,-1.0_dp,new_errors,error=error)

        DEALLOCATE (occ,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"occ")

      ELSE

        ! KC <- K*C
        CALL cp_fm_symm("L","U",nao,homo,maxocc,parameters,c,0.0_dp,kc(ispin)%matrix,error=error)

        IF (PRESENT(s_matrix)) THEN
          ! I guess that this copy can be avoided for LSD
          CALL copy_dbcsr_to_fm(s_matrix(1)%matrix,new_errors,error=error)
          ! sc <- S*C
          CALL cp_fm_symm("L","U",nao,homo,2.0_dp,new_errors,c,0.0_dp,sc,error=error)
          ! new_errors <- KC*(SC)^T - (SC)*(KC)^T = K*P*S - S*P*K
          CALL cp_gemm("N","T",nao,nao,homo,1.0_dp,sc,kc(ispin)%matrix, 0.0_dp,new_errors,error=error)
          CALL cp_gemm("N","T",nao,nao,homo,1.0_dp,kc(ispin)%matrix,sc,-1.0_dp,new_errors,error=error)
        ELSE
          ! new_errors <- KC*(C)^T - C*(KC)^T = K*P - P*K
          CALL cp_gemm("N","T",nao,nao,homo,1.0_dp,c,kc(ispin)%matrix, 0.0_dp,new_errors,error=error)
          CALL cp_gemm("N","T",nao,nao,homo,1.0_dp,kc(ispin)%matrix,c,-1.0_dp,new_errors,error=error)
        END IF

      END IF

      CALL cp_fm_maxabsval(new_errors,tmp,error=error)
      error_max = MAX(error_max,tmp)

    END DO

    ! Check, if a DIIS step is appropiate

    diis_step = ((diis_buffer%ncall >= my_nmixing).AND.(delta < eps_diis))

    output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%DIIS_INFO",&
                                       extension=".scfLog",error=error)
    IF (output_unit > 0) THEN
       WRITE (UNIT=output_unit,FMT="(/,T9,A,I4,/,(T9,A,ES12.3))")&
         "DIIS | Current SCF DIIS buffer size:         ",nb,&
         "DIIS | Maximum SCF DIIS error vector element:",error_max,&
         "DIIS | Current SCF convergence:              ",delta,&
         "DIIS | Threshold value for a DIIS step:      ",eps_diis
      IF (error_max < eps_diis) THEN
        WRITE (UNIT=output_unit,FMT="(T9,A)")&
          "DIIS | => The SCF DIIS buffer will be updated"
      ELSE
        WRITE (UNIT=output_unit,FMT="(T9,A)")&
          "DIIS | => No update of the SCF DIIS buffer"
      END IF
      IF (diis_step.AND.(error_max < eps_diis)) THEN
        WRITE (UNIT=output_unit,FMT="(T9,A,/)")&
          "DIIS | => A SCF DIIS step will be performed"
      ELSE
        WRITE (UNIT=output_unit,FMT="(T9,A,/)")&
          "DIIS | => No SCF DIIS step will be performed"
      END IF
    END IF

    ! Update the SCF DIIS buffer

    IF (error_max < eps_diis) THEN

      b => diis_buffer%b_matrix

      DO jb=1,nb
        b(jb,ib) = 0.0_dp
        DO ispin=1,nspin
          old_errors => diis_buffer%error(jb,ispin)%matrix
          new_errors => diis_buffer%error(ib,ispin)%matrix
          CALL cp_fm_trace(old_errors,new_errors,tmp,error=error)
          b(jb,ib) = b(jb,ib) + tmp
        END DO
        b(ib,jb) = b(jb,ib)
      END DO

    ELSE

      diis_step = .FALSE.

    END IF

    ! Perform DIIS step

    IF (diis_step) THEN

      nb1 = nb + 1

      ALLOCATE (a(nb1,nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "a",nb1*nb1*dp_size)
      ALLOCATE (b(nb1,nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "b",nb1*nb1*dp_size)
      ALLOCATE (ev(nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "ev",nb1*dp_size)

      ! Set up the linear DIIS equation system

      b(1:nb,1:nb) = diis_buffer%b_matrix(1:nb,1:nb)

      b(1:nb,nb1) = -1.0_dp
      b(nb1,1:nb) = -1.0_dp
      b(nb1,nb1) = 0.0_dp

      ! Solve the linear DIIS equation system

      CALL diamat_all(b(1:nb1,1:nb1),ev(1:nb1),error=error)

      a(1:nb1,1:nb1) = b(1:nb1,1:nb1)

      eigenvectors_discarded = .FALSE.

      DO jb=1,nb1
        IF (ABS(ev(jb)) < eigenvalue_threshold) THEN
          IF (output_unit > 0) THEN
            IF (.NOT.eigenvectors_discarded) THEN
              WRITE (UNIT=output_unit,FMT="(T9,A)")&
                "DIIS | Checking eigenvalues of the DIIS error matrix"
            END IF
            WRITE (UNIT=message,FMT="(T9,A,I6,A,ES10.1,A,ES10.1)")&
              "DIIS | Eigenvalue ",jb," = ",ev(jb)," is smaller than "//&
              "threshold ",eigenvalue_threshold
            CALL compress(message)
            WRITE (UNIT=output_unit,FMT="(T9,A)") TRIM(message)
            eigenvectors_discarded = .TRUE.
          END IF
          a(1:nb1,jb) = 0.0_dp
        ELSE
          a(1:nb1,jb) = a(1:nb1,jb)/ev(jb)
        END IF
      END DO

      IF ((output_unit > 0).AND.eigenvectors_discarded) THEN
        WRITE (UNIT=output_unit,FMT="(T9,A,/)")&
          "DIIS | The corresponding eigenvectors were discarded"
      END IF

      ev(1:nb) = MATMUL(a(1:nb,1:nb1),b(nb1,1:nb1))

      ! Update Kohn-Sham matrix

      DO ispin=1,nspin
        CALL cp_fm_set_all(kc(ispin)%matrix,0.0_dp,error=error)
        DO jb=1,nb
          parameters => diis_buffer%parameter(jb,ispin)%matrix
          CALL cp_fm_scale_and_add(1.0_dp,kc(ispin)%matrix,-ev(jb),parameters,error=error)
        END DO
      END DO

      DEALLOCATE (a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"a")
      DEALLOCATE (b,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"b")
      DEALLOCATE (ev,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ev")

    ELSE

      DO ispin=1,nspin
        parameters => diis_buffer%parameter(ib,ispin)%matrix
        CALL cp_fm_to_fm(parameters,kc(ispin)%matrix,error=error)
      END DO

    END IF

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
                                      "PRINT%DIIS_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE qs_diis_b_step

! *****************************************************************************
!> \brief clears the buffer
!> \param diis_buffer the buffer to clear
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE qs_diis_b_clear(diis_buffer,error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_clear', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
      diis_buffer%ncall = 0
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_diis_b_clear

! *****************************************************************************
!> \brief Update the SCF DIIS buffer in linear scaling SCF (LS-SCF),
!>        and if appropriate does a diis step.
!> \param diis_buffer ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param unit_nr ...
!> \param iscf ...
!> \param diis_step ...
!> \param eps_diis ...
!> \param nmixing ...
!> \param s_matrix ...
!> \param threshold ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      - Adapted for LS-SCF (10-11-14) from qs_diis_b_step
!> \author Fredy W. Aquino
! *****************************************************************************

  SUBROUTINE qs_diis_b_step_4lscf(diis_buffer,qs_env,ls_scf_env,unit_nr,iscf, &
                                  diis_step,eps_diis,nmixing,s_matrix,threshold,error)
! Note.- Input: ls_scf_env%matrix_p(ispin) , Density Matrix
!               matrix_ks (from qs_env)    , Kohn-Sham Matrix  (IN/OUT)            

    TYPE(qs_diis_buffer_type_sparse), &
      POINTER                                :: diis_buffer
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    INTEGER, INTENT(IN)                      :: unit_nr, iscf
    LOGICAL, INTENT(OUT)                     :: diis_step
    REAL(KIND=dp), INTENT(IN)                :: eps_diis
    INTEGER, INTENT(IN), OPTIONAL            :: nmixing
    TYPE(cp_dbcsr_type), OPTIONAL            :: s_matrix
    REAL(KIND=dp), INTENT(IN)                :: threshold
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_step_4lscf', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: eigenvalue_threshold = 1.0E-12_dp

    INTEGER                                  :: handle, ib, ispin, istat, jb, &
                                                my_nmixing, nb, nb1, nspin
    REAL(KIND=dp)                            :: error_max, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ev
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type)                      :: matrix_KSerr_t, matrix_tmp
    TYPE(cp_dbcsr_type), POINTER             :: new_errors, old_errors, &
                                                parameters
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)
    nspin=ls_scf_env%nspins
    diis_step = .FALSE.
    my_nmixing = 2
    IF (PRESENT(nmixing)) my_nmixing = nmixing
    NULLIFY (new_errors,old_errors,parameters,a,b)
    logger => cp_error_get_logger(error)
    ! Quick return, if no DIIS is requested
    IF (diis_buffer%nbuffer < 1) THEN
      CALL timestop(handle)
      RETURN
    END IF

! Getting current Kohn-Sham matrix from qs_env
    CALL get_qs_env(qs_env,              &
                    para_env=para_env,   &
                    matrix_ks=matrix_ks, &
                    error=error) 
    CALL qs_diis_b_check_i_alloc_sparse( &
                           diis_buffer,  &
                           ls_scf_env,   &
                           nspin,        &
                           error)
    error_max = 0.0_dp

    ib = MODULO(diis_buffer%ncall,diis_buffer%nbuffer) + 1
    diis_buffer%ncall = diis_buffer%ncall + 1
    nb = MIN(diis_buffer%ncall,diis_buffer%nbuffer)
! Create scratch arrays
    CALL cp_dbcsr_init(matrix_tmp,error=error)
    CALL cp_dbcsr_create(matrix_tmp, &
                         template=ls_scf_env%matrix_ks(1),&
                         matrix_type='N',error=error)
    CALL cp_dbcsr_set(matrix_tmp,0.0_dp,error=error)     ! reset matrix
    CALL cp_dbcsr_init(matrix_KSerr_t,error=error)
    CALL cp_dbcsr_create(matrix_KSerr_t,&
                         template=ls_scf_env%matrix_ks(1),&
                         matrix_type='N',error=error)
    CALL cp_dbcsr_set(matrix_KSerr_t,0.0_dp,error=error) ! reset matrix

    DO ispin=1,nspin ! ------ Loop-ispin----START

      new_errors => diis_buffer%error(ib,ispin)%matrix
      parameters => diis_buffer%parameter(ib,ispin)%matrix
      ! Copy the Kohn-Sham matrix K to the DIIS buffer
      CALL cp_dbcsr_copy(parameters,              & ! out
                         matrix_ks(ispin)%matrix, & ! in
                         error=error)

        IF (PRESENT(s_matrix)) THEN ! if-s_matrix ---------- START
! Calculate Kohn-Sham error (non-orthogonal)= K*P*S-(K*P*S)^T
! matrix_tmp = P*S
            CALL cp_dbcsr_multiply("N", "N",         &
                   1.0_dp,ls_scf_env%matrix_p(ispin),& 
                          s_matrix,                  &
                   0.0_dp,matrix_tmp,                &
                   filter_eps=threshold, error=error)    
! new_errors= K*P*S    
            CALL cp_dbcsr_multiply("N", "N",         &
                   1.0_dp,matrix_ks(ispin)%matrix,   & 
                          matrix_tmp,                &
                   0.0_dp,new_errors,                &
                   filter_eps=threshold,error=error) 
! matrix_KSerr_t= transpose(K*P*S)        
            CALL cp_dbcsr_transposed(matrix_KSerr_t, &
                                     new_errors,     &
                                     error=error)
! new_errors=K*P*S-transpose(K*P*S)
            CALL cp_dbcsr_add(new_errors,       &
                              matrix_KSerr_t,   &
                              1.0_dp,-1.0_dp,   &
                              error=error)
        ELSE   ! if-s_matrix ---------- MID
! Calculate Kohn-Sham error (orthogonal)= K*P - P*K
! new_errors=K*P
            CALL cp_dbcsr_multiply("N", "N",         &
                   1.0_dp,matrix_ks(ispin)%matrix,   & 
                          ls_scf_env%matrix_p(ispin),&
                   0.0_dp,new_errors,                &
                   filter_eps=threshold,error=error)  
! matrix_KSerr_t= transpose(K*P)   
            CALL cp_dbcsr_transposed(matrix_KSerr_t, & 
                                     new_errors,     & 
                                     error=error)
! new_errors=K*P-transpose(K*P)
            CALL cp_dbcsr_add(new_errors,            & 
                              matrix_KSerr_t,        &
                              1.0_dp,-1.0_dp,        &
                              error=error)
        END IF ! if-s_matrix ---------- END

       tmp=cp_dbcsr_maxabs(new_errors)
       error_max = MAX(error_max,tmp)

    END DO ! ------ Loop-ispin----END

    ! Check, if a DIIS step is appropiate

    diis_step = (diis_buffer%ncall >= my_nmixing)

    IF (unit_nr>0) THEN
         WRITE(unit_nr,'(A29,I3,A3,4(I3,A1))') &
            "DIIS: (ncall,nbuffer,ib,nb)=(",iscf,")=(", &
            diis_buffer%ncall,",",diis_buffer%nbuffer,",",ib,",",nb,")"
         WRITE(unit_nr,'(A57,I3,A3,L1,A1,F10.8,A1,F4.2,A1,L1,A1)')     &
                  "DIIS: (diis_step,error_max,eps_diis,error_max<eps_diis)=(",& 
                  iscf,")=(",diis_step,",",error_max,",",eps_diis ,",", &
                  (error_max<eps_diis),")"
         WRITE(unit_nr,'(A75)') &
         "DIIS: diis_step=T : Perform DIIS  error_max<eps_diis=T : Update DIIS buffer"
    ENDIF

    ! Update the SCF DIIS buffer
    IF (error_max < eps_diis) THEN
      b => diis_buffer%b_matrix
      DO jb=1,nb
        b(jb,ib) = 0.0_dp
        DO ispin=1,nspin
          old_errors => diis_buffer%error(jb,ispin)%matrix
          new_errors => diis_buffer%error(ib,ispin)%matrix
          CALL cp_dbcsr_trace(old_errors, &
                              new_errors, &
                              tmp,        & ! out : < f_i | f_j >
                              error=error)
          b(jb,ib) = b(jb,ib) + tmp
        END DO ! end-loop-ispin
        b(ib,jb) = b(jb,ib)
      END DO ! end-loop-jb
    ELSE
      diis_step = .FALSE.
    END IF

    ! Perform DIIS step
    IF (diis_step) THEN
      nb1 = nb + 1
      ALLOCATE (a(nb1,nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "a",nb1*nb1*dp_size)
      ALLOCATE (b(nb1,nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "b",nb1*nb1*dp_size)
      ALLOCATE (ev(nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "ev",nb1*dp_size)
      ! Set up the linear DIIS equation system
      b(1:nb,1:nb) = diis_buffer%b_matrix(1:nb,1:nb)
      b(1:nb, nb1) = -1.0_dp
      b(nb1 ,1:nb) = -1.0_dp
      b(nb1 , nb1) =  0.0_dp
      ! Solve the linear DIIS equation system
      CALL diamat_all(b(1:nb1,1:nb1),ev(1:nb1),error=error)
      a(1:nb1,1:nb1) = b(1:nb1,1:nb1)
      DO jb=1,nb1
        IF (ABS(ev(jb)) < eigenvalue_threshold) THEN
          a(1:nb1,jb) = 0.0_dp
        ELSE        
          a(1:nb1,jb) = a(1:nb1,jb)/ev(jb)
        END IF
      END DO ! end-loop-jb

      ev(1:nb) = MATMUL(a(1:nb,1:nb1),b(nb1,1:nb1))

      ! Update Kohn-Sham matrix
      IF (iscf.GE.ls_scf_env%iter_ini_diis) THEN  ! if-iscf-to-updateKS------ START

         IF (unit_nr>0) THEN
          WRITE(unit_nr,'(A40,I3)') 'DIIS: Updating Kohn-Sham matrix at iscf=',iscf
         ENDIF

        DO ispin=1,nspin
          CALL cp_dbcsr_set(matrix_ks(ispin)%matrix, & ! reset matrix
                            0.0_dp,                  &
                            error=error)  
          DO jb=1,nb
             parameters => diis_buffer%parameter(jb,ispin)%matrix
             CALL cp_dbcsr_add(matrix_ks(ispin)%matrix,parameters, &
                               1.0_dp,-ev(jb),                     &
                               error=error)
          END DO ! end-loop-jb
        END DO ! end-loop-ispin
      ENDIF ! if-iscf-to-updateKS------ END

      DEALLOCATE (a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"a")
      DEALLOCATE (b,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"b")
      DEALLOCATE (ev,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ev")

    ELSE
      DO ispin=1,nspin
       parameters => diis_buffer%parameter(ib,ispin)%matrix
       CALL cp_dbcsr_copy(parameters,              & ! out
                          matrix_ks(ispin)%matrix, & ! in
                          error=error)
      ENDDO ! end-loop-ispin
    END IF
      CALL cp_dbcsr_release(matrix_tmp    ,error=error)
      CALL cp_dbcsr_release(matrix_KSerr_t,error=error)
      CALL timestop(handle)

  END SUBROUTINE qs_diis_b_step_4lscf

! *****************************************************************************
!> \brief Allocate and initialize a DIIS buffer with a buffer size of nbuffer.
!> \param diis_buffer the buffer to initialize
!> \param ls_scf_env ...
!> \param nspin ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      - Adapted from qs_diis_b_check_i_alloc for sparse matrices and
!>        used in LS-SCF module (ls_scf_main) (10-11-14)
!> \author Fredy W. Aquino
!> \note
!>      check to allocate matrices only when needed
! *****************************************************************************

  SUBROUTINE qs_diis_b_check_i_alloc_sparse(diis_buffer,ls_scf_env, &
                                            nspin,error)

    TYPE(qs_diis_buffer_type_sparse), &
      POINTER                                :: diis_buffer
    TYPE(ls_scf_env_type)                    :: ls_scf_env
    INTEGER, INTENT(IN)                      :: nspin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'qs_diis_b_check_i_alloc_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ibuffer, ispin, &
                                                istat, nbuffer
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)

    nbuffer = diis_buffer%nbuffer

    IF (.NOT.ASSOCIATED(diis_buffer%error)) THEN
      ALLOCATE (diis_buffer%error(nbuffer,nspin),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "diis_buffer%error",0)

      DO ispin=1,nspin
        DO ibuffer=1,nbuffer
          ALLOCATE (diis_buffer%error(ibuffer,ispin)%matrix,STAT=istat)
          CALL cp_dbcsr_init(diis_buffer%error(ibuffer,ispin)%matrix, &
                             error=error)

          CALL cp_dbcsr_create(diis_buffer%error(ibuffer,ispin)%matrix, &
                               template=ls_scf_env%matrix_ks(1),&
                               matrix_type='N',error=error)
        END DO
      END DO
    END IF

    IF (.NOT.ASSOCIATED(diis_buffer%parameter)) THEN
      ALLOCATE (diis_buffer%parameter(nbuffer,nspin),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "diis_buffer%parameter",0)

      DO ispin=1,nspin
        DO ibuffer=1,nbuffer
          ALLOCATE (diis_buffer%parameter(ibuffer,ispin)%matrix,STAT=istat)
          CALL cp_dbcsr_init(diis_buffer%parameter(ibuffer,ispin)%matrix, &
                             error=error)
          CALL cp_dbcsr_create(diis_buffer%parameter(ibuffer,ispin)%matrix, &
                               template=ls_scf_env%matrix_ks(1),&
                               matrix_type='N',error=error)
        END DO
      END DO
    END IF

    IF (.NOT.ASSOCIATED(diis_buffer%b_matrix)) THEN
      ALLOCATE (diis_buffer%b_matrix(nbuffer+1,nbuffer+1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "diis_buffer%b_matrix",&
                                       (nbuffer + 1)**2*dp_size)

      diis_buffer%b_matrix = 0.0_dp
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_diis_b_check_i_alloc_sparse

! *****************************************************************************
!> \brief clears the DIIS buffer in LS-SCF calculation
!> \param diis_buffer the buffer to clear
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10-11-14 created [FA] modified from qs_diis_b_clear
!> \author Fredy W. Aquino
! *****************************************************************************

  SUBROUTINE qs_diis_b_clear_sparse(diis_buffer,error)

    TYPE(qs_diis_buffer_type_sparse), &
      POINTER                                :: diis_buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_clear_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
      diis_buffer%ncall = 0
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_diis_b_clear_sparse

! *****************************************************************************
!> \brief Allocates an SCF DIIS buffer for LS-SCF calculation
!> \param diis_buffer the buffer to create
!> \param nbuffer ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10-11-14 created [FA] modified from qs_diis_b_create
!> \author Fredy W. Aquino
! *****************************************************************************
  SUBROUTINE qs_diis_b_create_sparse(diis_buffer,nbuffer,error)

    TYPE(qs_diis_buffer_type_sparse), &
      POINTER                                :: diis_buffer
    INTEGER, INTENT(in)                      :: nbuffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_create_sparse', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    ALLOCATE(diis_buffer,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
      NULLIFY (diis_buffer%b_matrix)
      NULLIFY (diis_buffer%error)
      NULLIFY (diis_buffer%parameter)
      diis_buffer%nbuffer = nbuffer
      diis_buffer%ncall = 0
      last_diis_b_id = last_diis_b_id + 1
      diis_buffer%id_nr = last_diis_b_id
      diis_buffer%ref_count = 1
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_diis_b_create_sparse

END MODULE qs_diis
