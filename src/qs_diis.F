!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****** cp2k/qs_diis [1.0] *
!!
!!   NAME
!!     qs_diis
!!
!!   FUNCTION
!!     Apply the direct inversion in the iterative subspace (DIIS) of Pulay
!!     in the framework of an SCF iteration for convergence acceleration
!!
!!   LITERATURE
!!     - P. Pulay, Chem. Phys. Lett. 73, 393 (1980)
!!     - P. Pulay, J. Comput. Chem. 3, 556 (1982)
!!
!!   AUTHOR
!!     Matthias Krack (28.06.2000)
!!
!!   MODIFICATION HISTORY
!!     - Changed to BLACS matrix usage (08.06.2001,MK)
!!     - rewritten to include LSD (1st attempt) (01.2003, Joost VandeVondele)
!!     - DIIS for ROKS (05.04.06,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_diis

  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_column_scale,&
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
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dp_size
  USE mathlib,                         ONLY: diamat_all
  USE qs_diis_types,                   ONLY: qs_diis_buffer_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: moduleN = "qs_diis"
  INTEGER, SAVE               :: last_diis_b_id = 0

  ! Public subroutines

  PUBLIC :: qs_diis_b_clear,&
            qs_diis_b_create,&
            qs_diis_b_step

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

!!****f* qs_diis/qs_diis_b_create [1.0] *
!!
!!   NAME
!!     qs_diis_b_create
!!
!!   FUNCTION
!!     Allocates an SCF DIIS buffer
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - diis_buffer: the buffer to create
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************

  SUBROUTINE qs_diis_b_create(diis_buffer,nbuffer,error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    INTEGER, INTENT(in)                      :: nbuffer
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "qs_diis_b_create", &
                                   routineP = moduleN//":"//routineN

    INTEGER                                  :: handle,stat
    LOGICAL                                  :: failure

    ! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

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

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_diis_b_create

  !***************************************************************************

!!****f* qs_diis/qs_diis_b_check_i_alloc [1.0] *
!!
!!   NAME
!!     qs_diis_b_check_i_alloc
!!
!!   FUNCTION
!!     Allocate and initialize a DIIS buffer for nao*nao parameter
!!     variables and with a buffer size of nbuffer.
!!
!!   NOTES
!!     check to allocate matrixes only when needed, using a linked list?
!!
!!   ARGUMENTS
!!     - diis_buffer: the buffer to initialize
!!     - matrix_struct: the structure for the matrix of the buffer
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     - Creation (07.05.2001, Matthias Krack)
!!     - Changed to BLACS matrix usage (08.06.2001,MK)
!!     - DIIS for ROKS (05.04.06,MK)
!!
!!*** *********************************************************************

  SUBROUTINE qs_diis_b_check_i_alloc(diis_buffer,matrix_struct,nspin,&
                                     scf_section,error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    INTEGER, INTENT(IN)                      :: nspin
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "qs_diis_b_check_i_alloc",&
                                   routineP = moduleN//":"//routineN

    INTEGER                                  :: handle, ibuffer, ispin, istat, &
                                                nbuffer, output_unit
    LOGICAL                                  :: failure, ionode
    TYPE(cp_logger_type), POINTER            :: logger

    ! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.

    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)

    ionode = (logger%para_env%source == logger%para_env%mepos)

    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger,local=.FALSE.)
    ELSE
      output_unit= -1
    END IF

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
                            matrix_struct=matrix_struct)
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
                            matrix_struct=matrix_struct)
        END DO
      END DO
    END IF

    IF (.NOT.ASSOCIATED(diis_buffer%b_matrix)) THEN
      ALLOCATE (diis_buffer%b_matrix(nbuffer+1,nbuffer+1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "diis_buffer%b_matrix",&
                                       (nbuffer + 1)**2*dp_size)
      diis_buffer%b_matrix = 0.0_dp
      IF (ionode) THEN
        output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%DIIS_INFO",&
                                           extension=".scfLog",error=error)
      END IF
      IF (output_unit > 0) THEN
        WRITE (UNIT=output_unit,FMT="(/,T3,A)")&
          "DIIS | The SCF DIIS buffer was allocated and initialized"
      END IF
      CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
                                        "PRINT%DIIS_INFO",error=error)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_diis_b_check_i_alloc

  ! *****************************************************************************

!!****f* qs_diis/qs_diis_b_step [1.0] *
!!
!!   NAME
!!     qs_diis_b_step
!!
!!   FUNCTION
!!     Update the SCF DIIS buffer, and if appropriate does a diis step.
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     to do
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     - Creation (07.05.2001, Matthias Krack)
!!     - Changed to BLACS matrix usage (08.06.2001, MK)
!!     - 03.2003 rewamped [fawzi]
!!     - Adapted for high-spin ROKS (08.04.06,MK)
!!
!!*** *********************************************************************

  SUBROUTINE qs_diis_b_step(diis_buffer,mo_array,kc,sc,delta,error_max,&
                            diis_step,eps_diis,s_matrix,scf_section,roks,&
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
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: s_matrix
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(IN), OPTIONAL            :: roks
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "qs_diis_b_step",&
                                   routineP = moduleN//":"//routineN
    REAL(KIND=dp), PARAMETER    :: eigenvalue_threshold = 1.0E-12_dp

    CHARACTER(LEN=2*default_string_length)   :: message
    INTEGER                                  :: handle, homo, ib, imo, &
                                                ispin, istat, jb, &
                                                nao, nb, nb1, nmo, nspin, &
                                                output_unit
    LOGICAL                                  :: failure, ionode, mo_uocc, &
                                                my_roks, eigenvectors_discarded
    REAL(KIND=dp)                            :: maxocc, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ev
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occa,occb
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: occ
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_fm_type), POINTER                :: c, new_errors, old_errors, &
                                                parameters
    TYPE(cp_logger_type), POINTER            :: logger

    ! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    nspin = SIZE(mo_array)
    diis_step = .FALSE.
    failure = .FALSE.

    IF (PRESENT(roks)) THEN
      my_roks = .TRUE.
      nspin = 1
    ELSE
      my_roks = .FALSE.
    END IF

    NULLIFY (c,new_errors,old_errors,parameters,matrix_struct,a,b,occa,occb)
    logger => cp_error_get_logger(error)

    ionode = (logger%para_env%source == logger%para_env%mepos)

    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger,local=.FALSE.)
    ELSE
      output_unit= -1
    END IF

    ! Quick return, if no DIIS is requested

    IF (diis_buffer%nbuffer < 1) THEN
      CALL timestop(0.0_dp,handle)
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

      CALL cp_fm_to_fm(kc(ispin)%matrix,parameters)

      IF (my_roks) THEN

        ALLOCATE (occ(nmo),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "occ",nmo*dp_size)

        CALL get_mo_set(mo_set=mo_array(2)%mo_set,&
                        occupation_numbers=occb)

        DO imo=1,nmo
          occ(imo) = SQRT(occa(imo) + occb(imo))
        END DO

        CALL cp_fm_to_fm(c,sc)
        CALL cp_fm_column_scale(sc,occ(1:homo))

        ! KC <- K*C
        CALL cp_fm_symm("L","U",nao,homo,1.0_dp,parameters,sc,0.0_dp,kc(ispin)%matrix)

        IF (PRESENT(s_matrix)) THEN
          CALL copy_sm_to_fm(s_matrix(1)%matrix,new_errors)
          ! SC <- S*C
          CALL cp_fm_symm("L","U",nao,homo,1.0_dp,new_errors,c,0.0_dp,sc)
          CALL cp_fm_column_scale(sc,occ(1:homo))
        END IF

        ! new_errors <- KC*(SC)^T - (SC)*(KC)^T = K*P*S - S*P*K
        ! or for an orthogonal basis
        ! new_errors <- KC*C^T - C*(KC)^T = K*P - P*K with S = I
        CALL cp_fm_gemm("N","T",nao,nao,homo,1.0_dp,sc,kc(ispin)%matrix, 0.0_dp,new_errors)
        CALL cp_fm_gemm("N","T",nao,nao,homo,1.0_dp,kc(ispin)%matrix,sc,-1.0_dp,new_errors)

        DEALLOCATE (occ,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"occ")

      ELSE

!       Using the not integer occ. numbers makes the algorithm instable
!       IF(.NOT. mo_uocc) THEN
!         CALL cp_fm_to_fm(c,sc)
!         CALL cp_fm_column_scale(sc,occa(1:homo))
!         CALL cp_fm_symm("L","U",nao,homo,1.0_dp,parameters,sc,0.0_dp,kc(ispin)%matrix)
!       ELSE
          ! KC <- K*C
          CALL cp_fm_symm("L","U",nao,homo,maxocc,parameters,c,0.0_dp,kc(ispin)%matrix)
!       END IF

        IF (PRESENT(s_matrix)) THEN
          ! I guess that this copy can be avoided for LSD
          CALL copy_sm_to_fm(s_matrix(1)%matrix,new_errors)
          ! sc <- S*C
          CALL cp_fm_symm("L","U",nao,homo,2.0_dp,new_errors,c,0.0_dp,sc)
          ! new_errors <- KC*(SC)^T - (SC)*(KC)^T = K*P*S - S*P*K
          CALL cp_fm_gemm("N","T",nao,nao,homo,1.0_dp,sc,kc(ispin)%matrix, 0.0_dp,new_errors)
          CALL cp_fm_gemm("N","T",nao,nao,homo,1.0_dp,kc(ispin)%matrix,sc,-1.0_dp,new_errors)
        ELSE
          ! new_errors <- KC*(C)^T - C*(KC)^T = K*P - P*K
          CALL cp_fm_gemm("N","T",nao,nao,homo,1.0_dp,c,kc(ispin)%matrix, 0.0_dp,new_errors)
          CALL cp_fm_gemm("N","T",nao,nao,homo,1.0_dp,kc(ispin)%matrix,c,-1.0_dp,new_errors)
        END IF

      END IF

      CALL cp_fm_maxabsval(new_errors,tmp)
      error_max = MAX(error_max,tmp)

    END DO

    ! Check, if a DIIS step is appropiate

    diis_step = ((diis_buffer%ncall > 1).AND.(delta < eps_diis))

    IF (ionode) THEN
      output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%DIIS_INFO",&
                                         extension=".scfLog",error=error)
    END IF
    IF (output_unit > 0) THEN
       WRITE (UNIT=output_unit,FMT="(/,T3,A,I4,/,(T3,A,ES12.3))")&
         "DIIS | Current SCF DIIS buffer size:         ",nb,&
         "DIIS | Maximum SCF DIIS error vector element:",error_max,&
         "DIIS | Current SCF convergence:              ",delta,&
         "DIIS | Threshold value for a DIIS step:      ",eps_diis
      IF (error_max < eps_diis) THEN
        WRITE (UNIT=output_unit,FMT="(T3,A)")&
          "DIIS | => The SCF DIIS buffer will be updated"
      ELSE
        WRITE (UNIT=output_unit,FMT="(T3,A)")&
          "DIIS | => No update of the SCF DIIS buffer"
      END IF
      IF (diis_step.AND.(error_max < eps_diis)) THEN
        WRITE (UNIT=output_unit,FMT="(T3,A,/)")&
          "DIIS | => A SCF DIIS step will be performed"
      ELSE
        WRITE (UNIT=output_unit,FMT="(T3,A,/)")&
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
          CALL cp_fm_trace(old_errors,new_errors,tmp)
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

      CALL diamat_all(b(1:nb1,1:nb1),ev(1:nb1))

      a(1:nb1,1:nb1) = b(1:nb1,1:nb1)

      eigenvectors_discarded = .FALSE.

      DO jb=1,nb1
        IF (ABS(ev(jb)) < eigenvalue_threshold) THEN
          IF (output_unit > 0) THEN
            IF (.NOT.eigenvectors_discarded) THEN
              WRITE (UNIT=output_unit,FMT="(T3,A)")&
                "DIIS | Checking eigenvalues of the DIIS error matrix"
            END IF
            WRITE (UNIT=message,FMT="(A,I6,A,ES10.1,A,ES10.1)")&
              "DIIS | Eigenvalue ",jb," = ",ev(jb)," is smaller than "//&
              "the requested threshold ",eigenvalue_threshold
            CALL compress(message)
            WRITE (UNIT=output_unit,FMT="(T3,A)") TRIM(message)
            eigenvectors_discarded = .TRUE.
          END IF
          a(1:nb1,jb) = 0.0_dp
        ELSE
          a(1:nb1,jb) = a(1:nb1,jb)/ev(jb)
        END IF
      END DO

      IF ((output_unit > 0).AND.eigenvectors_discarded) THEN
        WRITE (UNIT=output_unit,FMT="(T3,A,/)")&
          "DIIS | The corresponding eigenvectors were discarded"
      END IF

      ev(1:nb) = MATMUL(a(1:nb,1:nb1),b(nb1,1:nb1))

      ! Update Kohn-Sham matrix

      DO ispin=1,nspin
        CALL cp_fm_set_all(kc(ispin)%matrix,0.0_dp)
        DO jb=1,nb
          parameters => diis_buffer%parameter(jb,ispin)%matrix
          CALL cp_fm_scale_and_add(1.0_dp,kc(ispin)%matrix,-ev(jb),parameters)
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
        CALL cp_fm_to_fm(parameters,kc(ispin)%matrix)
      END DO

    END IF

    IF (output_unit > 0) THEN
      CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
                                        "PRINT%DIIS_INFO",error=error)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_diis_b_step

  ! ***************************************************************************

!!****f* qs_diis/qs_diis_b_clear [1.0] *
!!
!!   NAME
!!     qs_diis_b_clear
!!
!!   FUNCTION
!!     clears the buffer
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - diis_buffer: the buffer to clear
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************

  SUBROUTINE qs_diis_b_clear(diis_buffer,error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "qs_diis_b_clear", &
                                   routineP = moduleN//":"//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    ! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
      diis_buffer%ncall = 0
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_diis_b_clear

  ! ***************************************************************************

END MODULE qs_diis
