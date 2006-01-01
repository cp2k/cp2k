!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!


!!****** cp2k/qs_diis [1.0] *
!!
!!   NAME
!!     qs_diis
!!
!!   FUNCTION
!!     Perform a direct inversion in the iterative subspace (DIIS).
!!
!!   AUTHOR
!!     Matthias Krack (28.06.2000)
!!
!!   MODIFICATION HISTORY
!!     - Changed to BLACS matrix usage (08.06.2001,MK)
!!     - rewritten to include LSD (1st attempt) (01.2003, Joost VandeVondele)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_diis
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
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
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathlib,                         ONLY: diamat_all
  USE qs_diis_types,                   ONLY: qs_diis_buffer_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="qs_diis"
  INTEGER, SAVE, PRIVATE :: last_diis_b_id=0

  ! *** Public subroutines ***

  PUBLIC :: qs_diis_b_step, qs_diis_b_create, qs_diis_b_clear
!!***
  ! *****************************************************************************

CONTAINS

!!****f* qs_diis/qs_diis_b_create [1.0] *
!!
!!   NAME
!!     qs_diis_b_create
!!
!!   SYNOPSIS
!!     Subroutine qs_diis_b_create(diis_buffer, nbuffer, error)
!!       Type(qs_diis_buffer_type), Pointer:: diis_buffer
!!       Integer, Intent (IN):: nbuffer
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_diis_b_create
!!
!!   FUNCTION
!!     allocates an scf buffer
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

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(diis_buffer,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       diis_buffer%nbuffer=nbuffer
       diis_buffer%ncall=0
       last_diis_b_id=last_diis_b_id+1
       diis_buffer%id_nr=last_diis_b_id
       diis_buffer%ref_count=1
       NULLIFY(diis_buffer%error,diis_buffer%parameter, diis_buffer%b_matrix)
    END IF
  END SUBROUTINE qs_diis_b_create
  !***************************************************************************

!!****f* qs_diis/qs_diis_b_check_i_alloc [1.0] *
!!
!!   NAME
!!     qs_diis_b_check_i_alloc
!!
!!   FUNCTION
!!     Allocate and initialize a DIIS buffer for "nao*nao" parameter
!!     variables and with a buffer size of "nbuffer".
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
!!     - Changed to BLACS matrix usage (08.06.2001, MK)
!!
!!*** *********************************************************************
  SUBROUTINE qs_diis_b_check_i_alloc(diis_buffer,&
       matrix_struct,nspin,error)

    TYPE(qs_diis_buffer_type), POINTER       :: diis_buffer
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    INTEGER, INTENT(IN)                      :: nspin
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_check_i_alloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ibuffer, ispin, istat, nbuffer
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)

    nbuffer=diis_buffer%nbuffer

    IF (.NOT.ASSOCIATED(diis_buffer%error)) THEN
       ALLOCATE (diis_buffer%error(nbuffer,nspin),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"diis_buffer%error",0)

       DO ispin=1,nspin
          DO ibuffer=1,nbuffer
             NULLIFY (diis_buffer%error(ibuffer,ispin)%matrix)
             CALL cp_fm_create(&
                  diis_buffer%error(ibuffer,ispin)%matrix,&
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
       IF (istat /= 0) CALL stop_memory(routineP,"diis_buffer%parameter",0)

       DO ispin=1,nspin
          DO ibuffer=1,nbuffer
             NULLIFY (diis_buffer%parameter(ibuffer,ispin)%matrix)
             CALL cp_fm_create(&
                  diis_buffer%parameter(ibuffer,ispin)%matrix,&
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
       IF (istat /= 0) CALL stop_memory(routineP,"diis_buffer%b_matrix",&
            (nbuffer + 1)**2*dp_size)
       diis_buffer%b_matrix=0.0_dp
       IF (logger%para_env%source==logger%para_env%mepos.AND.&
            logger%print_keys%diis_information) THEN
          WRITE (UNIT=cp_logger_get_default_unit_nr(logger),FMT="(/,T2,A)")&
               "The SCF DIIS buffer was allocated and initialised"
       END IF
    END IF

  END SUBROUTINE qs_diis_b_check_i_alloc
  ! *****************************************************************************

!!****f* qs_diis/qs_diis_b_step [1.0] *
!!
!!   NAME
!!     qs_diis_b_step
!!
!!   SYNOPSIS
!!     Subroutine qs_diis_b_step(diis_buffer, mo_array, kc, sc, delta,&
!!         error_max, diis_step, qs_env, error)
!!       Type(qs_diis_buffer_type), Pointer:: diis_buffer
!!       Type(cp_fm_type), Pointer:: sc
!!       Type(cp_fm_p_type), Dimension(:), Pointer:: kc
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!       Type(mo_set_p_type), Dimension(:), Pointer:: mo_array
!!       Real(KIND = dp), Intent (IN):: delta
!!       Real(KIND = dp), Intent (OUT):: error_max
!!       Logical, Intent (OUT):: diis_step
!!     End Subroutine qs_diis_b_step
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
!!
!!*** *********************************************************************
  SUBROUTINE qs_diis_b_step(diis_buffer,mo_array,kc,sc,delta,error_max,&
       diis_step,eps_diis,s_matrix,error)

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
      POINTER, OPTIONAL                      :: s_matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_diis_b_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ib, ispin, istat, jb, &
                                                nao, nb, nb1, nmo, nspin, &
                                                output_unit
    LOGICAL                                  :: failure, ionode, mo_uocc
    REAL(KIND=dp)                            :: tmp, maxocc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ev
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a, b
    TYPE(cp_fm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_fm_type), POINTER                :: c, new_errors, old_errors, &
                                                parameters
    TYPE(cp_logger_type), POINTER            :: logger

!   ---------------------------------------------------------------------------

    CALL timeset("scf_diis","I","",handle)

    nspin=SIZE(mo_array)
    failure=.FALSE.
    NULLIFY(c,new_errors,old_errors,parameters,matrix_struct,a,b,occupation)
    logger => cp_error_get_logger(error)

    diis_step = .FALSE.

    ionode = logger%para_env%source==logger%para_env%mepos
    IF (ionode) THEN
        output_unit= cp_logger_get_default_unit_nr(logger,local=.FALSE.)
    ELSE
        output_unit= 0
    ENDIF

    !   *** Quick return, if no DIIS is requested ***
    IF ( diis_buffer%nbuffer < 1) THEN
       CALL timestop(0.0_dp,handle)
       RETURN
    END IF

    CALL cp_fm_get_info(kc(1)%matrix,matrix_struct=matrix_struct,error=error)
    CALL qs_diis_b_check_i_alloc(diis_buffer,matrix_struct=matrix_struct,&
         nspin=nspin)

    error_max=0.0_dp

    ib = MODULO(diis_buffer%ncall,diis_buffer%nbuffer) + 1
    diis_buffer%ncall = diis_buffer%ncall + 1
    nb = MIN(diis_buffer%ncall,diis_buffer%nbuffer)

    DO ispin=1,nspin

       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
            nao=nao,&
            homo=nmo,&
            mo_coeff=c,&
            occupation_numbers=occupation,& 
            uniform_occupation= mo_uocc,&
            maxocc = maxocc)


       parameters => diis_buffer%parameter(ib,ispin)%matrix
       CALL cp_fm_to_fm(kc(ispin)%matrix,parameters)
!       Using the not integer occ. numbers makes the algorithm instable
!       IF(.NOT. mo_uocc) THEN
!         CALL cp_fm_to_fm(c,sc)
!         CALL cp_fm_column_scale(sc,occupation(1:nmo))
!         CALL cp_fm_symm("L","U",nao,nmo,1.0_dp,parameters,sc,0.0_dp,kc(ispin)%matrix)
!       ELSE
         CALL cp_fm_symm("L","U",nao,nmo,maxocc,parameters,c,0.0_dp,kc(ispin)%matrix)
!       END IF

       IF (PRESENT(s_matrix)) THEN
         new_errors => diis_buffer%error(ib,ispin)%matrix
         ! I guess that this copy can be avoided for LSD
         CALL copy_sm_to_fm(s_matrix(1)%matrix,new_errors)
         CALL cp_fm_symm("L","U",nao,nmo,2.0_dp,new_errors,c,0.0_dp,sc)

         CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_dp,sc,kc(ispin)%matrix, 0.0_dp,new_errors)
         CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_dp,kc(ispin)%matrix,sc,-1.0_dp,new_errors)
       ELSE
         new_errors => diis_buffer%error(ib,ispin)%matrix
         CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_dp,c,kc(ispin)%matrix, 0.0_dp,new_errors)
         CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_dp,kc(ispin)%matrix,c,-1.0_dp,new_errors)
       END IF

       CALL cp_fm_maxabsval(new_errors,tmp)
       error_max=MAX(error_max,tmp)

    ENDDO

    IF (ionode.AND.logger%print_keys%diis_information) THEN
       WRITE (UNIT=output_unit,FMT="(/,T2,A,E12.3)")&
            "Maximum SCF DIIS error vector element:",error_max
    END IF

    !   *** Check, if a DIIS step is appropiate ***

    diis_step = ((diis_buffer%ncall > 1).AND.(delta < eps_diis))

    !   *** Update the SCF DIIS buffer ***

    IF (error_max < eps_diis) THEN

       b => diis_buffer%b_matrix

       DO jb=1,nb
          b(jb,ib)=0.0_dp
          DO ispin=1,nspin
             old_errors => diis_buffer%error(jb,ispin)%matrix
             new_errors => diis_buffer%error(ib,ispin)%matrix
             CALL cp_fm_trace(old_errors,new_errors,tmp)
             b(jb,ib)=b(jb,ib)+tmp
          ENDDO
          b(ib,jb) = b(jb,ib)
       END DO

       IF (ionode.AND.logger%print_keys%diis_information) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
               "The SCF DIIS buffer was updated"
       END IF

    ELSE

       diis_step = .FALSE.

    END IF

    !   *** Perform DIIS update ***

    IF (diis_step) THEN

       nb1 = nb + 1

       ALLOCATE (a(nb1,nb1),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"a",nb1*nb1*dp_size)
       ALLOCATE (b(nb1,nb1),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"b",nb1*nb1*dp_size)
       ALLOCATE (ev(nb1),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"ev",nb1*dp_size)

       !     *** Set up the linear DIIS equation system ***

       b(1:nb,1:nb) = diis_buffer%b_matrix(1:nb,1:nb)

       b(1:nb,nb1) = -1.0_dp
       b(nb1,1:nb) = -1.0_dp
       b(nb1,nb1) = 0.0_dp

       !     *** Solve the linear DIIS equation system ***

       CALL diamat_all(b(1:nb1,1:nb1),ev(1:nb1))

       a(1:nb1,1:nb1) = b(1:nb1,1:nb1)

       DO jb=1,nb1
          IF (ABS(ev(jb)) < 1.0E-12_dp) THEN
             a(1:nb1,jb) = 0.0_dp
          ELSE
             a(1:nb1,jb) = a(1:nb1,jb)/ev(jb)
          END IF
       END DO

       ev(1:nb) = MATMUL(a(1:nb,1:nb1),b(nb1,1:nb1))

       !     *** Update Kohn-Sham matrix ***

       DO ispin=1,nspin
          CALL cp_fm_set_all(kc(ispin)%matrix,0.0_dp)
          DO jb=1,nb
             parameters => diis_buffer%parameter(jb,ispin)%matrix
             CALL cp_fm_scale_and_add(1.0_dp,kc(ispin)%matrix,-ev(jb),parameters)
          END DO
       ENDDO

       DEALLOCATE (a,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"a")
       DEALLOCATE (b,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"b")
       DEALLOCATE (ev,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineP,"ev")

    ELSE

       DO ispin=1,nspin
          parameters => diis_buffer%parameter(ib,ispin)%matrix
          CALL cp_fm_to_fm(parameters,kc(ispin)%matrix)
       ENDDO

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_diis_b_step
  ! *****************************************************************************

!!****f* qs_diis/qs_diis_b_clear [1.0] *
!!
!!   NAME
!!     qs_diis_b_clear
!!
!!   SYNOPSIS
!!     Subroutine qs_diis_b_clear(diis_buffer, error)
!!       Type(qs_diis_buffer_type), Pointer:: diis_buffer
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_diis_b_clear
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

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_diis_b_clear', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(diis_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(diis_buffer%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       diis_buffer%ncall=0
    END IF
  END SUBROUTINE qs_diis_b_clear
  !***************************************************************************

END MODULE qs_diis
