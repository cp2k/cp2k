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
!!
!!   SOURCE
!******************************************************************************

MODULE qs_diis

  USE kinds, ONLY: wp => dp,&
                   wp_size => dp_size

  USE qs_blacs,             ONLY: blacs_add,&
                                  blacs_gemm,&
                                  blacs_matrix_p_type,&
                                  blacs_matrix_type,&
                                  blacs_maxval,&
                                  blacs_set_all,&
                                  blacs_symm,&
                                  blacs_trace,&
                                  copy_blacs_to_blacs_matrix,&
                                  copy_sparse_to_blacs_matrix,&
                                  get_blacs_matrix_info,&
                                  replicate_blacs_matrix
  USE qs_environment_types, ONLY: get_qs_env,&
                                  qs_environment_type
  USE global_types,         ONLY: global_environment_type
  USE mathlib,              ONLY: diagonalize_matrix,&
                                  symmetrize_matrix
  USE qs_mo_types,          ONLY: get_mo_set,&
                                  mo_set_type
  USE sparse_matrix_types,  ONLY: real_matrix_p_type
  USE termination,          ONLY: stop_memory
  USE timings,              ONLY: timeset,&
                                  timestop

  IMPLICIT NONE

  PRIVATE

  TYPE diis_buffer_type
    PRIVATE
    INTEGER                                          :: nbuffer,ncall
    TYPE(blacs_matrix_p_type), DIMENSION(:), POINTER :: error,parameter
    REAL(wp), DIMENSION(:,:), POINTER                :: b_matrix
  END TYPE diis_buffer_type

  TYPE(blacs_matrix_type), POINTER :: new_errors,old_errors,parameters

  TYPE(diis_buffer_type) :: scf_diis_buffer
  REAL(wp)               :: eps_diis = 0.1_wp
  INTEGER                :: max_diis = 4
  LOGICAL                :: init_scf_diis_buffer_done = .FALSE.

! *** Public variables ***

  PUBLIC :: eps_diis,max_diis

! *** Public subroutines ***

  PUBLIC :: scf_diis

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_scf_diis_buffer(nbuffer,nao,prototype_matrix)

!   Purpose: Allocate and initialize a DIIS buffer for "nao*nao" parameter
!            variables and with a buffer size of "nbuffer".

!   History: - Creation (07.05.2001, Matthias Krack)
!            - Changed to BLACS matrix usage (08.06.2001, MK)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER :: prototype_matrix
    INTEGER, INTENT(IN)              :: nbuffer,nao

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_scf_diis_buffer (MODULE qs_diis)"

!   *** Local variables ***

    INTEGER :: ibuffer,istat

!   ---------------------------------------------------------------------------

    scf_diis_buffer%nbuffer = nbuffer
    scf_diis_buffer%ncall = 0

    ALLOCATE (scf_diis_buffer%error(nbuffer),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scf_diis_buffer%error",0)
    ALLOCATE (scf_diis_buffer%parameter(nbuffer),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scf_diis_buffer%parameter",0)

    DO ibuffer=1,nbuffer
      CALL replicate_blacs_matrix(&
        prototype_matrix=prototype_matrix,&
        new_matrix=scf_diis_buffer%error(ibuffer)%blacs_matrix,&
        name="SCF DIIS ERROR MATRIX")
      CALL replicate_blacs_matrix(&
        prototype_matrix=prototype_matrix,&
        new_matrix=scf_diis_buffer%parameter(ibuffer)%blacs_matrix,&
        name="SCF DIIS PARAMETER MATRIX")
    END DO

    ALLOCATE (scf_diis_buffer%b_matrix(nbuffer+1,nbuffer+1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scf_diis_buffer%b_matrix",&
                                     (nbuffer + 1)**2*wp_size)

    init_scf_diis_buffer_done = .TRUE.

  END SUBROUTINE allocate_scf_diis_buffer

! *****************************************************************************

  SUBROUTINE scf_diis(mo_set,kc,sc,delta,error_max,diis_step,qs_env,globenv)

!   Purpose: Update the SCF DIIS buffer.

!   History: - Creation (07.05.2001, Matthias Krack)
!            - Changed to BLACS matrix usage (08.06.2001, MK)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: kc,sc
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    REAL(wp), INTENT(IN)                      :: delta
    REAL(wp), INTENT(OUT)                     :: error_max
    LOGICAL, INTENT(OUT)                      :: diis_step

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE scf_diis (MODULE qs_diis)"

!   *** Local variables ***

    TYPE(blacs_matrix_type), POINTER  :: c

    INTEGER :: handle,ib,istat,jb,nao,nb,nb1,nmo,nparameter,output_unit
    LOGICAL :: ionode

    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s

    REAL(wp), DIMENSION(:), ALLOCATABLE :: ev
    REAL(wp), DIMENSION(:,:), POINTER   :: a,b

!   ---------------------------------------------------------------------------

    CALL timeset("scf_diis","I","",handle)

    diis_step = .FALSE.

!   *** Quick return, if no DIIS is requested ***

    IF (max_diis < 1) RETURN

    ionode = globenv%ionode
    output_unit= globenv%scr

    CALL get_qs_env(qs_env=qs_env,s=s)

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    homo=nmo,&
                    eigenvectors=c)

    nparameter = nao*nao

    IF (.NOT.init_scf_diis_buffer_done) THEN
      CALL allocate_scf_diis_buffer(max_diis,nao,kc)
      IF (ionode.AND.globenv%print%diis_information) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "The SCF DIIS buffer was allocated and initialised"
      END IF
    END IF

    ib = MODULO(scf_diis_buffer%ncall,scf_diis_buffer%nbuffer) + 1
    scf_diis_buffer%ncall = scf_diis_buffer%ncall + 1
    nb = MIN(scf_diis_buffer%ncall,scf_diis_buffer%nbuffer)

    parameters => scf_diis_buffer%parameter(ib)%blacs_matrix
    CALL copy_blacs_to_blacs_matrix(kc,parameters)
    CALL blacs_symm("L","U",nao,nmo,2.0_wp,parameters,c,0.0_wp,kc,globenv)

    new_errors => scf_diis_buffer%error(ib)%blacs_matrix
    CALL copy_sparse_to_blacs_matrix(s(1)%matrix,new_errors,globenv)
    CALL blacs_symm("L","U",nao,nmo,2.0_wp,new_errors,c,0.0_wp,sc,globenv)

    CALL blacs_gemm("N","T",nao,nao,nmo,1.0_wp,sc,kc, 0.0_wp,new_errors,&
                    globenv)
    CALL blacs_gemm("N","T",nao,nao,nmo,1.0_wp,kc,sc,-1.0_wp,new_errors,&
                    globenv)

!   *** Get maximum error ***

    CALL blacs_maxval(new_errors,error_max,globenv)

    IF (ionode.AND.globenv%print%diis_information) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,E12.3)")&
       "Maximum SCF DIIS error vector element:",error_max
    END IF

!   *** Check, if a DIIS step is appropiate ***

    diis_step = ((scf_diis_buffer%ncall > 1).AND.(delta < eps_diis))

!   *** Update the SCF DIIS buffer ***

    IF (error_max < eps_diis) THEN

      b => scf_diis_buffer%b_matrix

      DO jb=1,nb
        old_errors => scf_diis_buffer%error(jb)%blacs_matrix
        CALL blacs_trace(old_errors,new_errors,b(jb,ib),globenv)
        b(ib,jb) = b(jb,ib)
      END DO

      IF (ionode.AND.globenv%print%diis_information) THEN
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
      IF (istat /= 0) CALL stop_memory(routine,"a",nb1*nb1*wp_size)
      ALLOCATE (b(nb1,nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"b",nb1*nb1*wp_size)
      ALLOCATE (ev(nb1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ev",nb1*wp_size)

!     *** Set up the linear DIIS equation system ***

      b(1:nb,1:nb) = scf_diis_buffer%b_matrix(1:nb,1:nb)

      b(1:nb,nb1) = -1.0_wp
      b(nb1,1:nb) = -1.0_wp
      b(nb1,nb1) = 0.0_wp

!     *** Solve the linear DIIS equation system ***

      CALL diagonalize_matrix(b(1:nb1,1:nb1),ev(1:nb1))

      a(1:nb1,1:nb1) = b(1:nb1,1:nb1)

      DO jb=1,nb1
        IF (ABS(ev(jb)) < 1.0E-12_wp) THEN
          a(1:nb1,jb) = 0.0_wp
        ELSE
          a(1:nb1,jb) = a(1:nb1,jb)/ev(jb)
        END IF
      END DO

      ev(1:nb) = MATMUL(a(1:nb,1:nb1),b(nb1,1:nb1))

!     *** Update Kohn-Sham matrix ***

      CALL blacs_set_all(kc,0.0_wp,globenv)

      DO jb=1,nb
        parameters => scf_diis_buffer%parameter(jb)%blacs_matrix
        CALL blacs_add(1.0_wp,kc,-ev(jb),parameters,globenv)
      END DO

      DEALLOCATE (a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"a")
      DEALLOCATE (b,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"b")
      DEALLOCATE (ev,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ev")

    ELSE

      CALL copy_blacs_to_blacs_matrix(parameters,kc)

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf_diis

! *****************************************************************************

END MODULE qs_diis
