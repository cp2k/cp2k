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

  USE kinds, ONLY: wp => dp,&
                   wp_size => dp_size

  USE qs_blacs,             ONLY: cp_fm_add,&
                                  cp_fm_gemm,&
                                  cp_full_matrix_p_type,&
                                  cp_full_matrix_type,&
                                  cp_fm_maxval,&
                                  cp_fm_set_all,&
                                  cp_fm_symm,&
                                  cp_fm_trace,&
                                  copy_blacs_to_blacs_matrix,&
                                  copy_real_matrix_to_full_matrix,&
                                  cp_fm_get_info,&
                                  replicate_blacs_matrix
  USE qs_environment_types, ONLY: get_qs_env,&
                                  qs_environment_type
  USE global_types,         ONLY: global_environment_type
  USE mathlib,              ONLY: diagonalize_matrix,&
                                  symmetrize_matrix
  USE qs_mo_types,          ONLY: get_mo_set,&
                                  mo_set_p_type, &
                                  mo_set_type
  USE sparse_matrix_types,  ONLY: real_matrix_p_type
  USE termination,          ONLY: stop_memory
  USE timings,              ONLY: timeset,&
                                  timestop
  USE scf_control_types,    ONLY: scf_control_type

  IMPLICIT NONE

  PRIVATE

  TYPE diis_buffer_type
    PRIVATE
    INTEGER                                          :: nbuffer,ncall
    TYPE(cp_full_matrix_p_type), DIMENSION(:,:), POINTER :: error,PARAMETER
    REAL(wp), DIMENSION(:,:), POINTER                :: b_matrix
  END TYPE diis_buffer_type

  TYPE(cp_full_matrix_type), POINTER :: new_errors,old_errors,parameters

  TYPE(diis_buffer_type) :: scf_diis_buffer
  LOGICAL                :: init_scf_diis_buffer_done = .FALSE.

! *** Public subroutines ***

  PUBLIC :: scf_diis

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_scf_diis_buffer(nbuffer,prototype_matrix,nspin)

!   Purpose: Allocate and initialize a DIIS buffer for "nao*nao" parameter
!            variables and with a buffer size of "nbuffer".

!   History: - Creation (07.05.2001, Matthias Krack)
!            - Changed to BLACS matrix usage (08.06.2001, MK)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER :: prototype_matrix
    INTEGER, INTENT(IN)              :: nbuffer,nspin

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_scf_diis_buffer (MODULE qs_diis)"

!   *** Local variables ***

    INTEGER :: ibuffer,istat,ispin

!   ---------------------------------------------------------------------------

    scf_diis_buffer%nbuffer = nbuffer
    scf_diis_buffer%ncall = 0

    ALLOCATE (scf_diis_buffer%error(nbuffer,nspin),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scf_diis_buffer%error",0)
    ALLOCATE (scf_diis_buffer%parameter(nbuffer,nspin),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scf_diis_buffer%parameter",0)

    DO ispin=1,nspin
    DO ibuffer=1,nbuffer
      CALL replicate_blacs_matrix(&
        prototype_matrix=prototype_matrix,&
        new_matrix=scf_diis_buffer%error(ibuffer,ispin)%matrix,&
        name="SCF DIIS ERROR MATRIX")
      CALL replicate_blacs_matrix(&
        prototype_matrix=prototype_matrix,&
        new_matrix=scf_diis_buffer%parameter(ibuffer,ispin)%matrix,&
        name="SCF DIIS PARAMETER MATRIX")
    END DO
    END DO

    ALLOCATE (scf_diis_buffer%b_matrix(nbuffer+1,nbuffer+1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"scf_diis_buffer%b_matrix",&
                                     (nbuffer + 1)**2*wp_size)

    init_scf_diis_buffer_done = .TRUE.

  END SUBROUTINE allocate_scf_diis_buffer

! *****************************************************************************

  SUBROUTINE scf_diis(mo_array,kc,sc,delta,error_max,diis_step,qs_env,globenv)

!   Purpose: Update the SCF DIIS buffer.

!   History: - Creation (07.05.2001, Matthias Krack)
!            - Changed to BLACS matrix usage (08.06.2001, MK)

!   ***************************************************************************

    TYPE(cp_full_matrix_type), POINTER                 :: sc
    TYPE(cp_full_matrix_p_type),DIMENSION(:), POINTER  :: kc
    TYPE(global_environment_type), INTENT(IN)          :: globenv
    TYPE(mo_set_type), POINTER                         :: mo_set
    TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
    TYPE(qs_environment_type), INTENT(IN)              :: qs_env
    REAL(wp), INTENT(IN)                               :: delta
    REAL(wp), INTENT(OUT)                              :: error_max
    LOGICAL, INTENT(OUT)                               :: diis_step

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE scf_diis (MODULE qs_diis)"

!   *** Local variables ***

    TYPE(cp_full_matrix_type), POINTER  :: c
    TYPE(scf_control_type), POINTER     :: scf_control

    INTEGER :: handle,ib,istat,jb,nao,nb,nb1,nmo,output_unit,&
         max_diis,nspin,ispin
    LOGICAL :: ionode
    REAL(wp) :: eps_diis,tmp

    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s

    REAL(wp), DIMENSION(:), ALLOCATABLE :: ev
    REAL(wp), DIMENSION(:,:), POINTER   :: a,b

!   ---------------------------------------------------------------------------

    CALL timeset("scf_diis","I","",handle)

    nspin=SIZE(mo_array)

    NULLIFY(c,scf_control,s,a,b)

    diis_step = .FALSE.

    ionode = globenv%ionode
    output_unit= globenv%scr

    CALL get_qs_env(qs_env=qs_env,s=s,scf_control=scf_control)

    max_diis=scf_control%max_diis
    eps_diis=scf_control%eps_diis

!   *** Quick return, if no DIIS is requested ***
    IF (max_diis < 1) THEN
       CALL timestop(0.0_wp,handle)
       RETURN
    END IF

    IF (.NOT.init_scf_diis_buffer_done) THEN
      CALL allocate_scf_diis_buffer(max_diis,kc(1)%matrix,nspin)
      IF (ionode.AND.globenv%print%diis_information) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "The SCF DIIS buffer was allocated and initialised"
      END IF
    END IF

    error_max=0.0_wp

    ib = MODULO(scf_diis_buffer%ncall,scf_diis_buffer%nbuffer) + 1
    scf_diis_buffer%ncall = scf_diis_buffer%ncall + 1
    nb = MIN(scf_diis_buffer%ncall,scf_diis_buffer%nbuffer)

    DO ispin=1,nspin

       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                       nao=nao,&
                       homo=nmo,&
                       eigenvectors=c)


       parameters => scf_diis_buffer%parameter(ib,ispin)%matrix
       CALL copy_blacs_to_blacs_matrix(kc(ispin)%matrix,parameters)
       CALL cp_fm_symm("L","U",nao,nmo,2.0_wp,parameters,c,0.0_wp,kc(ispin)%matrix)

       new_errors => scf_diis_buffer%error(ib,ispin)%matrix
       ! I guess that this copy can be avoided for LSD
       CALL copy_real_matrix_to_full_matrix(s(1)%matrix,new_errors)
       CALL cp_fm_symm("L","U",nao,nmo,2.0_wp,new_errors,c,0.0_wp,sc)

       CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_wp,sc,kc(ispin)%matrix, 0.0_wp,new_errors)
       CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_wp,kc(ispin)%matrix,sc,-1.0_wp,new_errors)

       CALL cp_fm_maxval(new_errors,tmp)
       error_max=MAX(error_max,tmp)

    ENDDO

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
        b(jb,ib)=0.0_wp
        DO ispin=1,nspin
           old_errors => scf_diis_buffer%error(jb,ispin)%matrix
           new_errors => scf_diis_buffer%error(ib,ispin)%matrix
           CALL cp_fm_trace(old_errors,new_errors,tmp)
           b(jb,ib)=b(jb,ib)+tmp
        ENDDO
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

      DO ispin=1,nspin
         CALL cp_fm_set_all(kc(ispin)%matrix,0.0_wp)
         DO jb=1,nb
              parameters => scf_diis_buffer%parameter(jb,ispin)%matrix
              CALL cp_fm_add(1.0_wp,kc(ispin)%matrix,-ev(jb),parameters)
         END DO
      ENDDO

      DEALLOCATE (a,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"a")
      DEALLOCATE (b,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"b")
      DEALLOCATE (ev,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"ev")

    ELSE

      DO ispin=1,nspin
         parameters => scf_diis_buffer%parameter(ib,ispin)%matrix
         CALL copy_blacs_to_blacs_matrix(parameters,kc(ispin)%matrix)
      ENDDO

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf_diis

! *****************************************************************************

END MODULE qs_diis
