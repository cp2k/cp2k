!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/diis [1.0] *
!!
!!   NAME
!!     diis
!!
!!   FUNCTION
!!     Perform a direct inversion in the iterative subspace (DIIS).
!!
!!   AUTHOR
!!     Matthias Krack (28.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE diis

  USE kinds, ONLY: wp => dp

  USE core_hamiltonian,      ONLY: s
  USE global_types,          ONLY: global_environment_type
  USE mathlib,               ONLY: diagonalize_matrix,&
                                   symmetrize_matrix
  USE matrix_types,          ONLY: allocate_matrix,&
                                   copy_matrix,&
                                   get_block_node,&
                                   get_matrix_info,&
                                   real_matrix_set_type
  USE memory_utilities,      ONLY: reallocate
  USE mo_types,              ONLY: mo_set_type
  USE timings,               ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

  TYPE diis_buffer_type
    PRIVATE
    INTEGER                                           :: nbuffer,ncall
    TYPE(real_matrix_set_type), DIMENSION(:), POINTER :: error,parameter
    REAL(wp), DIMENSION(:,:), POINTER                 :: b_matrix
  END TYPE diis_buffer_type

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

  SUBROUTINE allocate_scf_diis_buffer(nbuffer,nao)

!   Purpose: Allocate and initialize a DIIS buffer for "nao*nao" parameter
!            variables and with a buffer size of "nbuffer".

!   History: - Creation (07.05.2001, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: nbuffer,nao

!   *** Local variables ***

    INTEGER :: ibuffer

!   ---------------------------------------------------------------------------

    scf_diis_buffer%nbuffer = nbuffer
    scf_diis_buffer%ncall = 0

    ALLOCATE (scf_diis_buffer%error(nbuffer))
    ALLOCATE (scf_diis_buffer%parameter(nbuffer))

    DO ibuffer=1,nbuffer
      NULLIFY (scf_diis_buffer%error(ibuffer)%matrix)
      CALL allocate_matrix(matrix=scf_diis_buffer%error(ibuffer)%matrix,&
                           nrow=nao,&
                           ncol=nao,&
                           matrix_name="SCF DIIS ERROR MATRIX",&
                           matrix_symmetry="symmetric")
      NULLIFY (scf_diis_buffer%parameter(ibuffer)%matrix)
      CALL allocate_matrix(matrix=scf_diis_buffer%parameter(ibuffer)%matrix,&
                           nrow=nao,&
                           ncol=nao,&
                           matrix_name="SCF DIIS PARAMETER MATRIX",&
                           matrix_symmetry="symmetric")
    END DO

    NULLIFY (scf_diis_buffer%b_matrix)
    scf_diis_buffer%b_matrix => reallocate(scf_diis_buffer%b_matrix,&
                                           1,nbuffer+1,&
                                           1,nbuffer+1)

    init_scf_diis_buffer_done = .TRUE.

  END SUBROUTINE allocate_scf_diis_buffer

! *****************************************************************************

  SUBROUTINE scf_diis(globenv,c,kc,sc,delta,error_max,diis_step)

!   Purpose: Update the SCF DIIS buffer.

!   History: - Creation (07.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(IN)             :: c
    TYPE(real_matrix_set_type), INTENT(INOUT) :: kc,sc
    REAL(wp), INTENT(IN)                      :: delta
    REAL(wp), INTENT(OUT)                     :: error_max
    LOGICAL, INTENT(OUT)                      :: diis_step

!   *** Local variables ***

    INTEGER :: handle,ib,jb,nao,nb,nb1,nmo,nparameter,output_unit

    REAL(wp), DIMENSION(:), POINTER   :: ev
    REAL(wp), DIMENSION(:,:), POINTER :: a,b,c_matrix,kc_matrix,&
                                         new_error_matrix,old_error_matrix,&
                                         parameter_matrix,sc_matrix

!   ---------------------------------------------------------------------------

    CALL timeset("scf_diis","I","",handle)

    diis_step = .FALSE.

!   *** Quick return, if no DIIS is requested ***

    IF (max_diis < 1) RETURN

    output_unit= globenv%scr

    CALL get_matrix_info(matrix=c%eigenvectors%matrix,nrow=nao)

    nmo = c%homo
    nparameter = nao*nao

    IF (.NOT.init_scf_diis_buffer_done) THEN
      CALL allocate_scf_diis_buffer(max_diis,nao)
      IF (globenv%print%diis_information) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "The SCF DIIS buffer was allocated and initialised"
      END IF
    END IF

    ib = MODULO(scf_diis_buffer%ncall,scf_diis_buffer%nbuffer) + 1
    scf_diis_buffer%ncall = scf_diis_buffer%ncall + 1
    nb = MIN(scf_diis_buffer%ncall,scf_diis_buffer%nbuffer)

    CALL get_block_node(matrix=scf_diis_buffer%parameter(ib)%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=parameter_matrix)

    CALL get_block_node(matrix=scf_diis_buffer%error(ib)%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=new_error_matrix)

    CALL get_block_node(matrix=c%eigenvectors%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=c_matrix)

    CALL get_block_node(matrix=kc%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=kc_matrix)

    CALL get_block_node(matrix=sc%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=sc_matrix)

    parameter_matrix(:,:) = kc_matrix(:,:)

    new_error_matrix(:,:) = 0.0_wp
    CALL copy_matrix(s%matrix,scf_diis_buffer%error(ib)%matrix)
    CALL symmetrize_matrix(new_error_matrix,"upper_to_lower")

    CALL dgemm("T","N",nao,nmo,nao,2.0_wp,parameter_matrix,nao,c_matrix,nao,&
               0.0_wp,kc_matrix,nao)
    CALL dgemm("T","N",nao,nmo,nao,2.0_wp,new_error_matrix,nao,c_matrix,nao,&
               0.0_wp,sc_matrix,nao)

    CALL dgemm("N","T",nao,nao,nmo,1.0_wp,sc_matrix,nao,kc_matrix,nao,0.0_wp,&
               new_error_matrix,nao)
    CALL dgemm("N","T",nao,nao,nmo,1.0_wp,kc_matrix,nao,sc_matrix,nao,-1.0_wp,&
               new_error_matrix,nao)

!   *** Get maximum error ***

    error_max = MAXVAL(ABS(new_error_matrix))

    IF (globenv%print%diis_information) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,E12.3)")&
       "Maximum SCF DIIS error vector element:",error_max
    END IF

!   *** Check, if a DIIS step is appropiate ***

    diis_step = ((scf_diis_buffer%ncall > 1).AND.(delta < eps_diis))

!   *** Update the SCF DIIS buffer ***

    IF (error_max < eps_diis) THEN

      DO jb=1,nb
        CALL get_block_node(matrix=scf_diis_buffer%error(jb)%matrix,&
                            block_row=1,&
                            block_col=1,&
                            block=old_error_matrix)
        CALL dgemm("T","N",1,1,nparameter,1.0_wp,old_error_matrix,nparameter,&
                   new_error_matrix,nparameter,0.0_wp,&
                   scf_diis_buffer%b_matrix(jb,ib),&
                   SIZE(scf_diis_buffer%b_matrix,1))
        scf_diis_buffer%b_matrix(ib,jb) = scf_diis_buffer%b_matrix(jb,ib)
      END DO

      IF (globenv%print%diis_information) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
          "The SCF DIIS buffer was updated"
      END IF

    ELSE

      diis_step = .FALSE.

    END IF

!   *** Perform DIIS update ***

    IF (diis_step) THEN

      nb1 = nb + 1

      NULLIFY (a,b,ev)

      a => reallocate(a,1,nb1,1,nb1)
      b => reallocate(b,1,nb1,1,nb1)
      ev => reallocate(ev,1,nb1)

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

      kc_matrix(:,:) = 0.0_wp

      DO jb=1,nb
        CALL get_block_node(matrix=scf_diis_buffer%parameter(jb)%matrix,&
                            block_row=1,&
                            block_col=1,&
                            block=parameter_matrix)
        kc_matrix(:,:) = kc_matrix(:,:) - ev(jb)*parameter_matrix(:,:)
      END DO

      DEALLOCATE (a,b,ev)

    ELSE

      kc_matrix(:,:) = parameter_matrix(:,:)

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf_diis

! *****************************************************************************

END MODULE diis
