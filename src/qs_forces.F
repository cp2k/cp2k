!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_forces [1.0] *
!!
!!   NAME
!!     qs_forces
!!
!!   FUNCTION
!!     Quickstep force driver routine
!!
!!   AUTHOR
!!     MK (12.06.2002)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_forces

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE dft_types,            ONLY: dft_control_type
  USE global_types,         ONLY: global_environment_type
  USE qs_blacs,             ONLY: allocate_blacs_matrix,&
                                  blacs_matrix_type,&
                                  copy_blacs_to_sparse_matrix,&
                                  deallocate_blacs_matrix,&
                                  get_blacs_matrix_info
  USE qs_build_ks_matrix,   ONLY: build_kohn_sham_matrix,qs_ks_env_type
  USE qs_core_energies,     ONLY: calculate_ecore_overlap
  USE qs_core_hamiltonian,  ONLY: build_core_hamiltonian_matrix
  USE qs_environment_types, ONLY: get_qs_env,&
                                  qs_environment_type
  USE qs_force_types,       ONLY: qs_force_type
  USE qs_mo_types,          ONLY: allocate_mo_set,&
                                  calculate_w_matrix,&
                                  get_mo_set,&
                                  mo_set_p_type,&
                                  mo_set_type
  USE qs_overlap,           ONLY: write_sparse_matrix
  USE sparse_matrix_types,  ONLY: real_matrix_p_type,&
                                  real_matrix_type,&
                                  replicate_matrix_structure,&
                                  symmetrise_diagonal_blocks
  USE termination,          ONLY: stop_memory,&
                                  stop_program
  USE timings,              ONLY: timeset,&
                                  timestop

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: update_forces,&
            write_forces

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE update_forces(ks_env,qs_env,globenv)

!   Purpose: Update the Quickstep forces.

!   History: - Creation (13.06.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type)                 :: qs_env
    TYPE(qs_ks_env_type)                      :: ks_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE update_forces (MODULE qs_forces)"

!   *** Local variables ***

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors,w_blacs
    TYPE(dft_control_type), POINTER  :: dft_control
    TYPE(mo_set_type), POINTER       :: mo_set
    TYPE(real_matrix_type), POINTER  :: w

    INTEGER :: handle,iatom,ikind,istat,nao,ncol_block,nrow_block,output_unit

    TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: c
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: k,p,p_old,s

!   ---------------------------------------------------------------------------

    CALL timeset("update_forces","I","",handle)
   
    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    s=s,&
                    c=c,&
                    k=k,&
                    p=p,&
                    p_old=p_old,&
                    force=force)

    IF (.NOT.dft_control%forces) RETURN

!   *** Force component from the core charge distribution overlap ***

    CALL calculate_ecore_overlap(qs_env=qs_env,&
                                 globenv=globenv,&
                                 calculate_forces=.TRUE.)

!   *** Load MO set (QS wavefunction) ***

    mo_set => c(1)%mo_set

    CALL get_mo_set(mo_set=mo_set,&
                    eigenvectors=mo_eigenvectors)

    CALL get_blacs_matrix_info(matrix=mo_eigenvectors,&
                               nrow_global=nao,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)

    CALL allocate_blacs_matrix(new_matrix=w_blacs,&
                               nrow_global=nao,&
                               ncol_global=nao,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="W MATRIX",&
                               globenv=globenv)

    NULLIFY (w)
    CALL replicate_matrix_structure(s(1)%matrix,w,"W MATRIX")
    CALL calculate_w_matrix(mo_set,w_blacs,globenv)
    CALL copy_blacs_to_sparse_matrix(w_blacs,w,globenv)
    CALL symmetrise_diagonal_blocks(w)
    IF (globenv%ionode.AND.globenv%print%w_matrix) THEN
      CALL write_sparse_matrix(w,4,6,qs_env,globenv)
    END IF

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.TRUE.,&
                                       p=p(1)%matrix,&
                                       w=w)

    CALL build_kohn_sham_matrix(ks_env=ks_env,&
                                qs_env=qs_env,&
                                globenv=globenv,&
                                ks=k,&
                                p=p,&
                                p_old=p_old,&
                                rebuild=.TRUE.,&
                                nrebuild=1,&
                                calculate_forces=.TRUE.)

    CALL write_forces(force,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE update_forces

! *****************************************************************************

  SUBROUTINE write_forces(qs_force,globenv)

!   Purpose: Write a Quickstep force data structure to 

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN)  :: globenv
    TYPE(qs_force_type), DIMENSION(:), POINTER :: qs_force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_forces (MODULE qs_force_types)"

!   *** Local variables ***

    INTEGER :: iatom,ikind,istat,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (.NOT.globenv%print%forces) RETURN

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routine,"The qs_force pointer is not associated "//&
                                "and cannot be printed")
    END IF

    output_unit = globenv%scr

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") "FORCES"

    DO ikind=1,SIZE(qs_force)
      WRITE (UNIT=output_unit,FMT="(/,T2,A,I5,/)") "ATOMIC KIND",ikind
      IF (ASSOCIATED(qs_force(ikind)%core_overlap)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%core_overlap,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "core_overlap:",iatom,qs_force(ikind)%core_overlap(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%gth_ppl)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%gth_ppl,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "gth_ppl:",iatom,qs_force(ikind)%gth_ppl(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%gth_ppnl)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%gth_ppnl,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "gth_ppnl:",iatom,qs_force(ikind)%gth_ppnl(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%kinetic)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%kinetic,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "kinetic:",iatom,qs_force(ikind)%kinetic(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%overlap)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%overlap,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "overlap:",iatom,qs_force(ikind)%overlap(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%rho_core)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%rho_core,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "rho_core:",iatom,qs_force(ikind)%rho_core(1:3,iatom)
        END DO
      END IF
      IF (ASSOCIATED(qs_force(ikind)%rho_elec)) THEN
        DO iatom=1,SIZE(qs_force(ikind)%rho_elec,2)
          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
            "rho_elec:",iatom,qs_force(ikind)%rho_elec(1:3,iatom)
        END DO
      END IF
!      IF (ASSOCIATED(qs_force(ikind)%total)) THEN
!        DO iatom=1,SIZE(qs_force(ikind)%total,2)
!          WRITE (UNIT=output_unit,FMT="(T2,A,T30,I6,3F12.6)")&
!            "total:",iatom,qs_force(ikind)%rho_total(1:3,iatom)
!        END DO
!      END IF
    END DO

  END SUBROUTINE write_forces

! *****************************************************************************

END MODULE qs_forces
