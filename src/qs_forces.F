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

  USE atomic_kind_types,    ONLY: atomic_kind_type,&
                                  get_atomic_kind_set
  USE dft_types,            ONLY: dft_control_type
  USE global_types,         ONLY: global_environment_type
  USE message_passing,      ONLY: mp_sum
  USE qs_blacs,             ONLY: allocate_blacs_matrix,&
                                  blacs_matrix_type,&
                                  copy_blacs_to_sparse_matrix,&
                                  deallocate_blacs_matrix,&
                                  get_blacs_matrix_info
  USE qs_build_ks_matrix,   ONLY: qs_ks_update_qs_env,&
                                  qs_ks_env_type,&
                                  qs_ks_did_change
  USE qs_core_energies,     ONLY: calculate_ecore_overlap
  USE qs_core_hamiltonian,  ONLY: build_core_hamiltonian_matrix
  USE qs_environment_types, ONLY: get_qs_env,&
                                  qs_environment_type
  USE qs_force_types,       ONLY: init_qs_force,&
                                  qs_force_type
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

    TYPE(dft_control_type), POINTER  :: dft_control
    TYPE(mo_set_type), POINTER       :: mo_set
    TYPE(real_matrix_type), POINTER  :: w

    INTEGER :: handle,ikind,istat,output_unit

    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
    TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: c
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: k,p,s

!   ---------------------------------------------------------------------------

    CALL timeset("update_forces","I","",handle)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    s=s,&
                    c=c,&
                    k=k,&
                    p_old=p,&
                    force=force)

    IF (.NOT.dft_control%forces) RETURN

    CALL init_qs_force(force)

!   *** Force component from the core charge distribution overlap ***

    CALL calculate_ecore_overlap(qs_env=qs_env,&
                                 globenv=globenv,&
                                 calculate_forces=.TRUE.)

!   *** Load MO set (QS wavefunction) ***

    mo_set => c(1)%mo_set

    NULLIFY (w)
    CALL replicate_matrix_structure(s(1)%matrix,w,"W MATRIX")
    ! k(1) is ks ?
    CALL calculate_w_matrix(mo_set,k(1)%matrix,w,globenv)
    IF (globenv%ionode.AND.globenv%print%w_matrix) THEN
      CALL write_sparse_matrix(w,4,6,qs_env,globenv)
    END IF

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                       globenv=globenv,&
                                       calculate_forces=.TRUE.,&
                                       p=p(1)%matrix,&
                                       w=w)
    CALL qs_ks_did_change(ks_env,qs_env,globenv,p_changed=.true.)
    CALL qs_ks_update_qs_env(ks_env=ks_env,&
                             qs_env=qs_env,&
                             global_env=globenv,&
                             rebuild=.TRUE.,&
                             calculate_forces=.TRUE.)

    DO ikind=1,SIZE(force)
      CALL mp_sum(force(ikind)%gth_ppl,globenv%group)
      CALL mp_sum(force(ikind)%gth_ppnl,globenv%group)
      CALL mp_sum(force(ikind)%kinetic,globenv%group)
      CALL mp_sum(force(ikind)%overlap,globenv%group)
      CALL mp_sum(force(ikind)%rho_elec,globenv%group)
      force(ikind)%total(:,:) = force(ikind)%total(:,:) +&
                                force(ikind)%core_overlap(:,:) +&
                                force(ikind)%gth_ppl(:,:) +&
                                force(ikind)%gth_ppnl(:,:) +&
                                force(ikind)%kinetic(:,:) +&
                                force(ikind)%overlap(:,:) +&
                                force(ikind)%rho_core(:,:) +&
                                force(ikind)%rho_elec(:,:)
    END DO

    CALL write_forces(force,atomic_kind_set,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE update_forces

! *****************************************************************************

  SUBROUTINE write_forces(qs_force,atomic_kind_set,globenv)

!   Purpose: Write a Quickstep force data structure to 

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(qs_force_type), DIMENSION(:), POINTER    :: qs_force

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_forces (MODULE qs_force_types)"

!   *** Local variables ***

    INTEGER :: i,iatom,ikind,istat,natom,output_unit

    INTEGER, DIMENSION(:), ALLOCATABLE :: atom_of_kind,kind_of
    REAL(wp) :: grand_total(3)

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (.NOT.globenv%print%forces) RETURN

    IF (.NOT.ASSOCIATED(qs_force)) THEN
      CALL stop_program(routine,"The qs_force pointer is not associated "//&
                                "and cannot be printed")
    END IF

    output_unit = globenv%scr

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom)
    ALLOCATE (atom_of_kind(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)
    ALLOCATE (kind_of(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of",natom*int_size)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)

    WRITE (UNIT=output_unit,&
           FMT="(/,/,T2,A,/,/,T3,A,T11,A,T23,A,T43,A1,2(13X,A1))")&
      "FORCES","Atom","Kind","Component","X","Y","Z"

    grand_total(:)=0.0_wp
    DO iatom=1,natom
      ikind = kind_of(iatom)
      i = atom_of_kind(iatom)
      WRITE (UNIT=output_unit,FMT="(/,(T2,I5,4X,I4,T18,A,T34,3F14.6))")&
        iatom,ikind,"       overlap",qs_force(ikind)%overlap(1:3,i),&
        iatom,ikind,"       kinetic",qs_force(ikind)%kinetic(1:3,i),&
        iatom,ikind,"       gth_ppl",qs_force(ikind)%gth_ppl(1:3,i),&
        iatom,ikind,"      gth_ppnl",qs_force(ikind)%gth_ppnl(1:3,i),&
        iatom,ikind,"  core_overlap",qs_force(ikind)%core_overlap(1:3,i),&
        iatom,ikind,"      rho_core",qs_force(ikind)%rho_core(1:3,i),&
        iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
        iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
        grand_total(1:3)=grand_total(1:3)+qs_force(ikind)%total(1:3,i)
    END DO
    WRITE(UNIT=output_unit,FMT="(/,(T2,A,T34,3F14.6))") "Sum of total", &
        grand_total(1:3)

  END SUBROUTINE write_forces

! *****************************************************************************

END MODULE qs_forces
