!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_force [1.0] *
!!
!!   NAME
!!     qs_force
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

MODULE qs_force

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_restraint
  USE particle_types,                  ONLY: particle_type
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: zero_qs_force,&
                                             qs_force_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env,&
                                             calculate_w_matrix
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: mo_set_p_type,&
                                             mo_set_type
  USE qs_overlap,                      ONLY: write_sparse_matrix
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_semi_empirical_hamiltonian,   ONLY: build_se_core_matrix
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure, &
                                             set_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_force"

! *** Public subroutines ***

  PUBLIC :: qs_forces

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_forces(qs_env,globenv)

!   Purpose: Calculate the Quickstep forces.

!   History: - Creation (29.10.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_forces", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                ispin, istat, natom, nspin
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, matrix_w
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: gapw, obasis
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: mo_derivs 
!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    NULLIFY (atomic_kind_set)
    NULLIFY (mos)
    NULLIFY (dft_control)
    NULLIFY (force)
    NULLIFY (matrix_ks)
    NULLIFY (ks_env)
    NULLIFY (particle_set)
    NULLIFY (matrix_s)
    NULLIFY (scf_control)
    NULLIFY (matrix_w)
    NULLIFY (rho)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    force=force,&
                    particle_set=particle_set,&
                    scf_control=scf_control)

    natom = SIZE(particle_set)

    ! zero out the forces
    DO iatom=1,natom
       particle_set(iatom)%f=0.0_dp
    END DO

    gapw = (dft_control%qs_control%method=="GAPW")
    IF(dft_control%qs_control%semi_empirical .AND. &
       dft_control%qs_control%se_control%orthogonal_basis) THEN
       obasis = .TRUE.
    ELSE
       obasis = .FALSE.
    END IF

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)

    ALLOCATE (kind_of(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of",natom*int_size)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)

    CALL zero_qs_force(force)

    CALL qs_energies(qs_env,globenv)

    CALL get_qs_env(qs_env=qs_env,mos=mos,mo_derivs=mo_derivs, &
                    matrix_ks=matrix_ks, ks_env=ks_env, matrix_s=matrix_s, rho=rho)


    IF ( .NOT. obasis ) THEN
      !   *** Build W matrix ***
      !   *** (from the orthonormality constraint lambda C^T S C)
      nspin = SIZE(mos)

      CALL allocate_matrix_set(matrix_w,nspin)

      DO ispin=1,nspin
        CALL replicate_matrix_structure(matrix_s(1)%matrix,&
             matrix_w(ispin)%matrix,"W MATRIX")
        IF (qs_env%dft_control%restricted .AND. ispin>1) THEN
          ! not very elegant, indeed ...
          CALL set_matrix(matrix_w(ispin)%matrix,0.0_dp)
        ELSE 
          mo_set => mos(ispin)%mo_set
          IF (scf_control%use_ot) THEN
            CALL calculate_w_matrix(mo_set,mo_derivs(ispin)%matrix,&
                 matrix_w(ispin)%matrix)
          ELSE
            CALL calculate_w_matrix(mo_set,matrix_w(ispin)%matrix)
          END IF
        ENDIF
        IF (globenv%print%w_matrix) THEN
          CALL write_sparse_matrix(matrix_w(ispin)%matrix,4,6,qs_env,globenv)
        END IF
      END DO
    END IF

!   *** from an eventual mulliken restraint
    IF (dft_control%qs_control%mulliken_restraint) THEN
        CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control, &
                                globenv%para_env,matrix_s(1)%matrix, rho%rho_ao,w_matrix=matrix_w)
    ENDIF
    CALL set_qs_env(qs_env=qs_env,matrix_w=matrix_w)

!   *** compute core forces (also overwrites matrix_w) ***
    IF (dft_control%qs_control%semi_empirical) THEN
      CALL build_se_core_matrix(qs_env=qs_env,globenv=globenv,&
                                calculate_forces=.TRUE.)
    ELSE
      CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                         globenv=globenv,&
                                         calculate_forces=.TRUE.)
    END IF

!   *** compute grid-based forces ***
    CALL qs_ks_update_qs_env(ks_env=ks_env,&
                             qs_env=qs_env,&
                             calculate_forces=.TRUE.)

!   *** distribute forces ***
    DO ikind=1,SIZE(force)
      CALL mp_sum(force(ikind)%overlap,globenv%group)
      CALL mp_sum(force(ikind)%kinetic,globenv%group)
      CALL mp_sum(force(ikind)%gth_ppl,globenv%group)
      CALL mp_sum(force(ikind)%gth_ppnl,globenv%group)
      CALL mp_sum(force(ikind)%all_potential,globenv%group)
      CALL mp_sum(force(ikind)%core_overlap,globenv%group)
      CALL mp_sum(force(ikind)%rho_core,globenv%group)
      CALL mp_sum(force(ikind)%rho_elec,globenv%group)
      CALL mp_sum(force(ikind)%vhxc_atom,globenv%group)
      CALL mp_sum(force(ikind)%g0s_Vh_elec,globenv%group)
      CALL mp_sum(force(ikind)%hartree_2c,globenv%group)
      CALL mp_sum(force(ikind)%hartree_3c,globenv%group)
      force(ikind)%total(:,:) = force(ikind)%total(:,:) +&
                                force(ikind)%core_overlap(:,:) +&
                                force(ikind)%gth_ppl(:,:) +&
                                force(ikind)%gth_ppnl(:,:) +&
                                force(ikind)%all_potential(:,:) +&
                                force(ikind)%kinetic(:,:) +&
                                force(ikind)%overlap(:,:) +&
                                force(ikind)%rho_core(:,:) +&
                                force(ikind)%rho_elec(:,:) +&
                                force(ikind)%vhxc_atom(:,:) +&
                                force(ikind)%g0s_Vh_elec(:,:) +&
                                force(ikind)%hartree_2c(:,:) +&
                                force(ikind)%hartree_3c(:,:)
    END DO

    DO iatom=1,natom
      ikind = kind_of(iatom)
      i = atom_of_kind(iatom)
      ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      ! the force is - dE/dR, what is called force is actually the gradient
      ! Things should have the right name
      ! The minus sign below is a hack
      ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      force(ikind)%other(1:3,i)=-particle_set(iatom)%f(1:3)
      force(ikind)%total(1:3,i)=force(ikind)%total(1:3,i)+force(ikind)%other(1:3,i)
      IF (particle_set(iatom)%is_fixed) THEN
        particle_set(iatom)%f = 0.0_dp
      ELSE
        particle_set(iatom)%f = -force(ikind)%total(1:3,i)
      END IF
    END DO

    IF(dft_control%qs_control%semi_empirical .AND. obasis) THEN
      CALL write_forces(force,atomic_kind_set,2,globenv)
    ELSEIF(gapw) THEN
      CALL write_forces(force,atomic_kind_set,1,globenv)
    ELSE
      CALL write_forces(force,atomic_kind_set,0,globenv)
    END IF

    CALL deallocate_matrix_set(matrix_w)
    CALL set_qs_env(qs_env=qs_env,matrix_w=matrix_w)

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    DEALLOCATE (kind_of,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of")

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_forces

! *****************************************************************************

  SUBROUTINE write_forces(qs_force,atomic_kind_set,ftype,globenv)

!   Purpose: Write a Quickstep force data structure to

!   History: - Creation (05.06.2002,MK)

!   ***************************************************************************

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(IN)                      :: ftype

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE write_forces (MODULE qs_force_types)"

    CHARACTER(LEN=20)                        :: fmtstr3
    CHARACTER(LEN=35)                        :: fmtstr2
    CHARACTER(LEN=48)                        :: fmtstr1
    INTEGER                                  :: i, iatom, ikind, istat, &
                                                natom, ndigits, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    REAL(KIND=dp), DIMENSION(3)              :: grand_total

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

!   *** Variable precision output of the forces ***

    ndigits = MAX(6 + globenv%print%precision,1)

    fmtstr1 = "(/,/,T2,A,/,/,T3,A,T11,A,T23,A,T40,A1,2(  X,A1))"
    WRITE (UNIT=fmtstr1(41:42),FMT="(I2)") ndigits + 5

    fmtstr2 = "(/,(T2,I5,4X,I4,T18,A,T34,3F  .  ))"
    WRITE (UNIT=fmtstr2(32:33),FMT="(I2)") ndigits
    WRITE (UNIT=fmtstr2(29:30),FMT="(I2)") ndigits + 6

    fmtstr3 = "(/,T3,A,T34,3F  .  )"
    WRITE (UNIT=fmtstr3(18:19),FMT="(I2)") ndigits
    WRITE (UNIT=fmtstr3(15:16),FMT="(I2)") ndigits + 6

    WRITE (UNIT=output_unit,FMT=fmtstr1)&
      "FORCES","Atom","Kind","Component","X","Y","Z"

    grand_total(:) = 0.0_dp

    SELECT CASE (ftype)
      CASE DEFAULT
        DO iatom=1,natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          WRITE (UNIT=output_unit,FMT=fmtstr2)&
            iatom,ikind,"       overlap",qs_force(ikind)%overlap(1:3,i),&
            iatom,ikind,"       kinetic",qs_force(ikind)%kinetic(1:3,i),&
            iatom,ikind,"       gth_ppl",qs_force(ikind)%gth_ppl(1:3,i),&
            iatom,ikind,"      gth_ppnl",qs_force(ikind)%gth_ppnl(1:3,i),&
            iatom,ikind,"  core_overlap",qs_force(ikind)%core_overlap(1:3,i),&
            iatom,ikind,"      rho_core",qs_force(ikind)%rho_core(1:3,i),&
            iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
            iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
            grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
        END DO
      CASE (1)
        DO iatom=1,natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          WRITE (UNIT=output_unit,FMT=fmtstr2)&
            iatom,ikind,"       overlap",qs_force(ikind)%overlap(1:3,i),&
            iatom,ikind,"       kinetic",qs_force(ikind)%kinetic(1:3,i),&
            iatom,ikind,"       gth_ppl",qs_force(ikind)%gth_ppl(1:3,i),&
            iatom,ikind,"      gth_ppnl",qs_force(ikind)%gth_ppnl(1:3,i),&
            iatom,ikind," all_potential",qs_force(ikind)%all_potential(1:3,i),&
            iatom,ikind,"  core_overlap",qs_force(ikind)%core_overlap(1:3,i),&
            iatom,ikind,"      rho_core",qs_force(ikind)%rho_core(1:3,i),&
            iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
            iatom,ikind,"      vhxc_atom",qs_force(ikind)%vhxc_atom(1:3,i),&
            iatom,ikind,"    g0s_Vh_elec",qs_force(ikind)%g0s_Vh_elec(1:3,i),& 
            iatom,ikind,"     hartree_2c",qs_force(ikind)%hartree_2c(1:3,i),& 
            iatom,ikind,"     hartree_3c",qs_force(ikind)%hartree_3c(1:3,i),&
            iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
            grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
        END DO
      CASE (2)
        DO iatom=1,natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          WRITE (UNIT=output_unit,FMT=fmtstr2)&
            iatom,ikind," all_potential",qs_force(ikind)%all_potential(1:3,i),&
            iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
            iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
            grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
        END DO
    END SELECT

    WRITE (UNIT=output_unit,FMT=fmtstr3) "Sum of total",grand_total(1:3)

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind")

    DEALLOCATE (kind_of,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"kind_of")

  END SUBROUTINE write_forces

! *****************************************************************************

END MODULE qs_force
