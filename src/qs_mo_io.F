!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition and initialisation of the mo data type.
!> \par History
!>      - adapted to the new QS environment data structure (02.04.2002,MK)
!>      - set_mo_occupation added (17.04.02,MK)
!>      - correct_mo_eigenvalues added (18.04.02,MK)
!>      - calculate_density_matrix moved from qs_scf to here (22.04.02,MK)
!>      - mo_set_p_type added (23.04.02,MK)
!>      - PRIVATE attribute set for TYPE mo_set_type (23.04.02,MK)
!>      - started conversion to LSD (1.2003, Joost VandeVondele)
!>      - Split of from qs_mo_types (07.2014, JGH)
!> \author Matthias Krack (09.05.2001,MK)
! *****************************************************************************
MODULE qs_mo_io

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_read_unformatted,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_write_unformatted
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE kpoint_types,                    ONLY: kpoint_type
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_irecv,&
                                             mp_isend,&
                                             mp_sum,&
                                             mp_sync,&
                                             mp_wait
  USE orbital_pointers,                ONLY: indco,&
                                             nco,&
                                             nso
  USE orbital_symbols,                 ONLY: cgf_symbol,&
                                             sgf_symbol
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_mo_occupation,                ONLY: set_mo_occupation
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             duplicate_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mo_io'

  PUBLIC :: wfn_restart_file_name,&
            read_mo_set,&
            write_mo_set,&
            write_rt_mos_to_restart,&
            read_rt_mos_from_restart

  INTERFACE read_mo_set
    MODULE PROCEDURE read_mo_set_from_restart
  END INTERFACE

  INTERFACE write_mo_set
    MODULE PROCEDURE write_mo_set_to_output_unit,write_mo_set_to_restart
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param mo_array ...
!> \param particle_set ...
!> \param dft_section ...
!> \param kpoints ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_mo_set_to_restart(mo_array,particle_set,dft_section,&
                                     kpoints,atomic_kind_set,qs_kind_set,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(kpoint_type), OPTIONAL, POINTER     :: kpoints
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_to_restart', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=30), DIMENSION(2), PARAMETER :: &
      keys = (/"SCF%PRINT%RESTART_HISTORY","SCF%PRINT%RESTART        "/)

    INTEGER                                  :: handle, ikey, ires, ispin
    LOGICAL                                  :: failure, write_kp
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(1),error=error),cp_p_file) .OR.  &
         BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(2),error=error),cp_p_file) ) THEN

       IF(PRESENT(kpoints)) THEN
          write_kp = (kpoints%nkp > 0)
       ELSE
          write_kp = .FALSE.
       END IF

       IF(mo_array(1)%mo_set%use_mo_coeff_b) THEN
          ! we are using the dbcsr mo_coeff
          ! we copy it to the fm for anycase
          DO ispin=1,SIZE(mo_array)
             IF(.not.ASSOCIATED(mo_array(ispin)%mo_set%mo_coeff_b)) THEN
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF
             CALL copy_dbcsr_to_fm(mo_array(ispin)%mo_set%mo_coeff_b,&
                  mo_array(ispin)%mo_set%mo_coeff,error=error)!fm->dbcsr
          ENDDO
       ENDIF

       DO ikey=1,SIZE(keys)

          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               dft_section,keys(ikey),error=error),cp_p_file)) THEN

             ires = cp_print_key_unit_nr(logger,dft_section,keys(ikey),&
                  extension=".wfn", file_status="REPLACE", file_action="WRITE",&
                  do_backup=.TRUE., file_form="UNFORMATTED", error=error)

             CALL write_mo_set_low(mo_array, particle_set=particle_set, &
                                   qs_kind_set=qs_kind_set,ires=ires, error=error)

             IF(write_kp) CALL write_mos_kp(kpoints, ires=ires, error=error)

             CALL cp_print_key_finished_output(ires,logger,dft_section,TRIM(keys(ikey)), error=error)
          END IF
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE write_mo_set_to_restart

! *****************************************************************************
!> \brief ...
!> \param kpoints ...
!> \param ires ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_mos_kp(kpoints,ires,error)

    TYPE(kpoint_type), OPTIONAL, POINTER     :: kpoints
    INTEGER, INTENT(IN)                      :: ires
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mos_kp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ic, ik, ikk, iogroup, &
                                                is, mosource, nc, ns, &
                                                request, tag
    LOGICAL                                  :: failure, is_source, moalloc
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_inter_kp
    TYPE(mo_set_p_type), DIMENSION(:, :), &
      POINTER                                :: mos
    TYPE(mo_set_type), POINTER               :: mo_set, mop

    CALL timeset(routineN,handle)
    failure = .FALSE.

    ! first we dump the kpoint info
    IF (ires > 0) THEN
       WRITE (ires) "KPOINT INFO"
       WRITE (ires) kpoints%nkp
       DO ik=1,kpoints%nkp
          WRITE (ires) kpoints%wkp(ik),kpoints%xkp(1:3,ik)
       END DO
    END IF

    ! para groups
    para_env => kpoints%para_env
    para_env_inter_kp => kpoints%para_env_inter_kp
    iogroup = 0
    IF(kpoints%iogrp) THEN
       iogroup = para_env_inter_kp%mepos
    END IF
    CALL mp_sum(iogroup,para_env_inter_kp%group)

    moalloc = .FALSE.
    ! as each proc owns at least one kpoint and all kpoint mo structures 
    ! are the same, we can get the array dimensions this way
    mos => kpoints%kp_env(1)%kpoint_env%mos
    ns = SIZE(mos,2)
    nc = SIZE(mos,1)
    ! now come the mos at each kpoint
    DO ik=1,kpoints%nkp
       IF (ires > 0) WRITE (ires) ik,kpoints%wkp(ik),kpoints%xkp(1:3,ik)
       is_source = .FALSE.
       IF (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2)) THEN
          ! this is our kpoint
          ikk = ik - kpoints%kp_range(1) + 1
          is_source = .TRUE.
          mos => kpoints%kp_env(ikk)%kpoint_env%mos
       END IF
       IF(is_source .AND. kpoints%iogrp) THEN
          ! this is a 'local' MO, print right away
          DO is=1,SIZE(mos,2)
             DO ic=1,SIZE(mos,1)
                IF (ires > 0) WRITE (ires) ic,is
                mo_set => mos(ic,is)%mo_set
                CALL write_mo_set_basic(mo_set, ires, error)
             END DO
          END DO
          ! this is for synchronization
          mosource = 0
          CALL mp_sum(mosource,para_env_inter_kp%group)
       ELSE IF(is_source) THEN
          mosource = para_env_inter_kp%mepos
          CALL mp_sum(mosource,para_env_inter_kp%group)
          tag = ik
          DO is=1,SIZE(mos,2)
             DO ic=1,SIZE(mos,1)
                request = ic*100 + is
                mo_set => mos(ic,is)%mo_set
                CALL mp_isend(mo_set%mo_coeff%local_data,iogroup,para_env_inter_kp%group,request,tag)
                CALL mp_wait(request)
             END DO
          END DO
       ELSE IF(kpoints%iogrp) THEN
          IF(.NOT.moalloc) THEN
             NULLIFY(mop)
             mos => kpoints%kp_env(1)%kpoint_env%mos
             mo_set => mos(1,1)%mo_set
             CALL duplicate_mo_set(mop,mo_set,error)
             moalloc = .TRUE.
          END IF
          mosource = 0
          CALL mp_sum(mosource,para_env_inter_kp%group)
          tag = ik
          DO is=1,SIZE(mos,2)
             DO ic=1,SIZE(mos,1)
                request = ic*100 + is
                CALL mp_irecv(mop%mo_coeff%local_data,mosource,para_env_inter_kp%group,request,tag)
                CALL mp_wait(request)
                IF (ires > 0) WRITE (ires) ic,is
                CALL write_mo_set_basic(mop, ires, error)
             END DO
          END DO
       ELSE
          mosource = 0
          CALL mp_sum(mosource,para_env_inter_kp%group)
       END IF
       CALL mp_sync(para_env_inter_kp%group)
    END DO

    IF(moalloc) THEN
       CALL deallocate_mo_set(mop,error)
       moalloc = .FALSE.
    END IF
    
    CALL timestop(handle)

  END SUBROUTINE write_mos_kp

! *****************************************************************************
!> \brief ...
!> \param mo_array ...
!> \param rt_mos ...
!> \param particle_set ...
!> \param dft_section ...
!> \param qs_kind_set ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_rt_mos_to_restart(mo_array,rt_mos,particle_set,dft_section,&
                                     qs_kind_set,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rt_mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_rt_mos_to_restart', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=43), DIMENSION(2), PARAMETER :: keys = (/&
      "REAL_TIME_PROPAGATION%PRINT%RESTART_HISTORY",&
      "REAL_TIME_PROPAGATION%PRINT%RESTART        "/)

    INTEGER                                  :: handle, ikey, ires
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(1),error=error),cp_p_file) .OR.  &
         BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(2),error=error),cp_p_file) ) THEN

       DO ikey=1,SIZE(keys)

          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               dft_section,keys(ikey),error=error),cp_p_file)) THEN

             ires = cp_print_key_unit_nr(logger,dft_section,keys(ikey),&
                  extension=".rtpwfn", file_status="REPLACE", file_action="WRITE",&
                  do_backup=.TRUE., file_form="UNFORMATTED", error=error)

             CALL write_mo_set_low(mo_array, rt_mos=rt_mos, qs_kind_set=qs_kind_set,&
                                   particle_set=particle_set, ires=ires, error=error)
             CALL cp_print_key_finished_output(ires,logger,dft_section,TRIM(keys(ikey)), error=error)
          END IF
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE write_rt_mos_to_restart

! *****************************************************************************
!> \brief ...
!> \param mo_array ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param ires ...
!> \param rt_mos ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_mo_set_low(mo_array, qs_kind_set, particle_set, ires, rt_mos, error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: ires
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rt_mos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, ikind, imat, iset, ishell, ispin, istat, lmax, &
      lshell, max_block, nao, natom, nmo, nset, nset_max, nshell_max, nspin
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    CALL timeset(routineN,handle)
    nspin = SIZE(mo_array)
    nao = mo_array(1)%mo_set%nao

    IF (ires>0) THEN
       !     *** create some info about the basis set first ***
       natom = SIZE(particle_set,1)
       nset_max = 0
       nshell_max = 0

       DO iatom=1,natom
          NULLIFY(orb_basis_set,dftb_parameter)
          CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind),&
               basis_set=orb_basis_set,dftb_parameter=dftb_parameter,error=error)
          IF (ASSOCIATED(orb_basis_set)) THEN
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  nset=nset,&
                  nshell=nshell,&
                  l=l)
             nset_max = MAX(nset_max,nset)
             DO iset=1,nset
                nshell_max = MAX(nshell_max,nshell(iset))
             END DO
          ELSEIF (ASSOCIATED(dftb_parameter)) THEN
             CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
             nset_max = MAX(nset_max,1)
             nshell_max = MAX(nshell_max,lmax+1)
          ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "Unknown basis type. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          END IF
       END DO

       ALLOCATE (nso_info(nshell_max,nset_max,natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nso_info(:,:,:) = 0

       ALLOCATE (nshell_info(nset_max,natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nshell_info(:,:) = 0

       ALLOCATE (nset_info(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nset_info(:) = 0

       DO iatom=1,natom
          NULLIFY(orb_basis_set,dftb_parameter)
          CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind),&
               basis_set=orb_basis_set,dftb_parameter=dftb_parameter,error=error)
          IF (ASSOCIATED(orb_basis_set)) THEN
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  nset=nset,&
                  nshell=nshell,&
                  l=l)
             nset_info(iatom) = nset
             DO iset=1,nset
                nshell_info(iset,iatom)=nshell(iset)
                DO ishell=1,nshell(iset)
                   lshell = l(ishell,iset)
                   nso_info(ishell,iset,iatom) = nso(lshell)
                END DO
             END DO
          ELSEIF (ASSOCIATED(dftb_parameter)) THEN
             CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
             nset_info(iatom) = 1
             nshell_info(1,iatom)=lmax+1
             DO ishell=1,lmax+1
                lshell = ishell-1
                nso_info(ishell,1,iatom) = nso(lshell)
             END DO
          ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Unknown basis set type. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          END IF
       END DO

       WRITE (ires) natom,nspin,nao,nset_max,nshell_max
       WRITE (ires) nset_info
       WRITE (ires) nshell_info
       WRITE (ires) nso_info

       DEALLOCATE (nset_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DEALLOCATE (nshell_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DEALLOCATE (nso_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! use the scalapack block size as a default for buffering columns
    CALL cp_fm_get_info(mo_array(1)%mo_set%mo_coeff,ncol_block=max_block,error=error)

    DO ispin=1,nspin
       nmo=mo_array(ispin)%mo_set%nmo
       IF ((ires>0).AND.(nmo > 0)) THEN
          WRITE (ires) nmo,&
               mo_array(ispin)%mo_set%homo,&
               mo_array(ispin)%mo_set%lfomo,&
               mo_array(ispin)%mo_set%nelectron
          WRITE (ires) mo_array(ispin)%mo_set%eigenvalues(1:nmo),&
               mo_array(ispin)%mo_set%occupation_numbers(1:nmo)
       END IF
       IF(PRESENT(rt_mos))THEN
          DO imat=2*ispin-1,2*ispin
             CALL cp_fm_write_unformatted(rt_mos(imat)%matrix,ires,error)
          END DO
       ELSE
          CALL cp_fm_write_unformatted(mo_array(ispin)%mo_set%mo_coeff,ires,error)
       END IF
    END DO

    CALL timestop(handle)

  END SUBROUTINE write_mo_set_low
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param mo_set ...
!> \param ires ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE write_mo_set_basic(mo_set, ires, error)

    TYPE(mo_set_type), POINTER               :: mo_set
    INTEGER                                  :: ires
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_basic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, nmo

    CALL timeset(routineN,handle)

    nao = mo_set%nao
    nmo=mo_set%nmo

    IF ((ires>0).AND.(nmo > 0)) THEN
       WRITE (ires) nao,nmo,mo_set%homo,mo_set%lfomo,mo_set%nelectron
       WRITE (ires) mo_set%eigenvalues(1:nmo),mo_set%occupation_numbers(1:nmo)
    END IF
    CALL cp_fm_write_unformatted(mo_set%mo_coeff,ires,error)

    CALL timestop(handle)

  END SUBROUTINE write_mo_set_basic

! *****************************************************************************
!> \brief ...
!> \param filename ...
!> \param exist ...
!> \param section ...
!> \param logger ...
!> \param xas ...
!> \param rtp ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE wfn_restart_file_name(filename,exist,section,logger,xas,rtp,error)
    CHARACTER(LEN=default_path_length), &
      INTENT(OUT)                            :: filename
    LOGICAL, INTENT(OUT)                     :: exist
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL, INTENT(IN), OPTIONAL            :: xas, rtp
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n_rep_val
    LOGICAL                                  :: my_rtp, my_xas
    TYPE(section_vals_type), POINTER         :: print_key

    my_xas = .FALSE.
    my_rtp = .FALSE.
    IF(PRESENT(xas)) my_xas = xas
    IF(PRESENT(rtp)) my_rtp = rtp

    exist = .FALSE.
    CALL section_vals_val_get(section,"WFN_RESTART_FILE_NAME",n_rep_val=n_rep_val,error=error)
    IF (n_rep_val>0) THEN
      CALL section_vals_val_get(section,"WFN_RESTART_FILE_NAME",c_val=filename,error=error)
    ELSE
      IF(my_xas) THEN
       ! try to read from the filename that is generated automatically from the printkey
        print_key => section_vals_get_subs_vals(section,"PRINT%RESTART",error=error)
        filename = cp_print_key_generate_filename(logger,print_key, &
                    extension="",my_local=.FALSE., error=error)
      ELSE IF (my_rtp)THEN
       ! try to read from the filename that is generated automatically from the printkey
        print_key => section_vals_get_subs_vals(section,"REAL_TIME_PROPAGATION%PRINT%RESTART",error=error)
        filename = cp_print_key_generate_filename(logger,print_key, &
                    extension=".rtpwfn",my_local=.FALSE., error=error)
      ELSE
        ! try to read from the filename that is generated automatically from the printkey
        print_key => section_vals_get_subs_vals(section,"SCF%PRINT%RESTART",error=error)
        filename = cp_print_key_generate_filename(logger,print_key, &
                    extension=".wfn", my_local=.FALSE., error=error)
      END IF
    ENDIF
    IF(.NOT.my_xas) THEN
      INQUIRE(FILE=filename,exist=exist)
    END IF

  END SUBROUTINE wfn_restart_file_name

! *****************************************************************************
!> \brief ...
!> \param mo_array ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param para_env ...
!> \param id_nr ...
!> \param multiplicity ...
!> \param dft_section ...
!> \param natom_mismatch ...
!> \param kpoints ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_mo_set_from_restart(mo_array,atomic_kind_set,qs_kind_set,particle_set,&
       para_env,id_nr,multiplicity,dft_section,natom_mismatch,kpoints,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: id_nr, multiplicity
    TYPE(section_vals_type), POINTER         :: dft_section
    LOGICAL, INTENT(OUT), OPTIONAL           :: natom_mismatch
    TYPE(kpoint_type), OPTIONAL, POINTER     :: kpoints
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mo_set_from_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: group, handle, ispin, natom, &
                                                nspin, restart_unit, source
    LOGICAL                                  :: exist, failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    failure = .FALSE.

    nspin = SIZE(mo_array)
    restart_unit = -1

    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN

      natom = SIZE(particle_set,1)
      CALL wfn_restart_file_name(file_name,exist,dft_section,logger,error=error)
      IF (id_nr/=0) THEN
         ! Is it one of the backup files?
         file_name = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

    END IF

    CALL read_mos_restart_low (mo_array,para_env=para_env, qs_kind_set=qs_kind_set,&
         particle_set=particle_set, natom=natom,&
         rst_unit=restart_unit, multiplicity=multiplicity, natom_mismatch=natom_mismatch, &
         error=error)
    IF (PRESENT(natom_mismatch)) THEN
       ! read_mos_restart_low only the io_node returns natom_mismatch, must broadcast it
       CALL mp_bcast(natom_mismatch,source,group)
       IF (natom_mismatch) THEN
          IF (para_env%ionode) CALL close_file(unit_number=restart_unit)
          CALL timestop(handle)
          RETURN
       ENDIF
    ENDIF

    IF(PRESENT(kpoints)) THEN
       IF(kpoints%nkp > 0) THEN
          CALL read_mos_kp(kpoints,restart_unit,error)
       END IF
    END IF

    ! Close restart file
    IF (para_env%ionode) CALL close_file(unit_number=restart_unit)

    DO ispin = 1,nspin
      CALL write_mo_set(mo_array(ispin)%mo_set,atomic_kind_set,qs_kind_set,&
                        particle_set,4,dft_section,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE read_mo_set_from_restart
! *****************************************************************************
!> \brief ...
!> \param mo_array ...
!> \param rt_mos ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param para_env ...
!> \param id_nr ...
!> \param multiplicity ...
!> \param dft_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_rt_mos_from_restart(mo_array,rt_mos,atomic_kind_set,qs_kind_set, &
                 particle_set,para_env,id_nr,multiplicity,dft_section,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rt_mos
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: id_nr, multiplicity
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_rt_mos_from_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: group, handle, ispin, natom, &
                                                nspin, restart_unit, source
    LOGICAL                                  :: exist, failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    failure = .FALSE.

    nspin = SIZE(mo_array)
    restart_unit = -1

    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN

      natom = SIZE(particle_set,1)
      CALL wfn_restart_file_name(file_name,exist,dft_section,logger,rtp=.TRUE.,error=error)
      IF (id_nr/=0) THEN
         ! Is it one of the backup files?
         file_name = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

    END IF

    CALL read_mos_restart_low (mo_array, rt_mos=rt_mos,para_env= para_env,&
         particle_set=particle_set, qs_kind_set=qs_kind_set, natom=natom,&
         rst_unit=restart_unit, multiplicity=multiplicity, error=error)

    ! Close restart file
    IF (para_env%ionode) CALL close_file(unit_number=restart_unit)

    DO ispin = 1,nspin
      CALL write_mo_set(mo_array(ispin)%mo_set,atomic_kind_set,qs_kind_set, &
                        particle_set,4,dft_section,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE read_rt_mos_from_restart
! *****************************************************************************
!> \brief ...
!> \param kpoints ...
!> \param ires ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_mos_kp(kpoints,ires,error)

    TYPE(kpoint_type), OPTIONAL, POINTER     :: kpoints
    INTEGER, INTENT(IN)                      :: ires
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mos_kp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ic, ik, ikk, is, jc, &
                                                jk, js, nc, nkp, ns
    LOGICAL                                  :: failure, is_source, moalloc
    REAL(KIND=dp)                            :: dd, wkp, xkp(3)
    TYPE(mo_set_p_type), DIMENSION(:, :), &
      POINTER                                :: mos
    TYPE(mo_set_type), POINTER               :: mo_set, mop

    CALL timeset(routineN,handle)
    failure = .FALSE.

    ! first we read the kpoint info
    IF (ires > 0) THEN
       READ (ires) 
       READ (ires) nkp
       CALL cp_assert(nkp==kpoints%nkp,cp_fatal_level,cp_assertion_failed,routineP,&
                      "Different number of kpoints",error,failure)
       DO ik=1,kpoints%nkp
          READ (ires) wkp,xkp(1:3)
          dd = ABS(wkp-kpoints%wkp(ik)) + SUM(ABS(xkp(1:3)-kpoints%xkp(1:3,ik)))
          CALL cp_assert(dd > 1.0e-12_dp,cp_fatal_level,cp_assertion_failed,routineP,&
                         "Different order/position of kpoints",error,failure)
       END DO
    END IF

    moalloc = .FALSE.
    ! as each proc owns at least one kpoint and all kpoint mo structures 
    ! are the same, we can get the array dimensions this way
    mos => kpoints%kp_env(1)%kpoint_env%mos
    ns = SIZE(mos,2)
    nc = SIZE(mos,1)
    ! now come the mos at each kpoint
    DO ik=1,kpoints%nkp
       IF (ires > 0) THEN
          READ (ires) jk,wkp,xkp(1:3)
          dd = ABS(ik-jk)
          dd = dd + ABS(wkp-kpoints%wkp(ik)) + SUM(ABS(xkp(1:3)-kpoints%xkp(1:3,ik)))
          CALL cp_assert(dd > 1.0e-12_dp,cp_fatal_level,cp_assertion_failed,routineP,&
                         "Different order/position of kpoints",error,failure)
       END IF
       is_source = .FALSE.
       IF (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2)) THEN
          ! this is our kpoint
          ikk = ik - kpoints%kp_range(1) + 1
          is_source = .TRUE.
          mos => kpoints%kp_env(ikk)%kpoint_env%mos
       END IF
       IF(is_source .AND. kpoints%iogrp) THEN
          ! this is a 'local' MO, print right away
          DO is=1,SIZE(mos,2)
             DO ic=1,SIZE(mos,1)
                IF (ires > 0) THEN
                   READ (ires) jc,js
                   CALL cp_assert((ic==jc .AND. is==js),cp_fatal_level,cp_assertion_failed,routineP,&
                         "Different order/position/type of kpoints",error,failure)
                END IF
                mo_set => mos(ic,is)%mo_set
                CALL read_mo_set_basic(mo_set, ires, error)
             END DO
          END DO
       ELSE IF(is_source) THEN
          !  not yet implemented
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ELSE IF(kpoints%iogrp) THEN
          IF(.NOT.moalloc) THEN
             NULLIFY(mop)
             mos => kpoints%kp_env(1)%kpoint_env%mos
             mo_set => mos(1,1)%mo_set
             CALL duplicate_mo_set(mop,mo_set,error)
             moalloc = .TRUE.
          END IF
          !  not yet implemented
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    IF(moalloc) THEN
       CALL deallocate_mo_set(mop,error)
       moalloc = .FALSE.
    END IF
    
    CALL timestop(handle)

  END SUBROUTINE read_mos_kp

! *****************************************************************************
!> \brief Reading the mos from apreviously defined restart file
!> \param mos ...
!> \param para_env ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param natom ...
!> \param rst_unit ...
!> \param multiplicity ...
!> \param rt_mos ...
!> \param natom_mismatch ...
!> \param error ...
!> \par History
!>      12.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE read_mos_restart_low (mos,  para_env, qs_kind_set, particle_set, natom, rst_unit, &
        multiplicity, rt_mos, natom_mismatch, error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: natom, rst_unit
    INTEGER, INTENT(in), OPTIONAL            :: multiplicity
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rt_mos
    LOGICAL, INTENT(OUT), OPTIONAL           :: natom_mismatch
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mos_restart_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: group, homo, homo_read, i, iatom, ikind, imat, irow, iset, &
      iset_read, ishell, ishell_read, iso, ispin, istat, lfomo_read, lmax, &
      lshell, my_mult, nao, nao_read, natom_read, nelectron, nelectron_read, &
      nmo, nmo_read, nnshell, nset, nset_max, nshell_max, nspin, nspin_read, &
      offset_read, source
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info, offset_info
    LOGICAL                                  :: failure, minbas, natom_match, &
                                                use_this
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig_read, occ_read
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer, vecbuffer_read
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    logger => cp_error_get_logger(error)

    nspin = SIZE(mos)
    nao = mos(1)%mo_set%nao
    my_mult = 0
    IF(PRESENT(multiplicity)) my_mult = multiplicity
    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN
       READ (rst_unit) natom_read,nspin_read,nao_read,nset_max,nshell_max
       IF(PRESENT(rt_mos)) THEN
          IF (nspin_read /= nspin) THEN
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                 "To change nspin is not possible. "//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
         END IF
      ELSE
         ! we should allow for restarting with different spin settings
         IF (nspin_read /= nspin) THEN
            WRITE(cp_logger_get_default_unit_nr(logger),*)  &
                 "READ RESTART : WARNING : nspin is not equal "
         END IF
         ! this case needs fixing of homo/lfomo/nelec/occupations ...
         IF (nspin_read > nspin) THEN
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                 "Reducing nspin is not possible. "//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
         ENDIF
      END IF

      natom_match = (natom_read == natom)

      IF (natom_match) THEN ! actually do the read read

         ! Let's make it possible to change the basis set
         ALLOCATE (nso_info(nshell_max,nset_max,natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE (nshell_info(nset_max,natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE (nset_info(natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE (offset_info(nshell_max,nset_max,natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

         IF (nao_read /= nao) THEN
            WRITE(cp_logger_get_default_unit_nr(logger),*) &
                 " READ RESTART : WARNING : DIFFERENT # AOs ",nao,nao_read
            IF(PRESENT(rt_mos))&
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                    "To change basis is not possible. "//&
                    CPSourceFileRef,&
                    only_ionode=.TRUE.)
         END IF

         READ (rst_unit) nset_info
         READ (rst_unit) nshell_info
         READ (rst_unit) nso_info

         i=1
         DO iatom=1,natom
            DO iset=1,nset_info(iatom)
               DO ishell=1,nshell_info(iset,iatom)
                  offset_info(ishell,iset,iatom) = i
                  i=i+nso_info(ishell,iset,iatom)
               END DO
            END DO
         END DO

         ALLOCATE(vecbuffer_read(1,nao_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      END IF ! natom_match
   END IF ! ionode

   ! make natom_match and natom_mismatch uniform across all nodes
   CALL mp_bcast(natom_match,source,group)
   IF (PRESENT(natom_mismatch)) natom_mismatch = .NOT. natom_match
   ! handle natom_match false
   IF (.NOT. natom_match) THEN
      IF (PRESENT(natom_mismatch)) THEN
         WRITE(cp_logger_get_default_unit_nr(logger),*) &
              " READ RESTART : WARNING : DIFFERENT natom, returning ",natom,natom_read
         RETURN
      ELSE
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
              "Incorrect number of atoms in restart file. "//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
      ENDIF
   ENDIF

   CALL mp_bcast(nspin_read,source,group)

   ALLOCATE (vecbuffer(1,nao),STAT=istat)
   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

   DO ispin=1,nspin

      nmo=mos(ispin)%mo_set%nmo
      homo=mos(ispin)%mo_set%homo
      mos(ispin)%mo_set%eigenvalues(:) = 0.0_dp
      mos(ispin)%mo_set%occupation_numbers(:) = 0.0_dp
      CALL cp_fm_set_all(mos(ispin)%mo_set%mo_coeff,0.0_dp,error=error)

      IF (para_env%ionode.AND.(nmo > 0)) THEN
         READ (rst_unit) nmo_read, homo_read, lfomo_read, nelectron_read
         ALLOCATE(eig_read(nmo_read), occ_read(nmo_read), STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         eig_read = 0.0_dp
         occ_read = 0.0_dp

         nmo = MIN(nmo,nmo_read)
         CALL cp_assert((nmo_read >= nmo),cp_warning_level,cp_assertion_failed,routineP,&
              "The number of MOs on the restart unit is smaller than the number of "//&
              "the allocated MOs. The MO set will be padded with zeros!"//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
        CALL cp_assert((nmo_read<=nmo),cp_warning_level,cp_assertion_failed,routineP,&
             "The number of MOs on the restart unit is greater than the number of "//&
             "the allocated MOs. The read MO set will be truncated!"//&
             CPSourceFileRef,&
             only_ionode=.TRUE.)

        READ (rst_unit) eig_read(1:nmo_read), occ_read(1:nmo_read)
        mos(ispin)%mo_set%eigenvalues(1:nmo) = eig_read(1:nmo)
        mos(ispin)%mo_set%occupation_numbers(1:nmo) = occ_read(1:nmo)
        DEALLOCATE(eig_read, occ_read, STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        mos(ispin)%mo_set%homo = homo_read
        mos(ispin)%mo_set%lfomo = lfomo_read
        IF (homo_read > nmo) THEN
          IF(nelectron_read==mos(ispin)%mo_set%nelectron) THEN
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                 "The number of occupied MOs on the restart unit is larger than "//&
                 "the allocated MOs. The read MO set will be truncated and the occupation numbers recalculated!"//&
CPSourceFileRef,&
                 only_ionode=.TRUE.)
             CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,error=error)
          ELSE
              ! can not make this a warning i.e. homo must be smaller than nmo
              ! otherwise e.g. set_mo_occupation will go out of bounds
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                "Number of occupied MOs on restart unit larger than allocated MOs. "//&
CPSourceFileRef,&
                 only_ionode=.TRUE.)
          END IF
        END IF
      END IF

      CALL mp_bcast(nmo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%homo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%lfomo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%nelectron,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%eigenvalues,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%occupation_numbers,source,group)
      IF(PRESENT(rt_mos))THEN
         DO imat=2*ispin-1,2*ispin
            DO i=1,nmo
               IF (para_env%ionode) THEN
                  READ (rst_unit) vecbuffer
               ELSE
                  vecbuffer(1,:) = 0.0_dp
               END IF
               CALL mp_bcast(vecbuffer,source,group)
               CALL cp_fm_set_submatrix(rt_mos(imat)%matrix,&
                    vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)
            END DO
         END DO
      ELSE
         DO i=1,nmo
            IF (para_env%ionode) THEN
               READ (rst_unit) vecbuffer_read
               ! now, try to assign the read to the real vector
               ! in case the basis set changed this involves some guessing
               irow=1
               DO iatom=1,natom
                  NULLIFY(orb_basis_set,dftb_parameter,l,nshell)
                  CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                  CALL get_qs_kind(qs_kind_set(ikind),&
                       basis_set=orb_basis_set,dftb_parameter=dftb_parameter,error=error)
                  IF (ASSOCIATED(orb_basis_set)) THEN
                     CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                          nset=nset,&
                          nshell=nshell,&
                          l=l)
                     minbas=.FALSE.
                  ELSEIF (ASSOCIATED(dftb_parameter)) THEN
                     CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                     nset = 1
                     minbas=.TRUE.
                  ELSE
                     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                          "Unknown basis set type. "//&
                          CPSourceFileRef,&
                          only_ionode=.TRUE.)
                  END IF

                  use_this = .TRUE.
                  iset_read = 1
                  DO iset=1,nset
                     ishell_read = 1
                     IF(minbas) THEN
                        nnshell = lmax+1
                     ELSE
                        nnshell = nshell(iset)
                     END IF
                     DO ishell=1,nnshell
                        IF(minbas) THEN
                           lshell = ishell-1
                        ELSE
                           lshell = l(ishell,iset)
                        END IF
                        IF (iset_read > nset_info(iatom)) use_this = .FALSE.
                        IF (use_this) THEN ! avoids out of bound access of the lower line if false
                           IF (nso(lshell) == nso_info(ishell_read,iset_read,iatom)) THEN
                              offset_read=offset_info(ishell_read,iset_read,iatom)
                              ishell_read=ishell_read+1
                              IF (ishell_read > nshell_info(iset,iatom)) THEN
                                 ishell_read = 1
                                 iset_read = iset_read+1
                              END IF
                           ELSE
                              use_this = .FALSE.
                           END IF
                        END IF
                        DO iso=1,nso(lshell)
                           IF (use_this) THEN
                              IF (offset_read-1+iso.LT.1 .OR. offset_read-1+iso.GT.nao_read) THEN
                                 vecbuffer(1,irow)=0.0_dp
                              ELSE
                                 vecbuffer(1,irow)=vecbuffer_read(1,offset_read-1+iso)
                              END IF
                           ELSE
                              vecbuffer(1,irow) = 0.0_dp
                           END IF
                           irow = irow + 1
                        END DO
                        use_this = .TRUE.
                     END DO
                  END DO
               END DO

            ELSE

               vecbuffer(1,:) = 0.0_dp

            END IF

            CALL mp_bcast(vecbuffer,source,group)
            CALL cp_fm_set_submatrix(mos(ispin)%mo_set%mo_coeff,&
                 vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)
         END DO
      END IF
      ! Skip extra MOs if there any
      IF (para_env%ionode) THEN
        !ignore nmo = 0
        IF(nmo>0) THEN
          DO i=nmo+1,nmo_read
            READ (rst_unit) vecbuffer_read
          END DO
        END IF
      END IF

      IF(.NOT.PRESENT(rt_mos)) THEN
        IF(ispin==1 .AND. nspin_read<nspin) THEN

          mos(ispin+1)%mo_set%homo = mos(ispin)%mo_set%homo
          mos(ispin+1)%mo_set%lfomo = mos(ispin)%mo_set%lfomo
          nelectron = mos(ispin)%mo_set%nelectron
          IF (my_mult.NE.1) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Restarting an LSD calculation from an LDA wfn only works for multiplicity=1 (singlets)."//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
          IF(mos(ispin+1)%mo_set%nelectron < 0)THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "LSD: too few electrons for this multiplisity. "//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
          mos(ispin+1)%mo_set%eigenvalues = mos(ispin)%mo_set%eigenvalues
          mos(ispin)%mo_set%occupation_numbers = mos(ispin)%mo_set%occupation_numbers/2.0_dp
          mos(ispin+1)%mo_set%occupation_numbers = mos(ispin)%mo_set%occupation_numbers
          CALL cp_fm_to_fm(mos(ispin)%mo_set%mo_coeff,mos(ispin+1)%mo_set%mo_coeff,error=error)
          EXIT
       END IF
    END IF
 END DO   ! ispin

    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (para_env%ionode) THEN
       DEALLOCATE(vecbuffer_read,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(offset_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(nso_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(nshell_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(nset_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE read_mos_restart_low

! *****************************************************************************
!> \brief ...
!> \param mo_set ...
!> \param ires ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE read_mo_set_basic (mo_set, ires, error)

    TYPE(mo_set_type), POINTER               :: mo_set
    INTEGER, INTENT(in)                      :: ires
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mo_set_basic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: mao, mmo, nao, nmo
    LOGICAL                                  :: failure

    mao = mo_set%nao
    mmo = mo_set%nmo

    IF ((ires>0).AND.(nmo > 0)) THEN
       READ (ires) nao,nmo,mo_set%homo,mo_set%lfomo,mo_set%nelectron
       CALL cp_assert(mao==nao,cp_fatal_level,cp_assertion_failed,routineP,&
                         "Different size of basis set",error,failure)
       CALL cp_assert(mmo==nmo,cp_fatal_level,cp_assertion_failed,routineP,&
                         "Different number of MOs",error,failure)
       READ (ires) mo_set%eigenvalues(1:nmo),mo_set%occupation_numbers(1:nmo)
    END IF
    CALL cp_fm_read_unformatted(mo_set%mo_coeff,ires,error)

  END SUBROUTINE read_mo_set_basic

! *****************************************************************************
!> \brief   Write the MO eigenvalues, MO occupation numbers and
!>          MO mo_coeff.
!> \param mo_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param before ...
!> \param dft_section ...
!> \param spin ...
!> \param last ...
!> \param error ...
!> \date    15.05.2001
!> \par History:
!>       - Optionally print Cartesian MOs (20.04.2005,MK)
!> \par Variables
!>       - after : Number of digits after point.
!>       - before: Number of digits before point.
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_mo_set_to_output_unit(mo_set,atomic_kind_set,qs_kind_set,particle_set,&
       before,dft_section,spin,last,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: before
    TYPE(section_vals_type), POINTER         :: dft_section
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: spin
    LOGICAL, INTENT(IN), OPTIONAL            :: last
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_to_output_unit', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: symbol
    CHARACTER(LEN=12), DIMENSION(:), POINTER :: bcgf_symbol
    CHARACTER(LEN=16)                        :: fmtstr5, fmtstr6
    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=2*default_string_length)   :: name
    CHARACTER(LEN=22)                        :: fmtstr2
    CHARACTER(LEN=25)                        :: fmtstr1, fmtstr7
    CHARACTER(LEN=27)                        :: fmtstr4
    CHARACTER(LEN=38)                        :: fmtstr3
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: bsgf_symbol
    INTEGER :: after, first_mo, from, iatom, icgf, ico, icol, ikind, imo, &
      irow, iset, isgf, ishell, iso, istat, iw, jcol, last_mo, left, lmax, &
      lshell, natom, ncgf, ncol, ncol_global, nrow_global, nset, nsgf, right, &
      scf_step, to, width
    INTEGER, DIMENSION(:), POINTER           :: mo_index_range, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    LOGICAL                                  :: failure, ionode, my_last, &
                                                p_cart, p_eval, p_evec, &
                                                p_occ, should_output
    REAL(KIND=dp)                            :: gap
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cmatrix, smatrix
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    NULLIFY (bcgf_symbol)
    NULLIFY (bsgf_symbol)
    NULLIFY (logger)
    NULLIFY (mo_index_range)
    NULLIFY (nshell)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    failure = .FALSE.
    CALL section_vals_val_get(dft_section,"PRINT%MO%EIGENVALUES",l_val=p_eval,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%EIGENVECTORS",l_val=p_evec,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%OCCUPATION_NUMBERS",l_val=p_occ,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%CARTESIAN",l_val=p_cart,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%MO_INDEX_RANGE",i_vals=mo_index_range,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%NDIGITS",i_val=after,error=error)
    after = MIN(MAX(after,1),16)
    should_output = BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
         "PRINT%MO",error=error),cp_p_file)

    IF ((.NOT.should_output).OR.(.NOT.(p_eval.OR.p_evec.OR.p_occ))) RETURN

    IF (PRESENT(last)) THEN
       my_last = last
    ELSE
       my_last = .FALSE.
    END IF

    scf_step = logger%iter_info%iteration(logger%iter_info%n_rlevel) - 1

    IF (p_evec) THEN
       CALL cp_fm_get_info(mo_set%mo_coeff,&
            nrow_global=nrow_global,&
            ncol_global=ncol_global,error=error)
       ALLOCATE(smatrix(nrow_global,ncol_global),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL cp_fm_get_submatrix(mo_set%mo_coeff,smatrix,error=error)
       IF (.NOT.ionode) THEN
          DEALLOCATE(smatrix,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

    iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%MO",&
         ignore_should_output=should_output,&
         extension=".MOLog",error=error)

    IF (iw > 0) THEN

       CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
       CALL get_qs_kind_set(qs_kind_set, ncgf=ncgf, nsgf=nsgf,error=error)

       ! Definition of the variable formats

       fmtstr1 = "(/,T2,21X,  (  X,I5,  X))"
       fmtstr2 = "(T2,21X,  (1X,F  .  ))"
       fmtstr3 = "(T2,I5,1X,I5,1X,A,1X,A6,  (1X,F  .  ))"

       width = before + after + 3
       ncol = INT(56/width)

       right = MAX((after-2),1)
       left =  width - right - 5

       WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
       WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
       WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

       WRITE (UNIT=fmtstr2(9:10),FMT="(I2)") ncol
       WRITE (UNIT=fmtstr2(16:17),FMT="(I2)") width - 1
       WRITE (UNIT=fmtstr2(19:20),FMT="(I2)") after

       WRITE (UNIT=fmtstr3(25:26),FMT="(I2)") ncol
       WRITE (UNIT=fmtstr3(32:33),FMT="(I2)") width - 1
       WRITE (UNIT=fmtstr3(35:36),FMT="(I2)") after

       IF (p_evec) THEN

          IF (p_cart) THEN

             IF (my_last) THEN
                name = "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "CARTESIAN MO EIGENVECTORS"
             ELSE
                WRITE (UNIT=name,FMT="(A,I6)")&
                     "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "CARTESIAN MO EIGENVECTORS AFTER SCF STEP",scf_step
             END IF

             ALLOCATE(cmatrix(ncgf,ncgf),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

             cmatrix = 0.0_dp

             ! Transform spherical MOs to Cartesian MOs

             icgf = 1
             isgf = 1
             DO iatom=1,natom
                NULLIFY (orb_basis_set,dftb_parameter)
                CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
                CALL get_qs_kind(qs_kind_set(ikind),&
                     basis_set=orb_basis_set,&
                     dftb_parameter=dftb_parameter,error=error)
                IF (ASSOCIATED(orb_basis_set)) THEN
                   CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                        nset=nset,&
                        nshell=nshell,&
                        l=l)
                   DO iset=1,nset
                      DO ishell=1,nshell(iset)
                         lshell = l(ishell,iset)
                         CALL dgemm("T","N",nco(lshell),nsgf,nso(lshell),1.0_dp,&
                              orbtramat(lshell)%s2c,nso(lshell),&
                              smatrix(isgf,1),nsgf,0.0_dp,&
                              cmatrix(icgf,1),ncgf)
                         icgf = icgf + nco(lshell)
                         isgf = isgf + nso(lshell)
                      END DO
                   END DO
                ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                   CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                   DO ishell=1,lmax+1
                      lshell = ishell-1
                      CALL dgemm("T","N",nco(lshell),nsgf,nso(lshell),1.0_dp,&
                           orbtramat(lshell)%s2c,nso(lshell),&
                           smatrix(isgf,1),nsgf,0.0_dp,&
                           cmatrix(icgf,1),ncgf)
                      icgf = icgf + nco(lshell)
                      isgf = isgf + nso(lshell)
                   END DO
                ELSE
                   CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                        "Unknown basis set type. "//&
                        CPSourceFileRef,&
                        only_ionode=.TRUE.)
                END IF
             END DO ! iatom

          ELSE

             IF (my_last) THEN
                name = "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "SPHERICAL MO EIGENVECTORS"
             ELSE
                WRITE (UNIT=name,FMT="(A,I6)")&
                     "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "SPHERICAL MO EIGENVECTORS AFTER SCF STEP",scf_step
             END IF

          END IF ! p_cart

       ELSE IF (p_occ.OR.p_eval) THEN

          IF (my_last) THEN
             name = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
          ELSE
             WRITE (UNIT=name,FMT="(A,I6)")&
                  "MO EIGENVALUES AND MO OCCUPATION NUMBERS AFTER "//&
                  "SCF STEP",scf_step
          END IF

       END IF ! p_evec

       CALL compress(name)

       ! Print headline

       IF (PRESENT(spin)) THEN
          WRITE (UNIT=iw,FMT="(/,/,T2,A)") spin//" "//TRIM(name)
       ELSE
          WRITE (UNIT=iw,FMT="(/,/,T2,A)") TRIM(name)
       END IF

       ! Check if only a subset of the MOs has to be printed
       IF ((mo_index_range(1) > 0).AND.&
            (mo_index_range(2) > 0).AND.&
            (mo_index_range(2) >= mo_index_range(1))) THEN
          first_mo = MAX(1,mo_index_range(1))
          last_mo = MIN(mo_set%nmo,mo_index_range(2))
       ELSE
          first_mo = 1
          last_mo = mo_set%nmo
       END IF

       IF (p_evec) THEN

          ! Print full MO information

          DO icol=first_mo,last_mo,ncol

             from = icol
             to = MIN((from+ncol-1),last_mo)

             WRITE (UNIT=iw,FMT=fmtstr1) (jcol,jcol=from,to)
             WRITE (UNIT=iw,FMT=fmtstr2) (mo_set%eigenvalues(jcol),jcol=from,to)
             WRITE (UNIT=iw,FMT="(A)") ""

             WRITE (UNIT=iw,FMT=fmtstr2) (mo_set%occupation_numbers(jcol),jcol=from,to)
             WRITE (UNIT=iw,FMT="(A)") ""

             irow = 1

             DO iatom=1,natom

                IF (iatom /= 1) WRITE (UNIT=iw,FMT="(A)") ""

                NULLIFY(orb_basis_set,dftb_parameter)
                CALL get_atomic_kind(particle_set(iatom)%atomic_kind,&
                          element_symbol=element_symbol, kind_number=ikind)
                CALL get_qs_kind(qs_kind_set(ikind),&
                     basis_set=orb_basis_set,&
                     dftb_parameter=dftb_parameter,error=error)

                IF (p_cart) THEN

                   IF (ASSOCIATED(orb_basis_set)) THEN
                      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                           nset=nset,&
                           nshell=nshell,&
                           l=l,&
                           cgf_symbol=bcgf_symbol)

                      icgf = 1
                      DO iset=1,nset
                         DO ishell=1,nshell(iset)
                            lshell = l(ishell,iset)
                            DO ico=1,nco(lshell)
                               WRITE (UNIT=iw,FMT=fmtstr3)&
                                    irow,iatom,ADJUSTR(element_symbol),bcgf_symbol(icgf),&
                                    (cmatrix(irow,jcol),jcol=from,to)
                               icgf = icgf + 1
                               irow = irow + 1
                            END DO
                         END DO
                      END DO
                   ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                      CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                      icgf = 1
                      DO ishell=1,lmax+1
                         lshell = ishell-1
                         DO ico=1,nco(lshell)
                            symbol = cgf_symbol(1,indco(1:3,icgf))
                            symbol(1:2) = "  "
                            WRITE (UNIT=iw,FMT=fmtstr3)&
                                 irow,iatom,ADJUSTR(element_symbol),symbol,&
                                 (cmatrix(irow,jcol),jcol=from,to)
                            icgf = icgf + 1
                            irow = irow + 1
                         END DO
                      END DO
                   ELSE
                      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                           "Unknown basis set type. "//&
                           CPSourceFileRef,&
                           only_ionode=.TRUE.)
                   END IF

                ELSE

                   IF (ASSOCIATED(orb_basis_set)) THEN
                      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                           nset=nset,&
                           nshell=nshell,&
                           l=l,&
                           sgf_symbol=bsgf_symbol)
                      isgf = 1
                      DO iset=1,nset
                         DO ishell=1,nshell(iset)
                            lshell = l(ishell,iset)
                            DO iso=1,nso(lshell)
                               WRITE (UNIT=iw,FMT=fmtstr3)&
                                    irow,iatom,ADJUSTR(element_symbol),bsgf_symbol(isgf),&
                                    (smatrix(irow,jcol),jcol=from,to)
                               isgf = isgf + 1
                               irow = irow + 1
                            END DO
                         END DO
                      END DO
                   ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                      CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                      isgf = 1
                      DO ishell=1,lmax+1
                         lshell = ishell-1
                         DO iso=1,nso(lshell)
                            symbol = sgf_symbol(1,lshell,-lshell+iso-1)
                            symbol(1:2) = "  "
                            WRITE (UNIT=iw,FMT=fmtstr3)&
                                 irow,iatom,ADJUSTR(element_symbol),symbol,&
                                 (smatrix(irow,jcol),jcol=from,to)
                            isgf = isgf + 1
                            irow = irow + 1
                         END DO
                      END DO
                   ELSE
                      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                           "Unknown basis type. "//&
                           CPSourceFileRef,&
                           only_ionode=.TRUE.)
                   END IF

                END IF ! p_cart

             END DO ! iatom

          END DO ! icol

          WRITE (UNIT=iw,FMT="(/)")

          ! Release work storage

          DEALLOCATE (smatrix,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          IF (p_cart) THEN
             DEALLOCATE (cmatrix,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF

       ELSE IF (p_occ.OR.p_eval) THEN

          fmtstr4 = "(T2,I9,2X,F28.  ,1X,F24.  )"
          WRITE (UNIT=fmtstr4(15:16),FMT="(I2)") after
          WRITE (UNIT=fmtstr4(25:26),FMT="(I2)") after
          WRITE (UNIT=iw,FMT="(/,A)")&
               "# MO index          MO eigenvalue [a.u.]            MO occupation"
          DO imo=first_mo,last_mo
             WRITE (UNIT=iw,FMT=fmtstr4)&
                  imo,mo_set%eigenvalues(imo),mo_set%occupation_numbers(imo)
          END DO
          fmtstr5 = "(A,T42,F24.  ,/)"
          WRITE (UNIT=fmtstr5(12:13),FMT="(I2)") after
          WRITE (UNIT=iw,FMT=fmtstr5)&
               "# Sum",accurate_sum(mo_set%occupation_numbers(:))

       END IF ! p_evec

       fmtstr6 = "(A,T17,F24.  ,/)"
       WRITE (UNIT=fmtstr6(12:13),FMT="(I2)") after
       WRITE (UNIT=iw,FMT=fmtstr6) "  Fermi energy:",mo_set%mu
       IF ((mo_set%homo>0)) THEN
         IF ((mo_set%occupation_numbers(mo_set%homo)==mo_set%maxocc).AND.(last_mo > mo_set%homo)) THEN
          gap = mo_set%eigenvalues(mo_set%homo+1) -&
                mo_set%eigenvalues(mo_set%homo)
          fmtstr7 = "(A,T17,F24.  ,A,F6.2,A,/)"
          WRITE (UNIT=fmtstr7(12:13),FMT="(I2)") after
          WRITE (UNIT=iw,FMT=fmtstr7)&
               "  HOMO-LUMO gap:",gap," = ",gap*evolt," eV"
         END IF
       END IF

    END IF ! iw

    CALL cp_print_key_finished_output(iw,logger,dft_section,"PRINT%MO",&
         ignore_should_output=should_output,&
         error=error)

  END SUBROUTINE write_mo_set_to_output_unit

END MODULE qs_mo_io
