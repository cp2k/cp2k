!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/xas_methods   [1.0] *
!!
!!   NAME
!!      xas_methods
!!
!!   FUNCTION
!!      Initialize the XAS orbitals for specific core excitations
!!      Either the GS orbitals are used as initial guess, or the
!!      xas mos are read from a previous calculation.
!!      In the latter case, the core-hole potetial should be the same.
!!
!!   NOTE
!!      The restart with the same core-hole potential should be checked
!!      and a wrong restart should stop the program
!!
!!   AUTHOR
!!     MI (09.2006)
!!
!!   MODIFICATION HISTORY
!!     created 09.2006
!!
!!   SOURCE
!******************************************************************************

MODULE xas_restart


  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_types,                     ONLY: cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast
  USE orbital_pointers,                ONLY: nso
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             make_basis_sm
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             write_mo_set_low
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE string_utilities,                ONLY: xstring
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_restart'

! *** Public subroutines ***

  PUBLIC ::  xas_read_restart, xas_write_restart, xas_initialize_rho


!****************************************************************************

CONTAINS

!!****f* xas_restart/read_mo_set_xas [1.0] *
!!
!!   NAME
!!     read_mo_set_xas
!!
!!   FUNCTION
!!
!!
!!   NOTES
!!     this routine should notbe public, it reads only the mos
!!     for an xas calculation, i.e. for each core excitation from a different
!!     restart file.
!!
!!   ARGUMENTS
!!     - mos: store the molecular orbitalsread fromthe restart file
!!     - qs_env: contains the para_env and info about basis set and atomic kinds
!!     - rst_unit: where to read the mos
!!     - error
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE read_mo_set_restart_xas(mos,qs_env,rst_unit,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: rst_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mo_set_restart_xas', &
      routineP = moduleN//':'//routineN

    INTEGER :: group, homo, i, iatom, irow, iset, iset_read, ishell, &
      ishell_read, iso, ispin, istat, lshell, nao, nao_read, natom, &
      natom_read, nmo, nmo_read, nset, nset_max, nshell_max, nspin, &
      nspin_read, offset_read, source
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info, offset_info
    LOGICAL                                  :: failure, use_this
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer, vecbuffer_read
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!   ---------------------------------------------------------------------------

    NULLIFY (logger, para_env, particle_set)
    NULLIFY (nset_info, nshell,l, nshell_info,nso_info, offset_info)

    CALL get_qs_env(qs_env=qs_env, para_env=para_env, particle_set=particle_set,error=error)
    group = para_env%group
    source = para_env%source
    logger => cp_error_get_logger()

    nspin = SIZE(mos,1)
    nao = mos(1)%mo_set%nao
    natom = SIZE(particle_set,1)

    IF (para_env%ionode) THEN

      READ (rst_unit) natom_read,nspin_read,nao_read,nset_max,nshell_max

      IF (nspin_read /= nspin) THEN
        CALL stop_program(routineP,"changed nspin, not possible")
      END IF
      IF (natom_read /= natom) THEN
        CALL stop_program(routineP,"changed natom, not possible")
      END IF

      ! Let's make it possible tochange the basis set
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
             " READ RESTART XAS : WARNING, DIFFERENT # AOs ",nao,nao_read
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

    END IF

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
        READ (rst_unit) nmo_read,&
                            mos(ispin)%mo_set%homo,&
                            mos(ispin)%mo_set%lfomo,&
                            mos(ispin)%mo_set%nelectron
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
        CALL cp_assert((homo==mos(ispin)%mo_set%homo),cp_warning_level,cp_assertion_failed,routineP,&
             "The homo has changed,  you are probably using a different core-hole potential."//&
CPSourceFileRef,&
             only_ionode=.TRUE.)
        READ (rst_unit) mos(ispin)%mo_set%eigenvalues(1:nmo),&
                            mos(ispin)%mo_set%occupation_numbers(1:nmo)
        IF (mos(ispin)%mo_set%homo > nmo) THEN
          ! can not make this a warning i.e. homo must be smaller than nmo
          ! otherwise e.g. set_mo_occupation will go out of bounds
          CALL stop_program(routineP,"The number of the occupied MOs on the "//&
                                     "restart unit is greater than the "//&
                                     "number of the allocated MOs.")
        END IF
      END IF

      CALL mp_bcast(nmo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%homo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%lfomo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%nelectron,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%eigenvalues,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%occupation_numbers,source,group)

      DO i=1,nmo
        IF (para_env%ionode) THEN
          READ (rst_unit) vecbuffer_read
          ! now, try to assign the read to the real vector
          ! in case the basis set changed this involves some guessing
          irow=1
          DO iatom=1,natom
            NULLIFY(orb_basis_set,l,nshell)
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                                 orb_basis_set=orb_basis_set)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   nshell=nshell,&
                                   l=l)

            use_this = .TRUE.
            iset_read = 1
            DO iset=1,nset
              ishell_read = 1
              DO ishell=1,nshell(iset)
                lshell = l(ishell,iset)
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
        CALL cp_fm_set_submatrix(mos(ispin)%mo_set%mo_coeff,vecbuffer,1,&
             i,nao,1,transpose=.TRUE.,error=error)

      END DO

      ! Skip extra MOs if there any
      IF (para_env%ionode) THEN
        DO i=nmo+1,nmo_read
          READ (rst_unit) vecbuffer_read
        END DO
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

  END SUBROUTINE read_mo_set_restart_xas

  ! ***************************************************************************

!!****f* xas_restart/xas_read_restart [1.0] *
!!
!!   NAME
!!     xas_read_restart
!!
!!   FUNCTION
!!     Set up for reading the restart
!!     corresponing to the excitation of iatom
!!     If the corresponding restart file does not exist
!!     the GS orbitals are used as initial guess
!!
!!   ARGUMENTS
!!     xas_control: control variables
!!     xas_section: input section for XAS calculations
!!     qs_env:
!!     iatom: index of the absorbing atom
!!     estate: index of the core-hole orbital
!!     error:
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     09.2006 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE xas_read_restart(xas_control,xas_section,qs_env,iatom,estate,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(section_vals_type), POINTER         :: xas_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    INTEGER, INTENT(INOUT)                   :: estate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_read_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: group, handle, homo, ia, ie, &
                                                ispin, istate, nspin, &
                                                output_unit, rst_unit, source
    LOGICAL                                  :: failure, file_exists, ionode
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.
    file_exists = .FALSE.

    NULLIFY(matrix_s,mos,occupation_numbers,para_env)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
      output_unit = -1
    END IF

    output_unit = cp_print_key_unit_nr(logger,xas_section,&
         "PRINT%PROGRAM_RUN_INFO",&
         extension=".Log",error=error)

    CALL get_qs_env( qs_env=qs_env, para_env=para_env ,error=error)
    group = para_env%group
    source = para_env%source

    IF(ionode) THEN
      ! check if the restart file exists
      CALL xstring (xas_control%xas_rootfname, ia, ie )
      filename = xas_control%xas_rootfname(ia:ie)//'_at'//&
                 TRIM(ADJUSTL(cp_to_string(iatom)))//'.rst'

      INQUIRE (FILE=filename,EXIST=file_exists)

      ! open file
      IF(file_exists) THEN

         CALL open_file(file_name=TRIM(filename),&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_position="REWIND",&
                     file_status="OLD",&
                     unit_number=rst_unit)

        IF(output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,T20,A,I5,/)")&
             "Read restart file for atom ", iatom

      ELSE
        IF(output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,T10,A,I5,A,/)")&
             "Restart file for atom ", iatom,&
             " not available. Initialization done with GS orbitals"
      END IF

    END IF

    CALL mp_bcast(file_exists,source,group)

    IF(file_exists) THEN
      CALL get_qs_env(qs_env=qs_env, mos=mos, matrix_s=matrix_s,error=error)
      CALL read_mo_set_restart_xas(mos,qs_env,rst_unit,error=error)

      ! close file
      IF(ionode) THEN
        CALL close_file(unit_number=rst_unit)
      END IF

      CALL get_mo_set(mo_set=mos(1)%mo_set,occupation_numbers=occupation_numbers,homo=homo)
      estate = 0
      DO istate = 1,xas_control%nexc_search
        IF(occupation_numbers(istate)==xas_control%occ_estate) THEN
          estate = istate
          EXIT
        END IF
      END DO
      IF(estate==0) THEN
        CALL stop_program(routineP,"Core-hole orbital not found,"//&
        " probably the core-hole potential of the restart file is different")
      END IF
      occupation_numbers(homo) =  xas_control%occ_homo

      nspin = SIZE(mos,1)
      DO ispin = 1,nspin
      ! ortho so that one can restart for different positions (basis sets?)
         NULLIFY(mo_coeff)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff,homo=homo)
         CALL make_basis_sm(mo_coeff,homo,matrix_s(1)%matrix,error=error)
      END DO
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xas_read_restart

  !*****************************************************************************

  SUBROUTINE xas_write_restart(xas_control,xas_section,scf_section,qs_env,iatom,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(section_vals_type), POINTER         :: xas_section, scf_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_write_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, ia, ie, output_unit, &
                                                rst_unit, source
    LOGICAL                                  :: failure, ionode
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    NULLIFY(mos,para_env,particle_set,logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         scf_section,"PRINT%RESTART",error=error),cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,xas_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".Log",error=error)
       
       CALL get_qs_env( qs_env=qs_env, mos=mos, &
            para_env=para_env, particle_set=particle_set ,error=error)
       
       ! Open file
       rst_unit = -1
       IF(ionode) THEN
          ! Check if the restart file exists
          CALL xstring (xas_control%xas_rootfname, ia, ie )
          filename = xas_control%xas_rootfname(ia:ie)//'_at'//&
               TRIM(ADJUSTL(cp_to_string(iatom)))//'.rst'
          
          CALL open_file(file_name=TRIM(filename),&
               file_action="WRITE",&
               file_form="UNFORMATTED",&
               file_status="REPLACE",&
               unit_number=rst_unit)
          IF(output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(/,T10,A,I5,A,A,/)")&
                  "Xas orbitals  for the absorbing atom ", iatom,&
                  " are written in ", filename
             
          END IF
       END IF
       ! Write mos
       CALL write_mo_set_low(mos, particle_set, para_env, rst_unit, error)
       ! Close file
       IF(ionode)  CALL close_file(unit_number=rst_unit)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xas_write_restart

  !*****************************************************************************
!****f* xas_restart/xas_initialize_rho [1.0] *
!!
!!   NAME
!!      xas_initialize_rho
!!
!!   FUNCTION
!!     Once the mos and the occupation numbers are initialized
!!     the electronic density of the excited state can be calclated
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     09-2006 MI created
!!
!!*** *********************************************************************

  SUBROUTINE xas_initialize_rho(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_initialize_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)
    failure=.FALSE.

    NULLIFY(dft_control,mos)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         mos=mos,rho=rho,error=error)

    DO ispin=1,SIZE(mos)
       CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
            density_matrix=rho%rho_ao(ispin)%matrix,error=error)
    END DO

    CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xas_initialize_rho

! *****************************************************************************

END MODULE  xas_restart
