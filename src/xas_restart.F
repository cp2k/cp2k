!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize the XAS orbitals for specific core excitations
!>       Either the GS orbitals are used as initial guess, or the
!>       xas mos are read from a previous calculation.
!>       In the latter case, the core-hole potetial should be the same.
!> \note
!>       The restart with the same core-hole potential should be checked
!>       and a wrong restart should stop the program
!> \par History
!>      created 09.2006
!> \author MI (09.2006)
! *****************************************************************************
MODULE xas_restart

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type,&
                                             cp_fm_write_unformatted
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gspace_mixing,                ONLY: mixing_init
  USE qs_mo_io,                        ONLY: wfn_restart_file_name
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE string_utilities,                ONLY: xstring
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             set_xas_env,&
                                             xas_environment_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_restart'

! *** Public subroutines ***

  PUBLIC ::  xas_read_restart, xas_write_restart, xas_initialize_rho, find_excited_core_orbital

CONTAINS

! *****************************************************************************
!> \brief Set up for reading the restart
!>      corresponing to the excitation of iatom
!>      If the corresponding restart file does not exist
!>      the GS orbitals are used as initial guess
!> \param xas_env ...
!> \param xas_section input section for XAS calculations
!>      qs_env:
!> \param qs_env ...
!> \param xas_method ...
!> \param iatom index of the absorbing atom
!> \param estate index of the core-hole orbital
!>      error:
!> \param error ...
!> \par History
!>      09.2006 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE xas_read_restart(xas_env,xas_section,qs_env, xas_method,iatom,estate,error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(section_vals_type), POINTER         :: xas_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: xas_method, iatom
    INTEGER, INTENT(OUT)                     :: estate
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_read_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER :: group, handle, i, ia, ie, ispin, istat, nao, nao_read, &
      nelectron, nexc_atoms, nexc_atoms_read, nexc_search, nexc_search_read, &
      nmo, nmo_read, output_unit, rst_unit, source, xas_estate, &
      xas_estate_read, xas_method_read
    LOGICAL                                  :: failure, file_exists
    REAL(dp)                                 :: occ_estate, occ_estate_read, &
                                                xas_nelectron, &
                                                xas_nelectron_read
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, &
                                                occupation_numbers
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig_read, occ_read
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)

    failure = .FALSE.
    file_exists = .FALSE.
    rst_unit =  -1

    NULLIFY(eigenvalues, matrix_s, mos, occupation_numbers, vecbuffer)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    output_unit = cp_print_key_unit_nr(logger,xas_section,&
         "PRINT%PROGRAM_RUN_INFO", extension=".Log",error=error)

    CALL get_qs_env( qs_env=qs_env, para_env=para_env ,error=error)
    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN
      CALL wfn_restart_file_name(filename,file_exists,xas_section,logger,&
           xas=.TRUE.,error=error)

      CALL xstring (filename, ia, ie )
      filename = filename(ia:ie)//'-at'//&
                 TRIM(ADJUSTL(cp_to_string(iatom)))//'.rst'

      INQUIRE (FILE=filename,EXIST=file_exists)
      ! open file
      IF (file_exists) THEN

         CALL open_file(file_name=TRIM(filename),&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_position="REWIND",&
                     file_status="OLD",&
                     unit_number=rst_unit)

         IF (output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,T20,A,I5,/)")&
             "Read restart file for atom ", iatom

      ELSE IF (.NOT.file_exists) THEN
          IF(output_unit > 0) WRITE (UNIT=output_unit,FMT="(/,T10,A,I5,A,/)")&
               "Restart file for atom ", iatom,&
               " not available. Initialization done with GS orbitals"
      END IF
    END IF
    CALL mp_bcast(file_exists,source,group)

    CALL get_xas_env(xas_env=xas_env,occ_estate=occ_estate,xas_estate=xas_estate, &
            xas_nelectron=xas_nelectron, nexc_search=nexc_search, nexc_atoms=nexc_atoms,error=error)

    IF(file_exists) THEN
      CALL get_qs_env(qs_env=qs_env, mos=mos, matrix_s=matrix_s,&
            error=error)

      IF(rst_unit>0) THEN
         READ(rst_unit) xas_method_read
         READ(rst_unit) nexc_search_read, nexc_atoms_read, occ_estate_read, xas_nelectron_read
         READ(rst_unit) xas_estate_read

         CALL cp_assert(xas_method_read==xas_method,cp_failure_level,cp_assertion_failed,routineP,&
                 " READ XAS RESTART: restart with different XAS method is not possible. "//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
         CALL cp_assert(nexc_atoms_read==nexc_atoms,cp_failure_level,cp_assertion_failed,routineP,&
                 " READ XAS RESTART: restart with different excited atoms "//&
                 " is not possible. Start instead a new XAS run with the new set of atoms. "//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
      ENDIF

      CALL mp_bcast(xas_estate_read,source,group)
      CALL set_xas_env(xas_env=xas_env,xas_estate=xas_estate_read,error=error)
      estate = xas_estate_read

      CALL get_mo_set(mo_set=mos(1)%mo_set,nao=nao)
      ALLOCATE (vecbuffer(1,nao),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      DO ispin = 1, SIZE(mos)
        CALL get_mo_set(mo_set=mos(ispin)%mo_set, nmo=nmo, eigenvalues=eigenvalues,&
             occupation_numbers=occupation_numbers,mo_coeff=mo_coeff,nelectron=nelectron)
        eigenvalues = 0.0_dp
        occupation_numbers = 0.0_dp
        CALL cp_fm_set_all(mo_coeff,0.0_dp,error=error)
        IF (para_env%ionode) THEN
          READ (rst_unit) nao_read, nmo_read
          CALL cp_assert(nao==nao_read,cp_failure_level,cp_assertion_failed,routineP,&
                 "To change basis is not possible. "//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
          ALLOCATE(eig_read(nmo_read), occ_read(nmo_read), STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          eig_read = 0.0_dp
          occ_read = 0.0_dp
          nmo = MIN(nmo,nmo_read)
          READ (rst_unit) eig_read(1:nmo_read), occ_read(1:nmo_read)
          eigenvalues(1:nmo) = eig_read(1:nmo)
          occupation_numbers(1:nmo) = occ_read(1:nmo)
          IF(nmo_read>nmo) THEN
            CALL cp_assert(occupation_numbers(nmo)<EPSILON(0.0_dp),&
                 cp_warning_level,cp_assertion_failed,routineP,&
                "The number of occupied MOs on the restart unit is larger than "//&
                 "the allocated MOs."//&
CPSourceFileRef,&
                 only_ionode=.TRUE.)

          END IF
          DEALLOCATE(eig_read, occ_read, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        ENDIF
        CALL mp_bcast(eigenvalues,source,group)
        CALL mp_bcast(occupation_numbers,source,group)

         DO i=1,nmo
            IF (para_env%ionode) THEN
               READ (rst_unit) vecbuffer
            ELSE
               vecbuffer(1,:) = 0.0_dp
            END IF
            CALL mp_bcast(vecbuffer,source,group)
            CALL cp_fm_set_submatrix(mo_coeff,&
                 vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)
         END DO
         ! Skip extra MOs if there any
         IF (para_env%ionode) THEN
           DO i=nmo+1,nmo_read
             READ (rst_unit) vecbuffer
           END DO
         END IF

      END DO ! ispin

      DEALLOCATE(vecbuffer,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

!      nspin = SIZE(mos,1)
!      DO ispin = 1,nspin
!      ! ortho so that one can restart for different positions (basis sets?)
!         NULLIFY(mo_coeff)
!         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff,homo=homo)
!         CALL make_basis_sm(mo_coeff,homo,matrix_s(1)%matrix,error=error)
!      END DO
    END IF !file_exist

    IF (para_env%ionode) THEN
      IF (file_exists) CALL close_file(unit_number=rst_unit)
    END IF

    CALL timestop(handle)

  END SUBROUTINE xas_read_restart

! *****************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param xas_section ...
!> \param qs_env ...
!> \param xas_method ...
!> \param iatom ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE xas_write_restart(xas_env,xas_section,qs_env,xas_method,iatom,error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(section_vals_type), POINTER         :: xas_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: xas_method, iatom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_write_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: my_middle
    INTEGER                                  :: handle, ispin, nao, &
                                                nexc_atoms, nexc_search, nmo, &
                                                output_unit, rst_unit, &
                                                xas_estate
    LOGICAL                                  :: failure
    REAL(dp)                                 :: occ_estate, xas_nelectron
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, &
                                                occupation_numbers
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(mos,logger,print_key)
    logger => cp_error_get_logger(error)

    CALL get_xas_env(xas_env=xas_env,occ_estate=occ_estate,xas_estate=xas_estate, &
            xas_nelectron=xas_nelectron, nexc_search=nexc_search, nexc_atoms=nexc_atoms, error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              xas_section,"PRINT%RESTART",used_print_key=print_key,error=error),&
              cp_p_file)) THEN

       output_unit = cp_print_key_unit_nr(logger,xas_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".Log",error=error)

       CALL get_qs_env( qs_env=qs_env, mos=mos, error=error)

       ! Open file
       rst_unit = -1
       my_middle= 'at'//TRIM(ADJUSTL(cp_to_string(iatom)))
       rst_unit = cp_print_key_unit_nr(logger,xas_section,"PRINT%RESTART",&
                  extension=".rst", file_status="REPLACE", file_action="WRITE",&
                  file_form="UNFORMATTED",middle_name=TRIM(my_middle), error=error)

       filename = cp_print_key_generate_filename(logger,print_key,&
                  middle_name=TRIM(my_middle),extension=".rst",&
                  my_local=.FALSE.,error=error)

       IF(output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T10,A,I5,A,A,/)")&
               "Xas orbitals  for the absorbing atom ", iatom,&
               " are written in ", TRIM(filename)

       END IF

       ! Write mos
       IF(rst_unit>0) THEN
         WRITE(rst_unit)  xas_method
         WRITE(rst_unit)  nexc_search, nexc_atoms, occ_estate, xas_nelectron
         WRITE(rst_unit)  xas_estate
       ENDIF
       DO ispin=1,SIZE(mos)
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo,&
              eigenvalues=eigenvalues,occupation_numbers=occupation_numbers)
         IF ((rst_unit>0)) THEN
           WRITE (rst_unit) nao, nmo
           WRITE (rst_unit) eigenvalues(1:nmo),&
                           occupation_numbers(1:nmo)
         END IF
         CALL cp_fm_write_unformatted(mo_coeff,rst_unit,error)
       END DO

! Close file
       CALL cp_print_key_finished_output(rst_unit,logger,xas_section,&
                  "PRINT%RESTART", error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE xas_write_restart

!****f* xas_restart/xas_initialize_rho [1.0] *

! *****************************************************************************
!> \brief Once the mos and the occupation numbers are initialized
!>      the electronic density of the excited state can be calclated
!> \param qs_env ...
!> \param error ...
!> \par History
!>      09-2006 MI created
!> \author MI
! *****************************************************************************
  SUBROUTINE xas_initialize_rho(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_initialize_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom
    TYPE(xas_environment_type), POINTER      :: xas_env

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(mos,rho,xas_env,para_env)

    CALL get_qs_env(qs_env,&
                    mos=mos,&
                    rho=rho,&
                    xas_env=xas_env,&
                    para_env=para_env,&
                    error=error)

    DO ispin=1,SIZE(mos)
       IF(ispin==1) THEN
         CALL set_mo_occupation(mo_set=qs_env%mos(ispin)%mo_set, smear=xas_env%smear,&
              xas_env=xas_env, error=error)
       ELSE
         CALL set_mo_occupation(mo_set=qs_env%mos(ispin)%mo_set, smear=xas_env%smear,&
              error=error)
       END IF
       CALL calculate_density_matrix(mo_set=mos(ispin)%mo_set,&
            density_matrix=rho%rho_ao(ispin)%matrix,error=error)

    END DO

    CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)

    IF(xas_env%mixing_method>1) THEN
       CALL get_qs_env(qs_env=qs_env,rho_atom_set=rho_atom,error=error)
       CALL mixing_init(xas_env%mixing_method,rho,xas_env%mixing_store,&
                           para_env,rho_atom=rho_atom,error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE xas_initialize_rho

! *****************************************************************************
!> \brief Find the index of the core orbital that has been excited by XAS
!> \param xas_env ...
!> \param mos ...
!> \param matrix_s ...
!> \param error ...
!> \par History
!>      03-2010 MI created
!> \author MI
! *****************************************************************************

  SUBROUTINE find_excited_core_orbital(xas_env,mos,matrix_s,error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'find_excited_core_orbital', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ic_max, ir_max, istat, m, &
                                                n, nao, nexc_search, nmo, &
                                                xas_estate
    INTEGER, DIMENSION(:), POINTER           :: col_indices
    LOGICAL                                  :: failure
    REAL(dp)                                 :: a_max, b_max, ip_energy, &
                                                occ_estate
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues, &
                                                occupation_numbers
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer, vecbuffer2
    TYPE(cp_fm_type), POINTER                :: excvec_coeff, excvec_overlap, &
                                                fm_work, mo_coeff

    failure = .FALSE.
    NULLIFY(excvec_coeff,excvec_overlap,fm_work,mo_coeff)
    ! Some elements from the xas_env
    CALL get_xas_env(xas_env=xas_env,excvec_coeff=excvec_coeff,&
         excvec_overlap=excvec_overlap,fm_work=fm_work,nexc_search=nexc_search,&
         xas_estate=xas_estate,occ_estate=occ_estate,error=error)
    CPPrecondition(ASSOCIATED(excvec_overlap),cp_failure_level,routineP,error,failure)

    CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
         eigenvalues=eigenvalues,occupation_numbers=occupation_numbers)
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer = 0.0_dp
    ALLOCATE(vecbuffer2(1,nexc_search),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer2 = 0.0_dp

  ! ** use the maximum overlap criterion to find the index of the excited orbital
    CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,mo_coeff,fm_work,ncol=nmo,error=error)
    CALL cp_gemm("T","N",1,xas_env%nexc_search,nao,1.0_dp,excvec_coeff,&
         fm_work,0.0_dp,excvec_overlap,b_first_col=1,error=error)
    CALL cp_fm_get_info ( matrix=excvec_overlap, col_indices = col_indices,&
         nrow_global = m, ncol_global = n ,error=error)
    CALL cp_fm_get_submatrix(excvec_overlap,vecbuffer2,1,1,&
         1,nexc_search,transpose=.FALSE.,error=error)

    b_max = 0.0_dp
    ic_max = xas_estate
    DO i = 1,nexc_search
       a_max = ABS(vecbuffer2(1,i))
       IF(a_max > b_max)THEN
          ic_max = i
          b_max = a_max
       ENDIF
    END DO

    IF(ic_max /= xas_estate) THEN
       ir_max =  xas_estate
       xas_estate = ic_max
       occupation_numbers(xas_estate) = occ_estate
       occupation_numbers(ir_max) = 1.0_dp
    END IF

    ! Ionization Potential
    iP_energy = eigenvalues(xas_estate)
    CALL set_xas_env(xas_env=xas_env, xas_estate=xas_estate, ip_energy=ip_energy,error=error)

    CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,xas_estate,&
         nao,1,transpose=.TRUE.,error=error)
    CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
         nao,1,transpose=.TRUE.,error=error)

    DEALLOCATE(vecbuffer,vecbuffer2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

 END SUBROUTINE find_excited_core_orbital

END MODULE xas_restart
