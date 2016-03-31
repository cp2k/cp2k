!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for reading and writing restart files.
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_io
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE cell_types,                      ONLY: cell_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_col_block_sizes,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_row_block_sizes
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE pao_input,                       ONLY: id2str
  USE pao_param,                       ONLY: pao_param_count
  USE pao_types,                       ONLY: pao_env_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_io'

  PUBLIC :: pao_read_restart, pao_write_restart

CONTAINS

! *****************************************************************************
!> \brief Reads restart file
!> \param pao ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE pao_read_restart(pao, qs_env)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pao_read_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: label
    INTEGER                                  :: i, iatom, natoms, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, buffer
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input

    CALL get_qs_env(qs_env, input=input, para_env=para_env)

    CPASSERT(LEN_TRIM(pao%restart_file)>0)
    IF(pao%iw>0) WRITE(pao%iw,'(A,A)') " PAO| Reading restart file: ", TRIM(pao%restart_file)

    row_blk_sizes => cp_dbcsr_row_block_sizes(pao%matrix_X)
    col_blk_sizes => cp_dbcsr_col_block_sizes(pao%matrix_X)
    IF(SIZE(row_blk_sizes) /= SIZE(col_blk_sizes)) CPABORT("matrix_X not squared")
    natoms = SIZE(col_blk_sizes)

    unit_nr = -1
    IF(para_env%mepos == para_env%source) THEN
       CALL open_file(file_name=pao%restart_file, file_status="OLD", file_form="FORMATTED",&
                      file_action="READ", unit_number=unit_nr)

       CALL read_restart_header(pao, qs_env, unit_nr)
    ENDIF

    !TODO: this is a serial algorithm
    DO i=1, natoms
       IF(unit_nr>0) THEN
          READ(unit_nr,fmt=*) label, iatom
          CPASSERT(TRIM(label)=="Xblock")
       ENDIF
       CALL mp_bcast(iatom, para_env%source, para_env%group)
       IF(pao%iw>0)  WRITE(pao%iw,*) "PAO| Restart found iatom: ", iatom
       CPASSERT(iatom <= natoms)

       ALLOCATE(buffer(row_blk_sizes(iatom), col_blk_sizes(iatom)))
       IF(unit_nr>0) THEN
          BACKSPACE(unit_nr)
          READ(unit_nr,fmt=*) label, iatom, buffer
       ENDIF
       CALL mp_bcast(buffer, para_env%source, para_env%group)
       CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=iatom, col=iatom, block=block_X, found=found)
       IF(ASSOCIATED(block_X)) THEN
          block_X = buffer
       ENDIF
       DEALLOCATE(buffer)
    ENDDO

    IF(unit_nr>0) CALL close_file(unit_number=unit_nr)
  END SUBROUTINE pao_read_restart


! *****************************************************************************
!> \brief Read and check header of restart file
!> \param pao ...
!> \param qs_env ...
!> \param unit_nr ...
! *****************************************************************************
  SUBROUTINE read_restart_header(pao, qs_env, unit_nr)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'read_restart_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kindname, label, str_in
    INTEGER                                  :: i1, i2, i3, iatom, ikind, &
                                                pao_basis_size, pao_pot_maxl, &
                                                pao_pot_neighbors, z
    REAL(dp)                                 :: diff, pao_pot_beta, r1
    REAL(dp), DIMENSION(3)                   :: pos_in
    REAL(dp), DIMENSION(3, 3)                :: hmat_in
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

   CALL get_qs_env(qs_env,&
                   cell=cell,&
                   particle_set=particle_set,&
                   atomic_kind_set=atomic_kind_set,&
                   qs_kind_set=qs_kind_set)

   DO WHILE(.TRUE.)
       READ(unit_nr, fmt=*) label
       BACKSPACE(unit_nr)

       IF(TRIM(label) == "Parametrization") THEN
          READ(unit_nr,fmt=*) label, str_in
          IF(TRIM(str_in) .NE. TRIM(ADJUSTL(id2str(pao%parameterization))))&
             CPABORT("Restart PAO parametrization does not match")

       ELSE IF(TRIM(label) == "Cell") THEN
          READ(unit_nr,fmt=*) label, hmat_in
          diff = MAXVAL(ABS(hmat_in - cell%hmat * angstrom))
          IF(diff > 1e-14)&
             CPWARN("Restaring from different cell dimension")

       ELSE IF(TRIM(label) == "Kind") THEN
          READ(unit_nr, fmt=*) label, ikind, str_in, i1
          CALL get_atomic_kind(atomic_kind_set(ikind), name=kindname, z=z)
          IF(str_in .NE. kindname)&
             CPABORT("Kind names do not match")
          IF(i1 /= z)&
             CPABORT("Atomic numbers do not match")

       ELSE IF(TRIM(label) == "PrimBasis") THEN
          READ(unit_nr, fmt=*) label, ikind, str_in
          CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set)
          IF(TRIM(str_in) .NE. TRIM(basis_set%name))&
             CPABORT("Primary Basis set does not match")

       ELSE IF(TRIM(label) == "PaoBasis") THEN
          READ(unit_nr, fmt=*) label, ikind, i1, i2, i3, r1
          CALL get_qs_kind(qs_kind_set(ikind),&
                           pao_basis_size=pao_basis_size,&
                           pao_potential_maxl=pao_pot_maxl,&
                           pao_potential_neighbors=pao_pot_neighbors,&
                           pao_potential_beta=pao_pot_beta)
          IF(i1 /= pao_basis_size)&
             CPABORT("PAO_BASIS_SIZE does not match")
          IF(i2 /= pao_pot_maxl)&
             CPABORT("PAO_POT_MAXL does not match")
          IF(i3 /= pao_pot_neighbors)&
             CPABORT("PAO_POT_NEIGHBORS does not match")
          IF(r1 /= pao_pot_beta)&
             CPABORT("PAO_POT_BETA does not match")

       ELSE IF(TRIM(label) == "Atom") THEN
          READ(unit_nr,fmt=*) label, iatom, str_in, pos_in
          IF(str_in .NE. particle_set(iatom)%atomic_kind%name)&
             CPABORT("Restart atomic kinds do not match.")
          diff = MAXVAL(ABS(pos_in - particle_set(iatom)%r * angstrom))
          IF(diff > 1e-10)&
             CPWARN("Restarting from different atom positions")

       ELSE IF(TRIM(label) == "Xblock") THEN
          EXIT ! end of header

       ELSE
          !CPWARN("Skipping restart header with label: "//TRIM(label))
          READ(unit_nr,fmt=*) label  ! just read again and ignore
       ENDIF
   ENDDO
 END SUBROUTINE read_restart_header


! *****************************************************************************
!> \brief Writes restart file
!> \param pao ...
!> \param qs_env ...
!> \param energy ...
! *****************************************************************************
  SUBROUTINE pao_write_restart(pao, qs_env, energy)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: energy

    CHARACTER(len=*), PARAMETER :: &
      printkey_section = 'DFT%LS_SCF%PAO%PRINT%RESTART', &
      routineN = 'pao_write_restart', routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, natoms, &
                                                unit_max, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: block_X, buffer
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()

    CALL get_qs_env(qs_env,&
                    input=input,&
                    natom=natoms,&
                    para_env=para_env)

    ! open file
    unit_nr = cp_print_key_unit_nr(logger,&
                              input,&
                              printkey_section,&
                              extension=".pao",&
                              file_action="WRITE",&
                              file_position="REWIND",&
                              file_status="UNKNOWN",&
                              do_backup=.TRUE.)

    ! although just rank-0 writes the trajectory it requires collective MPI calls
    unit_max = unit_nr
    CALL mp_max(unit_max, para_env%group)
    IF(unit_max>0) THEN
       IF(pao%iw>0) WRITE(pao%iw,'(A,A)') " PAO| Writing restart file."
       IF(unit_nr>0)&
          CALL write_restart_header(pao, qs_env, energy, unit_nr)

       !TODO: this is a serial algorithm
       ! write matrix_X
       row_blk_sizes => cp_dbcsr_row_block_sizes(pao%matrix_X)
       col_blk_sizes => cp_dbcsr_col_block_sizes(pao%matrix_X)
       DO iatom=1, natoms
          ALLOCATE(buffer(row_blk_sizes(iatom), col_blk_sizes(iatom)))
          NULLIFY(block_X)
          CALL cp_dbcsr_get_block_p(matrix=pao%matrix_X, row=iatom, col=iatom, block=block_X, found=found)
          IF(ASSOCIATED(block_X)) THEN
             IF(SIZE(block_X) > 0) & ! corner-case of zero pao parameters
               buffer(:,:) = block_X(:,:)
          ELSE
             buffer(:,:) = 0.0_dp
          ENDIF
          CALL mp_sum(buffer, para_env%group)
          IF(unit_nr>0) THEN
             WRITE(unit_nr, fmt="(A,I10,1X)",advance='no') "Xblock ", iatom
             WRITE(unit_nr, *) buffer
          ENDIF
          DEALLOCATE(buffer)
       ENDDO

       ! flush
       IF(unit_nr>0) FLUSH(unit_nr)

    ENDIF

    ! close file
    CALL cp_print_key_finished_output(unit_nr, logger, input, printkey_section)

    CALL timestop(handle)
  END SUBROUTINE pao_write_restart


! *****************************************************************************
!> \brief Writes header of restart file
!> \param pao ...
!> \param qs_env ...
!> \param energy ...
!> \param unit_nr ...
! *****************************************************************************
  SUBROUTINE write_restart_header(pao, qs_env, energy, unit_nr)
    TYPE(pao_env_type), POINTER              :: pao
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: energy
    INTEGER, INTENT(IN)                      :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'write_restart_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kindname
    INTEGER                                  :: iatom, ikind, nparams, &
                                                pao_basis_size, pao_pot_maxl, &
                                                pao_pot_neighbors, z
    REAL(dp)                                 :: pao_pot_beta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_type), POINTER        :: basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL get_qs_env(qs_env,&
                    cell=cell,&
                    particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set)


    WRITE(unit_nr, "(A,5X,F20.10)") "Energy", energy
    WRITE(unit_nr, "(A,5X,I0)") "Step", pao%istep
    WRITE(unit_nr,"(A,5X,A10)") "Parametrization", id2str(pao%parameterization)

    ! write kinds
    WRITE(unit_nr, "(A,5X,I0)") "Nkinds", SIZE(atomic_kind_set)
    DO ikind=1, SIZE(atomic_kind_set)
       CALL get_atomic_kind(atomic_kind_set(ikind), name=kindname, z=z)
       CALL get_qs_kind(qs_kind_set(ikind),&
                        pao_basis_size=pao_basis_size,&
                        pao_potential_maxl=pao_pot_maxl,&
                        pao_potential_neighbors=pao_pot_neighbors,&
                        pao_potential_beta=pao_pot_beta,&
                        basis_set=basis_set)
       CALL pao_param_count(pao, qs_env, ikind, nparams)
       WRITE(unit_nr,"(A,5X,I10,1X,A,1X,I3)") "Kind", ikind, TRIM(kindname), z
       WRITE(unit_nr,"(A,5X,I10,1X,A)") "PrimBasis", ikind, TRIM(basis_set%name)
       WRITE(unit_nr,"(A,5X,I10,1X,I3,1X,I3)",advance='no') "PaoBasis", ikind, pao_basis_size, pao_pot_maxl
       WRITE(unit_nr,"(1X,I3,1X,F20.16)") pao_pot_neighbors, pao_pot_beta
       WRITE(unit_nr,"(A,5X,I10,1X,I3)") "NParams", ikind, nparams
    ENDDO

    ! write cell
    WRITE(unit_nr,fmt="(A,5X)",advance='no') "Cell"
    WRITE(unit_nr,*) cell%hmat * angstrom

    ! write atoms
    WRITE(unit_nr, "(A,5X,I0)") "Natoms", SIZE(particle_set)
    DO iatom=1, SIZE(particle_set)
       kindname = particle_set(iatom)%atomic_kind%name
       WRITE(unit_nr,fmt="(A,5X,I10,5X,A,1X)",advance='no') "Atom ", iatom, TRIM(kindname)
       WRITE(unit_nr,*) particle_set(iatom)%r * angstrom
    ENDDO

  END SUBROUTINE write_restart_header

END MODULE pao_io
