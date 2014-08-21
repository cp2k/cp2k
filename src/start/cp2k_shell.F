!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interactive shell of CP2K
!> \note
!>       sample of a simple runner that uses the f77_interface
!>       it can beused to connect c programs, communicating through std-in/ std-out
!>
!>       positions are in angstrom, energies in evolt
!>
!>       commands:
!>       load filename:   loads the filename, returns the env_id, or -1 in case of error
!>       natom [env_id]:  returns the number of atoms in the environement env_id
!>                        (defaults to the last loaded)
!>       setpos [env_id]: sets the positions of the atoms, should be followed
!>                        by natom*3 (on a line) and then all the positions [angstrom]
!>       getpos [env_id]: gets the positions of the atoms, returns
!>                        natom*3 (on a line) and then all the positions [angstrom]
!>       calcE [env_id]:  calculate the energy and returns it (in eV)
!>       calcEF [env_id]: calculate the energy and forces and returns it,
!>                        first the energy on a line (in eV), then the natom*3 (on a line)
!>                        and finally all the values (in eV/angstrom)
!> \author Fawzi Mohamed
! *****************************************************************************
PROGRAM cp2k_shell
  USE cp2k_info,                       ONLY: compile_arch,&
                                             compile_date,&
                                             compile_host,&
                                             compile_revision,&
                                             cp2k_home,&
                                             cp2k_version,&
                                             id_cp2k_version,&
                                             print_cp2k_license
  USE cp_files,                        ONLY: close_file,&
                                             get_unit_number
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_interface,                   ONLY: &
       calc_energy_force, calc_force, create_force_env, destroy_force_env, &
       finalize_cp2k, get_energy, get_force, get_natom, get_pos, init_cp2k, &
       set_pos
  USE input_constants,                 ONLY: id_development_version
  USE input_cp2k,                      ONLY: create_cp2k_root_section
  USE input_section_types,             ONLY: section_release,&
                                             section_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_chdir,&
                                             m_flush,&
                                             m_getcwd,&
                                             m_getlog,&
                                             m_getpid,&
                                             m_hostnm
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_bcast
  USE physcon,                         ONLY: angstrom,&
                                             evolt
  USE string_utilities,                ONLY: compress,&
                                             uppercase
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  LOGICAL                                  :: eof, harsh,failure
  INTEGER                                  :: ierr,i,iostat,shift,&
                                              shift2,env_id,last_env_id,&
                                              n_atom,stat,n_atom2,pid
  INTEGER, PARAMETER                       :: sout=6 ! ,sin=5
  TYPE(cp_error_type)                      :: error
  TYPE(cp_para_env_type), POINTER          :: para_env
  TYPE(cp_logger_type), POINTER            :: logger
  CHARACTER(LEN=default_path_length)       :: cmdStr, inp_filename, out_filename
  CHARACTER(LEN=default_string_length)     :: cmd,units
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: pos,pos2
  REAL(KIND=dp)                            :: e_pot,err
  REAL(KIND=dp)                            :: e_fact,pos_fact

  LOGICAL, SAVE                            :: did_init = .FALSE.
  INTEGER, SAVE                            :: eof_stat
  TYPE(section_type), POINTER              :: input_declaration

  CHARACTER(LEN=default_string_length), DIMENSION(37) :: helpMsg

  helpMsg=(/ &
  ll('Commands'),&
  ll(' '),&
  ll(' If there is [env_id] it means that an optional env_id can be given,'),&
  ll(' if none is given it defaults to the last env_id loaded'),&
  ll(' All commands are case insensitive.'),&
  ll(' '),&
  ll(' INFO: returns some information about cp2k.'),&
  ll(' LOAD inp-filename [out-filename]: loads the filename, returns the env_id, or -1 in case of error'),&
  ll('   out-filename is optional and defaults to <inp-filename>.out'),&
  ll('   use "__STD_OUT__" for printing to the screen'),&
  ll(' BG_LOAD filename: loads the filename, without returning the env_id'),&
  ll(' LAST_ENV_ID: returns the env_id of the last environment loaded'),&
  ll(' DESTROY [env_id]: destroys the given environment (last and default env'),&
  ll('   might become invalid)'),&
  ll(' NATOM [env_id]: returns the number of atoms in the environement env_id'),&
  ll(' SET_POS [env_id]: sets the positions of the atoms, should be followed'),&
  ll('   by natom*3 (on a line) and then all the positions. Returns the max'),&
  ll('   change of the coordinates (useful to avoid extra calculations).'),&
  ll(' GET_POS [env_id]: gets the positions of the atoms, returns'),&
  ll('   natom*3 (on a line) and then all the positions then "* END" '),&
  ll('   (alone on a line)'),&
  ll(' GET_E [env_id]: gets the energy of the last calculation on env_id'),&
  ll(' GET_F [env_id]: gets the forces on the atoms,of the last calculation on '),&
  ll('   env_id, if only the energy was calculated the content is undefined. Returns'),&
  ll('   natom*3 (on a line) and then all the forces then "* END" (alone on'),&
  ll('   a line)'),&
  ll(' CALC_E [env_id]: calculate the energy and returns it'),&
  ll(' EVAL_E [env_id]: calculate the energy (without returning it)'),&
  ll(' CALC_EF [env_id]: calculate energy and forces and returns them,'),&
  ll('   first the energy on a line, then the natom*3 (on a line)'),&
  ll('   and finally all the values and "* END" (alone on a line)'),&
  ll(' EVAL_EF [env_id]: calculate the energy and forces (without returning them)'),&
  ll(' HARSH: stops on any error'),&
  ll(' PERMISSIVE: stops only on serious errors'),&
  ll(' UNITS: returns the units used for energy and position'),&
  ll(' UNITS_EV_A: sets the units to electron volt (energy)  and Angstrom (positions)'),&
  ll(' HELP: writes the present help') /)

  CALL init_cp2k(init_mpi=.TRUE.,ierr=ierr)
  IF (ierr/=0) CALL mp_abort('init_cp2k failure')

  NULLIFY(input_declaration)
  CALL create_cp2k_root_section(input_declaration, error)

  pos_fact=1.0_dp
  e_fact=1.0_dp
  units='au'
  harsh=.FALSE.
  failure=.FALSE.

  CALL cp_error_init(error)
  logger => cp_error_get_logger(error)
  para_env => logger%para_env

  eof=.FALSE.
  last_env_id=-1
  DO
     IF (para_env%mepos==para_env%source) THEN
        WRITE (sout,'("* READY")')
     END IF
     CALL m_flush(sout)
     IF (para_env%mepos==para_env%source) THEN
        READ(*,'(a)',iostat=iostat)cmdStr
        IF (iostat/=0) THEN
           IF (iostat==get_eof_stat()) THEN
              WRITE(sout,'(a)',iostat=iostat) '* EOF'
              ! ignoring iostat on purpose
           END IF
           eof=.TRUE.
        END IF
     END IF
     CALL mp_bcast(eof,para_env%source,para_env%group)
     IF (eof) EXIT
     CALL mp_bcast(cmdStr,para_env%source,para_env%group)
     shift=1
     DO
        IF(cmdStr(shift:shift)/=' ') EXIT
        shift=shift+1
        IF (shift>=LEN(cmdStr)) EXIT
     END DO
     shift2=shift
     DO
        IF(cmdStr(shift2:shift2)==' ') EXIT
        shift2=shift2+1
        IF (shift2>=LEN(cmdStr)) EXIT
     END DO
     shift2=shift2-1
     CALL my_assert(shift<=shift2,'* No command',failure)
     IF (failure) GOTO 10
     cmd=cmdStr(shift:shift2)
     CALL uppercase(cmd)
     SELECT CASE(cmd)
     CASE('HELP')
        IF (para_env%source==para_env%mepos) THEN
            DO i=1,SIZE(helpMsg)
                WRITE(sout,'(a)',iostat=iostat) TRIM(helpMsg(i))
                IF (iostat/=0) CALL mp_abort('writing help')
            END DO
            CALL m_flush(sout)
        END IF
     CASE('INFO','INFORMATION','LICENSE')
        IF (para_env%source==para_env%mepos) THEN
            CALL print_v_info(sout)
            CALL print_cp2k_license(sout)
        END IF
     CASE('LAST_ENV_ID')
        IF (para_env%source==para_env%mepos) THEN
           WRITE (sout,'(i10)',iostat=iostat) last_env_id
           IF (iostat/=0) CALL mp_abort(cmd//' failed writing last_env_id')
           CALL m_flush(sout)
        END IF
     CASE('BG_LOAD','BGLOAD')
        shift=shift2+1
        CALL my_assert(shift<=LEN(cmdStr),'load without file',failure)
        IF (failure) GOTO 10
        DO
           IF (cmdStr(shift:shift)/=' ') EXIT
           shift=shift+1
           IF (shift>=LEN(cmdStr)) EXIT
        END DO
        CALL my_assert(cmdStr(shift:shift)/=' ','load without file2',failure)
        IF (failure) GOTO 10
        CALL create_force_env(new_env_id=last_env_id,&
             input_declaration=input_declaration,&
             input_path=cmdStr(shift:LEN_TRIM(cmdStr)),&
             output_path=cmdStr(shift:LEN_TRIM(cmdStr))//'.out',&
             owns_out_unit=.TRUE.,ierr=ierr)
        IF (ierr/=0) last_env_id=-1
        IF (last_env_id<=0) failure=.TRUE.
     CASE('LOAD')
        ! find beginning of inp-filename
        shift=shift2+1
        CALL my_assert(shift<=LEN(cmdStr),'load without file',failure)
        IF (failure) GOTO 10
        DO
           IF (cmdStr(shift:shift)/=' ') EXIT
           shift=shift+1
           IF (shift>=LEN(cmdStr)) EXIT
        END DO
        CALL my_assert(cmdStr(shift:shift)/=' ','load without file2',failure)

        ! find end of inp-filename
        shift2 = shift+1
        DO WHILE(shift2 < LEN(cmdStr))
           shift2 = shift2+1
           IF (cmdStr(shift2:shift2)==' ') EXIT
        END DO
        inp_filename = TRIM(cmdStr(shift:shift2))
        out_filename = TRIM(inp_filename)//'.out'
        IF(shift2 < LEN(cmdStr) .AND. LEN_TRIM(cmdStr(shift2:)) > 0) &
           out_filename = ADJUSTL(cmdStr(shift2:))

        CALL create_force_env(new_env_id=last_env_id,&
             input_declaration=input_declaration,&
             input_path=inp_filename,&
             output_path=out_filename,&
             owns_out_unit=.TRUE.,ierr=ierr)
        IF (ierr/=0) last_env_id=-1
        IF (para_env%source==para_env%mepos) THEN
           WRITE (sout,'(i10)',iostat=iostat) last_env_id
           IF (iostat/=0) CALL mp_abort(cmd//' failed writing last_env_id')
           CALL m_flush(sout)
        END IF
        IF (last_env_id<=0) failure=.TRUE.
     CASE ('DESTROY')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        IF (failure) GOTO 10
        CALL destroy_force_env(env_id,ierr)
        CALL my_assert(ierr==0,cmd//' failed destroy_force_env',failure)
        IF (last_env_id==env_id) last_env_id=-1
     CASE('NATOM','N_ATOM')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        IF (failure) GOTO 10
        CALL get_natom(env_id, n_atom, ierr)
        CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
        IF (.NOT. failure .AND. para_env%mepos==para_env%source) THEN
           WRITE (sout,'(i10)',iostat=iostat) n_atom
           CALL m_flush(sout)
        END IF
     CASE('SETPOS','SET_POS')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        IF (failure) GOTO 10
        CALL get_natom(env_id, n_atom, ierr)
        CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
        IF (failure) GOTO 10
        ALLOCATE(pos(3*n_atom),pos2(3*n_atom),stat=stat)
        CALL my_assert(stat==0,cmd//' failed alloc',failure)
        IF (para_env%mepos==para_env%source) THEN
           READ (*,*,iostat=iostat) n_atom2
           IF (iostat/=0) CALL mp_abort('setpos read n_atom')
           IF (n_atom2/=SIZE(pos)) THEN
              CALL my_assert(.FALSE.,'setpos invalid number of atoms',failure)
              DO i=1,n_atom
                 READ(*,'(a)',iostat=iostat) cmdStr
                 CALL compress(cmdStr,full=.TRUE.)
                 CALL uppercase(cmdStr)
                 IF (cmdStr=='*END') EXIT
              END DO
              GOTO 10
           END IF
           READ (*,*,iostat=iostat) pos
           IF (iostat/=0) CALL mp_abort('setpos read coord')
           pos=pos/pos_fact
           READ(*,'(a)',iostat=iostat) cmdStr
           CALL compress(cmdStr,full=.TRUE.)
           CALL uppercase(cmdStr)
           CALL my_assert(cmdStr=='*END',' missing *END',failure)
        END IF
        CALL mp_bcast(pos,para_env%source,para_env%group)
        CALL get_pos(env_id, pos2,n_el=3*n_atom,ierr=ierr)
        CALL my_assert(ierr==0,'get_pos error',failure)
        CALL set_pos(env_id, new_pos=pos, n_el=3*n_atom, ierr=ierr)
        CALL my_assert(ierr==0,'set_pos error',failure)
        err=0.0_dp
        DO i=1,SIZE(pos)
           err=MAX(err,ABS(pos(i)-pos2(i)))
        END DO
        DEALLOCATE(pos,pos2,stat=stat)
        CALL my_assert(stat==0,'dealloc error',failure)
        IF (para_env%mepos==para_env%source) THEN
           WRITE (sout,'(ES22.13)') err*pos_fact
        END IF
        IF (stat/=0) CALL mp_abort('setpos failed dealloc')
     CASE('GETPOS','GET_POS')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        IF (failure) GOTO 10
        CALL get_natom(env_id, n_atom, ierr)
        CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
        IF (failure) GOTO 10
        ALLOCATE(pos(3*n_atom),stat=stat)
        IF (stat/=0) CALL mp_abort(cmd//' failed alloc')
        CALL get_pos(env_id, pos=pos, n_el=3*n_atom, ierr=ierr)
        CALL my_assert(ierr==0,'get_pos error',failure)
        pos=pos*pos_fact
        IF (.not.failure.and.para_env%mepos==para_env%source) THEN
           WRITE (sout,'(i10)',iostat=iostat) 3*n_atom
           IF (iostat/=0) CALL mp_abort(cmd//' write ncoord')
           WRITE (sout,'(3ES22.13)',iostat=iostat) pos
           IF (iostat/=0) CALL mp_abort(cmd//' write coord')
           WRITE (sout,'(a)',iostat=iostat) "* END"
           IF (iostat/=0) CALL mp_abort(cmd//' write coord')
           CALL m_flush(sout)
        END IF
        DEALLOCATE(pos,stat=stat)
        IF (stat/=0) CALL mp_abort(cmd//' failed dealloc')
     CASE('GETE','GET_E')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        CALL get_energy(env_id,e_pot,ierr)
        CALL my_assert(ierr==0,cmd//' failed get_energy',failure)
        e_pot=e_pot*e_fact
        IF (para_env%mepos==para_env%source) THEN
           WRITE (sout,'(ES22.13)',iostat=iostat) e_pot
           IF (iostat/=0) CALL mp_abort(cmd//' write e_pot')
           CALL m_flush(sout)
        END IF
     CASE('EVALE','EVAL_E')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        CALL calc_energy_force(env_id,calc_force=.FALSE.,ierr=ierr)
        CALL my_assert(ierr==0,cmd//' calc_energy_force failed',failure)
     CASE('CALCE','CALC_E')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        CALL calc_energy_force(env_id,calc_force=.FALSE.,ierr=ierr)
        CALL my_assert(ierr==0,cmd//' calc_energy_force failed',failure)
        CALL get_energy(env_id,e_pot,ierr)
        CALL my_assert(ierr==0,cmd//' failed get_energy',failure)
        e_pot=e_pot*e_fact
        IF (para_env%mepos==para_env%source) THEN
           WRITE (sout,'(ES22.13)',iostat=iostat) e_pot
           IF (iostat/=0) CALL mp_abort(cmd//' write e_pot')
           CALL m_flush(sout)
        END IF
     CASE('EVALEF','EVAL_EF')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        IF (failure) GOTO 10
        CALL calc_energy_force(env_id,calc_force=.TRUE.,ierr=ierr)
        CALL my_assert(ierr==0,cmd//' calc_energy_force failed',failure)
     CASE('GETF','GET_F')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        CALL get_natom(env_id, n_atom, ierr)
        CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
        IF (failure) GOTO 10
        ALLOCATE(pos(3*n_atom),stat=stat)
        IF (stat/=0) CALL mp_abort(cmd//' failed alloc')
        CALL get_force(env_id, frc=pos, n_el=3*n_atom, ierr=ierr)
        CALL my_assert(ierr==0,'get_force error',failure)
        pos=pos*(e_fact/pos_fact)
        IF (.not.failure.and.para_env%mepos==para_env%source) THEN
           WRITE (sout,'(i10)',iostat=iostat) 3*n_atom
           IF (iostat/=0) CALL mp_abort(cmd//' write nforce')
           WRITE (sout,'(3ES22.13)',iostat=iostat) pos
           IF (iostat/=0) CALL mp_abort(cmd//' write force')
           WRITE (sout,'("* END")',iostat=iostat)
           IF (iostat/=0) CALL mp_abort(cmd//' write coord')
           CALL m_flush(sout)
        END IF
        DEALLOCATE(pos,stat=stat)
        IF (stat/=0) CALL mp_abort(cmd//' failed dealloc')
     CASE('CALCEF','CALC_EF')
        env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
        CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
        IF (failure) GOTO 10
        CALL calc_energy_force(env_id,calc_force=.TRUE.,ierr=ierr)
        CALL my_assert(ierr==0,cmd//' calc_energy_force failed',failure)
        IF (failure) GOTO 10
        CALL get_energy(env_id,e_pot,ierr)
        CALL my_assert(ierr==0,cmd//' failed get_energy',failure)
        IF (failure) GOTO 10
        e_pot=e_pot*e_fact
        CALL get_natom(env_id, n_atom, ierr)
        CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
        IF (failure) GOTO 10
        ALLOCATE(pos(3*n_atom),stat=stat)
        IF (stat/=0) CALL mp_abort(cmd//' failed alloc')
        CALL get_force(env_id, frc=pos, n_el=3*n_atom, ierr=ierr)
        CALL my_assert(ierr==0,'get_force error',failure)
        pos=pos*(e_fact/pos_fact)
        IF (.not.failure.and.para_env%mepos==para_env%source) THEN
           WRITE (sout,'(ES22.13)',iostat=iostat) e_pot
           IF (iostat/=0) CALL mp_abort(cmd//' write e_pot')
           WRITE (sout,'(i10)',iostat=iostat) 3*n_atom
           IF (iostat/=0) CALL mp_abort(cmd//' write nforce')
           WRITE (sout,'(3ES22.13)',iostat=iostat) pos
           IF (iostat/=0) CALL mp_abort(cmd//' write force')
           WRITE (sout,'("* END")',iostat=iostat)
           IF (iostat/=0) CALL mp_abort(cmd//' write coord')
           CALL m_flush(sout)
        END IF
        DEALLOCATE(pos,stat=stat)
        IF (stat/=0) CALL mp_abort(cmd//' failed dealloc')
     CASE('UNITS_EVA','UNITS_EV_A')
        e_fact=evolt
        pos_fact=angstrom
        units='eV_A'
     CASE('UNITS_AU')
        e_fact=1.0_dp
        pos_fact=1.0_dp
        units='au'
     CASE('UNITS')
        IF (para_env%source==para_env%mepos) THEN
           WRITE (sout,'(a)',iostat=iostat) TRIM(units)
           IF (iostat/=0) CALL mp_abort(cmd//' writing units')
           CALL m_flush(sout)
        END IF
     CASE('HARSH')
        harsh=.TRUE.
     CASE('PERMISSIVE')
        harsh=.FALSE.
     CASE('CD','CHDIR')
        shift=shift2+1
        CALL my_assert(shift<=LEN(cmdStr),'cd without directory not accepted',failure)
        IF (failure) GOTO 10
        DO
           IF (cmdStr(shift:shift)/=' ') EXIT
           shift=shift+1
           IF (shift>=LEN(cmdStr)) EXIT
        END DO
        CALL my_assert(cmdStr(shift:shift)/=' ',&
             'cd without directory not accepted2',failure)
        CALL m_chdir(cmdStr(shift:LEN_TRIM(cmdStr)),ierr)
        CALL my_assert(ierr==0,'error changing directory',failure)
     CASE('PWD','CWD')
        IF (para_env%source==para_env%mepos) THEN
           CALL m_getcwd(cmdStr)
           WRITE (sout,'(a)')TRIM(cmdStr)
        END IF
     CASE default
        CALL my_assert(.FALSE.,'unknown command '//cmd,failure)
     END SELECT
10   IF (failure) THEN
        IF (harsh) CALL mp_abort('failure')
        IF (para_env%mepos==para_env%source) THEN
           WRITE(sout,'("* ERROR, type HELP for help")',iostat=iostat)
           IF (iostat/=0) CALL mp_abort('writing error')
           CALL m_flush(sout)
        END IF
        failure=.FALSE.
     END IF
  END DO

  CALL section_release(input_declaration,error=error)

  CALL finalize_cp2k(finalize_mpi=.TRUE.,ierr=ierr)
  IF (ierr/=0) CALL mp_abort('finalize_cp2k failure')

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param tst ...
!> \param msg ...
!> \param failure ...
! *****************************************************************************
  SUBROUTINE my_assert(tst,msg,failure)
    LOGICAL, INTENT(in)                      :: tst
    CHARACTER(LEN=*), INTENT(in)             :: msg
    LOGICAL, INTENT(inout)                   :: failure

    IF (.not.tst) THEN
       IF (harsh) CALL mp_abort(msg)
       IF (para_env%mepos==para_env%source) THEN
          WRITE(sout,'("* ERROR ",a)') msg
       END IF
       failure=.TRUE.
    END IF
  END SUBROUTINE my_assert

! *****************************************************************************
!> \brief ...
!> \param str ...
!> \param startI ...
!> \param default_val ...
!> \retval res ...
! *****************************************************************************
  FUNCTION parse_env_id(str,startI,default_val) RESULT(res)
    CHARACTER(LEN=*), INTENT(in)             :: str
    INTEGER, INTENT(in)                      :: startI, default_val
    INTEGER                                  :: res

    INTEGER                                  :: iostat, pos

    res=default_val
    pos=startI
    IF (pos<=LEN(str)) THEN
       DO
          IF (str(pos:pos)/=' ') EXIT
          pos=pos+1
          IF (pos>=LEN(str)) EXIT
       END DO
       IF (str(pos:pos)/=' ') THEN
          READ(str(pos:LEN(str)),*,iostat=iostat) res
          IF (iostat/=0) res=-1
       END IF
    END IF
  END FUNCTION parse_env_id

! *****************************************************************************
!> \brief ...
!> \param str ...
!> \retval res ...
! *****************************************************************************
  FUNCTION ll(str) RESULT(res)
    CHARACTER(LEN=*), INTENT(in)             :: str
    CHARACTER(LEN=default_string_length)     :: res

    res=str
  END FUNCTION

! *****************************************************************************
!> \brief ...
!> \param unit ...
! *****************************************************************************
 SUBROUTINE print_v_info(unit)
    INTEGER, INTENT(in)                      :: unit

    CHARACTER(LEN=default_path_length)       :: cwd
    CHARACTER(LEN=default_string_length)     :: dev_flag, host_name, user_name

    CALL m_getcwd(cwd)
    CALL m_getpid(pid)
    CALL m_getlog(user_name)
    CALL m_hostnm(host_name)
    dev_flag = ""
    IF (id_cp2k_version==id_development_version) dev_flag=" (Development Version)"
    WRITE (UNIT=unit,FMT="(A,A)")&
         " PROGRAM STARTED ON ",TRIM(host_name)
    WRITE (UNIT=unit,FMT="(A,A)")&
         " PROGRAM STARTED BY ",TRIM(user_name)
    WRITE (UNIT=unit,FMT="(A,i10)")&
         " PROGRAM PROCESS ID ",pid
    WRITE (UNIT=unit,FMT="(A,A)")&
         " PROGRAM STARTED IN ",TRIM(cwd)
    WRITE (UNIT=unit,FMT="(/,T2,A,T31,A50)")&
         "CP2K| version string: ",&
         ADJUSTR(TRIM(cp2k_version)//TRIM(dev_flag))
    WRITE (UNIT=unit,FMT="(T2,A,T41,A40)")&
         "CP2K| source code revision number:",&
         ADJUSTR(compile_revision)
    WRITE (UNIT=unit,FMT="(T2,A,T41,A40)")&
         "CP2K| is freely available from ",&
         ADJUSTR(TRIM(cp2k_home))
    WRITE (UNIT=unit,FMT="(T2,A,T31,A50)")&
         "CP2K| Program compiled at",&
         ADJUSTR(compile_date(1:MIN(50,LEN(compile_date))))
    WRITE (UNIT=unit,FMT="(T2,A,T31,A50)")&
         "CP2K| Program compiled on",&
         ADJUSTR(compile_host(1:MIN(50,LEN(compile_host))))
    WRITE (UNIT=unit,FMT="(T2,A,T31,A50)")&
         "CP2K| Program compiled for",&
         ADJUSTR(compile_arch(1:MIN(50,LEN(compile_arch))))
 END SUBROUTINE

! *****************************************************************************
!> \brief returns the number returned by iostat on eof
!> \retval res ...
! *****************************************************************************
  FUNCTION get_eof_stat() RESULT(res)
    INTEGER                                  :: res

    CHARACTER(LEN=5)                         :: rdBuf
    INTEGER                                  :: stat_n, unit

    IF (.NOT.did_init) THEN
       unit=get_unit_number()
       OPEN (UNIT=unit, FILE='.tmp_eof_code',FORM='FORMATTED', STATUS='NEW')

       WRITE(unit,fmt=*)'z'
       REWIND(unit)

       eof_stat=0
       stat_n=0
       DO
          READ (unit,'(a1)', iostat=stat_n)rdBuf
          IF (stat_n/=0) EXIT
       END DO
       eof_stat=stat_n
       IF (eof_stat >= 0) THEN
          STOP 'ERROR,cp_files:get_eof_stat failed identification of EOF IOSTAT'
          eof_stat=-1
       END IF
       CALL close_file(unit,'delete')
       did_init=.TRUE.
    END IF
    res=eof_stat
  END FUNCTION get_eof_stat

END PROGRAM cp2k_shell
