! sample of a simple runner that uses the f77_interface
! it can beused to connect c programs, communicating through std-in/ std-out
!
! positions are in angstrom, energies in evolt
!
! commands:
! load filename: loads the filename, returns the env_id, or -1 in case of error
! natom [env_id]: returns the number of atoms in the environement env_id
!  (defaults to the last loaded)
! setpos [env_id]: sets the positions of the atoms, should be followed 
!  by natom*3 (on a line) and then all the positions [angstrom]
! getpos [env_id]: gets the positions of the atoms, returns
!  natom*3 (on a line) and then all the positions [angstrom]
! calcE [env_id]: calculate the energy and returns it (in eV)
! calcEF [env_id]: calculate the energy and forces and returns it,
!  first the energy on a line (in eV), then the natom*3 (on a line)
!  and finally all the values (in eV/angstrom)
!
! fawzi
PROGRAM cp2k_runner
USE kinds, only: dp
use f77_interface
USE cp_para_types, only: cp_para_env_type
USE physcon, ONLY: angstrom,evolt
USE message_passing, ONLY: mp_bcast,mp_abort
USE cp_files, only: get_eof_stat
USE string_utilities, only: uppercase,compress
use machine, only: m_flush  
#include "cp_common_uses.h"
implicit none

INTEGER :: ierr,i,iostat,shift,shift2,env_id,last_env_id,n_atom,stat,n_atom2
INTEGER, parameter :: sout=6 ! ,sin=5
type(cp_error_type) :: error
TYPE(cp_para_env_type), pointer :: para_env
TYPE(cp_logger_type), pointer :: logger
character(len=80) :: cmdStr
character(len=20) :: cmd,units
logical :: eof
REAL(dp), DIMENSION(:), allocatable :: pos
real(dp) :: e_pot
REAL(dp) :: e_fact,pos_fact
LOGICAL :: harsh,failure

call init_cp2k(init_mpi=.true.,ierr=ierr)
IF (ierr/=0) CALL mp_abort('init_cp2k failure')

pos_fact=1.0_dp
e_fact=1.0_dp
units='au'
harsh=.false.
failure=.false.

call cp_error_init(error)
logger => cp_error_get_logger(error)
para_env => logger%para_env

eof=.false.
last_env_id=-1
do
   WRITE (sout,'("* READY")')
   call m_flush(sout)
   IF (para_env%mepos==para_env%source) THEN
      READ(*,'(a)',iostat=iostat)cmdStr
      IF (iostat/=0) THEN
         IF (iostat==get_eof_stat()) THEN
            write(sout,'(a)',iostat=iostat) '* EOF'
            ! ignoring iostat on purpose
         END IF
         eof=.true.
      END IF
   END IF
   CALL mp_bcast(eof,para_env%source,para_env%group)
   if (eof) exit
   CALL mp_bcast(cmdStr,para_env%source,para_env%group)
   shift=1
   DO
      IF(cmdStr(shift:shift)/=' ') EXIT
      shift=shift+1
      if (shift>=len(cmdStr)) EXIT
   END DO
   shift2=shift
   DO
      IF(cmdStr(shift2:shift2)==' ') EXIT
      shift2=shift2+1
      IF (shift2>=LEN(cmdStr)) EXIT
   END DO
   shift2=shift2-1
   CALL my_assert(shift<=shift2,'* No command',failure)
   if (failure) goto 10
   cmd=cmdStr(shift:shift2)
   call uppercase(cmd)
   SELECT CASE(cmd)
   CASE('HELP')
    write(sout,'(a)',iostat=iostat) 'Commands'
    if (iostat/=0) call mp_abort('writing help')
    write(sout,'(a)',iostat=iostat) &
    ' (if there is [env_id] it means that an optional env_id can be given,'
    if (iostat/=0) call mp_abort('writing help')
    write(sout,'(a)',iostat=iostat) &
    ' if none is given it defaults to the last env_id loaded)'
    if (iostat/=0) call mp_abort('writing help')
    write(sout,'(a)',iostat=iostat) ' '
    if (iostat/=0) call mp_abort('writing help')
    write(sout,'(a)',iostat=iostat) &
   ' load filename: loads the filename, returns the env_id, or -1 in case of error'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' natom [env_id]: returns the number of atoms in the environement env_id'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' setpos [env_id]: sets the positions of the atoms, should be followed'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' by natom*3 (on a line) and then all the positions'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' getpos [env_id]: gets the positions of the atoms, returns'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   '  natom*3 (on a line) and then all the positions then "* END" (alone on a line)'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' calcE [env_id]: calculate the energy and returns it'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' calcEF [env_id]: calculate the energy and forces and returns them,'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   '  first the energy on a line, then the natom*3 (on a line)'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   '  and finally all the values and "* END" (alone on a line)'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' harsh: stops on any error'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' permissive: stops only on serious errors'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' units: returns the units used for energy and position'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' units_ev_a: sets the units to electron volt (energy)  and Angstrom (positions)'
   if (iostat/=0) call mp_abort('writing help')
   write(sout,'(a)',iostat=iostat) &
   ' help: writes the present help'
   if (iostat/=0) call mp_abort('writing help')
   call m_flush(sout)
   CASE('LOAD')
      shift=shift2+1
      CALL my_assert(shift<=LEN(cmdStr),'load without file',failure)
      IF (failure) GOTO 10
      DO
         IF (cmdStr(shift:shift)/=' ') exit
         shift=shift+1
         IF (shift>=len(cmdStr)) EXIT
      END DO
      CALL my_assert(cmdStr(shift:shift)/=' ','load without file2',failure)
      IF (failure) GOTO 10
      CALL create_force_env(new_env_id=last_env_id,&
            input_path=cmdStr(shift:LEN_TRIM(cmdStr)),&
            output_path=cmdStr(shift:LEN_TRIM(cmdStr))//'.out',&
            owns_out_unit=.true.,ierr=ierr)
      IF (ierr/=0) last_env_id=-1
      IF (para_env%source==para_env%mepos) THEN
         WRITE (sout,'(i10)',iostat=iostat) last_env_id
         IF (iostat/=0) CALL mp_abort(cmd//' failed writing last_env_id')
         call m_flush(sout)
      END IF
      if (last_env_id<=0) failure=.true.
   CASE('NATOM','N_ATOM')
      env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
      CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
      if (failure) goto 10
      CALL get_natom(env_id, n_atom, ierr)
      CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
      IF (.NOT. failure .and. para_env%mepos==para_env%source) THEN
         WRITE (sout,'(i10)',iostat=iostat) n_atom
         call m_flush(sout)
      END IF
   CASE('SETPOS','SET_POS')
      env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
      CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
      if (failure) goto 10
      CALL get_natom(env_id, n_atom, ierr)
      CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
      if (failure) goto 10
      ALLOCATE(pos(3*n_atom),stat=stat)
      CALL my_assert(stat==0,cmd//' failed alloc',failure)
      IF (para_env%mepos==para_env%source) THEN
         READ (*,*,iostat=iostat) n_atom2
         IF (iostat/=0) CALL mp_abort('setpos read n_atom')
         IF (n_atom2/=SIZE(pos)) THEN
            CALL my_assert(.false.,'setpos invalid number of atoms',failure)
            DO i=1,n_atom
               READ(*,'(a)',iostat=iostat) cmdStr
               CALL compress(cmdStr,full=.true.)
               call uppercase(cmdStr)
               if (cmdStr=='*END') EXIT
            end do
            goto 10
         END IF
         READ (*,*,iostat=iostat) pos
         IF (iostat/=0) CALL mp_abort('setpos read coord')
         pos=pos/pos_fact
         READ(*,'(a)',iostat=iostat) cmdStr
         CALL compress(cmdStr,full=.true.)
         call uppercase(cmdStr)
         CALL my_assert(cmdStr=='*END',' missing *END',failure)
      END IF
      CALL mp_bcast(pos,para_env%source,para_env%group)
      CALL set_pos(env_id, new_pos=pos, n_el=3*n_atom, ierr=ierr)
      CALL my_assert(ierr==0,'set_pos error',failure)
      DEALLOCATE(pos,stat=stat)
      IF (stat/=0) CALL mp_abort('setpos failed dealloc')      
   CASE('GETPOS','GET_POS')
      env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
      CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
      if (failure) goto 10
      CALL get_natom(env_id, n_atom, ierr)
      CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
      if (failure) goto 10
      ALLOCATE(pos(3*n_atom),stat=stat)
      IF (stat/=0) CALL mp_abort(cmd//' failed alloc')
      CALL get_pos(env_id, pos=pos, n_el=3*n_atom, ierr=ierr)
      CALL my_assert(ierr==0,'get_pos error',failure)
      pos=pos*pos_fact
      IF (.not.failure.and.para_env%mepos==para_env%source) THEN
         WRITE (sout,'(i10)',iostat=iostat) 3*n_atom
         IF (iostat/=0) CALL mp_abort(cmd//' write ncoord')
         WRITE (sout,'(3ES24.15)',iostat=iostat) pos
         IF (iostat/=0) CALL mp_abort(cmd//' write coord')
         WRITE (sout,'(a)',iostat=iostat) "* END"
         IF (iostat/=0) CALL mp_abort(cmd//' write coord')
         call m_flush(sout)
      END IF
      DEALLOCATE(pos,stat=stat)
      IF (stat/=0) CALL mp_abort(cmd//' failed dealloc')      
   CASE('CALCE','CALC_E')
      env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
      CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
      CALL calc_energy_force(env_id,calc_force=.false.,ierr=ierr)
      CALL my_assert(ierr==0,cmd//' calc_energy_force failed',failure)
      CALL get_energy(env_id,e_pot,ierr)
      CALL my_assert(ierr==0,cmd//' failed get_energy',failure)
      e_pot=e_pot*e_fact
      IF (para_env%mepos==para_env%source) THEN
         WRITE (sout,'(ES24.15)',iostat=iostat) e_pot
         IF (iostat/=0) CALL mp_abort(cmd//' write e_pot')
         call m_flush(sout)
      END IF
   CASE('CALCEF','CALC_EF')
      env_id=parse_env_id(str=cmdStr,startI=shift2+1,default_val=last_env_id)
      CALL my_assert(env_id>0,cmd//' invalid env_id',failure)
      if (failure) goto 10
      CALL calc_energy_force(env_id,calc_force=.true.,ierr=ierr)
      CALL my_assert(ierr==0,cmd//' calc_energy_force failed',failure)
      if (failure) goto 10
      CALL get_energy(env_id,e_pot,ierr)
      CALL my_assert(ierr==0,cmd//' failed get_energy',failure)
      if (failure) goto 10
      e_pot=e_pot*e_fact
      CALL get_natom(env_id, n_atom, ierr)
      CALL my_assert(ierr==0,cmd//' failed get_natom',failure)
      if (failure) goto 10
      ALLOCATE(pos(3*n_atom),stat=stat)
      IF (stat/=0) CALL mp_abort(cmd//' failed alloc')
      CALL get_force(env_id, force=pos, n_el=3*n_atom, ierr=ierr)
      CALL my_assert(ierr==0,'get_force error',failure)
      pos=pos*(e_fact/pos_fact)
      IF (.not.failure.and.para_env%mepos==para_env%source) THEN
         WRITE (sout,'(ES24.15)',iostat=iostat) e_pot
         IF (iostat/=0) CALL mp_abort(cmd//' write e_pot')
         WRITE (sout,'(i10)',iostat=iostat) 3*n_atom
         IF (iostat/=0) CALL mp_abort(cmd//' write nforce')
         WRITE (sout,'(3ES24.15)',iostat=iostat) pos
         IF (iostat/=0) CALL mp_abort(cmd//' write force')
         WRITE (sout,'("* END")',iostat=iostat)         
         IF (iostat/=0) CALL mp_abort(cmd//' write coord')
         call m_flush(sout)
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
            WRITE (sout,'(a)',iostat=iostat) trim(units)
            IF (iostat/=0) CALL mp_abort(cmd//' writing units')
            call m_flush(sout)
         END IF
      CASE('HARSH')
         harsh=.true.
      CASE('PERMISSIVE')
         harsh=.false.
   CASE default
      CALL my_assert(.false.,'unknown command '//cmd,failure)
   END SELECT
10 IF (failure) THEN
      if (harsh) call mp_abort('failure')
      if (para_env%mepos==para_env%source) then
          WRITE(sout,'("* ERROR, type HELP for help")',iostat=iostat)
          if (iostat/=0) call mp_abort('writing error')
          call m_flush(sout)
      end if
      failure=.false.
   END IF
end do
call finalize_cp2k(finalize_mpi=.true.,ierr=ierr)
IF (ierr/=0) CALL mp_abort('finalize_cp2k failure')

contains
 
subroutine my_assert(tst,msg,failure)
    logical, intent(in) :: tst
    character(len=*), intent(in) :: msg
    logical, intent(inout) :: failure
    
    if (.not.tst) then
        if (harsh) call mp_abort(msg)
        if (para_env%mepos==para_env%source) then
            write(sout,'("* ERROR ",a)') msg
        end if
        failure=.true.
    end if
end subroutine

FUNCTION parse_env_id(str,startI,default_val) RESULT(res)
  CHARACTER(len=*), intent(in) :: str
  INTEGER, intent(in) :: startI,default_val
  INTEGER :: res
  
  INTEGER :: pos,iostat

  res=default_val
  pos=startI
  IF (pos<=len(str)) THEN
     DO
        IF (str(pos:pos)/=' ') EXIT
        pos=pos+1
        IF (pos>=len(str)) EXIT
     END DO
     IF (str(pos:pos)/=' ') THEN
        READ(str(pos:len(str)),*,iostat=iostat) res
        IF (iostat/=0) res=-1
     END IF
  END IF
END FUNCTION parse_env_id

END PROGRAM cp2k_runner
