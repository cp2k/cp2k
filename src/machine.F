!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000,2004  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/machine [1.0] *
!!
!!   NAME
!!     machine
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     APSI & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (05.07.2001) : added G95 interface
!!     - m_flush added (12.06.2002,MK)
!!     - Missing print_memory added (24.09.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE machine

#if defined ( __AIX )
  USE machine_aix, ONLY    : & 
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __IRIX )
  USE machine_irix, ONLY   : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __G95 )
  USE machine_g95, ONLY    : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __INTEL )
  USE machine_intel, ONLY  : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __PGI )
  USE machine_pgi, ONLY    : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __T3E )
  USE machine_t3e, ONLY    : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __DEC )
  USE machine_dec, ONLY    : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __ABSOFT )
  USE machine_absoft, ONLY : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __SUN )
  USE machine_sun, ONLY    : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __NAG )
  USE machine_nag, ONLY    : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#elif defined ( __SX )
  USE machine_sx, ONLY    : &
       m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: m_walltime, m_cputime, m_datum, m_flush, &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg, &
       print_memory, m_iargc, m_abort, m_chdir, m_loc_r, m_loc_c

END MODULE machine
!!*****
