!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE machine

#if defined ( __AIX )
  USE machine_aix, ONLY : m_walltime, m_cputime, m_datum,  &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
#elif defined ( __IRIX )
  USE machine_irix, ONLY : m_walltime, m_cputime, m_datum,  &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
#elif defined ( __PGI )
  USE machine_pgi, ONLY : m_walltime, m_cputime, m_datum,  &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
#elif defined ( __T3E )
  USE machine_t3e, ONLY : m_walltime, m_cputime, m_datum,  &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
#elif defined ( __DEC )
  USE machine_dec, ONLY : m_walltime, m_cputime, m_datum,  &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
#elif defined ( __ABSOFT )
  USE machine_absoft, ONLY : m_walltime, m_cputime, m_datum,  &
       m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: m_walltime, m_cputime, m_datum
  PUBLIC :: m_hostnm, m_getcwd, m_getlog, m_getuid, m_getpid, m_getarg
  
END MODULE machine
