!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/f_diff_high [1.0] *
!!
!!   NAME
!!     f_diff_high
!!
!!   FUNCTION
!!     This an example of a program that uses the f77_interface module.
!!     This program reads the file "input.inp" that should a valid
!!     cp2k input file for something that calculates the force 
!!     (MD for example).
!!     Then this program calculates the forces with finite differences and
!!     compares them to the analytic ones. The error and the percentual error
!!     with respect to the mean quadratic force are printed.
!!
!!   NOTES
!!     The easiest way to use it is to replace cp2k/src/cp2k.F with this
!!     file and recompile cp2k.
!!
!!     Using this method you could also get full access to cp2k
!!     by making the subroutines f77_interface:f_env_add_defaults
!!     and f77_interface:f_env_rm_defaults public.
!!     Very useful to mess around, but obviously it should be done with
!!     *extreme* care, that is the reason this two routines are not public
!!     at the moment.
!!     If you need it in production code think about extending the
!!     f77_interface.
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [fawzi]
!!
!!*** ***********************************************************************
PROGRAM f_diff_high
  USE f77_interface, ONLY: init_cp2k,create_force_env,get_natom,calc_energy_force,&
       get_energy,get_force,get_pos,set_pos,finalize_cp2k, destroy_force_env,&
       calc_energy
  IMPLICIT NONE
  INTEGER :: ierr, f_env_id, natom, i
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
  REAL(dp), DIMENSION(:), POINTER :: pos,f
  REAL(dp) :: e0,ep,em,p0,f_m,err_max,err
  REAL(dp), PARAMETER :: small=5.e-3

  CALL init_cp2k(.TRUE.,ierr)
  IF (ierr/=0) STOP "init_cp2k"
  
  CALL create_force_env(f_env_id,"input.inp","out.out",ierr=ierr)
  IF (ierr/=0) STOP "create_force_env"
  CALL get_natom(f_env_id,natom,ierr)
  IF (ierr/=0) STOP "get_natom"
  
  ALLOCATE(pos(natom*3),f(natom*3),stat=ierr)
  IF (ierr/=0) STOP "alloc"
  
  CALL calc_energy_force(f_env_id,calc_force=.TRUE.,ierr=ierr)
  IF (ierr/=0) STOP "calc_force1"
  CALL get_energy(f_env_id,e0,ierr)
  IF (ierr/=0) STOP "get_energy1"
  WRITE(*,*) "e0= ",e0
  CALL get_force(f_env_id,f,SIZE(f),ierr)
  IF (ierr/=0) STOP "get_force"
  CALL get_pos(f_env_id,pos,SIZE(pos),ierr)
  IF (ierr/=0) STOP "get_pos"
  f_m=0._dp
  DO i=1,SIZE(f)
     f_m=f_m+f(i)**2
  END DO
  f_m=SQRT(f_m/natom)
  WRITE(*,*)"f_m",f_m
  WRITE(*,*)
  err_max=0._dp
  DO i=1,SIZE(pos)
     p0=pos(i)
     pos(i)=p0+small
     CALL set_pos(f_env_id,pos,SIZE(pos),ierr)
     IF (ierr/=0) STOP "set_pos1"
     CALL calc_energy_force(f_env_id,calc_force=.FALSE.,ierr=ierr)
     IF (ierr/=0) STOP "calc_energy_force2"
     CALL get_energy(f_env_id,ep,ierr)
     IF (ierr/=0) STOP "get_energy2"
     pos(i)=p0-small
     CALL calc_energy(f_env_id,pos,SIZE(pos),em,ierr)
     IF (ierr/=0) STOP "calc_energy"
     pos(i)=p0
     err=(em-ep)/(2._dp*small)-f(i)
     WRITE(*,*) "err= ",err
     if (abs(err)>abs(err_max)) err_max=err
     WRITE (*,*) "f= ",f(i)," f_fdiff= ",(em-ep)/(2._dp*small),&
          " err= ",err," err_m= ",err/f_m*100._dp
  END DO

  WRITE (*,*)
  WRITE (*,*) "err_max ",err_max," err_max_m ",err_max/f_m*100._dp
!FM  CALL destroy_force_env(f_env_id,ierr)
!FM  IF (ierr/=0) STOP "destroy_force_env"

  CALL finalize_cp2k(.TRUE.,ierr)
  IF (ierr/=0) STOP "finalize_cp2k"

END PROGRAM f_diff_high
