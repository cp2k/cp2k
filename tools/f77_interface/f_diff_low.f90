!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/f_diff_low [1.0] *
!!
!!   NAME
!!     f_diff_low
!!
!!   FUNCTION
!!     This an example of a program that uses the f77 interface.
!!     This program reads the file "input.inp" that should a valid
!!     cp2k input file for something that calculates the force 
!!     (MD for example).
!!     Then this program calculates the forces with finite differences and
!!     compares them to the analytic ones. The error and the percentual error
!!     with respect to the mean quadratic force are printed.
!!
!!   NOTES
!!     You should be able to link it against the cp2k library
!!     Something like
!!       cd ../../makefiles
!!       make libs VERSION=sopt
!!       cd ../tools/f77_interface
!!       f95 f_diff_low.f90 -L${HOME}/cp2k/lib/Linux-i686-nag/sdbg \
!!           -lcp2k_lib -lcp2k_base_lib -L${HOME}/lib -llapack -lg2c
!!     should work
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [fawzi]
!!
!!*** ***********************************************************************
PROGRAM f_diff_low
  IMPLICIT NONE
  INTEGER :: ierr, f_env_id, natom, i
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
  REAL(dp), DIMENSION(:), POINTER :: pos,f
  REAL(dp) :: e0,ep,em,p0,f_m,err_max,err
  REAL(dp), PARAMETER :: small=5.e-3

  CALL cp_init_cp2k(1,ierr)
  IF (ierr/=0) STOP "init_cp2k"
  
  CALL cp_create_fenv(f_env_id,"input.inp","out.out",ierr)
  IF (ierr/=0) STOP "create_force_env"
  CALL cp_get_natom(f_env_id,natom,ierr)
  IF (ierr/=0) STOP "get_natom"
  
  ALLOCATE(pos(natom*3),f(natom*3),stat=ierr)
  IF (ierr/=0) STOP "alloc"
  
  CALL cp_calc_energy_force(f_env_id,1,ierr)
  IF (ierr/=0) STOP "calc_force1"
  CALL cp_get_energy(f_env_id,e0,ierr)
  IF (ierr/=0) STOP "get_energy1"
  WRITE(*,*) "e0= ",e0
  CALL cp_get_force(f_env_id,f,SIZE(f),ierr)
  IF (ierr/=0) STOP "get_force"
  CALL cp_get_pos(f_env_id,pos,SIZE(pos),ierr)
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
     CALL cp_set_pos(f_env_id,pos,SIZE(pos),ierr)
     IF (ierr/=0) STOP "set_pos1"
     CALL cp_calc_energy_force(f_env_id,1,ierr)
     IF (ierr/=0) STOP "calc_energy_force2"
     CALL cp_get_energy(f_env_id,ep,ierr)
     IF (ierr/=0) STOP "get_energy2"
     pos(i)=p0-small
     CALL cp_calc_energy(f_env_id,pos,SIZE(pos),em,ierr)
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
!FM  CALL cp_destroy_force_env(f_env_id,ierr)
!FM  IF (ierr/=0) STOP "destroy_force_env"

  CALL cp_finalize_cp2k(1,ierr)
  IF (ierr/=0) STOP "finalize_cp2k"

END PROGRAM f_diff_low
