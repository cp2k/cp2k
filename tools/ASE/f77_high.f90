PROGRAM f77_high
  USE f77_interface,     ONLY : default_para_env,&
       init_cp2k, finalize_cp2k,&
       create_force_env, destroy_force_env, set_pos, get_pos,&
       get_force, calc_energy_force, run_input, get_energy, &
       calc_energy, calc_force, check_input, get_natom
  IMPLICIT NONE
  INTEGER :: ierr, f_env_id, natom, i
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
  REAL(dp), DIMENSION(:), POINTER :: pos,f
  REAL(dp) :: e0,ep,em,p0,f_m,err_max,err
  REAL(dp), PARAMETER :: small=5.e-3
  CALL init_cp2k(.TRUE.,ierr)
  
  CALL create_force_env(f_env_id,"input.inp","out.out",ierr=ierr)
  CALL get_natom(f_env_id,natom,ierr)
  
  ALLOCATE(pos(natom*3),f(natom*3),stat=ierr)
  
  CALL calc_energy_force(f_env_id,calc_force=.TRUE.,ierr=ierr)
  CALL get_energy(f_env_id,e0,ierr)
  WRITE(*,*) "e0= ",e0
  CALL get_force(f_env_id,f,SIZE(f),ierr)
  CALL get_pos(f_env_id,pos,SIZE(pos),ierr)
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
     CALL calc_energy_force(f_env_id,calc_force=.FALSE.,ierr=ierr)
     CALL get_energy(f_env_id,ep,ierr)
     pos(i)=p0-small
     CALL set_pos(f_env_id,pos,SIZE(pos),ierr)
     CALL calc_energy_force(f_env_id,calc_force=.FALSE.,ierr=ierr)
     CALL get_energy(f_env_id,em,ierr)
     pos(i)=p0
     err=(em-ep)/(2._dp*small)-f(i)
     WRITE(*,*) "err= ",err
     if (abs(err)>abs(err_max)) err_max=err
     WRITE (*,*) "f= ",f(i)," f_fdiff= ",(em-ep)/(2._dp*small),&
          " err= ",err," err_m= ",err/f_m*100._dp
  END DO

  WRITE (*,*)
  WRITE (*,*) "err_max ",err_max," err_max_m ",err_max/f_m*100._dp

END PROGRAM f77_high
