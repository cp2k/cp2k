!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/f77_int_low *
!!
!!   NAME
!!     f77_int_low
!!
!!   FUNCTION
!!     interface to use cp2k as library (these are f77 wrappers to the 
!!     f77_interface functions, for a description see in that module)
!!
!!   NOTES
!!     I tried to use just integer, real(db), and CHARACTER(len=*) to
!!     make the c binding/wrapping as easy as possible.
!!     The absence of intent from the arguments is also done on purpose.
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [fawzi]
!!     11.2004 parallel version [fawzi]
!!
!!   SOURCE
!****************************************************************************
SUBROUTINE cp_init_cp2k(init_mpi,ierr)
  USE f77_interface, ONLY: icp => init_cp2k
  IMPLICIT NONE
  INTEGER :: init_mpi,ierr
  LOGICAL ::i_mpi
  i_mpi=init_mpi/=0
  CALL icp(i_mpi,ierr)
END SUBROUTINE cp_init_cp2k

SUBROUTINE cp_finalize_cp2k(finalize_mpi,ierr)
  USE f77_interface, ONLY: kcp => finalize_cp2k
  IMPLICIT NONE
  INTEGER :: finalize_mpi, ierr
  LOGICAL ::f_mpi
  f_mpi=finalize_mpi/=0
  CALL kcp(f_mpi,ierr)
END SUBROUTINE cp_finalize_cp2k

SUBROUTINE cp_create_fenv(new_env_id,input_file_path,output_file_path,&
     ierr)
  USE f77_interface, ONLY: cfe => create_force_env
  IMPLICIT NONE
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: new_env_id, ierr

  CALL cfe(new_env_id,input_file_path,output_file_path,ierr=ierr)
END SUBROUTINE cp_create_fenv

SUBROUTINE cp_create_fenv_comm(new_env_id,input_file_path,output_file_path,&
     mpi_comm,ierr)
  USE f77_interface, ONLY: cfe => create_force_env
  IMPLICIT NONE
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: new_env_id, mpi_comm, ierr

  CALL cfe(new_env_id,input_file_path,output_file_path,ierr,mpi_comm)
END SUBROUTINE cp_create_fenv_comm

SUBROUTINE cp_destroy_fenv(env_id, ierr)
  USE f77_interface, ONLY: dfe => destroy_force_env
  IMPLICIT NONE
  INTEGER :: env_id, ierr

  CALL dfe(env_id,ierr)
END SUBROUTINE cp_destroy_fenv

SUBROUTINE cp_set_pos(env_id, new_pos, n_el, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: sp => set_pos
  IMPLICIT NONE
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: new_pos

  CALL sp(env_id,new_pos,n_el,ierr)
END SUBROUTINE cp_set_pos

SUBROUTINE cp_get_natom(env_id, natom, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: gna => get_natom
  IMPLICIT NONE
  INTEGER :: natom, env_id, ierr

  CALL gna(env_id,natom,ierr)
END SUBROUTINE cp_get_natom

SUBROUTINE cp_get_pos(env_id, pos, n_el, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: gp => get_pos
  IMPLICIT NONE
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: pos

  CALL gp(env_id,pos,n_el,ierr)
END SUBROUTINE cp_get_pos

SUBROUTINE cp_get_force(env_id,force, n_el, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: gf => get_force
  IMPLICIT NONE
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: force

  CALL gf(env_id,force,n_el,ierr)
END SUBROUTINE cp_get_force

SUBROUTINE cp_get_energy(env_id, e_pot, ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: ge => get_energy
  IMPLICIT NONE
  INTEGER :: env_id, ierr
  REAL(kind=dp) :: e_pot

  CALL ge(env_id,e_pot,ierr)
END SUBROUTINE cp_get_energy

SUBROUTINE cp_calc_energy_force(env_id,calc_force,ierr)
  USE f77_interface, ONLY: cef => calc_energy_force
  IMPLICIT NONE
  INTEGER :: calc_force, env_id, ierr
  LOGICAL :: cf

  cf=(calc_force/=0)
  CALL cef(env_id,cf,ierr)
END SUBROUTINE cp_calc_energy_force

SUBROUTINE cp_calc_energy(env_id,pos,n_el,e_pot,ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: ce => calc_energy
  IMPLICIT NONE
  INTEGER :: calc_force, env_id, ierr, n_el
  REAL(dp) :: e_pot
  REAL(dp), DIMENSION(1:n_el) :: pos

  CALL ce(env_id,pos,n_el,e_pot,ierr)
END SUBROUTINE cp_calc_energy

SUBROUTINE cp_calc_force(env_id,pos,n_el_pos,e_pot,force,n_el_force,ierr)
  USE kinds, ONLY: dp
  USE f77_interface, ONLY: cf => calc_force
  IMPLICIT NONE
  INTEGER :: env_id, ierr, n_el_pos, n_el_force
  REAL(dp) :: e_pot
  REAL(dp), DIMENSION(1:n_el_pos) :: pos
  REAL(dp), DIMENSION(1:n_el_force) :: force

  CALL cf(env_id,pos,n_el,e_pot,force,n_el_force,ierr)
END SUBROUTINE cp_calc_force

SUBROUTINE cp_run_input(input_file_path,output_file_path,ierr)
  USE cp2k_runs, ONLY: rin => run_input
  IMPLICIT NONE
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: ierr

  CALL rin(input_file_path,output_file_path,ierr=ierr)
END SUBROUTINE cp_run_input

SUBROUTINE cp_run_input_comm(input_file_path,output_file_path,&
     mpi_comm,ierr)
  USE f77_interface, ONLY: rin => run_input
  IMPLICIT NONE
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: mpi_comm, ierr

  CALL rin(input_file_path,output_file_path,ierr,mpi_comm)
END SUBROUTINE cp_run_input_comm
