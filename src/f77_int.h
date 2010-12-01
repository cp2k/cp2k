INTERFACE

SUBROUTINE cp_init_cp2k(init_mpi,ierr)
  INTEGER :: init_mpi,ierr
END SUBROUTINE cp_init_cp2k

SUBROUTINE cp_finalize_cp2k(finalize_mpi,ierr)
  INTEGER :: finalize_mpi, ierr
END SUBROUTINE cp_finalize_cp2k

SUBROUTINE cp_create_fenv(new_env_id,input_file_path,output_file_path,&
     ierr)
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: new_env_id, ierr
END SUBROUTINE cp_create_fenv

SUBROUTINE cp_create_fenv_comm(new_env_id,input_file_path,output_file_path,&
     mpi_comm,ierr)
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: new_env_id, mpi_comm, ierr
END SUBROUTINE cp_create_fenv_comm

RECURSIVE SUBROUTINE cp_destroy_fenv(env_id, ierr)
  INTEGER :: env_id, ierr
END SUBROUTINE cp_destroy_fenv

SUBROUTINE cp_set_pos(env_id, new_pos, n_el, ierr)
  use kinds, only: dp
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: new_pos
END SUBROUTINE cp_set_pos

SUBROUTINE cp_set_vel(env_id, new_vel, n_el, ierr)
  use kinds, only: dp
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: new_vel
END SUBROUTINE cp_set_vel

SUBROUTINE cp_get_natom(env_id, natom, ierr)
  INTEGER :: natom, env_id, ierr
END SUBROUTINE cp_get_natom

SUBROUTINE cp_get_nparticle(env_id, nparticle, ierr)
  INTEGER :: nparticle, env_id, ierr
END SUBROUTINE cp_get_nparticle

SUBROUTINE cp_get_pos(env_id, pos, n_el, ierr)
  use kinds, only: dp
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: pos
END SUBROUTINE cp_get_pos

SUBROUTINE cp_get_force(env_id,force, n_el, ierr)
  use kinds, only: dp
  INTEGER :: n_el, env_id, ierr
  REAL(kind=dp), DIMENSION(1:n_el) :: force
END SUBROUTINE cp_get_force

RECURSIVE SUBROUTINE cp_get_energy(env_id, e_pot, ierr)
  use kinds, only: dp
  INTEGER :: env_id, ierr
  REAL(kind=dp) :: e_pot
END SUBROUTINE cp_get_energy

RECURSIVE SUBROUTINE cp_calc_energy_force(env_id,calc_force,ierr)
  INTEGER :: calc_force, env_id, ierr
END SUBROUTINE cp_calc_energy_force

RECURSIVE SUBROUTINE cp_calc_energy(env_id,pos,n_el,e_pot,ierr)
  use kinds, only: dp
  INTEGER :: env_id, ierr, n_el
  REAL(dp) :: e_pot
  REAL(dp), DIMENSION(1:n_el) :: pos
END SUBROUTINE cp_calc_energy

SUBROUTINE cp_calc_force(env_id,pos,n_el_pos,e_pot,force,n_el_force,ierr)
  use kinds, only: dp
  INTEGER :: env_id, ierr, n_el_pos, n_el_force
  REAL(dp) :: e_pot
  REAL(dp), DIMENSION(1:n_el_pos) :: pos
  REAL(dp), DIMENSION(1:n_el_force) :: force
END SUBROUTINE cp_calc_force

SUBROUTINE cp_run_input(input_file_path,output_file_path,ierr)
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: ierr
END SUBROUTINE cp_run_input

RECURSIVE SUBROUTINE cp_run_input_comm(input_file_path,output_file_path,&
     mpi_comm,ierr)
  CHARACTER(len=*) :: input_file_path, output_file_path
  INTEGER :: mpi_comm, ierr
END SUBROUTINE cp_run_input_comm

RECURSIVE SUBROUTINE cp_rep_init(rep_env_id,ierr)
  INTEGER :: rep_env_id,ierr
END SUBROUTINE cp_rep_init

RECURSIVE SUBROUTINE cp_rep_calc_e_f(rep_env_id,calc_f,ierr)
  INTEGER :: rep_env_id,calc_f,ierr
END SUBROUTINE cp_rep_calc_e_f

RECURSIVE SUBROUTINE cp_ep_init(ep_env_id,ierr)
  INTEGER :: ep_env_id, ierr
END SUBROUTINE cp_ep_init

RECURSIVE SUBROUTINE cp_ep_calc_e_f(ep_env_id,calc_f,ierr)
  INTEGER :: ep_env_id, calc_f, ierr
END SUBROUTINE cp_ep_calc_e_f

SUBROUTINE cp_do_shake(f_env_id,dt,shake_tol,ierr)
  use kinds, only: dp
  INTEGER :: f_env_id
  REAL(kind=dp) :: dt, shake_tol
  integer :: ierr
END SUBROUTINE cp_do_shake

END INTERFACE
